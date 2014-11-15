/*****************************************************************************
 *                                                                           *
 * Copyright 2014 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Dual Group (Gray UV + X-ray) Flux-Limited Diffusion solver class
/  Evolve Routine 
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Performs a time step update of the UV and/or Xray radiation 
/           fields.  Evolution is performed in an operator-split fashion:
/              We take radiation time steps to catch up with the 
/                hydrodynamic time, i.e. (dt_rad <= dt_hydro). 
/              Both the UV and X-ray radiation are computed with the same
/                step size.  
/              Since chemistry and gas cooling are handled externally to
/                this module, we generally prefer (dt_rad == dt_hydro), 
/                but subcycling is allowed for robustness when handling
/                fast transients (e.g. the first step after a source 
/                turns on).  
/              As a result, prior to completion, this routine updates the
/                maximum time step the overall Grid module can take so 
/                that subsequent steps are taken with the same sizes.
/
************************************************************************/
#ifdef TRANSFER
 
#include "DualFLD.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"


//#define FAIL_ON_NAN
#define NO_FAIL_ON_NAN

#define SOLUTION_DIAGNOSTICS


/* function prototypes */
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int RadiationGetUnits(float *RadiationUnits, FLOAT Time);



int DualFLD::Evolve(HierarchyEntry *ThisGrid, float dthydro)
{

  // Only continue if we own this grid
  if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber())
    return SUCCESS;

#ifdef USE_MPI
  //  check that MyProcessorNumber agrees with MPI process ID
  MPI_Arg MPI_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_id);
  if (MyProcessorNumber != MPI_id) {
    ENZO_VFAIL("ERROR: Enzo PID %"ISYM" doesn't match MPI ID %"ISYM"\n", MyProcessorNumber, int(MPI_id))
  }
#endif

  // in case MPI is not included
#ifndef MPI_INT
  int MPI_COMM_WORLD = 0;
#endif

  // start MPI timer for overall solver
#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif

  ////////////////////////////////////
  // Problem Setup Phase

  if (debug)  printf("\n DualFLD::Evolve:\n");

  // Set pointers to relevant fields
  float *XrRadiation = ThisGrid->GridData->AccessRadiationFrequency0();
  float *UVRadiation = ThisGrid->GridData->AccessRadiationFrequency1();
  if (UseXray && XrRadiation == NULL) 
    ENZO_FAIL("DualFLD::Evolve: could not obtain X-ray radiation field");
  if (UseUV && UVRadiation == NULL) 
    ENZO_FAIL("DualFLD::Evolve: could not obtain UV radiation field");

  // Get current time from Grid
  tnew = ThisGrid->GridData->ReturnTime();

  // set shortcut for size of BaryonFields
  int i, dim, size=1;
  for (dim=0; dim<rank; dim++)  size *= ArrDims[dim];

  // point U EnzoVector fields at existing radiation fields
  // (for ghost zone communication)
  int nrad=0;
  if (UseXray)  U->SetData(nrad++, XrRadiation);
  if (UseUV)    U->SetData(nrad++, UVRadiation);

  // rescale radiation fields to non-dimensionalize within solver
  nrad=0;
  if (UseXray)  U->scale_component(nrad++, 1.0/XrScale);
  if (UseUV)    U->scale_component(nrad++, 1.0/UVScale);

  // enforce boundary conditions on current radiation fields
  if (EnforceBoundary(ThisGrid) != SUCCESS)
    ENZO_FAIL("DualFLD::Evolve: EnforceBoundary failure!!");
  
  // copy current radiation fields into solution vector
  sol->copy(U);
  
  // begin communication of ghost zone information
  if (sol->exchange_start() != SUCCESS) 
    ENZO_FAIL("DualFLD::Evolve: vector exchange_start error");

  // update internal Enzo units for current times
  float TempUnits, RadUnits;
  double MassUnits;
  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, told) != SUCCESS) 
    ENZO_FAIL("DualFLD::Evolve: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, told) != SUCCESS) 
    ENZO_FAIL("DualFLD::Evolve: Error in RadiationGetUnits.");
  UVUnits0 = RadUnits*UVScale;
  XrUnits0 = RadUnits*XrScale;
  NiUnits0 = DenUnits/mh;

  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, tnew) != SUCCESS) 
    ENZO_FAIL("DualFLD::Evolve: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, tnew) != SUCCESS) 
    ENZO_FAIL("DualFLD::Evolve: Error in RadiationGetUnits.");
  UVUnits = RadUnits*UVScale;
  XrUnits = RadUnits*XrScale;
  NiUnits = DenUnits/mh;

  // output typical/maximum values
  float UTypVals[2], UMaxVals[2];
#ifdef SOLUTION_DIAGNOSTICS
  nrad=0;
  if (UseXray) {
    UTypVals[0] = U->rmsnorm_component(nrad);
    UMaxVals[0] = U->infnorm_component(nrad++);
  }
  if (UseUV) {
    UTypVals[1] = U->rmsnorm_component(nrad);
    UMaxVals[1] = U->infnorm_component(nrad++);
  }

  if (debug) {
    printf("   current internal (physical) quantities:\n");
    if (UseXray)
      printf("    Xray rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	     UTypVals[0], UTypVals[0]*XrUnits, 
	     UMaxVals[0], UMaxVals[0]*XrUnits);
    if (UseUV)
      printf("      UV rms = %13.7e (%8.2e), max = %13.7e (%8.2e)\n",
	     UTypVals[1], UTypVals[1]*UVUnits, 
	     UMaxVals[1], UMaxVals[1]*UVUnits);
  }
#endif

  // if autoScale enabled, determine scaling factor updates here
  float ScaleCorrTol = 1.e-2;
  float XrScaleCorr = 1.0;
  float UVScaleCorr = 1.0;
  if (StartAutoScale && autoScale) {
#ifndef SOLUTION_DIAGNOSTICS
    nrad=0;
    if (UseXray) {
      UTypVals[0] = U->rmsnorm_component(nrad);
      UMaxVals[0] = U->infnorm_component(nrad++);
    }
    if (UseUV) {
      UTypVals[1] = U->rmsnorm_component(nrad);
      UMaxVals[1] = U->infnorm_component(nrad++);
    }
#endif
    if ((UMaxVals[0] - UTypVals[0]) > ScaleCorrTol*UMaxVals[0])
      XrScaleCorr = UMaxVals[0];
    if ((UMaxVals[1] - UTypVals[1]) > ScaleCorrTol*UMaxVals[1])
      UVScaleCorr = UMaxVals[1];
  }

  // finish communication of neighbor information
  if (sol->exchange_end() != SUCCESS) 
    ENZO_FAIL("DualFLD::Evolve: vector exchange_end error\n");


  ////////////////////////////////////
  // Problem Solve Phase

  // internal time-stepping loop to catch up with Hydro time
  float stime2, ftime2;   // radiation timer
  int radstep, radstop;   // subcycle iterators
  float end_time = tnew + dthydro;
  radstop = 0;
  int maxtries = 100;
  for (radstep=0; radstep<maxtries; radstep++) {
      
    // start MPI timer for radiation solver
#ifdef USE_MPI
    stime2 = MPI_Wtime();
#else
    stime2 = 0.0;
#endif

    // update time-step information
    told = tnew;

    // keep trying time steps until radiation solver succeeds. 
    // Note: if we reach the minimum time step size, RadStep will call ENZO_FAIL
    int recompute_step = 0;
    do {

      // update time-step information.  Note: dtrad was set on previous 
      // iteration of solver, or by user input for first iteration
      tnew = told + dtrad;
      if ((tnew - end_time)/end_time > -1.0e-14) {   // do not exceed synchronization time
	tnew = end_time;
	radstop = 1;
      }
      dt = tnew - told;
      if (debug) {
	printf("\n rad attempt %"ISYM": dt=%7.1e, t=%7.1e (hydro dt=%7.1e, t=%7.1e)\n",
	       radstep, dt, tnew, dthydro, end_time);
	printf(" ----------------------------------------------------------------------\n");
      }

      // Xray evolution
      if (UseXray) {
	if (!XrStatic) {
	  if (debug)  printf("  Xray:");
	  recompute_step = this->RadStep(ThisGrid, 0);
	}
      }

      // UV evolution
      if (!recompute_step && UseUV) {
	if (!UVStatic) {
	  if (debug)  printf("  UV:");
	  recompute_step = this->RadStep(ThisGrid, 1);
	}
      }
      
      if (debug)  
	printf(" ======================================================================\n\n");

      // if either radiation step was unsuccessful, back-track to previous 
      // step and pull back on dtrad
      if (recompute_step) {
	dtrad = max(dtrad*0.5, mindt);
	tnew = told;
	radstop = 0;
      }

    } while (recompute_step);
    
    // stop MPI timer for radiation solver, increment total
#ifdef USE_MPI
    ftime2 = MPI_Wtime();
#else
    ftime2 = 0.0;
#endif
    HYPREtime += ftime2-stime2;
    
    // update the radiation time step size for next time step
    //   (limit growth at each cycle)
    float dt_est = this->ComputeTimeStep(U, sol);
    dtrad = min(dt_est, dtgrowth*dtrad);

    // update Enzo radiation fields with new values
    U->copy(sol);
    
    // // have U communicate neighbor information
    // if (U->exchange() != SUCCESS) 
    //   ENZO_FAIL("DualFLD::Evolve: vector exchange error");

    // break out of time-stepping loop if we've reached the end
    if (radstop)  break;
	
  } // end outer radiation time-stepping loop


  ////////////////////////////////////
  // Problem Cleanup and Preparation for Next Call Phase

  // if chemistry and cooling handled elsewhere, fill the rates
  if (RadiativeCooling) 
    this->FillRates(ThisGrid);

#ifdef SOLUTION_DIAGNOSTICS
  // output the total photo-heating/photo-ionization rates
  if (debug)  printf("  Photo-heating/photo-ionization rates:\n");

  float *phHI = ThisGrid->GridData->AccessKPhHI();
  EtaVec->SetData(0, phHI);
  float phHI_norm = EtaVec->rmsnorm();
  float phHI_max  = EtaVec->infnorm();
  if (debug)  printf("          phHI norm = %.2e,  max = %.2e\n", 
		     phHI_norm, phHI_max);

  float *photogamma = ThisGrid->GridData->AccessPhotoGamma();
  EtaVec->SetData(0, photogamma);
  float photogamma_norm = EtaVec->rmsnorm();
  float photogamma_max  = EtaVec->infnorm();
  if (debug)  printf("    photogamma norm = %.2e,  max = %.2e\n", 
		     photogamma_norm, photogamma_max);

  if (!RadiativeTransferHydrogenOnly) {
    float *phHeI = ThisGrid->GridData->AccessKPhHeI();
    EtaVec->SetData(0, phHeI);
    float phHeI_norm = EtaVec->rmsnorm();
    float phHeI_max  = EtaVec->infnorm();
    if (debug)  printf("         phHeI norm = %.2e,  max = %.2e\n", 
		       phHeI_norm, phHeI_max);

    float *phHeII = ThisGrid->GridData->AccessKPhHeII();
    EtaVec->SetData(0, phHeII);
    float phHeII_norm = EtaVec->rmsnorm();
    float phHeII_max  = EtaVec->infnorm();
    if (debug)  printf("        phHeII norm = %.2e,  max = %.2e\n", 
		       phHeII_norm, phHeII_max);
  }
#endif

  // update the radiation time step size for next time step
  if (dtrad != huge_number)
    ThisGrid->GridData->SetMaxRadiationDt(dtrad);

  // scale back to Enzo units
  nrad=0;
  if (UseXray)  U->scale_component(nrad++, XrScale);
  if (UseUV)    U->scale_component(nrad++, UVScale);

  // update scaling factors to account for new values
  if (StartAutoScale && autoScale) {
    if (UseXray)  XrScale *= XrScaleCorr;
    if (UseUV)    UVScale *= UVScaleCorr;
  }

  // stop MPI timer, add to cumulative clock, output to stdout
#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  RTtime += ftime-stime;
  if (debug)  printf("RadHydro cumulative time = %g (HYPRE = %g)\n\n",
		     RTtime, HYPREtime);

  return SUCCESS;
}



// This routine evolves the radiation subsystem within the DualFLD module.  
// This is performed in a robust manner; if the current radiation subsystem 
// cannot be solved to the desired tolerance (meaning the time step is 
// likely much too large), it will return a flag to the calling routine to 
// have the step recomputed with a smaller time step size.
//
// Note:  XrUv=0 -> Xray solve;  XrUv=1 -> UV solve
int DualFLD::RadStep(HierarchyEntry *ThisGrid, int XrUv)
{
   
  // initialize return value
  int recompute_step = 0;

  // update internal Enzo units for current times
  float TempUnits, RadUnits;
  double MassUnits;
  DenUnits0 = LenUnits0 = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits0, &LenUnits0, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, told) != SUCCESS) 
    ENZO_FAIL("DualFLD_RadStep: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, told) != SUCCESS) 
    ENZO_FAIL("DualFLD_RadStep: Error in RadiationGetUnits.");
  UVUnits0 = RadUnits*UVScale;
  XrUnits0 = RadUnits*XrScale;
  NiUnits0 = DenUnits0/mh;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(told, &a0, &adot0) != SUCCESS) 
      ENZO_FAIL("DualFLD_RadStep: Error in CosmologyComputeExpansionFactor.");
  adot0 /= TimeUnits;  // rescale to physical units

  DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = 1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, &TimeUnits, 
	       &VelUnits, &MassUnits, tnew) != SUCCESS) 
    ENZO_FAIL("DualFLD_RadStep: Error in GetUnits.");
  if (RadiationGetUnits(&RadUnits, tnew) != SUCCESS) 
    ENZO_FAIL("DualFLD_RadStep: Error in RadiationGetUnits.");
  UVUnits = RadUnits*UVScale;
  XrUnits = RadUnits*XrScale;
  NiUnits = DenUnits/mh;
  if (ComovingCoordinates) 
    if (CosmologyComputeExpansionFactor(tnew, &a, &adot) != SUCCESS) 
      ENZO_FAIL("DualFLD_RadStep: Error in CosmologyComputeExpansionFactor.");
  adot /= TimeUnits;  // rescale to physical units
    
  // rescale dt, told, tnew, adot to physical values for use within solver
  dt   *= TimeUnits;
  told *= TimeUnits;
  tnew *= TimeUnits;


  // set shortcut pointers for opacity, emissivity, current solution
  // (kphHI and photogamma must exist, and are currently unused)
  float *Opacity = ThisGrid->GridData->AccessKPhHI();
  if (Opacity == NULL)
    ENZO_FAIL("DualFLD_RadStep: could not obtain KPhHI field");
  float *Eta = ThisGrid->GridData->AccessPhotoGamma();
  if (Eta == NULL)
    ENZO_FAIL("DualFLD_RadStep: could not obtain photogamma field");
  int rad_idx = (UseXray && UseUV && XrUv) ? 1 : 0;
  float *Enew = sol->GetData(rad_idx);

  // set shortcut for size of BaryonFields
  int i, dim, size=1;
  for (dim=0; dim<rank; dim++)  size *= ArrDims[dim];

  // initialize emissivity sources to 0
  for (i=0; i<size; i++)  Eta[i] = 0.0;

  // copy externally supplied emissivity field into Eta, if applicable
#ifdef EMISSIVITY
  if (StarMakerEmissivityField > 0) {
    float *EtaSM = ThisGrid->GridData->AccessEmissivity0();
    if (EtaSM == NULL)
      ENZO_FAIL("DualFLD_RadStep: could not obtain Emissivity0 field");
    float factor = (XrUv) ? 1.0 : Xr_parameter/UV_parameter;
    for (i=0; i<size; i++)  Eta[i] += (EtaSM[i] * factor);
  }
#endif

  // update emissivity field with locally-defined sources
  if (this->RadiationSource(ThisGrid, XrUv, Eta) != SUCCESS)
    ENZO_FAIL("DualFLD_RadStep: Error in RadiationSource routine\n");

  // output emissivity information, if applicable
  EtaVec->SetData(0, Eta);
  float eta_norm = EtaVec->rmsnorm();
  float eta_max  = EtaVec->infnorm();
  if (debug)  printf("  emissivity norm = %g,  max = %g\n", eta_norm, eta_max);

  // turn on automatic scaling for next step if this step has nontrivial emissivity
  float ScaleCorrTol = 1.e-2;
  if ((eta_max - eta_norm) > ScaleCorrTol*eta_max)
    StartAutoScale = true;

  // compute updated opacities
  if (this->ComputeOpacity(ThisGrid, XrUv, Opacity) != SUCCESS) 
    ENZO_FAIL("DualFLD_RadStep: Error in Opacity routine");
    
   
#ifdef USE_HYPRE
  // HYPRE stuff

  // set solver parameters based on Xray vs UV physics
  float sol_tolerance;
  Eint32 sol_MGmaxit, sol_KryMaxit, sol_rlxtype, sol_npre, sol_npost, Krylov_method;
  int *totIters;
  if (XrUv == 0) {   // Xray
    sol_tolerance = sol_tolerance_Xr;
    sol_MGmaxit   = sol_MGmaxit_Xr;
    sol_KryMaxit  = sol_KryMaxit_Xr;
    sol_rlxtype   = sol_rlxtype_Xr;
    sol_npre      = sol_npre_Xr;
    sol_npost     = sol_npost_Xr;
    Krylov_method = sol_Krylov_Xr;
    totIters      = &totIters_Xr;
  } else {           // UV
    sol_tolerance = sol_tolerance_UV;
    sol_MGmaxit   = sol_MGmaxit_UV;
    sol_KryMaxit  = sol_KryMaxit_UV;
    sol_rlxtype   = sol_rlxtype_UV;
    sol_npre      = sol_npre_UV;
    sol_npost     = sol_npost_UV;
    Krylov_method = sol_Krylov_UV;
    totIters      = &totIters_UV;
  }
  
  // set up radiation linear system
  float rhsnorm;      // used for setting HYPRE solver tolerance
  if (this->SetupSystem(ThisGrid, XrUv, Enew, Opacity, Eta, rhsnorm) != SUCCESS) 
    ENZO_FAIL("DualFLD_RadStep: Error in SetupSystem routine");


#define NO_PRINT_HYPRE_MATRICES
#ifdef PRINT_HYPRE_MATRICES
  // dump HYPRE matrices to disk
  HYPRE_StructMatrixPrint("A-hypre",P,0);
  HYPRE_StructVectorPrint("B-hypre",rhsvec,0);
#endif  
 


  // skip solve if ||rhs|| < sol_tolerance  (i.e. old solution is fine)
  //  if (rhsnorm < sol_tolerance) {
  if (rhsnorm < 0.01*sol_tolerance) {
    if (debug) {
      printf("   no solve required: |rhs| = %.1e  <  tol = %.1e\n", rhsnorm, sol_tolerance);
    }
    // rescale dt, told, tnew, adot back to normalized values
    dt   /= TimeUnits;
    told /= TimeUnits;
    tnew /= TimeUnits;
    return 0;
  }

  // set linear solver tolerance (rescale to relative residual and not actual)
  bool use_abs_resid = (rhsnorm < 1.e-8);
  Eflt64 delta = (use_abs_resid) ? sol_tolerance : sol_tolerance/rhsnorm;
  delta = min(delta, 1.0e-2);

  // set up the solver and preconditioner [PFMG]
  //    create the solver & preconditioner
  HYPRE_StructSolver solver;            // HYPRE solver structure
  HYPRE_StructSolver preconditioner;    // HYPRE preconditioner structure
  switch (Krylov_method) {
  case 0:   // PCG
    HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
    break;
  case 2:   // GMRES
    HYPRE_StructGMRESCreate(MPI_COMM_WORLD, &solver);
    break;
  default:  // BiCGStab
    HYPRE_StructBiCGSTABCreate(MPI_COMM_WORLD, &solver);
    break;
  }
  HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &preconditioner);
  
  // Multigrid solver: for periodic dims, only coarsen until grid no longer divisible by 2
  Eint32 max_levels, level=-1;
  int Ndir;
  if (XrBdryType[0][0] == 0) {
    level = 0;
    Ndir = GlobDims[0];
    while ( Ndir%2 == 0 ) {
      level++;
      Ndir /= 2;
    }
  }
  max_levels = level;
  if (rank > 1) {
    if (XrBdryType[1][0] == 0) {
      level = 0;
      Ndir = GlobDims[1];
      while ( Ndir%2 == 0 ) {
	level++;
	Ndir /= 2;
      }
    }
    max_levels = min(level,max_levels);
  }
  if (rank > 2) {
    if (XrBdryType[2][0] == 0) {
      level = 0;
      Ndir = GlobDims[2];
      while ( Ndir%2 == 0 ) {
	level++;
	Ndir /= 2;
      }
    }
    max_levels = min(level,max_levels);
  }

  //    set preconditioner options
  if (max_levels > -1) 
    HYPRE_StructPFMGSetMaxLevels(preconditioner, max_levels);
  HYPRE_StructPFMGSetMaxIter(preconditioner, sol_MGmaxit);
  HYPRE_StructPFMGSetRelaxType(preconditioner, sol_rlxtype);
  HYPRE_StructPFMGSetNumPreRelax(preconditioner, sol_npre);
  HYPRE_StructPFMGSetNumPostRelax(preconditioner, sol_npost);
  
  //    set solver options
  switch (Krylov_method) {
  case 0:   // PCG
    HYPRE_StructPCGSetPrintLevel(solver, sol_printl);
    HYPRE_StructPCGSetLogging(solver, sol_log);
    HYPRE_StructPCGSetRelChange(solver, 1);
    if (rank > 1) {
      HYPRE_StructPCGSetMaxIter(solver, sol_KryMaxit);
      HYPRE_StructPCGSetPrecond(solver, 
				(HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
				(HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
				preconditioner);
    } else {    // ignore preconditioner for 1D tests (bug); increase CG its
      HYPRE_StructPCGSetMaxIter(solver, sol_KryMaxit*500);
    }
    if (use_abs_resid) {
      HYPRE_StructPCGSetTol(solver, 0.0);
      HYPRE_StructPCGSetAbsoluteTol(solver, delta);
    } else {
      HYPRE_StructPCGSetTol(solver, delta);
    }
    HYPRE_StructPCGSetup(solver, P, rhsvec, solvec);
    break;
  case 2:   // GMRES
    //  HYPRE_StructGMRESSetPrintLevel(solver, sol_printl);
    HYPRE_StructGMRESSetLogging(solver, sol_log);
    //  HYPRE_StructGMRESSetRelChange(solver, 1);
    if (rank > 1) {
      HYPRE_StructGMRESSetMaxIter(solver, sol_KryMaxit);
      HYPRE_StructGMRESSetKDim(solver, sol_KryMaxit);
      HYPRE_StructGMRESSetPrecond(solver, 
				  (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
				  (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
				  preconditioner);
    } else {    // ignore preconditioner for 1D tests (bug); increase its
      HYPRE_StructGMRESSetMaxIter(solver, sol_KryMaxit*50);
      HYPRE_StructGMRESSetKDim(solver, sol_KryMaxit*50);
    }
    if (use_abs_resid) {
      HYPRE_StructGMRESSetTol(solver, 0.0);
      HYPRE_StructGMRESSetAbsoluteTol(solver, delta);
    } else {
      HYPRE_StructGMRESSetTol(solver, delta);
    }
    HYPRE_StructGMRESSetup(solver, P, rhsvec, solvec);
    break;
  default:  // BiCGStab
    //  HYPRE_StructBiCGSTABSetPrintLevel(solver, sol_printl);
    HYPRE_StructBiCGSTABSetLogging(solver, sol_log);
    if (rank > 1) {
      HYPRE_StructBiCGSTABSetMaxIter(solver, sol_KryMaxit);
      HYPRE_StructBiCGSTABSetPrecond(solver, 
				     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
				     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
				     preconditioner);
    } else {    // ignore preconditioner for 1D tests (bug); increase its
      HYPRE_StructBiCGSTABSetMaxIter(solver, sol_KryMaxit*500);
    }
    if (use_abs_resid) {
      HYPRE_StructBiCGSTABSetTol(solver, 0.0);
      HYPRE_StructBiCGSTABSetAbsoluteTol(solver, delta);
    } else {
      HYPRE_StructBiCGSTABSetTol(solver, delta);
    }
    HYPRE_StructBiCGSTABSetup(solver, P, rhsvec, solvec);
    break;
  }
  
  // solve the linear system
  switch (Krylov_method) {
  case 0:   // PCG
    HYPRE_StructPCGSolve(solver, P, rhsvec, solvec);
    break;
  case 2:   // GMRES
    HYPRE_StructGMRESSolve(solver, P, rhsvec, solvec);
    break;
  default:  // BiCGStab
    HYPRE_StructBiCGSTABSolve(solver, P, rhsvec, solvec);
    break;
  }
  
#ifdef PRINT_HYPRE_MATRICES
  // dump HYPRE matrices to disk
  HYPRE_StructVectorPrint("X-hypre",solvec,0);
#endif  

  // extract solver & preconditioner statistics
  Eflt64 finalresid=1.0;  // HYPRE solver statistics
  Eint32 Sits=0, Pits=0;  // HYPRE solver statistics
  switch (Krylov_method) {
  case 0:   // PCG
    if (use_abs_resid)
      finalresid = 0.0;
    else
      HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &finalresid);
    HYPRE_StructPCGGetNumIterations(solver, &Sits);
    break;
  case 2:   // GMRES
    if (use_abs_resid)
      finalresid = 0.0;
    else
      HYPRE_StructGMRESGetFinalRelativeResidualNorm(solver, &finalresid);
    HYPRE_StructGMRESGetNumIterations(solver, &Sits);
    break;
  default:  // BiCGStab
    if (use_abs_resid)
      finalresid = 0.0;
    else
      HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(solver, &finalresid);
    HYPRE_StructBiCGSTABGetNumIterations(solver, &Sits);
    break;
  }
  HYPRE_StructPFMGGetNumIterations(preconditioner, &Pits);
  *totIters += Sits;
  if (debug) printf("   lin resid = %.1e (tol = %.1e, |rhs| = %.1e), its = (%i,%i)\n",
		    finalresid*rhsnorm, delta, rhsnorm, Sits, Pits);

  // check if residual is NaN
  if (finalresid != finalresid) {

#ifndef FAIL_ON_NAN
    if (dt > mindt*TimeUnits*1.00000001) {
      // allow remainder of function to complete (to reset units, etc.), 
      // but have calling routine update dt and compute step again.
      recompute_step = 1;
    } else {
#endif
      fprintf(stderr,"DualFLD_RadStep: could not solve problem at minimum step size!\n");
      
      // output linear system to disk
      if (debug)  printf("Writing out matrix to file P.mat\n");
      HYPRE_StructMatrixPrint("P.mat",P,0);
      if (debug)  printf("Writing out rhs to file b.vec\n");
      HYPRE_StructVectorPrint("b.vec",rhsvec,0);
      if (debug)  printf("Writing out current solution to file x.vec\n");
      HYPRE_StructVectorPrint("x.vec",solvec,0);
      
      ENZO_FAIL("Error in DualFLD_RadStep");
#ifndef FAIL_ON_NAN
    }
#endif
  }
  
  // if solve was successful: extract values and add to current solution
  if (!recompute_step) {
    Eint32 ilower[3] = {SolvIndices[0][0],SolvIndices[1][0],SolvIndices[2][0]};
    Eint32 iupper[3] = {SolvIndices[0][1],SolvIndices[1][1],SolvIndices[2][1]};
    int xBuff = GhDims[0][0]-SolvOff[0];
    int yBuff = (GhDims[1][0]-SolvOff[1])-SolvIndices[1][0];
    int zBuff = (GhDims[2][0]-SolvOff[2])-SolvIndices[2][0];
    for (int iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
      int Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];
      ilower[2] = iz;  iupper[2] = iz;
      for (int iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
	int Ybl = (iy+yBuff)*ArrDims[0];
	ilower[1] = iy;  iupper [1] = iy;
	HYPRE_StructVectorGetBoxValues(solvec, ilower, iupper, HYPREbuff);
	for (int ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	  Enew[Zbl+Ybl+xBuff+ix] += HYPREbuff[ix];
      }
    }
  }

  // destroy HYPRE solver & preconditioner structures
  switch (Krylov_method) {
  case 0:   // PCG
    HYPRE_StructPCGDestroy(solver);
    break;
  case 2:   // GMRES
    HYPRE_StructGMRESDestroy(solver);
    break;
  default:  // BiCGStab
    HYPRE_StructBiCGSTABDestroy(solver);
    break;
  }
  HYPRE_StructPFMGDestroy(preconditioner);
  
  // enforce a solution floor on radiation
  float epsilon=1.0;      // radiation floor
  while (epsilon*0.25 > 0.0)  epsilon*=0.5;
  for (int i=0; i<ArrDims[0]*ArrDims[1]*ArrDims[2]; i++)  
    Enew[i] = max(Enew[i], epsilon);

  // rescale dt, told, tnew, adot back to normalized values
  dt   /= TimeUnits;
  told /= TimeUnits;
  tnew /= TimeUnits;

  return recompute_step;

#else
  ENZO_FAIL("DualFLD_RadStep ERROR: module requires USE_HYPRE to be set!");
#endif

}




#endif   // TRANSFER
