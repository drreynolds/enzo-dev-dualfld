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
/  Initialization routine
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Allocates all necessary internal memory for problem 
/           definition and associated linear solver.  This begins the
/           interface between Enzo and the FLD solver module, so any
/           and all grid/index transformations must be performed and 
/           stored here.
/
************************************************************************/
#ifdef TRANSFER
#include "DualFLD.h"
#include "CosmologyParameters.h"
// #ifdef _OPENMP
// #include <omp.h>
// #endif

// function prototypes
int InitializeRateData(FLOAT Time);
int FreezeRateData(FLOAT Time, HierarchyEntry &TopGrid);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

// Function prototypes
int DualCosmoIonizationInitialize(FILE *fptr, FILE *Outfptr,
				  HierarchyEntry &TopGrid,
				  TopGridData &MetaData, int local);
int DualRHIonizationTestInitialize(FILE *fptr, FILE *Outfptr,
				   HierarchyEntry &TopGrid,
				   TopGridData &MetaData, int local);
int DualRHIonizationClumpInitialize(FILE *fptr, FILE *Outfptr,
				    HierarchyEntry &TopGrid,
				    TopGridData &MetaData, int local);
int DualRHIonizationSteepInitialize(FILE *fptr, FILE *Outfptr,
				    HierarchyEntry &TopGrid,
				    TopGridData &MetaData, int local);
int DualRadStreamTestInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);
int DualRadConstTestInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid,
			       TopGridData &MetaData, int local);

// character strings
EXTERN char outfilename[];



int DualFLD::Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData) {

  if (debug)  printf("Entering DualFLD::Initialize routine\n");


  // find the grid corresponding to this process from the Hierarcy
  HierarchyEntry *ThisGrid = &TopGrid;
  int i, dim, face, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) {
    ENZO_VFAIL("DualFLD::Initialize ERROR: p%"ISYM" could not locate his grid\n",
	       MyProcessorNumber)
  }

  
  //// Set parallelism information ////

// #ifdef _OPENMP
//   // output number of OpenMP threads that will be used in this run
//   int nthreads = omp_get_max_threads();
//   if (debug)
//     printf("DualFLD::Initialize: tasks have %"ISYM" available OpenMP threads\n",
//            nthreads);
// #endif


#ifndef MPI_INT
  // in case MPI is not included
  int MPI_PROC_NULL = -3;
  int MPI_COMM_WORLD = 0;
#endif

  //// Set local grid information ////

  // store problem rank 
  rank = MetaData.TopGridRank;

  // get processor layout from Grid
  for (dim=0; dim<rank; dim++) 
    layout[dim] = ThisGrid->GridData->GetProcessorLayout(dim);
  
  // get processor location in MPI grid
  for (dim=0; dim<rank; dim++) 
    location[dim] = ThisGrid->GridData->GetProcessorLocation(dim);

  // get neighbor information from grid
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) 
      NBors[dim][face] = ThisGrid->GridData->GetProcessorNeighbors(dim,face);

  // ensure that new BdryVals array pointers are set to NULL
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) {
      if (XrBdryVals[dim][face] != NULL) 
	delete[] XrBdryVals[dim][face];
      XrBdryVals[dim][face] = NULL;
      if (UVBdryVals[dim][face] != NULL) 
	delete[] UVBdryVals[dim][face];
      UVBdryVals[dim][face] = NULL;
    }

  // set up subdomain information
  //   EdgeVals gives the location of the left/right edge of the
  //      domain (no bdry) -- start with Enzo grid size
  for (dim=0; dim<rank; dim++) {
    EdgeVals[dim][0] = ThisGrid->GridData->GetGridLeftEdge(dim);
    EdgeVals[dim][1] = ThisGrid->GridData->GetGridRightEdge(dim);
  }

  //   LocDims holds the dimensions of the local domain, 
  //   active cells only (no ghost or boundary cells)
  for (dim=0; dim<rank; dim++)
    LocDims[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
      - ThisGrid->GridData->GetGridStartIndex(dim) + 1;

  // set flags denoting if this processor is on the external boundary
  for (dim=0; dim<rank; dim++) {
    if (layout[dim]==0) {
      OnBdry[dim][0] = OnBdry[dim][1] = true;
    }
    else {
      OnBdry[dim][0] = (location[dim] == 0);
      OnBdry[dim][1] = (location[dim] == layout[dim]-1);
    }
  }
  if (debug){
    printf("DualFLD::Initialize p%"ISYM": rank = %"ISYM"\n", MyProcessorNumber, rank);
    printf("DualFLD::Initialize p%"ISYM": layout = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,layout[0],layout[1],layout[2]);
    printf("DualFLD::Initialize p%"ISYM": location = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,location[0],location[1],location[2]);
  }

  // compute global dimension information
  for (dim=0; dim<rank; dim++)
    GlobDims[dim] = MetaData.TopGridDims[dim];

  // dx gives grid cell size (comoving, normalized units)
  for (dim=0; dim<rank; dim++)
    dx[dim] = (EdgeVals[dim][1]-EdgeVals[dim][0])/LocDims[dim];

  // compute global index information for this subdomain
  float fCellsLeft;
  for (dim=0; dim<rank; dim++) {

    // the global indexing is easy if we're at the left edge
    if (location[dim]==0)  SolvIndices[dim][0]=0;

    // otherwise we compute the number of intervening cells to left edge
    else {

      // get floating point value for number of cells
      fCellsLeft = (EdgeVals[dim][0] - DomainLeftEdge[dim])/dx[dim];

      // round floating point value to closest integer
      SolvIndices[dim][0] =  (long) (fCellsLeft >= 0.0) ?
	(trunc(fCellsLeft+0.5)) : (trunc(fCellsLeft-0.5));
    }

    // add on local size to obtain right edge indices
    SolvIndices[dim][1] = SolvIndices[dim][0] + LocDims[dim] - 1;
  }

  // store local array sizes (active + ghost)
  for (dim=0; dim<rank; dim++)
    ArrDims[dim] = LocDims[dim] + 2*NumberOfGhostZones;

  // set up ghost zone index information
  GhDims[0][0] = NumberOfGhostZones;
  GhDims[0][1] = NumberOfGhostZones;
  GhDims[1][0] = 0;
  GhDims[1][1] = 0;
  if (rank > 1) {
    GhDims[1][0] = NumberOfGhostZones;
    GhDims[1][1] = NumberOfGhostZones;
  }
  GhDims[2][0] = 0;
  GhDims[2][1] = 0;
  if (rank > 2) {
    GhDims[2][0] = NumberOfGhostZones;
    GhDims[2][1] = NumberOfGhostZones;
  }


  //// input/check solver parameters ////
  
  if (this->ReadParameters(MetaData, ThisGrid) != SUCCESS)
    ENZO_FAIL("DualFLD::Initialize: Error in ReadParameters.");


  //// General setup ////

  // a, adot give cosmological expansion & rate
  a = 1.0;  a0 = 1.0;  adot = 0.0;  adot0 = 0.0;

  //   for non-periodic domain, unset neighbor info.
  for (dim=0; dim<rank; dim++) {
    if ((OnBdry[dim][0]) && (XrBdryType[dim][0] != 0))
      NBors[dim][0] = MPI_PROC_NULL;
    if ((OnBdry[dim][1]) && (XrBdryType[dim][1] != 0))
      NBors[dim][1] = MPI_PROC_NULL;
  }

  // initialize the time values
  tnew = told = MetaData.Time;

  // dt* gives the time step sizes for each piece of physics
  dtrad  = initdt;                        // use the input value (scaled units)

  // set the global initial dt to match the radiation timestep
  dt = initdt;
  ThisGrid->GridData->SetMaxRadiationDt(dt);
  
  // get the current units values (used to help set the time step size)
  double MassUnits;
  float TempUnits;
  DenUnits=LenUnits=TempUnits=MassUnits=TimeUnits=VelUnits=aUnits=1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, 
	       &TimeUnits, &VelUnits, &MassUnits, MetaData.Time) == FAIL) 
    ENZO_FAIL("Error in GetUnits.");
  a = 1.0; adot = 0.0;
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(MetaData.Time, &a, &adot) == FAIL) 
      ENZO_FAIL("Error in CosmologyComputeExpansionFactor.");
    aUnits = 1.0/(1.0 + InitialRedshift);
  }

  // copy initial units values into "old" units
  a0 = a;
  adot0 = adot;
  NiUnits0 = NiUnits;
  DenUnits0 = DenUnits;
  LenUnits0 = LenUnits;

  // set up EnzoVector container for updated time step solution
  int nrad=0;
  if (UseXray) nrad++;
  if (UseUV)   nrad++;
  sol = new EnzoVector(LocDims[0], LocDims[1], LocDims[2], GhDims[0][0], 
		       GhDims[0][1], GhDims[1][0], GhDims[1][1], GhDims[2][0], 
		       GhDims[2][1], nrad, NBors[0][0], NBors[0][1], NBors[1][0], 
		       NBors[1][1], NBors[2][0], NBors[2][1], 1);

  // ensure that CoolData object has been set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) 
      ENZO_FAIL("DualFLD::Initialize Error in InitializeRateData.");

  // compute Radiation Energy spectrum integrals
  if (this->ComputeRadiationIntegrals() == FAIL) 
    ENZO_FAIL("DualFLD::Initialize Error in computing radiation spectrum integrals");

#ifdef USE_HYPRE

#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif
  // initialize HYPRE stuff
  //    initialize the diagnostic information
  totIters_Xr = 0;
  totIters_UV = 0;

  //    set up the grid
  //       create the grid object
  HYPRE_StructGridCreate(MPI_COMM_WORLD, rank, &grid);

  //       set my grid extents as if we have one part with multiple boxes.
  //       Have each processor describe it's own global extents
  Eint32 ilower[3] = {SolvIndices[0][0], SolvIndices[1][0], SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1], SolvIndices[1][1], SolvIndices[2][1]};
  HYPRE_StructGridSetExtents(grid, ilower, iupper);

  //       set grid periodicity
  Eint32 periodicity[3] = {0, 0, 0};
  if (XrBdryType[0][0] == 0)  periodicity[0] = GlobDims[0];
  if (XrBdryType[1][0] == 0)  periodicity[1] = GlobDims[1];
  if (XrBdryType[2][0] == 0)  periodicity[2] = GlobDims[2];
  HYPRE_StructGridSetPeriodic(grid, periodicity);
  
  //       assemble the grid
  HYPRE_StructGridAssemble(grid);

  //   set up the stencil
  if (rank == 1) 
    stSize = 3;
  else if (rank == 2)
    stSize = 5;
  else 
    stSize = 7;
  HYPRE_StructStencilCreate(rank, stSize, &stencil);

  //      set stencil entries
  Eint32 offset[3];
  Eint32 stentry=0;
  //         dependency to x2 left
  if (rank == 3) {
    offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //         dependency to x1 left
  if (rank >= 2) {
    offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //         dependency to x0 left
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //         dependency to self
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //         dependency to x0 right
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //         dependency to x1 right
  if (rank >= 2) {
    offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //         dependency to x2 right
  if (rank == 3) {
    offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }

  //   allocate temporary arrays
  int Nx = (SolvIndices[0][1]-SolvIndices[0][0]+1);
  int Ny = (SolvIndices[1][1]-SolvIndices[1][0]+1);
  int Nz = (SolvIndices[2][1]-SolvIndices[2][0]+1);
  HYPREbuff = new Eflt64[Nx];
  HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &P);
  HYPRE_StructMatrixInitialize(P);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &rhsvec);
  HYPRE_StructVectorInitialize(rhsvec);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &solvec);
  HYPRE_StructVectorInitialize(solvec);

#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  HYPREtime += ftime-stime;

#else  // ifdef USE_HYPRE

  ENZO_FAIL("DualFLD::Initialize ERROR: module requires USE_HYPRE to be set!");
  
#endif

  //   check HYPRE solver parameters
  if (sol_MGmaxit_Xr < 0) {
    fprintf(stderr,"Illegal DualFLDMaxMGItersXray = %i. Setting to 5\n",
	    sol_MGmaxit_Xr);
    sol_MGmaxit_Xr = 5;
  }
  if (sol_MGmaxit_UV < 0) {
    fprintf(stderr,"Illegal DualFLDMaxMGItersUV = %i. Setting to 3\n",
	    sol_MGmaxit_UV);
    sol_MGmaxit_UV = 3;
  }
  if (sol_KryMaxit_Xr < 0) {
    fprintf(stderr,"Illegal DualFLDMaxKryItersXray = %i. Setting to 3\n",
	    sol_KryMaxit_Xr);
    sol_KryMaxit_Xr = 3;
  }
  if (sol_KryMaxit_UV < 0) {
    fprintf(stderr,"Illegal DualFLDMaxKryItersUV = %i. Setting to 2\n",
	    sol_KryMaxit_UV);
    sol_KryMaxit_UV = 2;
  }
  if ((sol_rlxtype_Xr<0) || (sol_rlxtype_Xr>3)) {
    fprintf(stderr,"Illegal DualFLDMGRelaxTypeXray = %i. Setting to 2\n",
	    sol_rlxtype_Xr);
    sol_rlxtype_Xr = 2;
  }
  if ((sol_rlxtype_UV<0) || (sol_rlxtype_UV>3)) {
    fprintf(stderr,"Illegal DualFLDMGRelaxTypeUV = %i. Setting to 1\n",
	    sol_rlxtype_UV);
    sol_rlxtype_UV = 1;
  }
  if (sol_npre_Xr < 1) {
    fprintf(stderr,"Illegal DualFLDMGPreRelaxXray = %i. Setting to 3\n",
	    sol_npre_Xr);
    sol_npre_Xr = 3;
  }
  if (sol_npre_UV < 1) {
    fprintf(stderr,"Illegal DualFLDMGPreRelaxUV = %i. Setting to 1\n",
	    sol_npre_UV);
    sol_npre_UV = 1;
  }
  if (sol_npost_Xr < 1) {
    fprintf(stderr,"Illegal DualFLDMGPostRelaxXray = %i. Setting to 3\n",
	    sol_npost_Xr);
    sol_npost_Xr = 3;
  }
  if (sol_npost_UV < 1) {
    fprintf(stderr,"Illegal DualFLDMGPostRelaxUV = %i. Setting to 1\n",
	    sol_npost_UV);
    sol_npost_UV = 1;
  }
  if ((sol_tolerance_Xr < 1.0e-15) || (sol_tolerance_Xr > 1.0)) {
    fprintf(stderr,"Illegal DualFLDSolToleranceXray = %g. Setting to 1e-4\n",
	    sol_tolerance_Xr);
    sol_tolerance_Xr = 1.0e-4;
  }
  if ((sol_tolerance_UV < 1.0e-15) || (sol_tolerance_UV > 1.0)) {
    fprintf(stderr,"Illegal DualFLDSolToleranceUV = %g. Setting to 1e-4\n",
	    sol_tolerance_UV);
    sol_tolerance_UV = 1.0e-4;
  }
  if (sol_Krylov_Xr < 0 || sol_Krylov_Xr > 2) {
    fprintf(stderr,"Illegal DualFLDKrylovMethodXray = %"ISYM". Setting to 1 (BiCGStab)\n",
	    sol_Krylov_Xr);
    sol_Krylov_Xr = 1;
  }
  if (sol_Krylov_UV < 0 || sol_Krylov_UV > 2) {
    fprintf(stderr,"Illegal DualFLDKrylovMethodUV = %"ISYM". Setting to 1 (BiCGStab)\n",
	    sol_Krylov_UV);
    sol_Krylov_UV = 1;
  }


  ////////////////////////////////
  // set up the boundary conditions on the radiation field, 
  // depending on the ProblemType
  float ZERO = 0.0;
  float ONE  = 1.0;
  FILE *fptr = NULL;

  // set boundary conditions based on problem type
  // (default to homogeneous Dirichlet)
  switch (ProblemType) {
    
  // Ionization tests 0,1,7,8: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 410:
  case 411:
  case 417:
  case 418:
    // first call local problem initializer (to allocate/setup local data)
    if (DualRHIonizationTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("DualFLD::Initialize Error in RHIonizationTestInitialize.");
    
    // set BC on all non-periodic faces to zero-valued (Dirichlet or Neumann)
    for (dim=0; dim<MetaData.TopGridRank; dim++) {
      for (face=0; face<2; face++) {
	if (XrBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,0,1,&ZERO) == FAIL) {
	    ENZO_VFAIL("Error setting dim %"ISYM", face %"ISYM" Xray radiation BCs\n", dim, face)
	  }
	if (UVBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,1,1,&ZERO) == FAIL) {
	    ENZO_VFAIL("Error setting dim %"ISYM", face %"ISYM" UV radiation BCs\n", dim, face)
	  }
      }
    }
    break;
    
    
  // Ionization test 2: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all non-periodic faces.
  case 412:
    // first call local problem initializer (to allocate/setup local data)
    if (DualRHIonizationClumpInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("DualFLD::Initialize Error in RHIonizationClumpInitialize.");
    
    // set BC on all non-periodic faces to zero-valued (Dirichlet or Neumann)
    for (dim=0; dim<MetaData.TopGridRank; dim++) {
      for (face=0; face<2; face++) {
	if (XrBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,0,1,&ZERO) == FAIL) {
	    ENZO_VFAIL("Error setting dim %"ISYM", face %"ISYM" Xray radiation BCs\n", dim, face)
	  }
	if (UVBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,1,1,&ZERO) == FAIL) {
	    ENZO_VFAIL("Error setting dim %"ISYM", face %"ISYM" UV radiation BCs\n", dim, face)
	  }
      }
    }
    break;
    
    
  // Ionization test 13: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 413:
    // first call local problem initializer (to allocate/setup local data)
    if (DualRHIonizationSteepInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("DualFLD::Initialize Error in RHIonizationSteepInitialize.");
    
    // set BC on all non-periodic faces to zero-valued (Dirichlet or Neumann)
    for (dim=0; dim<MetaData.TopGridRank; dim++) {
      for (face=0; face<2; face++) {
	if (XrBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,0,1,&ZERO) == FAIL) {
	    ENZO_VFAIL("Error setting dim %"ISYM", face %"ISYM" Xray radiation BCs\n", dim, face)
	  }
	if (UVBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,1,1,&ZERO) == FAIL) {
	    ENZO_VFAIL("Error setting dim %"ISYM", face %"ISYM" UV radiation BCs\n", dim, face)
	  }
      }
    }
    break;
    
    
    
  // Ionization test 14: periodic boundary conditions on all faces (store no data).
  case 414:
    // first call local problem initializer (to allocate/setup local data)
    if (DualCosmoIonizationInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("DualFLD::Initialize Error in CosmoIonizationInitialize.");
    
    break;
    


  // Homogeneous test initializer
  case 416:
    // first call local problem initializer (to allocate/setup local data)
    if (DualRadConstTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("DualFLD::Initialize Error in DualRadHydroConstTestInitialize.");

    // set BC on all non-periodic faces to zero-valued (Dirichlet or Neumann)
    for (dim=0; dim<MetaData.TopGridRank; dim++) {
      for (face=0; face<2; face++) {
	if (XrBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,0,1,&ZERO) == FAIL) {
	    ENZO_VFAIL("DualFLD::Initialize Error setting dim %"ISYM
		       ", face %"ISYM" Xray radiation BCs\n", dim, face)
	  }
	if (UVBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,1,1,&ZERO) == FAIL) {
	    ENZO_VFAIL("DualFLD::Initialize Error setting dim %"ISYM
		       ", face %"ISYM" UV radiation BCs\n", dim, face)
	  }
      }
    }
    break;
    

    
  // Insert new problem intializers here...


  default:

    // set BC on all non-periodic faces to zero-valued (Dirichlet or Neumann)
    for (dim=0; dim<MetaData.TopGridRank; dim++) {
      for (face=0; face<2; face++) {
	if (XrBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,0,1,&ZERO) == FAIL) {
	    ENZO_VFAIL("DualFLD::Initialize Error setting dim %"ISYM
		       ", face %"ISYM" Xray radiation BCs\n", dim, face)
	  }
	if (UVBdryType[dim][face] != 0)
	  if (this->SetupBoundary(dim,face,1,1,&ZERO) == FAIL) {
	    ENZO_VFAIL("DualFLD::Initialize Error setting dim %"ISYM
		       ", face %"ISYM" UV radiation BCs\n", dim, face)
	  }
      }
    }
    break;

  }
  ////////////////////////////////

  // if using an isothermal model, freeze rate data, now that ICs exist
  if (isothermal) 
    if (FreezeRateData(MetaData.Time, TopGrid) == FAIL) 
      ENZO_FAIL("DualFLD::Initialize Error in FreezeRateData.");


  if (debug)  printf("DualFLD::Initialize: outputting parameters to log file\n");

  // initialize rate exchange fields to zero
  if (RadiativeCooling) {
    float *phHI       = ThisGrid->GridData->AccessKPhHI();
    float *phHeI      = ThisGrid->GridData->AccessKPhHeI();
    float *phHeII     = ThisGrid->GridData->AccessKPhHeII();
    float *photogamma = ThisGrid->GridData->AccessPhotoGamma();
    float *dissH2I    = ThisGrid->GridData->AccessKDissH2I();
    int i, size=1;
    for (int dim=0; dim<rank; dim++)  size *= ArrDims[dim];
    for (i=0; i<size; i++)  phHI[i] = 0.0;
    if (RadiativeTransferHydrogenOnly == FALSE) {
      for (i=0; i<size; i++)  phHeI[i]  = 0.0;
      for (i=0; i<size; i++)  phHeII[i] = 0.0;
    }
    for (i=0; i<size; i++)  photogamma[i] = 0.0;
    if (MultiSpecies > 1) 
      for (i=0; i<size; i++)  dissH2I[i] = 0.0;
  }

  ////////////////////////////////

  // output RadHydro solver parameters to output log file 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *outfptr;
    if ((outfptr = fopen(outfilename, "a")) == NULL) {
      ENZO_VFAIL("DualFLD::Initialize: Error opening parameter output file %s!!\n", outfilename)
    }
    else {
      // write parameters to log file and close
      this->WriteParameters(outfptr);
      fclose(outfptr);
    }
  }

  return SUCCESS;
}
#endif
