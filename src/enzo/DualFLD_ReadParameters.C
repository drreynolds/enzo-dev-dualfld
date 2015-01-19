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
/  ReadParameters routine
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Reads the DualFLD-specific solver parameters from the
/           input file.
/
************************************************************************/
#ifdef TRANSFER
#include "DualFLD.h"


int DualFLD::ReadParameters(TopGridData &MetaData, HierarchyEntry *ThisGrid) {

  if (debug)  printf("Entering DualFLD::ReadParameters routine\n");

  // local variables
  int dim, isrc, face;

  // set default module parameters
  int use_xray = 1;     // enable Xray propagation
  int use_uv   = 1;     // enable UV propagation
  int staticXr = 0;     // dynamic Xray radiation
  int staticUV = 0;     // dynamic UV radiation
  int xr_diff  = 0;     // normal Xray radiation equation
  UVSpectrum  = -1;     // monochromatic spectrum
  UVFrequency = 13.6;   // monochromatic spectrum frequency (eV)
  XrSpectrum  = -1;     // monochromatic spectrum
  XrFrequency = 500.0;  // monochromatic spectrum frequency (500 keV)
  theta  = 1.0;         // backwards euler implicit time discret.
  dtnorm = 2.0;         // use 2-norm for time step estimation
  dtgrowth = 1.1;       // 10% allowed growth in dt per step
  timeAccuracyXr = huge_number;  // no limit
  timeAccuracyUV = huge_number;  // no limit
  dt_control  = 2;      // PID controller
  UVScale     = 1.0;    // no radiation equation scaling
  XrScale     = 1.0;    // no radiation equation scaling
  int autoscale = 1;    // enable automatic variable scaling
  isothermal = 0;       // temperature-dependent chemistry
  for (dim=0; dim<rank; dim++)     // set default radiation boundaries to 
    for (face=0; face<2; face++) { //   periodic in each direction
      XrBdryType[dim][face] = 0;
      UVBdryType[dim][face] = 0;
    }

  // set default solver parameters
  sol_tolerance_Xr   = 1.0e-5;    // HYPRE solver tolerance
  sol_tolerance_UV   = 1.0e-5;
  sol_MGmaxit_Xr     = 5;         // HYPRE max multigrid iters
  sol_MGmaxit_UV     = 3;
  sol_KryMaxit_Xr    = 3;         // HYPRE max Krylov iters
  sol_KryMaxit_UV    = 2;
  sol_rlxtype_Xr     = 2;         // HYPRE relaxation type
  sol_rlxtype_UV     = 1;
  sol_npre_Xr        = 3;         // HYPRE num pre-smoothing steps
  sol_npre_UV        = 1;
  sol_npost_Xr       = 3;         // HYPRE num post-smoothing steps
  sol_npost_UV       = 1;
  sol_printl         = 1;         // HYPRE print level
  sol_log            = 1;         // HYPRE logging level
  sol_Krylov_Xr      = 1;         // BiCGStab outer solver
  sol_Krylov_UV      = 1;         // BiCGStab outer solver

  // set default limiter parameters
  FLD_Rmin = 1e-2;       // lower bound on R (normalized)
  FLD_Dmax = 1e-2;       // upper bound on D (normalized)
  FLD_Kmin = 1e-20;      // lower bound on opacity (normalized)
  FLD_Emin = 1e-30;      // lower bound on radiation (normalized)

  // set default chemistry constants
  hnu0_HI   = 13.6;      // ionization energy of HI   [eV]
  hnu0_HeI  = 24.6;      // ionization energy of HeI  [eV]
  hnu0_HeII = 54.4;      // ionization energy of HeII [eV]


  // if input file present, over-write defaults with module inputs
  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  int ret, ival;
  float fval1, fval2, fval3;
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;

  // check whether input file is non-null
  if (MetaData.RadHydroParameterFname != NULL) {
    if ((fptr = fopen(MetaData.RadHydroParameterFname, "r")) == NULL)
      fprintf(stderr,"Error opening RadHydro parameter file %s, using defaults\n",
	      MetaData.RadHydroParameterFname);
    else {

      // read until out of lines
      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
	ret = 0;
	ret += sscanf(line, "DualFLDIsothermal = %"ISYM, &isothermal);
	ret += sscanf(line, "DualFLDUseXray = %"ISYM, &use_xray);
	ret += sscanf(line, "DualFLDUseUV = %"ISYM, &use_uv);
	ret += sscanf(line, "DualFLDXrayStatic = %"ISYM, &staticXr);
	ret += sscanf(line, "DualFLDUVStatic = %"ISYM, &staticUV);
	ret += sscanf(line, "DualFLDXrayDiffusive = %"ISYM, &xr_diff);
	ret += sscanf(line, "DualFLDUVSpectrum = %"ISYM, &UVSpectrum);
	ret += sscanf(line, "DualFLDUVFrequency = %"FSYM, &UVFrequency);
	ret += sscanf(line, "DualFLDXraySpectrum = %"ISYM, &XrSpectrum);
	ret += sscanf(line, "DualFLDXrayFrequency = %"FSYM, &XrFrequency);
	ret += sscanf(line, "DualFLDMaxDt = %"FSYM, &maxdt);
	ret += sscanf(line, "DualFLDMinDt = %"FSYM, &mindt);
	ret += sscanf(line, "DualFLDInitDt = %"FSYM, &initdt);
	ret += sscanf(line, "DualFLDDtNorm = %"FSYM, &dtnorm);
	ret += sscanf(line, "DualFLDDtGrowth = %"FSYM, &dtgrowth);
	ret += sscanf(line, "DualFLDDtControl = %"ISYM, &dt_control);
	ret += sscanf(line, "DualFLDUVTimeAccuracy = %"FSYM, &timeAccuracyUV);
	ret += sscanf(line, "DualFLDXrayTimeAccuracy = %"FSYM, &timeAccuracyXr);
	ret += sscanf(line, "DualFLDUVScaling = %"FSYM, &UVScale);
	ret += sscanf(line, "DualFLDXrayScaling = %"FSYM, &XrScale);
	ret += sscanf(line, "DualFLDAutoScaling = %"ISYM, &autoscale);
	ret += sscanf(line, "DualFLDTheta = %"FSYM, &theta);
	ret += sscanf(line, "DualFLDUVBoundaryX0Faces = %"ISYM" %"ISYM, 
		      UVBdryType[0], UVBdryType[0]+1);
	if (rank > 1) {
	  ret += sscanf(line, "DualFLDUVBoundaryX1Faces = %"ISYM" %"ISYM,
			UVBdryType[1], UVBdryType[1]+1);
	  if (rank > 2) {
	    ret += sscanf(line, "DualFLDUVBoundaryX2Faces = %"ISYM" %"ISYM,
			  UVBdryType[2], UVBdryType[2]+1);
	  }
	}
	ret += sscanf(line, "DualFLDXrayBoundaryX0Faces = %"ISYM" %"ISYM, 
		      XrBdryType[0], XrBdryType[0]+1);
	if (rank > 1) {
	  ret += sscanf(line, "DualFLDXrayBoundaryX1Faces = %"ISYM" %"ISYM,
			XrBdryType[1], XrBdryType[1]+1);
	  if (rank > 2) {
	    ret += sscanf(line, "DualFLDXrayBoundaryX2Faces = %"ISYM" %"ISYM,
			  XrBdryType[2], XrBdryType[2]+1);
	  }
	}
	ret += sscanf(line, "DualFLDSolToleranceXray = %"FSYM, &sol_tolerance_Xr);
	ret += sscanf(line, "DualFLDSolToleranceUV = %"FSYM, &sol_tolerance_UV);
	ret += sscanf(line, "DualFLDMaxMGItersXray = %i", &sol_MGmaxit_Xr);
	ret += sscanf(line, "DualFLDMaxMGItersUV = %i", &sol_MGmaxit_UV);
	ret += sscanf(line, "DualFLDMaxKryItersXray = %i", &sol_KryMaxit_Xr);
	ret += sscanf(line, "DualFLDMaxKryItersUV = %i", &sol_KryMaxit_UV);
	ret += sscanf(line, "DualFLDMGRelaxTypeXray = %i", &sol_rlxtype_Xr);
	ret += sscanf(line, "DualFLDMGRelaxTypeUV = %i", &sol_rlxtype_UV);
	ret += sscanf(line, "DualFLDMGPreRelaxXray = %i", &sol_npre_Xr);
	ret += sscanf(line, "DualFLDMGPreRelaxUV = %i", &sol_npre_UV);
	ret += sscanf(line, "DualFLDMGPostRelaxXray = %i", &sol_npost_Xr);
	ret += sscanf(line, "DualFLDMGPostRelaxUV = %i", &sol_npost_UV);
	ret += sscanf(line, "DualFLDKrylovMethodXray = %"ISYM, &sol_Krylov_Xr);
	ret += sscanf(line, "DualFLDKrylovMethodUV = %"ISYM, &sol_Krylov_UV);
	ret += sscanf(line, "DualFLDLimiterRmin = %"FSYM, &FLD_Rmin);
	ret += sscanf(line, "DualFLDLimiterDmax = %"FSYM, &FLD_Dmax);
	ret += sscanf(line, "DualFLDLimiterKmin = %"FSYM, &FLD_Kmin);
	ret += sscanf(line, "DualFLDLimiterEmin = %"FSYM, &FLD_Emin);

	ret += sscanf(line, "DualFLDWeakScaling = %"ISYM, &WeakScaling);
	ret += sscanf(line, "DualFLDNumSourcesXray = %"ISYM, &NumSourcesXr);
	ret += sscanf(line, "DualFLDNumSourcesUV = %"ISYM, &NumSourcesUV);
	if (sscanf(line, "DualFLDSourceLocationXray[%"ISYM"] = %"FSYM" %"FSYM" %"FSYM, 
	 	   &ival, &fval1, &fval2, &fval3) == 4) {
	  ret++;  
	  if (ival >= MAX_FLD_SOURCES) {
	    ENZO_VFAIL("DualFLDSourceLocationXray source %"ISYM" > maximum allowed.\n", ival)
          }
	  SourceLocationXr[ival][0] = fval1;
	  SourceLocationXr[ival][1] = fval2;
	  SourceLocationXr[ival][2] = fval3;
	}
	if (sscanf(line, "DualFLDSourceLocationUV[%"ISYM"] = %"FSYM" %"FSYM" %"FSYM, 
	 	   &ival, &fval1, &fval2, &fval3) == 4) {
	  ret++;  
	  if (ival >= MAX_FLD_SOURCES) {
	    ENZO_VFAIL("DualFLDSourceLocationUV source %"ISYM" > maximum allowed.\n", ival)
          }
	  SourceLocationUV[ival][0] = fval1;
	  SourceLocationUV[ival][1] = fval2;
	  SourceLocationUV[ival][2] = fval3;
	}
	if (sscanf(line, "DualFLDSourceEnergyXray[%"ISYM"] = %"FSYM, &ival, &fval1) == 2) {
	  ret++;  
	  if (ival >= MAX_FLD_SOURCES) {
	    ENZO_VFAIL("DualFLDSourceEnergyXray source %"ISYM" > maximum allowed.\n", ival)
          }
	  SourceEnergyXr[ival] = fval1;
	}
	if (sscanf(line, "DualFLDSourceEnergyUV[%"ISYM"] = %"FSYM, &ival, &fval1) == 2) {
	  ret++;  
	  if (ival >= MAX_FLD_SOURCES) {
	    ENZO_VFAIL("DualFLDSourceEnergyUV source %"ISYM" > maximum allowed.\n", ival)
          }
	  SourceEnergyUV[ival] = fval1;
	}

      }  // end loop over file lines
    }  // end successful file open
  }  // end if file name exists
 
  // clean up
  delete[] dummy;
  rewind(fptr);
  fclose(fptr);

  // set Xray/UV flags based on input values
  UseXray = (use_xray == 1);
  UseUV = (use_uv == 1);

  // set XrayDiffusive flag based on input value; if set, disable temporal accuracy limit
  XrayDiffusive = (xr_diff == 1);
  if (XrayDiffusive)  timeAccuracyXr = huge_number;


  // set static radiation fields based on input values staticXr and staticUV
  XrStatic = (staticXr == 1);
  UVStatic = (staticUV == 1);
  
  // update static flags if a radiation field is unused
  if (!UseXray)  XrStatic = 1;
  if (!UseUV)    UVStatic = 1;

  //// Check input parameters ////

  // ensure that at least one radiation field is enabled
  if (!(UseXray || UseUV))
    ENZO_FAIL("DualFLD::ReadParameters Error: module enabled but all radiation physics disabled");

  // check that these give appropriate values, otherwise set dim to periodic
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) {
      if (UseXray)
	if ((XrBdryType[dim][face] < 0) || (XrBdryType[dim][face] >= NUM_FLD_BDRY_TYPES)) {
	  fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDXrayBoundaryX*Face = %"ISYM
		  "(dim %"ISYM", face %"ISYM"), re-setting Xray BC to periodic\n", dim, face,
		  XrBdryType[dim][face]);
	  XrBdryType[dim][0] = XrBdryType[dim][1] = 0;
	}
      if (UseUV)
	if ((UVBdryType[dim][face] < 0) || (UVBdryType[dim][face] >= NUM_FLD_BDRY_TYPES)) {
	  fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDUVBoundaryX*Face = %"ISYM
		  "(dim %"ISYM", face %"ISYM"), re-setting UV BC to periodic\n", dim, face,
		  UVBdryType[dim][face]);
	  UVBdryType[dim][0] = UVBdryType[dim][1] = 0;
	}
    }
  
  // check that periodic faces match
  for (dim=0; dim<rank; dim++) {
    if (UseXray)
      if ((XrBdryType[dim][0]*XrBdryType[dim][1] == 0) && 
	  (XrBdryType[dim][0]+XrBdryType[dim][1] != 0)) {
	fprintf(stderr,"DualFLD::ReadParameters Warning: non-matching periodic Xray BCs, dim %"ISYM"\n",dim);
	XrBdryType[dim][0] = 0;
	XrBdryType[dim][1] = 0;
      }
    if (UseUV)
      if ((UVBdryType[dim][0]*UVBdryType[dim][1] == 0) && 
	  (UVBdryType[dim][0]+UVBdryType[dim][1] != 0)) {
	fprintf(stderr,"DualFLD::ReadParameters Warning: non-matching periodic UV BCs, dim %"ISYM"\n",dim);
	UVBdryType[dim][0] = 0;
	UVBdryType[dim][1] = 0;
      }
  }

  // check that Xray and UV periodicity match
  if (UseUV && UseXray)
    for (dim=0; dim<rank; dim++) {
      if (((XrBdryType[dim][0]*XrBdryType[dim][1] == 0) && 
	   (UVBdryType[dim][0]*UVBdryType[dim][1] != 0)) ||
	  ((UVBdryType[dim][0]*UVBdryType[dim][1] == 0) && 
	   (XrBdryType[dim][0]*XrBdryType[dim][1] != 0))) {
	fprintf(stderr,"DualFLD::ReadParameters Warning: non-matching periodicity of Xray and UV fields, dim %"ISYM", Xr = (%"ISYM",%"ISYM"), UV = (%"ISYM",%"ISYM")\n", dim, XrBdryType[dim][0], XrBdryType[dim][1], UVBdryType[dim][0], UVBdryType[dim][1]);
	XrBdryType[dim][0] = 0;
	XrBdryType[dim][1] = 0;
	UVBdryType[dim][0] = 0;
	UVBdryType[dim][1] = 0;
      }
    }

  // if only one field is used, set other field to match BC type
  // (for simplified checking later on) 
  if (UseUV && !UseXray) 
    for (dim=0; dim<rank; dim++)
      for (face=0; face<2; face++)
	XrBdryType[dim][face] = UVBdryType[dim][face];
  if (UseXray && !UseUV) 
    for (dim=0; dim<rank; dim++)
      for (face=0; face<2; face++)
	UVBdryType[dim][face] = XrBdryType[dim][face];

  if (debug) {
    if (UseXray)
      printf("DualFLD::ReadParameters: XrBdryType = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
	     XrBdryType[0][0], XrBdryType[0][1], XrBdryType[1][0], 
	     XrBdryType[1][1], XrBdryType[2][0], XrBdryType[2][1]);
    if (UseUV)
      printf("DualFLD::ReadParameters: UVBdryType = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
	     UVBdryType[0][0], UVBdryType[0][1], UVBdryType[1][0], 
	     UVBdryType[1][1], UVBdryType[2][0], UVBdryType[2][1]);
  }

  // check for valid radiation spectra (these are numbered from -1 upwards)
  if (UseUV && ((UVSpectrum < -1) || (UVSpectrum > NUM_FLD_SED_TYPES-2))) {
    ENZO_VFAIL("DualFLD::ReadParameters: UV enabled, but DualFLDUVSpectrum = %"ISYM" is illegal\n",
	       UVSpectrum);
  }
  if (UseXray && ((XrSpectrum < -1) || (XrSpectrum > NUM_FLD_SED_TYPES-2))) {
    ENZO_VFAIL("DualFLD::ReadParameters: Xray enabled, but DualFLDXraySpectrum = %"ISYM" is illegal\n",
	       XrSpectrum);
  }

  // monochromatic radiation frequencies must be >0
  if ((UVSpectrum < 0) && (UVFrequency <= 0.0)) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal UVFrequency = %g\n",UVFrequency);
    fprintf(stderr,"   re-setting UVFrequency 13.6\n");
    UVFrequency = 13.6;   // default is the hydrogen ionization threshold
  }
  if ((XrSpectrum < 0) && (XrFrequency <= 0.0)) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal XrFrequency = %g\n",XrFrequency);
    fprintf(stderr,"   re-setting XrFrequency 500.0\n");
    XrFrequency = 500.0;  // default is 0.5 keV
  }
  
  // check for a legal number of sources
  if (NumSourcesXr > MAX_FLD_SOURCES) {
    ENZO_VFAIL("DualFLD::ReadParameters: too many Xray radiation sources requested (%"
	       ISYM" > %"ISYM"). Halting run\n", NumSourcesXr, MAX_FLD_SOURCES);
  }
  if (NumSourcesUV > MAX_FLD_SOURCES) {
    ENZO_VFAIL("DualFLD::ReadParameters: too many UV radiation sources requested (%"
	       ISYM" > %"ISYM"). Halting run\n", NumSourcesUV, MAX_FLD_SOURCES);
  }
  if ((NumSourcesXr > 0) && (!UseXray)) {
    fprintf(stderr,"DualFLD::ReadParameters: NumSourcesXr = %"ISYM
	    ", but X-ray radiation is disabled!  re-setting NumSourcesXr\n", NumSourcesXr);
    NumSourcesXr = 0;
  }
  if ((NumSourcesUV > 0) && (!UseUV)) {
    fprintf(stderr,"DualFLD::ReadParameters: NumSourcesUV = %"ISYM
	    ", but UV radiation is disabled!  re-setting NumSourcesUV\n", NumSourcesUV);
    NumSourcesUV = 0;
  }
  if (debug) {
    printf("DualFLD::ReadParameters NumSourcesXr = %"ISYM"\n", NumSourcesXr);
    printf("DualFLD::ReadParameters NumSourcesUV = %"ISYM"\n", NumSourcesUV);
  }

  // check for legal SourceLocationXr/SourceLocationUV
  for (isrc=0; isrc<NumSourcesXr; isrc++) {
    for (dim=0; dim<rank; dim++) 
      if ((SourceLocationXr[isrc][dim] < DomainLeftEdge[dim]) || 
	  (SourceLocationXr[isrc][dim] > DomainRightEdge[dim])) {
	ENZO_VFAIL("DualFLD::ReadParameters: Xray source %"ISYM
		   " is outside the computational domain. Halting run\n", isrc);
      }
    if (debug) {
      printf("DualFLD::ReadParameters Xray source %"ISYM": location =", isrc);
      for (dim=0; dim<rank; dim++)  printf(" %g", SourceLocationXr[isrc][dim]);
      printf("\n");
    }
  }
  for (isrc=0; isrc<NumSourcesUV; isrc++) {
    for (dim=0; dim<rank; dim++) 
      if ((SourceLocationUV[isrc][dim] < DomainLeftEdge[dim]) || 
	  (SourceLocationUV[isrc][dim] > DomainRightEdge[dim])) {
	ENZO_VFAIL("DualFLD::ReadParameters: UV source %"ISYM
		   " is outside the computational domain. Halting run\n", isrc);
      }
    if (debug) {
      printf("DualFLD::ReadParameters UV source %"ISYM": location =", isrc);
      for (dim=0; dim<rank; dim++)  printf(" %g", SourceLocationUV[isrc][dim]);
      printf("\n");
    }
  }

  // check for legal SourceEnergy
  for (isrc=0; isrc<NumSourcesXr; isrc++)
    if (SourceEnergyXr[isrc] < 0.0) {
      ENZO_VFAIL("DualFLD::ReadParameters: Xray source %"ISYM
		 " has illegal energy = %g. Halting run\n", isrc, SourceEnergyXr[isrc])
    }
  for (isrc=0; isrc<NumSourcesUV; isrc++)
    if (SourceEnergyUV[isrc] < 0.0) {
      ENZO_VFAIL("DualFLD::ReadParameters: UV source %"ISYM
		 " has illegal energy = %g. Halting run\n", isrc, SourceEnergyUV[isrc])
    }
  if (debug) {
    for (isrc=0; isrc<NumSourcesXr; isrc++)
      printf("DualFLD::ReadParameters Xray source %"ISYM": energy = %g photons/sec\n", 
	     isrc, SourceEnergyXr[isrc]);
    for (isrc=0; isrc<NumSourcesUV; isrc++)
      printf("DualFLD::ReadParameters UV source %"ISYM": energy = %g photons/sec\n", 
	     isrc, SourceEnergyUV[isrc]);
  }

  // if doing a weak-scaling run, adjust source locations appropriately
  if (WeakScaling) {
    for (isrc=0; isrc<NumSourcesXr; isrc++) 
      for (dim=0; dim<rank; dim++) {
	float frac = (SourceLocationXr[isrc][dim] - DomainLeftEdge[dim]) 
                   / (DomainRightEdge[dim] - DomainLeftEdge[dim]);
	OriginalSourceLocationXr[isrc][dim] = SourceLocationXr[isrc][dim];
	SourceLocationXr[isrc][dim] = ThisGrid->GridData->GetGridLeftEdge(dim) +
                                      frac*(ThisGrid->GridData->GetGridRightEdge(dim) -
					    ThisGrid->GridData->GetGridLeftEdge(dim));
      }
    for (isrc=0; isrc<NumSourcesUV; isrc++) 
      for (dim=0; dim<rank; dim++) {
	float frac = (SourceLocationUV[isrc][dim] - DomainLeftEdge[dim]) 
                   / (DomainRightEdge[dim] - DomainLeftEdge[dim]);
	OriginalSourceLocationUV[isrc][dim] = SourceLocationUV[isrc][dim];
	SourceLocationUV[isrc][dim] = ThisGrid->GridData->GetGridLeftEdge(dim) +
                                      frac*(ThisGrid->GridData->GetGridRightEdge(dim) -
					    ThisGrid->GridData->GetGridLeftEdge(dim));
      }
  }
  

  // maxdt gives the maximum radiation time step size
  if (maxdt <= 0.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDMaxDt = %g\n",maxdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    maxdt = huge_number;  // default is no limit
  }

  // mindt gives the minimum radiation time step size
  if (mindt < 0.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDMinDt = %g\n",mindt);
    fprintf(stderr,"   re-setting to %g\n",0.0);
    mindt = 0.0;  // default is 0.0
  }

  // initdt gives the initial time step size
  if (initdt <= 0.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDInitDt = %g\n",initdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    initdt = huge_number;  // default is no limit
  }

  // dt_control gives the time step controller algorithm (starts numbering at -1)
  if ((dt_control < -1) || (dt_control > NUM_FLD_DT_CONTROLLERS-2)) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDDtControl = %"ISYM"\n", dt_control);
    fprintf(stderr,"   re-setting to -1 (original controller)\n");
    dt_control = -1;
  }
  
  // *Scale give variable scalings for implicit solver
  if (UVScale <= 0.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDUVScaling = %g\n",UVScale);
    fprintf(stderr,"   re-setting to 1.0\n");
    UVScale = 1.0;  // default is no scaling
  }
  if (XrScale <= 0.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDXrayScaling = %g\n",XrScale);
    fprintf(stderr,"   re-setting to 1.0\n");
    XrScale = 1.0;  // default is no scaling
  }
  autoScale = (autoscale != 0);  // set bool based on integer input
  if (debug)
    printf("DualFLD::ReadParameters: UVScale = %g, XrScale = %g, autoScale = %i\n",
	   UVScale, XrScale, autoscale);


  // timeAccuracy* gives the desired percent change in values per step
  if (timeAccuracyUV <= 0.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDUVTimeAccuracy = %g\n",timeAccuracyUV);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    timeAccuracyUV = huge_number;  // default is no limit
  }
  if (timeAccuracyXr <= 0.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDXrayTimeAccuracy = %g\n",timeAccuracyXr);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    timeAccuracyXr = huge_number;  // default is no limit
  }
  if (!UseUV && (timeAccuracyUV != huge_number)) {
    fprintf(stderr,"DualFLD::ReadParameters: warning, DualFLDUVTimeAccuracy set but UV field is unused\n");
    timeAccuracyUV = huge_number;  // disable
  }
  if (!UseXray && (timeAccuracyXr != huge_number)) {
    fprintf(stderr,"DualFLD::ReadParameters: warning, DualFLDXrayTimeAccuracy set but Xray field is unused\n");
    timeAccuracyXr = huge_number;  // disable
  }

  // dtnorm gives the norm for calculating allowed relative change per step
  if (dtnorm < 0.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDDtNorm = %g\n",dtnorm);
    fprintf(stderr,"   re-setting to 2.0 (2-norm)\n");
    dtnorm = 2.0;  // default is 2-norm
  }

  // dtgrowth gives the maximum growth factor in dt per step
  if (dtgrowth < 1.0 || dtgrowth > 10.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDDtGrowth = %g\n",dtgrowth);
    fprintf(stderr,"   re-setting to 1.1\n");
    dtgrowth = 1.1;
  }

  // theta gives the implicit time-stepping method (1->BE, 0.5->Trap, 0->FE)
  if ((theta < 0.0) || (theta > 1.0)) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDTheta = %g\n", theta);
    fprintf(stderr,"   re-setting theta to 1.0 (Backwards Euler)\n");
    theta = 1.0;  // default is backwards Euler
  }

  // check limiter parameters
  if (FLD_Rmin <= 0.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDLimiterRmin = %g\n",FLD_Rmin);
    fprintf(stderr,"   re-setting to 0.01\n");
    FLD_Rmin = 0.01;
  }
  if (FLD_Dmax <= 0.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDLimiterDmax = %g\n",FLD_Dmax);
    fprintf(stderr,"   re-setting to 0.01\n");
    FLD_Dmax = 0.01;
  }
  if (FLD_Kmin <= 0.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDLimiterKmin = %g\n",FLD_Kmin);
    fprintf(stderr,"   re-setting to 1e-20\n");
    FLD_Kmin = 1e-20;
  }
  if (FLD_Emin <= 0.0) {
    fprintf(stderr,"DualFLD::ReadParameters: illegal DualFLDLimiterEmin = %g\n",FLD_Emin);
    fprintf(stderr,"   re-setting to 1e-30\n");
    FLD_Emin = 1e-30;
  }

  // ensure that Enzo was called with RadiativeCooling enabled 
  // (since AMRFLDSplit doesn't handle chemistry/cooling)
  if (!RadiativeCooling) 
    ENZO_FAIL("DualFLD::ReadParameters: RadiativeCooling must be on!  Halting run");

  //   check linear solver parameters
  if ((sol_Krylov_Xr < 0) || (sol_Krylov_Xr >= NUM_FLD_SOL_TYPES)) {
    fprintf(stderr,"Illegal DualFLDKrylovMethodXray = %"ISYM",  Setting to 1 (BiCGStab)\n", 
	    sol_Krylov_Xr);
    sol_Krylov_Xr = 1;
  }
  if ((sol_Krylov_UV < 0) || (sol_Krylov_UV >= NUM_FLD_SOL_TYPES)) {
    fprintf(stderr,"Illegal DualFLDKrylovMethodUV = %"ISYM",  Setting to 1 (BiCGStab)\n", 
	    sol_Krylov_UV);
    sol_Krylov_UV = 1;
  }

  if ((sol_tolerance_Xr < 1.0e-10) || (sol_tolerance_Xr > 1.0)) {
    fprintf(stderr,"Illegal DualFLDSolToleranceXray = %g. Setting to 1e-2\n",
	    sol_tolerance_Xr);
    sol_tolerance_Xr = 1e-2;
  }
  if ((sol_tolerance_UV < 1.0e-10) || (sol_tolerance_UV > 1.0)) {
    fprintf(stderr,"Illegal DualFLDSolToleranceUV = %g. Setting to 1e-2\n",
	    sol_tolerance_UV);
    sol_tolerance_UV = 1e-2;
  }

  if (sol_MGmaxit_Xr <= 0) {
    fprintf(stderr,"Illegal DualFLDMaxMGItersXray = %"ISYM". Setting to 10\n",
	    sol_MGmaxit_Xr);
    sol_MGmaxit_Xr = 10;
  }
  if (sol_MGmaxit_UV <= 0) {
    fprintf(stderr,"Illegal DualFLDMaxMGItersUV = %"ISYM". Setting to 10\n",
	    sol_MGmaxit_UV);
    sol_MGmaxit_UV = 10;
  }
  if (sol_KryMaxit_Xr <= 0) {
    fprintf(stderr,"Illegal DualFLDMaxKryItersXray = %"ISYM". Setting to 10\n",
	    sol_KryMaxit_Xr);
    sol_KryMaxit_Xr = 10;
  }
  if (sol_KryMaxit_UV <= 0) {
    fprintf(stderr,"Illegal DualFLDMaxKryItersUV = %"ISYM". Setting to 10\n",
	    sol_KryMaxit_UV);
    sol_KryMaxit_UV = 10;
  }

  if ((sol_rlxtype_Xr < 0) || (sol_rlxtype_Xr >= NUM_HYPRE_RLX_TYPES)) {
    fprintf(stderr,"Illegal DualFLDMGRelaxTypeXray = %"ISYM". Setting to 1 (weighted Jacobi)\n",
	    sol_rlxtype_Xr);
    sol_rlxtype_Xr = 1;
  }
  if ((sol_rlxtype_UV < 0) || (sol_rlxtype_UV >= NUM_HYPRE_RLX_TYPES)) {
    fprintf(stderr,"Illegal DualFLDMGRelaxTypeUV = %"ISYM". Setting to 1 (weighted Jacobi)\n",
	    sol_rlxtype_UV);
    sol_rlxtype_UV = 1;
  }

  if (sol_npre_Xr < 1) {
    fprintf(stderr,"Illegal DualFLDMGPreRelaxXray = %"ISYM", Setting to 1\n",
	    sol_npre_Xr);
    sol_npre_Xr = 1;
  }
  if (sol_npre_UV < 1) {
    fprintf(stderr,"Illegal DualFLDMGPreRelaxUV = %"ISYM", Setting to 1\n",
	    sol_npre_UV);
    sol_npre_UV = 1;
  }

  if (sol_npost_Xr < 1) {
    fprintf(stderr,"Illegal DualFLDMGPostRelaxXray = %"ISYM". Setting to 1\n",
	    sol_npost_Xr);
    sol_npost_Xr = 1;
  }
  if (sol_npost_UV < 1) {
    fprintf(stderr,"Illegal DualFLDMGPostRelaxUV = %"ISYM". Setting to 1\n",
	    sol_npost_UV);
    sol_npost_UV = 1;
  }

  
  return SUCCESS;

}
#endif   /* TRANSFER */
