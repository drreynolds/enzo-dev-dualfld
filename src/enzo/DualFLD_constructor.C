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
/  Constructor routine
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Initializes all values to illegal numbers, and sets all 
/           arrays to NULL;  Requires call to Initialize to actually 
/           set up these values.
/
************************************************************************/
#ifdef TRANSFER
#include "DualFLD.h"


DualFLD::DualFLD() {

//   if (debug)  printf("\nEntering DualFLD::constructor routine\n");
  int dim, face, i, src;
#ifndef MPI_INT
  int MPI_PROC_NULL = -3;
#endif

  // initialize total RT time to zero
  RTtime = 0.0;
  HYPREtime = 0.0;

  // initialize HYPRE values to -1/NULL
  mattype          = -1;
  stSize           = -1;
#ifdef USE_HYPRE
  grid             = NULL;
  stencil          = NULL;
#endif
  sol_tolerance_Xr = -1.0;
  sol_tolerance_UV = -1.0;
  sol_MGmaxit_Xr   = -1;
  sol_MGmaxit_UV   = -1;
  sol_KryMaxit_Xr  = -1;
  sol_KryMaxit_UV  = -1;
  sol_rlxtype_Xr   = -1;
  sol_rlxtype_UV   = -1;
  sol_npre_Xr      = -1;
  sol_npre_UV      = -1;
  sol_npost_Xr     = -1;
  sol_npost_UV     = -1;
  sol_printl       = -1;
  sol_log          = -1;
  sol_Krylov_Xr    = 1;
  sol_Krylov_UV    = 1;
  totIters_Xr      = 0;
  totIters_UV      = 0;
  for (dim=0; dim<3; dim++) {
    for (face=0; face<2; face++)
      SolvIndices[dim][face] = 0;
    SolvOff[dim] = 0;
  }

  // initialize HYPRE structures to NULL
#ifdef USE_HYPRE
  P      = NULL;
  rhsvec = NULL;
  solvec = NULL;
#endif

  // initialize HYPRE interface array to NULL
  HYPREbuff = NULL;

  // initialize problem grid information to -1/NULL
  rank = -1;
  for (dim=0; dim<3; dim++) {
    layout[dim]   = 1;  // initialize for single-processor run
    location[dim] = 0;  // initialize for single-processor run
    LocDims[dim]  = 1;
    ArrDims[dim]  = 1;
    GlobDims[dim] = 1;
    dx[dim]       = 1.0;
    for (face=0; face<2; face++) {
      OnBdry[dim][face]     = false;
      NBors[dim][face]      = MPI_PROC_NULL;
      GhDims[dim][face]     = 0;
      XrBdryType[dim][face] = -1;
      UVBdryType[dim][face] = -1;
      EdgeVals[dim][face]   = -1.0;
      XrBdryVals[dim][face] = NULL;
      UVBdryVals[dim][face] = NULL;
    }
  }

  // limiter parameters
  FLD_Rmin = -1.0;
  FLD_Dmax = -1.0;
  FLD_Kmin = -1.0;
  FLD_Emin = -1.0;
  
  // initialize time-stepping related data to -1/NULL
  initdt         = 1.0e20;
  maxdt          = 1.0e20;
  mindt          = 0.0;
  timeAccuracyXr = 1.0e20;
  timeAccuracyUV = 1.0e20;
  dt_control     = 2;
  Err_cur_Xr     = 1.0;
  Err_cur_UV     = 1.0;
  Err_new_Xr     = 1.0;
  Err_new_UV     = 1.0;
  dtnorm         = 0.0;
  dtgrowth       = 1.1;
  tnew           = -1.0;
  told           = -1.0;
  dt             = -1.0;
  dtrad          = -1.0;
  theta          = -1.0;

  // initialize problem defining data 
  isothermal   = 0;
  UseXray      = false;
  UseUV        = false;
  XrStatic     = false;
  UVStatic     = false;
  // NGammaDotXr  = 0.0;
  // NGammaDotUV  = 0.0;
  // EtaRadius    = 0.0;
  // EtaCenter[0] = 0.0;
  // EtaCenter[1] = 0.0;
  // EtaCenter[2] = 0.0;

  // initialize ionization source parameters
  NumSourcesXr = 0;
  NumSourcesUV = 0;
  for (src=0; src<MAX_FLD_SOURCES; src++) {
    for (dim=0; dim<3; dim++) {
      SourceLocationXr[src][dim] = 0.0;
      SourceLocationUV[src][dim] = 0.0;
      OriginalSourceLocationXr[src][dim] = 0.0;
      OriginalSourceLocationUV[src][dim] = 0.0;
    }
    SourceEnergyXr[src] = -1.0;
    SourceEnergyUV[src] = -1.0;
  }
  WeakScaling = 0;         // standard run, do not replicate input sources
  
  // initialize constants
  a           = 1.0;
  a0          = 1.0;
  adot        = 0.0;
  adot0       = 0.0;
  aUnits      = 1.0;
  autoScale   = true;
  StartAutoScale = false;
  UVScale     = 1.0;
  UVUnits     = 1.0;
  UVUnits0    = 1.0;
  XrScale     = 1.0;
  XrUnits     = 1.0;
  XrUnits0    = 1.0;
  NiUnits     = 1.0;
  NiUnits0    = 1.0;
  DenUnits    = 1.0;
  DenUnits0   = 1.0;
  LenUnits    = 1.0;
  LenUnits0   = 1.0;
  TimeUnits   = 1.0;
  VelUnits    = 1.0;

  // frequency-space data
  hnu0_HI   = -1.0;
  hnu0_HeI  = -1.0;
  hnu0_HeII = -1.0;
  UVSpectrum  = -1;
  UVFrequency = 0.0;
  XrSpectrum  = -1;
  XrFrequency = 0.0;
  for (i=0; i<7; i++)  RadIntXr[i] = 0.0;
  for (i=0; i<7; i++)  RadIntUV[i] = 0.0;
  Xr_SED = NULL;
  UV_SED = NULL;

  // private solver storage
  sol    = NULL;
  U      = NULL;
  EtaVec = NULL;

}
#endif
