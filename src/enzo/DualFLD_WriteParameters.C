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
/  Parameter output routine
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Writes all necessary internal parameters for problem 
/           restart.
/
************************************************************************/
#ifdef TRANSFER
#include "DualFLD.h"

int DualFLD::WriteParameters(FILE *fptr)
{

  // output all DualFLD parameters
  fprintf(fptr, "DualFLDIsothermal = %i\n", isothermal);

  if (UseXray) 
    fprintf(fptr, "DualFLDUseXray = 1\n");
  else
    fprintf(fptr, "DualFLDUseXray = 0\n");
  if (UseUV)
    fprintf(fptr, "DualFLDUseUV = 1\n");
  else
    fprintf(fptr, "DualFLDUseUV = 0\n");

  if (XrStatic)
    fprintf(fptr, "DualFLDXrayStatic = 1\n");
  else
    fprintf(fptr, "DualFLDXrayStatic = 0\n");
  if (UVStatic)
    fprintf(fptr, "DualFLDUVStatic = 1\n");
  else
    fprintf(fptr, "DualFLDUVStatic = 0\n");

  if (XrayDiffusive)
    fprintf(fptr, "DualFLDXrayDiffusive = 1\n");
  else
    fprintf(fptr, "DualFLDXrayDiffusive = 0\n");

  fprintf(fptr, "DualFLDUVSpectrum = %"ISYM"\n",    UVSpectrum);
  fprintf(fptr, "DualFLDXraySpectrum = %"ISYM"\n",  XrSpectrum);

  fprintf(fptr, "DualFLDUVFrequency = %22.16e\n",   UVFrequency);
  fprintf(fptr, "DualFLDXrayFrequency = %22.16e\n", XrFrequency);

  fprintf(fptr, "DualFLDMaxDt = %22.16e\n",         maxdt);
  fprintf(fptr, "DualFLDMinDt = %22.16e\n",         mindt);
  fprintf(fptr, "DualFLDInitDt = %22.16e\n",        dtrad);  // set restart initdt to current step
  fprintf(fptr, "DualFLDDtNorm = %22.16e\n",        dtnorm);
  fprintf(fptr, "DualFLDDtGrowth = %22.16e\n",      dtgrowth);
  fprintf(fptr, "DualFLDDtControl = %i\n",          dt_control);
  fprintf(fptr, "DualFLDUVTimeAccuracy = %22.16e\n",   timeAccuracyUV);
  fprintf(fptr, "DualFLDXrayTimeAccuracy = %22.16e\n", timeAccuracyXr);
  fprintf(fptr, "DualFLDUVScaling = %22.16e\n",     UVScale);
  fprintf(fptr, "DualFLDXrayScaling = %22.16e\n",   XrScale);
  if (autoScale) {
    fprintf(fptr, "DualFLDAutoScaling = 1\n");
  } else {
    fprintf(fptr, "DualFLDAutoScaling = 0\n");
  }
  fprintf(fptr, "DualFLDTheta = %22.16e\n",         theta);

  fprintf(fptr, "DualFLDUVBoundaryX0Faces = %"ISYM" %"ISYM"\n", 
	  UVBdryType[0][0], UVBdryType[0][1]);
  if (rank > 1) {
    fprintf(fptr, "DualFLDUVBoundaryX1Faces = %"ISYM" %"ISYM"\n", 
	    UVBdryType[1][0], UVBdryType[1][1]);
    if (rank > 2) {
      fprintf(fptr, "DualFLDUVBoundaryX2Faces = %"ISYM" %"ISYM"\n", 
	      UVBdryType[2][0], UVBdryType[2][1]);
    }
  }
  fprintf(fptr, "DualFLDXrayBoundaryX0Faces = %"ISYM" %"ISYM"\n", 
	  XrBdryType[0][0], XrBdryType[0][1]);
  if (rank > 1) {
    fprintf(fptr, "DualFLDXrayBoundaryX1Faces = %"ISYM" %"ISYM"\n", 
	    XrBdryType[1][0], XrBdryType[1][1]);
    if (rank > 2) {
      fprintf(fptr, "DualFLDXrayBoundaryX2Faces = %"ISYM" %"ISYM"\n", 
	      XrBdryType[2][0], XrBdryType[2][1]);
    }
  }

  fprintf(fptr, "DualFLDSolToleranceXray = %22.16e\n", sol_tolerance_Xr);
  fprintf(fptr, "DualFLDSolToleranceUV = %22.16e\n",   sol_tolerance_UV);

  fprintf(fptr, "DualFLDMaxMGItersXray = %i\n",        sol_MGmaxit_Xr);
  fprintf(fptr, "DualFLDMaxMGItersUV = %i\n",          sol_MGmaxit_UV);

  fprintf(fptr, "DualFLDMaxKryItersXray = %i\n",       sol_KryMaxit_Xr);
  fprintf(fptr, "DualFLDMaxKryItersUV = %i\n",         sol_KryMaxit_UV);

  fprintf(fptr, "DualFLDMGRelaxTypeXray = %i\n",       sol_rlxtype_Xr);    
  fprintf(fptr, "DualFLDMGRelaxTypeUV = %i\n",         sol_rlxtype_UV);

  fprintf(fptr, "DualFLDMGPreRelaxXray = %i\n",        sol_npre_Xr);
  fprintf(fptr, "DualFLDMGPreRelaxUV = %i\n",          sol_npre_UV);

  fprintf(fptr, "DualFLDMGPostRelaxXray = %i\n",       sol_npost_Xr);
  fprintf(fptr, "DualFLDMGPostRelaxUV = %i\n",         sol_npost_UV);

  fprintf(fptr, "DualFLDKrylovMethodXray = %i\n",      sol_Krylov_Xr);
  fprintf(fptr, "DualFLDKrylovMethodUV = %i\n",        sol_Krylov_UV);

  fprintf(fptr, "DualFLDLimiterRmin = %22.16e\n",      FLD_Rmin);
  fprintf(fptr, "DualFLDLimiterDmax = %22.16e\n",      FLD_Dmax);
  fprintf(fptr, "DualFLDLimiterKmin = %22.16e\n",      FLD_Kmin);
  fprintf(fptr, "DualFLDLimiterEmin = %22.16e\n",      FLD_Emin);

  fprintf(fptr, "DualFLDWeakScaling = %i\n",           WeakScaling);
  if (UseXray) {
    fprintf(fptr, "DualFLDNumSourcesXray = %i\n",      NumSourcesXr);
    if (WeakScaling) 
      for (int isrc=0; isrc<NumSourcesXr; isrc++) {
	fprintf(fptr, "DualFLDSourceLocationXray[%"ISYM"] = %22.16e %22.16e %22.16e\n", isrc, 
		OriginalSourceLocationXr[isrc][0], OriginalSourceLocationXr[isrc][1], 
		OriginalSourceLocationXr[isrc][2]);
	fprintf(fptr, "DualFLDSourceEnergyXray[%"ISYM"] = %22.16e\n", isrc, SourceEnergyXr[isrc]);
      }
    else
      for (int isrc=0; isrc<NumSourcesXr; isrc++) {
	fprintf(fptr, "DualFLDSourceLocationXray[%"ISYM"] = %22.16e %22.16e %22.16e\n", isrc, 
		SourceLocationXr[isrc][0], SourceLocationXr[isrc][1], SourceLocationXr[isrc][2]);
	fprintf(fptr, "DualFLDSourceEnergyXray[%"ISYM"] = %22.16e\n", isrc, SourceEnergyXr[isrc]);
      }
  }

  if (UseUV) {
    fprintf(fptr, "DualFLDNumSourcesUV = %i\n",        NumSourcesUV);
    if (WeakScaling) 
      for (int isrc=0; isrc<NumSourcesUV; isrc++) {
	fprintf(fptr, "DualFLDSourceLocationUV[%"ISYM"] = %22.16e %22.16e %22.16e\n", isrc, 
		OriginalSourceLocationUV[isrc][0], OriginalSourceLocationUV[isrc][1], 
		OriginalSourceLocationUV[isrc][2]);
	fprintf(fptr, "DualFLDSourceEnergyUV[%"ISYM"] = %22.16e\n", isrc, SourceEnergyUV[isrc]);
      }
    else
      for (int isrc=0; isrc<NumSourcesUV; isrc++) {
	fprintf(fptr, "DualFLDSourceLocationUV[%"ISYM"] = %22.16e %22.16e %22.16e\n", isrc, 
		SourceLocationUV[isrc][0], SourceLocationUV[isrc][1], SourceLocationUV[isrc][2]);
	fprintf(fptr, "DualFLDSourceEnergyUV[%"ISYM"] = %22.16e\n", isrc, SourceEnergyUV[isrc]);
      }
  }

  // // if doing an ionization problem (ProblemTypes 210-215),  
  // // output additional parameters 
  // if ((ProblemType >= 210) && (ProblemType <= 215)) {
  //   fprintf(fptr, "DualFLDNGammaDotXray = %22.16e\n", NGammaDotXr);
  //   fprintf(fptr, "DualFLDNGammaDotUV = %22.16e\n",   NGammaDotUV);
  //   fprintf(fptr, "DualFLDEtaRadius = %22.16e\n",     EtaRadius);
  //   fprintf(fptr, "DualFLDEtaCenter = %22.16e %22.16e %22.16e\n",
  // 	    EtaCenter[0], EtaCenter[1], EtaCenter[2]);
  // }

  // output relevant units: although these aren't required for restart, 
  // cosmology runs never output the units (why?), making data analysis tough
  fprintf(fptr, "DensityUnits = %22.16e\n", DenUnits);
  fprintf(fptr, "LengthUnits = %22.16e\n",  LenUnits);
  fprintf(fptr, "TimeUnits = %22.16e\n",    TimeUnits);

  return SUCCESS;
}
#endif
