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
/  SetupSources routine
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Handles user inputs to set up relevant desired emissivity
/           sources for each radiation frequency/group.
/
************************************************************************/
#ifdef TRANSFER
#include "DualFLD.h"
#include "SED.h"
#include "BlackbodySED.h"
#include "MonochromaticSED.h"


// function prototypes
float SED_integral(SED &sed, float a, float b, bool convertHz);


// SED "moment" class for computing radiation field contributions
class SED_SpectrumMoment : public virtual SED {
 private:
  SED *baseSED;     // base SED
 public:
  SED_SpectrumMoment(SED &base) { this->baseSED = &base; };
  bool monochromatic() { return this->baseSED->monochromatic(); };
  float lower_bound() { return this->baseSED->lower_bound(); };
  float upper_bound() { return this->baseSED->upper_bound(); };
  float value(float hnu) { return this->baseSED->value(hnu)*hnu0_HI/hnu; };
};


int DualFLD::SetupSources(HierarchyEntry *RootGrid,
			  int SourceType[MAX_FLD_SOURCES], 
			  float SourceEnergy[MAX_FLD_SOURCES]) {

  if (debug)  printf("Entering DualFLD::SetupSources routine\n");

  // if this is a weak scaling test, overwrite SourceLocation for each 
  // source to replicate source setup on each root grid tile
  int isrc, ibin, dim;
  if (WeakScaling) {
    if (UseXray)
      for (isrc=0; isrc<NumSourcesXr; isrc++) 
	for (dim=0; dim<rank; dim++) {
	  float frac = (SourceLocationXr[isrc][dim] - DomainLeftEdge[dim]) 
                     / (DomainRightEdge[dim] - DomainLeftEdge[dim]);
	  OriginalSourceLocationXr[isrc][dim] = SourceLocationXr[isrc][dim];
	  SourceLocationXr[isrc][dim] = RootGrid->GridData->GetGridLeftEdge(dim) +
                                        frac*(RootGrid->GridData->GetGridRightEdge(dim) -
			                      RootGrid->GridData->GetGridLeftEdge(dim));
	}
    if (UseUV)
      for (isrc=0; isrc<NumSourcesUV; isrc++) 
	for (dim=0; dim<rank; dim++) {
	  float frac = (SourceLocationUV[isrc][dim] - DomainLeftEdge[dim]) 
                     / (DomainRightEdge[dim] - DomainLeftEdge[dim]);
	  OriginalSourceLocationUV[isrc][dim] = SourceLocationUV[isrc][dim];
	  SourceLocationUV[isrc][dim] = RootGrid->GridData->GetGridLeftEdge(dim) +
                                        frac*(RootGrid->GridData->GetGridRightEdge(dim) -
			                      RootGrid->GridData->GetGridLeftEdge(dim));
	}
  }

  // fill SourceGroupEnergy for sources that are specified by type/energy
  for (isrc=0; isrc<NumSources; isrc++) {

    // skip sources with SourceType of -1 (indicates alternate source strategy)
    if (SourceType[isrc] == -1)  continue;

    // set up source SED
    SED *tmp_src;
    switch (SourceType[isrc]) {
    case 0:      // monochromatic source at hnu = 13.6 eV
      tmp_src = new MonochromaticSED(13.6);
      break;
    case 1:      // blackbody spectrum at T = 1e5 K
      tmp_src = new BlackbodySED(1.0e5);
      break;
    default:
      ENZO_VFAIL("DualFLD::SetupSources error, source type %"ISYM" does not exist\n", 
		 SourceType[isrc]);
    }

    // set up source SED spectrum moment
    SED_SpectrumMoment tmp_src_m(*tmp_src);

    // compute totalintegral of SED over all ionizing frequencies
    float total_integral = SED_integral(*tmp_src, hnu0_HI, -1.0, true);

    // first handle Xray contribution

    // second handle UV contribution
    for (ibin=0; ibin<NumRadiationFields; ibin++) {

      // skip monochromatic fields
      if (FieldMonochromatic[ibin]) {
	SourceGroupEnergy[isrc][ibin] = 0.0;
	continue;
      }

      // computing source contributions to this radiation field
      SourceGroupEnergy[isrc][ibin] = SourceEnergy[isrc] * hplanck / total_integral
	* SED_integral(tmp_src_m, FrequencyBand[ibin][0], FrequencyBand[ibin][1], true);
    }

    // report results to stdout
    if (debug) {
      printf("DualFLD::Initialize source %"ISYM" has group energies", isrc);
      for (ibin=0; ibin<NumRadiationFields; ibin++)  
	printf(" %g", SourceGroupEnergy[isrc][ibin]);
      printf("\n");
    }

    // clean up
    delete tmp_src;

  }   // end isrc loop

  return SUCCESS;
}
#endif // TRANSFER
