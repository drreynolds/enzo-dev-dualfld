/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Dual Flux-Limited Diffusion Solver, emissivity field calculation routine 
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Computes the emissivity field to enforce on the specified 
/           radiation energy equation.  The source terms are added to 
/           any sources already present in the provided emissivity field, 
/           Eta (assumed to be the size of a BaryonField).  All required 
/           non-DualFLD values are extracted from the provided 
/           HierarchyEntry, ThisGrid.
/
************************************************************************/
#ifdef TRANSFER
#include "DualFLD.h"
#include "phys_constants.h"
#include "BlackbodySED.h"


// function prototypes
float SED_integral(SED &sed, float a, float b, bool convertHz);
 

// SED "moment" class for computing radiation field contributions
class Source_Moment : public virtual SED {
 private:
  SED *baseSED;     // base SED
 public:
  Source_Moment(SED &base) { this->baseSED = &base; };
  bool monochromatic() { return this->baseSED->monochromatic(); };
  float lower_bound() { return this->baseSED->lower_bound(); };
  float upper_bound() { return this->baseSED->upper_bound(); };
  float value(float hnu) { return this->baseSED->value(hnu)*13.6/hnu; };
};


int DualFLD::RadiationSource(HierarchyEntry *ThisGrid, int XrUV, float *Eta) {

  // check that output field data exists
  if (Eta==NULL)
    ENZO_FAIL("DualFLD::Opacity ERROR: no Opacity array!");


  ////////////////////////////
  // compute emissivity for this radiation field 

  // initialize constants and reusable local variables
  int i, j, k, isrc, dim;
  float dV = 1.0;
  float lUn = (LenUnits + LenUnits0)*0.5;
  for (dim=0; dim<rank; dim++) 
    dV *= (dx[dim]*lUn);       // cell volume (proper)
  float h_nu0 = 13.6*ev2erg;   // ionization energy of HI [ergs]
  float h_nu1 = 24.6*ev2erg;   // ionization energy of HeI [ergs]
  float h_nu2 = 54.4*ev2erg;   // ionization energy of HeII [ergs]
  float cellZl, cellZr, cellYl, cellYr, cellXl, cellXr, cellXc, cellYc, cellZc;

  // copy source information to local arrays for simplified computatation
  int NumSources, Spectrum;
  float Frequency;
  float SourceLocation[MAX_FLD_SOURCES][3];
  float SourceEnergy[MAX_FLD_SOURCES];
  SED *RadSED;
  if (XrUV) {    // UV sources
    Spectrum = UVSpectrum;
    Frequency = UVFrequency;
    NumSources = NumSourcesUV;
    for (isrc=0; isrc<NumSources; isrc++) {
      SourceEnergy[isrc] = SourceEnergyUV[isrc];
      for (dim=0; dim<rank; dim++) 
	SourceLocation[isrc][dim] = SourceLocationUV[isrc][dim];
    }
    RadSED = UV_SED;
  } else {       // Xray sources
    Spectrum = XrSpectrum;
    Frequency = XrFrequency;
    NumSources = NumSourcesXr;
    for (isrc=0; isrc<NumSources; isrc++) {
      SourceEnergy[isrc] = SourceEnergyXr[isrc];
      for (dim=0; dim<rank; dim++) 
	SourceLocation[isrc][dim] = SourceLocationXr[isrc][dim];
    }
    RadSED = Xr_SED; 
 }

  // compute integral of SED over all ionizing frequencies
  float total_integral = SED_integral(*RadSED, hnu0_HI, -1.0, true);


  // iterate over sources, adding emissivity to this field at specified location
  for (isrc=0; isrc<NumSources; isrc++) {

    // compute SED "moment" for this source, to properly deposit energy
    Source_Moment tmp_src_m(*RadSED);
    float total_moment = SED_integral(tmp_src_m, 0.01, -1.0, true);

    // all sources have radius of one cell; count number of cells to receive source
    int num_cells = 0;
    for (k=0; k<ArrDims[2]; k++) {
      cellZc = EdgeVals[2][0] + (k-GhDims[2][0]+0.5)*dx[2];        // z-center (comoving) for this cell
      for (j=0; j<ArrDims[1]; j++) {
	cellYc = EdgeVals[1][0] + (j-GhDims[1][0]+0.5)*dx[1];      // y-center (comoving) for this cell
	for (i=0; i<ArrDims[0]; i++) {
	  cellXc = EdgeVals[0][0] + (i-GhDims[0][0]+0.5)*dx[0];    // x-center (comoving) for this cell
	  if ( (fabs(cellXc-SourceLocation[isrc][0]) < dx[0]) &&
	       (fabs(cellYc-SourceLocation[isrc][1]) < dx[1]) &&
	       (fabs(cellZc-SourceLocation[isrc][2]) < dx[2]) )
	    num_cells++;                        // cell is within source region
	} // x-loop
      } // y-loop
    } // z-loop

    // equi-partition energy among affected cells
    float cell_energy = SourceEnergy[isrc] * total_integral * hnu0_HI * ev2erg / total_moment / num_cells / dV;
    for (k=GhDims[2][0]; k<LocDims[2]+GhDims[2][0]; k++) {
      cellZc = EdgeVals[2][0] + (k-GhDims[2][0]+0.5)*dx[2];        // z-center (comoving) for this cell
      for (j=GhDims[1][0]; j<LocDims[1]+GhDims[1][0]; j++) {
	cellYc = EdgeVals[1][0] + (j-GhDims[1][0]+0.5)*dx[1];      // y-center (comoving) for this cell
	for (i=GhDims[0][0]; i<LocDims[0]+GhDims[0][0]; i++) {
	  cellXc = EdgeVals[0][0] + (i-GhDims[0][0]+0.5)*dx[0];    // x-center (comoving) for this cell
	  if ( (fabs(cellXc-SourceLocation[isrc][0]) < dx[0]) &&
	       (fabs(cellYc-SourceLocation[isrc][1]) < dx[1]) &&
	       (fabs(cellZc-SourceLocation[isrc][2]) < dx[2]) )
	    Eta[(k*ArrDims[1] + j)*ArrDims[0] + i] += cell_energy;
	} // x-loop
      } // y-loop
    } // z-loop

  } // end loop over sources from input parameters


	
  /////////////
  // add emissivity based on ProblemType, if desired
	
  switch (ProblemType) {

  case 432:    // emissivity flux along x=0 wall (NGammaDot photons/s/cm^2)
	  
    // only relevant for UV radiation
    if (XrUV) {

      // place ionization sources along left wall (if on this subdomain)
      if (EdgeVals[0][0] == 0.0) {

	// set up source SED and the "moment" SED
	BlackbodySED BB_src(1.0e5);
	Source_Moment tmp_src_m(BB_src);

	// compute integral of SED over all ionizing frequencies
	float tot_int = SED_integral(BB_src, hnu0_HI, -1.0, true);

	// compute integral of moment SED over this band
	float tot_mom = SED_integral(tmp_src_m, 0.01, -1.0, true);

	// determine this group's portion of total blackbody emissivity
	float wall_energy = 1e6 * hnu0_HI * ev2erg * tot_int / tot_mom / dx[0] / lUn;

	// place along wall (i=GhDims[0][0])
	for (k=GhDims[2][0]; k<LocDims[2]+GhDims[2][0]; k++)
	  for (j=GhDims[1][0]; j<LocDims[1]+GhDims[1][0]; j++)
	    Eta[(k*ArrDims[1] + j)*ArrDims[0] + GhDims[0][0]] = wall_energy;
      }
    }
    break;
	  
  }  // switch (ProblemType)
	
  return SUCCESS;
}

#endif
