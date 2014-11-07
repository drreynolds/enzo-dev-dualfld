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
class BB_Moment : public virtual SED {
 private:
  SED *baseSED;     // base SED
 public:
  BB_Moment(SED &base) { this->baseSED = &base; };
  bool monochromatic() { return this->baseSED->monochromatic(); };
  float lower_bound() { return this->baseSED->lower_bound(); };
  float upper_bound() { return this->baseSED->upper_bound(); };
  float value(float hnu) { return this->baseSED->value(hnu)*hnu*ev2erg/hplanck; };
};


int DualFLD::RadiationSource(HierarchyEntry *ThisGrid, int XrUv, float *Eta)
{

  // set dimension information
  int ghZl = (rank > 2) ? DEFAULT_GHOST_ZONES : 0;
  int ghYl = (rank > 1) ? DEFAULT_GHOST_ZONES : 0;
  int ghXl = DEFAULT_GHOST_ZONES;
  int n3[] = {1, 1, 1};
  int dx[] = {1.0, 1.0, 1.0};
  for (int dim=0; dim<rank; dim++) {
    n3[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
            - ThisGrid->GridData->GetGridStartIndex(dim) + 1;
    dx[dim] = ( ThisGrid->GridData->GetGridRightEdge(dim) -
		ThisGrid->GridData->GetGridLeftEdge(dim) ) / n3[dim];
  }
  int x0len = n3[0] + 2*ghXl;
  int x1len = n3[1] + 2*ghYl;
  int x2len = n3[2] + 2*ghZl;
  float x0L = Temp->GridData->GetGridLeftEdge(0);
  float x1L = Temp->GridData->GetGridLeftEdge(1);
  float x2L = Temp->GridData->GetGridLeftEdge(2);
  float lUn = (LenUnits + LenUnits0)*0.5;
  
  // access chemistry fields
  float *HI=NULL, *HeI=NULL, *HeII=NULL;
  HI   = ThisGrid->GridData->AccessHIDensity();
  HeI  = ThisGrid->GridData->AccessHeIDensity();
  HeII = ThisGrid->GridData->AccessHeIIDensity();

  // check that required field data exists
  if (HI==NULL)
    ENZO_FAIL("DualFLD::Opacity ERROR: no HI array!");
  if (Opacity==NULL)
    ENZO_FAIL("DualFLD::Opacity ERROR: no Opacity array!");
  if (RadiativeTransferHydrogenOnly == FALSE) {
    if (HeI==NULL)
      ENZO_FAIL("DualFLD::Opacity ERROR: no HeI array!");
    if (HeII==NULL)
      ENZO_FAIL("DualFLD::Opacity ERROR: no HeII array!");
  }


  ////////////////////////////
  // compute emissivity field

  // initialize field to zero
  for (i=0; i<x0len*x1len*x2len; i++)  Eta[i] = 0.0;
  
  // initialize constants and reusable local variables
  float dV = 1.0;
  for (int dim=0; dim<rank; dim++) 
    dV *= (dx[dim]*lUn);       // cell volume (proper)
  float h_nu0 = 13.6*ev2erg;   // ionization energy of HI [ergs]
  float h_nu1 = 24.6*ev2erg;   // ionization energy of HeI [ergs]
  float h_nu2 = 54.4*ev2erg;   // ionization energy of HeII [ergs]
  int i, j, k;
  float cellZl, cellZr, cellYl, cellYr, cellXl, cellXr, cellXc, cellYc, cellZc;


  // iterate over sources, adding emissivity to this field at specified location
  for (int isrc=0; isrc<NumSources; isrc++) {
    
    // all sources have radius of one cell; count number of cells to receive source
    int num_cells = 0;
    for (k=ghZl; k<n3[2]+ghZl; k++) {
      cellZc = x2L + (k-ghZl+0.5)*dx[2];        // z-center (comoving) for this cell
      for (j=ghYl; j<n3[1]+ghYl; j++) {
	cellYc = x1L + (j-ghYl+0.5)*dx[1];      // y-center (comoving) for this cell
	for (i=ghXl; i<n3[0]+ghXl; i++) {
	  cellXc = x0L + (i-ghXl+0.5)*dx[0];    // x-center (comoving) for this cell
	  if ( (fabs(cellXc-SourceLocation[ibin][0]) < dx[0]) &&
	       (fabs(cellYc-SourceLocation[ibin][1]) < dx[1]) &&
	       (fabs(cellZc-SourceLocation[ibin][2]) < dx[2]) )
	    num_cells++;                        // cell is within source region
	} // x-loop
      } // y-loop
    } // z-loop

    // equi-partition energy among affected cells
    float cell_energy = SourceGroupEnergy[isrc][XrUv] / num_cells / dV;
    for (k=ghZl; k<n3[2]+ghZl; k++) {
      cellZc = x2L + (k-ghZl+0.5)*dx[2];        // z-center (comoving) for this cell
      for (j=ghYl; j<n3[1]+ghYl; j++) {
	cellYc = x1L + (j-ghYl+0.5)*dx[1];      // y-center (comoving) for this cell
	for (i=ghXl; i<n3[0]+ghXl; i++) {
	  cellXc = x0L + (i-ghXl+0.5)*dx[0];    // x-center (comoving) for this cell
	  if ( (fabs(cellXc-SourceLocation[ibin][0]) < dx[0]) &&
	       (fabs(cellYc-SourceLocation[ibin][1]) < dx[1]) &&
	       (fabs(cellZc-SourceLocation[ibin][2]) < dx[2]) )
	    Eta[(k*x1len + j)*x0len + i] += cell_energy;
	} // x-loop
      } // y-loop
    } // z-loop

  } // end loop over sources from input parameters


	
  /////////////
  // add emissivity based on ProblemType, if desired
	
  switch (ProblemType) {

  case 412:    // emissivity flux along x=0 wall (NGammaDot photons/s/cm^2)
	  
    // only relevant for radiation groups (not individual frequencies)
    if (!FieldMonochromatic[XrUv]) {

      // place ionization sources along left wall (if on this subdomain)
      if (x0L == 0.0) {

	// set up source SED and the "moment" SED
	BlackbodySED tmp_src(1.0e5);
	BB_Moment tmp_src_m(tmp_src);

	// compute integral of SED over all ionizing frequencies
	float total_integral = SED_integral(tmp_src, 13.6, -1.0, true);

	// compute integral of moment SED over this band
	float band_integral = SED_integral(tmp_src_m, FrequencyBand[XrUv][0], 
					   FrequencyBand[XrUv][1], true);

	// determine this group's portion of total blackbody emissivity
	float wall_energy = 1e6 * hplanck * band_integral 
	  / total_integral / dx[0] / lUn;

	// place along wall (i=ghXl)
	for (k=ghZl; k<n3[2]+ghZl; k++)
	  for (j=ghYl; j<n3[1]+ghYl; j++)
	    Eta[(k*x1len + j)*x0len + ghXl] = wall_energy;
      }
    }
    break;
	  
  }  // switch (ProblemType)
	
  return SUCCESS;
}

#endif
