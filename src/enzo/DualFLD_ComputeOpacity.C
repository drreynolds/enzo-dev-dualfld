/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Dual Flux-Limited Diffusion Solver, opacity field calculation routine 
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Computes the opacity field throughout the domain, storing it 
/           in the provided array, Opacity (assumed to be the size of a 
/           BaryonField).  All non-DualFLD required values are extracted 
/           from the provided HierarchyEntry, ThisGrid.
/
************************************************************************/
#ifdef TRANSFER
#include "DualFLD.h"


int DualFLD::ComputeOpacity(HierarchyEntry *ThisGrid, int XrUv, float *Opacity)
{

  // set dimension information
  int ghZl = (rank > 2) ? NumberOfGhostZones : 0;
  int ghYl = (rank > 1) ? NumberOfGhostZones : 0;
  int ghXl = NumberOfGhostZones;
  int n3[] = {1, 1, 1};
  for (int dim=0; dim<rank; dim++)
    n3[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
            - ThisGrid->GridData->GetGridStartIndex(dim) + 1;
  int x0len = n3[0] + 2*ghXl;
  int x1len = n3[1] + 2*ghYl;
  int x2len = n3[2] + 2*ghZl;
  
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
  // compute opacity field

  // local variables to be reused within loop
  int i;
  float HIconst, HeIconst, HeIIconst;
  if (XrUv) {   // UV opacity
    HIconst   = RadIntUV[1] / RadIntUV[0];
    HeIconst  = RadIntUV[3] / RadIntUV[0] / 4.0;
    HeIIconst = RadIntUV[5] / RadIntUV[0] / 4.0;
  } else {      // Xray opacity
    HIconst   = RadIntXr[1] / RadIntXr[0];
    HeIconst  = RadIntXr[3] / RadIntXr[0] / 4.0;
    HeIIconst = RadIntXr[5] / RadIntXr[0] / 4.0;
  }

  // Hydrogen-only calculation
  if (RadiativeTransferHydrogenOnly) {
    for (i=0; i<x0len*x1len*x2len; i++) 
      Opacity[i] = HI[i]*HIconst;
    
  // Hydrogen + Helium calculation
  } else {
    for (i=0; i<x0len*x1len*x2len; i++) 
      Opacity[i] = HI[i]*HIconst + HeI[i]*HeIconst + HeII[i]*HeIIconst;
  }
  
  return SUCCESS;
}

#endif  /* TRANSFER */
