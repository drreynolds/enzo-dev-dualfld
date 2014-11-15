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
/  FillRates routine.
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Fills the photo-ionization and photo-heating arrays used by
/           chemistry and cooling routines using time-averaged internal
/           values for the radiation.
/ 
/  NOTE: In order to save on memory, the photo-heating rates are 
/        combined into a single rate, and scaled by the current number 
/        density of HI, to be later unpacked by rescaling back with HI.  
/        This loses accuracy in the case that during chemistry 
/        subcycling the chemistry changes significantly, since we retain 
/        the initial rate scaling but use updated HI values in rescaling.
/
************************************************************************/
#ifdef TRANSFER
#include "DualFLD.h"
#include "phys_constants.h"


int DualFLD::FillRates(HierarchyEntry *ThisGrid) {

//   if (debug)
//     printf("Entering DualFLD::FillRates routine\n");

  // access relevant fields
  float *rho, *HI, *HeI, *HeII, *phHI, *phHeI, *phHeII, *photogamma, *dissH2I, *UV, *Xr;
  rho = HI = HeI = HeII = phHI = phHeI = phHeII = photogamma = dissH2I = NULL;
  rho  = ThisGrid->GridData->AccessDensity();
  HI   = ThisGrid->GridData->AccessHIDensity();
  HeI  = ThisGrid->GridData->AccessHeIDensity();
  HeII = ThisGrid->GridData->AccessHeIIDensity();
  phHI = ThisGrid->GridData->AccessKPhHI();
  phHeI = ThisGrid->GridData->AccessKPhHeI();
  phHeII = ThisGrid->GridData->AccessKPhHeII();
  photogamma = ThisGrid->GridData->AccessPhotoGamma();
  dissH2I = ThisGrid->GridData->AccessKDissH2I();
  Xr = ThisGrid->GridData->AccessRadiationFrequency0();
  UV = ThisGrid->GridData->AccessRadiationFrequency1();

  // check that required field data exists
  if (rho==NULL)
    ENZO_FAIL("DualFLD::FillRates ERROR: missing Density BaryonField");
  if (HI==NULL)
    ENZO_FAIL("DualFLD::FillRates ERROR: missing HIDensity BaryonField");
  if (phHI==NULL)
    ENZO_FAIL("DualFLD::FillRates ERROR: missing KPhHI BaryonField");
  if (photogamma==NULL)
    ENZO_FAIL("DualFLD::FillRates ERROR: missing PhotoGamma BaryonField");
  if (UseXray && Xr==NULL)
    ENZO_FAIL("DualFLD::FillRates ERROR: missing RadiationFrequency0 BaryonField");
  if (UseUV && UV==NULL)
    ENZO_FAIL("DualFLD::FillRates ERROR: missing RadiationFrequency1 BaryonField");
  if (RadiativeTransferHydrogenOnly == FALSE) {
    if (HeI==NULL)
      ENZO_FAIL("DualFLD::FillRates ERROR: missing HeIDensity BaryonField");
    if (HeII==NULL)
      ENZO_FAIL("DualFLD::FillRates ERROR: missing HeIIDensity BaryonField");
    if (phHeI==NULL)
      ENZO_FAIL("DualFLD::FillRates ERROR: missing KPhHeI BaryonField");
    if (phHeII==NULL)
      ENZO_FAIL("DualFLD::FillRates ERROR: missing KPhHeII BaryonField");
  }
  if (MultiSpecies > 1)
    if (dissH2I==NULL)
      ENZO_FAIL("DualFLD::FillRates ERROR: missing KDissH2I BaryonField");

  // set some physical constants
  float dom      = DenUnits*a*a*a/mh;
  float tbase1   = TimeUnits;
  float xbase1   = LenUnits/a/aUnits;
  float dbase1   = DenUnits*a*a*a*aUnits*aUnits*aUnits;
  float coolunit = aUnits*aUnits*aUnits*aUnits*aUnits * xbase1*xbase1
                 * mh*mh / tbase1/tbase1/tbase1 / dbase1;
  float rtunits  = ev2erg/TimeUnits/coolunit/dom;

  //float UVUn     = UVUnits;   // original
  //float XrUn     = XrUnits;
  float UVUn     = (UVUnits+UVUnits0)*0.5;   // arithmetic mean
  float XrUn     = (XrUnits+XrUnits0)*0.5;
  //float UVUn     = sqrt(UVUnits*UVUnits0);   // geometric mean
  //float XrUn     = sqrt(XrUnits*XrUnits0);
  //float UVUn     = 2.0*UVUnits*UVUnits0/(UVUnits+UVUnits0);   // harmonic mean
  //float XrUn     = 2.0*XrUnits*XrUnits0/(XrUnits+XrUnits0);

  // temporary vector 'sol' is not currently used, store electron fraction there
  float *xe = sol->GetData(0);

  // compute the size of the fields
  int i, dim, size=1;
  for (dim=0; dim<rank; dim++)  size *= ArrDims[dim];

  // fill electron fraction over grid
  float HFrac = CoolData.HydrogenFractionByMass;
  if (RadiativeTransferHydrogenOnly == 1)
    for (i=0; i<size; i++)
      xe[i] = max(1.0 - HI[i]/(rho[i]*HFrac), 1e-4);   // xe = nHII / nH
  else
    for (i=0; i<size; i++)   // xe = (nHII + nHeII/4 + nHeIII/2) / (nH + 2*nHe)
      xe[i] = max((rho[i]*(HFrac+1.0)*0.5 - HI[i] - 0.5*HeI[i] - 0.25*HeII[i])
		  / (rho[i]*(HFrac+1.0)*0.5), 1e-4);

  // fill HI photo-ionization rate
  for (i=0; i<size; i++)  phHI[i] = 0.0;
  if (UseXray) {
    float pHIconstXr = clight*TimeUnits*XrUn*RadIntXr[1]/RadIntXr[0]/(hnu0_HI*ev2erg);
    for (i=0; i<size; i++)
      phHI[i] += Xr[i]*pHIconstXr*0.3908*pow((1.0-pow(xe[i], 0.4092)), 1.7592);
  }
  if (UseUV) {
    float pHIconstUV = clight*TimeUnits*UVUn*RadIntUV[2]/hplanck/RadIntUV[0];
    for (i=0; i<size; i++)  
      phHI[i] += UV[i]*pHIconstUV;
  }

  // fill HeI and HeII photo-ionization rates
  if (!RadiativeTransferHydrogenOnly) {
    for (i=0; i<size; i++)  phHeI[i] = 0.0;
    for (i=0; i<size; i++)  phHeII[i] = 0.0;
    if (UseXray) {
      float pHeIconstXr  = clight*TimeUnits*XrUn*RadIntXr[3]/RadIntXr[0]/(hnu0_HeI*ev2erg);
      for (i=0; i<size; i++)  
	phHeI[i] += Xr[i]*pHeIconstXr*0.0554*pow((1.0-pow(xe[i], 0.4614)), 1.6660);
    }
    if (UseUV) {
      float pHeIconstUV  = clight*TimeUnits*UVUn*RadIntUV[4]/hplanck/RadIntUV[0];
      float pHeIIconstUV = clight*TimeUnits*UVUn*RadIntUV[6]/hplanck/RadIntUV[0];
      for (i=0; i<size; i++)  phHeI[i] += UV[i]*pHeIconstUV;
      for (i=0; i<size; i++)  phHeII[i] += UV[i]*pHeIIconstUV;
    }
  }

  // fill photo-heating rate
  for (i=0; i<size; i++)  photogamma[i] = 0.0;
  if (UseXray) {
    float phScaleXr    = clight*TimeUnits*XrUn/RadIntXr[0]/VelUnits/VelUnits/mh/rtunits*0.9971;
    float GHIconstXr   = phScaleXr*(RadIntXr[1] - hnu0_HI*ev2erg/hplanck*RadIntXr[2]);
    float GHeIconstXr  = phScaleXr*(RadIntXr[3] - hnu0_HeI*ev2erg/hplanck*RadIntXr[4]);
    float GHeIIconstXr = phScaleXr*(RadIntXr[5] - hnu0_HeII*ev2erg/hplanck*RadIntXr[6]);
    if (RadiativeTransferHydrogenOnly) {
      for (i=0; i<size; i++)
	photogamma[i] += Xr[i]*GHIconstXr*(1.0 - pow(1.0 - pow(xe[i], 0.2663), 1.3163));
    } else {
      for (i=0; i<size; i++) {
	photogamma[i] += Xr[i]/HI[i]*(GHIconstXr*HI[i] + GHeIconstXr*HeI[i] + GHeIIconstXr*HeII[i])
	  *(1.0 - pow(1.0 - pow(xe[i], 0.2663), 1.3163));
      }
    }
  }
  if (UseUV) {
    float phScaleUV    = clight*TimeUnits*UVUn/RadIntUV[0]/VelUnits/VelUnits/mh/rtunits;
    float GHIconstUV   = phScaleUV*(RadIntUV[1] - hnu0_HI*ev2erg/hplanck*RadIntUV[2]);
    float GHeIconstUV  = phScaleUV*(RadIntUV[3] - hnu0_HeI*ev2erg/hplanck*RadIntUV[4]);
    float GHeIIconstUV = phScaleUV*(RadIntUV[5] - hnu0_HeII*ev2erg/hplanck*RadIntUV[6]);
    if (RadiativeTransferHydrogenOnly) {
      for (i=0; i<size; i++)  photogamma[i] += UV[i]*GHIconstUV;
    } else {
      for (i=0; i<size; i++)  
	photogamma[i] += UV[i]/HI[i]*(GHIconstUV*HI[i] + GHeIconstUV*HeI[i] + GHeIIconstUV*HeII[i]);
    }
  }

  // fill H2 dissociation rate (zero for grey FLD problems)
  if (MultiSpecies > 1) 
    for (i=0; i<size; i++)  dissH2I[i] = 0.0;

  // return success
  return SUCCESS;

}
#endif
