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
/  Radiation Spectrum Integration routine 
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Computes the UV and Xray radiation spectrum integrals
/              int_{nu0}^{inf} chiUV(nu) dnu
/              int_{nu0}^{inf} chiUV(nu)*sigHI(nu) dnu
/              int_{nu0}^{inf} chiUV(nu)*sigHI(nu)/nu dnu
/              int_{nu1}^{inf} chiUV(nu)*sigHeI(nu) dnu
/              int_{nu1}^{inf} chiUV(nu)*sigHeI(nu)/nu dnu
/              int_{nu2}^{inf} chiUV(nu)*sigHeII(nu) dnu
/              int_{nu2}^{inf} chiUV(nu)*sigHeII(nu)/nu dnu
/           where nu0_* is the ionization threshold of the relevant 
/           chemical species, chiUV(nu) is the spectrum of the relevant 
/           radiation energy density, sigHI(nu) is the ionization cross 
/           section of HI, sigHeI(nu) is the ionization cross section of 
/           HeI, and sigHeII(nu) is the ionization cross section of HeII.
/           The Xray integrals are defined similarly.
/
************************************************************************/
#ifdef TRANSFER
#include "DualFLD.h"
#include "phys_constants.h"
#include "CrossSectionSED.h"
#include "BlackbodySED.h"
#include "MonochromaticSED.h"
#include "PowerLawSED.h"
#include "PopIISED.h"
 
// function prototypes
float SED_integral(SED &sed, float a, float b, bool convertHz);

// SED types for computing cross-section integrals 
class CrossSection_type2 : public virtual SED {
 private:
  SED *radSED;                  // radiation SED
  CrossSectionSED *chemSED;     // cross section SED
 public:
  CrossSection_type2(SED &rSED, CrossSectionSED &cSED) { 
    radSED = &rSED; chemSED = &cSED; };
  bool monochromatic() { return chemSED->monochromatic(); };
  float lower_bound() { return max(chemSED->lower_bound(),radSED->lower_bound()); };
  float upper_bound() { return min(chemSED->upper_bound(),radSED->upper_bound()); };
  float value(float hnu) { 
    return ((radSED->value(hnu))*(chemSED->value(hnu))); 
  };
};

class CrossSection_type3 : public virtual SED {
 private:
  SED *radSED;                  // radiation SED
  CrossSectionSED *chemSED;     // cross section SED
 public:
  CrossSection_type3(SED &rSED, CrossSectionSED &cSED) { 
    radSED = &rSED; chemSED = &cSED; };
  bool monochromatic() { return chemSED->monochromatic(); };
  float lower_bound() { return max(chemSED->lower_bound(),radSED->lower_bound()); };
  float upper_bound() { return min(chemSED->upper_bound(),radSED->upper_bound()); };
  float value(float hnu) { 
    return ((radSED->value(hnu))*(chemSED->value(hnu))/(hnu*ev2erg/hplanck)); 
  };
};




int DualFLD::ComputeRadiationIntegrals() {

  if (debug)  printf("Entering DualFLD::ComputeRadiationIntegrals\n");

  // local variables
  float epsilon = 1.0;                   // floating point roundoff
  while ((1.0 + epsilon) > 1.0)  epsilon*=0.5;

  // first handle Xray radiation integrals
  if (UseXray) {

    // create SED object for each field
    switch (XrSpectrum) {
    case -1:     // monochromatic SED at hnu = Frequency [eV]
      Xr_SED = new MonochromaticSED(XrFrequency);
      break;
    case 0:      // power law SED, with exponent stored in Frequency (base is HI ionization threshold)
      Xr_SED = new PowerLawSED(XrFrequency, hnu0_HI);
      break;
    case 1:      // blackbody spectrum at T = 1e5 K
      Xr_SED = new BlackbodySED(1.0e5);
      break;
    case 2:      // PopII SED
      Xr_SED = new PopIISED();
      break;
    }

    // create SED objects for each chemical species and integrand
    CrossSectionSED    sHI(  0);
    CrossSection_type2 sHI2(  *Xr_SED, sHI);
    CrossSection_type3 sHI3(  *Xr_SED, sHI);
    CrossSectionSED    sHeI( 1);
    CrossSection_type2 sHeI2( *Xr_SED, sHeI);
    CrossSection_type3 sHeI3( *Xr_SED, sHeI);
    CrossSectionSED    sHeII(2);
    CrossSection_type2 sHeII2(*Xr_SED, sHeII);
    CrossSection_type3 sHeII3(*Xr_SED, sHeII);

    // perform integration
    if (XrSpectrum == -1) {   // monochromatic; integrals are just evaluation points

      // spectrum integral
      RadIntXr[0] = 1.0;

      // HI integrals
      float hnu_eval = XrFrequency*(1.0 + 2.0*epsilon);
      RadIntXr[1] = sHI2.value(hnu_eval);
      RadIntXr[2] = sHI3.value(hnu_eval);

      // HeI integrals
      RadIntXr[3] = sHeI2.value(hnu_eval);
      RadIntXr[4] = sHeI3.value(hnu_eval);

      // HeII integrals
      RadIntXr[5] = sHeII2.value(hnu_eval);
      RadIntXr[6] = sHeII3.value(hnu_eval);

    } else {                 // non-monochromatic

      // spectrum integral
      RadIntXr[0] = SED_integral(*Xr_SED, hnu0_HI, -1.0, true);
     
      // HI integrals
      RadIntXr[1] = SED_integral(sHI2, hnu0_HI, -1.0, true);
      RadIntXr[2] = SED_integral(sHI3, hnu0_HI, -1.0, true);

      // HeI integrals
      RadIntXr[3] = SED_integral(sHeI2, hnu0_HeI, -1.0, true);
      RadIntXr[4] = SED_integral(sHeI3, hnu0_HeI, -1.0, true);

      // HeI integrals
      RadIntXr[5] = SED_integral(sHeII2, hnu0_HeII, -1.0, true);
      RadIntXr[6] = SED_integral(sHeII3, hnu0_HeII, -1.0, true);

    }

    if (debug) {
      printf("  Computed X-ray Radiation Integrals:\n");
      printf("    intChiE          = %22.16e\n", RadIntXr[0]);
      printf("    intChiESigHI     = %22.16e\n", RadIntXr[1]);
      printf("    intChiESigHInu   = %22.16e\n", RadIntXr[2]);
      printf("    intChiESigHeI    = %22.16e\n", RadIntXr[3]);
      printf("    intChiESigHeInu  = %22.16e\n", RadIntXr[4]);
      printf("    intChiESigHeII   = %22.16e\n", RadIntXr[5]);
      printf("    intChiESigHeIInu = %22.16e\n", RadIntXr[6]);
    }
  }


  // second handle UV radiation integrals
  if (UseUV) {

    // create SED object for each field
    switch (UVSpectrum) {
    case -1:     // monochromatic SED at hnu = Frequency [eV]
      UV_SED = new MonochromaticSED(UVFrequency);
      break;
    case 0:      // power law SED, with exponent stored in Frequency (base is HI ionization threshold)
      UV_SED = new PowerLawSED(UVFrequency, hnu0_HI);
      break;
    case 1:      // blackbody spectrum at T = 1e5 K
      UV_SED = new BlackbodySED(1.0e5);
      break;
    case 2:      // PopII SED
      UV_SED = new PopIISED();
      break;
    }

    // create SED objects for each chemical species and integrand
    CrossSectionSED    sHI(  0);
    CrossSection_type2 sHI2(  *UV_SED, sHI);
    CrossSection_type3 sHI3(  *UV_SED, sHI);
    CrossSectionSED    sHeI( 1);
    CrossSection_type2 sHeI2( *UV_SED, sHeI);
    CrossSection_type3 sHeI3( *UV_SED, sHeI);
    CrossSectionSED    sHeII(2);
    CrossSection_type2 sHeII2(*UV_SED, sHeII);
    CrossSection_type3 sHeII3(*UV_SED, sHeII);

    // perform integration
    if (UVSpectrum == -1) {   // monochromatic; integrals are just evaluation points

      // spectrum integral
      RadIntUV[0] = 1.0;

      // HI integrals
      float hnu_eval = UVFrequency*(1.0 + 2.0*epsilon);
      RadIntUV[1] = sHI2.value(hnu_eval);
      RadIntUV[2] = sHI3.value(hnu_eval);

      // HeI integrals
      RadIntUV[3] = sHeI2.value(hnu_eval);
      RadIntUV[4] = sHeI3.value(hnu_eval);

      // HeII integrals
      RadIntUV[5] = sHeII2.value(hnu_eval);
      RadIntUV[6] = sHeII3.value(hnu_eval);

    } else {                 // non-monochromatic

      // spectrum integral
      RadIntUV[0] = SED_integral(*UV_SED, hnu0_HI, -1.0, true);
     
      // HI integrals
      RadIntUV[1] = SED_integral(sHI2, hnu0_HI, -1.0, true);
      RadIntUV[2] = SED_integral(sHI3, hnu0_HI, -1.0, true);

      // HeI integrals
      RadIntUV[3] = SED_integral(sHeI2, hnu0_HeI, -1.0, true);
      RadIntUV[4] = SED_integral(sHeI3, hnu0_HeI, -1.0, true);

      // HeI integrals
      RadIntUV[5] = SED_integral(sHeII2, hnu0_HeII, -1.0, true);
      RadIntUV[6] = SED_integral(sHeII3, hnu0_HeII, -1.0, true);

    }

    if (debug) {
      printf("  Computed UV Radiation Integrals:\n");
      printf("    intChiE          = %22.16e\n", RadIntUV[0]);
      printf("    intChiESigHI     = %22.16e\n", RadIntUV[1]);
      printf("    intChiESigHInu   = %22.16e\n", RadIntUV[2]);
      printf("    intChiESigHeI    = %22.16e\n", RadIntUV[3]);
      printf("    intChiESigHeInu  = %22.16e\n", RadIntUV[4]);
      printf("    intChiESigHeII   = %22.16e\n", RadIntUV[5]);
      printf("    intChiESigHeIInu = %22.16e\n", RadIntUV[6]);
    }
  }

  return SUCCESS;
}

#endif
