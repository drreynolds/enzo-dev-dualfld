/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Pop II radiation SED instantiation of the SED base class.
/
/  written by: Daniel Reynolds
/  date:       November, 2014
/  modified1:  
/
/  PURPOSE: This implements a PopII SED, as described in Ricotti etal, 
/           ApJ 575:33-48, 2002.
/
************************************************************************/
#include "PopIISED.h"
#include "phys_constants.h"

// monochromatic return function
bool PopIISED::monochromatic() { 
  return false;  // not monochromatic
}

// lower bound function
float PopIISED::lower_bound() { 
  return 0.0;    // lower bound of hnu=0
}

// upper bound function
float PopIISED::upper_bound() { 
  return 54.4*ev2erg/hplanck;   // HeII ionization threshold
}

// main SED function
float PopIISED::value(float hnu) {
  float nu = hnu*ev2erg/hplanck;       // convert frequency to Hz
  float nu0 = 13.6*ev2erg/hplanck;     // ionization threshold of Hydrogen (Hz)
  float nu1 = 2.5*nu0;                 // ionization of Wolf Reyet stars + HeliumI
  float nu2 = 4.0*nu0;                 // ionization threshold of HeliumII
  float GBbrk = 2.5;                   // parameters in Ricotti 2002 fig.4 fit

  // SED with photons above 4 Ryd truncated
  if (nu < nu1)
    return (1.0/nu0/POW(nu/nu0, 1.8));
  else if (nu < nu2)
    return (1.0/nu0/GBbrk/POW(nu/nu0, 1.8));
  else
    return 0.0;
}
