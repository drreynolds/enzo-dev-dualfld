/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Power law radiation SED instantiation of the SED base class.
/
/  written by: Daniel Reynolds
/  date:       November, 2014
/  modified1:  
/
/  PURPOSE: This is an example SED class implementation, that is used 
/           for test problems.
/
************************************************************************/

#include "PowerLawSED.h"

// constructor
PowerLawSED::PowerLawSED(float exponent_, float hnu0_) { 
  this->exponent = exponent_; 
  this->hnu0 = hnu0_;
};

// monochromatic return function
bool PowerLawSED::monochromatic() { 
  return false;  // not monochromatic
}

// lower bound function
float PowerLawSED::lower_bound() { 
  return 0.0;    // lower bound of hnu=0
}

// upper bound function
float PowerLawSED::upper_bound() { 
  return -1.0;   // no upper bound
}

// main SED function
float PowerLawSED::value(float hnu) {
  return (POW((hnu/hnu0),exponent));
}
