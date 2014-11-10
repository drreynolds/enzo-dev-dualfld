/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Power Law SED instantiation of the SED base class.
/
/  written by: Daniel Reynolds
/  date:       November, 2014
/  modified1:  
/
/  PURPOSE: This is an example SED class implementation, that is used 
/           for test problems.
/
************************************************************************/

#ifndef POWERLAW_SED_DEFINED__
#define POWERLAW_SED_DEFINED__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "phys_constants.h"
#include "SED.h"

class PowerLawSED : public virtual SED {

 private:

  // SED-defining data
  float exponent;      // exponent for power-law
  float hnu0;          // base frequency for power-law

 public:

  // constructor
  PowerLawSED(float exponent_, float hnu0_);

  // required functions
  bool monochromatic();
  float lower_bound();
  float upper_bound();
  float value(float hnu);

};
  
#endif
