/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  PopII radiation SED instantiation of the SED base class.
/
/  written by: Daniel Reynolds
/  date:       November, 2014
/  modified1:  
/
/  PURPOSE: This implements a PopII SED, as described in Ricotti etal, 
/           ApJ 575:33-48, 2002.
/
************************************************************************/

#ifndef POPII_SED_DEFINED__
#define POPII_SED_DEFINED__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "phys_constants.h"
#include "SED.h"

class PopIISED : public virtual SED {

 public:

  // required functions
  bool monochromatic();
  float lower_bound();
  float upper_bound();
  float value(float hnu);

};
  
#endif
