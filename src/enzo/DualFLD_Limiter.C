/*****************************************************************************
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Dual Flux-Limited Diffusion Solver, flux-limiter evaluation routine 
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Computes the flux limiter for a given cell face.
/
************************************************************************/
#ifdef TRANSFER
#include "phys_constants.h"
#include "DualFLD.h"


float DualFLD::Limiter(float E1, float E2, float k1, float k2, 
		       float nUn, float lUn, float dxi)
{

  // limiter bounds
  float Rmin = MIN(FLD_Rmin/lUn, 1e-20);            // 1st is astro/cosmo, 2nd is lab frame
  float Dmax = MAX(FLD_Dmax * clight * lUn, 1e20);  // 1st is astro/cosmo, 2nd is lab frame
  float Eavg = MAX((E1 + E2)*0.5, FLD_Emin);        // arithmetic mean, bound from zero
  k1 = MAX(k1, FLD_Kmin);                           // bound away from zero since we'll divide
  k2 = MAX(k2, FLD_Kmin);
  Scalar kap = 2.0*k1*k2/(k1+k2)*nUn;               // harmonic mean
  Scalar R = MAX(dxi*fabs(E1-E2)/Eavg, Rmin);

  // compute limiter
  return (MIN(c_light/sqrt(9.0*kap*kap + R*R), Dmax));

}

#endif
