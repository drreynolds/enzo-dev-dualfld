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
/  Destructor routine
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Frees all memory allocated for the implicit FLD problem.
/
************************************************************************/
#ifdef TRANSFER
#include "DualFLD.h"


DualFLD::~DualFLD() {

//   if (debug)  printf("Entering DualFLD::destructor routine\n");

  // delete HYPRE objects
#ifdef USE_HYPRE
  HYPRE_StructStencilDestroy(stencil);
  HYPRE_StructGridDestroy(grid);
#endif

  // delete EnzoVectors and other internal arrays
  //   EnzoVectors require deleting the structure
  delete sol;

  //   arrays require deleting the array
#ifdef USE_HYPRE
  HYPRE_StructVectorDestroy(rhsvec);
  HYPRE_StructVectorDestroy(solvec);
  HYPRE_StructMatrixDestroy(P);
#endif
  delete[] HYPREbuff;

  // delete SED objects
  delete Xr_SED;
  delete UV_SED;

  // delete boundary condition arrays
  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++) {
      if (XrBdryVals[i][j] != NULL)  
	delete[] XrBdryVals[i][j];
      if (UVBdryVals[i][j] != NULL)  
	delete[] UVBdryVals[i][j];
    }

}
#endif
