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
/  EnforceBoundary routine
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Enforces boundary conditions on a FLD problem vector.
/
/           Note: Neumann values are enforced on the first 
/                 layer of ghost zones using a first-order central 
/                 difference approximation to the first (outward-normal) 
/                 derivative.
/           Note: Since the internal radiation variables are comoving 
/                 and normalized, we renormalize the boundary conditions 
/                 as they are enforced to match the internal units.
/
************************************************************************/
#ifdef TRANSFER
#include "DualFLD.h"


int DualFLD::EnforceBoundary(HierarchyEntry *ThisGrid) {

  // access relevant fields
  float *UV, *Xr;
  Xr = ThisGrid->GridData->AccessRadiationFrequency0();
  UV = ThisGrid->GridData->AccessRadiationFrequency1();

  // check that required field data exists
  if (UseXray && Xr==NULL)
    ENZO_FAIL("DualFLD::FillRates ERROR: missing RadiationFrequency0 BaryonField");
  if (UseUV && UV==NULL)
    ENZO_FAIL("DualFLD::FillRates ERROR: missing RadiationFrequency1 BaryonField");

  // temporary variables
  int i, j, k, idx, i2, j2, k2, idx2, idxbc;

  // enforce boundary conditions on Xray radiation field
  if (UseXray) {
    float *udata = Xr;
    float EUnits = XrUnits;

    float dxa = dx[0]*LenUnits;
    // x0 left boundary
    //   Dirichlet
    if (OnBdry[0][0] && (XrBdryType[0][0]==1)) {
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<GhDims[0][0]; i++) {
	    idxbc = k*LocDims[1] + j;
	    idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i;
	    udata[idx] = XrBdryVals[0][0][idxbc]/EUnits;
	  }
    }
    //   Neumann
    if (OnBdry[0][0] && (XrBdryType[0][0]==2)) {
      i = -1;  i2 = i+1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	  idx2 = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i2+GhDims[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dxa*XrBdryVals[0][0][idxbc]/EUnits;
	}
    }
    
    // x0 right boundary
    //   Dirichlet
    if (OnBdry[0][1] && (XrBdryType[0][1]==1)) {
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=ArrDims[0]-GhDims[0][0]; i<ArrDims[0]; i++) {
	    idxbc = k*LocDims[1] + j;
	    idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i;
	    udata[idx] = XrBdryVals[0][1][idxbc]/EUnits;
	  }
    }
    //   Neumann
    if (OnBdry[0][1] && (XrBdryType[0][1]==2)) {
      i = LocDims[0];  i2 = i-1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	  idx2 = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i2+GhDims[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dxa*XrBdryVals[0][1][idxbc]/EUnits;
	}
    }

    if (rank > 1) {
      float dya = dx[1]*LenUnits;
      // x1 left boundary
      //   Dirichlet
      if (OnBdry[1][0] && (XrBdryType[1][0]==1)) {
	for (k=0; k<LocDims[2]; k++)
	  for (j=0; j<GhDims[1][0]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = ((k+GhDims[2][0])*ArrDims[1] + j)*ArrDims[0] + i+GhDims[0][0];
	      idxbc = i*LocDims[2] + k;
	      udata[idx] = XrBdryVals[1][0][idxbc]/EUnits;
	    }
      }
      //   Neumann
      if (OnBdry[1][0] && (XrBdryType[1][0]==2)) {
	j = -1;  j2 = j+1;
	for (k=0; k<LocDims[2]; k++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idx2 = ((k+GhDims[2][0])*ArrDims[1] + j2+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = udata[idx2] + dya*XrBdryVals[1][0][idxbc]/EUnits;
	  }
      }
      
      // x1 right boundary
      //   Dirichlet
      if (OnBdry[1][1] && (XrBdryType[1][1]==1)) {
	for (k=0; k<LocDims[2]; k++)
	  for (j=ArrDims[1]-GhDims[1][0]; j<ArrDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = ((k+GhDims[2][0])*ArrDims[1] + j)*ArrDims[0] + i+GhDims[0][0];
	      idxbc = i*LocDims[2] + k;
	      udata[idx] = XrBdryVals[1][1][idxbc]/EUnits;
	    }
      }
      //   Neumann
      if (OnBdry[1][1] && (XrBdryType[1][1]==2)) {
	j = LocDims[1];  j2 = j-1;
	for (k=0; k<LocDims[2]; k++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idx2 = ((k+GhDims[2][0])*ArrDims[1] + j2+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = udata[idx2] + dya*XrBdryVals[1][1][idxbc]/EUnits;
	  }
      }
    }  // end if rank > 1
     
    if (rank > 2) {
      float dza = dx[2]*LenUnits;
      // x2 left boundary
      //   Dirichlet
      if (OnBdry[2][0] && (XrBdryType[2][0]==1)) {
	for (k=0; k<GhDims[2][0]; k++)
	  for (j=0; j<LocDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = (k*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	      idxbc = j*LocDims[0] + i;
	      udata[idx] = XrBdryVals[2][0][idxbc]/EUnits;
	    }
      }
      //   Neumann
      if (OnBdry[2][0] && (XrBdryType[2][0]==2)) {
	k = -1;  k2 = k+1;
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idx2 = ((k2+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = udata[idx2] + dza*XrBdryVals[2][0][idxbc]/EUnits;
	  }
      }
      
      // x2 right boundary
      //   Dirichlet
      if (OnBdry[2][1] && (XrBdryType[2][1]==1)) {
	for (k=ArrDims[2]-GhDims[2][0]; k<ArrDims[2]; k++)
	  for (j=0; j<LocDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = (k*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	      idxbc = j*LocDims[0] + i;
	      udata[idx] = XrBdryVals[2][1][idxbc]/EUnits;
	    }
      }
      //   Neumann
      if (OnBdry[2][1] && (XrBdryType[2][1]==2)) {
	k = LocDims[2];  k2 = k-1;
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idx2 = ((k2+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = udata[idx2] + dza*XrBdryVals[2][1][idxbc]/EUnits;
	  }
      }
    }  // end if rank > 2

  }  // if (UseXray)



  // enforce boundary conditions on UV radiation field
  if (UseUV) {
    float *udata = UV;
    float EUnits = UVUnits;

    float dxa = dx[0]*LenUnits;
    // x0 left boundary
    //   Dirichlet
    if (OnBdry[0][0] && (UVBdryType[0][0]==1)) {
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<GhDims[0][0]; i++) {
	    idxbc = k*LocDims[1] + j;
	    idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i;
	    udata[idx] = UVBdryVals[0][0][idxbc]/EUnits;
	  }
    }
    //   Neumann
    if (OnBdry[0][0] && (UVBdryType[0][0]==2)) {
      i = -1;  i2 = i+1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	  idx2 = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i2+GhDims[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dxa*UVBdryVals[0][0][idxbc]/EUnits;
	}
    }
    
    // x0 right boundary
    //   Dirichlet
    if (OnBdry[0][1] && (UVBdryType[0][1]==1)) {
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++)
	  for (i=ArrDims[0]-GhDims[0][0]; i<ArrDims[0]; i++) {
	    idxbc = k*LocDims[1] + j;
	    idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i;
	    udata[idx] = UVBdryVals[0][1][idxbc]/EUnits;
	  }
    }
    //   Neumann
    if (OnBdry[0][1] && (UVBdryType[0][1]==2)) {
      i = LocDims[0];  i2 = i-1;
      for (k=0; k<LocDims[2]; k++)
	for (j=0; j<LocDims[1]; j++) {
	  idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	  idx2 = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i2+GhDims[0][0];
	  idxbc = k*LocDims[1] + j;
	  udata[idx] = udata[idx2] + dxa*UVBdryVals[0][1][idxbc]/EUnits;
	}
    }

    if (rank > 1) {
      float dya = dx[1]*LenUnits;
      // x1 left boundary
      //   Dirichlet
      if (OnBdry[1][0] && (UVBdryType[1][0]==1)) {
	for (k=0; k<LocDims[2]; k++)
	  for (j=0; j<GhDims[1][0]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = ((k+GhDims[2][0])*ArrDims[1] + j)*ArrDims[0] + i+GhDims[0][0];
	      idxbc = i*LocDims[2] + k;
	      udata[idx] = UVBdryVals[1][0][idxbc]/EUnits;
	    }
      }
      //   Neumann
      if (OnBdry[1][0] && (UVBdryType[1][0]==2)) {
	j = -1;  j2 = j+1;
	for (k=0; k<LocDims[2]; k++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idx2 = ((k+GhDims[2][0])*ArrDims[1] + j2+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = udata[idx2] + dya*UVBdryVals[1][0][idxbc]/EUnits;
	  }
      }
      
      // x1 right boundary
      //   Dirichlet
      if (OnBdry[1][1] && (UVBdryType[1][1]==1)) {
	for (k=0; k<LocDims[2]; k++)
	  for (j=ArrDims[1]-GhDims[1][0]; j<ArrDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = ((k+GhDims[2][0])*ArrDims[1] + j)*ArrDims[0] + i+GhDims[0][0];
	      idxbc = i*LocDims[2] + k;
	      udata[idx] = UVBdryVals[1][1][idxbc]/EUnits;
	    }
      }
      //   Neumann
      if (OnBdry[1][1] && (UVBdryType[1][1]==2)) {
	j = LocDims[1];  j2 = j-1;
	for (k=0; k<LocDims[2]; k++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idx2 = ((k+GhDims[2][0])*ArrDims[1] + j2+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idxbc = i*LocDims[2] + k;
	    udata[idx] = udata[idx2] + dya*UVBdryVals[1][1][idxbc]/EUnits;
	  }
      }
    }  // end if rank > 1
     
    if (rank > 2) {
      float dza = dx[2]*LenUnits;
      // x2 left boundary
      //   Dirichlet
      if (OnBdry[2][0] && (UVBdryType[2][0]==1)) {
	for (k=0; k<GhDims[2][0]; k++)
	  for (j=0; j<LocDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = (k*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	      idxbc = j*LocDims[0] + i;
	      udata[idx] = UVBdryVals[2][0][idxbc]/EUnits;
	    }
      }
      //   Neumann
      if (OnBdry[2][0] && (UVBdryType[2][0]==2)) {
	k = -1;  k2 = k+1;
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idx2 = ((k2+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = udata[idx2] + dza*UVBdryVals[2][0][idxbc]/EUnits;
	  }
      }
      
      // x2 right boundary
      //   Dirichlet
      if (OnBdry[2][1] && (UVBdryType[2][1]==1)) {
	for (k=ArrDims[2]-GhDims[2][0]; k<ArrDims[2]; k++)
	  for (j=0; j<LocDims[1]; j++)
	    for (i=0; i<LocDims[0]; i++) {
	      idx = (k*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	      idxbc = j*LocDims[0] + i;
	      udata[idx] = UVBdryVals[2][1][idxbc]/EUnits;
	    }
      }
      //   Neumann
      if (OnBdry[2][1] && (UVBdryType[2][1]==2)) {
	k = LocDims[2];  k2 = k-1;
	for (j=0; j<LocDims[1]; j++)
	  for (i=0; i<LocDims[0]; i++) {
	    idx = ((k+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idx2 = ((k2+GhDims[2][0])*ArrDims[1] + j+GhDims[1][0])*ArrDims[0] + i+GhDims[0][0];
	    idxbc = j*LocDims[0] + i;
	    udata[idx] = udata[idx2] + dza*UVBdryVals[2][1][idxbc]/EUnits;
	  }
      }
    }  // end if rank > 2

  }  // if (UseUV)

  return SUCCESS;
}

#endif
