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
/  Time Step Computation Routine
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: Computes the rad-hydro time step size.  We note that this 
/           value affects the global hydrodynamics time step: 
/                 dt = min(dt_hydro,dt_radiation).
/           This routine is called with scaled arguments.
/
************************************************************************/
#ifdef TRANSFER
#include "DualFLD.h"
 

float DualFLD::ComputeTimeStep(EnzoVector *uold, EnzoVector *unew) {

  // If timeAccuracyXr or timeAccuracyUV is set, compute maximum time 
  // step as estimate allowing timeAccuracy* relative error
  float dt_est = huge_number;    // max time step (normalized units)
  float dt_est_Xr = huge_number;
  float dt_est_UV = huge_number;
  if (!UseXray)  timeAccuracyXr = huge_number;  // double-check that these are valid
  if (!UseUV)    timeAccuracyUV = huge_number;
  if ((timeAccuracyXr != huge_number) || (timeAccuracyUV != huge_number)) {

    // local variables
    float diff, w, tmp;
    int i, j, k, l;
    float safety = 1.0;
    float Err_floor = 1e-8;
    float den = ((theta < 0.51) && (theta > 0.49)) ? 2.0 : 1.0;

    // compute volume factor to account for fact that error estimate 
    // integrates over domain volume
    float Vol = 1.0;
    float dV  = 1.0;
    for (i=0; i<rank; i++) {
      Vol *= (DomainRightEdge[i]-DomainLeftEdge[i]);
      dV  *= dx[i];
    }
    float VolFac = (dtnorm > 0.0) ? pow(Vol,1.0/dtnorm) : 1.0;
      
    // update old error estimates
    float Err_old_Xr = Err_cur_Xr;
    float Err_old_UV = Err_cur_UV;
    Err_cur_Xr = Err_new_Xr;
    Err_cur_UV = Err_new_UV;
    
    // compute estimate of the change in each field:
    //       relerr_fac = || (unew - uold) / w ||_p
    //    with the scaling vector w given by
    //       w = dtfactor*[sqrt(|unew*uold|) + atol]
    //    and where we have the following parameters:
    //       p - norm choice (input), 0->max norm, otherwise the p-norm
    //       atol - 0.1 (assumes units are all normalized)
    float atol = 0.1;       // assumes variables are nearly normalized
    float *Eold, *Enew;

    // initialize variables
    int x0len = ArrDims[0];
    int x1len = ArrDims[1];
    Err_new_Xr = 0.0;
    Err_new_UV = 0.0;

    // perform local estimates for the Xray relative change
    if (timeAccuracyXr != huge_number) {
      Eold = uold->GetData(0);
      Enew = unew->GetData(0);
      if (dtnorm > 0.0) {
	for (k=GhDims[2][0]; k<LocDims[2]+GhDims[2][0]; k++) 
	  for (j=GhDims[1][0]; j<LocDims[1]+GhDims[1][0]; j++)
	    for (i=GhDims[0][0]; i<LocDims[0]+GhDims[0][0]; i++) {
	      w = sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
		          * Eold[(k*x1len + j)*x0len + i])) + atol;
	      diff = Enew[(k*x1len + j)*x0len + i] 
		   - Eold[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      Err_new_Xr += POW(tmp,dtnorm)*dV;
	    }
      } else {
	for (k=GhDims[2][0]; k<LocDims[2]+GhDims[2][0]; k++) 
	  for (j=GhDims[1][0]; j<LocDims[1]+GhDims[1][0]; j++)
	    for (i=GhDims[0][0]; i<LocDims[0]+GhDims[0][0]; i++) {
	      w = sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
			  * Eold[(k*x1len + j)*x0len + i])) + atol;
	      diff = Enew[(k*x1len + j)*x0len + i] 
  		   - Eold[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      Err_new_Xr = max(Err_new_Xr, tmp);
	    }
      }
    }

    // perform local estimates for UV relative change
    if (timeAccuracyUV != huge_number) {
      int UVint = (UseXray) ? 1 : 0;
      Eold = uold->GetData(UVint);
      Enew = unew->GetData(UVint);
      if (dtnorm > 0.0) {
	for (k=GhDims[2][0]; k<LocDims[2]+GhDims[2][0]; k++) 
	  for (j=GhDims[1][0]; j<LocDims[1]+GhDims[1][0]; j++)
	    for (i=GhDims[0][0]; i<LocDims[0]+GhDims[0][0]; i++) {
	      w = sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
			  * Eold[(k*x1len + j)*x0len + i])) + atol;
	      diff = Enew[(k*x1len + j)*x0len + i] 
		   - Eold[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      Err_new_UV += POW(tmp,dtnorm)*dV;
	    }
      } else {
	for (k=GhDims[2][0]; k<LocDims[2]+GhDims[2][0]; k++) 
	  for (j=GhDims[1][0]; j<LocDims[1]+GhDims[1][0]; j++)
	    for (i=GhDims[0][0]; i<LocDims[0]+GhDims[0][0]; i++) {
	      w = sqrt(fabs(Enew[(k*x1len + j)*x0len + i]
			  * Eold[(k*x1len + j)*x0len + i])) + atol;
	      diff = Enew[(k*x1len + j)*x0len + i] 
  		   - Eold[(k*x1len + j)*x0len + i];
	      tmp = fabs(diff/w);
	      Err_new_UV = max(Err_new_UV, tmp);
	    }
      }
    }

    // communicate to obtain overall error estimates
    float loc_err[] = {Err_new_Xr, Err_new_UV};
    float glob_err[2];
#ifdef USE_MPI
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg vars = 2;
    if (dtnorm > 0.0) 
      MPI_Allreduce(loc_err,glob_err,vars,DataType,MPI_SUM,MPI_COMM_WORLD);
    else
      MPI_Allreduce(loc_err,glob_err,vars,DataType,MPI_MAX,MPI_COMM_WORLD);
#else
    for (l=0; l<2; l++)  glob_err[l] = loc_err[l];
#endif

    // finish off error estimate norms
    if (dtnorm > 0.0) {
      Err_new_Xr = max(pow(glob_err[0], 1.0/dtnorm)/timeAccuracyXr/VolFac, Err_floor);
      Err_new_UV = max(pow(glob_err[1], 1.0/dtnorm)/timeAccuracyUV/VolFac, Err_floor);
    } else {
      Err_new_Xr = max(glob_err[0]/timeAccuracyXr, Err_floor);
      Err_new_UV = max(glob_err[1]/timeAccuracyUV, Err_floor);
    }

    // Set time step depending on how it has been set up by the user.
    //    dt_control determines the time adaptivity algorithm:
    //             0 -> I controller
    //             1 -> PI controller
    //             2 -> PID controller
    //          else -> original time controller
    if (dt_control == 0) {

      float k1 = -1.0/den;
      if (timeAccuracyXr != huge_number)
	dt_est_Xr = safety * dt * pow(Err_new_Xr,k1);
      if (timeAccuracyUV != huge_number)
	dt_est_UV = safety * dt * pow(Err_new_UV,k1);

      if (debug) {
	printf("  DualFLD_ComputeTimestep (I):\n");
	if (UseXray && (timeAccuracyXr != huge_number))
	  printf("     Xray Err = %8.2e, dt_est = %8.2e (%.1fx change)\n", 
		 Err_new_Xr, dt_est_Xr, dt_est_Xr/dt);
	if (UseUV && (timeAccuracyUV != huge_number))
	  printf("       UV Err = %8.2e, dt_est = %8.2e (%.1fx change)\n", 
		 Err_new_UV, dt_est_UV, dt_est_UV/dt);
      }

    } else if (dt_control == 1) {

      float k1 = -0.7/den;
      float k2 =  0.4/den;
      if (timeAccuracyXr != huge_number)
	dt_est_Xr = safety * dt * pow(Err_new_Xr,k1) * pow(Err_cur_Xr,k2);
      if (timeAccuracyUV != huge_number)
	dt_est_UV = safety * dt * pow(Err_new_UV,k1) * pow(Err_cur_UV,k2);
      if (debug) {
	printf("  DualFLD_ComputeTimestep (PI):\n");
	if (UseXray && (timeAccuracyXr != huge_number))
	  printf("     Xray Errs = %8.2e %8.2e, dt_est = %8.2e (%.1fx change)\n", 
		 Err_new_Xr, Err_cur_Xr, dt_est_Xr, dt_est_Xr/dt);
	if (UseUV && (timeAccuracyUV != huge_number))
	  printf("       UV Errs = %8.2e %8.2e, dt_est = %8.2e (%.1fx change)\n", 
		 Err_new_UV, Err_cur_UV, dt_est_UV, dt_est_UV/dt);
      }

    } else if (dt_control == 2) {

      float k1 = -0.49/den;
      float k2 =  0.34/den;
      float k3 = -0.1/den;
      if (timeAccuracyXr != huge_number)
	dt_est_Xr = safety * dt * pow(Err_new_Xr,k1) * pow(Err_cur_Xr,k2) * pow(Err_old_Xr,k3);
      if (timeAccuracyUV != huge_number)
	dt_est_UV = safety * dt * pow(Err_new_UV,k1) * pow(Err_cur_UV,k2) * pow(Err_old_UV,k3);
      if (debug) {
	printf("  DualFLD_ComputeTimestep (PID):\n");
	if (UseXray && (timeAccuracyXr != huge_number))
	  printf("     Xray Errs = %8.2e %8.2e %8.2e, dt_est = %8.2e (%.1fx change)\n", 
		 Err_new_Xr, Err_cur_Xr, Err_old_Xr, dt_est_Xr, dt_est_Xr/dt);
	if (UseUV && (timeAccuracyUV != huge_number))
	  printf("       UV Errs = %8.2e %8.2e %8.2e, dt_est = %8.2e (%.1fx change)\n", 
		 Err_new_UV, Err_cur_UV, Err_old_UV, dt_est_UV, dt_est_UV/dt);
      }

    } else {

      if (timeAccuracyXr != huge_number)
	dt_est_Xr = safety * dt / Err_new_Xr;
      if (timeAccuracyUV != huge_number)
	dt_est_UV = safety * dt / Err_new_UV;
      if (debug) {
	printf("  DualFLD_ComputeTimestep (orig):\n");
	if (UseXray && (timeAccuracyXr != huge_number))
	  printf("     Xray Err = %8.2e, dt_est = %8.2e (%.1fx change)\n", 
		 Err_new_Xr, dt_est_Xr, dt_est_Xr/dt);
	if (UseUV && (timeAccuracyUV != huge_number))
	  printf("       UV Err = %8.2e, dt_est = %8.2e (%.1fx change)\n", 
		 Err_new_UV, dt_est_UV, dt_est_UV/dt);
      }

    }

    // set dtrad as the minimum of the two estimates
    dt_est = min(dt_est_Xr, dt_est_UV);

  }

  // account for min/max time step size (according to user)
  dt_est = max(dt_est, mindt);
  dt_est = min(dt_est, maxdt);

  return dt_est;
}

#endif
