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
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE: This class defines problem-specific functions for a pair of 
/           split (UV & X-ray) implicit flux-limited diffusion solves.
/
************************************************************************/
#ifdef TRANSFER
#ifndef DUALFLD_DEFINED__
#define DUALFLD_DEFINED__

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "EnzoVector.h"
#include "ImplicitProblemABC.h"
#include "SED.h"


// increase the maximum number of FLD sources as needed
#define MAX_FLD_SOURCES       100


// DualFLD is hard-coded with four SED types:
//    -1 -> monochromatic at specified frequency
//     0 -> power law with specified slope
//     1 -> blackbody with temperature 1e5 K
//     2 -> PopII SED (like Ricotti etal, ApJ 575:33-48, 2002)
#define NUM_FLD_SED_TYPES  4

// DualFLD is hard-coded with three boundary condition types:
//     0 -> periodic
//     1 -> dirichlet
//     2 -> neumann
#define NUM_FLD_BDRY_TYPES  3

// DualFLD is hard-coded with three Krylov linear solver types:
//     0 -> PCG
//     1 -> BiCGStab
//     2 -> GMRES
#define NUM_FLD_SOL_TYPES  5

// DualFLD is hard-coded with four time step adaptivity algorithms:
//     -1 -> original time controller
//      0 -> I controller
//      1 -> PI controller
//      2 -> PID controller
#define NUM_FLD_DT_CONTROLLERS  4

// HYPRE is hard-coded with four relaxation algorithms:
//    0 -> Jacobi
//    1 -> weighted Jacobi
//    2 -> red-black Gauss-Seidel (symmetric)
//    3 -> red-black Gauss-Seidel (nonsymmetric)
#define NUM_HYPRE_RLX_TYPES  4



class DualFLD : public virtual ImplicitProblemABC {

 private:
  
  // overall time spent in solver and components
  float RTtime;
  float HYPREtime;
  
  // HYPRE Struct-specific data
  Eint32 mattype;                // HYPRE matrix type for solve
  Eint32 stSize;                 // stencil size
#ifdef USE_HYPRE
  HYPRE_StructGrid grid;         // HYPRE grid object for setup
  HYPRE_StructStencil stencil;   // stencil object
#endif

  // HYPRE Solver-specific data
  float  sol_tolerance_Xr;       // desired solver tolerance
  float  sol_tolerance_UV;
  Eint32 sol_MGmaxit_Xr;         // maximum number of MG iterations
  Eint32 sol_MGmaxit_UV;
  Eint32 sol_KryMaxit_Xr;        // maximum number of PCG iterations
  Eint32 sol_KryMaxit_UV;
  Eint32 sol_rlxtype_Xr;         // relaxation type:
                                 //    0,1 -> weighted Jacobi
                                 //    2,3 -> red-black Gauss-Seidel
  Eint32 sol_rlxtype_UV;
  Eint32 sol_npre_Xr;            // num. pre-relaxation sweeps
  Eint32 sol_npre_UV;
  Eint32 sol_npost_Xr;           // num. post-relaxation sweeps
  Eint32 sol_npost_UV;
  Eint32 sol_printl;             // print output level
  Eint32 sol_log;                // amount of logging
  int    sol_Krylov_Xr;          // type of outer solver to use for each field
  int    sol_Krylov_UV;
  Eint32 SolvIndices[3][2];      // L/R edge indices of subdomain in global mesh
                                 // Note: these INCLUDE Dirichlet zones, even 
                                 //   though those are not included as active 
                                 //   data in the vectors or physics routines.
  int SolvOff[3];                // offset between HYPRE mesh and active mesh; 
                                 //   typically 0, but will be 1 for inclusion 
                                 //   of Dirichlet zones in HYPRE grid.

  // HYPRE interface temporary data
#ifdef USE_HYPRE
  HYPRE_StructMatrix P;          // holds radiation matrix
  HYPRE_StructVector rhsvec;     // holds radiation rhs vector
  HYPRE_StructVector solvec;     // holds radiation solution vector
#endif
  Eflt64 *HYPREbuff;             // holds contiguous sections of rhs/sol

  // HYPRE solver diagnostics
  int totIters_Xr;               // total linear iterations for solves
  int totIters_UV;               // total linear iterations for solves

  // General problem grid information
  bool OnBdry[3][2]; // denotes if proc owns piece of boundary
  int rank;          // Rank of self-gravity problem
  int layout[3];     // number of procs in each dim (1-based)
  int location[3];   // location of this proc in each dim (0-based)
  int NBors[3][2];   // process IDs of L/R neighbors in each dim
  int LocDims[3];    // implicit problem local dims (no ghost or bdry cells)
  int ArrDims[3];    // local array sizes (includes ghost and bdry cells)
  int GhDims[3][2];  // ghost cells at each face (includes Dirichlet bdry zones)
  int GlobDims[3];   // implicit problem global dimensions (active cells only)
  float dx[3];             // mesh size in each dimension
  float EdgeVals[3][2];    // L/R edges of this proc's subdomain
  float *UVBdryVals[3][2];   // boundary values for UV BCs
  float *XrBdryVals[3][2];   // boundary values for Xray BCs

  // limiter parameters
  float FLD_Rmin;    // lower bound on R (normalized)
  float FLD_Dmax;    // upper bound on D (normalized)
  float FLD_Kmin;    // lower bound on opacity (normalized)
  float FLD_Emin;    // lower bound on radiation (normalized)

  // time-stepping related data
  float initdt;        // initial radiation time step size
  float maxdt;         // maximum radiation/chemistry/heating time step size
  float mindt;         // minimum radiation/chemistry/heating time step size
  float timeAccuracyXr; // desired relative change in Xray field per step
  float timeAccuracyUV; // desired relative change in UV field per step
  int   dt_control;    // time step controller algorithm:
  float Err_cur_Xr;    // storage for error estimates in adaptive time stepper
  float Err_cur_UV;
  float Err_new_Xr;
  float Err_new_UV;
  float dtnorm;        // norm choice for computing relative change:
                       //    0 -> max pointwise norm (default)
                       //   >0 -> rms p-norm over entire domain
  float dtgrowth;      // time step growth factor (1 < dtgrowth < 10)
  float tnew;          // new time
  float told;          // old time
  float dt;            // time step size
  float dtrad;         // radiation time step size (subcycled)
  float theta;         // implicitness parameter (1->BE, 0.5->CN, 0->FE)
  
  // problem defining data
  int   isothermal;    // flag to force an isothermal chemistry model
  bool  UseXray;       // flag to enable/disable Xray radiation
  bool  UseUV;         // flag to enable/disable UV radiation
  bool  XrStatic;      // flag to denote a static Xray radiation field
  bool  UVStatic;      // flag to denote a static UV radiation field

  // ionization source parameters
  int   NumSourcesXr;                                  // number of Xray input sources
  int   NumSourcesUV;                                  // number of UV input sources
  float SourceLocationXr[MAX_FLD_SOURCES][3];          // Xray source locations
  float SourceLocationUV[MAX_FLD_SOURCES][3];          // UV source locations
  float SourceEnergyXr[MAX_FLD_SOURCES];               // energy per Xray source [photons/sec]
  float SourceEnergyUV[MAX_FLD_SOURCES];               // energy per UV source [photons/sec]
  int   WeakScaling;                                   // whether to replicate inputs on each process
  float OriginalSourceLocationXr[MAX_FLD_SOURCES][3];  // necessary to restart WeakScaling runs
  float OriginalSourceLocationUV[MAX_FLD_SOURCES][3];

  // cosmology and scaling constants
  FLOAT a;             // cosmology expansion coefficient
  FLOAT a0;            // cosmology expansion coefficient (old time)
  FLOAT adot;          // time-derivative of a
  FLOAT adot0;         // time-derivative of a (old time)
  float aUnits;        // expansion parameter scaling
  bool  autoScale;     // flag to enable/disable automatic scaling factors
  bool  StartAutoScale;  // flag to turn begin automatic scaling in a run
  float UVScale;       // scaling factor for UV radiation energy density
  float UVUnits;       // UV radiation energy density unit conversion factor
  float UVUnits0;      // UV radiation energy density unit conversion factor
  float XrScale;       // scaling factor for X-ray radiation energy density
  float XrUnits;       // X-ray radiation energy density unit conversion factor
  float XrUnits0;      // X-ray radiation energy density unit conversion factor
  float NiUnits;       // species density unit conversion factor
  float NiUnits0;      // species density unit conversion factor

  float DenUnits;      // density scaling factor
  float LenUnits;      // length scaling factor
  float TimeUnits;     // time scaling factor
  float VelUnits;      // velocity scaling factor
  float DenUnits0;     // density scaling factor
  float LenUnits0;     // length scaling factor

  // storage for integrals over radiation spectrum (set during initialization)
  float hnu0_HI;       // HI ionization threshold (eV)
  float hnu0_HeI;      // HeI ionization threshold (eV)
  float hnu0_HeII;     // HeII ionization threshold (eV)
  int   UVSpectrum;    // UV radiation spectrum choice:
                       //   2 -> PopII SED (like Ricotti etal, ApJ 575:33–48, 2002)
                       //   1 -> 1e5 black body spectrum
                       //   0 -> simple power law spectrum (exponent -1.5)
                       //  -1 -> monochromatic spectrum at UVFrequency
  float UVFrequency;   // UV radiation frequency (if UVSpectrum = -1)
  int   XrSpectrum;    // X-ray spectrum choice (same values as above)
  float XrFrequency;   // X-ray radiation frequency (if XrSpectrum = -1)
  float RadIntUV[7];   // UV spectrum integrals:
                       //   0: int_{nu0}^{inf} chiUV(nu) dnu
                       //   1: int_{nu0}^{inf} chiUV(nu)*sigHI(nu) dnu
                       //   2: int_{nu0}^{inf} chiUV(nu)*sigHI(nu)/nu dnu
                       //   3: int_{nu1}^{inf} chiUV(nu)*sigHeI(nu) dnu
                       //   4: int_{nu1}^{inf} chiUV(nu)*sigHeI(nu)/nu dnu
                       //   5: int_{nu2}^{inf} chiUV(nu)*sigHeII(nu) dnu
                       //   6: int_{nu2}^{inf} chiUV(nu)*sigHeII(nu)/nu dnu
  float RadIntXr[7];   // X-ray spectrum integrals (same form as above)
  SED  *Xr_SED;
  SED  *UV_SED;

  // private solver storage
  EnzoVector *sol;     // solution vector

  // private computation routines
  int ReadParameters(TopGridData &MetaData, HierarchyEntry *ThisGrid);
  int ComputeOpacity(HierarchyEntry *ThisGrid, int XrUv, float *Opacity);
  int ComputeRadiationIntegrals();
  float ComputeTimeStep(EnzoVector *uold, EnzoVector *unew);
  int EnforceBoundary(HierarchyEntry *ThisGrid);
  int FillRates(HierarchyEntry *ThisGrid);
  float Limiter(float E1, float E2, float k1, float k2, float nUn, float lUn, float dxi);
  int RadiationSource(HierarchyEntry *ThisGrid, int XrUv, float *Eta);
  int SetupSystem(HierarchyEntry *ThisGrid, int XrUv, float *E,
		  float *Opacity, float *Eta, float &rhsnorm);
  int RadStep(HierarchyEntry *ThisGrid, int XrUV);


 public:

  // boundary type in each dimension, face:
  //    0->periodic
  //    1->dirichlet
  //    2->neumann
  int XrBdryType[3][2];
  int UVBdryType[3][2];

  ///////////////////////////////////////
  // FLD-Specific Routines

  // Constructor
  DualFLD();
  
  // Destructor
  ~DualFLD();

  // Problem Initializer
  int Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData);
  
  // Problem Evolver
  int Evolve(HierarchyEntry *ThisGrid, float deltat);
  
  // Write module parameters to file
  int WriteParameters(FILE *fptr);

  // Problem Boundary Condition setup (called once or at each time step, 
  //    must be called for each locally-owned external face separately)
  int SetupBoundary(int Dimension, int Face, int XuUv, 
		    int BdryConst, float *BdryData);

};


#endif   /* DUALFLD_DEFINED__ */
#endif   /* TRANSFER */
