/***********************************************************************
/
/  INITIALIZE DUAL FLD TEST -- CLUMP IONIZATION TEST
/
/  written by: Daniel Reynolds
/  date:       November 2014
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "phys_constants.h"


/* default constants */
#define DEFAULT_MU 0.6       // mean molecular mass
#define MIN_TEMP 1.0         // minimum temperature [K]


// function prototypes
int InitializeRateData(FLOAT Time);


int DualRHIonizationClumpInitialize(FILE *fptr, FILE *Outfptr, 
				    HierarchyEntry &TopGrid,
				    TopGridData &MetaData, int local)
{
#ifdef TRANSFER
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("Entering DualRHIonizationClumpInitialize routine\n");

  char *kphHIName    = "HI_kph";
  char *kphHeIName   = "HeI_kph";
  char *kphHeIIName  = "HeII_kph";
  char *gammaName    = "PhotoGamma";
  char *kdissH2IName = "H2I_kdiss";
  char *DensName     = "Density";
  char *TEName       = "TotalEnergy";
  char *IEName       = "Internal_Energy";
  char *Vel0Name     = "x-velocity";
  char *Vel1Name     = "y-velocity";
  char *Vel2Name     = "z-velocity";
  char *RadName0     = "Xray_Radiation";
  char *RadName1     = "UV_Radiation";
  char *HIName       = "HI_Density";
  char *HIIName      = "HII_Density";
  char *HeIName      = "HeI_Density";
  char *HeIIName     = "HeII_Density";
  char *HeIIIName    = "HeIII_Density";
  char *DeName       = "Electron_Density";

  // local declarations
  char line[MAX_LINE_LENGTH];
  int  dim, ret;

  // Setup and parameters:
  float X0Velocity           = 0.0;
  float X1Velocity           = 0.0;
  float X2Velocity           = 0.0;
  float NumDensityIn         = 0.04;
  float NumDensityOut        = 0.0002;
  float TemperatureIn        = 40.0;     // K
  float TemperatureOut       = 8000.0;   // K
  float UVRadiation          = -1.0;
  float XrRadiation          = -1.0;
  float InitialFractionHII   = 0.0;
  float InitialFractionHeII  = 0.0;
  float InitialFractionHeIII = 0.0;
  float ClumpCenterX         = 5.0;  // normalized units
  float ClumpCenterY         = 3.3;
  float ClumpCenterZ         = 3.3;
  float ClumpRadius          = 0.8;
  int   use_xray             = 0;
  int   use_uv               = 0;

  // overwrite input from RadHydroParamFile file, if it exists
  if (MetaData.RadHydroParameterFname != NULL) {
    FILE *RHfptr;
    if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
      while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	ret = 0;
	// read relevant problem parameters
	ret += sscanf(line, "DualFLD_Velocity = %"FSYM" %"FSYM" %"FSYM,
		      &X0Velocity, &X1Velocity, &X2Velocity);
	ret += sscanf(line, "DualFLD_NumDensityIn = %"FSYM, &NumDensityIn);
	ret += sscanf(line, "DualFLD_NumDensityOut = %"FSYM, &NumDensityOut);
	ret += sscanf(line, "DualFLD_TemperatureIn = %"FSYM, &TemperatureIn);
	ret += sscanf(line, "DualFLD_TemperatureOut = %"FSYM, &TemperatureOut);
	ret += sscanf(line, "DualFLD_UVRadiationEnergy = %"FSYM, &UVRadiation);
	ret += sscanf(line, "DualFLD_XrayRadiationEnergy = %"FSYM, &XrRadiation);
	ret += sscanf(line, "DualFLD_InitialFractionHII = %"FSYM, &InitialFractionHII);
	if (!RadiativeTransferHydrogenOnly) {
	  ret += sscanf(line, "DualFLD_InitialFractionHeII = %"FSYM, &InitialFractionHeII);
	  ret += sscanf(line, "DualFLD_InitialFractionHeIII = %"FSYM, &InitialFractionHeIII);
	}
	ret += sscanf(line, "DualFLDUseXray = %"ISYM, &use_xray);
	ret += sscanf(line, "DualFLDUseUV = %"ISYM, &use_uv);
	
	ret += sscanf(line, "ClumpCenter = %"FSYM" %"FSYM" %"FSYM,
		      &ClumpCenterX, &ClumpCenterY, &ClumpCenterZ);
	ret += sscanf(line, "ClumpRadius = %"FSYM, &ClumpRadius);

      } // end input from parameter file
      fclose(RHfptr);
    }
  }

  /* error checking */
  if (use_xray && (XrRadiation < 0.0)) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: x-ray radiation enabled, but specified value %f is illegal; overriding input.\n", 
	      XrRadiation);
    XrRadiation = 1.e-30;
  }
  if (!use_xray && (XrRadiation >= 0.0)) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: x-ray radiation value specified, but field not engaged; disabling.\n");
    XrRadiation = -1.0;
  }
  if (use_uv && (UVRadiation < 0.0)) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: uv radiation enabled, but specified value %f is illegal; overriding input.\n", 
	      UVRadiation);
    UVRadiation = 1.e-30;
  }
  if (!use_uv && (UVRadiation >= 0.0)) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: uv radiation value specified, but field not engaged; disabling.\n");
    UVRadiation = -1.0;
  }

  // set up CoolData object if not already set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) {
      fprintf(stderr,"Error in InitializeRateData.\n");
      return FAIL;
    }

  // if temperature specified and not internal energy, perform conversion here
  TemperatureIn  = max(TemperatureIn, MIN_TEMP); // enforce minimum
  TemperatureOut = max(TemperatureOut,MIN_TEMP); // enforce minimum
  float nH, HI, HII, nHe, HeI, HeII, HeIII, ne, num_dens, mu;
  if (RadiativeTransferHydrogenOnly) {
    HI = 1.0 - InitialFractionHII;
    HII = InitialFractionHII;
    ne = HII;
    num_dens = HI + HII + ne;
    mu = 1.0/num_dens;
  } else {
    nH = CoolData.HydrogenFractionByMass;
    nHe = (1.0 - CoolData.HydrogenFractionByMass);
    HI = nH*(1.0 - InitialFractionHII);
    HII = nH*InitialFractionHII;
    HeII = nHe*InitialFractionHeII;
    HeIII = nHe*InitialFractionHeIII;
    HeI = nHe - HeII - HeIII;
    ne = HII + HeII/4.0 + HeIII/2.0;
    num_dens = 0.25*(HeI + HeII + HeIII) + HI + HII + ne;
    mu = 1.0/num_dens;
  }
  // compute the internal energy
  float IEnergyIn  = kboltz*TemperatureIn/mu/mh/(Gamma-1.0);
  float IEnergyOut = kboltz*TemperatureOut/mu/mh/(Gamma-1.0);
  
  // set up the grid(s) on this level
  HierarchyEntry *Temp = &TopGrid;
  while (Temp != NULL) {
    if (Temp->GridData->DualRHIonizationClumpInitializeGrid(
			NumDensityIn, NumDensityOut, 
			X0Velocity, X1Velocity, X2Velocity, 
			IEnergyIn, IEnergyOut, XrRadiation, UVRadiation, 
			InitialFractionHII, InitialFractionHeII, 
			InitialFractionHeIII, ClumpCenterX, ClumpCenterY, 
			ClumpCenterZ, ClumpRadius, local) == FAIL) {
      fprintf(stderr, "Error in DualRHIonizationClumpInitializeGrid.\n");
      return FAIL;
    }
    Temp = Temp->NextGridThisLevel;
  }


  // set up field names and units
  int BaryonField = 0;
  DataLabel[BaryonField++] = DensName;
  DataLabel[BaryonField++] = TEName;
  if (DualEnergyFormalism) 
    DataLabel[BaryonField++] = IEName;
  DataLabel[BaryonField++] = Vel0Name;
  DataLabel[BaryonField++] = Vel1Name;
  DataLabel[BaryonField++] = Vel2Name;
  if (XrRadiation >= 0.0)
    DataLabel[BaryonField++] = RadName0;
  if (UVRadiation >= 0.0)
    DataLabel[BaryonField++] = RadName1;
  DataLabel[BaryonField++] = DeName;
  DataLabel[BaryonField++] = HIName;
  DataLabel[BaryonField++] = HIIName;
  if (!RadiativeTransferHydrogenOnly) {
    DataLabel[BaryonField++] = HeIName;
    DataLabel[BaryonField++] = HeIIName;
    DataLabel[BaryonField++] = HeIIIName;
  }

  // if using external chemistry/cooling, set rate labels
  DataLabel[BaryonField++] = kphHIName;
  DataLabel[BaryonField++] = gammaName;
  if (!RadiativeTransferHydrogenOnly) {
    DataLabel[BaryonField++] = kphHeIName;
    DataLabel[BaryonField++] = kphHeIIName;
  }
  if (MultiSpecies > 1)
    DataLabel[BaryonField++] = kdissH2IName;

  for (int i=0; i<BaryonField; i++) 
    DataUnits[i] = NULL;

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;
 
#endif

}
