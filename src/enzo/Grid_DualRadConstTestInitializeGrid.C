/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE RAD-HYDRO IONIZATION TEST) 
/
/  written by: Daniel Reynolds
/  date:       July 2007
/
/  PURPOSE: 
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"


// function prototypes
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);


int grid::DualRadConstTestInitializeGrid(float DensityConstant, 
					 float VxConstant, 
					 float VyConstant, 
					 float VzConstant, 
					 float IEConstant, 
					 float XrConstant, 
					 float UVConstant, 
					 float InitialFractionHII, 
					 float InitialFractionHeII, 
					 float InitialFractionHeIII, 
					 int   local) {

#ifdef TRANSFER
//   if (debug)
//     fprintf(stdout,"Entering grid::DualRadConstTestInitializeGrid routine\n");

  // determine whether data should be allocated/initialized or not
  int NewData = TRUE;
  if ((ParallelRootGridIO == TRUE) && (local == 0))
    NewData = FALSE;

  // if grids allocated and already set up (i.e. restart), return
  if ((NumberOfBaryonFields > 5) && (BaryonField[5] != NULL))
    return SUCCESS;

  // create necessary baryon fields
  int RhoNum, TENum, IENum, V0Num, V1Num, V2Num, XrNum, UVNum, DeNum, 
    HINum, HIINum, HeINum, HeIINum, HeIIINum, kphHINum, kphHeINum, 
    kphHeIINum, gammaNum, kdissH2INum;
  NumberOfBaryonFields = 0;
  FieldType[RhoNum = NumberOfBaryonFields++] = Density;
  FieldType[TENum = NumberOfBaryonFields++]  = TotalEnergy;
  if (DualEnergyFormalism) 
    FieldType[IENum = NumberOfBaryonFields++]  = InternalEnergy;
  FieldType[V0Num = NumberOfBaryonFields++]    = Velocity1;
  FieldType[V1Num = NumberOfBaryonFields++]    = Velocity2;
  FieldType[V2Num = NumberOfBaryonFields++]    = Velocity3;
  if (XrConstant >= 0.0)
    FieldType[XrNum = NumberOfBaryonFields++]  = RadiationFreq0;
  if (UVConstant >= 0.0)
    FieldType[UVNum = NumberOfBaryonFields++]  = RadiationFreq1;
  FieldType[DeNum = NumberOfBaryonFields++]    = ElectronDensity;
  FieldType[HINum = NumberOfBaryonFields++]    = HIDensity;
  FieldType[HIINum = NumberOfBaryonFields++]   = HIIDensity;
  if (!RadiativeTransferHydrogenOnly) {
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;    
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
  }
  // set external chemistry/cooling rate fields
  FieldType[kphHINum = NumberOfBaryonFields++] = kphHI;
  FieldType[gammaNum = NumberOfBaryonFields++] = PhotoGamma;
  if (!RadiativeTransferHydrogenOnly) {
    FieldType[kphHeINum  = NumberOfBaryonFields++] = kphHeI;
    FieldType[kphHeIINum = NumberOfBaryonFields++] = kphHeII;
  }
  if (MultiSpecies > 1)
    FieldType[kdissH2INum = NumberOfBaryonFields++] = kdissH2I;


  // set the subgrid static flag (necessary??)
  SubgridsAreStatic = FALSE;

  // Return if this doesn't concern us.
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  // Get various units
  double MassUnits=1;
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr,"Error in GetUnits.\n");
    return FAIL;
  }
  if (debug && NewData) {
    fprintf(stdout,"  Internal Unit Conversion Factors:\n");
    fprintf(stdout,"         length = %g\n",LengthUnits);
    fprintf(stdout,"           mass = %lg\n",MassUnits);
    fprintf(stdout,"           time = %g\n",TimeUnits);
  }

  // compute size of fields
  int dim, i;
  int size = 1;
  for (dim=0; dim<GridRank; dim++)  size *= GridDimension[dim];
 
  // allocate fields
  if (NewData == TRUE) {
    for (int field=0; field<NumberOfBaryonFields; field++)
      if (BaryonField[field] == NULL)
	BaryonField[field] = new float[size];
      for (i=0; i<size; i++)
	BaryonField[RhoNum][i] = DensityConstant/DensityUnits;
    
    // set fluid density, total energy, [internal energy,] velocities, 
    // radiation energy, electron density, chemical species
    float TEConstant = (IEConstant + 0.5*(VxConstant*VxConstant + 
					  VyConstant*VyConstant + 
					  VzConstant*VzConstant));
    float HIIConstant = InitialFractionHII*CoolData.HydrogenFractionByMass*DensityConstant;
    float HIConstant = CoolData.HydrogenFractionByMass*DensityConstant - HIIConstant;
    float HeIIConstant = InitialFractionHeII*DensityConstant*(1.0-CoolData.HydrogenFractionByMass);
    float HeIIIConstant = InitialFractionHeIII*DensityConstant*(1.0-CoolData.HydrogenFractionByMass);
    float HeIConstant = (1.0-CoolData.HydrogenFractionByMass)*DensityConstant 
      - HeIIConstant - HeIIIConstant;
    float DeConstant = HIIConstant + 0.25*HeIIConstant + 0.5*HeIIIConstant;
    float eUnits = VelocityUnits*VelocityUnits;
    float RUnits = DensityUnits*eUnits;

    for (i=0; i<size; i++) {
      BaryonField[RhoNum][i]   = DensityConstant/DensityUnits;
      BaryonField[TENum][i]    = TEConstant/eUnits;
      BaryonField[V0Num][i]    = VxConstant/VelocityUnits;
      BaryonField[V1Num][i]    = VyConstant/VelocityUnits;
      BaryonField[V2Num][i]    = VzConstant/VelocityUnits;
      if (XrConstant >= 0.0)
	BaryonField[XrNum][i]  = XrConstant/RUnits;
      if (UVConstant >= 0.0)
	BaryonField[UVNum][i]  = UVConstant/RUnits;
      BaryonField[DeNum][i]    = DeConstant/DensityUnits;
      BaryonField[HINum][i]    = HIConstant/DensityUnits;
      BaryonField[HIINum][i]   = HIIConstant/DensityUnits;
      if (!RadiativeTransferHydrogenOnly) {
	BaryonField[HeINum][i]   = HeIConstant/DensityUnits;
	BaryonField[HeIINum][i]  = HeIIConstant/DensityUnits;
	BaryonField[HeIIINum][i] = HeIIIConstant/DensityUnits;
      }
    }
    if (DualEnergyFormalism)
      for (i=0; i<size; i++)
	BaryonField[IENum][i] = IEConstant/eUnits;

    // set external chemistry/cooling rate fields
    for (i=0; i<size; i++)  BaryonField[kphHINum][i] = 0.0;
    for (i=0; i<size; i++)  BaryonField[gammaNum][i] = 0.0;
    if (!RadiativeTransferHydrogenOnly) {
      for (i=0; i<size; i++)  BaryonField[kphHeINum][i]  = 0.0;
      for (i=0; i<size; i++)  BaryonField[kphHeIINum][i] = 0.0;
    }
    if (MultiSpecies > 1)
      for (i=0; i<size; i++)  BaryonField[kdissH2INum][i] = 0.0;


    // output some information on the test problem
    if (debug && NewData) {
      printf("\n  Initializing constant fields using CGS values:\n");
      printf("        density = %g\n",DensityConstant);
      printf("   total energy = %g\n",TEConstant);
      if (DualEnergyFormalism)
	printf("   internal energy = %g\n",IEConstant);
      printf("     x-velocity = %g\n",VxConstant);
      printf("     y-velocity = %g\n",VyConstant);
      printf("     z-velocity = %g\n",VzConstant);
      if (XrConstant >= 0.0)
	printf(" Xray radiation = %g\n",XrConstant);
      if (UVConstant >= 0.0)
	printf("   UV radiation = %g\n",UVConstant);
      printf("      electrons = %g\n",DeConstant);
      printf("            nHI = %g\n",HIConstant);
      printf("           nHII = %g\n",HIIConstant);
      if (!RadiativeTransferHydrogenOnly) {
	printf("           nHeI = %g\n",HeIConstant);
	printf("          nHeII = %g\n",HeIIConstant);
	printf("         nHeIII = %g\n",HeIIIConstant);
      }
      
      printf("\n  Corresponding Enzo internal values:\n");
      printf("        density = %g\n",BaryonField[RhoNum][1]);
      printf("   total energy = %g\n",BaryonField[TENum][1]);
      if (DualEnergyFormalism)
	printf("   internal energy = %g\n",BaryonField[IENum][1]);
      printf("     x-velocity = %g\n",BaryonField[V0Num][1]);
      printf("     y-velocity = %g\n",BaryonField[V1Num][1]);
      printf("     z-velocity = %g\n",BaryonField[V2Num][1]);
      if (XrConstant >= 0.0)
	printf(" Xray radiation = %g\n",BaryonField[XrNum][1]);
      if (UVConstant >= 0.0)
	printf("   UV radiation = %g\n",BaryonField[UVNum][1]);
      printf("      electrons = %g\n",BaryonField[DeNum][1]);
      printf("            nHI = %g\n",BaryonField[HINum][1]);
      printf("           nHII = %g\n",BaryonField[HIINum][1]);
      if (!RadiativeTransferHydrogenOnly) {
	printf("           nHeI = %g\n",BaryonField[HeINum][1]);
	printf("          nHeII = %g\n",BaryonField[HeIINum][1]);
	printf("         nHeIII = %g\n",BaryonField[HeIIINum][1]);
      }
    }

  } // end if NewData == TRUE

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;

#endif

}
