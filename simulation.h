/*
* Name : simulation.h
* Auto :
* Brief: Contains the simualtion headers
*/

#ifndef SIMULATUION_H_
#define SIMULATUION_H_

// external libaries
#include <iostream>
#include <Eigen/Dense>
#include <EigenRand/EigenRand>
#include <time.h>
#include <cmath>
#include <fstream>

// internal headers files
#include "typedefs.h"
#include "mathFun.h"
#include "settings.h"
#include "constants.h"
#include "postprocessing.h"
#include "interactions.h"

// namespace definition
using namespace Eigen;
using namespace std;

// random number generator
Rand::Vmt19937_64 generator;

/* function definitions */

// Brief: Runs the complete simulation
// Return: program status
int main(void);

////////////////////////*initialization functions*///////////////////////////////////////

// Brief: This function initializes the simulation
// return: void
void InitSimulation(void);

// Brief: This function sets the particle properties
// param[in/out]: PartProps_S* pPartProps - pointer to a PartProps_S object
// return: void
void InitParticleProps(PartProps_S* pPartProps);

// Brief: This function sets the initial particle positions
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Params_S* pParams - pointer to params_S object
// return: void
void InitCoordsFromSketch(WorkingVar_S* pWorkVar, Params_S* pParams);

// Brief: This function sets the initial particle positions from last simulation step
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Params_S* pParams - pointer to params_S object
// return: void
void InitCoordsFromSim(WorkingVar_S* pWorkVar, Params_S* pParams);

// Brief: This function initializes the integration coefficant buffers
// param[in/out]: IntCoefficant_S* pIntCoeff - pointer to a IntCoefficant_S object
// return: void
void InitIntCoeff(IntCoefficant_S* pIntCoeff);

// Brief: This function initializes the buffers
// param[in/out]: Buffer_S* pBuffer - pointer to a Buffer_S object
// return: void
void InitBuffers(Buffer_S* pBuffer);

////////////////////////*simulation functions*///////////////////////////////////////

// Brief: This function performs the simulation
// return: void
void RunSimulation(void);

// Brief: This function sets the external field depending on time in z direction
// param[in/out]: ArrayXXd* pExtFluxDens - pointer to the external field object
// param[in]	: double tInt -  current integration time in s
// return: void
void DefExtField(ArrayXXd* pExtFluxDens, double tInt);

// Brief: This function creates random thermal fluctuations
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]	: Buffer1d_S* pBuffer1d - pointer to Buffer1d_S object
// param[in]	: double actual timestep
// return: void
void DefThermFluct(WorkingVar_S * pWorkVar, Buffer1d_S * pBuffer, double deltaT);

// Brief: This function applies periodic boundary conditions
// param[in/out]: ArrayXXd* pArr - pointer to an array of coords to check
// param[in]	: Buffer3d_S* pBuffer3d - pointer to Buffer3d_S object
// param[in]: Params_S* pParams - pointer to params_S object
// return: void
void ApplyBoundaryCondArr(ArrayXXd* pArr, Buffer3d_S* pBuffer3d, Params_S* pParams);

// Brief: This function computes the forces acting on a particle
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]	: Buffer_S* pBuffer - pointer to Buffer_S object
// param[in]: Params_S* pParams - pointer to params_S object
// return: void
void CompInteractions(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

// Brief: This function computes the demagnetizing field if the interactions are enabled
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]	: Buffer_S* pBuffer - pointer to Buffer_S object
// return: void
void CompDemagField(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams, CoordsRot_S* pCoordsRot);

// Brief: This function performs the numerical integration with the selected integration method
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out]: IntCoefficant_S* pIntCoeff - pointer to a IntCoefficant_S object
// param[in/out]: Buffer_S* pBuffer - pointer to Buffer_S* object 
// param[in]	: double* pTDelta - timestep for integration (adjusted by "AdjustTimestep")
// param[out]	: double* pTDeltaUsed - pointer to the actual used timestep by the integration routine
// return: void
void IntegrationRot(WorkingVar_S* pWorkVar, IntCoefficant_S* pIntCoeff, Buffer_S* pBuffer, Params_S* pParams, double* pTDelta, double* pTInt);

// Brief: Calculate the runge-kutta coefficients for magnetic moment and easy axes movement
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out]: Coord_S* pCoord - pointer to the current Coord_S object
// param[in/out]: RkCoeff_S* pRkCoeffs - pointer to the current RK coefficants
// param[in]: Buffer_S* pBuffer - pointer to Buffer_S* object 
// return: void
void RKCoeffs(WorkingVar_S* pWorkVar, CoordsRot_S* pCoord, RkCoeff_S* pRkCoeffs, Buffer_S* pBuffer);

// Brief: Adjusting timeStep to meet tolerance requirements 
// param[in] : WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : double error - error value
// param[in] : double power - power for timestep adjustment
// param[out]: double* pTDelta - new time delate
// param[out]: bool* pCond - new loop condidtion
// param[out]: bool* pNoFailed - new failed condidtion
// return: void
void AdjustTimeStep(WorkingVar_S* pWorkVar, double error, double power, double* pTDelta, bool* pCond, bool* pNoFailed);

// Brief: Integration of translational movement
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: deltaT timestep used in the rotational integration
// return: void
void IntegrationTrans(WorkingVar_S* pWorkVar, double deltaT);

#endif  /* SIMULATUION_H_ */

