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
#include "postprocessing.h"
#include "interactions.h"

// namespace definition
using namespace Eigen;
using namespace std;

int main(void);

////////////////////////*initialization functions*///////////////////////////////////////

void InitSimulation(void);

void InitParticleProps(PartProps_S* pPartProps);

void InitCoordsFromSketch(WorkingVar_S* pWorkVar, Params_S* pParams);

void InitCoordsFromSim(WorkingVar_S* pWorkVar, Params_S* pParams);

void InitIntCoeff(IntCoefficient_S* pIntCoeff);

void InitBuffers(Buffer_S* pBuffer);

void InitTwoStateApprox(WorkingVar_S* pWorkVar);



////////////////////////*simulation functions*///////////////////////////////////////

void RunSimulation(void);

void DefExtField(ArrayXXd* pExtFluxDens, double tInt);

void DefThermFluct(WorkingVar_S * pWorkVar, Buffer1d_S * pBuffer, double deltaT);

void ApplyBoundaryCondArr(ArrayXXd* pArr, Buffer3d_S* pBuffer3d, Params_S* pParams);

void CompInteractions(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams, OutputVar_S* pOutputVar);

void TwoStateApprox(WorkingVar_S* pWorkVar, double deltaT);

void IntegrationRot(WorkingVar_S* pWorkVar, IntCoefficient_S* pIntCoeff, Buffer_S* pBuffer, Params_S* pParams, double* pTDelta, double* pTInt);

void RKCoeffs(WorkingVar_S* pWorkVar, CoordsRot_S* pCoord, RkCoeff_S* pRkCoeffs, Buffer_S* pBuffer);

void AdjustTimeStep(WorkingVar_S* pWorkVar, double error, double power, double* pTDelta, bool* pCond, bool* pNoFailed);

void IntegrationTrans(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, double deltaTused, Params_S* pParams);

#endif  /* SIMULATUION_H_ */

