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

// internal headers files
#include "mathFun.h"
#include "postprocessing.h"
#include "typedefs.h"
#include "configuration.h"

// namespace definition
using namespace Eigen;
using namespace std;

// Brief: Runs the complete simulation
// Return: program status

int main(void);

////////////////////////*initialization functions*///////////////////////////////////////

void InitSimulation(void);

void InitParticleProps(PartProps_S* pPartProps);

void InitCoords(Coord_S* pCoords);

void InitIntCoeff(IntCoefficient_S* pIntCoeff);

void InitBuffers(Buffer3d_S* pBuffer3d, Buffer1d_S* pBuffer1d);



////////////////////////*simulation functions*///////////////////////////////////////

void RunSimulation(void);

void DefExtField(ArrayXXd* pExtFluxDens, double tInt);

void DefThermFluct(WorkingVar_S* pWorkVar, double tDelta);

void Integration(WorkingVar_S* pWorkVar, IntCoefficient_S* pIntCoeff, double* pTDelta, double* pTInt);

void RKCoeffs(WorkingVar_S* pWorkVar, Coord_S* pCoord, RkCoeff_S* pRkCoeffs);

void AdjustTimeStep(WorkingVar_S* pWorkVar, double error, double power, double* pTDelta, bool* pCond, bool* pNoFailed);

#endif  /* SIMULATUION_H_ */
