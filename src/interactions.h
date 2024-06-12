/*
* Name : interaction.h
* Autor:
* Brief: Contain functions for the interaction calculations of the dipoles
*/

#ifndef INTERACTION_H_
#define INTERACTION_H_

// external libraries
#include <iostream>
#include <Eigen/Dense>
#include <limits>
#include <cmath>
#include <iomanip>

// internal libaries
#include "typedefs.h"
#include "mathFun.h"
#include "postprocessing.h"

using namespace std;
using namespace Eigen;

/////////////////////// init stuff for short range interactions //////////////////

void InitCutOffSR(Params_S* pParams, PartProps_S* pPartProp);

double FunCutOffError(double x, double c1, double delta, double r, bool derivation);

//////////////////// stuff interactions //////////////////////////////

void CompInteractVecs(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

void ApplyBoundaryCondVec(Array<double, 1, 3>* pVec, Array<bool, 1, 3>* pBufferVec, Params_S* pParams);


/////////// find optimal ewald parameters ///////////////

void CompEwaldParams(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

double FunRealError(double x, double M, double fAl, Params_S* pParams);

double FunRezError(double x, double M, double fAl, Params_S* pParams);

void InitEwaldConst(Params_S* pParams, Buffer_S* pBuffer);



/////////// perform ewald summation /////////////////////

void CompDipInteractEwaldR(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams, OutputVar_S* pOutputVar);

void CompDipInteractEwaldF(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

void EvalSinCos(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

void CompDipFieldEwaldR(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

void CompDipFieldEwaldF(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

void CheckEwaldSum(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams, OutputVar_S* pOutputVar);



/////////////////////// short range interactions//////////////////////////////////////////

void CompSterRep(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, int idx1, int idx2, int i);

void CompElectrostatRep(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, int idx1, int idx2, int i);

void CompForceVdW(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, int idx1, int idx2, int i);

#endif  /* INTERACTION_H_ */
