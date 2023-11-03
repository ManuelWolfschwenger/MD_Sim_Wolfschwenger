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
#include <Eigen/CXX11/Tensor>
#include <limits>
#include <cmath>
#include <iomanip>

// internal libaries
#include "typedefs.h"
#include "constants.h"
#include "settings.h"
#include "mathFun.h"
#include "postprocessing.h"

using namespace std;
using namespace Eigen;

/////////////// neighbour list functions ////////////////

// Brief: This function initializes the parameters for the cell and neighbour list
// param[in]: Params_S* pParams - pointer to a Params_S object
// param[in]: PartProps_S* pPartProp - pointer to a PartProps_S object
// param[in/out]: NebrList_S* pLists - pointer to a NebrList_S object
// return: void
void InitNebrList(Params_S* pParams, NebrList_S* pLists, PartProps_S* pPartProp);

// Brief: calculation ofF_vdw/F_dip
// param[in]: double x dependent variable
// param[in]: bool derivation => true first derivation
// return: void
double FunCutOffError(double x, double c1, double delta, double r, bool derivation);

// Brief: This function builds the neighbour list
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Params_S* pParams - pointer to a Params_S object
// param[in/out]: NebrList_S* pLists - pointer to a NebrList_S object
// return: void
void BuildNebrList(WorkingVar_S* pWorkVar, Params_S* pParams, NebrList_S* pLists);

// Brief: cell wrapping, applies periodic boundary conditions to cells
// param[in] :Params_S* pParams - pointer to Params_S object
// param[in/out]: NebrList_S* pLists - pointer to a NebrList_S object
// param[in] :int col - index of vector column
// param[in/out] : pM2Vec - pointer to m2Vec
// param[in/out] : pShift - pointer to shifVec
// return: void
void CellWrapAround(Array<int, 1, 3>* pM2Vec, Array<double, 1, 3>* pShiftVec, Params_S* pParams, NebrList_S* pLists, int col);

// Brief: This function checks if the neighbour list has to be refreshed
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Params_S* pParams - pointer to a Params_S object
// param[in]: double deltaTused - used timestep
// return: void
void RefreshNebrList(WorkingVar_S* pWorkVar, Params_S* pParams, double deltaTused);


//////////////////////// interaction functions ///////////////////////////

// Brief: This function computes the interaction vector for one interaction pair and writes the data to an array
// // This function is only suitable for the all pairs method!!!
// This function is only suitable for the all pairs method!!!
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: Params_S* pParams - pointer to params_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompInteractVecsLR(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

// Brief: This function computes the interaction vector for one interaction pair and writes the data to an array 
// This function is only suitable for the neighbour list method!!!!!!!
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: Params_S* pParams - pointer to params_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompInteractVecNL(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams, int idx1, int idx2, int i);

// Brief: This function applies periodic boundary conditions
// param[in/out]: (Array<double, 1, 3>* pVec - pointer to a vector of coords to checked
// param[in/out]: (Array<bool, 1, 3>* pBoolVec - pointer to a vector of bools
// param[in]: Params_S* pParams - pointer to params_S object
// return: void
void ApplyBoundaryCondVec(Array<double, 1, 3>* pVec, Array<bool, 1, 3>* pBufferVec, Params_S* pParams);

/////////// find optimal ewald parameters ///////////////

// Brief: calculation of optimal Ewald summation parameters
// param[in]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
// param[in/out]: Params_S* pParams - pointer to Params_S object
// return: void
void CompEwaldParams(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

// Brief: calculation of rms error of real space forces
// param[in]:EwaldParams_S* pEwaldParams - pointer to a EwaldParams_S object
// param[in]: double x dependent variable
// param[in]: bool derivation => true first derivation
// return: void
double FunRealError(double x, double M, double fAl, Params_S* pParams);

// Brief: calculation of rms error of fourier space forces
// param[in]:EwaldParams_S* pEwaldParams - pointer to a EwaldParams_S object
// param[in]: double x dependent variable
// param[in]: bool derivation => true first derivation
// return: void
double FunRezError(double x, double M, double fAl, Params_S* pParams);

// Brief: precalculation of Constants for ewald summation
// param[in/out]: Params_S* pParams - pointer to Params_S object
// param[in/out]: Buffer_S* pBuffer - pointer to buffer_S object
// return: void
void InitEwaldConst(Params_S* pParams, Buffer_S* pBuffer);

/////////// perform ewald summation /////////////////////

// Brief: calculation of dipole interaction force in real space / must be called in all pairs loop
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompDipInteractEwaldR(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

// Brief: calculation of dipole interaction force in fourier space 
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompDipInteractEwaldF(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

// Brief: precalculation of sin and cos for fourier space
// param[in]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out]: Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: Params_S* pParams - pointer to Params_S object
// return: void
void EvalSinCos(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

// Brief: calculation of dipole interaction force in real space / must be called in all pairs loop
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: Params_S* pParams - pointer to Params_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompDipFieldEwaldR(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

// Brief: calculation of dipole interaction field in fourier space 
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompDipFieldEwaldF(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);

// Brief: this function calculates the exact force on the particles to check the ewald summation
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
// return: void
void CheckEwaldSum(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams);



/////////////////////// short range interactions//////////////////////////////////////////

// Brief: calculation of attractive and repulsive force (Rosensweig
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompSterRep(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, int idx1, int idx2, int i);

// Brief: calculation of electrostatic force 
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompElectrostatRep(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, int idx1, int idx2, int i);

// Brief: calculation of attractive van der Waals force
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompForceVdW(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, int idx1, int idx2, int i);

// Brief: check if particles overlap + correcting
// param[in/out] : WorkingVar_S* pWorkVar - pointer to workVar_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void HardSphereRep(WorkingVar_S* pWorkVar, int idx1, int idx2, int i);

// Brief: check if system is stable + correcting (instable if particle cores would overlap)
// param[in/out] : WorkingVar_S* pWorkVar - pointer to workVar_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CheckStability(WorkingVar_S* pWorkVar, int idx1, int idx2, int i);

#endif  /* INTERACTION_H_ */
