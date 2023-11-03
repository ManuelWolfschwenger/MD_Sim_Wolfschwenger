/*
* Name : mathFun.h
* Autor:
* Brief: Contains the simualtion headers
*/

#ifndef MATHFUN_H_
#define MATHFUN_H_

// external libaries
#include <Eigen/Dense>

//internal header files
#include "constants.h"

using namespace Eigen;

// Brief: Rowwise cross product of arrays
// param[in] : ArrayXXd* pA - pointer to a ArrayXXd object
// param[in] : ArrayXXd* pB - pointer to a ArrayXXd object
// param[out]: ArrayXXd* pRes - pointer to a ArrayXXd object
// return: void
void RowWiseCrossProd(ArrayXXd* pA, ArrayXXd* pB, ArrayXXd* pRes);

// Brief: This function calculates the integral if a function with the simpsonregel (https://de.wikipedia.org/wiki/Simpsonregel)
// param[in] : double dMag...magnetic radius of particle
// return: double
double IntegralEffDiameter(double dMag);

// Brief: calculates value for integratal to evaluate equivalent hard sphere diameter of rosensweig potential
// param[in] : double dMag...magnetic radius of particle
// param[in] : double r...distance between particles
// return: double
double funEffDiameter(double dMag, double r);

#endif  /* MATHFUN_H_ */
