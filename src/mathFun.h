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
#include "typedefs.h"

using namespace Eigen;


void RowWiseCrossProd(ArrayXXd* pA, ArrayXXd* pB, ArrayXXd* pRes);

double IntegralEffDiameter(double dMag);

double funEffDiameter(double dMag, double r);

void VectorProjection(WorkingVar_S* pWorkVar);

ArrayXXd FindLocalMinMax(WorkingVar_S* pWorkVar, int i);

#endif  /* MATHFUN_H_ */
