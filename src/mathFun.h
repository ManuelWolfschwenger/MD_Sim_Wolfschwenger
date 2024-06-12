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

using namespace Eigen;


void RowWiseCrossProd(ArrayXXd* pA, ArrayXXd* pB, ArrayXXd* pRes);

double IntegralEffDiameter(double dMag);

double funEffDiameter(double dMag, double r);

#endif  /* MATHFUN_H_ */
