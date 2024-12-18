/*
* Name : mathFun.h
* Autor:
* Brief: Contains the simualtion headers
*/

#ifndef MATHFUN_H_
#define MATHFUN_H_

// external libaries
#include <Eigen/Dense>

void RowWiseCrossProd(Eigen::ArrayXXd* pA, Eigen::ArrayXXd* pB, Eigen::ArrayXXd* pRes);

#endif  /* MATHFUN_H_ */