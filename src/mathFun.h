/*
* Name : mathFun.h
* Autor:
* Brief: Contains the simualtion headers
*/

#ifndef MATHFUN_H_
#define MATHFUN_H_

// external libaries
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <bits/stdc++.h>

using namespace Eigen;
using namespace std;

void RowWiseCrossProd(ArrayXXd* pA, ArrayXXd* pB, ArrayXXd* pRes);

#endif  /* MATHFUN_H_ */