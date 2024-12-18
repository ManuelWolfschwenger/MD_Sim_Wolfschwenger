/*
* Name : mathFun.cpp
* Auto :
* Brief: Simulation source file
*/

// own header file
#include "mathFun.h"

// Brief: Rowwise cross product
// param[in] : ArrayXXd* pA - pointer to a ArrayXXd object
// param[in] : ArrayXXd* pB - pointer to a ArrayXXd object
// param[out]: ArrayXXd* pRes - pointer to a ArrayXXd object
// return: void
void RowWiseCrossProd(ArrayXXd* pA, ArrayXXd* pB, ArrayXXd* pRes)
{
	pRes->col(0) = pA->col(1) * pB->col(2) - pA->col(2) * pB->col(1);
	pRes->col(1) = pA->col(2) * pB->col(0) - pA->col(0) * pB->col(2);
	pRes->col(2) = pA->col(0) * pB->col(1) - pA->col(1) * pB->col(0);
}
