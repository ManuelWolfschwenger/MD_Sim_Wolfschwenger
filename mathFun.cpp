/*
* Name : mathFun.cpp
* Auto :
* Brief: Mathematical functions 
*/

// own header file
#include "mathFun.h"

// Brief: Rowwise cross product of arrays
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


// Brief: This function calculates the integral if a function with the simpsonregel (https://de.wikipedia.org/wiki/Simpsonregel)
// param[in] : double dMag...magnetic radius of particle
// return: double
double IntegralEffDiameter(double dMag)
{
	double area = 0;

	double lowBound = 0;
	double highBound = dMag + 2 * effLenSurfMol;

	int steps = 1000;
	double stepSize = (highBound - lowBound) / steps;

	for (int i = 0; i < steps; i++)
	{
		area += stepSize / 6 * (funEffDiameter(dMag, i * stepSize) + 4 * funEffDiameter(dMag, (i + 0.5) * stepSize) + funEffDiameter(dMag, (i + 1) * stepSize));
	}

	return area;
}


// Brief: calculates value for integratal to evaluate equivalent hard sphere diameter of rosensweig potential
// param[in] : double dMag...magnetic radius of particle
// param[in] : double r...distance between particles
// return: double
double funEffDiameter(double dMag, double r)
{
	double dDelta = dMag + 2 * effLenSurfMol;
	return 1 - exp(-my_pi * dMag * dMag * densSurfMol / (2 * effLenSurfMol) * (dDelta - r * (log(dDelta / (r + 1e-20)) + 1)));  //1e-20 to avoid division by zero
}
