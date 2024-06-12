/*
* Name : mathFun.cpp
* Auto :
* Brief: Mathematical functions 
*/

// own header file
#include "mathFun.h"
#include "configuration.h"

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


// Brief: This function calculates the integral if a function with the simpson rule (https://de.wikipedia.org/wiki/Simpsonregel)
// param[in] : double dMag...magnetic radius of particle
// return: double
double IntegralEffDiameter(double dMag)
{
	double area = 0;

	double lowBound = 0;
	double highBound = dMag + 2 * config->getEffLenSurfMol();

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
// return: double effective diameter
double funEffDiameter(double dMag, double r)
{
	double dDelta = dMag + 2 * config->getEffLenSurfMol();
	return 1 - exp(-config->getMyPi() * dMag * dMag * config->getDensSurfMol() / (2 * config->getEffLenSurfMol()) * (dDelta - r * (log(dDelta / (r + 1e-20)) + 1)));  //1e-20 to avoid division by zero
}

// Brief: normalize vectors and project m onto the plane spanned by B and n / could be probably parallelized
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// return: void
void VectorProjection(WorkingVar_S* pWorkVar)
{
	RowWiseCrossProd(&pWorkVar->coordsRot.Bvec, &pWorkVar->coordsRot.posEa, &pWorkVar->coordsRot.vecNormal);
	pWorkVar->coordsRot.vecNormal.colwise() /= pWorkVar->coordsRot.vecNormal.rowwise().norm();

	pWorkVar->coordsRot.mProj = pWorkVar->coordsRot.posMm - pWorkVar->coordsRot.vecNormal.colwise()*(pWorkVar->coordsRot.posMm*pWorkVar->coordsRot.vecNormal).rowwise().sum();
	pWorkVar->coordsRot.posMm = pWorkVar->coordsRot.mProj.colwise() / pWorkVar->coordsRot.mProj.rowwise().norm();
}

// Brief: find local extrema of potential energy with Newton mehtod
// param[in]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: int i - index of particle
// return: void
 ArrayXXd FindLocalMinMax(WorkingVar_S* pWorkVar, int i)
 {
	double errTol = 1e-23; //error tolerance
	double x,f,df,Beff;

	vector<double>  minPsi, minEn, maxPsi, maxEn;

	Beff = sqrt((pWorkVar->demagFluxDens.row(i)*pWorkVar->demagFluxDens.row(i)).sum());
	pWorkVar->partProp.dPotEn = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*pWorkVar->partProp.sinPsi*pWorkVar->partProp.cosPsi +
			config->getSatMag()*pWorkVar->partProp.volMag(i)*Beff*(pWorkVar->partProp.psi - pWorkVar->coordsRot.phiIs(i)).sin();

	for (int j = 1; j < pWorkVar->partProp.dPotEn.rows(); j++)
	{
		//finding local minima
		if (pWorkVar->partProp.dPotEn(j-1) < 0 && pWorkVar->partProp.dPotEn(j) > 0) //=> Minimum crossed
		{
			x = pWorkVar->partProp.psi(j);

			//Newton iteration
			f = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*sin(x)*cos(x) +
					config->getSatMag()*pWorkVar->partProp.volMag(i)*Beff*sin(x - pWorkVar->coordsRot.phiIs(i));

			while (abs(f) > errTol)
			{
				df = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*cos(2*x) +
					config->getSatMag()*pWorkVar->partProp.volMag(i)*Beff*cos(x - pWorkVar->coordsRot.phiIs(i));
				
				x -= f/df;

				f = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*sin(x)*cos(x) +
					config->getSatMag()*pWorkVar->partProp.volMag(i)*Beff*sin(x - pWorkVar->coordsRot.phiIs(i));
			}	

			minPsi.push_back(x);
			minEn.push_back(config->getAnisEn()*pWorkVar->partProp.volMag(i)*pow(sin(x),2.0) - 
			config->getSatMag()*pWorkVar->partProp.volMag(i)*Beff*cos(x - pWorkVar->coordsRot.phiIs(i)));
		}	

		//finding local maxima
		if (pWorkVar->partProp.dPotEn(j-1) > 0 && pWorkVar->partProp.dPotEn(j) < 0) //=> Maximum crossed
		{
			x = pWorkVar->partProp.psi(j-1);

			//Newton iteration
			f = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*sin(x)*cos(x) +
					config->getSatMag()*pWorkVar->partProp.volMag(i)*Beff*sin(x - pWorkVar->coordsRot.phiIs(i));

			while (abs(f) > errTol)
			{
				df = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*cos(2*x) +
					config->getSatMag()*pWorkVar->partProp.volMag(i)*Beff*cos(x - pWorkVar->coordsRot.phiIs(i));
				
				x -= f/df;

				f = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*sin(x)*cos(x) +
					config->getSatMag()*pWorkVar->partProp.volMag(i)*Beff*sin(x - pWorkVar->coordsRot.phiIs(i));
			}	

			maxPsi.push_back(x);
			maxEn.push_back(config->getAnisEn()*pWorkVar->partProp.volMag(i)*pow(sin(x),2.0) - 
			config->getSatMag()*pWorkVar->partProp.volMag(i)*Beff*cos(x - pWorkVar->coordsRot.phiIs(i)));
		}			
	}
	
	// store data in array
	int len = minPsi.size();
	ArrayXXd extrema;
	extrema.resize(4,len);
	extrema.row(0) = Map<VectorXd, Unaligned>(minPsi.data(), minPsi.size());
	extrema.row(1) = Map<VectorXd, Unaligned>(minEn.data(), minEn.size());
	extrema.row(2) = Map<VectorXd, Unaligned>(maxPsi.data(), maxPsi.size());
	extrema.row(3) = Map<VectorXd, Unaligned>(maxEn.data(), maxEn.size());

	return extrema;
 }