/*
* Name : interaction.cpp
* Auto :
* Brief: Contain functions for the interaction calculations of the particles
*/

// own header file
#include "interactions.h"
#include "configuration.h"

//////////////////////// interaction functions ///////////////////////////

/////////////////////// init stuff for short range interactions //////////////////

// Brief: This function initializes the parameters for the cell and neighbour list
// param[in]: Params_S* pParams - pointer to a Params_S object
// param[in]: PartProps_S* pPartProp - pointer to a PartProps_S object
// param[in/out]: NebrList_S* pLists - pointer to a NebrList_S object
// return: void
void InitCutOffSR(Params_S* pParams, PartProps_S* pPartProp)
{
	double mMax;       //maximum magnetic moment
	double rMax;       //maximum magnetic radius
	double c1; // constant for F_vdw/F_dip
	double errNewton; // error for finding root of function from above

	mMax = pPartProp->magMom.maxCoeff(); // biggest particle
	rMax = pPartProp->rMag.maxCoeff();

	if (config->getEnableCoatingPot())
	{
		// find value where F_vdw/F_dip = errTol to cut off short range forces
		c1 = 128 * config->getHamaker() * pow(rMax, 6.0) * config->getMyPi() / (18 * config->getMu0() * mMax * mMax);
		errNewton = 100;
		pParams->rCutSR = 1.1 * 2 * rMax; //start value for iteration

		// Newton iteration for finding root of function
		do
		{
			pParams->rCutSR -= FunCutOffError(pParams->rCutSR, c1, config->getErrTolSR(), rMax, false) / FunCutOffError(pParams->rCutSR, c1, config->getErrTolSR(), rMax, true);
			errNewton = abs(FunCutOffError(pParams->rCutSR, c1, config->getErrTolSR(), rMax, false));

		} while (errNewton > 1e-7);

		////////////// now we know the cut off radius to make a desired error of errtolSR in Fvdw/Fdip/////////////////

		if (pParams->rCutSR < 2 * pPartProp->rHydr.maxCoeff())
		{
			cout << "error particles can only interact if they overlap => edit cut off radius" << endl;
			exit(-1);
		}
	}
	else // if the short range interactions are treated with the hard sphere potential
	{
		pParams->rCutSR = 2 * rMax;
	}
}

// Brief: calculation of F_vdw/F_dip
// param[in]: double x dependent variable
// param[in]: bool derivation => true first derivation
// return: void
double FunCutOffError(double x, double c1, double delta, double r, bool derivation)
{
	double val;

	if (derivation)
	{
		val = c1 * (4 * r * r + 3 * x * x) / pow(4 * r * r - x * x, 3.0);
	}
	else
	{
		val = c1 * x / pow(x * x - 4 * r * r, 2.0) - delta;
	}

	return val;
}

//////////////////// stuff for interactions //////////////////////////////

// Brief: This function computes the interaction vector for all interaction pairs and writes the data to an array (vectors could be possible also used for short range interactions)
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: Params_S* pParams - pointer to params_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompInteractVecs(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams)
{
	double r, dCont, d;
	int pairs = 0;

	// have to go through all pairs because rcCutLR > L/3
	for (int idx1 = 0; idx1 < config->getNumbPart() - 1; idx1++)
	{
		for (int idx2 = idx1 + 1; idx2 < config->getNumbPart(); idx2++)
		{
			//vector pointing from dipole i to dipole j
			pBuffer->bufferVec.vec_r = pWorkVar->coordsTrans.pos.row(idx2) - pWorkVar->coordsTrans.pos.row(idx1);

			// evaluate shortest distance (across boundaries or not)
			ApplyBoundaryCondVec(&pBuffer->bufferVec.vec_r, &pBuffer->bufferVec.vec_b, pParams);

			r = pBuffer->bufferVec.vec_r.matrix().norm(); //distance between centerpoints of particle i and j

			if (r < pParams->rCutLR)
			{
				//store vectorlength and vector coordinates in interaction arrays for vectors
				pWorkVar->intLists.intVecsLen.row(pairs) = r;
				pWorkVar->intLists.intVecsCoords.row(pairs) = pBuffer->bufferVec.vec_r / r;
				pWorkVar->intLists.intPairs(2 * pairs) = idx1;
				pWorkVar->intLists.intPairs(2 * pairs + 1) = idx2;
				//everything is stored in lists because the vectors are used in every partial step of the solver

				if (r < pParams->rCutSR)
				{
					pWorkVar->intLists.intRangeSR(pairs) = true;
				}

				dCont = pWorkVar->partProp.rMag(idx1) + pWorkVar->partProp.rMag(idx2) + 1e-10; //1e-10 to avoid singularity in vdw attraction
				// too correct unphysical overlaps / especially for hard sphere potential
				if (r < dCont)
				{
					d = dCont - r; //overlapping distance

					pWorkVar->coordsTrans.pos.row(idx1) -= 0.5 * d * pWorkVar->intLists.intVecsCoords.row(pairs); //move back half the overlapping distance
					pWorkVar->coordsTrans.pos.row(idx2) += 0.5 * d * pWorkVar->intLists.intVecsCoords.row(pairs);
					pWorkVar->intLists.intVecsLen.row(pairs) = dCont;
				}

				pairs++;
			}
		}
	}

	pWorkVar->intLists.pairs = pairs;
}


// Brief: This function applies periodic boundary conditions for the interaction vector
// param[in/out]: (Array<double, 1, 3>* pVec - pointer to a vector of coords to checked
// param[in/out]: (Array<bool, 1, 3>* pBoolVec - pointer to a vector of bools
// param[in]: Params_S* pParams - pointer to params_S object
// return: void
void ApplyBoundaryCondVec(Array<double, 1, 3>* pVec, Array<bool, 1, 3>* pBoolVec, Params_S* pParams)
{
	//use boolean bufferVec as idx array
	//if particle position is greater than lbox/2 => moove to other side of box
	*pBoolVec = *pVec > pParams->lBox / 2;
	*pVec -= pParams->lBox * (*pBoolVec).cast<double>(); //multiply arrays elementwise, cast to convert bool to double

	//if particle position is smaller than -lbox/2 => moove to other side of box
	*pBoolVec = *pVec < -pParams->lBox / 2;
	*pVec += pParams->lBox * (*pBoolVec).cast<double>();
}

/////////// find optimal ewald parameters ///////////////

// Brief: calculation of optimal Ewald summation parameters
// param[in]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
// param[in/out]: Params_S* pParams - pointer to Params_S object
// return: void
void CompEwaldParams(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams)
{
	double maxVal, t, tMean;
	double M, fAl, tAl;
	double nCutLRtest;
	int idx1, idx2, minIdx, maxIdx, tIter;
	const int steps = 10000;
	const int rcSteps = 30;
	tIter = 10; //change it to small values for debugging!!!!

	ArrayXXd funVals, paramArr; // should be stored on heap
	funVals.resize(steps, 2);
	paramArr.resize(rcSteps, 4); // rc, time, alpha, kc

	// constants from particle properties
	M = (pWorkVar->partProp.magMom.pow(2.0)).sum();
	fAl = 3 * config->getMu0() * pow(pWorkVar->partProp.magMom.minCoeff(), 2.0) / (2.0 * config->getMyPi() * pow(2.0 * pWorkVar->partProp.rHydr.minCoeff(), 4.0));  //force between aligned dipoles
	tAl = pow(pWorkVar->partProp.magMom.minCoeff(), 2.0) * config->getMu0() / (16.0 * config->getMyPi() * pow(pWorkVar->partProp.rHydr.minCoeff(), 3.0)); //torque on dip1 in field from dip2 when in contract and magMoments 90degrees to each other

	paramArr.col(0).setLinSpaced(rcSteps, 0.2, 0.5); //variation of rcut to try

	for (int j = 0; j < rcSteps; j++)
	{
		// set rCut for real space
		pParams->rCutLR = pParams->lBox * paramArr(j,0);

		////////////// finding alpha for specific errTol ///////////////////
		
		funVals.col(0).setLinSpaced(steps, 0, 2e8);

		idx1 = 0;
		idx2 = 0;
		for (int i = 0; i < steps; i++)
		{
			funVals(i, 1) = FunRealError(funVals(i, 0), M, fAl, pParams);  // fill array with function values

			if (funVals(i, 1) < 0)
			{
				idx2 = i;
				idx1 = i - 1;
				break;
			}
		}

		if (idx1 == 0 || idx2 == 0)
		{
			cout << "error: alpha for error tolerance not found!" << endl;
		}

		paramArr(j, 2) = -funVals(idx1, 1) / (funVals(idx2, 1) - funVals(idx1, 1)) * (funVals(idx2, 0) - funVals(idx1, 0)) + funVals(idx1, 0); // interpolate to get root
		pParams->alpha = paramArr(j, 2);
		// ===>>>> now we know alpha for a given error tolerance!!!!!

		////////////// finding kc for specific errTol ///////////////////
		// 
		funVals.col(0).setLinSpaced(steps, 0, 30);

		for (int i = 0; i < steps; i++)
		{
			funVals(i, 1) = FunRezError(funVals(i, 0), M, fAl, pParams);  // fill array with function values
		}

		funVals.col(1).maxCoeff(&maxIdx);
		maxVal = funVals.col(1).maxCoeff();

		// check if max value is greater than 0
		if (maxVal < 0)
		{
			cout << "attention: maximum value of error < 0 in kc computation => can't find a kc value to reach the desired error" << endl;
			cout << "the reziprocal error made now is " << abs(maxVal) / (config->getErrTolEwald() / sqrt(2.0)) * 100.0 << "% smaller than the desired one" <<  endl;
			cout << "rc was L*" << paramArr(j, 0) << "\n" << endl;
			exit(-1);
		}

		// start at max val => decaying curve
		for (int i = maxIdx; i < steps; i++)
		{
			if (funVals(i, 1) < 0)
			{
				idx2 = i;
				idx1 = i - 1;
				break;
			}
		}

		nCutLRtest = -funVals(idx1, 1) / (funVals(idx2, 1) - funVals(idx1, 1)) * (funVals(idx2, 0) - funVals(idx1, 0)) + funVals(idx1, 0);

		if (nCutLRtest < 2.0) // important because otherwise function EvalSinCos doesn't work anymore
		{
			nCutLRtest = 2.0; // if it would be smaller and we increase the cut off => smaller error => no problem
		}
		
		paramArr(j, 3) = nCutLRtest;
		pParams->nCutLR = paramArr(j, 3);
		
		// ===>>>> now we know nCutLR for a given error tolerance!!!!!

		//////////////////////// measuring simulation time of ewald summation ////////////////////////////////////

		pWorkVar->intForce.setZero();         //force array, initialized to zero
		pWorkVar->demagFluxDens.setZero();    //demagnetizing array, initialized to zero
		pWorkVar->partProp.magMomVec = pWorkVar->coordsRot.posMm.colwise() * pWorkVar->partProp.magMom; //create magnetic moment vectors
		tMean = 0;

		// evaluate interaction vectors
		CompInteractVecs(pWorkVar, pBuffer, pParams);

		// init constants for Ewald summation
		InitEwaldConst(pParams, pBuffer);
		EvalSinCos(pWorkVar, pBuffer, pParams);

		// measuring time and calculate mean of it
		for (int k = 0; k < tIter; k++)
		{
			t = clock(); // measure time for real space calculation
		
			CompDipFieldEwaldR(pWorkVar, pBuffer, pParams);
			CompDipFieldEwaldF(pWorkVar, pBuffer, pParams);

			t = (clock() - t) / CLOCKS_PER_SEC;
			tMean += t;
		}
		tMean /= tIter; // mean simulation time 
		paramArr(j, 1) = tMean;

		cout << j + 1 << ", ";		
	}

	paramArr.col(1).minCoeff(&minIdx);
	pParams->rCutLR = paramArr(minIdx, 0) * pParams->lBox;
	pParams->alpha = paramArr(minIdx, 2);
	pParams->nCutLR = paramArr(minIdx, 3);

	cout << "\n" << endl;
	cout << std::setprecision(5);
	//cout << "        rc        time        alpha          kc" << endl;
	//cout << paramArr << endl;
	//cout << "\n" << endl;
	cout << "the optimal paramters for the  ewald-summation are: rc = L*" << paramArr(minIdx, 0) << " alpha*L = " << pParams->alpha * pParams->lBox << " kc = " << pParams->nCutLR << endl;

	pParams->nCutLR = floor(paramArr(minIdx, 3) + 0.5); // round up kc for grid vector
	cout << "kc is rounded to " << pParams->nCutLR << endl;
	
	if (pParams->alpha * pParams->rCutLR < 2.0)
	{
		cout << "attention: error estimates are made under the assumption: alpha*rc >> 1, but alpha*rc = " << pParams->alpha * pParams->rCutLR << endl;
	}

	// calculation of errors in real and reziprocal space
	double cc = 4.0 * pow(pParams->alpha * pParams->rCutLR, 4.0) +  6.0 * pow(pParams->alpha * pParams->rCutLR, 2.0) + 3.0;
	double dc = 8.0 * pow(pParams->alpha * pParams->rCutLR, 6.0) + 20.0 * pow(pParams->alpha * pParams->rCutLR, 4.0) + 30.0 * pow(pParams->alpha * pParams->rCutLR, 2.0) + 15.0;
	double bc = 2.0 * pow(pParams->alpha * pParams->rCutLR, 2.0) + 1.0;

	double errRealForce = config->getFacEwald() * M / sqrt(pow(pParams->lBox, 3.0) * pow(pParams->alpha, 4.0) * pow(pParams->rCutLR, 9.0) * config->getNumbPart()) * sqrt(13.0 / 6.0 * cc * cc + 2.0 / 15.0 * dc * dc - 13.0 / 15.0 * cc * dc) * exp(-pow(pParams->alpha * pParams->rCutLR, 2.0));
	double errRezForce = config->getFacEwald() * 8.0 * config->getMyPi() * M / pow(pParams->lBox, 3.0) * pParams->alpha * sqrt(2.0 * config->getMyPi() * pow(pParams->nCutLR, 3.0) / (15.0 * config->getNumbPart())) * exp(-pow(config->getMyPi() * pParams->nCutLR / (pParams->alpha * pParams->lBox), 2.0));
	double errRealTorque = config->getFacEwald() * M / sqrt(pow(pParams->lBox, 3.0) * pow(pParams->alpha, 4.0) * pow(pParams->rCutLR, 7.0) * config->getNumbPart()) * sqrt(0.5 * bc * bc + 0.2 * cc * cc) * exp(-pow(pParams->alpha * pParams->rCutLR, 2.0));
	double errRezTorque = config->getFacEwald() * 4.0 * M / pow(pParams->lBox, 2.0) * pParams->alpha * sqrt(config->getMyPi() * pParams->nCutLR / (5.0 * config->getNumbPart())) * exp(-pow(config->getMyPi() * pParams->nCutLR / (pParams->alpha * pParams->lBox), 2.0));
	
	cout << "the resulting force error should be: " << sqrt(errRezForce * errRezForce + errRealForce * errRealForce) / fAl << "*F_contact" << endl;
	cout << "the resulting torque error should be: " << sqrt(errRezTorque * errRezTorque + errRealTorque * errRealTorque) / tAl << "*T_contact" << "\n" << endl;
}

// Brief: calculation of rms error of real space forces
// param[in]: double x - dependent variable
// param[in]: double M - sum of squared magnetic moments
// param[in]: double fAl - force between aligned dipoles in contact to relate to
// param[in]: Params_S* pParams - pointer to Params_S object
// return: double realError - error for real space
double FunRealError(double x, double M, double fAl, Params_S* pParams)
{
	double a = x * pParams->rCutLR;
	double b = M / sqrt(pow(pParams->lBox, 3.0) * pow(pParams->rCutLR, 9.0) * config->getNumbPart() * pow(x, 4.0));
	double c = 4.0 * pow(a, 4.0) + 6.0 * pow(a, 2.0) + 3.0;
	double d = 8.0 * pow(a, 6.0) + 20.0 * pow(a, 4.0) + 30.0 * pow(a, 2.0) + 15.0;

	return config->getFacEwald() * b * sqrt(13.0 / 6.0 * c * c + 2.0 / 15.0 * d * d - 13.0 / 15.0 * c * d) * exp(-a * a) / fAl - config->getErrTolEwald() / sqrt(2.0);
}

// Brief: calculation of rms error of real space forces
// param[in]: double x - dependent variable
// param[in]: double M - sum of squared magnetic moments
// param[in]: double fAl - force between aligned dipoles in contact to relate to
// param[in]: Params_S* pParams - pointer to Params_S object
// return: double realError - error for real space
double FunRezError(double x, double M, double fAl, Params_S* pParams)
{
	return  config->getFacEwald() * 8.0 * M * config->getMyPi() / pow(pParams->lBox, 3.0) * pParams->alpha * sqrt(2.0 * config->getMyPi() * x * x * x / (15.0 * config->getNumbPart())) * exp(-pow(config->getMyPi() * x / (pParams->alpha * pParams->lBox), 2.0)) / fAl - config->getErrTolEwald() / sqrt(2.0);
}

// Brief: precalculation of Constants for ewald summation
// param[in/out]: Params_S* pParams - pointer to Params_S object
// param[in/out]: Buffer_S* pBuffer - pointer to buffer_S object
// return: void
void InitEwaldConst(Params_S* pParams, Buffer_S* pBuffer)
{
	// initialize size of Ewald arrays in buffer structure
	pBuffer->tSin.resize(pParams->nCutLR + 1, config->getNumbPart());
	pBuffer->tCos.resize(pParams->nCutLR + 1, config->getNumbPart());

	// constants for fourier space
	pParams->nc = (int)pParams->nCutLR;
	pParams->nc2 = pParams->nCutLR * pParams->nCutLR;
	pParams->gu = 2 * config->getMyPi() * (pParams->lBox * pParams->lBox * pParams->lBox);
	pParams->gr = 4 * config->getMyPi() * pParams->gu / pParams->lBox;
	pParams->gs = 2 * pParams->gu;
	pParams->w = pow(config->getMyPi() / (pParams->lBox * pParams->alpha), 2.0);
	pParams->cFluxDens = 4.0 * config->getMyPi() * config->getFacEwald() / pow(pParams->lBox, 3.0);
	pParams->cForce = 8.0 * config->getMyPi() * config->getMyPi() * config->getFacEwald() / pow(pParams->lBox, 4.0);
	pParams->cSurf = config->getFacEwald() * 4.0 * config->getMyPi() / (3 * pow(pParams->lBox, 3.0));   //(2*epsilon +1) = 3 for box surrounded by vacuum 
	pParams->cSelf = config->getFacEwald() * 4.0 * pow(pParams->alpha, 3.0) / (3.0 * sqrt(config->getMyPi()));

	// constants for eal space
	pParams-> alpha2 = pParams->alpha * pParams->alpha;
	pParams->alpha4 = pParams->alpha2 * pParams->alpha2;
}

/////////// perform ewald summation /////////////////////

// Brief: calculation of dipole interaction force in real space and evaluation of pressure tensor components
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to Buffer_S object
// param[in]: Params_S* pParams - pointer to Params_S object
// return: void
void CompDipInteractEwaldR(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams, OutputVar_S* pOutputVar)
{
	int idx1, idx2;
	double r, r2, t, a1, a2, a3, dot1, dot2;
	double alpha = pParams->alpha;
	double alpha2 = pParams->alpha2;

	// loop for long range interaction / real part of ewald sum
	for (int i = 0; i < pWorkVar->intLists.pairs; i++)
	{
		idx1 = pWorkVar->intLists.intPairs(2 * i);
		idx2 = pWorkVar->intLists.intPairs(2 * i + 1);

		r = pWorkVar->intLists.intVecsLen(i);
		r2 = r * r;
		t = 2.0 * alpha * exp(-alpha2 * r2) * config->getIrPi() / r2;
		a1 = erfc(alpha * r) / (r2 * r) + t;
		a2 = 3.0 * a1 / r2 + 2.0 * alpha2 * t;
		a3 = 5.0 * a2 / r2 + 4.0 * pParams->alpha4 * t;

		dot1 = pWorkVar->partProp.magMomVec.row(idx1).matrix().dot(pWorkVar->intLists.intVecsCoords.row(i).matrix());
		dot2 = pWorkVar->partProp.magMomVec.row(idx2).matrix().dot(pWorkVar->intLists.intVecsCoords.row(i).matrix());

		// force calculation
		pBuffer->bufferVec.vec1 = config->getFacEwald() * ((a2 * pWorkVar->partProp.magMomVec.row(idx1).matrix().dot(pWorkVar->partProp.magMomVec.row(idx2).matrix()) -
			a3 * dot1 * dot2 * r2) * pWorkVar->intLists.intVecsCoords.row(i) * r +
			a2 * r * (dot2 * pWorkVar->partProp.magMomVec.row(idx1) + dot1 * pWorkVar->partProp.magMomVec.row(idx2)));

		// flux density calculation
		pWorkVar->demagFluxDens.row(idx1) += config->getFacEwald() * (a2 * dot2 * pWorkVar->intLists.intVecsCoords.row(i) * r2 - a1 * pWorkVar->partProp.magMomVec.row(idx2)); //field on particle i;
		pWorkVar->demagFluxDens.row(idx2) += config->getFacEwald() * (a2 * dot1 * pWorkVar->intLists.intVecsCoords.row(i) * r2 - a1 * pWorkVar->partProp.magMomVec.row(idx1)); //field on particle j;

		// add to forces
		pWorkVar->intForce.row(idx1) -= pBuffer->bufferVec.vec1;
		pWorkVar->intForce.row(idx2) += pBuffer->bufferVec.vec1;

		// calculate pressure tensor for viscosity
		EvalPressTensLR(pWorkVar, pOutputVar, pBuffer, i);
	}
}

// Brief: calculation of dipole interaction force in fourier space 
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to Buffer_S object
// param[in]: Params_S* pParams - pointer to Params_S object
// return: void
void CompDipInteractEwaldF(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams)
{
	double e, lenSqr, sumC, sumS, preFac;
	double pc, ps;
	int nc = pParams->nc;

	Vector3d vc, vs;

	// loops over upper part of cube
	for (int nz = 0; nz <= nc; nz++)
	{
		for (int ny = -nc; ny <= nc; ny++)
		{
			for (int nx = -nc; nx <= nc; nx++)
			{
				pBuffer->bufferVec.vec1 << nx, ny, nz;
				lenSqr = pBuffer->bufferVec.vec1.matrix().dot(pBuffer->bufferVec.vec1.matrix());

				if (lenSqr == 0 || lenSqr > pParams->nc2)
				{
					continue;
				}

				e = 2.0 * exp(-pParams->w * lenSqr) / lenSqr;
				
				if (nz == 0)
				{
					e *= 0.5; //symmetric cube of cells in z direction
				}

				// sums are the same for all particles, but vary with cells
				sumC = 0;
				sumS = 0;
		
				for (int i = 0; i < config->getNumbPart(); i++)
				{
					vc << pBuffer->tCos(abs(nx), i).x(), pBuffer->tCos(abs(ny), i).y(), pBuffer->tCos(abs(nz), i).z(); //vc
					vs << pBuffer->tSin(abs(nx), i).x(), pBuffer->tSin(abs(ny), i).y(), pBuffer->tSin(abs(nz), i).z(); //vs

					if (nx < 0)
					{
						vs.x() = -vs.x();
					}
					if (ny < 0)
					{
						vs.y() = -vs.y();
					}

					pc = vc.x() * vc.y() * vc.z() - vc.x() * vs.y() * vs.z() - vs.x() * vc.y() * vs.z() - vs.x() * vs.y() * vc.z();
					ps = vs.x() * vc.y() * vc.z() + vc.x() * vs.y() * vc.z() + vc.x() * vc.y() * vs.z() - vs.x() * vs.y() * vs.z();

					preFac = pBuffer->bufferVec.vec1.matrix().dot(pWorkVar->partProp.magMomVec.row(i).matrix());
					sumC += preFac * pc;
					sumS += preFac * ps;
				}

				for (int i = 0; i < config->getNumbPart(); i++)
				{
					vc << pBuffer->tCos(abs(nx), i).x(), pBuffer->tCos(abs(ny), i).y(), pBuffer->tCos(abs(nz), i).z(); //vc
					vs << pBuffer->tSin(abs(nx), i).x(), pBuffer->tSin(abs(ny), i).y(), pBuffer->tSin(abs(nz), i).z(); //vs

					if (nx < 0)
					{
						vs.x() = -vs.x();
					}
					if (ny < 0)
					{
						vs.y() = -vs.y();
					}

					pc = vc.x() * vc.y() * vc.z() - vc.x() * vs.y() * vs.z() - vs.x() * vc.y() * vs.z() - vs.x() * vs.y() * vc.z();
					ps = vs.x() * vc.y() * vc.z() + vc.x() * vs.y() * vc.z() + vc.x() * vc.y() * vs.z() - vs.x() * vs.y() * vs.z();
					
					pWorkVar->intForce.row(i) += pParams->cForce * e * (pBuffer->bufferVec.vec1.matrix().dot(pWorkVar->partProp.magMomVec.row(i).matrix())) * (sumC * ps - sumS * pc) * pBuffer->bufferVec.vec1;
					pWorkVar->demagFluxDens.row(i) -= pParams->cFluxDens * e * (sumC * pc + sumS * ps) * pBuffer->bufferVec.vec1;
				}			
			}
		}
	}

	//surface term and self term,
	for (int i = 0; i < config->getNumbPart(); i++)
	{
		pWorkVar->demagFluxDens.row(i) += pParams->cSelf * pWorkVar->partProp.magMomVec.row(i) - pParams->cSurf * pWorkVar->partProp.magMomSum;
	}
}


// Brief: precalculation of sin and cos for fourier space => faster calculation (from Rapaport "The Art of molecular dynamics simulation")
// param[in]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out]: Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: Params_S* pParams - pointer to Params_S object
// return: void
void EvalSinCos(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams)
{
	Array<double,1,3> t, tt, u, w;

	t.setConstant(2 * config->getMyPi() / pParams->lBox); //t...buffervec 1

	for (int i = 0; i < config->getNumbPart(); i++)
	{
		tt = t * pWorkVar->coordsTrans.pos.row(i); //tt...buffervec 2
		pBuffer->tCos(0, i).setConstant(1);
		pBuffer->tSin(0, i).setZero();
		pBuffer->tCos(1, i) = tt.cos();
		pBuffer->tSin(1, i) = tt.sin();
		u = 2.0 * pBuffer->tCos(1, i); //u...buffervec 3
		pBuffer->tCos(2, i) = u * pBuffer->tCos(1, i); //only have 2 rows!!!!
		pBuffer->tSin(2, i) = u * pBuffer->tSin(1, i);
		tt.setConstant(1);
		pBuffer->tCos(2, i) -= tt;

		for (int j = 3; j <= pParams->nCutLR; j++)
		{
			w = u * pBuffer->tCos(j - 1, i); //w...buffervec_r
			pBuffer->tCos(j, i) = w - pBuffer->tCos(j - 2, i);
			w = u * pBuffer->tSin(j - 1, i);
			pBuffer->tSin(j, i) = w - pBuffer->tSin(j - 2, i);
		}
	}
}


// Brief: calculation of dipole interaction field in real space / for partial solver steps
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: Params_S* pParams - pointer to Params_S object
// return: void
void CompDipFieldEwaldR(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams)
{
	int idx1, idx2;
	double r, r2, t, a1, a2, a3, dot1, dot2;
	double alpha = pParams->alpha;
	double alpha2 = pParams->alpha2;

	// loop for long range interaction / real part of ewald sum
	for (int i = 0; i < pWorkVar->intLists.pairs; i++)
	{
		idx1 = pWorkVar->intLists.intPairs(2 * i);
		idx2 = pWorkVar->intLists.intPairs(2 * i + 1);

		r = pWorkVar->intLists.intVecsLen(i);
		r2 = r * r;
		t = 2.0 * alpha * exp(-alpha2 * r2) * config->getIrPi() / r2;
		a1 = erfc(alpha * r) / (r2 * r) + t;
		a2 = 3.0 * a1 / r2 + 2.0 * alpha2 * t;
		a3 = 5.0 * a2 / r2 + 4.0 * pParams->alpha4 * t;

		dot1 = pWorkVar->partProp.magMomVec.row(idx1).matrix().dot(pWorkVar->intLists.intVecsCoords.row(i).matrix());
		dot2 = pWorkVar->partProp.magMomVec.row(idx2).matrix().dot(pWorkVar->intLists.intVecsCoords.row(i).matrix());

		// flux density calculation
		pWorkVar->demagFluxDens.row(idx1) += config->getFacEwald() * (a2 * dot2 * pWorkVar->intLists.intVecsCoords.row(i) * r2 - a1 * pWorkVar->partProp.magMomVec.row(idx2)); //field on particle i;
		pWorkVar->demagFluxDens.row(idx2) += config->getFacEwald() * (a2 * dot1 * pWorkVar->intLists.intVecsCoords.row(i) * r2 - a1 * pWorkVar->partProp.magMomVec.row(idx1)); //field on particle j;
	}
}

// Brief: calculation of dipole interaction field in fourier space / for partial solver steps
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
// return: void
void CompDipFieldEwaldF(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams)
{
	double e, lenSqr, sumC, sumS, preFac;
	double pc, ps;
	int nc = pParams->nc;

	Vector3d vc, vs;

	// loops over upper part of cube
	for (int nz = 0; nz <= nc; nz++)
	{
		for (int ny = -nc; ny <= nc; ny++)
		{
			for (int nx = -nc; nx <= nc; nx++)
			{
				pBuffer->bufferVec.vec1 << nx, ny, nz;
				lenSqr = pBuffer->bufferVec.vec1.matrix().dot(pBuffer->bufferVec.vec1.matrix());

				if (lenSqr == 0 || lenSqr > pParams->nc2)
				{
					continue;
				}

				e = 2.0 * exp(-pParams->w * lenSqr) / lenSqr;

				if (nz == 0)
				{
					e *= 0.5; //symmetric cube of cells in z direction
				}

				// sums are the same for all particles, but vary with cells
				sumC = 0;
				sumS = 0;

				for (int i = 0; i < config->getNumbPart(); i++)
				{
					vc << pBuffer->tCos(abs(nx), i).x(), pBuffer->tCos(abs(ny), i).y(), pBuffer->tCos(abs(nz), i).z(); //vc
					vs << pBuffer->tSin(abs(nx), i).x(), pBuffer->tSin(abs(ny), i).y(), pBuffer->tSin(abs(nz), i).z(); //vs

					if (nx < 0)
					{
						vs.x() = -vs.x();
					}
					if (ny < 0)
					{
						vs.y() = -vs.y();
					}

					pc = vc.x() * vc.y() * vc.z() - vc.x() * vs.y() * vs.z() - vs.x() * vc.y() * vs.z() - vs.x() * vs.y() * vc.z();
					ps = vs.x() * vc.y() * vc.z() + vc.x() * vs.y() * vc.z() + vc.x() * vc.y() * vs.z() - vs.x() * vs.y() * vs.z();

					preFac = pBuffer->bufferVec.vec1.matrix().dot(pWorkVar->partProp.magMomVec.row(i).matrix());
					sumC += preFac * pc;
					sumS += preFac * ps;
				}

				for (int i = 0; i < config->getNumbPart(); i++)
				{
					vc << pBuffer->tCos(abs(nx), i).x(), pBuffer->tCos(abs(ny), i).y(), pBuffer->tCos(abs(nz), i).z(); //vc
					vs << pBuffer->tSin(abs(nx), i).x(), pBuffer->tSin(abs(ny), i).y(), pBuffer->tSin(abs(nz), i).z(); //vs

					if (nx < 0)
					{
						vs.x() = -vs.x();
					}
					if (ny < 0)
					{
						vs.y() = -vs.y();
					}

					pc = vc.x() * vc.y() * vc.z() - vc.x() * vs.y() * vs.z() - vs.x() * vc.y() * vs.z() - vs.x() * vs.y() * vc.z();
					ps = vs.x() * vc.y() * vc.z() + vc.x() * vs.y() * vc.z() + vc.x() * vc.y() * vs.z() - vs.x() * vs.y() * vs.z();

					pWorkVar->demagFluxDens.row(i) -= pParams->cFluxDens * e * (sumC * pc + sumS * ps) * pBuffer->bufferVec.vec1;
				}
			}
		}
	}

	//surface term and self term,
	for (int i = 0; i < config->getNumbPart(); i++)
	{
		pWorkVar->demagFluxDens.row(i) += pParams->cSelf * pWorkVar->partProp.magMomVec.row(i) - pParams->cSurf * pWorkVar->partProp.magMomSum;
	}
}

// Brief: this function calculates the exact force on the particles to check the ewald summation
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
// return: void
void CheckEwaldSum(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams, OutputVar_S* pOutputVar)
{
    ArrayXXd  exactForces, ewaldForces, exactFields, ewaldFields, exactTorques, ewaldTorques, diff;
	double errIs, fAl, tAl, bAl;

	exactForces.resize(config->getNumbPart(), 3);
	ewaldForces.resize(config->getNumbPart(), 3);
	exactFields.resize(config->getNumbPart(), 3);
	ewaldFields.resize(config->getNumbPart(), 3);
	exactTorques.resize(config->getNumbPart(), 3);
	ewaldTorques.resize(config->getNumbPart(), 3);
	diff.resize(config->getNumbPart(), 3);

	exactForces.setZero();
	ewaldForces.setZero();
	exactFields.setZero();
	ewaldFields.setZero();
	exactTorques.setZero();
	ewaldTorques.setZero();
	diff.setZero();

	// store solution in arrays
	ewaldForces = pWorkVar->intForce - pWorkVar->thermForce;
	ewaldFields = pWorkVar->demagFluxDens;
	RowWiseCrossProd(&pWorkVar->partProp.magMomVec, &ewaldFields, &ewaldTorques);

	//////////////////////////////// compute exact solution //////////////////////////////////////////////
	pParams->alpha = 16.0 / pParams->lBox;
	pParams->nCutLR = 40;

	// set workvar to zero
	pWorkVar->intForce.setZero();
	pWorkVar->demagFluxDens.setZero();

	InitEwaldConst(pParams, pBuffer); //comp constants depending on new parameter
	cout << "attention: restart simulation after checking Ewald sum and disable it, because arrays then have the wrong size" << "\n" << endl;
	EvalSinCos(pWorkVar, pBuffer, pParams);

	//Ewald Summation
	CompDipInteractEwaldR(pWorkVar, pBuffer, pParams, pOutputVar);
	CompDipInteractEwaldF(pWorkVar, pBuffer, pParams);

	exactForces = pWorkVar->intForce;
	exactFields = pWorkVar->demagFluxDens;
	RowWiseCrossProd(&pWorkVar->partProp.magMomVec, &exactFields, &exactTorques);

	cout << std::setprecision(20);
	cout << "forces on particle 1 in N" << endl;
	cout << exactForces.row(1) << " exact" << endl;
	cout << ewaldForces.row(1) << " Ewald" << endl;

	cout << "\n" << endl;
	cout << "magnetic fields on particle 1 in T" << endl;
	cout << exactFields.row(1) << endl;
	cout << ewaldFields.row(1) << endl;

	cout << "\n" << endl;
	cout << "torques on particle 1 in Nm" << endl;
	cout << exactTorques.row(1) << endl;
	cout << ewaldTorques.row(1) << "\n" << endl;

	//////////////////////////// compare to ewald summation /////////////////////////////////////////

	diff = exactForces - ewaldForces;
	errIs = sqrt(((diff * diff).rowwise().sum()).sum() / config->getNumbPart()); // rms error made by ewald summation
	fAl = 3 * config->getMu0() * pow(pWorkVar->partProp.magMom.minCoeff(), 2.0) / (2 * config->getMyPi() * pow(2 * pWorkVar->partProp.rHydr.minCoeff(), 4.0));  //force between aligned dipoles in contact
	cout << "the actual rms force error is " << errIs / fAl<< "*F_contact" << endl;

	diff = exactFields - ewaldFields;
	errIs = sqrt(((diff * diff).rowwise().sum()).sum() / config->getNumbPart()); // rms error made by ewald summation
	bAl =  config->getMu0() * pWorkVar->partProp.magMom.minCoeff() / (2 * config->getMyPi() * pow(2 * pWorkVar->partProp.rHydr.minCoeff(), 3.0));  //field of dipole i at dipole j when in contact
	cout << "the actual rms field error is " << errIs / bAl << "*B_contact" << endl;

	diff = exactTorques - ewaldTorques;
	errIs = sqrt(((diff * diff).rowwise().sum()).sum() / config->getNumbPart()); // rms error made by ewald summation
	tAl = pow(pWorkVar->partProp.magMom.minCoeff(), 2.0) * config->getMu0() / (16.0 * config->getMyPi() * pow(pWorkVar->partProp.rHydr.minCoeff(), 3.0));
	cout << "the actual rms torque error is " << errIs / tAl << "*T_contact" << endl;
}


/////////////////////// short range interactions //////////////////////////////////////////

// Brief: calculation of repulsive force (Rosensweig) + adding pair interaction force to buffervec2 for pressure tensor calculation
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompSterRep(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, int idx1, int idx2, int i)
{
	double dDelta = 2* (pWorkVar->partProp.rMag(idx1) + config->getEffLenSurfMol());
	double r = pWorkVar->intLists.intVecsLen(i);

	// repulsive potential
	if (r <= dDelta)
	{
		pBuffer->bufferVec.vec1 = (2 * pWorkVar->partProp.rMag(idx2) / (pWorkVar->partProp.rMag(idx1) + pWorkVar->partProp.rMag(idx2))) * //correction for spheres of differnt size
			config->getSterRepFac() * pWorkVar->partProp.rMag(idx1) * pWorkVar->partProp.rMag(idx1) * log(dDelta / r) * pWorkVar->intLists.intVecsCoords.row(i);
	}
	else
	{
		pBuffer->bufferVec.vec1.setZero();
	}

	//add to forces
	pWorkVar->intForce.row(idx1) -= pBuffer->bufferVec.vec1; 
	pWorkVar->intForce.row(idx2) += pBuffer->bufferVec.vec1;

	// add force to buffervec 2 for pressure tensor calculation
	pBuffer->bufferVec.vec2 += pBuffer->bufferVec.vec1;
}

// Brief: calculation of electrostatic force + adding pair interaction force to buffervec2 for pressure tensor calculation
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompElectrostatRep(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, int idx1, int idx2, int i)
{
	double r = pWorkVar->intLists.intVecsLen(i);
	
	pBuffer->bufferVec.vec1 = (pWorkVar->partProp.rMag(idx1) * pWorkVar->partProp.rMag(idx2) / (pWorkVar->partProp.rMag(idx1) + pWorkVar->partProp.rMag(idx2))) * config->getElectrostatRepFac()
		* exp(-config->getKappa() * (r - pWorkVar->partProp.rMag(idx1) - pWorkVar->partProp.rMag(idx2))) * pWorkVar->intLists.intVecsCoords.row(i);
	
	//add to forces
	pWorkVar->intForce.row(idx1) -= pBuffer->bufferVec.vec1;
	pWorkVar->intForce.row(idx2) += pBuffer->bufferVec.vec1;

	// add force to buffervec 2 for pressure tensor calculation
	pBuffer->bufferVec.vec2 += pBuffer->bufferVec.vec1;
}

// Brief: calculation of attractive van der Waals force + adding pair interaction force to buffervec2 for pressure tensor calculation
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompForceVdW(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, int idx1, int idx2, int i)
{
	double r = pWorkVar->intLists.intVecsLen(i);

	pBuffer->bufferVec.vec1 = -config->getHamaker() * r * pow(4 * pWorkVar->partProp.rMag(idx1) * pWorkVar->partProp.rMag(idx2), 3.0) / (6 *
		pow((r - pWorkVar->partProp.rMag(idx1) - pWorkVar->partProp.rMag(idx2)) * (r + pWorkVar->partProp.rMag(idx1) + pWorkVar->partProp.rMag(idx2)), 2.0) *
		pow(r * r - pow(pWorkVar->partProp.rMag(idx1) - pWorkVar->partProp.rMag(idx2), 2.0), 2.0)
		) * pWorkVar->intLists.intVecsCoords.row(i);
	 
	//add to forces
	pWorkVar->intForce.row(idx1) -= pBuffer->bufferVec.vec1; 
	pWorkVar->intForce.row(idx2) += pBuffer->bufferVec.vec1;

	// add force to buffervec 2 for pressure tensor calculation
	pBuffer->bufferVec.vec2 += pBuffer->bufferVec.vec1;
}
