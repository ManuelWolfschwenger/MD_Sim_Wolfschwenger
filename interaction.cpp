/*
* Name : interaction.cpp
* Auto :
* Brief: Contain functions for the interaction calculations of the dipoles
*/

// own header file
#include "interactions.h"

/////////////// neighbour list functions ////////////////

// Brief: This function initializes the parameters for the cell and neighbour list
// param[in]: Params_S* pParams - pointer to a Params_S object
// param[in]: PartProps_S* pPartProp - pointer to a PartProps_S object
// param[in/out]: NebrList_S* pLists - pointer to a NebrList_S object
// return: void
void InitNebrList(Params_S* pParams, NebrList_S* pLists, PartProps_S* pPartProp)
{
	double mMax;       //maximum magnetic moment
	double rMax;       //maximum magnetic radius
	double c1; // constant for F_vdw/F_dip
	double errNewton; // error for finding root of function from above
	
	mMax = pPartProp->magMom.maxCoeff(); // biggest particle
	rMax = pPartProp->rMag.maxCoeff();

	if (ENABLE_COATING_POT)
	{
		// find value where F_vdw/F_dip = errTol to cut off short range forces
		c1 = 128 * hamaker * pow(rMax, 6.0) * my_pi / (18 * mu0 * mMax * mMax);
		errNewton = 100;
		pLists->rCutSR = 1.1 * 2 * rMax; //start value for iteration

		// Newton iteration for finding root of function
		do
		{
			pLists->rCutSR -= FunCutOffError(pLists->rCutSR, c1, errTolSR, rMax, false) / FunCutOffError(pLists->rCutSR, c1, errTolSR, rMax, true);
			errNewton = abs(FunCutOffError(pLists->rCutSR, c1, errTolSR, rMax, false));

		} while (errNewton > 1e-7);
		
		////////////// now we know the cut off radius to make a desired error of errtolSR in Fvdw/Fdip/////////////////

		if (pLists->rCutSR < 2 * pPartProp->rHydr.maxCoeff())
		{
			cout << "error particles can only interact if they overlap => edit cut off radius" << endl;
			exit(-1);
		}
	}
	else // if the short range interactions are treated with the hard sphere potential
	{
		pLists->rCutSR = 2*rMax;
	}

	//thickness of neighbour shell
	pLists->dNebrShell = nebrShellFac * rMagMean;

	//maximum radius for possible neighbours
	pLists->rNebr = pLists->rCutSR + pLists->dNebrShell;

	//calculate cell geometry
	pLists->nCells = (int)floor(pParams->lBox / pLists->rNebr);

	if (pLists->nCells < 3)
	{
		cout << "number of cells too small, increase number of particles or decrease neighbour circle!" << "\n" << endl;
		exit(EXIT_FAILURE); //stops the simulation
	}

	pLists->numbCells = (int)pow(pLists->nCells, 3);
	pLists->lCell = pParams->lBox / pLists->nCells; //real length of cell

	//initialize offset array
	pLists->offset.row(0) << 0, 0, 0;
	pLists->offset.row(1) << 1, 0, 0;
	pLists->offset.row(2) << 1, 1, 0;
	pLists->offset.row(3) << 0, 1, 0;
	pLists->offset.row(4) << -1, 1, 0;
	pLists->offset.row(5) << 0, 0, 1;
	pLists->offset.row(6) << 1, 0, 1;
	pLists->offset.row(7) << 1, 1, 1;
	pLists->offset.row(8) << 0, 1, 1;
	pLists->offset.row(9) << -1, 1, 1;
	pLists->offset.row(10) << -1, 0, 1;
	pLists->offset.row(11) << -1, -1, 1;
	pLists->offset.row(12) << 0, -1, 1;
	pLists->offset.row(13) << 1, -1, 1;

	// resize cell and neighbour list
	pLists->cellList.resize(numbPart + pLists->numbCells);

	pLists->dispSum = 0;    // init sum of displacements to zero
}

// Brief: calculation ofF_vdw/F_dip
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

// Brief: This function builds the neighbour list
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Params_S* pParams - pointer to a Params_S object
// param[in/out]: NebrList_S* pLists - pointer to a NebrList_S object
// return: void
void BuildNebrList(WorkingVar_S* pWorkVar, Params_S* pParams, NebrList_S* pLists)
{
	// create cell list
	// linked list:header is numbPart long following by numbCells
	pLists->cellList.setConstant(-2);                       //just to mark that it can not be a particle if something goes wrong later
	pLists->cellList.bottomRows(pLists->numbCells).setConstant(-1); //fill all cells with -1

	pLists->cx = ((pWorkVar->coordsTrans.pos.col(0) + 0.5 * pParams->lBox) / pLists->lCell).cast<int>(); //cells x indizes
	pLists->cy = ((pWorkVar->coordsTrans.pos.col(1) + 0.5 * pParams->lBox) / pLists->lCell).cast<int>(); //cells y indizes
	pLists->cz = ((pWorkVar->coordsTrans.pos.col(2) + 0.5 * pParams->lBox) / pLists->lCell).cast<int>(); //cells z indizes
	pLists->c = (pLists->cz * pLists->nCells + pLists->cy) * pLists->nCells + pLists->cx + numbPart; //cell indizes

	// linked list:header is numbPart long following by numbCells
	for (int i = 0; i < numbPart; i++)
	{
		pLists->cellList(i) = pLists->cellList(pLists->c(i));
		pLists->cellList(pLists->c(i)) = i;
	}

	// create neighbour-list
	Array<int, 1, 3> m1Vec;
	Array<int, 1, 3> m2Vec;
	Array<double, 1, 3> shiftVec;

	int m1, m2;
	int j1, j2;
	int nebrIdx = 0;

	double dr;

	//init maximum size of neighbour list
	pLists->nebrList.resize((numbPart - 1) * numbPart);

	for (int mz = 0; mz < pLists->nCells; mz++)
	{
		for (int my = 0; my < pLists->nCells; my++)
		{
			for (int mx = 0; mx < pLists->nCells; mx++)
			{
				m1Vec << mx, my, mz;
				m1 = (mz * pLists->nCells + my) * pLists->nCells + mx + numbPart; //index of cell 1

				for (int off = 0; off < 14; off++)
				{
					m2Vec = m1Vec + pLists->offset.row(off);

					shiftVec.setZero();

					CellWrapAround(&m2Vec, &shiftVec, pParams, pLists, 0);
					CellWrapAround(&m2Vec, &shiftVec, pParams, pLists, 1);
					CellWrapAround(&m2Vec, &shiftVec, pParams, pLists, 2);

					//now we know the index of the second cell
					m2 = (m2Vec(2) * pLists->nCells + m2Vec(1)) * pLists->nCells + m2Vec(0) + numbPart; //index of cell 2

					j1 = pLists->cellList(m1);

					while (j1 >= 0) //if its -1 there is no further particle in the cell
					{
						j2 = pLists->cellList(m2);

						while (j2 >= 0)
						{
							if (m1 != m2 || j2 < j1)
							{
								dr = (pWorkVar->coordsTrans.pos.row(j1) - pWorkVar->coordsTrans.pos.row(j2) - shiftVec).matrix().norm(); //distance between two particles

								if (dr < pLists->rNebr)
								{
									//write idx of neighbours which are in the interaction range to list
									pLists->nebrList(2 * nebrIdx) = j1;
									pLists->nebrList(2 * nebrIdx + 1) = j2;
									nebrIdx++;
								}
							}

							j2 = pLists->cellList(j2);
						}

						j1 = pLists->cellList(j1);
					}
				}
			}
		}
	}

	pLists->pairsSR = nebrIdx; //pairs contained in nebrlist

	//resize list to the actual size of pairs
	pLists->nebrList.conservativeResize(nebrIdx * 2);

	// resize array for interaction vectors
	pLists->intVecsCoordsSR.resize(pLists->pairsSR, 3);   // resize to number of pairs rows and 3 colums for x,y,z and 1 column for the vector length
	pLists->intVecsLenSR.resize(pLists->pairsSR);         // length of vectors
	pLists->intRangeSR.resize(pLists->pairsSR);           // is r < rCut?
}

// Brief: cell wrapping, applies periodic boundary conditions to cells
// param[in] :Params_S* pParams - pointer to Params_S object
// param[in/out]: NebrList_S* pLists - pointer to a NebrList_S object
// param[in] :int col - index of vector column
// param[in/out] : pM2Vec - pointer to m2Vec
// param[in/out] : pShift - pointer to shifVec
// return: void
void CellWrapAround(Array<int, 1, 3>* pM2Vec, Array<double, 1, 3>* pShiftVec, Params_S* pParams, NebrList_S* pLists, int col)
{
	// checks if m2Vec has a component exceeding the cell grid bounds 
	if ((*pM2Vec)(col) >= pLists->nCells)
	{
		(*pM2Vec)(col) = 0;
		(*pShiftVec)(col) = pParams->lBox;
	}
	if ((*pM2Vec)(col) < 0)
	{
		(*pM2Vec)(col) = pLists->nCells - 1;
		(*pShiftVec)(col) = -pParams->lBox;
	}
}

// Brief: This function checks if the neighbour list has to be refreshed
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Params_S* pParams - pointer to a Params_S object
// param[in]: double deltaTused - used timestep
// return: void
void RefreshNebrList(WorkingVar_S* pWorkVar, Params_S* pParams, double deltaTused)
{
	double vMax;

	vMax = (pWorkVar->coordsTrans.vel.rowwise().norm()).maxCoeff(); //maximum velocity of particles

	pWorkVar->nebrListData.dispSum += vMax * deltaTused; //maximum displacement

	// build a new list 
	if (pWorkVar->nebrListData.dispSum > 0.5 * pWorkVar->nebrListData.dNebrShell)
	{
		BuildNebrList(pWorkVar, pParams, &pWorkVar->nebrListData);
		pWorkVar->nebrListData.dispSum = 0;

		//cout << "\n" << "neighbour list refreshed" << "\n" << endl;
	}
}


//////////////////////// interaction functions ///////////////////////////

// Brief: This function computes the interaction vector for all interaction pairs and writes the data to an array
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: Params_S* pParams - pointer to params_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompInteractVecsLR(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams)
{
	double r;
	int pairsLR = 0;

	// have to go through all pairs because rcCutLR > L/3
	for (int idx1 = 0; idx1 < numbPart - 1; idx1++)
	{
		for (int idx2 = idx1 + 1; idx2 < numbPart; idx2++)
		{
			//vector pointing from dipole i to dipole j
			pBuffer->bufferVec.vec_r = pWorkVar->coordsTrans.pos.row(idx2) - pWorkVar->coordsTrans.pos.row(idx1);

			// evaluate shortest distance (across boundaries or not)
			ApplyBoundaryCondVec(&pBuffer->bufferVec.vec_r, &pBuffer->bufferVec.vec_b, pParams);

			r = pBuffer->bufferVec.vec_r.matrix().norm(); //distance between centerpoints of particle i and j

			if (r < pParams->rCutLR)
			{
				//store vectorlength and vector coordinates in interaction arrays for vectors
				pWorkVar->intListsLR.intVecsLenLR.row(pairsLR) = r;
				pWorkVar->intListsLR.intVecsCoordsLR.row(pairsLR) = pBuffer->bufferVec.vec_r / r;
				pWorkVar->intListsLR.intPairsLR(2 * pairsLR) = idx1;
				pWorkVar->intListsLR.intPairsLR(2 * pairsLR + 1) = idx2;
				pairsLR++;
			}
		}
	}

	pWorkVar->intListsLR.pairsLR = pairsLR;
}

// Brief: This function computes the interaction vector for one interaction pair and writes the data to an array 
// This function is only suitable for the neighbour list method, can be possible included in distance calculation for long range interactions!!!!!!!!!!!!!!!!!!!!!!!!!
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: Params_S* pParams - pointer to params_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompInteractVecNL(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams, int idx1, int idx2, int i) //This function is only suitable for the neighbour list method!!!
{
	double r;

	//vector pointing from dipole i to dipole j
	pBuffer->bufferVec.vec_r = pWorkVar->coordsTrans.pos.row(idx2) - pWorkVar->coordsTrans.pos.row(idx1);

	// evaluate shortest distance (across boundaries or not)
	ApplyBoundaryCondVec(&pBuffer->bufferVec.vec_r, &pBuffer->bufferVec.vec_b, pParams);

	r = pBuffer->bufferVec.vec_r.matrix().norm(); //distance between centerpoints of particle i and j

	// check if distance is in interaction range
	if (r <= pWorkVar->nebrListData.rCutSR)
	{
		//store vectorlength and vector coordinates in interaction arrays for vectors
		pWorkVar->nebrListData.intRangeSR.row(i) = true;
		pWorkVar->nebrListData.intVecsLenSR.row(i) = r;
		pWorkVar->nebrListData.intVecsCoordsSR.row(i) = pBuffer->bufferVec.vec_r / r;
	}
	else
	{
		pWorkVar->nebrListData.intRangeSR.row(i) = false;
	}
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
	int idx1, idx2, minIdx, maxIdx, tIter;
	const int steps = 10000;
	const int rcSteps = 30;
	tIter = 20;

	ArrayXXd funVals, paramArr; // should be stored on heap
	funVals.resize(steps, 2);
	paramArr.resize(rcSteps, 4); // rc, time, alpha, kc

	// constants from particle properties
	M = (pWorkVar->partProp.magMom.pow(2.0)).sum();
	fAl = 3 * mu0 * pow(pWorkVar->partProp.magMom.minCoeff(), 2.0) / (2.0 * my_pi * pow(2.0 * pWorkVar->partProp.rHydr.minCoeff(), 4.0));  //force between aligned dipoles
	tAl = pow(pWorkVar->partProp.magMom.minCoeff(), 2.0) * mu0 / (16.0 * my_pi * pow(pWorkVar->partProp.rHydr.minCoeff(), 3.0)); //torque on dip1 in field from dip2 when in contract and magMoments 90degrees to each other

	paramArr.col(0).setLinSpaced(rcSteps, 0.2, 0.5); //variation of rcut to try

	for (int j = 0; j < rcSteps; j++)
	{
		// set rCut for real space
		pParams->rCutLR = pParams->lBox * paramArr(j,0);

		////////////// finding alpha for specific errTol ///////////////////
		//
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

		paramArr(j, 2) = -funVals(idx1, 1) / (funVals(idx2, 1) - funVals(idx1, 1)) * (funVals(idx2, 0) - funVals(idx1, 0)) + funVals(idx1, 0);
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
			cout << "error: maximum value of error < 0 in kc computation => decrease error tolerance!" << endl;
			cout << "rc was L*" << paramArr(j, 0) << endl;
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

		paramArr(j, 3) = -funVals(idx1, 1) / (funVals(idx2, 1) - funVals(idx1, 1)) * (funVals(idx2, 0) - funVals(idx1, 0)) + funVals(idx1, 0);
		pParams->nCutLR = paramArr(j, 3);

		// ===>>>> now we know nCutLR for a given error tolerance!!!!!

		//////////////////////// measuring simulation time of ewald summation ////////////////////////////////////

		pWorkVar->intForce.setZero();         //force array, initialized to zero
		pWorkVar->demagFluxDens.setZero();    //demagnetizing array, initialized to zero
		pWorkVar->partProp.magMomVec = pWorkVar->coordsRot.posMm.colwise() * pWorkVar->partProp.magMom; //create magnetic moment vectors
		tMean = 0;

		// evaluate interaction vectors
		CompInteractVecsLR(pWorkVar, pBuffer, pParams);

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
	cout << std::setprecision(10);
	cout << "        rc        time        alpha          kc" << endl;
	cout << paramArr << endl;
	cout << "\n" << endl;
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

	double errRealForce = facEwald * M / sqrt(pow(pParams->lBox, 3.0) * pow(pParams->alpha, 4.0) * pow(pParams->rCutLR, 9.0) * numbPart) * sqrt(13.0 / 6.0 * cc * cc + 2.0 / 15.0 * dc * dc - 13.0 / 15.0 * cc * dc) * exp(-pow(pParams->alpha * pParams->rCutLR, 2.0));
	double errRezForce = facEwald * 8.0 * my_pi * M / pow(pParams->lBox, 3.0) * pParams->alpha * sqrt(2.0 * my_pi * pow(pParams->nCutLR, 3.0) / (15.0 * numbPart)) * exp(-pow(my_pi * pParams->nCutLR / (pParams->alpha * pParams->lBox), 2.0));
	double errRealTorque = facEwald * M / sqrt(pow(pParams->lBox, 3.0) * pow(pParams->alpha, 4.0) * pow(pParams->rCutLR, 7.0) * numbPart) * sqrt(0.5 * bc * bc + 0.2 * cc * cc) * exp(-pow(pParams->alpha * pParams->rCutLR, 2.0));
	double errRezTorque = facEwald * 4.0 * M / pow(pParams->lBox, 2.0) * pParams->alpha * sqrt(my_pi * pParams->nCutLR / (5.0 * numbPart)) * exp(-pow(my_pi * pParams->nCutLR / (pParams->alpha * pParams->lBox), 2.0));
	
	cout << "the resulting force error should be: " << sqrt(errRezForce * errRezForce + errRealForce * errRealForce) / fAl << "*F_contact" << endl;
	cout << "the resulting torque error should be: " << sqrt(errRezTorque * errRezTorque + errRealTorque * errRealTorque) / tAl << "*T_contact" << "\n" << endl;
}

// Brief: calculation of rms error of real space forces
// param[in]:EwaldParams_S* pEwaldParams - pointer to a EwaldParams_S object
// param[in]: double x dependent variable
// return: void
double FunRealError(double x, double M, double fAl, Params_S* pParams)
{
	double a = x * pParams->rCutLR;
	double b = M / sqrt(pow(pParams->lBox, 3.0) * pow(pParams->rCutLR, 9.0) * numbPart * pow(x, 4.0));
	double c = 4.0 * pow(a, 4.0) + 6.0 * pow(a, 2.0) + 3.0;
	double d = 8.0 * pow(a, 6.0) + 20.0 * pow(a, 4.0) + 30.0 * pow(a, 2.0) + 15.0;

	return facEwald * b * sqrt(13.0 / 6.0 * c * c + 2.0 / 15.0 * d * d - 13.0 / 15.0 * c * d) * exp(-a * a) / fAl - errTolEwald / sqrt(2.0);
}

// Brief: calculation of rms error of fourier space forces
// param[in]:EwaldParams_S* pEwaldParams - pointer to a EwaldParams_S object
// param[in]: double x dependent variable
// return: void
double FunRezError(double x, double M, double fAl, Params_S* pParams)
{
	return  facEwald * 8.0 * M * my_pi / pow(pParams->lBox, 3.0) * pParams->alpha * sqrt(2.0 * my_pi * x * x * x / (15.0 * numbPart)) * exp(-pow(my_pi * x / (pParams->alpha * pParams->lBox), 2.0)) / fAl - errTolEwald / sqrt(2.0);
}

// Brief: precalculation of Constants for ewald summation
// param[in/out]: Params_S* pParams - pointer to Params_S object
// param[in/out]: Buffer_S* pBuffer - pointer to buffer_S object
// return: void
void InitEwaldConst(Params_S* pParams, Buffer_S* pBuffer)
{
	// initialize size of Ewald arrays in buffer structure
	pBuffer->tSin.resize(pParams->nCutLR + 1, numbPart);
	pBuffer->tCos.resize(pParams->nCutLR + 1, numbPart);

	// constants for fourier space
	pParams->nc = (int)pParams->nCutLR;
	pParams->nc2 = pParams->nCutLR * pParams->nCutLR;
	pParams->gu = 2 * my_pi * (pParams->lBox * pParams->lBox * pParams->lBox);
	pParams->gr = 4 * my_pi * pParams->gu / pParams->lBox;
	pParams->gs = 2 * pParams->gu;
	pParams->w = pow(my_pi / (pParams->lBox * pParams->alpha), 2.0);
	pParams->cFluxDens = 4.0 * my_pi * facEwald / pow(pParams->lBox, 3.0);
	pParams->cForce = 8.0 * my_pi * my_pi * facEwald / pow(pParams->lBox, 4.0);
	pParams->cSurf = facEwald * 4.0 * my_pi / (3 * pow(pParams->lBox, 3.0));   //(2*epsilon +1) = 3 for box surrounded by vacuum 
	pParams->cSelf = facEwald * 4.0 * pow(pParams->alpha, 3.0) / (3.0 * sqrt(my_pi));

	// constants for eal space
	pParams-> alpha2 = pParams->alpha * pParams->alpha;
	pParams->alpha4 = pParams->alpha2 * pParams->alpha2;
}

/////////// perform ewald summation /////////////////////

// Brief: calculation of dipole interaction force in real space / must be called in all pairs loop /also surface
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
// return: void
void CompDipInteractEwaldR(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams)
{
	int idx1, idx2;
	double r, r2, t, a1, a2, a3, dot1, dot2;
	double alpha = pParams->alpha;
	double alpha2 = pParams->alpha2;

	// loop for long range interaction / real part of ewald sum
	for (int i = 0; i < pWorkVar->intListsLR.pairsLR; i++)
	{
		idx1 = pWorkVar->intListsLR.intPairsLR(2 * i);
		idx2 = pWorkVar->intListsLR.intPairsLR(2 * i + 1);

		r = pWorkVar->intListsLR.intVecsLenLR(i);
		r2 = r * r;
		t = 2.0 * alpha * exp(-alpha2 * r2) * irPi / r2;
		a1 = erfc(alpha * r) / (r2 * r) + t;
		a2 = 3.0 * a1 / r2 + 2.0 * alpha2 * t;
		a3 = 5.0 * a2 / r2 + 4.0 * pParams->alpha4 * t;

		dot1 = pWorkVar->partProp.magMomVec.row(idx1).matrix().dot(pWorkVar->intListsLR.intVecsCoordsLR.row(i).matrix());
		dot2 = pWorkVar->partProp.magMomVec.row(idx2).matrix().dot(pWorkVar->intListsLR.intVecsCoordsLR.row(i).matrix());

		// force calculation
		pBuffer->bufferVec.vec1 = facEwald * ((a2 * pWorkVar->partProp.magMomVec.row(idx1).matrix().dot(pWorkVar->partProp.magMomVec.row(idx2).matrix()) -
			a3 * dot1 * dot2 * r2) * pWorkVar->intListsLR.intVecsCoordsLR.row(i) * r +
			a2 * r * (dot2 * pWorkVar->partProp.magMomVec.row(idx1) + dot1 * pWorkVar->partProp.magMomVec.row(idx2)));

		// flux density calculation
		pWorkVar->demagFluxDens.row(idx1) += facEwald * (a2 * dot2 * pWorkVar->intListsLR.intVecsCoordsLR.row(i) * r2 - a1 * pWorkVar->partProp.magMomVec.row(idx2)); //field on particle i;
		pWorkVar->demagFluxDens.row(idx2) += facEwald * (a2 * dot1 * pWorkVar->intListsLR.intVecsCoordsLR.row(i) * r2 - a1 * pWorkVar->partProp.magMomVec.row(idx1)); //field on particle j;

		// add to forces
		pWorkVar->intForce.row(idx1) -= pBuffer->bufferVec.vec1;
		pWorkVar->intForce.row(idx2) += pBuffer->bufferVec.vec1;
	}
}

// Brief: calculation of dipole interaction force in fourier space 
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
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
		
				for (int i = 0; i < numbPart; i++)
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

				for (int i = 0; i < numbPart; i++)
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
	for (int i = 0; i < numbPart; i++)
	{
		pWorkVar->demagFluxDens.row(i) += pParams->cSelf * pWorkVar->partProp.magMomVec.row(i) - pParams->cSurf * pWorkVar->partProp.magMomSum;
	}
}

// Brief: precalculation of sin and cos for fourier space
// param[in]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out]: Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: Params_S* pParams - pointer to Params_S object
// return: void
void EvalSinCos(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams)
{
	Array<double,1,3> t, tt, u, w;

	t.setConstant(2 * my_pi / pParams->lBox); //t...buffervec 1

	for (int i = 0; i < numbPart; i++)
	{
		tt = t * pWorkVar->coordsTrans.pos.row(i); //tt...buffervec 2
		pBuffer->tCos(0, i).setConstant(1);
		pBuffer->tSin(0, i).setZero();
		pBuffer->tCos(1, i) = tt.cos();
		pBuffer->tSin(1, i) = tt.sin();
		u = 2.0 * pBuffer->tCos(1, i); //u...buffervec 3
		pBuffer->tCos(2, i) = u * pBuffer->tCos(1, i);
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

// Brief: calculation of dipole interaction force in real space / must be called in all pairs loop
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
	for (int i = 0; i < pWorkVar->intListsLR.pairsLR; i++)
	{
		idx1 = pWorkVar->intListsLR.intPairsLR(2 * i);
		idx2 = pWorkVar->intListsLR.intPairsLR(2 * i + 1);

		r = pWorkVar->intListsLR.intVecsLenLR(i);
		r2 = r * r;
		t = 2.0 * alpha * exp(-alpha2 * r2) * irPi / r2;
		a1 = erfc(alpha * r) / (r2 * r) + t;
		a2 = 3.0 * a1 / r2 + 2.0 * alpha2 * t;
		a3 = 5.0 * a2 / r2 + 4.0 * pParams->alpha4 * t;

		dot1 = pWorkVar->partProp.magMomVec.row(idx1).matrix().dot(pWorkVar->intListsLR.intVecsCoordsLR.row(i).matrix());
		dot2 = pWorkVar->partProp.magMomVec.row(idx2).matrix().dot(pWorkVar->intListsLR.intVecsCoordsLR.row(i).matrix());

		// flux density calculation
		pWorkVar->demagFluxDens.row(idx1) += facEwald * (a2 * dot2 * pWorkVar->intListsLR.intVecsCoordsLR.row(i) * r2 - a1 * pWorkVar->partProp.magMomVec.row(idx2)); //field on particle i;
		pWorkVar->demagFluxDens.row(idx2) += facEwald * (a2 * dot1 * pWorkVar->intListsLR.intVecsCoordsLR.row(i) * r2 - a1 * pWorkVar->partProp.magMomVec.row(idx1)); //field on particle j;
	}
}

// Brief: calculation of dipole interaction force in fourier space 
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

				for (int i = 0; i < numbPart; i++)
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

				for (int i = 0; i < numbPart; i++)
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
	for (int i = 0; i < numbPart; i++)
	{
		pWorkVar->demagFluxDens.row(i) += pParams->cSelf * pWorkVar->partProp.magMomVec.row(i) - pParams->cSurf * pWorkVar->partProp.magMomSum;
	}
}

// Brief: this function calculates the exact force on the particles to check the ewald summation
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Buffer_S* pBuffer - pointer to buffer_S object
// return: void
void CheckEwaldSum(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams)
{
    ArrayXXd  exactForces, ewaldForces, exactFields, ewaldFields, exactTorques, ewaldTorques, diff;
	double errIs, fAl, tAl, bAl;

	exactForces.resize(numbPart, 3);
	ewaldForces.resize(numbPart, 3);
	exactFields.resize(numbPart, 3);
	ewaldFields.resize(numbPart, 3);
	exactTorques.resize(numbPart, 3);
	ewaldTorques.resize(numbPart, 3);
	diff.resize(numbPart, 3);

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
	cout << "attention: restart simumlation after checking Ewald sum and disable it, because arrays then have the wrong size" << endl;
	EvalSinCos(pWorkVar, pBuffer, pParams);

	//Ewald Summation
	CompDipInteractEwaldR(pWorkVar, pBuffer, pParams);
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
	errIs = sqrt(((diff * diff).rowwise().sum()).sum() / numbPart); // rms error made by ewald summation
	fAl = 3 * mu0 * pow(pWorkVar->partProp.magMom.minCoeff(), 2.0) / (2 * my_pi * pow(2 * pWorkVar->partProp.rHydr.minCoeff(), 4.0));  //force between aligned dipoles in contact
	cout << "the actual rms force error is " << errIs / fAl<< "*F_contact" << endl;

	diff = exactFields - ewaldFields;
	errIs = sqrt(((diff * diff).rowwise().sum()).sum() / numbPart); // rms error made by ewald summation
	bAl =  mu0 * pWorkVar->partProp.magMom.minCoeff() / (2 * my_pi * pow(2 * pWorkVar->partProp.rHydr.minCoeff(), 3.0));  //field of dipole i at dipole j when in contact
	cout << "the actual rms field error is " << errIs / bAl << "*B_contact" << endl;

	diff = exactTorques - ewaldTorques;
	errIs = sqrt(((diff * diff).rowwise().sum()).sum() / numbPart); // rms error made by ewald summation
	tAl = pow(pWorkVar->partProp.magMom.minCoeff(), 2.0) * mu0 / (16.0 * my_pi * pow(pWorkVar->partProp.rHydr.minCoeff(), 3.0));
	cout << "the actual rms torque error is " << errIs / tAl << "*T_contact" << endl;
}


/////////////////////// short range interactions //////////////////////////////////////////

// Brief: calculation of repulsive force (Rosensweig)
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompSterRep(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, int idx1, int idx2, int i)
{
	double dDelta = 2* (pWorkVar->partProp.rMag(idx1) + effLenSurfMol);
	double r = pWorkVar->nebrListData.intVecsLenSR(i);

	// repulsive potential
	if (r <= dDelta)
	{
		pBuffer->bufferVec.vec1 = (2 * pWorkVar->partProp.rMag(idx2) / (pWorkVar->partProp.rMag(idx1) + pWorkVar->partProp.rMag(idx2))) * //correction for spheres of differnt size
			sterRepFac * pWorkVar->partProp.rMag(idx1) * pWorkVar->partProp.rMag(idx1) * log(dDelta / r) * pWorkVar->nebrListData.intVecsCoordsSR.row(i);
	}
	else
	{
		pBuffer->bufferVec.vec1.setZero();
	}

	//add to forces
	pWorkVar->intForce.row(idx1) -= pBuffer->bufferVec.vec1; 
	pWorkVar->intForce.row(idx2) += pBuffer->bufferVec.vec1;
}

// Brief: calculation of electrostatic force 
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompElectrostatRep(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, int idx1, int idx2, int i)
{
	double r = pWorkVar->nebrListData.intVecsLenSR(i);
	
	pBuffer->bufferVec.vec1 = (pWorkVar->partProp.rMag(idx1) * pWorkVar->partProp.rMag(idx2) / (pWorkVar->partProp.rMag(idx1) + pWorkVar->partProp.rMag(idx2))) * electrostatRepFac
		* exp(-kappa * (r - pWorkVar->partProp.rMag(idx1) - pWorkVar->partProp.rMag(idx2))) * pWorkVar->nebrListData.intVecsCoordsSR.row(i);
	
	//add to forces
	pWorkVar->intForce.row(idx1) -= pBuffer->bufferVec.vec1;
	pWorkVar->intForce.row(idx2) += pBuffer->bufferVec.vec1;
}

// Brief: calculation of attractive van der Waals force
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CompForceVdW(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, int idx1, int idx2, int i)
{
	double r = pWorkVar->nebrListData.intVecsLenSR(i);

	pBuffer->bufferVec.vec1 = -hamaker * r * pow(4 * pWorkVar->partProp.rMag(idx1) * pWorkVar->partProp.rMag(idx2), 3.0) / (6 *
		pow((r - pWorkVar->partProp.rMag(idx1) - pWorkVar->partProp.rMag(idx2)) * (r + pWorkVar->partProp.rMag(idx1) + pWorkVar->partProp.rMag(idx2)), 2.0) *
		pow(r * r - pow(pWorkVar->partProp.rMag(idx1) - pWorkVar->partProp.rMag(idx2), 2.0), 2.0)
		) * pWorkVar->nebrListData.intVecsCoordsSR.row(i);
	 
	//add to forces
	pWorkVar->intForce.row(idx1) -= pBuffer->bufferVec.vec1; 
	pWorkVar->intForce.row(idx2) += pBuffer->bufferVec.vec1;
}

// Brief: check if particles overlap + correcting
// param[in/out] : WorkingVar_S* pWorkVar - pointer to workVar_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void HardSphereRep(WorkingVar_S* pWorkVar, int idx1, int idx2, int i)
{
	double r = pWorkVar->nebrListData.intVecsLenSR(i);

	// checks if particles overlaps
	if (r <= (pWorkVar->partProp.rHydr(idx1) + pWorkVar->partProp.rHydr(idx2)))
	{
		double d = pWorkVar->partProp.rHydr(idx1) + pWorkVar->partProp.rHydr(idx2) - r; //overlapping distance

		pWorkVar->coordsTrans.pos.row(idx1) -= 0.5 * d * pWorkVar->nebrListData.intVecsCoordsSR.row(i); //move back half the overlapping distance
		pWorkVar->coordsTrans.pos.row(idx2) += 0.5 * d * pWorkVar->nebrListData.intVecsCoordsSR.row(i);
	}
}

// Brief: check if system is stable + correcting (instable if particle cores would overlap)
// param[in/out] : WorkingVar_S* pWorkVar - pointer to workVar_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void CheckStability(WorkingVar_S* pWorkVar, int idx1, int idx2, int i)
{
	double d =  pWorkVar->partProp.rMag(idx1) + pWorkVar->partProp.rMag(idx2) + 1e-10 - pWorkVar->nebrListData.intVecsLenSR(i); //overlapping distance + safety for singularity in vdw attraction

	// checks if particles overlap
	if (d > 0) 
	{
		pWorkVar->nebrListData.intVecsLenSR(i) = pWorkVar->partProp.rMag(idx1) + pWorkVar->partProp.rMag(idx2) + 2e-10;
		pWorkVar->coordsTrans.pos.row(idx1) -= 0.5 * d * pWorkVar->nebrListData.intVecsCoordsSR.row(i); //move back half the overlapping distance
		pWorkVar->coordsTrans.pos.row(idx2) += 0.5 * d * pWorkVar->nebrListData.intVecsCoordsSR.row(i);
	}
}
