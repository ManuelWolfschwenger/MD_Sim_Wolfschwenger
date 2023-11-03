/*
* Name : simulation.cpp
* Auto :
* Brief: Simulation source file
*/

// own header file
#include "simulation.h"

// global variables
static WorkingVar_S     workVar;    // simulation working varibales
static Buffer_S         buffer;     // buffer 
static IntCoefficant_S  intCoeff;   // integration coefficant varibales
static OutputVar_S      outputVar;  // output variables
static Params_S         params;     // global parameters not known at compile time

// Brief: Runs the complete simulation
// Return: program status
int main(void)
{
	// initialize the simulation varibales
	InitSimulation();

	// run the simulation
	RunSimulation();
}

////////////////////////*initialization functions*///////////////////////////////////////

// Brief: This function initializes the simulation
// return: void
void InitSimulation(void)
{
	WorkingVar_S* pWorkVar;      // pointer to the simulation working variables
	Buffer_S* pBuffer;           // pointer to buffer
	IntCoefficant_S* pIntCoeff;  // pointer to the integration coefficant variables
	Params_S* pParams;           // pointer to parameters not known at compile time

	// initialize the pointers
	pWorkVar = &workVar;
	pBuffer = &buffer;
	pIntCoeff = &intCoeff;
	pParams = &params;
	
	// initialize the particle properties
	InitParticleProps(&pWorkVar->partProp);

	// set Params
	pParams->lBox = pow(pWorkVar->partProp.volMag.sum() / volFrac, 1.0 / 3.0);  //box volume^1/3

	if (READ_FROM_FILE)
	{
		// continue the simulation from before
		InitCoordsFromSim(pWorkVar, pParams);
	} 
	else
	{
		// initialize the particle starting positions from random numbers
		InitCoordsFromSketch(pWorkVar, pParams);
	}
	
	// initialize the integration coefficient buffers
	InitIntCoeff(pIntCoeff);

	// initialize the buffers
	InitBuffers(pBuffer);

	// initialize neighbour list
	InitNebrList(pParams, &pWorkVar->nebrListData, &pWorkVar->partProp);

	// build neighbour list
	BuildNebrList(pWorkVar, pParams, &pWorkVar->nebrListData);
	
	// resizing interaction arrays
	int maxPairs = (numbPart * (numbPart - 1)) / 2;             // all pairs
	pWorkVar->intListsLR.intVecsCoordsLR.resize(maxPairs, 3);   // resize to number of pairs rows and 3 colums for x,y,z
	pWorkVar->intListsLR.intVecsLenLR.resize(maxPairs);         // length of vectors
	pWorkVar->intListsLR.intPairsLR.resize(2 * maxPairs);       // indizes of pairs

	// initialize the external field 
	pWorkVar->extFluxDens.resize(numbPart, 3);
	pWorkVar->extFluxDens.setZero();  //in case we start with a relaxation, dont need to set to zero every iteration in the loop

	// initialize demag field
	pWorkVar->demagFluxDens.resize(numbPart, 3);
	pWorkVar->demagFluxDens.setZero();

	// initialize forces
	pWorkVar->intForce.resize(numbPart, 3);
	pWorkVar->intForce.setZero();

	// initialize the thermal fluctuations
	pWorkVar->thermTorque.resize(numbPart, 3);
	pWorkVar->thermField.resize(numbPart, 3);
	pWorkVar->thermForce.resize(numbPart, 3);

	// compute optimal Ewald-Parameters
	CompEwaldParams(pWorkVar, pBuffer, pParams);
	InitEwaldConst(pParams, pBuffer);

	//store data in txt
	WriteData2TXT(pWorkVar, pParams, "data.txt");
}

// Brief: This function sets the particle properties
// param[in/out]: PartProps_S* pPartProps - pointer to a PartProps_S object
// return: void
void InitParticleProps(PartProps_S* pPartProps)
{
	pPartProps->rMag.resize(numbPart);
	pPartProps->rHydr.resize(numbPart);
	pPartProps->volMag.resize(numbPart);
	pPartProps->volHydr.resize(numbPart);
	pPartProps->magMom.resize(numbPart);
	pPartProps->zetaRot.resize(numbPart);
	pPartProps->zetaTrans.resize(numbPart);
	pPartProps->velEaConst1.resize(numbPart);
	pPartProps->velEaConst2.resize(numbPart);
	pPartProps->magMomVec.resize(numbPart, 3); //array with magnetic moment vectors
	pPartProps->magMomSum.resize(1, 3);        //array with summed components

	// switch depending on the size distribution
	switch (SIZE_DIST)
	{
	case SIZE_DIST_EQUAL:
	{
		// set all particle radii to same value
		pPartProps->rMag.setConstant(rMagMean);

		if (ENABLE_COATING_POT)
		{
			// calculate equivalent hard sphere diameter from repulsive potential	
			//(integration works till dHydr = dMag, when the steric potential is further decreased the equivalent hard sphere diameter is equal to the magnetic diameter (Rosensweig)
			double dHydr = IntegralEffDiameter(2 * rMagMean);

			if (dHydr >= 2*rMagMean )
			{
				pPartProps->rHydr.setConstant(dHydr / 2.0);
			}
			else
			{
				pPartProps->rHydr.setConstant(rMagMean);
			}		
		}
		else
		{
			// use thickness of shell to add it to magnetic radius
			pPartProps->rHydr.setConstant(rMagMean + dShell);
		}

		pPartProps->volMag.setConstant(4.0 / 3.0 * my_pi * rMagMean* rMagMean* rMagMean);                                           // mag volume
		pPartProps->volHydr.setConstant(4.0 / 3.0 * my_pi * pPartProps->rHydr(0)* pPartProps->rHydr(0)* pPartProps->rHydr(0));      // hydrdyn volume
		pPartProps->magMom.setConstant(satMag * pPartProps->volMag(0));                                                             // mag Moment
		pPartProps->zetaRot.setConstant(8.0 * my_pi * vis * pPartProps->rHydr(0)* pPartProps->rHydr(0)* pPartProps->rHydr(0));      // rotational friction coefficient
		pPartProps->zetaTrans.setConstant(6.0 * my_pi * vis * pPartProps->rHydr(0));                                                // translational friction coefficient
		pPartProps->velEaConst1.setConstant(2.0 * anisEn * (pPartProps->volMag(0) / pPartProps->zetaRot(0)));                       // constants for later DGL solving
		pPartProps->velEaConst2.setConstant(1/pPartProps->zetaRot(0));                                                              // constants for later DGL solving
	}
	break;

	case SIZE_DIST_LOG:
	{
		// log-normal distributed sizes
		double mu = log(pow(rMagMean, 2.0) / sqrt(pow(rMagMean, 2.0) + pow(sigma, 2.0)));
		double sigmaLog = sqrt(log(1.0 + pow(sigma / rMagMean , 2.0)));

		// lognormal distribution with mean = 0, stdev = 1.0
		pPartProps->rMag = Rand::lognormal<MatrixXd>(numbPart, 1, generator, mu, sigmaLog);

		if (ENABLE_COATING_POT)
		{
			double dHydr;

			// calculate equivalent hard sphere diameter from repulsive potential
			for (int i = 0; i < numbPart; i++)
			{
				dHydr = IntegralEffDiameter(2 * pPartProps->rMag(i));

				if (dHydr >= 2 * pPartProps->rMag(i))
				{
					pPartProps->rHydr(i) = dHydr / 2.0;
				}
				else
				{
					pPartProps->rHydr(i) = pPartProps->rMag(i);
				}
			}
		}
		else
		{
			// use thickness of shell to add it on magnetic radius
			pPartProps->rHydr = pPartProps->rMag + dShell;
		}

		pPartProps->volMag = 4.0 / 3.0 * my_pi * pPartProps->rMag.pow(3.0);                                       // mag volume
		pPartProps->volHydr = 4.0 / 3.0 * my_pi * pPartProps->rHydr.pow(3.0);                                     // hydrdyn volume
		pPartProps->magMom = satMag * pPartProps->volMag;                                                         // mag Moment
		pPartProps->zetaRot = 8.0 * my_pi * vis * pPartProps->rHydr.pow(3.0);                                     // rotational friction coefficient
		pPartProps->zetaTrans = 6.0 * my_pi * vis * pPartProps->rHydr;                                            // translational friction coefficient
		pPartProps->velEaConst1 = 2.0 * anisEn * (pPartProps->volMag / pPartProps->zetaRot);                      // constants for later DGL solving
		pPartProps->velEaConst2 = pPartProps->zetaRot.cwiseInverse();                                             // constants for later DGL solving
	}
	break;
	}
}         

// Brief: This function sets the initial particle positions
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Params_S* pParams - pointer to params_S object
// return: void
void InitCoordsFromSketch(WorkingVar_S* pWorkVar, Params_S* pParams)
{
	/////////////////////translational positions/////////////////////////////////

	//define size of the dynamic arrays
	pWorkVar->coordsTrans.pos.resize(numbPart, 3);
	pWorkVar->coordsTrans.vel.resize(numbPart, 3);
    pWorkVar->coordsTrans.vel.setZero();   //set initial velocities to zero

	switch (INIT_CONFIG_TRANS)
	{
	case 1: //random positions
	{
		Array<double, 1, 3> xyzVec;
		Array<double, 1, 3> xyzDiff;

		double rDiff;

		for (int i = 0; i < numbPart; i++)
		{
		setRandomCoords: //goto statement

			xyzVec = xyzVec.setRandom() * (pParams->lBox / 2.0 - pWorkVar->partProp.rHydr(i)); 

			if (i != 0) // first particle can be everywhere
			{
				for (int j = 0; j < i ; j++) //all particles exceptional the one with number ip
				{
					xyzDiff = pWorkVar->coordsTrans.pos.row(j) - xyzVec;
					rDiff = xyzDiff.row(0).matrix().norm(); //distance of center points

					if (rDiff < pWorkVar->partProp.rHydr(i) + pWorkVar->partProp.rHydr(j))
					{
						goto setRandomCoords;
					}
				}
			}

			pWorkVar->coordsTrans.pos.row(i) = xyzVec;

		}

		break;
	}

	case 2: //positions in a grid
	{
		// check if number of particles per dimension (x,y,z) is integer

		double npDim = pow(numbPart, 1.0 / 3.0);
		if (floor(npDim) == npDim)
		{
			cout << "grid can be created" << endl;
		}
		else
		{
			cout << "grid can not be created, edit number of particles!" << endl;
		}

		double gap = pParams->lBox / npDim; //gap between particles
		int ip = 0;

		for (int ix = 0; ix < npDim; ix++)
		{
			for (int iy = 0; iy < npDim; iy++)
			{
				for (int iz = 0; iz < npDim; iz++)
				{
					pWorkVar->coordsTrans.pos.row(ip).col(0) = (ix + 0.5) * gap - pParams->lBox / 2.0;
					pWorkVar->coordsTrans.pos.row(ip).col(1) = (iy + 0.5) * gap - pParams->lBox / 2.0;
					pWorkVar->coordsTrans.pos.row(ip).col(2) = (iz + 0.5) * gap - pParams->lBox / 2.0;
					ip++;
				}
			}
		}

		break;
	}
	}

	/////////////////////rotational positions/////////////////////////////////

	//azimuth and elevation angle
	ArrayXd azimuth(numbPart);
	ArrayXd elevation(numbPart);

	//define size of the dynamic arrays
	pWorkVar->coordsRot.posEa.resize(numbPart, 3);
	pWorkVar->coordsRot.posMm.resize(numbPart, 3);

	// create uniformly random distributed vector components in a sphere
	azimuth.setRandom(); //random numbers between -1 and 1
	azimuth = (azimuth + 1.0) * my_pi;//uniformly distributed random numbers between 0 and 2*pi

	elevation.setRandom(); //uniformly distributed random numbers between 0 and 1
	elevation = (elevation + 1.0) / 2.0; 	
	elevation = elevation.cwiseSqrt().asin() * 2.0;

	// switch depending on the initial configuration
	switch (INIT_CONFIG_ROT)
	{
	case 1:
	{
		pWorkVar->coordsRot.posEa.col(0) = elevation.sin() * azimuth.cos();
		pWorkVar->coordsRot.posEa.col(1) = elevation.sin() * azimuth.sin();
		pWorkVar->coordsRot.posEa.col(2) = elevation.cos();

		// all moments point in z direction / easy axis random see above
		pWorkVar->coordsRot.posMm.col(0).setZero(); //x-direction
		pWorkVar->coordsRot.posMm.col(1).setZero(); //y-direction
		pWorkVar->coordsRot.posMm.col(2).setOnes(); //z-direction
		break;
	}
	case 2:
	{
		pWorkVar->coordsRot.posEa.col(0) = elevation.sin() * azimuth.cos();
		pWorkVar->coordsRot.posEa.col(1) = elevation.sin() * azimuth.sin();
		pWorkVar->coordsRot.posEa.col(2) = elevation.cos();

		// all moments lie on their easy axes
		pWorkVar->coordsRot.posMm = pWorkVar->coordsRot.posEa;
		break;
	}	
	case 3:
	{
		pWorkVar->coordsRot.posEa.col(0).setZero(); //x-direction
		pWorkVar->coordsRot.posEa.col(1).setZero(); //y-direction
		pWorkVar->coordsRot.posEa.col(2).setOnes(); //z-direction

		// all moments lie on their easy axes
		pWorkVar->coordsRot.posMm = pWorkVar->coordsRot.posEa;
		break;
	}
	}
}

// Brief: This function sets the initial particle positions from last simulation step
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Params_S* pParams - pointer to params_S object
// return: void
void InitCoordsFromSim(WorkingVar_S* pWorkVar, Params_S* pParams)
{
	MatrixXd initCoords;
	initCoords = openData("coords_final.txt");

	cout << initCoords << "\n" << endl;
	cout << "the array has the dimensions " << initCoords.rows() << "x" << initCoords.cols() << "\n" << endl;

	//define size of the dynamic arrays
	pWorkVar->coordsTrans.pos.resize(numbPart, 3);
	pWorkVar->coordsTrans.vel.resize(numbPart, 3);

	pWorkVar->coordsRot.posEa.resize(numbPart, 3);
	pWorkVar->coordsRot.posMm.resize(numbPart, 3);

	pWorkVar->coordsTrans.pos.col(0) = initCoords.col(0);
	pWorkVar->coordsTrans.pos.col(1) = initCoords.col(1);
	pWorkVar->coordsTrans.pos.col(2) = initCoords.col(2);

	pWorkVar->coordsTrans.vel.col(0) = initCoords.col(3);
	pWorkVar->coordsTrans.vel.col(1) = initCoords.col(4);
	pWorkVar->coordsTrans.vel.col(2) = initCoords.col(5);

	pWorkVar->coordsRot.posMm.col(0) = initCoords.col(6);
	pWorkVar->coordsRot.posMm.col(1) = initCoords.col(7);
	pWorkVar->coordsRot.posMm.col(2) = initCoords.col(8);

	pWorkVar->coordsRot.posEa.col(0) = initCoords.col(9);
	pWorkVar->coordsRot.posEa.col(1) = initCoords.col(10);
	pWorkVar->coordsRot.posEa.col(2) = initCoords.col(11);
}

// Brief: This function initializes the integration coefficant buffers
// param[in/out]: IntCoefficant_S* pIntCoeff - pointer to a IntCoefficant_S object
// return: void
void InitIntCoeff(IntCoefficant_S* pIntCoeff)
{
	pIntCoeff->coords2.posEa.resize(numbPart, 3);
	pIntCoeff->coords2.posMm.resize(numbPart, 3);
	pIntCoeff->coords3.posEa.resize(numbPart, 3);
	pIntCoeff->coords3.posMm.resize(numbPart, 3);
	pIntCoeff->coords4.posEa.resize(numbPart, 3);
	pIntCoeff->coords4.posMm.resize(numbPart, 3);
	pIntCoeff->coords5.posEa.resize(numbPart, 3);
	pIntCoeff->coords5.posMm.resize(numbPart, 3);
	pIntCoeff->coords6.posEa.resize(numbPart, 3);
	pIntCoeff->coords6.posMm.resize(numbPart, 3);
	pIntCoeff->coords7.posEa.resize(numbPart, 3);
	pIntCoeff->coords7.posMm.resize(numbPart, 3);

	pIntCoeff->k1.knMm.resize(numbPart, 3);
	pIntCoeff->k1.knEa.resize(numbPart, 3);
	pIntCoeff->k2.knMm.resize(numbPart, 3);
	pIntCoeff->k2.knEa.resize(numbPart, 3);
	pIntCoeff->k3.knMm.resize(numbPart, 3);
	pIntCoeff->k3.knEa.resize(numbPart, 3);
	pIntCoeff->k4.knMm.resize(numbPart, 3);
	pIntCoeff->k4.knEa.resize(numbPart, 3);
	pIntCoeff->k5.knMm.resize(numbPart, 3);
	pIntCoeff->k5.knEa.resize(numbPart, 3);
	pIntCoeff->k6.knMm.resize(numbPart, 3);
	pIntCoeff->k6.knEa.resize(numbPart, 3);
	pIntCoeff->k7.knMm.resize(numbPart, 3);
	pIntCoeff->k7.knEa.resize(numbPart, 3);
}

// Brief: This function initializes the buffers
// param[in/out]: Buffer_S* pBuffer - pointer to a Buffer_S object
// return: void
void InitBuffers(Buffer_S* pBuffer)
{
	pBuffer->buffer3d.buffer1.resize(numbPart, 3);
	pBuffer->buffer3d.buffer2.resize(numbPart, 3);
	pBuffer->buffer3d.buffer3.resize(numbPart, 3);
	pBuffer->buffer3d.buffer_b.resize(numbPart, 3);

	pBuffer->buffer1d.buffer1.resize(numbPart, 1);
}

////////////////////////*simulation functions*///////////////////////////////////////

// Brief: This function performs the simulation
// return: void
void RunSimulation(void)
{
	WorkingVar_S*     pWorkVar;     // pointer to the simulation working varibales
	Buffer_S*          pBuffer;     // pointer to buffer
	IntCoefficant_S* pIntCoeff;     // pointer to the integration coefficant varibales
	OutputVar_S*    pOutputVar;     // pointer to output variables
	Params_S*          pParams;     // pointer to parameter structure
	double              tStart;     // start time to measure
	double             runTime;     // runtime in s
	double                tInt;     // current integration time in s
	double              deltaT;     // integration timestep
	double          deltaTused;     // used integration timestep from solver

	// initialize the pointers
	pWorkVar = &workVar;
	pBuffer = &buffer;
	pIntCoeff = &intCoeff;
	pOutputVar = &outputVar;
	pParams = &params;

	// use always different seed for random number generator	
	srand((unsigned int)time(NULL));

	// initialize time settings
	tInt = 0;              // set integration time to 0
	deltaT = deltaTinit;   // set init timestep
	int tCount = 0;        // progress counter
	int steps = 0;

	// start measurement
	tStart = clock();

	// simulation loop
	while (tInt <= tEnd) //dont hit end time exactly => does not matter for my applications
	{
		// displaying simulation progress in % 0,1,2,3,4,5....
		if (floor(tInt / tEnd * dataPoints) == tCount)
		{
			runTime = (clock() - tStart) / CLOCKS_PER_SEC;

			cout << "actual timestep: " << deltaT << "s" << endl;
			cout << "simulation progress: " << tCount / dataPoints * 100.0 << "%" << endl;
			cout << "simulation will end in about " << (tEnd - tInt) * (runTime / tInt) / 3600.0 << "h" << "\n" << endl;

			tCount++;

			//store coords in txt
			WriteCoords2TXT(&workVar, "coords.txt", tInt);
		}
		
		// set external magnetic field depending on time
		DefExtField(&pWorkVar->extFluxDens, tInt);

		// set thermal fluctuations
		DefThermFluct(pWorkVar, &pBuffer->buffer1d, deltaT);
		
		// compute forces
		EvalSinCos(pWorkVar, pBuffer, pParams);
		CompInteractions(pWorkVar, pBuffer, pParams);

		// check if ewald summation works correct
		//CheckEwaldSum(pWorkVar, pBuffer, pParams);
		
		// numerical integration routine for rotational movement
  		IntegrationRot(pWorkVar, pIntCoeff, pBuffer, pParams, &deltaT, &deltaTused); 

		//numerical integration for translational movement
		IntegrationTrans(pWorkVar, deltaTused);

		// Apply periodic boundary conditions
		ApplyBoundaryCondArr(&pWorkVar->coordsTrans.pos, &pBuffer->buffer3d, pParams);

		// refresh neigbour list
		RefreshNebrList(pWorkVar, pParams, deltaTused);

		//time counting
		tInt += deltaTused;
		steps++;
	}

	// stop measurment
	runTime = (clock() - tStart) / CLOCKS_PER_SEC;
	cout << "simulation progress: 100%" << endl;
	cout << "  simulationtime = " << runTime/3600. << " h" << endl;
	cout << steps << " timesteps needed" << endl;

	// store end coords in txt
	WriteCoords2TXT(&workVar, "coords.txt", tInt);
	WriteCoords2TXT(&workVar, "coords_final.txt", 0);
}

// Brief: This function sets the external field depending on time in z direction
// param[in/out]: ArrayXXd* pExtFluxDens - pointer to the external field object
// param[in]	: double tInt - current integration time in s
// return: void
void DefExtField(ArrayXXd* pExtFluxDens, double tInt)
{
	if (tInt < tMag)
	{
		pExtFluxDens->col(2).setConstant(magFluxDens);
	}
	else
	{
		pExtFluxDens->col(2).setZero();
	}
}

// Brief: This function creates random thermal fluctuations
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]	: Buffer1d_S* pBuffer1d - pointer to Buffer1d_S object
// param[in]	: double actual timestep
// return: void
void DefThermFluct(WorkingVar_S* pWorkVar, Buffer1d_S* pBuffer1d, double deltaT)
{
		if (ENABLE_THERM_TORQUE)
		{
			pBuffer1d->buffer1 = ((2.0 * kB * temp / deltaT) * pWorkVar->partProp.zetaRot).cwiseSqrt();

			// normal distribution with mean = 0, stdev = 1.0
			pWorkVar->thermTorque = Rand::normal<ArrayXXd>(numbPart, 3, generator, 0, 1.0);

			//multiplication of nx3 matrix with nx1 vector / each column of matrix elementwise with vector
			pWorkVar->thermTorque = pWorkVar->thermTorque.colwise() * pBuffer1d->buffer1;  //works like in matlab	
		}
		else
		{
			pWorkVar->thermTorque.setZero();
		}

		if (ENABLE_THERM_FIELD)
		{
			pBuffer1d->buffer1 = (2.0 * kB * temp * magDamp / (gyroMr * satMag * deltaT) * pWorkVar->partProp.volMag.cwiseInverse()).cwiseSqrt();

			pWorkVar->thermField = Rand::normal<ArrayXXd>(numbPart, 3, generator, 0, 1.0);
			pWorkVar->thermField = pWorkVar->thermField.colwise() * pBuffer1d->buffer1;
		}
		else
		{
			pWorkVar->thermField.setZero();
		}

		if (ENABLE_THERM_FORCE)
		{
			pBuffer1d->buffer1 = ((2.0 * kB * temp / deltaT) * pWorkVar->partProp.zetaTrans).cwiseSqrt();

			pWorkVar->thermForce = Rand::normal<ArrayXXd>(numbPart, 3, generator, 0, 1.0);
			pWorkVar->thermForce = pWorkVar->thermForce.colwise() * pBuffer1d->buffer1;
		}
		else
		{
			pWorkVar->thermForce.setZero();
		}
}

// Brief: This function applies periodic boundary conditions
// param[in/out]: ArrayXXd* pArr - pointer to an array of coords to check
// param[in]	: Buffer3d_S* pBuffer3d - pointer to Buffer3d_S object
// param[in]: Params_S* pParams - pointer to params_S object
// return: void
void ApplyBoundaryCondArr(ArrayXXd* pArr, Buffer3d_S* pBuffer3d, Params_S* pParams)
{
	//use buffer_b ... boolean as idx array
	//if particle position is greater than lbox/2 => moove to other side of box
	pBuffer3d->buffer_b = *pArr > pParams->lBox / 2;
	*pArr -= pParams->lBox * pBuffer3d->buffer_b.cast<double>(); //multiply arrays elementwise, cast to convert bool to double

	//if particle position is smaller than -lbox/2 => moove to other side of box
	pBuffer3d->buffer_b = *pArr < -pParams->lBox / 2;
	*pArr += pParams->lBox * pBuffer3d->buffer_b.cast<double>();
}

// Brief: This function computes the forces acting on a particle
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]	: Buffer_S* pBuffer - pointer to Buffer_S object
// param[in]: Params_S* pParams - pointer to params_S object
// return: void
void CompInteractions(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams)
{
	pWorkVar->intForce.setZero();         //force array, initialized to zero
	pWorkVar->demagFluxDens.setZero();    //demagnetizing array, initialized to zero

	pWorkVar->partProp.magMomVec = pWorkVar->coordsRot.posMm.colwise() * pWorkVar->partProp.magMom; //create magnetic moment vectors
	pWorkVar->partProp.magMomSum = pWorkVar->partProp.magMomVec.colwise().sum(); // term for surface contribution in Ewald summation

	// compute interaction vectors and create lists for long and short range interactions
	CompInteractVecsLR(pWorkVar, pBuffer, pParams);

	//Ewald summation
	CompDipInteractEwaldR(pWorkVar, pBuffer, pParams);
	CompDipInteractEwaldF(pWorkVar, pBuffer, pParams);

	// loop for short range interactions
	int idx1, idx2;
	for (int i = 0; i < pWorkVar->nebrListData.pairsSR; i++)
	{
		idx1 = pWorkVar->nebrListData.nebrList(2 * i);
		idx2 = pWorkVar->nebrListData.nebrList(2 * i + 1);

		CompInteractVecNL(pWorkVar, pBuffer, pParams, idx1, idx2, i);

		// only executes interaction computation if particles are within interaction range
		if (pWorkVar->nebrListData.intRangeSR(i) == true)
		{
			if (ENABLE_COATING_POT)
			{
				//check if core surfaces are too close + manually correct it to avoid singularity of vdw attraction
				CheckStability(pWorkVar, idx1, idx2, i);

				//steric repulsion
				CompSterRep(pWorkVar, pBuffer, idx1, idx2, i);

				//electrostatic repulsion
				CompElectrostatRep(pWorkVar, pBuffer, idx1, idx2, i);

				// vdw forces
				CompForceVdW(pWorkVar, pBuffer, idx1, idx2, i);
			}
			else
			{
				//check and correct overlapping of particles manually
				HardSphereRep(pWorkVar, idx1, idx2, i);
			}
		}
	}

	//add thermal fluctations to force array
	pWorkVar->intForce += pWorkVar->thermForce; 
}

// Brief: This function computes the demagnetizing field
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]	: Buffer_S* pBuffer - pointer to Buffer_S object
// param[in]    : Params_S* pParams - pointer to params_S object
// param[in]    : CoordsRot_S* pCoordsRot - pointer to CoordsRot object - intermediate coords from solver
// return: void
void CompDemagField(WorkingVar_S * pWorkVar, Buffer_S * pBuffer, Params_S * pParams, CoordsRot_S * pCoordsRot)
{
	pWorkVar->demagFluxDens.setZero();    //demagnetizing array, set to zero
	pWorkVar->partProp.magMomVec = pCoordsRot->posMm.colwise() * pWorkVar->partProp.magMom; //create magnetic moment vectors with intermediate coords
	pWorkVar->partProp.magMomSum = pWorkVar->partProp.magMomVec.colwise().sum(); // term for surface contribution in Ewald summation

	CompDipFieldEwaldR(pWorkVar, pBuffer, pParams);
	CompDipFieldEwaldF(pWorkVar, pBuffer, pParams);
}

// Brief: This function performs the numerical integration with the selected integration method
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out]: IntCoefficant_S* pIntCoeff - pointer to a IntCoefficant_S object
// param[in]: Buffer_S* pBuffer - pointer to Buffer_S* object 
// param[in]	: double* pTDelta - timestep for integration (adjusted by "AdjustTimestep")
// param[out]	: double* pTDeltaUsed - pointer to the actual used timestep by the integration routine
// return: void
void IntegrationRot(WorkingVar_S* pWorkVar, IntCoefficant_S* pIntCoeff, Buffer_S* pBuffer, Params_S* pParams, double* pDeltaT, double* pDeltaTused)
{
	bool    cond;       // condition for while loop, set by AdjustTimestep
	bool    noFailed;   // condition for while loop, set by AdjustTimestep
	double  errorMm;    // error Mm
	double 	deltaT;		
	double  power;      // error convergence

	// set initial parameters
	cond = true;
	noFailed = true;

	switch (SOLVER)
	{
	case HEUN_EULER:
	{
		power = 1.0 / 2.0;

		while (cond)
		{		
			// pre load deltaT for faster execution
			deltaT = (*pDeltaT);
			RKCoeffs(pWorkVar, &pWorkVar->coordsRot, &pIntCoeff->k1, pBuffer);

			pIntCoeff->coords2.posEa = pWorkVar->coordsRot.posEa + deltaT * pIntCoeff->k1.knEa;
			pIntCoeff->coords2.posMm = pWorkVar->coordsRot.posMm + deltaT * pIntCoeff->k1.knMm;
			
			CompDemagField(pWorkVar, pBuffer, pParams, &pIntCoeff->coords2);
			RKCoeffs(pWorkVar, &pIntCoeff->coords2, &pIntCoeff->k2, pBuffer);

			// error calculation
			pBuffer->buffer3d.buffer1 = 0.5 * pIntCoeff->k1.knMm - 0.5 * pIntCoeff->k2.knMm; //error of solver
			errorMm = deltaT * (pBuffer->buffer3d.buffer1.rowwise().norm()).maxCoeff(); //same as normcontrol in matlab

			// save old deltaT and adjust timestep
			*pDeltaTused = deltaT;		
			AdjustTimeStep(pWorkVar, errorMm, power, pDeltaT, &cond, &noFailed);			
		}

		pWorkVar->coordsRot.posEa += *pDeltaTused * (0.5 * pIntCoeff->k1.knEa + 0.5 * pIntCoeff->k2.knEa); //2.order solution
		pWorkVar->coordsRot.posMm += *pDeltaTused * (0.5 * pIntCoeff->k1.knMm + 0.5 * pIntCoeff->k2.knMm);

		break;
	}
	case BOG_SHAMP:
	{
		power = 1.0 / 3.0;

		while (cond)
		{
			// pre load deltaT for faster execution
			deltaT = (*pDeltaT);

			RKCoeffs(pWorkVar, &pWorkVar->coordsRot, &pIntCoeff->k1, pBuffer);

			pIntCoeff->coords2.posEa = pWorkVar->coordsRot.posEa + deltaT * 1.0 / 2.0 * pIntCoeff->k1.knEa;
			pIntCoeff->coords2.posMm = pWorkVar->coordsRot.posMm + deltaT * 1.0 / 2.0 * pIntCoeff->k1.knMm;

			// calculate demag field based on intermediate coords
			CompDemagField(pWorkVar, pBuffer, pParams, &pIntCoeff->coords2);
			RKCoeffs(pWorkVar, &pIntCoeff->coords2, &pIntCoeff->k2, pBuffer);

			pIntCoeff->coords3.posEa = pWorkVar->coordsRot.posEa + deltaT *  3.0 / 4.0 * pIntCoeff->k2.knEa;
			pIntCoeff->coords3.posMm = pWorkVar->coordsRot.posMm + deltaT *  3.0 / 4.0 * pIntCoeff->k2.knMm;

			// calculate demag field based on intermediate coords
			CompDemagField(pWorkVar, pBuffer, pParams, &pIntCoeff->coords3);
			RKCoeffs(pWorkVar, &pIntCoeff->coords3, &pIntCoeff->k3, pBuffer);

			pIntCoeff->coords4.posEa = pWorkVar->coordsRot.posEa + deltaT * (2.0 / 9.0 * pIntCoeff->k1.knEa + 1.0 / 3.0 * pIntCoeff->k2.knEa + 4.0 / 9.0 * pIntCoeff->k3.knEa);
			pIntCoeff->coords4.posMm = pWorkVar->coordsRot.posMm + deltaT * (2.0 / 9.0 * pIntCoeff->k1.knMm + 1.0 / 3.0 * pIntCoeff->k2.knMm + 4.0 / 9.0 * pIntCoeff->k3.knMm);	

			// calculate demag field based on intermediate coords
			CompDemagField(pWorkVar, pBuffer, pParams, &pIntCoeff->coords4);
			RKCoeffs(pWorkVar, &pIntCoeff->coords4, &pIntCoeff->k4, pBuffer);

			pBuffer->buffer3d.buffer1 = -5.0 / 72.0 * pIntCoeff->k1.knMm + 1.0 / 12.0 * pIntCoeff->k2.knMm + 1.0 / 9.0 * pIntCoeff->k3.knMm - 1.0 / 8.0 * pIntCoeff->k4.knMm;
			errorMm = deltaT * (pBuffer->buffer3d.buffer1.rowwise().norm()).maxCoeff();

			// save old deltaT and adjust timestep
			*pDeltaTused = deltaT;
			AdjustTimeStep(pWorkVar, errorMm, power, pDeltaT, &cond, &noFailed);
		}

		pWorkVar->coordsRot.posEa = pIntCoeff->coords4.posEa;
		pWorkVar->coordsRot.posMm = pIntCoeff->coords4.posMm;

		break;
	}
	case DOR_PRI:
	{
		power = 1.0 / 5.0;

		while (cond)
		{
			// pre load deltaT for faster execution
			deltaT = (*pDeltaT);

			RKCoeffs(pWorkVar, &pWorkVar->coordsRot, &pIntCoeff->k1, pBuffer);

			pIntCoeff->coords2.posEa = pWorkVar->coordsRot.posEa + deltaT * 1.0 / 5.0 * pIntCoeff->k1.knEa;
			pIntCoeff->coords2.posMm = pWorkVar->coordsRot.posMm + deltaT * 1.0 / 5.0 * pIntCoeff->k1.knMm;

			// calculate demag field based on intermediate coords
			CompDemagField(pWorkVar, pBuffer, pParams, &pIntCoeff->coords2);
			RKCoeffs(pWorkVar, &pIntCoeff->coords2, &pIntCoeff->k2, pBuffer);

			pIntCoeff->coords3.posEa = pWorkVar->coordsRot.posEa + deltaT * (3.0 / 40.0 * pIntCoeff->k1.knEa + 9.0 / 40.0 * pIntCoeff->k2.knEa);
			pIntCoeff->coords3.posMm = pWorkVar->coordsRot.posMm + deltaT * (3.0 / 40.0 * pIntCoeff->k1.knMm + 9.0 / 40.0 * pIntCoeff->k2.knMm);

			// calculate demag field based on intermediate coords
			CompDemagField(pWorkVar, pBuffer, pParams, &pIntCoeff->coords3);
			RKCoeffs(pWorkVar, &pIntCoeff->coords3, &pIntCoeff->k3, pBuffer);

			pIntCoeff->coords4.posEa = pWorkVar->coordsRot.posEa + deltaT * (44.0 / 45.0 * pIntCoeff->k1.knEa - 56.0 / 15.0 * pIntCoeff->k2.knEa + 32.0 / 9.0 * pIntCoeff->k3.knEa);
			pIntCoeff->coords4.posMm = pWorkVar->coordsRot.posMm + deltaT * (44.0 / 45.0 * pIntCoeff->k1.knMm - 56.0 / 15.0 * pIntCoeff->k2.knMm + 32.0 / 9.0 * pIntCoeff->k3.knMm);

			// calculate demag field based on intermediate coords
			CompDemagField(pWorkVar, pBuffer, pParams, &pIntCoeff->coords4);
			RKCoeffs(pWorkVar, &pIntCoeff->coords4, &pIntCoeff->k4, pBuffer);

			pIntCoeff->coords5.posEa = pWorkVar->coordsRot.posEa + deltaT * (19372.0 / 6561.0 * pIntCoeff->k1.knEa - 25360.0 / 2187.0 * pIntCoeff->k2.knEa + 64448.0 / 6561.0 * pIntCoeff->k3.knEa - 212.0 / 729.0 * pIntCoeff->k4.knEa);
			pIntCoeff->coords5.posMm = pWorkVar->coordsRot.posMm + deltaT * (19372.0 / 6561.0 * pIntCoeff->k1.knMm - 25360.0 / 2187.0 * pIntCoeff->k2.knMm + 64448.0 / 6561.0 * pIntCoeff->k3.knMm - 212.0 / 729.0 * pIntCoeff->k4.knMm);

			// calculate demag field based on intermediate coords
			CompDemagField(pWorkVar, pBuffer, pParams, &pIntCoeff->coords5);
			RKCoeffs(pWorkVar, &pIntCoeff->coords5, &pIntCoeff->k5, pBuffer);

			pIntCoeff->coords6.posEa = pWorkVar->coordsRot.posEa + deltaT * (9017.0 / 3168.0 * pIntCoeff->k1.knEa - 355.0 / 33.0 * pIntCoeff->k2.knEa + 46732.0 / 5247.0 * pIntCoeff->k3.knEa + 49.0 / 176.0 * pIntCoeff->k4.knEa - 5103.0 / 18656.0 * pIntCoeff->k5.knEa);
			pIntCoeff->coords6.posMm = pWorkVar->coordsRot.posMm + deltaT * (9017.0 / 3168.0 * pIntCoeff->k1.knMm - 355.0 / 33.0 * pIntCoeff->k2.knMm + 46732.0 / 5247.0 * pIntCoeff->k3.knMm + 49.0 / 176.0 * pIntCoeff->k4.knMm - 5103.0 / 18656.0 * pIntCoeff->k5.knMm);

			// calculate demag field based on intermediate coords
			CompDemagField(pWorkVar, pBuffer, pParams, &pIntCoeff->coords6);
			RKCoeffs(pWorkVar, &pIntCoeff->coords6, &pIntCoeff->k6, pBuffer);

			pIntCoeff->coords7.posEa = pWorkVar->coordsRot.posEa + deltaT * (35.0 / 384.0 * pIntCoeff->k1.knEa + 500.0 / 1113.0 * pIntCoeff->k3.knEa + 125.0 / 192.0 * pIntCoeff->k4.knEa - 2187.0 / 6784.0 * pIntCoeff->k5.knEa + 11.0 / 84.0 * pIntCoeff->k6.knEa);
			pIntCoeff->coords7.posMm = pWorkVar->coordsRot.posMm + deltaT * (35.0 / 384.0 * pIntCoeff->k1.knMm + 500.0 / 1113.0 * pIntCoeff->k3.knMm + 125.0 / 192.0 * pIntCoeff->k4.knMm - 2187.0 / 6784.0 * pIntCoeff->k5.knMm + 11.0 / 84.0 * pIntCoeff->k6.knMm);

			// calculate demag field based on intermediate coords
			CompDemagField(pWorkVar, pBuffer, pParams, &pIntCoeff->coords7);
			RKCoeffs(pWorkVar, &pIntCoeff->coords7, &pIntCoeff->k7, pBuffer);

			// error calculation
			pBuffer->buffer3d.buffer1 = 71.0 / 57600.0 * pIntCoeff->k1.knMm - 71.0 / 16695.0 * pIntCoeff->k3.knMm + 71.0 / 1920.0 * pIntCoeff->k4.knMm - 17253.0 / 339200.0 * pIntCoeff->k5.knMm + 22.0 / 525.0 * pIntCoeff->k6.knMm - 1.0 / 40.0 * pIntCoeff->k7.knMm;
			errorMm = deltaT * (pBuffer->buffer3d.buffer1.rowwise().norm()).maxCoeff();
			
			// save old deltaT and adjust timestep
			*pDeltaTused = deltaT;
			AdjustTimeStep(pWorkVar, errorMm, power, pDeltaT, &cond, &noFailed);
		}

		pWorkVar->coordsRot.posEa += *pDeltaTused * (35.0 / 384.0 * pIntCoeff->k1.knEa + 500.0 / 1113.0 * pIntCoeff->k3.knEa + 125.0 / 192.0 * pIntCoeff->k4.knEa - 2187.0 / 6784 * pIntCoeff->k5.knEa + 11.0 / 84.0 * pIntCoeff->k6.knEa); //5.order solution
		pWorkVar->coordsRot.posMm += *pDeltaTused * (35.0 / 384.0 * pIntCoeff->k1.knMm + 500.0 / 1113.0 * pIntCoeff->k3.knMm + 125.0 / 192.0 * pIntCoeff->k4.knMm - 2187.0 / 6784 * pIntCoeff->k5.knMm + 11.0 / 84.0 * pIntCoeff->k6.knMm);

		break;
	}
	}

	// make sure that the length of the unit vectors stay 1
	pWorkVar->coordsRot.posMm.colwise() /= pWorkVar->coordsRot.posMm.rowwise().norm();
	pWorkVar->coordsRot.posEa.colwise() /= pWorkVar->coordsRot.posEa.rowwise().norm();
}

// Brief: Calculate the runge-kutta coefficients for magnetic moment and easy axes movement
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out]: Coord_S* pCoord - pointer to the current Coord_S object
// param[in/out]: RkCoeff_S* pRkCoeffs - pointer to the current RK coefficants
// param[in]: Buffer_S* pBuffer - pointer to Buffer_S* object 
// return: void
void RKCoeffs(WorkingVar_S* pWorkVar, CoordsRot_S* pCoord, RkCoeff_S* pRkCoeffs, Buffer_S* pBuffer)
{
	// matlab code: fluxDensEff = extFluxDens + anisConst*dot(posMm,posEa).*posEa + fieldTherm;

	pBuffer->buffer1d.buffer1 = (pCoord->posMm * pCoord->posEa).rowwise().sum(); //dot Product dot(posMm,posEa)
	pBuffer->buffer3d.buffer1 = pWorkVar->extFluxDens + pWorkVar->demagFluxDens + anisConst * (pCoord->posEa.colwise() * pBuffer->buffer1d.buffer1) + pWorkVar->thermField; //effektive flux density

	// check if particles can rotate free
	if (true == ENABLE_MOBILIZATION)
	{
		//matlab code knEa = velEaConst1.*dot(posMm,posEa).*(posMm - dot(posMm,posEa).*posEa) + velEaConst2.*cross(torqueTherm, posEa);

		RowWiseCrossProd(&pWorkVar->thermTorque, &pCoord->posEa, &pBuffer->buffer3d.buffer2);

		pRkCoeffs->knEa = (pCoord->posMm - pCoord->posEa.colwise() * pBuffer->buffer1d.buffer1).colwise() * (pWorkVar->partProp.velEaConst1 * pBuffer->buffer1d.buffer1) + pBuffer->buffer3d.buffer2.colwise() * pWorkVar->partProp.velEaConst2;
	}
	else
	{
		pRkCoeffs->knEa.setZero();
	}

	//velMm = -velMmConst*(cross(posMm,fluxDensEff) + magDamp*cross(posMm,cross(posMm,fluxDensEff)));
	RowWiseCrossProd(&pCoord->posMm, &pBuffer->buffer3d.buffer1, &pBuffer->buffer3d.buffer2);
	RowWiseCrossProd(&pCoord->posMm, &pBuffer->buffer3d.buffer2, &pBuffer->buffer3d.buffer3);

	pRkCoeffs->knMm = -velMmConst * (pBuffer->buffer3d.buffer2 + magDamp * pBuffer->buffer3d.buffer3);
}

// Brief: Adjusting timeStep to meet tolerance requirements 
// param[in] : WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : double error - error value
// param[in] : double power - power for timestep adjustment
// param[out]: double* pTDelta - new time delate
// param[out]: bool* pCond - new loop condidtion
// param[out]: bool* pNoFailed - new failed condidtion
// return: void
void AdjustTimeStep(WorkingVar_S* pWorkVar, double error, double power, double* pDeltaT, bool* pCond, bool* pNoFailed)
{
	double tDeltaNew;

	// norm(e(i)) <= max(RelTol*norm(y(i)),AbsTol(i)) from matlab, norm(y(i)) is 1 in the case of unit vectors
	if (error > absTol) //the following code is from Matlab Ode45 solver
	{
		if (*pNoFailed == true)
		{
			*pNoFailed = false;
			tDeltaNew = max(deltaTmin, (*pDeltaT) * max(0.1, 0.9 * pow(absTol / error, power)));
		}
		else
		{
			tDeltaNew = max(deltaTmin, ((*pDeltaT) * 0.5));
		}

		pWorkVar->thermField *= sqrt(tDeltaNew / (*pDeltaT));
		pWorkVar->thermTorque *= sqrt(tDeltaNew / (*pDeltaT));
		pWorkVar->thermForce *= sqrt(tDeltaNew / (*pDeltaT));
	}
	else
	{
		if (*pNoFailed == true)
		{
			tDeltaNew = (*pDeltaT) * min(max(0.9 * pow(absTol / error, power), 0.2), 5.0);
		}
		else
		{
			tDeltaNew = (*pDeltaT);
		}

		*pCond = false;
	}

	//make sure the timestep stays below tDeltaMax
	*pDeltaT = min(tDeltaNew, deltaTmax);
}

// Brief: Integration of translational movement
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: deltaT timestep used in the rotational integration
// return: void
void IntegrationTrans(WorkingVar_S* pWorkVar, double deltaTused)
{
	pWorkVar->coordsTrans.pos += deltaTused * (pWorkVar->intForce.colwise() / pWorkVar->partProp.zetaTrans);
}