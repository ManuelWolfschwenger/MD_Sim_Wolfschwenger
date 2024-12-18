/*
* Name : simulation.cpp
* Auto :
* Brief: Simulation source file
*/

// own header file
#include "simulation.h"

// global variables
static WorkingVar_S     workVar;    // simulation working varibales
static IntCoefficient_S  intCoeff;   // integration coefficant varibales
static OutputVar_S      outputVar;  // output variables
Configuration* config = Configuration::getInstance(); // load configuration from configuration.ini

// random number generator
Rand::Vmt19937_64 generator{ (uint64_t)(config->getSeed()) };

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

// Brief: This function initializes the simulation and the simulation arrays
// return: void
void InitSimulation(void)
{
	WorkingVar_S* pWorkVar;      // pointer to the simulation working variables
	IntCoefficient_S* pIntCoeff;  // pointer to the integration coefficant variables

	// initialize the pointers
	pWorkVar = &workVar;
	pIntCoeff = &intCoeff;

	// initialize the particle properties for this simulation
	InitParticleProps(&pWorkVar->partProp);

	// initialize the particle starting rotational positions 
	InitCoords(&pWorkVar->coords);
	
	// initialize the integration coefficant buffers
	InitIntCoeff(pIntCoeff);

	// initialize the buffers
	InitBuffers(&pWorkVar->buffer3d, &pWorkVar->buffer1d);

	// initialize the external field 
	pWorkVar->extFluxDens.resize(config->getNumbPart(), 3); 
	pWorkVar->extFluxDens.setZero();                     //in case we start with a relaxation, dont need to set to zero every iteration in the loop

	// initialize the thermal fluctuations
	pWorkVar->thermTorque.resize(config->getNumbPart(), 3);
	pWorkVar->thermField.resize(config->getNumbPart(), 3);
	pWorkVar->thermTorque.setZero();
	pWorkVar->thermField.setZero();

	// store data in txt
	WriteData2TXT(pWorkVar, "data.txt");
}

// Brief: This function sets the particle properties like size, volume, magnetic moment...
// param[in/out]: PartProps_S* pPartProps - pointer to a PartProps_S object
// return: void
void InitParticleProps(PartProps_S* pPartProps)
{
	ArrayXd dHydr(config->getNumbPart(), 1); //hydrodynamic diameters of particles

	// switch depending on the size distribution
	switch (config->getSizeDist())
	{
	case SIZE_DIST_EQUAL:
	{
		// set all particle diameters to same value
		dHydr.setConstant(2.0 * config->getRHydrMean());
	}
	break;

	case SIZE_DIST_LOG:
	{
		// log-normal distributed sizes
		double dHydrMean = 2.0 * config->getRHydrMean();
		double mu = log(pow(dHydrMean, 2.0) / sqrt(pow(dHydrMean, 2.0) + pow(config->getSigma(), 2.0)));
		double sigmaLog = sqrt(log(1.0 + pow(config->getSigma() / dHydrMean, 2.0)));

		// lognormal distribution with mean = 0, stdev = 1.0
		dHydr = Rand::lognormal<MatrixXd>(config->getNumbPart(), 1, generator, mu, sigmaLog);
	}
	break;
	}

	const double rRatio = config->getRMagMean() / config->getRHydrMean();  // ratio to scale the magnetic radius depending on the distribution of the hydrodynamic radius

	// assign values to structure the PartProps_S object
	pPartProps->rHydr = dHydr / 2.0;                                                                         // hydrodyn radius
	pPartProps->rMag = pPartProps->rHydr * rRatio;                                                           // mag radius
	pPartProps->volMag = 4.0 / 3.0 * config->getMyPi() * pPartProps->rMag.pow(3.0);                          // mag volume
	pPartProps->volHydr = 4.0 / 3.0 * config->getMyPi() * pPartProps->rHydr.pow(3.0);                        // hydr. volume
	pPartProps->magMom = config->getSatMag() * pPartProps->volMag;                                           // mag Moment
	pPartProps->zetaRot = 8.0 * config->getMyPi() * config->getVis() * pPartProps->rHydr.pow(3.0);           // rotational friction coefficient
	pPartProps->velEaConst1 = 2.0 * config->getAnisEn() * (pPartProps->volMag / pPartProps->zetaRot);        // constant for later DGL solving
	pPartProps->velEaConst2 = pPartProps->zetaRot.cwiseInverse();                                            // constant for later DGL solving
}

// Brief: This function sets the initial particle rotational positions
// param[in/out]: Coord_S* pCoords - pointer to a Coord_S object
// return: void
void InitCoords(Coord_S* pCoords)
{
	//azimuth and elevation angle
	ArrayXd azimuth(config->getNumbPart());
	ArrayXd elevation(config->getNumbPart());

	//define size of the dynamic arrays
	pCoords->posEa.resize(config->getNumbPart(), 3);
	pCoords->posMm.resize(config->getNumbPart(), 3);

	// create uniformly random distributed vector components in a sphere => see also in my paper
	azimuth = Rand::balanced<ArrayXXd>(config->getNumbPart(), 1, generator); //random numbers between -1 and 1
	azimuth += 1; //random numbers between 0 and 2
	azimuth *= 0.5; //random numbers between 0 and 1
	azimuth *= 2 * config->getMyPi();//uniformly distributed random numbers between 0 and 2*pi

	elevation = Rand::balanced<ArrayXXd>(config->getNumbPart(), 1, generator); //random numbers between -1 and 1
	elevation += 1; //random numbers between 0 and 2
	elevation *= 0.5; //random numbers between 0 and 1
	elevation = elevation.cwiseSqrt().asin() * 2.0;

	// switch depending on the initial configuration
	switch (config->getInitConfigRot())
	{
	case 1:
	{
		pCoords->posEa.col(0) = elevation.sin() * azimuth.cos();
		pCoords->posEa.col(1) = elevation.sin() * azimuth.sin();
		pCoords->posEa.col(2) = elevation.cos();

		// all moments point in z direction / easy axis random see above
		pCoords->posMm.col(0).setZero(); //x-direction
		pCoords->posMm.col(1).setZero(); //y-direction
		pCoords->posMm.col(2).setOnes(); //z-direction
		break;
	}
	case 2:
	{
		pCoords->posEa.col(0) = elevation.sin() * azimuth.cos();
		pCoords->posEa.col(1) = elevation.sin() * azimuth.sin();
		pCoords->posEa.col(2) = elevation.cos();

		// all moments lie on their easy axes, / easy axis random see above
		pCoords->posMm = pCoords->posEa;
		break;
	}
	case 3:
	{
		// easy axis in z-direction
		pCoords->posEa.col(0).setZero(); //x-direction
		pCoords->posEa.col(1).setZero(); //y-direction
		pCoords->posEa.col(2).setOnes(); //z-direction

		// all moments lie on their easy axes
		pCoords->posMm = pCoords->posEa;
		break;
	}
	}
}

// Brief: This function initializes the integration coefficient buffers
// param[in/out]: IntCoefficient_S* pIntCoeff - pointer to a IntCoefficient_S object
// return: void
void InitIntCoeff(IntCoefficient_S* pIntCoeff)
{
	// buffers for Dormand Price Method (mostly used), for Heun-Euler and Bogacki Shampine some remain empty afterwards for lower order solvers
	pIntCoeff->coords2.posEa.resize(config->getNumbPart(), 3);
	pIntCoeff->coords2.posMm.resize(config->getNumbPart(), 3);
	pIntCoeff->coords3.posEa.resize(config->getNumbPart(), 3);
	pIntCoeff->coords3.posMm.resize(config->getNumbPart(), 3);
	pIntCoeff->coords4.posEa.resize(config->getNumbPart(), 3);
	pIntCoeff->coords4.posMm.resize(config->getNumbPart(), 3);
	pIntCoeff->coords5.posEa.resize(config->getNumbPart(), 3);
	pIntCoeff->coords5.posMm.resize(config->getNumbPart(), 3);
	pIntCoeff->coords6.posEa.resize(config->getNumbPart(), 3);
	pIntCoeff->coords6.posMm.resize(config->getNumbPart(), 3);
	pIntCoeff->coords7.posEa.resize(config->getNumbPart(), 3);
	pIntCoeff->coords7.posMm.resize(config->getNumbPart(), 3);

	pIntCoeff->k1.knMm.resize(config->getNumbPart(), 3);
	pIntCoeff->k1.knEa.resize(config->getNumbPart(), 3);
	pIntCoeff->k2.knMm.resize(config->getNumbPart(), 3);
	pIntCoeff->k2.knEa.resize(config->getNumbPart(), 3);
	pIntCoeff->k3.knMm.resize(config->getNumbPart(), 3);
	pIntCoeff->k3.knEa.resize(config->getNumbPart(), 3);
	pIntCoeff->k4.knMm.resize(config->getNumbPart(), 3);
	pIntCoeff->k4.knEa.resize(config->getNumbPart(), 3);
	pIntCoeff->k5.knMm.resize(config->getNumbPart(), 3);
	pIntCoeff->k5.knEa.resize(config->getNumbPart(), 3);
	pIntCoeff->k6.knMm.resize(config->getNumbPart(), 3);
	pIntCoeff->k6.knEa.resize(config->getNumbPart(), 3);
	pIntCoeff->k7.knMm.resize(config->getNumbPart(), 3);
	pIntCoeff->k7.knEa.resize(config->getNumbPart(), 3);
}

// Brief: This function initializes the buffers
// param[in/out]: Buffer3d_S* pBuffer3d - pointer to a Buffer3d_S object
// param[in/out]: Buffer1d_S* pBuffer1d - pointer to a Buffer1d_S object
// return: void
void InitBuffers(Buffer3d_S* pBuffer3d, Buffer1d_S* pBuffer1d)
{
	pBuffer3d->buffer1.resize(config->getNumbPart(), 3);
	pBuffer3d->buffer2.resize(config->getNumbPart(), 3);
	pBuffer3d->buffer3.resize(config->getNumbPart(), 3);

	pBuffer1d->buffer1.resize(config->getNumbPart(), 1);
}

////////////////////////*simulation functions*///////////////////////////////////////

// Brief: This function performs the simulation
// return: void
void RunSimulation(void)
{
	WorkingVar_S*       pWorkVar;   // pointer to the simulation working varibales
	IntCoefficient_S*   pIntCoeff;  // pointer to the integration coefficant varibales
	OutputVar_S*        pOutputVar; // pointer to output variables
	double              tStart;     // start timestamp
	double              runTime;    // runtime in s
	double              tInt;       // current integration time in s
	double              deltaT;     // timestep

	// initialize the pointers
	pWorkVar = &workVar;
	pIntCoeff = &intCoeff;
	pOutputVar = &outputVar;

	// initialize timesettings
	tInt = 0;
	deltaT = config->getDeltaTInit();
	int tCount = 0;
	int tCountData = 0;
	int steps = 0;

	// start measurment
	tStart = clock();
	
	// simulation loop
	while (tInt <= config->getTEnd())
	{
		if (floor(tInt / config->getTEnd() * 100) == tCount)
		{
			runTime = (clock() - tStart) / CLOCKS_PER_SEC;

			cout << "actual timestep: " << deltaT << "s" << endl;
			cout << "actual time: " << tInt << "s" << endl;
			cout << "simulation progress: " << tCount << "%" << endl;
			cout << "simulation will end in about " << (config->getTEnd() - tInt) * (runTime / tInt) / 3600.0 << "h" << "\n" << endl;

			tCount++;
		}

		// displaying simulation progress in % 0,1,2,3,4,5....
		if (floor(tInt / config->getTEnd() * config->getDataPoints()) == tCountData)
		{
			//magnetization vector
			pOutputVar->mZ.push_back(pWorkVar->coords.posMm.col(2).mean());
			pOutputVar->nZ.push_back(pWorkVar->coords.posEa.col(2).abs().mean());
			pOutputVar->t.push_back(tInt);

			//store coords in txt
			//WriteCoords2TXT(&workVar, "coords.txt", tInt);

			tCountData++;
		}

		// set external magnetic field depending on time
		DefExtField(&pWorkVar->extFluxDens, tInt);

		// set thermal fluctuations
		DefThermFluct(pWorkVar, deltaT);
		
		// numerical integration routine / tInt is set inside
		Integration(pWorkVar, pIntCoeff, &deltaT, &tInt);

		steps++;
	}

	// stop measurment
	runTime = (clock() - tStart) / CLOCKS_PER_SEC;
	cout << "simulation progress: 100%" << endl;
	cout << "  simulationtime = " << runTime / 3600. << " h" << endl;
	cout << steps << " timesteps needed" << endl;

	//postprocessing
	saveData("time.txt", Map<VectorXd, Unaligned>(pOutputVar->t.data(), pOutputVar->t.size()), 0);
	saveData("mz.txt", Map<VectorXd, Unaligned>(pOutputVar->mZ.data(), pOutputVar->mZ.size()), 0);
	saveData("nz.txt", Map<VectorXd, Unaligned>(pOutputVar->nZ.data(), pOutputVar->nZ.size()), 0);
	saveData("phi.txt", Map<VectorXd, Unaligned>(pOutputVar->phi.data(), pOutputVar->phi.size()), 0);

	cout << "\n simulation ended successfuly" << endl;
}

// Brief: This function sets the external field depending on time in z direction
// param[in/out]: ArrayXXd* pExtFluxDens - pointer to the external field object
// param[in]	: double tInt - current integration time in s
// return: void
void DefExtField(ArrayXXd* pExtFluxDens, double tInt)
{
	if (tInt < config->getTMag())
	{
		pExtFluxDens->col(2).setConstant(config->getMagFluxDens());
	}
	else
	{
		pExtFluxDens->col(2).setZero();
	}
}

// Brief: This function creates random thermal fluctuations depending on the timestep
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]	: double tDelta
// return: void
void DefThermFluct(WorkingVar_S* pWorkVar, double deltaT)
{
		if (config->getEnableThermTorque())
		{
			pWorkVar->buffer1d.buffer1 = ((2.0 * config->getKB() * config->getTemp() / deltaT) * pWorkVar->partProp.zetaRot).cwiseSqrt();

			// normal distribution with mean = 0, stdev = 1.0
			pWorkVar->thermTorque = Rand::normal<ArrayXXd>(config->getNumbPart(), 3, generator, 0, 1.0);

			//multiplication of nx3 matrix with nx1 vector / each column of matrix elementwise with vector
			pWorkVar->thermTorque = pWorkVar->thermTorque.colwise() * pWorkVar->buffer1d.buffer1;  //works like in matlab	
		}
		else
		{
			pWorkVar->thermTorque.setZero();
		}

		
		if (config->getEnableThermField())
		{
			pWorkVar->buffer1d.buffer1 = (2.0 * config->getKB() * config->getTemp() * config->getMagDamp() / (config->getGyroMr() * config->getSatMag() * deltaT) * pWorkVar->partProp.volMag.cwiseInverse()).cwiseSqrt();

			pWorkVar->thermField = Rand::normal<ArrayXXd>(config->getNumbPart(), 3, generator, 0, 1.0);
			pWorkVar->thermField = pWorkVar->thermField.colwise() * pWorkVar->buffer1d.buffer1;
		}
		else
		{
			pWorkVar->thermField.setZero();
		}
}

// Brief: This function performs the numerical integration with the selected integration method
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out]: IntCoefficant_S* pIntCoeff - pointer to a IntCoefficant_S object
// param[in]	: double* pTDelta - pointer to timestep
// param[in/out]: double* pTInt - pointer to simulation time
// return: void
void Integration(WorkingVar_S* pWorkVar, IntCoefficient_S* pIntCoeff, double* pDeltaT, double* pTInt)
{
	bool   cond;       // condition for while loop, set by AdjustTimestep
	bool   noFailed;   // condition for while loop, set by AdjustTimestep
	double deltaTold;  // old timeStep
	double error;      // error
	double deltaT;	   // tDelta buffer
	double power;      // error convergence

	// set initial parameters
	cond = true;
	noFailed = true;

	switch (config->getSolver())
	{
	case 0: //Heun-Euler
	{
		power = 0.5;

		while (cond)
		{
			// pre load tDelta for faster execution
			deltaT = (*pDeltaT);

			RKCoeffs(pWorkVar, &pWorkVar->coords, &pIntCoeff->k1);

			pIntCoeff->coords2.posEa = pWorkVar->coords.posEa + deltaT * pIntCoeff->k1.knEa;
			pIntCoeff->coords2.posMm = pWorkVar->coords.posMm + deltaT * pIntCoeff->k1.knMm;

			RKCoeffs(pWorkVar, &pIntCoeff->coords2, &pIntCoeff->k2);

			pWorkVar->buffer3d.buffer2 = 0.5 * pIntCoeff->k1.knMm - 0.5 * pIntCoeff->k2.knMm;
			error = deltaT * (pWorkVar->buffer3d.buffer2.rowwise().norm()).maxCoeff();

			// save old tDelta and adjust timestep
			deltaTold = deltaT;
			AdjustTimeStep(pWorkVar, error, power, pDeltaT, &cond, &noFailed);
		}

		pWorkVar->coords.posEa += deltaTold * (0.5 * pIntCoeff->k1.knEa + 0.5 * pIntCoeff->k2.knEa); //2.order solution
		pWorkVar->coords.posMm += deltaTold * (0.5 * pIntCoeff->k1.knMm + 0.5 * pIntCoeff->k2.knMm);

		break;
	}
	case 1: //Bogacki-Shampine
	{
		power = 1.0 / 3.0;

		while (cond)
		{
			// pre load tDelta for faster execution
			deltaT = (*pDeltaT);

			RKCoeffs(pWorkVar, &pWorkVar->coords, &pIntCoeff->k1);

			pIntCoeff->coords2.posEa = pWorkVar->coords.posEa + deltaT * 1.0 / 2.0 * pIntCoeff->k1.knEa;
			pIntCoeff->coords2.posMm = pWorkVar->coords.posMm + deltaT * 1.0 / 2.0 * pIntCoeff->k1.knMm;

			RKCoeffs(pWorkVar, &pIntCoeff->coords2, &pIntCoeff->k2);

			pIntCoeff->coords3.posEa = pWorkVar->coords.posEa + deltaT * 3.0 / 4.0 * pIntCoeff->k2.knEa;
			pIntCoeff->coords3.posMm = pWorkVar->coords.posMm + deltaT * 3.0 / 4.0 * pIntCoeff->k2.knMm;

			RKCoeffs(pWorkVar, &pIntCoeff->coords3, &pIntCoeff->k3);

			pIntCoeff->coords4.posEa = pWorkVar->coords.posEa + deltaT * (2.0 / 9.0 * pIntCoeff->k1.knEa + 1.0 / 3.0 * pIntCoeff->k2.knEa + 4.0 / 9.0 * pIntCoeff->k3.knEa);
			pIntCoeff->coords4.posMm = pWorkVar->coords.posMm + deltaT * (2.0 / 9.0 * pIntCoeff->k1.knMm + 1.0 / 3.0 * pIntCoeff->k2.knMm + 4.0 / 9.0 * pIntCoeff->k3.knMm);

			RKCoeffs(pWorkVar, &pIntCoeff->coords4, &pIntCoeff->k4);

			pWorkVar->buffer3d.buffer2 = -5.0 / 72.0 * pIntCoeff->k1.knMm + 1.0 / 12.0 * pIntCoeff->k2.knMm + 1.0 / 9.0 * pIntCoeff->k3.knMm - 1.0 / 8.0 * pIntCoeff->k4.knMm;
			error = deltaT * (pWorkVar->buffer3d.buffer2.rowwise().norm()).maxCoeff();

			// save old tDelta and adjust timestamp
			deltaTold = deltaT;
			AdjustTimeStep(pWorkVar, error, power, pDeltaT, &cond, &noFailed);
		}

		pWorkVar->coords.posEa = pIntCoeff->coords4.posEa;
		pWorkVar->coords.posMm = pIntCoeff->coords4.posMm;

		break;
	}
	case 2: //Dormand Prince
	{
		power = 0.2;

		while (cond)
		{
			// pre load tDelta for faster execution
			deltaT = (*pDeltaT);

			RKCoeffs(pWorkVar, &pWorkVar->coords, &pIntCoeff->k1);

			pIntCoeff->coords2.posEa = pWorkVar->coords.posEa + deltaT * 1.0 / 5.0 * pIntCoeff->k1.knEa;
			pIntCoeff->coords2.posMm = pWorkVar->coords.posMm + deltaT * 1.0 / 5.0 * pIntCoeff->k1.knMm;

			RKCoeffs(pWorkVar, &pIntCoeff->coords2, &pIntCoeff->k2);

			pIntCoeff->coords3.posEa = pWorkVar->coords.posEa + deltaT * (3.0 / 40.0 * pIntCoeff->k1.knEa + 9.0 / 40.0 * pIntCoeff->k2.knEa);
			pIntCoeff->coords3.posMm = pWorkVar->coords.posMm + deltaT * (3.0 / 40.0 * pIntCoeff->k1.knMm + 9.0 / 40.0 * pIntCoeff->k2.knMm);

			RKCoeffs(pWorkVar, &pIntCoeff->coords3, &pIntCoeff->k3);

			pIntCoeff->coords4.posEa = pWorkVar->coords.posEa + deltaT * (44.0 / 45.0 * pIntCoeff->k1.knEa - 56.0 / 15.0 * pIntCoeff->k2.knEa + 32.0 / 9.0 * pIntCoeff->k3.knEa);
			pIntCoeff->coords4.posMm = pWorkVar->coords.posMm + deltaT * (44.0 / 45.0 * pIntCoeff->k1.knMm - 56.0 / 15.0 * pIntCoeff->k2.knMm + 32.0 / 9.0 * pIntCoeff->k3.knMm);

			RKCoeffs(pWorkVar, &pIntCoeff->coords4, &pIntCoeff->k4);

			pIntCoeff->coords5.posEa = pWorkVar->coords.posEa + deltaT * (19372.0 / 6561.0 * pIntCoeff->k1.knEa - 25360.0 / 2187.0 * pIntCoeff->k2.knEa + 64448.0 / 6561.0 * pIntCoeff->k3.knEa - 212.0 / 729.0 * pIntCoeff->k4.knEa);
			pIntCoeff->coords5.posMm = pWorkVar->coords.posMm + deltaT * (19372.0 / 6561.0 * pIntCoeff->k1.knMm - 25360.0 / 2187.0 * pIntCoeff->k2.knMm + 64448.0 / 6561.0 * pIntCoeff->k3.knMm - 212.0 / 729.0 * pIntCoeff->k4.knMm);

			RKCoeffs(pWorkVar, &pIntCoeff->coords5, &pIntCoeff->k5);

			pIntCoeff->coords6.posEa = pWorkVar->coords.posEa + deltaT * (9017.0 / 3168.0 * pIntCoeff->k1.knEa - 355.0 / 33.0 * pIntCoeff->k2.knEa + 46732.0 / 5247.0 * pIntCoeff->k3.knEa + 49.0 / 176.0 * pIntCoeff->k4.knEa - 5103.0 / 18656.0 * pIntCoeff->k5.knEa);
			pIntCoeff->coords6.posMm = pWorkVar->coords.posMm + deltaT * (9017.0 / 3168.0 * pIntCoeff->k1.knMm - 355.0 / 33.0 * pIntCoeff->k2.knMm + 46732.0 / 5247.0 * pIntCoeff->k3.knMm + 49.0 / 176.0 * pIntCoeff->k4.knMm - 5103.0 / 18656.0 * pIntCoeff->k5.knMm);

			RKCoeffs(pWorkVar, &pIntCoeff->coords6, &pIntCoeff->k6);

			pIntCoeff->coords7.posEa = pWorkVar->coords.posEa + deltaT * (35.0 / 384.0 * pIntCoeff->k1.knEa + 500.0 / 1113.0 * pIntCoeff->k3.knEa + 125.0 / 192.0 * pIntCoeff->k4.knEa - 2187.0 / 6784.0 * pIntCoeff->k5.knEa + 11.0 / 84.0 * pIntCoeff->k6.knEa);
			pIntCoeff->coords7.posMm = pWorkVar->coords.posMm + deltaT * (35.0 / 384.0 * pIntCoeff->k1.knMm + 500.0 / 1113.0 * pIntCoeff->k3.knMm + 125.0 / 192.0 * pIntCoeff->k4.knMm - 2187.0 / 6784.0 * pIntCoeff->k5.knMm + 11.0 / 84.0 * pIntCoeff->k6.knMm);

			RKCoeffs(pWorkVar, &pIntCoeff->coords7, &pIntCoeff->k7);

			pWorkVar->buffer3d.buffer1 = 71.0 / 57600.0 * pIntCoeff->k1.knMm - 71.0 / 16695.0 * pIntCoeff->k3.knMm + 71.0 / 1920.0 * pIntCoeff->k4.knMm - 17253.0 / 339200.0 * pIntCoeff->k5.knMm + 22.0 / 525.0 * pIntCoeff->k6.knMm - 1.0 / 40.0 * pIntCoeff->k7.knMm;
			error = deltaT * (pWorkVar->buffer3d.buffer1.rowwise().norm()).maxCoeff();

			// save old tDelta and adjust timestamp
			deltaTold = deltaT;
			AdjustTimeStep(pWorkVar, error, power, pDeltaT, &cond, &noFailed);
		}

		pWorkVar->coords.posEa += deltaTold * (35.0 / 384.0 * pIntCoeff->k1.knEa + 500.0 / 1113.0 * pIntCoeff->k3.knEa + 125.0 / 192.0 * pIntCoeff->k4.knEa - 2187.0 / 6784 * pIntCoeff->k5.knEa + 11.0 / 84.0 * pIntCoeff->k6.knEa); //5.order solution
		pWorkVar->coords.posMm += deltaTold * (35.0 / 384.0 * pIntCoeff->k1.knMm + 500.0 / 1113.0 * pIntCoeff->k3.knMm + 125.0 / 192.0 * pIntCoeff->k4.knMm - 2187.0 / 6784 * pIntCoeff->k5.knMm + 11.0 / 84.0 * pIntCoeff->k6.knMm);

		break;
	}
	}

	//time counting
	*pTInt += deltaTold;

	// make sure that the length of the unit vectors stay 1
	pWorkVar->coords.posMm.colwise() /= pWorkVar->coords.posMm.rowwise().norm();
	pWorkVar->coords.posEa.colwise() /= pWorkVar->coords.posEa.rowwise().norm();
}

// Brief: Calculate the runge-kutta coefficients for magnetic moment and easy axes movement
// pCoords is necessary although coords can be adressed via pWorkVar, but the function also uses coords not stored in the workingVar
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out]: Coord_S* pCoord - pointer to the current Coord_S object
// param[in/out]: RkCoeff_S* pRkCoeffs - pointer to the current RK coefficants
// return: void
void RKCoeffs(WorkingVar_S* pWorkVar, Coord_S* pCoord, RkCoeff_S* pRkCoeffs)
{
	// matlab code: fluxDensEff = extFluxDens + anisConst*dot(posMm,posEa).*posEa + fieldTherm;

	pWorkVar->buffer1d.buffer1 = (pCoord->posMm * pCoord->posEa).rowwise().sum(); //dot Product
	pWorkVar->buffer3d.buffer1 = pWorkVar->extFluxDens + (pCoord->posEa.colwise() * pWorkVar->buffer1d.buffer1) * config->getAnisConst() + pWorkVar->thermField; //effektive flux density

	// check if particles can rotate free
	if (config->getEnableMobilization())
	{
		//matlab code knEa = velEaConst1.*dot(posMm,posEa).*(posMm - dot(posMm,posEa).*posEa) + velEaConst2.*cross(torqueTherm, posEa);

		RowWiseCrossProd(&pWorkVar->thermTorque, &pCoord->posEa, &pWorkVar->buffer3d.buffer2);
		pRkCoeffs->knEa = (pCoord->posMm - pCoord->posEa.colwise() * pWorkVar->buffer1d.buffer1).colwise() * (pWorkVar->partProp.velEaConst1 * pWorkVar->buffer1d.buffer1) + pWorkVar->buffer3d.buffer2.colwise() * pWorkVar->partProp.velEaConst2;
	}
	else
	{
		pRkCoeffs->knEa.setZero();
	}

	//velMm = -velMmConst*(cross(posMm,fluxDensEff) + magDamp*cross(posMm,cross(posMm,fluxDensEff)));
	RowWiseCrossProd(&pCoord->posMm, &pWorkVar->buffer3d.buffer1, &pWorkVar->buffer3d.buffer2);
	RowWiseCrossProd(&pCoord->posMm, &pWorkVar->buffer3d.buffer2, &pWorkVar->buffer3d.buffer3);

	pRkCoeffs->knMm = -config->getVelMmConst() * (pWorkVar->buffer3d.buffer2 + config->getMagDamp() * pWorkVar->buffer3d.buffer3);
}

// Brief: Adjusting timeStep to meet tolerance requirements 
// param[in] : WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : double error - error value
// param[in] : double power - power for timestep adjustment
// param[in/out]: double* pDeltaT - new timestep
// param[in/out]: bool* pCond - new loop condidtion
// param[in/out]: bool* pNoFailed - new failed condidtion
// return: void
void AdjustTimeStep(WorkingVar_S* pWorkVar, double error, double power, double* pDeltaT, bool* pCond, bool* pNoFailed)
{
	double tDeltaNew;

	if (error > config->getAbsTol()) //the following code is from Matlab Ode45 solver
	{
		if (*pNoFailed == true)
		{
			*pNoFailed = false;
			tDeltaNew = max(config->getDeltaTMin(), (*pDeltaT) * max(0.1, 0.9 * pow(config->getAbsTol() / error, power)));
		}
		else
		{
			tDeltaNew = max(config->getDeltaTMin(), ((*pDeltaT) * 0.5));
		}

		pWorkVar->thermField *= sqrt(tDeltaNew / (*pDeltaT));
		pWorkVar->thermTorque *= sqrt(tDeltaNew / (*pDeltaT));
	}
	else
	{
		if (*pNoFailed == true)
		{
			tDeltaNew = (*pDeltaT) * min(max(0.9 * pow(config->getAbsTol() / error, power), 0.2), 5.0); 
		}
		else
		{
			tDeltaNew = (*pDeltaT);
		}

		*pCond = false;
	}

	//make sure the timestep stays below tDeltaMax
	*pDeltaT = min(tDeltaNew, config->getDeltaTMax());
	*pDeltaT = max(tDeltaNew, config->getDeltaTMin());
}
