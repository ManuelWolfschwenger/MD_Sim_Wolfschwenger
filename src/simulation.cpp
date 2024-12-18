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
	pWorkVar->coords.Bvec.setZero();                     //set direction of field
	pWorkVar->coords.Bvec.col(2).setOnes();

	// initialize the thermal fluctuations
	pWorkVar->thermTorque.resize(config->getNumbPart(), 3);
	pWorkVar->thermTorque.setZero();

	// put magnetic moments in energy minima
	InitTwoStateApprox(pWorkVar);

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
		break;
	}
	
	case SIZE_DIST_LOG:
	{
		// log-normal distributed sizes
		double dHydrMean = 2.0 * config->getRHydrMean();
		double mu = log(pow(dHydrMean, 2.0) / sqrt(pow(dHydrMean, 2.0) + pow(config->getSigma(), 2.0)));
		double sigmaLog = sqrt(log(1.0 + pow(config->getSigma() / dHydrMean, 2.0)));

		// lognormal distribution with mean = 0, stdev = 1.0
		dHydr = Rand::lognormal<MatrixXd>(config->getNumbPart(), 1, generator, mu, sigmaLog);
		break;
	}	
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

	// stuff for two state approx
	double Hk = 2 * config->getAnisEn() / config->getSatMag();
	double gamma1 = config->getGyroMr()/(1+pow(config->getMagDamp(),2.0));
	pPartProps->tau0 = 1/(2*config->getMagDamp()*gamma1)*sqrt(2*config->getMyPi()*config->getKB()*config->getTemp()/(pow(Hk,3.0)*config->getSatMag()))*(pPartProps->volMag.cwiseInverse()).cwiseSqrt();
	pPartProps->psi.setLinSpaced(200, -config->getMyPi()/50.0,2*config->getMyPi()-0.001); //-0.001 => otherwise last Element of pot Energy would be positive
	pPartProps->cosPsi = pPartProps->psi.cos();
	pPartProps->sinPsi = pPartProps->psi.sin();
	pPartProps->sinPsi2 = (pPartProps->psi.sin()).pow(2.0);
	pPartProps->dPotEn.resize(pPartProps->psi.rows());
	pPartProps->ddPotEn.resize(pPartProps->psi.rows());
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
	pCoords->omegaEa.resize(config->getNumbPart(), 3);
	pCoords->Bvec.resize(config->getNumbPart(), 3);
	pCoords->mProj.resize(config->getNumbPart(), 3);
	pCoords->vecNormal.resize(config->getNumbPart(), 3);
	pCoords->psiIs.resize(config->getNumbPart(), 1);
	pCoords->phiIs.resize(config->getNumbPart(), 1);
	pCoords->state.resize(config->getNumbPart(), 1);

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
	pIntCoeff->coords3.posEa.resize(config->getNumbPart(), 3);
	pIntCoeff->coords4.posEa.resize(config->getNumbPart(), 3);
	pIntCoeff->coords5.posEa.resize(config->getNumbPart(), 3);
	pIntCoeff->coords6.posEa.resize(config->getNumbPart(), 3);
	pIntCoeff->coords7.posEa.resize(config->getNumbPart(), 3);

	pIntCoeff->k1.knEa.resize(config->getNumbPart(), 3);
	pIntCoeff->k2.knEa.resize(config->getNumbPart(), 3);
	pIntCoeff->k3.knEa.resize(config->getNumbPart(), 3);
	pIntCoeff->k4.knEa.resize(config->getNumbPart(), 3);
	pIntCoeff->k5.knEa.resize(config->getNumbPart(), 3);
	pIntCoeff->k6.knEa.resize(config->getNumbPart(), 3);
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

// Brief: Initializes Two state approximation / puts mag moment in energy mimima
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// return: void
void InitTwoStateApprox(WorkingVar_S* pWorkVar)
{
	ArrayXXd extrema;

	//normalize vectors and project m onto the plane spanned by B and n
	VectorProjection(pWorkVar);

	pWorkVar->coords.psiIs = ((pWorkVar->coords.posMm*pWorkVar->coords.posEa).rowwise().sum()).acos(); //angle between posEa and posMm after projection
	pWorkVar->coords.phiIs = ((pWorkVar->coords.posEa*pWorkVar->coords.Bvec).rowwise().sum()).acos();  //angle between posEa and magnetic field

    //set external flux density to determine initial positions
	if (config->getTMag() > 0)
	{
		pWorkVar->extFluxDens.col(2).setConstant(config->getMagFluxDens());
	}
	
	for (int i = 0; i < config->getNumbPart(); i++)
	{
		extrema = FindLocalMinMax(pWorkVar,i);

		////////////////// put mag moments in minima //////////////////////////////
		if (extrema.cols() == 2) //check if there are two minima
		{
			if (pWorkVar->coords.psiIs(i) > extrema(2,0) && pWorkVar->coords.psiIs(i) < extrema(2,1))
			{
				pWorkVar->coords.psiIs(i) = extrema(0,1);
				pWorkVar->coords.state(i) = 2;
			}
			else
			{
				pWorkVar->coords.psiIs(i) = extrema(0,0);
				pWorkVar->coords.state(i) = 1;
			}
		}
		else
		{
			pWorkVar->coords.psiIs(i) = extrema(0); //if there is only one minimum
			pWorkVar->coords.state(i) = 0;
		} 	
		
		pWorkVar->coords.posMm.row(i) = cos(pWorkVar->coords.psiIs(i))*pWorkVar->coords.posEa.row(i) +
			sin(pWorkVar->coords.psiIs(i))*(pWorkVar->coords.Bvec.row(i) - (pWorkVar->coords.Bvec.row(i)*pWorkVar->coords.posEa.row(i)).sum()*pWorkVar->coords.posEa.row(i));
	}
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
		if ((int)(floor(tInt / config->getTEnd() * 100)) == tCount)
		{
			runTime = (clock() - tStart) / CLOCKS_PER_SEC;

			cout << "actual timestep: " << deltaT << "s" << endl;
			cout << "actual time: " << tInt << "s" << endl;
			cout << "simulation progress: " << tCount << "%" << endl;
			cout << "simulation will end in about " << (config->getTEnd() - tInt) * (runTime / tInt) / 3600.0 << "h" << "\n" << endl;

			tCount++;
		}

		if ((int)(floor(tInt / config->getTEnd() * config->getDataPoints())) == tCountData)
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

		//two state approximation
		TwoStateApprox(pWorkVar, deltaT);

		//evaluation of angular velocity
		evalOmega(pWorkVar, pOutputVar);
		
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
	saveData("mZ.txt", Map<VectorXd, Unaligned>(pOutputVar->mZ.data(), pOutputVar->mZ.size()), 0);
	saveData("nZ.txt", Map<VectorXd, Unaligned>(pOutputVar->nZ.data(), pOutputVar->nZ.size()), 0);
	saveData("omegaDiffx.txt", Map<VectorXd, Unaligned>(pOutputVar->omegaX.data(), pOutputVar->omegaX.size()), 0);
	saveData("omegaDiffy.txt", Map<VectorXd, Unaligned>(pOutputVar->omegaY.data(), pOutputVar->omegaY.size()), 0);
	saveData("omegaDiffz.txt", Map<VectorXd, Unaligned>(pOutputVar->omegaZ.data(), pOutputVar->omegaZ.size()), 0);
	
	if (config->getDataPoints() > steps)
	{
		cout << "\n attention: number of datapoints is smaller than simulation steps" << endl;
		cout << "simulation ended" << endl;
	}
	else
	{
		cout << "\n simulation ended successfuly" << endl;
	}
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
}

// Brief: normalize vectors and project m onto the plane spanned by B and n / could be probably parallelized
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// return: void
void VectorProjection(WorkingVar_S* pWorkVar)
{
	RowWiseCrossProd(&pWorkVar->coords.Bvec, &pWorkVar->coords.posEa, &pWorkVar->coords.vecNormal);
	pWorkVar->coords.vecNormal.colwise() /= pWorkVar->coords.vecNormal.rowwise().norm();

	pWorkVar->coords.mProj = pWorkVar->coords.posMm - pWorkVar->coords.vecNormal.colwise()*(pWorkVar->coords.posMm*pWorkVar->coords.vecNormal).rowwise().sum();
	pWorkVar->coords.posMm = pWorkVar->coords.mProj.colwise() / pWorkVar->coords.mProj.rowwise().norm();
}

// Brief: normalize vectors and project m onto the plane spanned by B and n / could be probably parallelized
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// return: void
void VectorProjectionSolver(WorkingVar_S* pWorkVar, ArrayXXd* pPosEa)
{
	RowWiseCrossProd(&pWorkVar->coords.Bvec, pPosEa, &pWorkVar->coords.vecNormal);
	pWorkVar->coords.vecNormal.colwise() /= pWorkVar->coords.vecNormal.rowwise().norm();

	pWorkVar->coords.mProj = pWorkVar->coords.posMm - pWorkVar->coords.vecNormal.colwise()*(pWorkVar->coords.posMm*pWorkVar->coords.vecNormal).rowwise().sum();
	pWorkVar->coords.posMm = pWorkVar->coords.mProj.colwise() / pWorkVar->coords.mProj.rowwise().norm();
}

// Brief: find local extrema of potential energy with Newton method
// param[in]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: int i - index of particle
// return: void
 ArrayXXd FindLocalMinMax(WorkingVar_S* pWorkVar, int i)
 {
	double x,f,df;
	double B = pWorkVar->extFluxDens(0,2); //sets magfluxdensity according to time

	vector<double>  minPsi, minEn, maxPsi, maxEn;

	pWorkVar->partProp.dPotEn = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*pWorkVar->partProp.sinPsi*pWorkVar->partProp.cosPsi +
			config->getSatMag()*pWorkVar->partProp.volMag(i)*B*(pWorkVar->partProp.psi - pWorkVar->coords.phiIs(i)).sin();

	for (int j = 1; j < pWorkVar->partProp.dPotEn.rows(); j++)
	{
		//finding local minima
		if (pWorkVar->partProp.dPotEn(j-1) < 0 && pWorkVar->partProp.dPotEn(j) > 0) //=> Minimum crossed
		{
			x = pWorkVar->partProp.psi(j-1);

			//Newton iteration
			f = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*sin(x)*cos(x) +
					config->getSatMag()*pWorkVar->partProp.volMag(i)*B*sin(x - pWorkVar->coords.phiIs(i));

			while (abs(f) > config->getErrTolMinMax())
			{
				df = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*cos(2*x) +
					config->getSatMag()*pWorkVar->partProp.volMag(i)*B*cos(x - pWorkVar->coords.phiIs(i));
				
				x -= f/df;

				f = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*sin(x)*cos(x) +
					config->getSatMag()*pWorkVar->partProp.volMag(i)*B*sin(x - pWorkVar->coords.phiIs(i));
			}	

			minPsi.push_back(x);
			minEn.push_back(config->getAnisEn()*pWorkVar->partProp.volMag(i)*pow(sin(x),2.0) - 
			config->getSatMag()*pWorkVar->partProp.volMag(i)*B*cos(x - pWorkVar->coords.phiIs(i)));
		}	

		//finding local maxima
		if (pWorkVar->partProp.dPotEn(j-1) > 0 && pWorkVar->partProp.dPotEn(j) < 0) //=> Maximum crossed
		{
			x = pWorkVar->partProp.psi(j-1);

			//Newton iteration
			f = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*sin(x)*cos(x) +
					config->getSatMag()*pWorkVar->partProp.volMag(i)*B*sin(x - pWorkVar->coords.phiIs(i));

			while (abs(f) > config->getErrTolMinMax())
			{
				df = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*cos(2*x) +
					config->getSatMag()*pWorkVar->partProp.volMag(i)*B*cos(x - pWorkVar->coords.phiIs(i));
				
				x -= f/df;

				f = 2*config->getAnisEn()*pWorkVar->partProp.volMag(i)*sin(x)*cos(x) +
					config->getSatMag()*pWorkVar->partProp.volMag(i)*B*sin(x - pWorkVar->coords.phiIs(i));
			}	

			maxPsi.push_back(x);
			maxEn.push_back(config->getAnisEn()*pWorkVar->partProp.volMag(i)*pow(sin(x),2.0) - 
			config->getSatMag()*pWorkVar->partProp.volMag(i)*B*cos(x - pWorkVar->coords.phiIs(i)));
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

// Brief: Two state approximation for jumps of magnetic moment
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: double deltaT - current integration timestep
// return: void
void TwoStateApprox(WorkingVar_S* pWorkVar, double deltaT)
{
	ArrayXXd extrema;
	double dE1, dE2, r1, r2, numb, p;
	bool jump;

	//normalize vectors and project m onto the plane spanned by B and n
	VectorProjection(pWorkVar);

	pWorkVar->coords.phiIs = ((pWorkVar->coords.posEa*pWorkVar->coords.Bvec).rowwise().sum()).acos();  //angle between posEa and magnetic field

	for (int i = 0; i < config->getNumbPart(); i++)
	{
		extrema = FindLocalMinMax(pWorkVar,i);

		///////////////// put mag moments back in minima after easy axis movement //////////////////////////////
		if (extrema.cols() == 2) //check if there are two minima
		{
			if (pWorkVar->coords.state(i) == 1)
			{
				pWorkVar->coords.psiIs(i) = extrema(0,0);
			}
			else if (pWorkVar->coords.state(i) == 2)
			{
				pWorkVar->coords.psiIs(i) = extrema(0,1);
			}
			else if (pWorkVar->coords.state(i) == 0)
			{
				if (pWorkVar->coords.psiIs(i) > extrema(2,0) && pWorkVar->coords.psiIs(i) < extrema(2,1))
				{
					pWorkVar->coords.psiIs(i) = extrema(0,1);
					pWorkVar->coords.state(i) = 2;
				}
				else
				{
					pWorkVar->coords.psiIs(i) = extrema(0,0);
					pWorkVar->coords.state(i) = 1;
				}
			}
		}
		else
		{
			pWorkVar->coords.psiIs(i) = extrema(0); //if there is only one minimum
			pWorkVar->coords.state(i) = 0;
		} 	
		
		pWorkVar->coords.posMm.row(i) = cos(pWorkVar->coords.psiIs(i))*pWorkVar->coords.posEa.row(i) +
			sin(pWorkVar->coords.psiIs(i))*(pWorkVar->coords.Bvec.row(i) - (pWorkVar->coords.Bvec.row(i)*pWorkVar->coords.posEa.row(i)).sum()*pWorkVar->coords.posEa.row(i)); 

		////////////////calculate possible jumps of mag moments ////////////////////////////
		if (extrema.cols() == 2)
		{
			//energy differences
			dE1 = extrema.row(3).minCoeff() - extrema(1,0);
			dE2 = extrema.row(3).minCoeff() - extrema(1,1);
			r1 = 1.0/(2.0*pWorkVar->partProp.tau0(i))*exp(-dE1/(config->getKB()*config->getTemp())); //switching rate from minimum 1 to minimum 2
			r2 = 1.0/(2.0*pWorkVar->partProp.tau0(i))*exp(-dE2/(config->getKB()*config->getTemp())); //switching rate from minimum 2 to minimum 1

			numb = ((double) rand() / (RAND_MAX + 1.0)); //random number between 0 and 1
			jump = false;

			if (pWorkVar->coords.state(i) == 1)
			{
				p = r1/(r1+r2)*(1-exp(-(r1+r2)*deltaT));

				if (p > numb)
				{
					pWorkVar->coords.state(i) = 2;
					pWorkVar->coords.psiIs(i) = extrema(0,1);
					pWorkVar->coords.posMm.row(i) = cos(pWorkVar->coords.psiIs(i))*pWorkVar->coords.posEa.row(i) +
			 			sin(pWorkVar->coords.psiIs(i))*(pWorkVar->coords.Bvec.row(i) - (pWorkVar->coords.Bvec.row(i)*pWorkVar->coords.posEa.row(i)).sum()*pWorkVar->coords.posEa.row(i));
					jump = true;
				}
			}
			if (pWorkVar->coords.state(i) == 2 && jump == false)
			{
				p = r2/(r1+r2)*(1-exp(-(r1+r2)*deltaT));

				if (p > numb)
				{
					pWorkVar->coords.state(i) = 1;
					pWorkVar->coords.psiIs(i) = extrema(0,0);
					pWorkVar->coords.posMm.row(i) = cos(pWorkVar->coords.psiIs(i))*pWorkVar->coords.posEa.row(i) +
			 			sin(pWorkVar->coords.psiIs(i))*(pWorkVar->coords.Bvec.row(i) - (pWorkVar->coords.Bvec.row(i)*pWorkVar->coords.posEa.row(i)).sum()*pWorkVar->coords.posEa.row(i));
				}
			}
		} 
	}	
}

// Brief: Two state approximation for readjustment of magnetic moment positions in partial solver steps
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]	: Buffer_S* pBuffer - pointer to Buffer_S object
// param[in]: double deltaT - current integration timestep
// return: void
void TwoStateApproxSolver(WorkingVar_S* pWorkVar, ArrayXXd* pPosEa)
{
	ArrayXXd extrema;
	Array<double,1,3> vecN;
	double dE1, dE2, r1, r2, numb, p;

	//normalize vectors and project m onto the plane spanned by B and n
	VectorProjectionSolver(pWorkVar, pPosEa);

	pWorkVar->coords.phiIs = (((*pPosEa)*pWorkVar->coords.Bvec).rowwise().sum()).acos();  //angle between posEa and magnetic field

	for (int i = 0; i < config->getNumbPart(); i++)
	{
		extrema = FindLocalMinMax(pWorkVar,i);
		//cout << extrema << endl;

		////////////////// put mag moments back in minima after easy axis movement //////////////////////////////
		if (extrema.cols() == 2) //check if there are two minima
		{
			if (pWorkVar->coords.state(i) == 1)
			{
				pWorkVar->coords.psiIs(i) = extrema(0,0);
			}
			else if (pWorkVar->coords.state(i) == 2)
			{
				pWorkVar->coords.psiIs(i) = extrema(0,1);
			}
			else if (pWorkVar->coords.state(i) == 0)
			{
				if (pWorkVar->coords.psiIs(i) > extrema(2,0) && pWorkVar->coords.psiIs(i) < extrema(2,1))
				{
					pWorkVar->coords.psiIs(i) = extrema(0,1);
					pWorkVar->coords.state(i) = 2;
				}
				else
				{
					pWorkVar->coords.psiIs(i) = extrema(0,0);
					pWorkVar->coords.state(i) = 1;
				}
			}
		}
		else
		{
			pWorkVar->coords.psiIs(i) = extrema(0); //if there is only one minimum
			pWorkVar->coords.state(i) = 0;
		} 	
		
		vecN = pWorkVar->coords.Bvec.row(i) - (pWorkVar->coords.Bvec.row(i)*(*pPosEa).row(i)).sum()*(*pPosEa).row(i);
		vecN.colwise() /=vecN.rowwise().norm();
		pWorkVar->coords.posMm.row(i) = cos(pWorkVar->coords.psiIs(i))*(*pPosEa).row(i) + sin(pWorkVar->coords.psiIs(i))*vecN;
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
	case 0:
	{
		power = 0.5;

		while (cond)
		{
			deltaT = (*pDeltaT);

			RKCoeffs(pWorkVar, &pWorkVar->coords, &pIntCoeff->k1);

			pIntCoeff->coords2.posEa = pWorkVar->coords.posEa + deltaT * pIntCoeff->k1.knEa;
			pIntCoeff->coords2.posEa.rowwise().normalize(); //normalize vectors
			TwoStateApproxSolver(pWorkVar, &pIntCoeff->coords2.posEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords2, &pIntCoeff->k2);

			pWorkVar->buffer3d.buffer2 = 0.5 * pIntCoeff->k1.knEa - 0.5 * pIntCoeff->k2.knEa;
			error = deltaT * (pWorkVar->buffer3d.buffer2.rowwise().norm()).maxCoeff();
			
			// save old tDelta and adjust timestep
			deltaTold = deltaT;
			AdjustTimeStep(pWorkVar, error, power, pDeltaT, &cond, &noFailed);
		}

		pWorkVar->coords.posEa += deltaTold * (0.5 * pIntCoeff->k1.knEa + 0.5 * pIntCoeff->k2.knEa) ; //2.order solution

		break;
	}
	case 1:
	{
		power = 1.0 / 3.0;

		while (cond)
		{
			// pre load tDelta for faster execution
			deltaT = (*pDeltaT);

			RKCoeffs(pWorkVar, &pWorkVar->coords, &pIntCoeff->k1);

			pIntCoeff->coords2.posEa = pWorkVar->coords.posEa + deltaT * 1.0 / 2.0 * pIntCoeff->k1.knEa;
			pIntCoeff->coords2.posEa.rowwise().normalize(); //normalize vectors
			TwoStateApproxSolver(pWorkVar, &pIntCoeff->coords2.posEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords2, &pIntCoeff->k2);

			pIntCoeff->coords3.posEa = pWorkVar->coords.posEa + deltaT * 3.0 / 4.0 * pIntCoeff->k2.knEa;
			pIntCoeff->coords3.posEa.rowwise().normalize(); //normalize vectors
			TwoStateApproxSolver(pWorkVar, &pIntCoeff->coords3.posEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords3, &pIntCoeff->k3);

			pIntCoeff->coords4.posEa = pWorkVar->coords.posEa + deltaT * (2.0 / 9.0 * pIntCoeff->k1.knEa + 1.0 / 3.0 * pIntCoeff->k2.knEa + 4.0 / 9.0 * pIntCoeff->k3.knEa);
			pIntCoeff->coords4.posEa.rowwise().normalize(); //normalize vectors
			TwoStateApproxSolver(pWorkVar, &pIntCoeff->coords4.posEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords4, &pIntCoeff->k4);

			pWorkVar->buffer3d.buffer2 = -5.0 / 72.0 * pIntCoeff->k1.knEa + 1.0 / 12.0 * pIntCoeff->k2.knEa + 1.0 / 9.0 * pIntCoeff->k3.knEa - 1.0 / 8.0 * pIntCoeff->k4.knEa;
			error = deltaT * (pWorkVar->buffer3d.buffer2.rowwise().norm()).maxCoeff();

			// save old tDelta and adjust timestamp
			deltaTold = deltaT;
			AdjustTimeStep(pWorkVar, error, power, pDeltaT, &cond, &noFailed);
		}

		pWorkVar->coords.posEa = pIntCoeff->coords4.posEa;

		break;
	}
	case 2: 
	{
		power = 0.2;

		while (cond)
		{
			// pre load tDelta for faster execution
			deltaT = (*pDeltaT);

			RKCoeffs(pWorkVar, &pWorkVar->coords, &pIntCoeff->k1);

			pIntCoeff->coords2.posEa = pWorkVar->coords.posEa + deltaT * 1.0 / 5.0 * pIntCoeff->k1.knEa;
			pIntCoeff->coords2.posEa.rowwise().normalize(); //normalize vectors
			TwoStateApproxSolver(pWorkVar, &pIntCoeff->coords2.posEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords2, &pIntCoeff->k2);

			pIntCoeff->coords3.posEa = pWorkVar->coords.posEa + deltaT * (3.0 / 40.0 * pIntCoeff->k1.knEa + 9.0 / 40.0 * pIntCoeff->k2.knEa);
			pIntCoeff->coords3.posEa.rowwise().normalize(); //normalize vectors
			TwoStateApproxSolver(pWorkVar, &pIntCoeff->coords3.posEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords3, &pIntCoeff->k3);

			pIntCoeff->coords4.posEa = pWorkVar->coords.posEa + deltaT * (44.0 / 45.0 * pIntCoeff->k1.knEa - 56.0 / 15.0 * pIntCoeff->k2.knEa + 32.0 / 9.0 * pIntCoeff->k3.knEa);
			pIntCoeff->coords4.posEa.rowwise().normalize(); //normalize vectors
			TwoStateApproxSolver(pWorkVar, &pIntCoeff->coords4.posEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords4, &pIntCoeff->k4);

			pIntCoeff->coords5.posEa = pWorkVar->coords.posEa + deltaT * (19372.0 / 6561.0 * pIntCoeff->k1.knEa - 25360.0 / 2187.0 * pIntCoeff->k2.knEa + 64448.0 / 6561.0 * pIntCoeff->k3.knEa - 212.0 / 729.0 * pIntCoeff->k4.knEa);
			pIntCoeff->coords5.posEa.rowwise().normalize(); //normalize vectors
			TwoStateApproxSolver(pWorkVar, &pIntCoeff->coords5.posEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords5, &pIntCoeff->k5);

			pIntCoeff->coords6.posEa = pWorkVar->coords.posEa + deltaT * (9017.0 / 3168.0 * pIntCoeff->k1.knEa - 355.0 / 33.0 * pIntCoeff->k2.knEa + 46732.0 / 5247.0 * pIntCoeff->k3.knEa + 49.0 / 176.0 * pIntCoeff->k4.knEa - 5103.0 / 18656.0 * pIntCoeff->k5.knEa);
			pIntCoeff->coords6.posEa.rowwise().normalize();
			TwoStateApproxSolver(pWorkVar, &pIntCoeff->coords6.posEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords6, &pIntCoeff->k6);

			pIntCoeff->coords7.posEa = pWorkVar->coords.posEa + deltaT * (35.0 / 384.0 * pIntCoeff->k1.knEa + 500.0 / 1113.0 * pIntCoeff->k3.knEa + 125.0 / 192.0 * pIntCoeff->k4.knEa - 2187.0 / 6784.0 * pIntCoeff->k5.knEa + 11.0 / 84.0 * pIntCoeff->k6.knEa);
			pIntCoeff->coords7.posEa.rowwise().normalize();
			TwoStateApproxSolver(pWorkVar, &pIntCoeff->coords7.posEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords7, &pIntCoeff->k7);

			pWorkVar->buffer3d.buffer1 = 71.0 / 57600.0 * pIntCoeff->k1.knEa - 71.0 / 16695.0 * pIntCoeff->k3.knEa + 71.0 / 1920.0 * pIntCoeff->k4.knEa - 17253.0 / 339200.0 * pIntCoeff->k5.knEa + 22.0 / 525.0 * pIntCoeff->k6.knEa - 1.0 / 40.0 * pIntCoeff->k7.knEa;
			error = deltaT * (pWorkVar->buffer3d.buffer1.rowwise().norm()).maxCoeff();

			// save old tDelta and adjust timestamp
			deltaTold = deltaT;
			AdjustTimeStep(pWorkVar, error, power, pDeltaT, &cond, &noFailed);
		}

		pWorkVar->coords.posEa += deltaTold * (35.0 / 384.0 * pIntCoeff->k1.knEa + 500.0 / 1113.0 * pIntCoeff->k3.knEa + 125.0 / 192.0 * pIntCoeff->k4.knEa - 2187.0 / 6784 * pIntCoeff->k5.knEa + 11.0 / 84.0 * pIntCoeff->k6.knEa); //5.order solution;;
		
		break;
	}
	}

	//time counting
	*pTInt += deltaTold;

	pWorkVar->coords.posEa.rowwise().normalize();
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

	pWorkVar->buffer1d.buffer1 = (pWorkVar->coords.posMm * pCoord->posEa).rowwise().sum(); //dot Product

	// check if particles can rotate free
	if (config->getEnableMobilization())
	{
		//add vorticity in y-direction to torque/zetaRot = angular velocity
		pWorkVar->buffer3d.buffer1 = pWorkVar->thermTorque.colwise() * pWorkVar->partProp.velEaConst2;
		pWorkVar->buffer3d.buffer1.col(1) += config->getShearRate() * 0.5;

		RowWiseCrossProd(&pWorkVar->buffer3d.buffer1, &pCoord->posEa, &pWorkVar->buffer3d.buffer2);

		pRkCoeffs->knEa = (pWorkVar->coords.posMm - pCoord->posEa.colwise() * pWorkVar->buffer1d.buffer1).colwise() * (pWorkVar->partProp.velEaConst1 * pWorkVar->buffer1d.buffer1) + pWorkVar->buffer3d.buffer2;
	}
	else
	{
		pRkCoeffs->knEa.setZero();
	}
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
	double deltaTnew;

	if (error > config->getAbsTol()) //the following code is from Matlab Ode45 solver
	{
		if (*pNoFailed == true)
		{
			*pNoFailed = false;
			deltaTnew = max(config->getDeltaTMin(), (*pDeltaT) * max(0.1, 0.9 * pow(config->getAbsTol() / error, power)));
		}
		else
		{
			deltaTnew = max(config->getDeltaTMin(), ((*pDeltaT) * 0.5));
		}

		pWorkVar->thermTorque *= sqrt(deltaTnew / (*pDeltaT));
	}
	else
	{
		if (*pNoFailed == true)
		{
			deltaTnew = (*pDeltaT) * min(max(0.9 * pow(config->getAbsTol() / error, power), 0.2), 5.0); 
		}
		else
		{
			deltaTnew = (*pDeltaT);
		}

		*pCond = false;
	}

	//make sure the timestep stays below tDeltaMax
	*pDeltaT = min(max(deltaTnew, config->getDeltaTMin()), config->getDeltaTMax());
}
