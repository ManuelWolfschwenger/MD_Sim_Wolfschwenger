/*
* Name : simulation.cpp
* Auto :
* Brief: Simulation source file
*/

// own header file
#include "simulation.h"
#include "configuration.h"

// global variables
static WorkingVar_S     workVar;    // simulation working varibales
static Buffer_S         buffer;     // buffer 
static IntCoefficient_S  intCoeff;   // integration coefficant varibales
static OutputVar_S      outputVar;  // output variables
static Params_S         params;     // global parameters not known at compile time
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

// Brief: This function initializes the simulation
// return: void
void InitSimulation(void)
{
	WorkingVar_S* pWorkVar;      // pointer to the simulation working variables
	Buffer_S* pBuffer;           // pointer to buffer
	IntCoefficient_S* pIntCoeff;  // pointer to the integration coefficant variables
	Params_S* pParams;           // pointer to parameters not known at compile time
	OutputVar_S* pOutputVar;     // pointer to output variables

	// initialize the pointers
	pWorkVar = &workVar;
	pBuffer = &buffer;
	pIntCoeff = &intCoeff;
	pParams = &params;
	pOutputVar = &outputVar;
	
	// initialize the particle properties
	InitParticleProps(&pWorkVar->partProp);

	// set Params
	pParams->lBox = pow(pWorkVar->partProp.volHydr.sum() / config->getVolFrac(), 1.0 / 3.0);  //box volume^1/3
	pParams->volBox = pow(pParams->lBox, 3.0);
	pParams->transOutput = false;

	if (config->getReadFromFile())
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

	// init short range cut off radius
	InitCutOffSR(pParams, &pWorkVar->partProp);
	
	// resizing interaction arrays
	int maxPairs = (config->getNumbPart() * (config->getNumbPart() - 1)) / 2;             // all pairs
	pWorkVar->intLists.intVecsCoords.resize(maxPairs, 3);   // resize to number of pairs rows and 3 colums for x,y,z
	pWorkVar->intLists.intVecsLen.resize(maxPairs);         // length of vectors
	pWorkVar->intLists.intPairs.resize(2 * maxPairs);       // indizes of pairs
	pWorkVar->intLists.intRangeSR.resize(maxPairs);           // true/false if in interaction range

	// initialize the external field 
	pWorkVar->extFluxDens.resize(config->getNumbPart(), 3);
	pWorkVar->extFluxDens.setZero();  //in case we start with a relaxation, dont need to set to zero every iteration in the loop

	// initialize demag field
	pWorkVar->demagFluxDens.resize(config->getNumbPart(), 3);
	pWorkVar->demagFluxDens.setZero();

	// initialize forces
	pWorkVar->intForce.resize(config->getNumbPart(), 3);
	pWorkVar->intForce.setZero();

	// initialize the thermal fluctuations
	pWorkVar->thermTorque.resize(config->getNumbPart(), 3);
	pWorkVar->thermForce.resize(config->getNumbPart(), 3);

	// compute optimal Ewald-Parameters
	CompEwaldParams(pWorkVar, pBuffer, pParams);
	InitEwaldConst(pParams, pBuffer);

	//init calculation of transport coefficients
	pOutputVar->posRef = pWorkVar->coordsTrans.pos;
	pOutputVar->truePos = pWorkVar->coordsTrans.pos;

	//store data in txt
	WriteData2TXT(pWorkVar, pParams, "data.txt");

    //////////// stuff to set initial mag momnets in correct energy minima
	// set external magnetic field depending on time
	DefExtField(&pWorkVar->extFluxDens, 0);

	// compute forces
	EvalSinCos(pWorkVar, pBuffer, pParams);
	CompInteractions(pWorkVar, pBuffer, pParams, pOutputVar);

	// put magnetic moments in energy minima, needs interacting field to put mag moments in correct minima
	InitTwoStateApprox(pWorkVar);
}

// Brief: This function sets the particle properties like radius, mass, friction coefficients...
// param[in/out]: PartProps_S* pPartProps - pointer to a PartProps_S object
// return: void
void InitParticleProps(PartProps_S* pPartProps)
{
	pPartProps->rMag.resize(config->getNumbPart());         // magnetic radius
	pPartProps->rHydr.resize(config->getNumbPart());        // hydrodynamic radius

	// switch depending on the size distribution
	switch (config->getSizeDist())
	{
	case SIZE_DIST_EQUAL:
	{
		// set all particle radii to same value
		pPartProps->rMag.setConstant(config->getRMagMean());

		if (config->getEnableCoatingPot())
		{
			// calculate equivalent hard sphere diameter from repulsive potential	
			//(integration works till dHydr = dMag, when the steric potential is further decreased the equivalent hard sphere diameter is equal to the magnetic diameter (Rosensweig)
			double dHydr = IntegralEffDiameter(2 * config->getRMagMean());

			if (dHydr >= 2*config->getRMagMean() )
			{
				pPartProps->rHydr.setConstant(dHydr / 2.0);
			}
			else
			{
				pPartProps->rHydr.setConstant(config->getRMagMean());
			}		
		}
		else
		{
			// use thickness of shell to add it to magnetic radius
			pPartProps->rHydr.setConstant(config->getRMagMean() + config->getDShell());
		}
		break;
	}
	
	case SIZE_DIST_LOG:
	{
		// log-normal distributed sizes
		double mu = log(pow(config->getRMagMean(), 2.0) / sqrt(pow(config->getRMagMean(), 2.0) + pow(config->getSigma(), 2.0))); //log normal parameters
		double sigmaLog = sqrt(log(1.0 + pow(config->getSigma() / config->getRMagMean() , 2.0)));

		pPartProps->rMag = Rand::lognormal<MatrixXd>(config->getNumbPart(), 1, generator, mu, sigmaLog);

		if (config->getEnableCoatingPot())
		{
			double dHydr;

			// calculate equivalent hard sphere diameter from repulsive potential
			for (int i = 0; i < config->getNumbPart(); i++)
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
			pPartProps->rHydr = pPartProps->rMag + config->getDShell();
		}
	break;
	}	
	}

	pPartProps->volMag = 4.0 / 3.0 * config->getMyPi() * pPartProps->rMag.pow(3.0);														// mag volume
	pPartProps->volHydr = 4.0 / 3.0 * config->getMyPi() * pPartProps->rHydr.pow(3.0);													// hydrdyn volume
	pPartProps->magMom = config->getSatMag() * pPartProps->volMag;																		// mag Moment
	pPartProps->zetaRot = 8.0 * config->getMyPi() * config->getVis() * pPartProps->rHydr.pow(3.0);										// rotational friction coefficient
	pPartProps->zetaTrans = 6.0 * config->getMyPi() * config->getVis() * pPartProps->rHydr;												// translational friction coefficient
	pPartProps->velEaConst1 = 2.0 * config->getAnisEn() * (pPartProps->volMag / pPartProps->zetaRot);									// constants for later DGL solving
	pPartProps->velEaConst2 = pPartProps->zetaRot.cwiseInverse();
	pPartProps->magMomVec.resize(config->getNumbPart(), 3); //array with magnetic moment vectors
	pPartProps->magMomSum.resize(1, 3);        //array with summed components
	
	// stuff for two state approx
	double Hk = 2 * config->getAnisEn() / config->getSatMag();
	double gamma1 = config->getGyroMr()/(1+pow(config->getMagDamp(),2.0));
	pPartProps->tau0 = 1/(2*config->getMagDamp()*gamma1)*sqrt(2*config->getMyPi()*config->getKB()*config->getTemp()/(pow(Hk,3.0)*config->getSatMag()))*(pPartProps->volMag.cwiseInverse()).cwiseSqrt();
	pPartProps->psi.setLinSpaced(200, -config->getMyPi()/50.0,2*config->getMyPi()-0.001); //-0.001 => otherwise last Element of pot Energy would be positive
	pPartProps->cosPsi = pPartProps->psi.cos();
	pPartProps->sinPsi = pPartProps->psi.sin();
	pPartProps->sinPsi2 = (pPartProps->psi.sin()).pow(2.0);
	pPartProps->dPotEn.resize(pPartProps->psi.rows());
	pPartProps->ddPotEn.resize(pPartProps->psi.rows());																		// constants for later DGL solving
}         

// Brief: This function sets the initial particle positions (rotational and translational) according to the definition in settings
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Params_S* pParams - pointer to params_S object
// return: void
void InitCoordsFromSketch(WorkingVar_S* pWorkVar, Params_S* pParams)
{
	/////////////////////translational positions/////////////////////////////////

	//define size of the dynamic arrays
	pWorkVar->coordsTrans.pos.resize(config->getNumbPart(), 3);

	switch (config->getInitConfigTrans())
	{
	case 1: //random positions
	{
		Array<double, 1, 3> xyzVec;
		Array<double, 1, 3> xyzDiff;

		double rDiff;

		// finds random positions of particles in the box and makes sure that do not overlap
		for (int i = 0; i < config->getNumbPart(); i++)
		{
		setRandomCoords: //goto statement

			xyzVec = Rand::balanced<ArrayXXd>(1, 3, generator);
			xyzVec *= (pParams->lBox / 2.0 - pWorkVar->partProp.rHydr(i));

			if (i != 0) // first particle can be everywhere
			{
				for (int j = 0; j < i; j++) //all particles exceptional the one with number ip
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
		double gap = pParams->lBox / config->getNpDim(); //gap between particles
		int ip = 0;
		cout << "number of particles per dimension in the grid:" << config->getNpDim() << endl;

		for (int ix = 0; ix < config->getNpDim(); ix++)
		{
			for (int iy = 0; iy < config->getNpDim(); iy++)
			{
				for (int iz = 0; iz < config->getNpDim(); iz++)
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
	ArrayXd azimuth(config->getNumbPart());
	ArrayXd elevation(config->getNumbPart());

	//define size of the dynamic arrays
	pWorkVar->coordsRot.posEa.resize(config->getNumbPart(), 3);
	pWorkVar->coordsRot.posMm.resize(config->getNumbPart(), 3);
	pWorkVar->coordsRot.Bvec.resize(config->getNumbPart(), 3);
	pWorkVar->coordsRot.mProj.resize(config->getNumbPart(), 3);
	pWorkVar->coordsRot.vecNormal.resize(config->getNumbPart(), 3);
	pWorkVar->coordsRot.psiIs.resize(config->getNumbPart(), 1);
	pWorkVar->coordsRot.phiIs.resize(config->getNumbPart(), 1);
	pWorkVar->coordsRot.state.resize(config->getNumbPart(), 1);

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

		// all moments lie on their easy axes, / easy axis random see above
		pWorkVar->coordsRot.posMm = pWorkVar->coordsRot.posEa;
		break;
	}
	case 3:
	{
		// easy axis in z-direction
		pWorkVar->coordsRot.posEa.col(0).setZero(); //x-direction
		pWorkVar->coordsRot.posEa.col(1).setZero(); //y-direction
		pWorkVar->coordsRot.posEa.col(2).setOnes(); //z-direction

		// all moments lie on their easy axes
		pWorkVar->coordsRot.posMm = pWorkVar->coordsRot.posEa;
		break;
	}
	}
}

// Brief: This function sets the initial particle positions from last simulation step of the last simulation, not updated anymore => just keep it if i need it at some time!!!
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: Params_S* pParams - pointer to params_S object
// return: void
void InitCoordsFromSim(WorkingVar_S* pWorkVar, Params_S* pParams)
{
	MatrixXd initCoords = openData("coords_final.txt");

	cout << initCoords << "\n" << endl;
	// check if dimensions are the one that are expected Nx3
	cout << "the array has the dimensions " << initCoords.rows() << "x" << initCoords.cols() << "\n" << endl;

	//define size of the dynamic arrays
	pWorkVar->coordsTrans.pos.resize(config->getNumbPart(), 3);

	pWorkVar->coordsRot.posEa.resize(config->getNumbPart(), 3);
	pWorkVar->coordsRot.posMm.resize(config->getNumbPart(), 3);

	pWorkVar->coordsTrans.pos.col(0) = initCoords.col(0);
	pWorkVar->coordsTrans.pos.col(1) = initCoords.col(1);
	pWorkVar->coordsTrans.pos.col(2) = initCoords.col(2);

	pWorkVar->coordsRot.posMm.col(0) = initCoords.col(6);
	pWorkVar->coordsRot.posMm.col(1) = initCoords.col(7);
	pWorkVar->coordsRot.posMm.col(2) = initCoords.col(8);

	pWorkVar->coordsRot.posEa.col(0) = initCoords.col(9);
	pWorkVar->coordsRot.posEa.col(1) = initCoords.col(10);
	pWorkVar->coordsRot.posEa.col(2) = initCoords.col(11);
}

// Brief: This function initializes the integration coefficient buffers
// param[in/out]: IntCoefficient_S* pIntCoeff - pointer to a IntCoefficient_S object
// return: void
void InitIntCoeff(IntCoefficient_S* pIntCoeff)
{
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
// param[in/out]: Buffer_S* pBuffer - pointer to a Buffer_S object
// return: void
void InitBuffers(Buffer_S* pBuffer)
{
	pBuffer->buffer3d.buffer1.resize(config->getNumbPart(), 3);
	pBuffer->buffer3d.buffer2.resize(config->getNumbPart(), 3);
	pBuffer->buffer3d.buffer3.resize(config->getNumbPart(), 3);
	pBuffer->buffer3d.buffer_b.resize(config->getNumbPart(), 3);

	pBuffer->buffer1d.buffer1.resize(config->getNumbPart(), 1);
}

// Brief: Initializes Two state approximation / puts mag moment in energy mimima
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// return: void
void InitTwoStateApprox(WorkingVar_S* pWorkVar)
{
	ArrayXXd extrema;
	Array<double,1,3> vecN;

	//direction vector of magnetic flux density
	pWorkVar->coordsRot.Bvec = pWorkVar->demagFluxDens.colwise()/pWorkVar->demagFluxDens.rowwise().norm();

	//normalize vectors and project m onto the plane spanned by B and n
	VectorProjection(pWorkVar);

	pWorkVar->coordsRot.psiIs = ((pWorkVar->coordsRot.posMm*pWorkVar->coordsRot.posEa).rowwise().sum()).acos(); //angle between posEa and posMm after projection
	pWorkVar->coordsRot.phiIs = ((pWorkVar->coordsRot.posEa*pWorkVar->coordsRot.Bvec).rowwise().sum()).acos();  //angle between posEa and magnetic field

	for (int i = 0; i < config->getNumbPart(); i++)
	{
		extrema = FindLocalMinMax(pWorkVar,i);

		////////////////// put mag moments in minima //////////////////////////////
		if (extrema.cols() == 2) //check if there are two minima
		{
			if (pWorkVar->coordsRot.psiIs(i) > extrema(2,0) && pWorkVar->coordsRot.psiIs(i) < extrema(2,1))
			{
				pWorkVar->coordsRot.psiIs(i) = extrema(0,1);
				pWorkVar->coordsRot.state(i) = 2;
			}
			else
			{
				pWorkVar->coordsRot.psiIs(i) = extrema(0,0);
				pWorkVar->coordsRot.state(i) = 1;
			}
		}
		else
		{
			pWorkVar->coordsRot.psiIs(i) = extrema(0); //if there is only one minimum
			pWorkVar->coordsRot.state(i) = 0;
		} 	
		
		vecN = pWorkVar->coordsRot.Bvec.row(i) - (pWorkVar->coordsRot.Bvec.row(i)*pWorkVar->coordsRot.posEa.row(i)).sum()*pWorkVar->coordsRot.posEa.row(i);
		vecN.colwise() /=vecN.rowwise().norm();
		pWorkVar->coordsRot.posMm.row(i) = cos(pWorkVar->coordsRot.psiIs(i))*pWorkVar->coordsRot.posEa.row(i) + sin(pWorkVar->coordsRot.psiIs(i))*vecN;
	}
}

////////////////////////*simulation functions*///////////////////////////////////////

// Brief: This function performs the simulation
// return: void
void RunSimulation(void)
{
	WorkingVar_S*      pWorkVar;    // pointer to the simulation working varibales
	Buffer_S*           pBuffer;    // pointer to buffer
	IntCoefficient_S* pIntCoeff;    // pointer to the integration coefficient varibales
	OutputVar_S*     pOutputVar;    // pointer to output variables
	Params_S*           pParams;    // pointer to parameter structure

	double               tStart;    // start time to measure
	double              runTime;    // runtime in s
	double                 tInt;    // current integration time in s
	double               deltaT;    // integration timestep
	double           deltaTused;    // used integration timestep from solver

	// initialize the pointers
	pWorkVar = &workVar;
	pBuffer = &buffer;
	pIntCoeff = &intCoeff;
	pOutputVar = &outputVar;
	pParams = &params;

	// initialize time settings
	tInt = 0;              // set integration time to 0
	deltaT = config->getDeltaTInit();   // set init timestep
	int steps = 0;         // counts integration steps
	int countCoords = 0;     // progress counter for coords exportation

	// start measurement
	tStart = clock();

	// simulation loop
	while (tInt <= config->getTEnd()) //dont hit end time exactly => does not matter for my applications
	{		
		// set external magnetic field depending on time
		DefExtField(&pWorkVar->extFluxDens, tInt);

		// set thermal fluctuations
		DefThermFluct(pWorkVar, &pBuffer->buffer1d, deltaT);
		
		// compute forces
		EvalSinCos(pWorkVar, pBuffer, pParams);
		CompInteractions(pWorkVar, pBuffer, pParams, pOutputVar);

		// check if ewald summation works correct
		//CheckEwaldSum(pWorkVar, pBuffer, pParams, pOutputVar);
		
		// numerical integration routine for rotational movement
  		IntegrationRot(pWorkVar, pIntCoeff, pBuffer, pParams, &deltaT, &deltaTused); 

		//numerical integration for translational movement
		IntegrationTrans(pWorkVar, pBuffer, deltaTused, pParams);

		//two state approximation
		TwoStateApprox(pWorkVar, deltaT);

		// Apply periodic boundary conditions
		ApplyBoundaryCondArr(&pWorkVar->coordsTrans.pos, &pBuffer->buffer3d, pParams);

		///////////////////postprocessing stuff////////////////////////////
		// displaying simulation progress in % 0,1,2,3,4,5....
		if (floor(tInt / config->getTEnd() * config->getNPCoords()) == countCoords)
		{
			runTime = (clock() - tStart) / CLOCKS_PER_SEC;

			cout << "actual timestep: " << deltaT << "s" << endl;
			cout << "actual time: " << tInt << "s" << endl;
			cout << "simulation progress: " << countCoords / config->getNPCoords() * 100.0 << "%" << endl;
			cout << "simulation will end in about " << (config->getTEnd() - tInt) * (runTime / tInt) / 3600.0 << "h" << endl;

			countCoords++;

			//store coords in txt
			WriteCoords2TXT(&workVar, "coords.txt", tInt);

			//enable translational output
			pParams->transOutput = true;
		}

		EvalTransCoeffs(pWorkVar, pBuffer, pOutputVar, pParams);
	
		// magnetization vector
		pOutputVar->t.push_back(tInt);
		pOutputVar->nz.push_back(pWorkVar->coordsRot.posEa.col(2).abs().mean());
		pOutputVar->mz.push_back(pWorkVar->coordsRot.posMm.col(2).mean());

		//time counting
		tInt += deltaTused;
		steps++;
	}

	// stop measurment
	runTime = (clock() - tStart) / CLOCKS_PER_SEC;
	cout << "simulation progress: 100%" << endl;
	cout << "  simulationtime = " << runTime/3600. << " h" << " (" << runTime << " s)" << endl;
	cout << steps << " timesteps needed" << endl;

	// store end coords in txt
	WriteCoords2TXT(&workVar, "coords.txt", tInt);

	// export diffusion and viscosity data
	ExportDiffVis(pOutputVar, tInt, steps);

	//export magnetization vector
	saveData("mz.txt", Map<VectorXd, Unaligned>(pOutputVar->mz.data(), pOutputVar->mz.size()), 0);
	saveData("nz.txt", Map<VectorXd, Unaligned>(pOutputVar->nz.data(), pOutputVar->nz.size()), 0);
	saveData("t.txt", Map<VectorXd, Unaligned>(pOutputVar->t.data(), pOutputVar->t.size()), 0);

	cout << "simulation ended successfully \n" << endl;
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

// Brief: This function creates random thermal fluctuations, equations can be found in my paper
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]	: Buffer1d_S* pBuffer1d - pointer to Buffer1d_S object
// param[in]	: double actual timestep
// return: void
void DefThermFluct(WorkingVar_S* pWorkVar, Buffer1d_S* pBuffer1d, double deltaT)
{
		if (config->getEnableThermTorque()) //thermal torque
		{
			pBuffer1d->buffer1 = ((2.0 * config->getKB() * config->getTemp() / deltaT) * pWorkVar->partProp.zetaRot).cwiseSqrt();

			// normal distribution with mean = 0, stdev = 1.0
			pWorkVar->thermTorque = Rand::normal<ArrayXXd>(config->getNumbPart(), 3, generator, 0, 1.0);

			//multiplication of nx3 matrix with nx1 vector / each column of matrix elementwise with vector
			pWorkVar->thermTorque = pWorkVar->thermTorque.colwise() * pBuffer1d->buffer1;  //works like in matlab	
		}
		else
		{
			pWorkVar->thermTorque.setZero();
		}

		if (config->getEnableThermForce()) //thermal force
		{
			pBuffer1d->buffer1 = ((2.0 * config->getKB() * config->getTemp() / deltaT) * pWorkVar->partProp.zetaTrans).cwiseSqrt();

			pWorkVar->thermForce = Rand::normal<ArrayXXd>(config->getNumbPart(), 3, generator, 0, 1.0);
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

// Brief: This function computes the forces acting on the particles and computes the pressure tensor for viscosity evaluation
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]	: Buffer_S* pBuffer - pointer to Buffer_S object
// param[in]: Params_S* pParams - pointer to params_S object
// param[in/out]: OutputVar_S* pOutputVar - pointer to a OutputVar_S object
// return: void
void CompInteractions(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, Params_S* pParams, OutputVar_S* pOutputVar)
{
	pWorkVar->intForce.setZero();         //force array, initialized to zero
	pWorkVar->demagFluxDens.setZero();    //demagnetizing field array, initialized to zero

	//pressure tensor components
	pOutputVar->Pxy = 0; //set pressure tensor components zero before adding forces for viscosity
	pOutputVar->Pxz = 0;
	pOutputVar->Pyz = 0;

	pWorkVar->partProp.magMomVec = pWorkVar->coordsRot.posMm.colwise() * pWorkVar->partProp.magMom; //create magnetic moment vectors
	pWorkVar->partProp.magMomSum = pWorkVar->partProp.magMomVec.colwise().sum(); // term for surface contribution in Ewald summation

	// compute interaction vectors and create lists for long and short range interactions (used because distances are again used in partial solver steps)
	CompInteractVecs(pWorkVar, pBuffer, pParams);

	//Ewald summation, fourier and real space contributions
	CompDipInteractEwaldR(pWorkVar, pBuffer, pParams, pOutputVar);
	CompDipInteractEwaldF(pWorkVar, pBuffer, pParams);

	// loop for short range interactions
	int idx1, idx2;
	for (int i = 0; i < pWorkVar->intLists.pairs; i++) // looping through all pairs but list of pairs for short range interactions would be very small, is it better to create extra lists??
	{
		idx1 = pWorkVar->intLists.intPairs(2 * i);
		idx2 = pWorkVar->intLists.intPairs(2 * i + 1);

		// only executes interaction computation if particles are within interaction range
		if (pWorkVar->intLists.intRangeSR(i))
		{
			if (config->getEnableCoatingPot())
			{
				// set buffer vec to zero to sum up short range forces for pressure tensor calculation
				pBuffer->bufferVec.vec2.setZero();

				//steric repulsion
				CompSterRep(pWorkVar, pBuffer, idx1, idx2, i);

				// electrostatic repulsion
				CompElectrostatRep(pWorkVar, pBuffer, idx1, idx2, i);

				// vdw forces
				CompForceVdW(pWorkVar, pBuffer, idx1, idx2, i);

				// compute pressure tensor components for short range interactions
				EvalPressTensSR(pWorkVar, pOutputVar, pBuffer, i);
			}
		}
	}
	
	//add thermal fluctations to force array
	pWorkVar->intForce += pWorkVar->thermForce; 

	//add external flux dens to interacting field
	pWorkVar->demagFluxDens += pWorkVar->extFluxDens;
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

	//Ewald summation
	CompDipFieldEwaldR(pWorkVar, pBuffer, pParams);
	CompDipFieldEwaldF(pWorkVar, pBuffer, pParams);
} 

// Brief: Two state approximation for jumps of magnetic moment
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: double deltaT - current integration timestep
// return: void
void TwoStateApprox(WorkingVar_S* pWorkVar, double deltaT)
{
	ArrayXXd extrema;
	Array<double,1,3> vecN;
	double dE1, dE2, r1, r2, numb, p;
	bool jump;

	//direction vector of magnetic flux density
	pWorkVar->coordsRot.Bvec = pWorkVar->demagFluxDens.colwise()/pWorkVar->demagFluxDens.rowwise().norm();	

	//normalize vectors and project m onto the plane spanned by B and n
	VectorProjection(pWorkVar);

	pWorkVar->coordsRot.phiIs = ((pWorkVar->coordsRot.posEa*pWorkVar->coordsRot.Bvec).rowwise().sum()).acos();  //angle between posEa and magnetic field

	for (int i = 0; i < config->getNumbPart(); i++)
	{
		extrema = FindLocalMinMax(pWorkVar,i);

		////////////////// put mag moments back in minima after easy axis movement //////////////////////////////
		if (extrema.cols() == 2) //check if there are two minima
		{
			if (pWorkVar->coordsRot.state(i) == 1)
			{
				pWorkVar->coordsRot.psiIs(i) = extrema(0,0);
			}
			else if (pWorkVar->coordsRot.state(i) == 2)
			{
				pWorkVar->coordsRot.psiIs(i) = extrema(0,1);
			}
			else if (pWorkVar->coordsRot.state(i) == 0)
			{
				if (pWorkVar->coordsRot.psiIs(i) > extrema(2,0) && pWorkVar->coordsRot.psiIs(i) < extrema(2,1))
				{
					pWorkVar->coordsRot.psiIs(i) = extrema(0,1);
					pWorkVar->coordsRot.state(i) = 2;
				}
				else
				{
					pWorkVar->coordsRot.psiIs(i) = extrema(0,0);
					pWorkVar->coordsRot.state(i) = 1;
				}
			}
		}
		else
		{
			pWorkVar->coordsRot.psiIs(i) = extrema(0); //if there is only one minimum
			pWorkVar->coordsRot.state(i) = 0;
		} 	
		
		vecN = pWorkVar->coordsRot.Bvec.row(i) - (pWorkVar->coordsRot.Bvec.row(i)*pWorkVar->coordsRot.posEa.row(i)).sum()*pWorkVar->coordsRot.posEa.row(i);
		vecN.colwise() /=vecN.rowwise().norm();
		pWorkVar->coordsRot.posMm.row(i) = cos(pWorkVar->coordsRot.psiIs(i))*pWorkVar->coordsRot.posEa.row(i) + sin(pWorkVar->coordsRot.psiIs(i))*vecN;

		////////////////calculate possible jumps of mag moments ////////////////////////////
		if (extrema.cols() == 2)
		{
			//ernergy differences
			dE1 = extrema.row(3).minCoeff() - extrema(1,0);
			dE2 = extrema.row(3).minCoeff() - extrema(1,1);
			r1 = 1.0/(2.0*pWorkVar->partProp.tau0(i))*exp(-dE1/(config->getKB()*config->getTemp())); //switching rate from minimum 1 to minimum 2
			r2 = 1.0/(2.0*pWorkVar->partProp.tau0(i))*exp(-dE2/(config->getKB()*config->getTemp())); //switching rate from minimum 2 to minimum 1

			numb = ((double) rand() / (RAND_MAX + 1.0)); //random number between 0 and 1
			jump = false;

			if (pWorkVar->coordsRot.state(i) == 1)
			{
				p = r1/(r1+r2)*(1-exp(-(r1+r2)*deltaT));

				if (p > numb)
				{
					pWorkVar->coordsRot.state(i) = 2;
					pWorkVar->coordsRot.psiIs(i) = extrema(0,1);

					vecN = pWorkVar->coordsRot.Bvec.row(i) - (pWorkVar->coordsRot.Bvec.row(i)*pWorkVar->coordsRot.posEa.row(i)).sum()*pWorkVar->coordsRot.posEa.row(i);
					vecN.colwise() /=vecN.rowwise().norm();
					pWorkVar->coordsRot.posMm.row(i) = cos(pWorkVar->coordsRot.psiIs(i))*pWorkVar->coordsRot.posEa.row(i) + sin(pWorkVar->coordsRot.psiIs(i))*vecN;

					jump = true;
				}
			}
			if (pWorkVar->coordsRot.state(i) == 2 && jump == false)
			{
				p = r2/(r1+r2)*(1-exp(-(r1+r2)*deltaT));

				if (p > numb)
				{
					pWorkVar->coordsRot.state(i) = 1;
					pWorkVar->coordsRot.psiIs(i) = extrema(0,0);

					vecN = pWorkVar->coordsRot.Bvec.row(i) - (pWorkVar->coordsRot.Bvec.row(i)*pWorkVar->coordsRot.posEa.row(i)).sum()*pWorkVar->coordsRot.posEa.row(i);
					vecN.colwise() /=vecN.rowwise().norm();
					pWorkVar->coordsRot.posMm.row(i) = cos(pWorkVar->coordsRot.psiIs(i))*pWorkVar->coordsRot.posEa.row(i) + sin(pWorkVar->coordsRot.psiIs(i))*vecN;
				}
			}
		} 
	}	
}

// Brief: This function performs the numerical integration of the rotational equations of motion with the selected integration method
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out]: IntCoefficient_S* pIntCoeff - pointer to a IntCoefficient_S object
// param[in]    : Buffer_S* pBuffer - pointer to Buffer_S* object 
// param[in]	: double* pDeltaT - timestep for integration (adjusted by "AdjustTimestep")
// param[out]	: double* pDeltaTused - pointer to the actual used timestep by the integration routine
// return: void
void IntegrationRot(WorkingVar_S* pWorkVar, IntCoefficient_S* pIntCoeff, Buffer_S* pBuffer, Params_S* pParams, double* pDeltaT, double* pDeltaTused)
{
	bool    cond;       // condition for while loop, set by AdjustTimestep
	bool    noFailed;   // condition for while loop, set by AdjustTimestep
	double  errorMm;    // error for magnetic moment
	double 	deltaT;		
	double  power;      // error convergence

	// set initial parameters
	cond = true;
	noFailed = true;

	switch (config->getSolver())
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
			RKCoeffs(pWorkVar, &pIntCoeff->coords2, &pIntCoeff->k2, pBuffer);

			// error calculation
			pBuffer->buffer3d.buffer1 = 0.5 * pIntCoeff->k1.knEa - 0.5 * pIntCoeff->k2.knEa; //error of solver
			errorMm = deltaT * (pBuffer->buffer3d.buffer1.rowwise().norm()).maxCoeff(); //same as normcontrol in matlab

			// save old deltaT and adjust timestep
			*pDeltaTused = deltaT;		
			AdjustTimeStep(pWorkVar, errorMm, power, pDeltaT, &cond, &noFailed);			
		}

		pWorkVar->coordsRot.posEa += *pDeltaTused * (0.5 * pIntCoeff->k1.knEa + 0.5 * pIntCoeff->k2.knEa); //2.order solution
		break;
	}
	case BOGACKI_SHAMPINE:
	{
		power = 1.0 / 3.0;

		while (cond)
		{
			// pre load deltaT for faster execution
			deltaT = (*pDeltaT);

			RKCoeffs(pWorkVar, &pWorkVar->coordsRot, &pIntCoeff->k1, pBuffer);

			pIntCoeff->coords2.posEa = pWorkVar->coordsRot.posEa + deltaT * 1.0 / 2.0 * pIntCoeff->k1.knEa;
			RKCoeffs(pWorkVar, &pIntCoeff->coords2, &pIntCoeff->k2, pBuffer);

			pIntCoeff->coords3.posEa = pWorkVar->coordsRot.posEa + deltaT *  3.0 / 4.0 * pIntCoeff->k2.knEa;
			RKCoeffs(pWorkVar, &pIntCoeff->coords3, &pIntCoeff->k3, pBuffer);

			pIntCoeff->coords4.posEa = pWorkVar->coordsRot.posEa + deltaT * (2.0 / 9.0 * pIntCoeff->k1.knEa + 1.0 / 3.0 * pIntCoeff->k2.knEa + 4.0 / 9.0 * pIntCoeff->k3.knEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords4, &pIntCoeff->k4, pBuffer);

			pBuffer->buffer3d.buffer1 = -5.0 / 72.0 * pIntCoeff->k1.knEa + 1.0 / 12.0 * pIntCoeff->k2.knEa + 1.0 / 9.0 * pIntCoeff->k3.knEa - 1.0 / 8.0 * pIntCoeff->k4.knEa;
			errorMm = deltaT * (pBuffer->buffer3d.buffer1.rowwise().norm()).maxCoeff();

			// save old deltaT and adjust timestep
			*pDeltaTused = deltaT;
			AdjustTimeStep(pWorkVar, errorMm, power, pDeltaT, &cond, &noFailed);
		}

		pWorkVar->coordsRot.posEa = pIntCoeff->coords4.posEa;
		break;
	}
	case DORMAND_PRINCE:
	{
		power = 1.0 / 5.0;

		while (cond)
		{
			// pre load deltaT for faster execution
			deltaT = (*pDeltaT);

			RKCoeffs(pWorkVar, &pWorkVar->coordsRot, &pIntCoeff->k1, pBuffer);

			pIntCoeff->coords2.posEa = pWorkVar->coordsRot.posEa + deltaT * 1.0 / 5.0 * pIntCoeff->k1.knEa;
			RKCoeffs(pWorkVar, &pIntCoeff->coords2, &pIntCoeff->k2, pBuffer);

			pIntCoeff->coords3.posEa = pWorkVar->coordsRot.posEa + deltaT * (3.0 / 40.0 * pIntCoeff->k1.knEa + 9.0 / 40.0 * pIntCoeff->k2.knEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords3, &pIntCoeff->k3, pBuffer);

			pIntCoeff->coords4.posEa = pWorkVar->coordsRot.posEa + deltaT * (44.0 / 45.0 * pIntCoeff->k1.knEa - 56.0 / 15.0 * pIntCoeff->k2.knEa + 32.0 / 9.0 * pIntCoeff->k3.knEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords4, &pIntCoeff->k4, pBuffer);

			pIntCoeff->coords5.posEa = pWorkVar->coordsRot.posEa + deltaT * (19372.0 / 6561.0 * pIntCoeff->k1.knEa - 25360.0 / 2187.0 * pIntCoeff->k2.knEa + 64448.0 / 6561.0 * pIntCoeff->k3.knEa - 212.0 / 729.0 * pIntCoeff->k4.knEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords5, &pIntCoeff->k5, pBuffer);

			pIntCoeff->coords6.posEa = pWorkVar->coordsRot.posEa + deltaT * (9017.0 / 3168.0 * pIntCoeff->k1.knEa - 355.0 / 33.0 * pIntCoeff->k2.knEa + 46732.0 / 5247.0 * pIntCoeff->k3.knEa + 49.0 / 176.0 * pIntCoeff->k4.knEa - 5103.0 / 18656.0 * pIntCoeff->k5.knEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords6, &pIntCoeff->k6, pBuffer);

			pIntCoeff->coords7.posEa = pWorkVar->coordsRot.posEa + deltaT * (35.0 / 384.0 * pIntCoeff->k1.knEa + 500.0 / 1113.0 * pIntCoeff->k3.knEa + 125.0 / 192.0 * pIntCoeff->k4.knEa - 2187.0 / 6784.0 * pIntCoeff->k5.knEa + 11.0 / 84.0 * pIntCoeff->k6.knEa);
			RKCoeffs(pWorkVar, &pIntCoeff->coords7, &pIntCoeff->k7, pBuffer);

			// error calculation
			pBuffer->buffer3d.buffer1 = 71.0 / 57600.0 * pIntCoeff->k1.knEa - 71.0 / 16695.0 * pIntCoeff->k3.knEa + 71.0 / 1920.0 * pIntCoeff->k4.knEa - 17253.0 / 339200.0 * pIntCoeff->k5.knEa + 22.0 / 525.0 * pIntCoeff->k6.knEa - 1.0 / 40.0 * pIntCoeff->k7.knEa;
			errorMm = deltaT * (pBuffer->buffer3d.buffer1.rowwise().norm()).maxCoeff();
			
			// save old deltaT and adjust timestep
			*pDeltaTused = deltaT;
			AdjustTimeStep(pWorkVar, errorMm, power, pDeltaT, &cond, &noFailed);
		}

		pWorkVar->coordsRot.posEa += *pDeltaTused * (35.0 / 384.0 * pIntCoeff->k1.knEa + 500.0 / 1113.0 * pIntCoeff->k3.knEa + 125.0 / 192.0 * pIntCoeff->k4.knEa - 2187.0 / 6784 * pIntCoeff->k5.knEa + 11.0 / 84.0 * pIntCoeff->k6.knEa); //5.order solution
		break;
	}
	}

	// make sure that the length of the unit vectors stay 1
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

	pBuffer->buffer1d.buffer1 = (pWorkVar->coordsRot.posMm * pCoord->posEa).rowwise().sum(); //dot Product dot(posMm,posEa)
	
	// check if particles can rotate free
	if (config->getEnableMobilization())
	{
		//matlab code knEa = velEaConst1.*dot(posMm,posEa).*(posMm - dot(posMm,posEa).*posEa) + velEaConst2.*cross(torqueTherm, posEa);

		RowWiseCrossProd(&pWorkVar->thermTorque, &pCoord->posEa, &pBuffer->buffer3d.buffer2);

		pRkCoeffs->knEa = (pWorkVar->coordsRot.posMm - pCoord->posEa.colwise() * pBuffer->buffer1d.buffer1).colwise() * (pWorkVar->partProp.velEaConst1 * pBuffer->buffer1d.buffer1) + pBuffer->buffer3d.buffer2.colwise() * pWorkVar->partProp.velEaConst2;
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
// param[out]: double* pTDelta - new time delate
// param[out]: bool* pCond - new loop condidtion
// param[out]: bool* pNoFailed - new failed condidtion
// return: void
void AdjustTimeStep(WorkingVar_S* pWorkVar, double error, double power, double* pDeltaT, bool* pCond, bool* pNoFailed)
{
	double deltaTnew;

	// norm(e(i)) <= max(RelTol*norm(y(i)),AbsTol(i)) from matlab, norm(y(i)) is 1 in the case of unit vectors
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

		// rescale thermal fluctuations
		pWorkVar->thermTorque *= sqrt(deltaTnew / (*pDeltaT));
		pWorkVar->thermForce *= sqrt(deltaTnew / (*pDeltaT));
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

// Brief: Integration of translational movement
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: double deltaT - timestep used in the rotational integration
// return: void
void IntegrationTrans(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, double deltaTused, Params_S* pParams)
{
	pWorkVar->coordsTrans.pos += deltaTused * (pWorkVar->intForce.colwise() / pWorkVar->partProp.zetaTrans);

	if (pParams->transOutput)
	{
		pParams->transOutput = false;
		pBuffer->buffer3d.buffer2 = deltaTused * (pWorkVar->intForce.colwise() / pWorkVar->partProp.zetaTrans);
		cout << "maximum displacement:" << (pBuffer->buffer3d.buffer2.rowwise().norm()/pParams->lBox).maxCoeff() << " of lBox and " << (pBuffer->buffer3d.buffer2.rowwise().norm()/config->getRMagMean()).maxCoeff()<< " of rMag \n"<< endl;
	}
}