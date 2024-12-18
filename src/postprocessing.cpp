/*
* Name : postprocessing.cpp
* Auto :
* Brief: Simulation source file
*/

// own header file
#include "postprocessing.h"

// Brief: This function writes some data from a matrix into a txt file
// param[in]: string filename - name of file
// param[in]: MatrixXd matrix - matrix to store in file
// return: void
void saveData(string filename, MatrixXd  matrix, double tInt)
{
	//const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

	ofstream file(filename, ios::app); //ios::app...seek to the end of stream before each write

	if (file.is_open())
	{
		if (tInt == 0)
		{
			ofstream file(filename, std::ios::out | std::ios::trunc); //out...allows output (writing operations), ofstream should have this automatically set, trunc...old file contents are removed
		}

		//file << matrix.format(CSVFormat);
		file << matrix << endl;
		file.close();
	}
}

// Brief: This function writes data to a txt-file and saves it in the working directory
// param[in]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// return: void
void WriteCoords2TXT(WorkingVar_S* pWorkVar, string filename, double tInt)
{
	MatrixXd res(config->getNumbPart(), 7);

	res.col(0) = pWorkVar->coords.posMm.col(0); //posMm x
	res.col(1) = pWorkVar->coords.posMm.col(1); //posMm y
	res.col(2) = pWorkVar->coords.posMm.col(2); //posMm z

	res.col(3) = pWorkVar->coords.posEa.col(0); //posEa x
	res.col(4) = pWorkVar->coords.posEa.col(1); //posEa y
	res.col(5) = pWorkVar->coords.posEa.col(2); //posEa z

	res.col(6).setZero();
	res(0, 6) = tInt;

	saveData(filename, res, tInt);
}

// Brief: This function writes data to a txt-file and saves it in the working directory
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// return: void
void WriteData2TXT(WorkingVar_S* pWorkVar, string filename)
{
	ArrayXd res(10, 1);

	res.setZero();
	res(0) = config->getAnisEn();
	res(1) = config->getNumbPart();
	res(2) = config->getMagFluxDens();
	res(3) = config->getSatMag();
	res(4) = config->getTemp();
	res(5) = config->getRMagMean();

	saveData(filename, res, 0);
}