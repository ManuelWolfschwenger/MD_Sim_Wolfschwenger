/*
* Name : postprocessing.cpp
* Auto :
* Brief: Simulation source file
*/

// own header file
#include "postprocessing.h"
#include "configuration.h"

// Brief: This function writes some data from a matrix into a txt file and stores it
// param[in]: string filename - name of file
// param[in]: MatrixXd matrix - matrix to store in file
// return: void
void saveData(string filename, MatrixXd  matrix, double tInt)
{
	//const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

	ofstream file(filename, ios::app); //ios::app...seek to the end of stream before each write

	if (file.is_open())
	{
		if (tInt == 0) //to decide wheather the file content is removed are just added on
		{
			ofstream file(filename, std::ios::out | std::ios::trunc); //out...allows output (writing operations), ofstream should have this automatically set, trunc...old file contents are removed
		}

		//file << matrix.format(CSVFormat);
		file << matrix << "\n" << endl;
		file.close();
	}
}

// Brief: This function reads some data from a txt file and stores it in a Matrix
// param[in]: string fileToOpen - name of file
// return: MatrixXd matrix - matrix to store in file
MatrixXd openData(string fileToOpen)
{
	// the matrix entries are stored in this variable row-wise. For example if we have the matrix:
	// M=[a b c 
	//    d e f]
	// the entries are stored as matrixEntries=[a,b,c,d,e,f], that is the variable "matrixEntries" is a row vector
	// later on, this vector is mapped into the Eigen matrix format
	vector<double> matrixEntries;

	// in this object we store the data from the matrix
	ifstream matrixDataFile(fileToOpen);

	// this variable is used to store the row of the matrix that contains commas 
	string matrixRowString;

	// this variable is used to store the matrix entry;
	string matrixEntry;

	// this variable is used to track the number of rows
	int matrixRowNumber = 0;


	while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and store every line into the string variable matrixRowString
	{
		stringstream matrixRowStringStream(matrixRowString); //convert matrixRowString that is a string to a stream variable.

		while (matrixRowStringStream >> matrixEntry) //while (getline(matrixRowStringStream, matrixEntry, ',')) // here we read pieces of the stream matrixRowStringStream until every comma, and store the resulting character into the matrixEntry
		{
			matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector storing all the matrix entries
		}
		matrixRowNumber++; //update the column numbers
	}
	matrixRowNumber -= 1; //because one empty line is inserted

	// here we convet the vector variable into the matrix and return the resulting object, 
	// note that matrixEntries.data() is the pointer to the first memory location at which the entries of the vector matrixEntries are stored;
	return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);

}

// Brief: This function writes the coords to a txt-file and saves it in the working directory
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// return: void
void WriteCoords2TXT(WorkingVar_S* pWorkVar, string filename, double tInt)
{
	MatrixXd res(config->getNumbPart(), 10);

	res.col(0) = pWorkVar->coordsTrans.pos.col(0); //x
	res.col(1) = pWorkVar->coordsTrans.pos.col(1); //y
	res.col(2) = pWorkVar->coordsTrans.pos.col(2); //z

	res.col(3) = pWorkVar->coordsRot.posMm.col(0); //posMm x
	res.col(4) = pWorkVar->coordsRot.posMm.col(1); //posMm y
	res.col(5) = pWorkVar->coordsRot.posMm.col(2); //posMm z

	res.col(6) = pWorkVar->coordsRot.posEa.col(0);  //posEa x
	res.col(7) = pWorkVar->coordsRot.posEa.col(1); //posEa y
	res.col(8) = pWorkVar->coordsRot.posEa.col(2); //posEa z

	res.col(9).setZero();
	res(0, 9) = tInt;

	saveData(filename, res, tInt);
}

// Brief: This function writes data to a txt-file and saves it in the working directory
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// return: void
void WriteData2TXT(WorkingVar_S* pWorkVar, Params_S* pParams, string filename)
{
	MatrixXd res(config->getNumbPart(), 3);

	res.col(0) = pWorkVar->partProp.rHydr; //hydrodyn radius
	res.col(1) = pWorkVar->partProp.rMag;  //core radius

	res.col(2).setZero();
	res(0, 2) = pParams->lBox;
	res(1, 2) = config->getNumbPart();
	res(2, 2) = config->getMagFluxDens();
	res(3, 2) = config->getSatMag();
	res(4, 2) = config->getTemp();
	res(5, 2) = config->getTMag();
	res(6, 2) = config->getTRelax();
	res(7, 2) = config->getTEnd();
	res(8, 2) = config->getVolFrac();
	res(9, 2) = config->getInitConfigRot();
	res(10, 2) = config->getInitConfigTrans();
	res(11, 2) = config->getSeed();
	res(12, 2) = config->getAnisEn();
	res(13, 2) = config->getDensSurfMol();
	res(14, 2) = config->getEffLenSurfMol();
	res(15, 2) = config->getHamaker();
	res(16, 2) = config->getZValency();
	res(17, 2) = config->getEpsilonR();
	res(18, 2) = config->getGamma0();
	res(19, 2) = config->getAbsTol();
	res(20, 2) = config->getErrTolEwald();
	res(21, 2) = config->getErrTolSR();

	saveData(filename, res, 0);
}

// Brief: This function measures time that a function takes to execute
// param[in/out]: double* pTime, pointer to time variable
// param[in]:     bool start, begin or end measurement 
// param[in]:     string fun, function name
// return: void
void TimeMeasure(double* pTime, bool start, string fun)
{
	if (start)
	{
		*pTime = clock();
	}
	else
	{
		*pTime = (clock() - *pTime) / CLOCKS_PER_SEC;	
		cout << fun << ":   " << *pTime << "s" << endl;
	}
}



/////////////eval transport coefficients/////////////////////////////

// Brief: This function evaluates the diffusion coefficient and the viscosity via Einstein and GK relation
// param[in]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out]: OutputVar_S* pOutputVar - pointer to a OutputVar_S object
// param[in]: Params_S* pParams - pointer to params_S object
// param[in]: Buffer_S* pBuffer - pointer to Buffer_S* object 
// return: void
void EvalTransCoeffs(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, OutputVar_S* pOutputVar, Params_S* pParams)
{
	// Einstein diffusion
	//remove wrap around effects from particle position => Rapaport S122
	pOutputVar->truePos = pWorkVar->coordsTrans.pos + ((pOutputVar->truePos - pWorkVar->coordsTrans.pos) / pParams->lBox).round() * pParams->lBox;

	pBuffer->buffer3d.buffer1 = pOutputVar->truePos - pOutputVar->posRef;
	pOutputVar->MSDx.push_back((pBuffer->buffer3d.buffer1.col(0) * pBuffer->buffer3d.buffer1.col(0)).mean());
	pOutputVar->MSDy.push_back((pBuffer->buffer3d.buffer1.col(1) * pBuffer->buffer3d.buffer1.col(1)).mean());
	pOutputVar->MSDz.push_back((pBuffer->buffer3d.buffer1.col(2) * pBuffer->buffer3d.buffer1.col(2)).mean());

	//GK viscosity
	pOutputVar->Pxy /= pParams->volBox;
	pOutputVar->Pxz /= pParams->volBox;
	pOutputVar->Pyz /= pParams->volBox;

	//Einstein viscosity
	pOutputVar->PxyVec.push_back(pOutputVar->Pxy);
	pOutputVar->PxzVec.push_back(pOutputVar->Pxz);
	pOutputVar->PyzVec.push_back(pOutputVar->Pyz);
}

// Brief: This function calculates the pressure tensor components for long range interactions
// param[in]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in/out]: OutputVar_S* pOutputVar - pointer to a OutputVar_S object
// param[in]	: Buffer_S* pBuffer - pointer to Buffer_S object
// param[in]    : int i - index of interaction pair
// return: void
void EvalPressTensLR(WorkingVar_S* pWorkVar, OutputVar_S* pOutputVar, Buffer_S* pBuffer, int i)
{
	pOutputVar->Pxy += pWorkVar->intLists.intVecsLen(i) * pWorkVar->intLists.intVecsCoords.row(i)(0) * pBuffer->bufferVec.vec1(1);
	pOutputVar->Pxz += pWorkVar->intLists.intVecsLen(i) * pWorkVar->intLists.intVecsCoords.row(i)(0) * pBuffer->bufferVec.vec1(2);
	pOutputVar->Pyz += pWorkVar->intLists.intVecsLen(i) * pWorkVar->intLists.intVecsCoords.row(i)(1) * pBuffer->bufferVec.vec1(2);
}

// Brief: This function calculates the pressure tensor components for short range interactions
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in] : Buffer_S* pBuffer - pointer to buffer_S object
// param[in]: int idx1, int idx2, particle indizes
// param[in]: int i, pair indizes
// return: void
void EvalPressTensSR(WorkingVar_S* pWorkVar, OutputVar_S* pOutputVar, Buffer_S* pBuffer, int i)
{
	pOutputVar->Pxy += pWorkVar->intLists.intVecsLen(i) * pWorkVar->intLists.intVecsCoords.row(i)(0) * pBuffer->bufferVec.vec2(1);
	pOutputVar->Pxz += pWorkVar->intLists.intVecsLen(i) * pWorkVar->intLists.intVecsCoords.row(i)(0) * pBuffer->bufferVec.vec2(2);
	pOutputVar->Pyz += pWorkVar->intLists.intVecsLen(i) * pWorkVar->intLists.intVecsCoords.row(i)(1) * pBuffer->bufferVec.vec2(2);
}

// Brief: This function calculates the mean of the acf and resizes the arrays to export them as txt
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// param[in]	: Buffer_S* pBuffer - pointer to Buffer_S object
// return: void
void ExportDiffVis(OutputVar_S* pOutputVar, double tInt, int steps) //cut off zeros
{
	// mean timestep
	Array<double, 1, 1> deltaTmean;
	deltaTmean(0,0) = tInt / steps;

	saveData("dtMean.txt", deltaTmean, 0);
	saveData("Pxy.txt", Map<VectorXd, Unaligned>(pOutputVar->PxyVec.data(), pOutputVar->PxyVec.size()), 0);
	saveData("Pxz.txt", Map<VectorXd, Unaligned>(pOutputVar->PxzVec.data(), pOutputVar->PxzVec.size()), 0);
	saveData("Pyz.txt", Map<VectorXd, Unaligned>(pOutputVar->PyzVec.data(), pOutputVar->PyzVec.size()), 0);
	saveData("MSDx.txt", Map<VectorXd, Unaligned>(pOutputVar->MSDx.data(), pOutputVar->MSDx.size()), 0);
	saveData("MSDy.txt", Map<VectorXd, Unaligned>(pOutputVar->MSDy.data(), pOutputVar->MSDy.size()), 0);
	saveData("MSDz.txt", Map<VectorXd, Unaligned>(pOutputVar->MSDz.data(), pOutputVar->MSDz.size()), 0);
	saveData("MSDrot.txt", Map<VectorXd, Unaligned>(pOutputVar->MSDrot.data(), pOutputVar->MSDrot.size()), 0);
}