/*
* Name : postprocessing.cpp
* Auto :
* Brief: Simulation source file
*/

// own header file
#include "postprocessing.h"

// Brief: This function writes some data from a matrix into a csv file and stores it
// param[in]: string filename - name of file
// param[in]: MatrixXd matrix - matrix to store in file
// return: void
void saveData(string filename, MatrixXd  matrix, double tInt)
{
	//const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

	ofstream file(filename, ios::app);

	if (file.is_open())
	{
		if (tInt == 0)
		{
			ofstream file(filename, std::ios::out | std::ios::trunc);
		}

		//file << matrix.format(CSVFormat);
		file << matrix << "\n" << endl;
		file.close();
	}
}

// Brief: This function reads some data from a csv file and stores it in Matrix
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
	MatrixXd res(numbPart, 13);

	res.col(0) = pWorkVar->coordsTrans.pos.col(0); //x
	res.col(1) = pWorkVar->coordsTrans.pos.col(1); //y
	res.col(2) = pWorkVar->coordsTrans.pos.col(2); //z

	res.col(3) = pWorkVar->coordsTrans.vel.col(0); //vx
	res.col(4) = pWorkVar->coordsTrans.vel.col(1); //vy
	res.col(5) = pWorkVar->coordsTrans.vel.col(2); //vz

	res.col(6) = pWorkVar->coordsRot.posMm.col(0); //posMm x
	res.col(7) = pWorkVar->coordsRot.posMm.col(1); //posMm y
	res.col(8) = pWorkVar->coordsRot.posMm.col(2); //posMm z

	res.col(9) = pWorkVar->coordsRot.posEa.col(0);  //posEa x
	res.col(10) = pWorkVar->coordsRot.posEa.col(1); //posEa y
	res.col(11) = pWorkVar->coordsRot.posEa.col(2); //posEa z

	res.col(12).setZero();
	res(0, 12) = tInt;

	saveData(filename, res, tInt);
}

// Brief: This function writes data to a txt-file and saves it in the working directory
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// return: void
void WriteData2TXT(WorkingVar_S* pWorkVar, Params_S* pParams, string filename)
{
	MatrixXd res(numbPart, 3);

	res.col(0) = pWorkVar->partProp.rHydr; //hydrodyn radius
	res.col(1) = pWorkVar->partProp.rMag;  //core radius

	res.col(2).setZero();
	res(0, 2) = pParams->lBox;
	res(1, 2) = numbPart;
	res(2, 2) = magFluxDens;
	res(3, 2) = satMag;
	res(4, 2) = temp;

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