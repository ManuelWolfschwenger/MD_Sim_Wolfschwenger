/*
* Name : postprocessing.h
* Auto :
* Brief: contains postprocessing functions
*/

#ifndef POSTPROCESSING_H_
#define POSTPROCESSING_H_

// external libaries
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

//internal libaries
#include "typedefs.h"
#include "constants.h"
#include "settings.h"

// namespace definition
using namespace Eigen;
using namespace std;

// Brief: This function writes some data from a matrix into a csv file and stores it
// param[in]: string filename - name of file
// param[in]: MatrixXd matrix - matrix to store in file
// return: void
void saveData(string fileName, MatrixXd  matrix, double tInt);

// Brief: This function reads some data from a csv file and stores it in Matrix
// param[in]: string fileToOpen - name of file
// return: MatrixXd matrix - matrix to store in file
MatrixXd openData(string fileToOpen);

// Brief: This function writes the coords to a txt-file and saves it in the working directory
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// return: void
void WriteCoords2TXT(WorkingVar_S* pWorkVar, string filename, double tInt);

// Brief: This function writes data to a txt-file and saves it in the working directory
// param[in/out]: WorkingVar_S* pWorkVar - pointer to a WorkingVar_S object
// return: void
void WriteData2TXT(WorkingVar_S* pWorkVar, Params_S* pParams, string filename);

// Brief: This function measures time that a function takes to execute
// param[in/out]: double* pTime, pointer to time variable
// param[in]:     bool start, begin or end measurement 
// param[in]:     string fun, function name
// return: void
void TimeMeasure(double* pTime, bool start, string fun);

#endif  /* POSTPROCESSING_H_ */

