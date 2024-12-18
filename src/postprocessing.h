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

//internal libaries
#include "typedefs.h"
#include "configuration.h"

// namespace definition
using namespace Eigen;
using namespace std;

void saveData(string filename, MatrixXd  matrix, double tInt);

void WriteData2TXT(WorkingVar_S* pWorkVar, string filename);

void WriteCoords2TXT(WorkingVar_S* pWorkVar, string filename, double tInt);

#endif  /* POSTPROCESSING_H_ */
