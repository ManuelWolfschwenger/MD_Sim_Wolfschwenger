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

// namespace definition
using namespace Eigen;
using namespace std;


void saveData(string fileName, MatrixXd  matrix, double tInt);

MatrixXd openData(string fileToOpen);

void WriteCoords2TXT(WorkingVar_S* pWorkVar, string filename, double tInt);

void WriteData2TXT(WorkingVar_S* pWorkVar, Params_S* pParams, string filename);

void TimeMeasure(double* pTime, bool start, string fun);

/////////////eval transport coefficients/////////////////////////////

void EvalTransCoeffs(WorkingVar_S* pWorkVar, Buffer_S* pBuffer, OutputVar_S* pOutputVar, Params_S* pParams);

void EvalPressTensLR(WorkingVar_S* pWorkVar, OutputVar_S* pOutputVar, Buffer_S* pBuffer, int i);

void EvalPressTensSR(WorkingVar_S* pWorkVar, OutputVar_S* pOutputVar, Buffer_S* pBuffer, int i);

void ExportDiffVis(OutputVar_S* pOutputVar, double tInt, int steps);

#endif  /* POSTPROCESSING_H_ */

