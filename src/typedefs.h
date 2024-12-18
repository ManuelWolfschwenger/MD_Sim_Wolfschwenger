
/*
* Name : typedefs.h
* Auto :
* Brief: typedefs
*/

#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_

//external libraries
#include <Eigen/Dense>

// namespace definition
using namespace Eigen;
using namespace std;

/* structure definitions */

// particale properties structure
typedef struct
{
	ArrayXd 		rHydr;
	ArrayXd 		rMag;
	ArrayXd 		volMag;
	ArrayXd 		volHydr;
	ArrayXd 		magMom;
	ArrayXd 		zetaRot;
	ArrayXd 		velEaConst1;
	ArrayXd 		velEaConst2;
	ArrayXd 		tau0;
	ArrayXd 		psi;
	ArrayXd         sinPsi;
	ArrayXd         cosPsi;
	ArrayXd         sinPsi2;
	ArrayXd         dPotEn;
	ArrayXd         ddPotEn;
} PartProps_S;

// particale coordinates structure
typedef struct
{
	ArrayXXd 		posMm;
	ArrayXXd 		posEa;
	ArrayXXd        vecNormal;
	ArrayXXd        mProj;
	ArrayXXd        Bvec;
	ArrayXd         psiIs;
	ArrayXd         phiIs;
	ArrayXi         state;
	ArrayXXd 		omegaEa;
} Coord_S;

typedef struct
{
	ArrayXXd 		knEa;
	ArrayXXd 		knEaAss;
} RkCoeff_S;

// integration coefficant structure
typedef struct
{
	RkCoeff_S 		k1;
	RkCoeff_S 		k2;
	RkCoeff_S 		k3;
	RkCoeff_S 		k4;
	RkCoeff_S 		k5;
	RkCoeff_S 		k6;
	RkCoeff_S 		k7;
	Coord_S 		coords2;
	Coord_S 		coords3;
	Coord_S 		coords4;
	Coord_S 		coords5;
	Coord_S 		coords6;
	Coord_S 		coords7;
} IntCoefficient_S;

// temporal 3 dimensional buffers
typedef struct
{
	ArrayXXd 		buffer1;		//	3 dimensional buffer
	ArrayXXd 		buffer2;		//	3 dimensional buffer
	ArrayXXd 		buffer3;		//	3 dimensional buffer
} Buffer3d_S;

// temporal 1 dimensional buffers
typedef struct
{
	ArrayXd 		buffer1;		//	1 dimensional buffer
} Buffer1d_S;

// working varibales structure
typedef struct
{
	PartProps_S 	partProp;   	// PartProps_S object
	Coord_S     	coords;     	// Coord_S object
	ArrayXXd 		extFluxDens;	// 
	ArrayXXd 		thermTorque;	// 
	Buffer3d_S		buffer3d;		// 3 dimensional buffers
	Buffer1d_S		buffer1d;		// 1 dimensional buffers
} WorkingVar_S;

typedef struct
{
	vector<double>  t;         //time vector
	vector<double>  mZ;      //magnetization vector
	vector<double>  nZ;      //magnetization vector
	vector<double> omegaX, omegaY, omegaZ;
}OutputVar_S;

#endif 