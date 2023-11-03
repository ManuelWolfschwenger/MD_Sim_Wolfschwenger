#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_

#include <Eigen/Dense>
#include <vector>

//intern header files
#include "constants.h"

// namespace definition
using namespace Eigen;
using namespace std;

// typedef for boolean array
typedef Array<bool, Dynamic, 3> ArrayXXb;
typedef Array<bool, Dynamic, 1> ArrayXb;

/* structure definitions */

// particle properties structure
typedef struct
{
	ArrayXd 		rHydr;
	ArrayXd 		rMag;
	ArrayXd 		volMag;
	ArrayXd 		volHydr;
	ArrayXd 		magMom;
	ArrayXd 		zetaRot;
	ArrayXd         zetaTrans;
	ArrayXd 		velEaConst1;
	ArrayXd 		velEaConst2;
	ArrayXXd        magMomVec;
	ArrayXXd        magMomSum;
} PartProps_S;

// particle rotational coordinates structure
typedef struct
{
	ArrayXXd 		posMm;
	ArrayXXd 		posEa;
} CoordsRot_S;

// particale translational coordinates structure
typedef struct
{
	ArrayXXd 		pos;
	ArrayXXd 		vel;
} CoordsTrans_S;

typedef struct
{
	ArrayXXd 		knMm;
	ArrayXXd 		knEa;
} RkCoeff_S;

// integration coefficient structure
typedef struct
{
	RkCoeff_S 		k1;
	RkCoeff_S 		k2;
	RkCoeff_S 		k3;
	RkCoeff_S 		k4;
	RkCoeff_S 		k5;
	RkCoeff_S 		k6;
	RkCoeff_S 		k7;
	CoordsRot_S 		coords2;
	CoordsRot_S 		coords3;
	CoordsRot_S 		coords4;
	CoordsRot_S 		coords5;
	CoordsRot_S 		coords6;
	CoordsRot_S 		coords7;
} IntCoefficant_S;

// temporal 3 dimensional buffers
typedef struct
{
	ArrayXXd 		buffer1;		//	3 dimensional buffer
	ArrayXXd 		buffer2;		//	3 dimensional buffer
	ArrayXXd 		buffer3;		//	3 dimensional buffer
	ArrayXXb 		buffer_b;		//	3 dimensional buffer boolean
} Buffer3d_S;

// temporal 1 dimensional buffers
typedef struct
{
	ArrayXd 		buffer1;		//	1 dimensional buffer
} Buffer1d_S;

//temporal vector buffer
typedef struct
{
	Array<double, 1, 3> vec1; 
	Array<double, 1, 3> vec2;
	Array<double, 1, 3> vec3;
	Array<double, 1, 3> vec_r; //vector for interaction radi      
	Array<bool, 1, 3> vec_b;   //bool vector
}BufferVec_S;

typedef struct
{
	Buffer3d_S		buffer3d;		// 3 dimensional buffers
	Buffer1d_S		buffer1d;		// 1 dimensional buffers
	BufferVec_S     bufferVec;      // Vector buffer
	Array<Array<double, 1, 3>, Dynamic, Dynamic> tCos, tSin; //Arrays for Ewald summation
} Buffer_S;

typedef struct
{
	vector<double>  t;         //time vector
	vector<double>  magZ;      //magnetization vector
}OutputVar_S;

typedef struct
{
	double lBox;        // length of simulation box
	double alpha;       // Ewald parameter
	double rCutLR;      // cut off for real space sum
	double nCutLR;      // cut off for reziprokal space
	double nc2, gu, gr, gs, w, cFluxDens, cForce, cSurf, cSelf; //constants for fourier space of Ewaldsum
	double alpha2, alpha4;
	int nc;
}Params_S;

typedef struct
{
	ArrayXXd intVecsCoordsLR;      // unit interaction vectors for interaction pairs
	ArrayXd  intVecsLenLR;         // length of interaction vectors for interaction pairs
	ArrayXi  intPairsLR;            // indices of interaction pairs
	int pairsLR;
}IntListLR_S;

typedef struct
{
	Array<int, 14, 3> offset;      // offset array
	ArrayXi nebrList;              // neighbour list
	ArrayXi cellList;              // cell-list
	Array<int, numbPart, 1> cx;    // cell indizes x
	Array<int, numbPart, 1> cy;    // cell indizes y
	Array<int, numbPart, 1> cz;    // cell indizes z
	Array<int, numbPart, 1> c;     // cell indizes 
	ArrayXXd intVecsCoordsSR;        // interaction vectors for interaction pairs
	ArrayXd intVecsLenSR;            // length of interaction vectors for interaction pairs
	ArrayXb intRangeSR;              // true for r<rCut, false for r>rCut

	double rNebr;                  // maximum radius of possible neighbours
	double dNebrShell;             // thickness of neighbour shell
	double lCell;                  // length of cell
	int nCells;                    // number of cells per dimension
	int numbCells;                 // whole number of cells
	double dispSum;                // sum of maximum displacements
	int pairsSR;                     // pairs in nebrtab
	double rCutSR;                   // cut off radius
}NebrList_S;

// working variables structure
typedef struct
{
	PartProps_S 	partProp;   	
	CoordsRot_S     coordsRot;     	
	CoordsTrans_S   coordsTrans;    
	ArrayXXd 		extFluxDens;	
	ArrayXXd 		demagFluxDens;	
	ArrayXXd 		thermTorque;	 
	ArrayXXd 		thermField;		
	ArrayXXd 		thermForce;		
	ArrayXXd 		intForce;		
	IntListLR_S     intListsLR;
	NebrList_S      nebrListData;

} WorkingVar_S;

#endif  /* TYPEDEFS_H_ */
