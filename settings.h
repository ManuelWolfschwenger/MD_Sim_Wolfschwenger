/*
* Name : setting.h
* Auto :
* Brief: Contains the general settings for the simulation
*/

#ifndef SETTINGS_H_
#define SETTINGS_H_

// initial configuration
#define INIT_CONFIG_ROT         2       // 1...z - orientated mag moment, random easy axis, 2...random easy axis and mag moment(equal direction), 3...easy axis and mag moment in z direction
#define INIT_CONFIG_TRANS       1       // 1...random coordinates, 2...organized in a grid
#define READ_FROM_FILE          false    // starts simulatation where last simulations ended attention: change name of coords file before!!!!!!!!!!!!!

//size distribution
#define SIZE_DIST_EQUAL         0
#define SIZE_DIST_LOG           1
#define SIZE_DIST               SIZE_DIST_EQUAL       // 1...equal size, 2...lognormal size dist<

// mobilization               
#define ENABLE_MOBILIZATION     true   // false... immobilized, true...mobilized

//thermal fluctuations
#define ENABLE_THERM_FIELD      true    // true...enable thermal fluctuations, else false
#define ENABLE_THERM_TORQUE     true    // true...enable thermal torque, else false
#define ENABLE_THERM_FORCE      true    // true...enable thermal torque, else false

// interactions
#define ENABLE_COATING_POT      false    // true...interaction via coating potentials, false...manually correct overlap/ hard sphere

// adaptive timestepping solver
#define HEUN_EULER              0 //Heun-Euler
#define BOG_SHAMP               1 //Bogacki-Shampine
#define DOR_PRI                 2 //Dormand-Prince
#define SOLVER                  DOR_PRI

#endif  /* SETTINGS_H_ */
