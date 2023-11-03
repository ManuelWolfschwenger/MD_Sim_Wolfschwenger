/*
* Name : constant.h
* Auto :
* Brief: Contains the mathematical constants
*/

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

// external libaries
#include <cmath>

//number of particles
const int numbPart = 500;                     // number of particles

// particle properties constants
const double magDamp = 0.001;                      // magnetic damping cnstant
const double temp = 300.0;                        // temperature in K
const double satMag = 4.5e5;                      // saturation magnetization	
const double anisEn = 1e4;                        // anisotropy energy
const double vis = 853.8e-6;                     // viscosity of medium (water)

// assemble properties
const double rMagMean = 8e-9;                    // mean magnetic radius
const double dShell = 0;                       // shell thickness for hard sphere diameter (only used with manual hard sphere correction , when REPULSION is disabled)
const double volFrac = 0.1;                      // particle volume concentration

// steric repulsion
const double densSurfMol = 5e17;                  // surface density of surfactant molecules in molecules/m^2
const double effLenSurfMol = 0.2 * rMagMean;      // effective length of polymer molecules

// vdw attraction
const double hamaker = 33e-21;                    // hamaker constant in J

// electrostat repulsion
const double z_valency = 1;                       // valency of ions
const double epsilon_r = 78.5;                    // relative permittivity
const double gamma_0 = 0.04;                      // surface potential in V
const double conc = 0.05 * 1000;                   // concentration of ions from mol/l to mol/m^3
 
// time settings
const double tMag = 0;                             // magnetization time in s
const double tRelax = 3e-5;                        // relaxation time in s
const double deltaTinit = 1e-12;                   // initial guess of deltaT                 

// external magnetic field in T
const double magFluxDens = 0.003;                    // external magnetic flux density

// solver settings
const double absTol = 1e-4;                        // relative error tolerance for rotational motion
const double deltaTmin = 1e-14;                    // minimal deltaT for integration
const double deltaTmax = 1e-9;                    // maximal deltaT for integration

// treatment of long range interactions
const double errTolEwald = 1e-4;                               

// short range interactions cell list and neighbour list
const double nebrShellFac = 3;                     // nebrShellFac*rMagMean = dNebrShell (Neighbour-Shell thickness)
const double errTolSR = 0.01;                      // F_vdw/F_dip to cut off short range forces

// output settings
const double dataPoints = 2000.0;

// physical constants
const double my_pi = 3.14159265359;               // pi
const double mu0 = 4.0 * my_pi * pow(10, -7);     // magnetic field constant
const double kB = 1.380648e-23;                   // boltzmann constant
const double gyroMr = 1.7595e11;                  // gyromagnetic ratio
const double avogadro = 6.02214076e23;            // avogadro constant in mol^-1
const double epsilon_0 = 8.8541e-12;              // permittivity of vacuum in C^2N^-1m^-2
const double el = 1.6e-19;                        // charge of an electron in C

// constants for faster calculation
const double facEwald = 1e-7;                    //mu0/(4*pi)
const double velMmConst = gyroMr / (1.0 + pow(magDamp, 2.0));
const double anisConst = 2.0 * anisEn / satMag;
const double sterRepFac = 4 * kB * temp * my_pi * densSurfMol / (2 * effLenSurfMol);
const double gamma = tanh(z_valency * el * gamma_0 / (4 * kB * temp));
const double epsilon = epsilon_0 * epsilon_r;
const double kappa = sqrt(el * el * avogadro * conc * z_valency * z_valency / (epsilon * kB * temp));
const double electrostatRepFac = kappa * 64 * my_pi * epsilon * pow(kB * temp / (z_valency * el),2.0) * gamma * gamma;
const double irPi = 1 / sqrt(my_pi);

//calculations
const double sigma = sqrt(0.1) * 2.0 * rMagMean;   // standard deviation of size distribution
const double tEnd = tMag + tRelax;                 // end time in s

#endif  //CONSTANTS_H_ */

