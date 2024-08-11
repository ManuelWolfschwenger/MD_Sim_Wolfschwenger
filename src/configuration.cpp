#include <iostream>
#include <cmath>
#include <sstream>
#include "configuration.h"

Configuration::Configuration() : reader("../configuration.ini") {
	if (reader.ParseError() != 0) {
		std::cerr << "Error loading configuration file." << std::endl;
	}
	else {
		initConfigRot_ = reader.GetInteger("initial_configuration", "INIT_CONFIG_ROT", -1);
		initConfigTrans_ = reader.GetInteger("initial_configuration", "INIT_CONFIG_TRANS", -1);
		readFromFile_ = reader.GetBoolean("initial_configuration", "READ_FROM_FILE", false);

		sizeDist_ = reader.GetInteger("size_distribution", "SIZE_DIST", -1);

		enableMobilization_ = reader.GetBoolean("mobilization", "ENABLE_MOBILIZATION", false);

		enableThermField_ = reader.GetBoolean("thermal_fluctuations", "ENABLE_THERM_FIELD", false);
		enableThermTorque_ = reader.GetBoolean("thermal_fluctuations", "ENABLE_THERM_TORQUE", false);
		enableThermForce_ = reader.GetBoolean("thermal_fluctuations", "ENABLE_THERM_FORCE", false);
		seed_ = reader.GetInteger("thermal_fluctuations", "seed", -1);

		enableCoatingPot_ = reader.GetBoolean("interactions", "ENABLE_COATING_POT", false);

		solver_ = reader.GetInteger("adaptive_timestepping_solver", "SOLVER", -1);

		npDim_ = reader.GetInteger("particle_properties", "npDim", -1);
		magDamp_ = reader.GetReal("particle_properties", "magDamp", -1);
		temp_ = reader.GetReal("particle_properties", "temp", -1);
		satMag_ = reader.GetReal("particle_properties", "satMag", -1);
		anisEn_ = reader.GetReal("particle_properties", "anisEn", -1);

		rMagMean_ = reader.GetReal("assemble_properties", "rMagMean", -1);
		dShell_ = reader.GetReal("assemble_properties", "dShell", -1);
		volFrac_ = reader.GetReal("assemble_properties", "volFrac", -1);

		densSurfMol_ = reader.GetReal("steric_repulsion", "densSurfMol", -1);
		effLenSurfMol_ = 0.2 * rMagMean_;

		hamaker_ = reader.GetReal("vdw_attraction", "hamaker", -1);

		z_valency_ = reader.GetReal("electrostat_repulsion", "z_valency", -1);
		epsilon_r_ = reader.GetReal("electrostat_repulsion", "epsilon_r", -1);
		gamma_0_ = reader.GetReal("electrostat_repulsion", "gamma_0", -1);
		conc_ = reader.GetReal("electrostat_repulsion", "conc", -1);

		tMag_ = reader.GetReal("time_settings", "tMag", -1);
		tRelax_ = reader.GetReal("time_settings", "tRelax", -1);
		deltaTinit_ = reader.GetReal("time_settings", "deltaTinit", -1);

		magFluxDens_ = reader.GetReal("external_magnetic_field_t", "magFluxDens", -1);

		absTol_ = reader.GetReal("solver_settings", "absTol", -1);
		deltaTmin_ = reader.GetReal("solver_settings", "deltaTmin", -1);
		deltaTmax_ = reader.GetReal("solver_settings", "deltaTmax", -1);

		errTolEwald_ = reader.GetReal("long_range_interactions", "errTolEwald", -1);
		errTolSR_ = reader.GetReal("short_range_interactions", "errTolSR", -1);

		npCoords_ = reader.GetReal("output_settings", "npCoords", -1);

		my_pi_ = reader.GetReal("physical_constants", "my_pi", -1);
		mu0_ = 4.0 * my_pi_ * pow(10, -7);
		kB_ = reader.GetReal("physical_constants", "kB", -1);
		gyroMr_ = reader.GetReal("physical_constants", "gyroMr", -1);
		avogadro_ = reader.GetReal("physical_constants", "avogadro", -1);
		epsilon_0_ = reader.GetReal("physical_constants", "epsilon_0", -1);
		el_ = reader.GetReal("physical_constants", "el", -1);

		numbPart_ = npDim_ * npDim_ * npDim_;
		facEwald_ = reader.GetReal("faster_calculation", "facEwald", -1);
		velMmConst_ = gyroMr_ / (1.0 + pow(magDamp_, 2.0));
		anisConst_ = 2.0 * anisEn_ / satMag_;
		sterRepFac_ = 4 * kB_ * temp_ * my_pi_ * densSurfMol_ / (2 * effLenSurfMol_);
		gamma_ = tanh(z_valency_ * el_ * gamma_0_ / (4 * kB_ * temp_));
		epsilon_ = epsilon_0_ * epsilon_r_;
		kappa_ = sqrt(el_ * el_ * avogadro_ * conc_ * z_valency_ * z_valency_ / (epsilon_ * kB_ * temp_));
		electrostatRepFac_ = kappa_ * 64 * my_pi_ * epsilon_ * pow(kB_ * temp_ / (z_valency_ * el_), 2.0) * gamma_ * gamma_;
		irPi_ = 1 / sqrt(my_pi_);

		vis_ = (2.414e-5) * pow(10, 247.8 / (temp_ - 140));
		sigma_ = sqrt(0.1) * 2.0 * rMagMean_;
		tEnd_ = tMag_ + tRelax_;
	}
}

Configuration* Configuration::getInstance() {
	if (!instance) {
		instance = new Configuration();
	}
	return instance;
}

int Configuration::getInitConfigRot() const {
	return initConfigRot_;
}

int Configuration::getInitConfigTrans() const {
	return initConfigTrans_;
}

bool Configuration::getReadFromFile() const {
	return readFromFile_;
}

int Configuration::getSizeDist() const {
	return sizeDist_;
}

bool Configuration::getEnableMobilization() const {
	return enableMobilization_;
}

bool Configuration::getEnableThermField() const {
	return enableThermField_;
}

bool Configuration::getEnableThermTorque() const {
	return enableThermTorque_;
}

int Configuration::getSeed() const {
	return seed_;
}

bool Configuration::getEnableThermForce() const {
	return enableThermForce_;
}

bool Configuration::getEnableCoatingPot() const {
	return enableCoatingPot_;
}

int Configuration::getSolver() const {
	return solver_;
}

// Number of particles
int Configuration::getNpDim() const {
	return npDim_;
}

int Configuration::getNumbPart() const {
	return numbPart_;
}

// Magnetic damping constant
double Configuration::getMagDamp() const {
	return magDamp_;
}

// Temperature in K
double Configuration::getTemp() const {
	return temp_;
}

// Saturation magnetization
double Configuration::getSatMag() const {
	return satMag_;
}

// Anistropy energy
double Configuration::getAnisEn() const {
	return anisEn_;
}

// Viscosity of medium (water9
double Configuration::getVis() const {
	return vis_;
}

// Mean magnetic radius
double Configuration::getRMagMean() const {
	return rMagMean_;
}

// Shell thickness for hard sphere diameter (only used with manual hard sphere correction , when REPULSION is disabled)
double Configuration::getDShell() const {
	return dShell_;
}

// Particle volume concentration
double Configuration::getVolFrac() const {
	return volFrac_;
}

// Surface density of surfactant molecules in molecules/m^2
double Configuration::getDensSurfMol() const {
	return densSurfMol_;
}

// Effective length of polymer molecules
double Configuration::getEffLenSurfMol() const {
	return effLenSurfMol_;
}

// hamaker constant in J
double Configuration::getHamaker() const {
	return hamaker_;
}

// Valency of ions
double Configuration::getZValency() const {
	return z_valency_;
}
// Relative permittivity
double Configuration::getEpsilonR() const {
	return epsilon_r_;
}
// Surface potential in V
double Configuration::getGamma0() const {
	return gamma_0_;
}

// Concentration of ions from mol/l to mol/m^3
double Configuration::getConc() const {
	return conc_;
}

// Magnetization time in s
double Configuration::getTMag() const {
	return tMag_;
}
// Relaxation time in s
double Configuration::getTRelax() const {
	return tRelax_;
}

// Initial guess of deltaT 
double Configuration::getDeltaTInit() const {
	return deltaTinit_;
}

// External magnetic flux density
double Configuration::getMagFluxDens() const {
	return magFluxDens_;
}

// Relative error tolerance for rotational motion
double Configuration::getAbsTol() const {
	return absTol_;
}

// Minimal deltaT for integration
double Configuration::getDeltaTMin() const {
	return deltaTmin_;
}

// Maximal deltaT for integration
double Configuration::getDeltaTMax() const {
	return deltaTmax_;
}

// Treatment of long range interactions
double Configuration::getErrTolEwald() const {
	return errTolEwald_;
}

// F_vdw/F_dip to cut off short range forces
double Configuration::getErrTolSR() const {
	return errTolSR_;
}

double Configuration::getNPCoords() const {
	return npCoords_;
}

double Configuration::getMyPi() const {
	return my_pi_;
}

// Magnetic field constant
double Configuration::getMu0() const {
	return mu0_;
}

// Boltzmann constant
double Configuration::getKB() const {
	return kB_;
}

// Gyromagnetic ratio
double Configuration::getGyroMr() const {
	return gyroMr_;
}

// Avogadro constant in mol^-1
double Configuration::getAvogadro() const {
	return avogadro_;
}

// Permittivity of vacuum in C^2N^-1m^-2
double Configuration::getEpsilon0() const {
	return epsilon_0_;
}

// Charge of an electron in C
double Configuration::getEl() const {
	return el_;
}

//mu0/(4*pi)
double Configuration::getFacEwald() const {
	return facEwald_;
}

double Configuration::getVelMmConst() const {
	return velMmConst_;
}

double Configuration::getAnisConst() const {
	return anisConst_;
}

double Configuration::getSterRepFac() const {
	return sterRepFac_;
}

double Configuration::getGamma() const {
	return gamma_;
}

double Configuration::getEpsilon() const {
	return epsilon_;
}

double Configuration::getKappa() const {
	return kappa_;
}

double Configuration::getElectrostatRepFac() const {
	return electrostatRepFac_;
}

double Configuration::getIrPi() const {
	return irPi_;
}

// Standard deviation of size distribution
double Configuration::getSigma() const {
	return sigma_;
}

// End time in s
double Configuration::getTEnd() const {
	return tEnd_;
}


std::string Configuration::toString() const {
	std::stringstream ss;
	ss << "numbPart_: " << numbPart_ << std::endl;
	ss << "seed_: " << seed_ << std::endl;
	ss << "npDim_: " << npDim_ << std::endl;
	ss << "magDamp_: " << magDamp_ << std::endl;
	ss << "temp_: " << temp_ << std::endl;
	ss << "satMag_: " << satMag_ << std::endl;
	ss << "anisEn_: " << anisEn_ << std::endl;
	ss << "rMagMean_: " << rMagMean_ << std::endl;
	ss << "dShell_: " << dShell_ << std::endl;
	ss << "volFrac_: " << volFrac_ << std::endl;
	ss << "densSurfMol_: " << densSurfMol_ << std::endl;
	ss << "effLenSurfMol_: " << effLenSurfMol_ << std::endl;
	ss << "hamaker_: " << hamaker_ << std::endl;
	ss << "z_valency_: " << z_valency_ << std::endl;
	ss << "epsilon_r_: " << epsilon_r_ << std::endl;
	ss << "gamma_0_: " << gamma_0_ << std::endl;
	ss << "conc_: " << conc_ << std::endl;
	ss << "tMag_: " << tMag_ << std::endl;
	ss << "tRelax_: " << tRelax_ << std::endl;
	ss << "deltaTinit_: " << deltaTinit_ << std::endl;
	ss << "magFluxDens_: " << magFluxDens_ << std::endl;
	ss << "absTol_: " << absTol_ << std::endl;
	ss << "deltaTmin_: " << deltaTmin_ << std::endl;
	ss << "deltaTmax_: " << deltaTmax_ << std::endl;
	ss << "errTolEwald_: " << errTolEwald_ << std::endl;
	ss << "errTolSR_: " << errTolSR_ << std::endl;
	ss << "npCoords_: " << npCoords_ << std::endl;
	ss << "my_pi_: " << my_pi_ << std::endl;
	ss << "mu0_: " << mu0_ << std::endl;
	ss << "kB_: " << kB_ << std::endl;
	ss << "gyroMr_: " << gyroMr_ << std::endl;
	ss << "avogadro_: " << avogadro_ << std::endl;
	ss << "epsilon_0_: " << epsilon_0_ << std::endl;
	ss << "el_: " << el_ << std::endl;
	ss << "facEwald_: " << facEwald_ << std::endl;
	ss << "velMmConst_: " << velMmConst_ << std::endl;
	ss << "anisConst_: " << anisConst_ << std::endl;
	ss << "sterRepFac_: " << sterRepFac_ << std::endl;
	ss << "gamma_: " << gamma_ << std::endl;
	ss << "epsilon_: " << epsilon_ << std::endl;
	ss << "kappa_: " << kappa_ << std::endl;
	ss << "electrostatRepFac_: " << electrostatRepFac_ << std::endl;
	ss << "irPi_: " << irPi_ << std::endl;
	ss << "vis_: " << vis_ << std::endl;
	ss << "sigma_: " << sigma_ << std::endl;
	ss << "tEnd_: " << tEnd_ << std::endl;
	return ss.str();
}

Configuration::~Configuration() {
	delete instance;
}

Configuration* Configuration::instance = nullptr;
