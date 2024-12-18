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

		sizeDist_ = reader.GetInteger("size_distribution", "SIZE_DIST", -1);
		enableMobilization_ = reader.GetBoolean("mobilization", "ENABLE_MOBILIZATION", false);

		enableThermTorque_ = reader.GetBoolean("thermal_fluctuations", "ENABLE_THERM_TORQUE", false);
		seed_ = reader.GetInteger("thermal_fluctuations", "seed", -1);

		solver_ = reader.GetInteger("adaptive_timestepping_solver", "SOLVER", -1);

		numbPart_ = reader.GetInteger("particle_properties", "numbPart", -1);
		magDamp_ = reader.GetReal("particle_properties", "magDamp", -1);
		temp_ = reader.GetReal("particle_properties", "temp", -1);
		satMag_ = reader.GetReal("particle_properties", "satMag", -1);
		anisEn_ = reader.GetReal("particle_properties", "anisEn", -1);

		rMagMean_ = reader.GetReal("assemble_properties", "rMagMean", -1);
		rHydrMean_ = reader.GetReal("assemble_properties", "rHydrMean", -1);
		shearRate_ = reader.GetReal("assemble_properties", "shearRate", -1);

		tMag_ = reader.GetReal("time_settings", "tMag", -1);
		tRelax_ = reader.GetReal("time_settings", "tRelax", -1);
		deltaTinit_ = reader.GetReal("time_settings", "deltaTinit", -1);

		magFluxDens_ = reader.GetReal("external_magnetic_field_t", "magFluxDens", -1);

		absTol_ = reader.GetReal("solver_settings", "absTol", -1);
		deltaTmin_ = reader.GetReal("solver_settings", "deltaTmin", -1);
		deltaTmax_ = reader.GetReal("solver_settings", "deltaTmax", -1);
		errTolMinMax_ = reader.GetReal("solver_settings", "errTolMinMax", -1);

		dataPoints_ = reader.GetReal("output_settings", "dataPoints", -1);

		my_pi_ = reader.GetReal("physical_constants", "my_pi", -1);
		mu0_ = 4.0 * my_pi_ * pow(10, -7);
		kB_ = reader.GetReal("physical_constants", "kB", -1);
		gyroMr_ = reader.GetReal("physical_constants", "gyroMr", -1);
		anisConst_ = 2.0 * anisEn_ / satMag_;

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

int Configuration::getSizeDist() const {
	return sizeDist_;
}

bool Configuration::getEnableMobilization() const {
	return enableMobilization_;
}

bool Configuration::getEnableThermTorque() const {
	return enableThermTorque_;
}

int Configuration::getSeed() const {
	return seed_;
}

int Configuration::getSolver() const {
	return solver_;
}

// Number of particles
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

double Configuration::getRHydrMean() const {
	return rHydrMean_;
}

double Configuration::getShearRate() const {
	return shearRate_;
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

double Configuration::getErrTolMinMax() const {
	return errTolMinMax_;
}

double Configuration::getDataPoints() const {
	return dataPoints_;
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

double Configuration::getAnisConst() const {
	return anisConst_;
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
	ss << "magDamp_: " << magDamp_ << std::endl;
	ss << "temp_: " << temp_ << std::endl;
	ss << "satMag_: " << satMag_ << std::endl;
	ss << "anisEn_: " << anisEn_ << std::endl;
	ss << "rMagMean_: " << rMagMean_ << std::endl;
	ss << "rHydrMean_: " << rHydrMean_ << std::endl;
	ss << "shearRate_: " << shearRate_ << std::endl;
	ss << "tMag_: " << tMag_ << std::endl;
	ss << "tRelax_: " << tRelax_ << std::endl;
	ss << "deltaTinit_: " << deltaTinit_ << std::endl;
	ss << "magFluxDens_: " << magFluxDens_ << std::endl;
	ss << "absTol_: " << absTol_ << std::endl;
	ss << "deltaTmin_: " << deltaTmin_ << std::endl;
	ss << "deltaTmax_: " << deltaTmax_ << std::endl;
	ss << "errTolMinMax_: " << errTolMinMax_ << std::endl;
	ss << "dataPoints_: " << dataPoints_ << std::endl;
	ss << "my_pi_: " << my_pi_ << std::endl;
	ss << "mu0_: " << mu0_ << std::endl;
	ss << "kB_: " << kB_ << std::endl;
	ss << "gyroMr_: " << gyroMr_ << std::endl;
	ss << "anisConst_: " << anisConst_ << std::endl;
	ss << "vis_: " << vis_ << std::endl;
	ss << "sigma_: " << sigma_ << std::endl;
	ss << "tEnd_: " << tEnd_ << std::endl;
	return ss.str();
}

Configuration::~Configuration() {
	delete instance;
}

Configuration* Configuration::instance = nullptr;
