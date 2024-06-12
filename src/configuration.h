#ifndef CONFIGURATION_H
#define CONFIGURATION_H 
#include "INIReader.h"

#define SIZE_DIST_EQUAL 0   // equal size    
#define SIZE_DIST_LOG 1     // lognormal size distribution

#define HEUN_EULER 0 
#define BOGACKI_SHAMPINE 1
#define DORMAND_PRINCE 2 

class Configuration {
private:
    INIReader reader;

    int initConfigRot_ = 0; 
    int initConfigTrans_ = 0;
    bool readFromFile_ = false;

    int sizeDist_ = 0;

    bool enableMobilization_ = false;

    bool enableThermField_ = false;
    bool enableThermTorque_ = false;
    bool enableThermForce_ = false;
    int seed_ = 0;
    

    bool enableCoatingPot_ = false;

    int solver_ = 0;

    int numbPart_ = 0;
    int npDim_ = 0;
    double magDamp_ = 0;
    double temp_ = 0;
    double satMag_ = 0;
    double anisEn_ = 0;
    double vis_ = 0;
    double rhoMag_ = 0; 
    double rhoShell_ = 500; 

    double rMagMean_ = 0;
    double dShell_ = 0;
    double volFrac_ = 0;

    double densSurfMol_ = 0;
    double effLenSurfMol_ = 0; 

    double hamaker_ = 0;

    double z_valency_ = 0;
    double epsilon_r_ = 0;
    double gamma_0_ = 0;
    double conc_ = 0;

    double tMag_ = 0;
    double tRelax_ = 0;
    double deltaTinit_ = 0;

    double magFluxDens_ = 0;

    double absTol_ = 0;
    double deltaTmin_ = 0;
    double deltaTmax_ = 0;

    double errTolEwald_ = 0;

    double nebrShellFac_ = 0;
    double errTolSR_ = 0;

    double dataPoints_ = 0;
    double npCoords_ = 0;
    int lAcf_ = 0;
    double tEquil_ = 0;

    double my_pi_ = 0;
    double mu0_ = 0; 
    double kB_ = 0;
    double gyroMr_ = 0;
    double avogadro_ = 0;
    double epsilon_0_ = 0;
    double el_ = 0;

    double facEwald_ = 0;
    double velMmConst_ = 0;
    double anisConst_ = 0;
    double sterRepFac_ = 0;
    double gamma_ = 0;
    double epsilon_ = 0;
    double kappa_ = 0;
    double electrostatRepFac_ = 0;
    double irPi_ = 0;

    double sigma_ = 0; 
    double tEnd_ = 0; 

    // Singleton instance
    static Configuration* instance;

    Configuration();

public:
    static Configuration* getInstance();

    int getInitConfigRot() const;
    int getInitConfigTrans() const;
    bool getReadFromFile() const;

    int getSizeDist() const;

    bool getEnableMobilization() const;

    bool getEnableThermField() const;
    bool getEnableThermTorque() const;
    bool getEnableThermForce() const;
    int getSeed() const;

    bool getEnableCoatingPot() const;

    int getSolver() const;

    int getNpDim() const;
    int getNumbPart() const;
    double getMagDamp() const;
    double getTemp() const;
    double getSatMag() const;
    double getAnisEn() const;
    double getVis() const;
    double getRhoMag() const;
    double getRhoShell() const;

    double getRMagMean() const;
    double getDShell() const;
    double getVolFrac() const;

    double getDensSurfMol() const;
    double getEffLenSurfMol() const; 

    double getHamaker() const;

    double getZValency() const;
    double getEpsilonR() const;
    double getGamma0() const;
    double getConc() const;

    double getTMag() const;
    double getTRelax() const;
    double getDeltaTInit() const;

    double getMagFluxDens() const;

    double getAbsTol() const;
    double getDeltaTMin() const;
    double getDeltaTMax() const;

    double getErrTolEwald() const;

    double getNebrShellFac() const;
    double getErrTolSR() const;

    double getDataPoints() const;
    double getNPCoords() const;
    int getLAcf() const;
    double getTEquil() const;

    double getMyPi() const;
    double getMu0() const; 
    double getKB() const;
    double getGyroMr() const;
    double getAvogadro() const;
    double getEpsilon0() const;
    double getEl() const;

    double getFacEwald() const;
    double getVelMmConst() const;
    double getAnisConst() const;
    double getSterRepFac() const;
    double getGamma() const;
    double getEpsilon() const;
    double getKappa() const;
    double getElectrostatRepFac() const;
    double getIrPi() const;

    double getSigma() const;
    double getTEnd() const;

    std::string toString() const;

    ~Configuration();
};

extern Configuration* config; 

#endif // CONFIG_LOADER_H
