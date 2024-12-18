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

    int sizeDist_ = 0;

    bool enableMobilization_ = false;

    bool enableThermField_ = false;
    bool enableThermTorque_ = false;
    int seed_ = 0;

    int solver_ = 0;

    int numbPart_ = 0;
    double magDamp_ = 0;
    double temp_ = 0;
    double satMag_ = 0;
    double anisEn_ = 0;
    double vis_ = 0;

    double rMagMean_ = 0;
    double rHydrMean_ = 0;
    double shearRate_ = 0;

    double tMag_ = 0;
    double tRelax_ = 0;
    double deltaTinit_ = 0;

    double magFluxDens_ = 0;

    double absTol_ = 0;
    double deltaTmin_ = 0;
    double deltaTmax_ = 0;
    double errTolMinMax_ = 0;

    double dataPoints_ = 0;

    double my_pi_ = 0;
    double mu0_ = 0; 
    double kB_ = 0;
    double gyroMr_ = 0;

    double velMmConst_ = 0;
    double anisConst_ = 0;

    double sigma_ = 0; 
    double tEnd_ = 0; 

    // Singleton instance
    static Configuration* instance;

    Configuration();

public:
    static Configuration* getInstance();

    int getInitConfigRot() const;

    int getSizeDist() const;

    bool getEnableMobilization() const;
    bool getEnableThermTorque() const;
    int getSeed() const;

    int getSolver() const;

    int getNumbPart() const;
    double getMagDamp() const;
    double getTemp() const;
    double getSatMag() const;
    double getAnisEn() const;
    double getVis() const;

    double getRMagMean() const;
    double getRHydrMean() const;
    double getShearRate() const;

    double getTMag() const;
    double getTRelax() const;
    double getDeltaTInit() const;

    double getMagFluxDens() const;

    double getAbsTol() const;
    double getDeltaTMin() const;
    double getDeltaTMax() const;
    double getErrTolMinMax() const;

    double getDataPoints() const;

    double getMyPi() const;
    double getMu0() const; 
    double getKB() const;
    double getGyroMr() const;

    double getAnisConst() const;

    double getSigma() const;
    double getTEnd() const;

    std::string toString() const;

    ~Configuration();
};

extern Configuration* config; 

#endif // CONFIG_LOADER_H
