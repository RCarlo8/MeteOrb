// Class definitions for CoreSky objects
// and prototypes for precise orbit interpolation methods

#ifndef CORESKY_H
#define CORESKY_H

#include <fstream>
#include <sstream> 
#include <iostream>
#include <iomanip>  
#include <string>
#include <iterator>
#include <map>
#include <vector>
#include <armadillo>


#include "utility.h"


using namespace std;
using namespace arma;

const vector<char> defSatelliteSystems = {'G','E'};

class CoreSky
{
    public:

        CoreSky();
        CoreSky(string pathFilename, vector<char> inputSS, ios_base::openmode );
        ~CoreSky();

        void setSS(vector<char> inputSS);
        void setPathFilenameMode(string pathFilename,
                                    ios_base::openmode );

        vec getEpoch();
        double getStartTime();
        double getEndTime();
        vector<string> getSatIds();
        void printSatIds();
        int getNumOfsat();
        mat getClock_corr();
        cube getCoord();

        double gpsTimeSp3(vector<double> epochInfo);

        void clearCoord();
        void cleaSetPoly();
        void clearVariables();

        void initHeader();
        int readHeader();
        void readCoord();

        void computePolyCoef();
        void computeClockPolyCoef();
        cube coordInterpolate(vec time, vec satId);
        void printInterpolatedPos(cube X, vec time, vec sat);
            
        vector<double> SP3EpochTimeOrganizer(string line);
        void savePos(cube Pos);

    private:

        string               pathFilename;
        fstream              fileStream;
        ios_base::openmode   fileMode;
        bool                 isempty = true;
        vector<char> satelliteSystems;
        
        // Data info variables

        char                 formatVersion;
        char                 modeFlag;
        
        double               SP3StartTime;
        double               SP3EndTime;
        unsigned long        numberSP3epochs;

        string               dataUsed;
        string               coordFrame;
        string               orbitType;
        string               sourceAgency;

        unsigned long        gpsWeek;
        double               secsOfWeek;
        double               SP3interval;
        long                 SP3mjd;
        double               SP3fmjd;

        // Sat coord variables
        
        vector<string>       satIds;
        vector<double>       epochs;
        vec                  ref_epoch;
        
        // Coordinates variables

        cube                 coord;
        mat                  clock;
        unsigned short       numberSP3svs;
        vector<unsigned short> svAccu;
        vector<unsigned short> sp3PRN;
        unsigned short       numberSVparams;
        unsigned short       numberGoodPRNs;
        unsigned short       numberGoodACCURs;

        // Interpolation variables
        vector<cube>         coef_poly;
        vector<mat>          coef_poly_clock;

};

#endif




      