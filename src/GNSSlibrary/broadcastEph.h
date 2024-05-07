/*
    descrivi questa classe
*/
#ifndef BROADCASTEPH_H
#define BROADCASTEPH_H

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <map>
#include <vector>
#include <filesystem>
#include <sys/stat.h>
#include <armadillo>
#include "rinex3Nav.h"
#include "utility.h"


// Constants for position computations
const double mu_GPS = 398600500000000;	
const double mu_GAL = 398600441800000; 			    // m3/sec2
const double omega_dot_E = 7.2921151467 * 0.00001;  // rad/sec
const double c = 299792458;                         // m/sec
// Constant for epoch calculations
const int GPS_step = 2*3600;                        
const int GAL_step = 10*60;
const int step = 30;

// BROADCASTEPH class
class BroadcastEph 
{   
    public : 
        // RINEX data 
        int version;
        int rinexType;
        string rinexTypeC;
        double start;
        // Reading and writing streams
        ifstream fin_nav;
        ofstream fout_nav;

    public :
        // CONSTRUCTORS
        BroadcastEph();
        ~BroadcastEph();
        
        // METHODS
        /*
        * Function to set rinex type variable
        */
        void setRinexType(string type);
        /*
        * Function to check Rinex version and type
        */
        void checkRinexVersionAndType(int &rinex_version, int &rinex_type, std::ifstream &fin);
        /*
        * Function to compute position of all GPS satellites in the RINEX
        */
        void ComputePosGPS(Rinex3Nav NAV, std::ofstream &fout_nav);
        /*
        * Function to compute position of all GAL satellites in the RINEX
        */
        void ComputePosGAL(Rinex3Nav NAV, std::ofstream &fout_nav);
        /*
        * Function to read Rinex v3 Files
        */
        void ReadRinex3();
        /*
        * Function to read Rinex v2 Files
        */
        void ReadRinex2();
};
#endif // BROADCASTEPH_H
