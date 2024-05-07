/* 
* 
* 
*  
*  
*  
*  
*/
#ifndef RINEX3NAV_H
#define RINEX3NAV_H

#include <fstream>
#include <sstream> 
#include <iostream>
#include <iomanip>  
#include <string>
#include <iterator>
#include <map>
#include <vector>
#include <armadillo>
#include "rinexNav.h" 
#include "utility.h"

using namespace std;
using namespace arma;

class Rinex3Nav 
{
	
	public:

		struct DataGPS {
			bool isAvailable;
			int PRN;
			vector<double> epochInfo;
			double gpsTime;
			double clockBias; 
			double clockDrift; 
			double clockDriftRate; 
			double IODE; 
			double Crs;
			double Delta_n; 
			double Mo;
			double Cuc;
			double Eccentricity;  
			double Cus;  
			double Sqrt_a; 
			double TOE; 
			double Cic;  
			double OMEGA;  
			double CIS;  
			double Io;  
			double Crc;  
			double Omega;  
			double Omega_dot;  
			double IDOT;  
			double L2_codes_channel;  
			double GPS_week;  
			double L2_P_data_flag;  
			double svAccuracy;  
			double svHealth;  
			double TGD;  
			double IODC;  
			double transmission_time;  
			double fit_interval;
		}; struct dataGPS;

		struct DataGAL {
			bool isAvailable;
			int PRN;
			vector<double> epochInfo;
			double gpsTime;
			double clockBias; 
			double clockDrift; 
			double clockDriftRate;
			double IOD;
			double Crs;
			double Delta_n;
			double Mo;
			double Cuc;
			double Eccentricity;
			double Cus;
			double Sqrt_a;
			double TOE;
			double Cic;
			double OMEGA;
			double CIS;
			double Io;
			double Crc;
			double Omega;
			double Omega_dot;
			double IDOT;
			double GAL_week;
			double SISA;
			double svHealth;
			double BGD_E5a;
			double BGD_E5b;
			double transmission_time;
		}; struct dataGAL;

		struct HeaderGPS {
			// Ionospheric alpha and beta constants
			vector<double> ialpha;
			vector<double> ibeta;
			// Time System correction
			vector<double> GPUT;
		}; 

		struct HeaderGAL {
			// Ionospheric alpha constants
			double ialpha;
			// Time System correction
			double GAUT;
			double GPGA;
			double leapSec;
		};

		// Parameters containers :
		map<int, vector<Rinex3Nav::DataGPS>> _navGPS;
		map<int, vector<Rinex3Nav::DataGAL>> _navGAL;

		// Header data containers :
		HeaderGPS _headerGPS;
		HeaderGAL _headerGAL;

		double startTime;
		double endTime;

		Rinex3Nav();
		~Rinex3Nav();
		void readGPS(std::ifstream& inputfileGPS); 
		void readGAL(std::ifstream& inputfileGAL); 
		void readMixed(std::ifstream& inputfileMixed); // for mixed navigation files
		int closerMinorEpoch(vector<Rinex3Nav::DataGPS> vettore, int valoreDiRiferimento);
		int closerMinorEpoch(vector<Rinex3Nav::DataGAL> vettore, int valoreDiRiferimento); 

	 
		
};

#endif // RINEX3NAV_H 
