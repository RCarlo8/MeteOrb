#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <map>
#include <vector>
#include <filesystem>

#include <stdexcept>  // Libreria per gestione eccezioni

#include <sys/stat.h>
#include <armadillo>

#include "broadcastEph.h"
#include "rinexNav.h"
#include "rinex3Nav.h"
#include "utility.h"

using namespace std;
using namespace arma;

// CONSTRUCTORS
BroadcastEph::BroadcastEph(){};
BroadcastEph::~BroadcastEph(){};

// METHODS

void BroadcastEph::setRinexType(string type){
	rinexTypeC = type;
};

void BroadcastEph::checkRinexVersionAndType(int &rinex_version, int &rinex_type, std::ifstream &fin)
{
	const std::string sTokenVER = "RINEX VERSION / TYPE";
	std::string line;
	int nLines = 0;
	while (!fin.eof())
	{
		nLines++;
		line.clear();
		getline(fin, line, '\n');
		size_t found_VER = line.find(sTokenVER);

		if (found_VER != std::string::npos)
		{
			std::istringstream iss(line);
			std::vector<std::string> words{std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
			// Rinex Version
			if ((std::stod(words[0])) < 4 && (std::stod(words[0])) >= 3)
			{
				rinex_version = 3;
			}
			else if ((std::stod(words[0])) < 3 && (std::stod(words[0]))>= 2)
			{
				rinex_version = 2;
			}
			else
			{
				//rinex_version = NULL;
				throw runtime_error("ERROR: THE RINEX VERSION IS NOT SUPPORTED");
			}
			// Rinex Type
			std::string type = words[3].substr(0, 1);
			if (type == "G")
			{
				rinex_type = 0;
				
			}
			else if (type == "E")
			{
				rinex_type = 1;
			}
			else if (type == "M")
			{
				rinex_type = 2;
			}
			else
			{
				throw runtime_error("ERROR: NOT POSSIBLE TO DEFINE TYPE OF RINEX");
			}
			// Set rinex type
			setRinexType(type);
			
			words.clear();
			break;
		}
		// Fail safe
		if (nLines > 1)
		{
			// Rinex version usually provided in first line of file
			std::cout << "RINEX VERSION TYPE NOT FOUND IN FIRST LINE OF FILE! \n";
			break;
		}
	}
}

void BroadcastEph::ComputePosGPS(Rinex3Nav NAV, ofstream &fout_nav)
{
	// Navigation data storage variables
	vector<Rinex3Nav::DataGPS> vecNavGPS;
	//Rinex3Nav::DataGPS epochNavGPS;

	map<int, vector<Rinex3Nav::DataGPS>>::iterator itNAV;
	vector<Rinex3Nav::DataGPS>::iterator itObs;
	string sat_name;
	// Iterate through PRN's to find corresponding ephemeris data
	for (itNAV = NAV._navGPS.begin(); itNAV != NAV._navGPS.end(); itNAV++)
	{
		// Current Satellite PRN
		int prn = itNAV->first;
		// Save vector of navigatinal GPS data of a setellite
		vecNavGPS = itNAV->second;
		// Compute starting time of the observatio
		stringstream ss;
		ss << 'G' << setw(2) << setfill('0') << prn;
		sat_name = ss.str();
		
		// iterazione sulla giornata.
		int end = start + SECPERDAY;

		for (int t = start; t <= end; t+=step)
		{
			int index = NAV.closerMinorEpoch(vecNavGPS,t);
			if(index != -1)
			{
				if ((t - vecNavGPS[index].gpsTime) < GPS_step)
				{
					Rinex3Nav::DataGPS obs = vecNavGPS[index];

					double tk = t - obs.TOE;

				if (tk > 302400)
				{

					tk -= 604800;
				}
				else if (tk < -302400)
				{

					tk += 604800;
				}
				// dts
				double dts = obs.clockBias + tk * obs.clockDrift + tk * tk * obs.clockDriftRate;
				// tk += dts;
				// cout << "dts : " << dts << endl;
				//
				//  Semi-major axis
				double A = pow(obs.Sqrt_a, 2);
				// Computed mean motion
				double n0 = sqrt(mu_GPS / (pow(A, 3)));
				// Corrected mean motion
				double n = n0 + obs.Delta_n;
				// Compute the mean anomaly for tk
				double Mk = obs.Mo + n * tk;
				// Solve iteratively kepler equation for eccentricity anomaly
				double e = obs.Eccentricity;
				double Ek = Mk;
				int counter = 0;
				while (counter < 1000)
				{
					double Ek_prev = Ek;
					Ek = Ek_prev + (Mk - Ek_prev + e * sin(Ek_prev)) / (1 - e * cos(Ek_prev));
					counter++;
				}

				// Compute true anomaly
				double vk = atan2(sqrt(1 - e * e) * sin(Ek), cos(Ek) - e);
				// double vk = 2*atan(sqrt((1+e)/(1-e))*tan(Ek/2));
				//  Compute argument of latitude
				double fi_k = vk + obs.Omega;
				// Corrected argument of latitude
				double uk = fi_k + obs.Cuc * cos(2 * fi_k) + obs.Cus * sin(2 * fi_k);
				// Compute the radial distance rk, considering corrections crc and crs
				double rk = A * (1 - e * cos(Ek)) + obs.Crc * cos(2 * fi_k) + obs.Crs * sin(2 * fi_k);
				// Compute the inclination ik of the orbital plane
				double ik = obs.Io + obs.Cic * cos(2 * fi_k) + obs.CIS * sin(2 * fi_k) + obs.IDOT * tk;

				// POSITIONS IN OBITAL
				double xk_o = rk * cos(uk);
				double yk_o = rk * sin(uk);

				// Compute the longitude of the ascending node λk (with respect to Greenwich)
				double omega_k = obs.OMEGA + (obs.Omega_dot - omega_dot_E) * tk - omega_dot_E * obs.TOE;

				// EARTH FIXED POSITIONS
				double xk = xk_o * cos(omega_k) - yk_o * cos(ik) * sin(omega_k);
				double yk = xk_o * sin(omega_k) + yk_o * cos(ik) * cos(omega_k);
				double zk = yk_o * sin(ik);

				// SATELLITES VELOCITY COMPUTATION
				double Ek_dot = n / (1 - e * cos(Ek));
				double vk_dot = Ek_dot * sqrt(1 - e * e) / (1 - e * cos(Ek));
				double ik_dot = obs.IDOT + 2 * vk_dot * (obs.CIS * cos(2 * fi_k) - obs.Cic * sin(2 * fi_k));
				double uk_dot = vk_dot + 2 * vk_dot * (obs.Cus * cos(2 * fi_k) - obs.Cuc * sin(2 * fi_k));
				double rk_dot = e * A * Ek_dot * sin(Ek) + 2 * vk_dot * (obs.Crs * cos(2 * fi_k) - obs.Crc * sin(2 * fi_k));
				double omega_k_dot = obs.Omega_dot - omega_dot_E;

				// IN PLANE VELOCITIES
				double xk_o_dot = rk_dot * cos(uk) - rk * uk_dot * sin(uk);
				double yk_o_dot = rk_dot * sin(uk) + rk * uk_dot * cos(uk);

				// EARTH FIXED VELOCITIES
				double xk_dot = -xk_o * omega_k_dot * sin(omega_k) + xk_o_dot * cos(omega_k) - yk_o_dot * sin(omega_k) * cos(ik) - yk_o * (omega_k_dot * cos(omega_k) * cos(ik) - ik_dot * sin(omega_k) * sin(ik));
				double yk_dot = xk_o * omega_k_dot * cos(omega_k) + xk_o_dot * sin(omega_k) + yk_o_dot * cos(omega_k) * cos(ik) - yk_o * (omega_k_dot * sin(omega_k) * cos(ik) + ik_dot * cos(omega_k) * sin(ik));
				double zk_dot = yk_o_dot * sin(ik) + yk_o * ik_dot * cos(ik);

				writeCsvLine(fout_nav,sat_name,gpstime2Date(t),xk,yk,zk);

				}
			}
	
		}
		
		/*  
		
			THE FOLLOWING COMMENTED CODE PROCESS ALL THE OBS IN NAV FILE
			AND FOR EACH OBS COMPUTE POSITION IN THE NEXT TWO HOUR EVERY
			30 SECONDS.




		for (itObs = vecNavGPS.begin(); itObs != vecNavGPS.end(); itObs++)
		{
			
			double t_start = itObs->gpsTime - itObs->GPS_week * 604800 - 18;
			double t_step = 30;
			double t_end = t_start + 2 * 3600;

			for (int i = t_start; i <= t_end; i += t_step)
			{
				double tk = i - itObs->TOE;

				if (tk > 302400)
				{

					tk -= 604800;
				}
				else if (tk < -302400)
				{

					tk += 604800;
				}
				// dts
				double dts = itObs->clockBias + tk * itObs->clockDrift + tk * tk * itObs->clockDriftRate;
				// tk += dts;
				// cout << "dts : " << dts << endl;
				//
				//  Semi-major axis
				double A = pow(itObs->Sqrt_a, 2);
				// Computed mean motion
				double n0 = sqrt(mu_GPS / (pow(A, 3)));
				// Corrected mean motion
				double n = n0 + itObs->Delta_n;
				// Compute the mean anomaly for tk
				double Mk = itObs->Mo + n * tk;
				// Solve iteratively kepler equation for eccentricity anomaly
				double e = itObs->Eccentricity;
				double Ek = Mk;
				int counter = 0;
				while (counter < 1000)
				{
					double Ek_prev = Ek;
					Ek = Ek_prev + (Mk - Ek_prev + e * sin(Ek_prev)) / (1 - e * cos(Ek_prev));
					counter++;
				}

				// Compute true anomaly
				double vk = atan2(sqrt(1 - e * e) * sin(Ek), cos(Ek) - e);
				// double vk = 2*atan(sqrt((1+e)/(1-e))*tan(Ek/2));
				//  Compute argument of latitude
				double fi_k = vk + itObs->Omega;
				// Corrected argument of latitude
				double uk = fi_k + itObs->Cuc * cos(2 * fi_k) + itObs->Cus * sin(2 * fi_k);
				// Compute the radial distance rk, considering corrections crc and crs
				double rk = A * (1 - e * cos(Ek)) + itObs->Crc * cos(2 * fi_k) + itObs->Crs * sin(2 * fi_k);
				// Compute the inclination ik of the orbital plane
				double ik = itObs->Io + itObs->Cic * cos(2 * fi_k) + itObs->CIS * sin(2 * fi_k) + itObs->IDOT * tk;

				// POSITIONS IN OBITAL
				double xk_o = rk * cos(uk);
				double yk_o = rk * sin(uk);

				// Compute the longitude of the ascending node λk (with respect to Greenwich)
				double omega_k = itObs->OMEGA + (itObs->Omega_dot - omega_dot_E) * tk - omega_dot_E * itObs->TOE;

				// EARTH FIXED POSITIONS
				double xk = xk_o * cos(omega_k) - yk_o * cos(ik) * sin(omega_k);
				double yk = xk_o * sin(omega_k) + yk_o * cos(ik) * cos(omega_k);
				double zk = yk_o * sin(ik);

				// SATELLITES VELOCITY COMPUTATION
				double Ek_dot = n / (1 - e * cos(Ek));
				double vk_dot = Ek_dot * sqrt(1 - e * e) / (1 - e * cos(Ek));
				double ik_dot = epochNavGPS.IDOT + 2 * vk_dot * (epochNavGPS.CIS * cos(2 * fi_k) - epochNavGPS.Cic * sin(2 * fi_k));
				double uk_dot = vk_dot + 2 * vk_dot * (epochNavGPS.Cus * cos(2 * fi_k) - epochNavGPS.Cuc * sin(2 * fi_k));
				double rk_dot = e * A * Ek_dot * sin(Ek) + 2 * vk_dot * (epochNavGPS.Crs * cos(2 * fi_k) - epochNavGPS.Crc * sin(2 * fi_k));
				double omega_k_dot = epochNavGPS.Omega_dot - omega_dot_E;

				// IN PLANE VELOCITIES
				double xk_o_dot = rk_dot * cos(uk) - rk * uk_dot * sin(uk);
				double yk_o_dot = rk_dot * sin(uk) + rk * uk_dot * cos(uk);

				// EARTH FIXED VELOCITIES
				double xk_dot = -xk_o * omega_k_dot * sin(omega_k) + xk_o_dot * cos(omega_k) - yk_o_dot * sin(omega_k) * cos(ik) - yk_o * (omega_k_dot * cos(omega_k) * cos(ik) - ik_dot * sin(omega_k) * sin(ik));
				double yk_dot = xk_o * omega_k_dot * cos(omega_k) + xk_o_dot * sin(omega_k) + yk_o_dot * cos(omega_k) * cos(ik) - yk_o * (omega_k_dot * sin(omega_k) * cos(ik) + ik_dot * cos(omega_k) * sin(ik));
				double zk_dot = yk_o_dot * sin(ik) + yk_o * ik_dot * cos(ik);

				writeCsvLine(fout_nav,sat_name,gpstime2Date(itObs->GPS_week * 604800 + 18 + i),xk,yk,zk);
				//writeCsvLine(fout_nav,sat_name,gpstime2Date(itObs->GPS_week * 604800 + 18 + i),xk,yk,zk,xk_dot,yk_dot,zk_dot,dts);
			}
		}*/
	} 
}

// Compute positions of GAL satellites from NAV parameters
void BroadcastEph::ComputePosGAL(Rinex3Nav NAV, ofstream &fout_nav)
{
	// Navigation data storage variables
	vector<Rinex3Nav::DataGAL> vecNavGAL;
	//Rinex3Nav::DataGAL epochNavGAL;
	// ITERATORS
	map<int, vector<Rinex3Nav::DataGAL>>::iterator itNAV;
	
	//vector<Rinex3Nav::DataGAL>::iterator itObs;
	string sat_name;


	// Iterate through PRN's to find corresponding ephemeris data
	for (itNAV = NAV._navGAL.begin(); itNAV != NAV._navGAL.end(); itNAV++)
	{
		// Current Satellite PRN
		int prn = itNAV->first;
		// Save vector of navigatinal GPS data of a setellite
		vecNavGAL = itNAV->second;
		
		// Compute starting time of the observation
		stringstream ss;
		ss << 'E' << setw(2) << setfill('0') << prn;
		sat_name = ss.str();

				// iterazione sulla giornata.
		int end = start + SECPERDAY;

		for (int t = start; t <= end; t+=step)
		{
			int index = NAV.closerMinorEpoch(vecNavGAL,t);

			if(index != -1)
			{
				if ((t - vecNavGAL[index].gpsTime) < GAL_step)
				{
					Rinex3Nav::DataGAL obs = vecNavGAL[index];

					double tk = t - obs.TOE;

				if (tk > 302400)
				{
					tk -= 604800;
				}
				else if (tk < -302400)
				{

					tk += 604800;
				}
				// dts
				double dts = obs.clockBias + tk * obs.clockDrift + tk * tk * obs.clockDriftRate;
				// tk += dts;
				// cout << "dts : " << dts << endl;
				//
				//  Semi-major axis
				double A = pow(obs.Sqrt_a, 2);
				// Computed mean motion
				double n0 = sqrt(mu_GPS / (pow(A, 3)));
				// Corrected mean motion
				double n = n0 + obs.Delta_n;
				// Compute the mean anomaly for tk
				double Mk = obs.Mo + n * tk;
				// Solve iteratively kepler equation for eccentricity anomaly
				double e = obs.Eccentricity;
				double Ek = Mk;
				int counter = 0;
				while (counter < 1000)
				{
					double Ek_prev = Ek;
					Ek = Ek_prev + (Mk - Ek_prev + e * sin(Ek_prev)) / (1 - e * cos(Ek_prev));
					counter++;
				}

				// Compute true anomaly
				double vk = atan2(sqrt(1 - e * e) * sin(Ek), cos(Ek) - e);
				// double vk = 2*atan(sqrt((1+e)/(1-e))*tan(Ek/2));
				//  Compute argument of latitude
				double fi_k = vk + obs.Omega;
				// Corrected argument of latitude
				double uk = fi_k + obs.Cuc * cos(2 * fi_k) + obs.Cus * sin(2 * fi_k);
				// Compute the radial distance rk, considering corrections crc and crs
				double rk = A * (1 - e * cos(Ek)) + obs.Crc * cos(2 * fi_k) + obs.Crs * sin(2 * fi_k);
				// Compute the inclination ik of the orbital plane
				double ik = obs.Io + obs.Cic * cos(2 * fi_k) + obs.CIS * sin(2 * fi_k) + obs.IDOT * tk;

				// POSITIONS IN OBITAL
				double xk_o = rk * cos(uk);
				double yk_o = rk * sin(uk);

				// Compute the longitude of the ascending node λk (with respect to Greenwich)
				double omega_k = obs.OMEGA + (obs.Omega_dot - omega_dot_E) * tk - omega_dot_E * obs.TOE;

				// EARTH FIXED POSITIONS
				double xk = xk_o * cos(omega_k) - yk_o * cos(ik) * sin(omega_k);
				double yk = xk_o * sin(omega_k) + yk_o * cos(ik) * cos(omega_k);
				double zk = yk_o * sin(ik);

				// SATELLITES VELOCITY COMPUTATION
				double Ek_dot = n / (1 - e * cos(Ek));
				double vk_dot = Ek_dot * sqrt(1 - e * e) / (1 - e * cos(Ek));
				double ik_dot = obs.IDOT + 2 * vk_dot * (obs.CIS * cos(2 * fi_k) - obs.Cic * sin(2 * fi_k));
				double uk_dot = vk_dot + 2 * vk_dot * (obs.Cus * cos(2 * fi_k) - obs.Cuc * sin(2 * fi_k));
				double rk_dot = e * A * Ek_dot * sin(Ek) + 2 * vk_dot * (obs.Crs * cos(2 * fi_k) - obs.Crc * sin(2 * fi_k));
				double omega_k_dot = obs.Omega_dot - omega_dot_E;

				// IN PLANE VELOCITIES
				double xk_o_dot = rk_dot * cos(uk) - rk * uk_dot * sin(uk);
				double yk_o_dot = rk_dot * sin(uk) + rk * uk_dot * cos(uk);

				// EARTH FIXED VELOCITIES
				double xk_dot = -xk_o * omega_k_dot * sin(omega_k) + xk_o_dot * cos(omega_k) - yk_o_dot * sin(omega_k) * cos(ik) - yk_o * (omega_k_dot * cos(omega_k) * cos(ik) - ik_dot * sin(omega_k) * sin(ik));
				double yk_dot = xk_o * omega_k_dot * cos(omega_k) + xk_o_dot * sin(omega_k) + yk_o_dot * cos(omega_k) * cos(ik) - yk_o * (omega_k_dot * sin(omega_k) * cos(ik) + ik_dot * cos(omega_k) * sin(ik));
				double zk_dot = yk_o_dot * sin(ik) + yk_o * ik_dot * cos(ik);

				writeCsvLine(fout_nav,sat_name,gpstime2Date(t),xk,yk,zk);

				}
			}
	
		}

		/*
		for (itObs = vecNavGAL.begin(); itObs != vecNavGAL.end(); itObs++)
		{
			// TODO converter local seconds tow
			double t_start = itObs->gpsTime - itObs->GAL_week * SECPERWEEK - numLeapSeconds;
			double t_step = 30;
			double t_end = t_start + (10 * 60);
			
			for (int i = t_start; i < t_end; i += t_step)
			{
				double tk = i - itObs->TOE;
				if (tk > (SECPERWEEK/2))
				{
					tk -= SECPERWEEK;
				}
				else if (tk < -(SECPERWEEK/2))
				{
					tk += SECPERWEEK;
				}
				// dts
				double dts = itObs->clockBias + tk * itObs->clockDrift + tk * tk * itObs->clockDriftRate;
				//
				// cout << "dts : " << dts << endl;
				//
				//  Semi-major axis
				double A = pow(itObs->Sqrt_a, 2);
				// Computed mean motion
				double n0 = sqrt(mu_GAL / (pow(A, 3)));
				// Corrected mean motion
				double n = n0 + itObs->Delta_n;
				// Compute the mean anomaly for tk
				double Mk = itObs->Mo + n * tk;
				// Solve iteratively kepler equation for eccentricity anomaly
				double e = itObs->Eccentricity;
				double Ek = Mk;
				int counter = 0;
				while (counter < 1000)
				{
					double Ek_prev = Ek;
					Ek = Ek_prev + (Mk - Ek_prev + e * sin(Ek_prev)) / (1 - e * cos(Ek_prev));
					counter++;
				}

				// Compute true anomaly
				double vk = atan2(sqrt(1 - e * e) * sin(Ek), cos(Ek) - e);
				// double vk = 2*atan(sqrt((1+e)/(1-e))*tan(Ek/2));
				//  Compute argument of latitude
				double fi_k = vk + itObs->Omega;
				// Corrected argument of latitude
				double uk = fi_k + itObs->Cuc * cos(2 * fi_k) + itObs->Cus * sin(2 * fi_k);
				// Compute the radial distance rk, considering corrections crc and crs
				double rk = A * (1 - e * cos(Ek)) + itObs->Crc * cos(2 * fi_k) + itObs->Crs * sin(2 * fi_k);
				// Compute the inclination ik of the orbital plane
				double ik = itObs->Io + itObs->Cic * cos(2 * fi_k) + itObs->CIS * sin(2 * fi_k) + itObs->IDOT * tk;

				// POSITIONS IN OBITAL
				double xk_o = rk * cos(uk);
				double yk_o = rk * sin(uk);

				// Compute the longitude of the ascending node λk (with respect to Greenwich)
				double omega_k = itObs->OMEGA + (itObs->Omega_dot - omega_dot_E) * tk - omega_dot_E * itObs->TOE;

				// EARTH FIXED POSITIONS
				double xk = xk_o * cos(omega_k) - yk_o * cos(ik) * sin(omega_k);
				double yk = xk_o * sin(omega_k) + yk_o * cos(ik) * cos(omega_k);
				double zk = yk_o * sin(ik);
				
				// SATELLITES VELOCITY COMPUTATION
				double Ek_dot = n/(1-e*cos(Ek));
				double vk_dot = Ek_dot*sqrt(1-e*e)/(1-e*cos(Ek));
				double ik_dot = itObs->IDOT + 2*vk_dot*(itObs->CIS*cos(2*fi_k) - itObs->Cic*sin(2*fi_k));
				double uk_dot = vk_dot + 2*vk_dot*(itObs->Cus*cos(2*fi_k) - itObs->Cuc*sin(2*fi_k));
				double rk_dot = e*A*Ek_dot*sin(Ek) + 2*vk_dot*(itObs->Crs*cos(2*fi_k) - itObs->Crc*sin(2*fi_k));
				double omega_k_dot = itObs->Omega_dot - omega_dot_E;

				// IN PLANE VELOCITIES
				double xk_o_dot = rk_dot*cos(uk) - rk*uk_dot*sin(uk);
				double yk_o_dot = rk_dot*sin(uk) + rk*uk_dot*cos(uk);

				// EARTH FIXED VELOCITIES
				double xk_dot = -xk_o*omega_k_dot*sin(omega_k) + xk_o_dot*cos(omega_k) - yk_o_dot*sin(omega_k)*cos(ik) - yk_o*(omega_k_dot*cos(omega_k)*cos(ik) - ik_dot*sin(omega_k)*sin(ik));
				double yk_dot = xk_o*omega_k_dot*cos(omega_k) + xk_o_dot*sin(omega_k) + yk_o_dot*cos(omega_k)*cos(ik) - yk_o*(omega_k_dot*sin(omega_k)*cos(ik) + ik_dot*cos(omega_k)*sin(ik));
				double zk_dot = yk_o_dot*sin(ik) + yk_o*ik_dot*cos(ik);

				// Write result to file
				writeCsvLine(fout_nav,sat_name,gpstime2Date(itObs->GAL_week * 604800 + 18 + i),xk,yk,zk);
				//writeCsvLine(fout_nav,sat_name,gpstime2Date(itObs->GAL_week * 604800 + 18 + i),xk,yk,zk,xk_dot,yk_dot,zk_dot,dts);
				
			}
		}*/
	}
}


void BroadcastEph::ReadRinex3()
{	
	std::cout << "|" << setw(59) << "|" << endl;
	std::cout << "| - ................ Start reading file ................ - |" << endl;
	// NAVIGATION FILE DATA OBJECT
	Rinex3Nav NAV;

	// Check block operations for storing and computing positions
	if (rinexType == 0) // Rinex 3 only GPS
	{
		NAV.readGPS(fin_nav);
		this->ComputePosGPS(NAV,fout_nav);
	}
	else if(rinexType == 1) // Rinex 3 only GALILEO
	{
		NAV.readGAL(fin_nav);
		this->ComputePosGAL(NAV,fout_nav);
	}
	else if(rinexType ==2) // Rinex 3 MIXED
	{	
		NAV.readMixed(this->fin_nav);
		this->ComputePosGAL(NAV,this->fout_nav);
		this->ComputePosGPS(NAV,this->fout_nav);
	}
	
	fout_nav.close(); 
	fin_nav.close();
	
	std::cout << "| - ............... Reached End of file ................ - |" << endl;
}

int main(int argc, char const *argv[])
{
	// Path to the directory
	string path = "data/navigational/";
	// This structure would distinguish a file from a directory
	struct stat sb;
	// list of files
	vector<string> files;
	// BroadcastEph object
	BroadcastEph bcast;
	// Looping until all the items of the directory are exhausted
	for (const auto &entry : fs::directory_iterator(path))
	{
		// Converting the path to const char * in the
		// subsequent lines
		std::__fs::filesystem::path outfilename = entry.path();
		std::string outfilename_str = outfilename.string();
		const char *path = outfilename_str.c_str();

		// Testing whether the path points to a
		// non-directory or not If it does, displays path
		if (stat(path, &sb) == 0 && !(sb.st_mode & S_IFDIR))
			files.push_back(outfilename_str);
	}

	std::sort(files.begin(), files.end());

	// Iteration on all files in the directory
	for (size_t i = 1; i < files.size(); i++)
	{	
		std::cout << "------------------------------------------------------------" << endl;
		string pathRinexNav = files[i];
		std::cout << "| - *******************  RINEX FILE  ******************* - |" << endl;
		std::cout << "| - " << pathRinexNav << " - |" << endl;		
		string output = "output/NAV/" + pathRinexNav.substr(18, 30) + ".csv";
		// Operazioni per inizializzare file lettura e scrittura
		fileSafeOut(output, bcast.fout_nav);
		fileSafeIn(pathRinexNav, bcast.fin_nav);
		
		string rec = pathRinexNav.substr(18,34);
		cout << rec << endl;
		// Obtain start time from file name
		//string receiverID = rec.substr(0,9);
		//cout << receiverID << endl;

		int year = stoi(pathRinexNav.substr(30,4));
		int doy = stoi(pathRinexNav.substr(34,3));
		double start = gpsTime(dayOfYearToDateTime(doy,year));
		//
		std::cout << "| - NAV FILE reference datetime : " 
				<< std::setw(22) << gpstime2Date(start) << " - |" << std::endl;		
		std::cout << "| - Start time in GPS seconds : " 
				<< std::setw(24) << start << " - |" << std::endl;
		std::cout << "| - File doy : " 
				<< std::setw(41) << doy << " - |" << std::endl;		
		//
		bcast.version = 0;
		bcast.rinexType = 0;
		bcast.start = start;

		try
		{
			bcast.checkRinexVersionAndType(bcast.version, bcast.rinexType, bcast.fin_nav);

			std::cout << "| - Rinex version : " << std::setw(36)
				 << bcast.version << " - |" << endl		// | 2 | 3 |
				 << "| - Type of rinex : " << std::setw(36)
				 << bcast.rinexTypeC << " - |" << endl; // | G GPS | E GALILEO | M MIXED |
			
			// Local variable for the checks
			int version = bcast.version;
			int rinexType = bcast.rinexType;

			// Version 2 decision block
			if (version == 2)
			{
				// Rinex2 only supports GPS files
				if (rinexType == 0)
				{
					// ReadRinex2(fin_obs, fin_nav, fout_log);
				}
				else
				{
					throw runtime_error("ERROR: Rinex 2 File Reader only supports GPS files.\n");
				}
			}
			// Version 3 decision block
			if (version == 3)
			{
				bcast.ReadRinex3();
			}

			std::cout << "------------------------------------------------------------" << endl;
		}
		catch(const std::exception& e)
		{
			std::cerr << e.what() << '\n';
		}
		// End of operations
		std::cout << "------------------------------------------------------------" << endl;
	}
}
