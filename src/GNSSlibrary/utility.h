/*
*
*
*/

#include <iostream>
#include <iomanip>  
#include <string>
#include <fstream>
#include <sstream> 
#include <iterator>
#include <map>
#include <vector>
#include <iostream>
#include <ctime>

using namespace std;

#ifndef UTILITY_H_
#define UTILITY_H_

// constant to handle file system
namespace fs = std::__fs::filesystem;

const double   GPS_ZERO = 315964800.0;
const long     JAN61980 = 44244;
const long     JAN11901 = 15385;
const double   JD_GPS_epoch = 2444244.5;

const int      DAYINWEEK = 7;
const double   SECPERDAY = 86400.0;
const double   SECPERWEEK = 604800.0;

// Number of leap seconds as of September 2021
const int numLeapSeconds = 18;

// Functions
void fixTimeStamp(double* hr, double* min, double* sec);
double gpsTime(vector<double> epochInfo);
double computeGpsTime(int year, int month, int day, int hour, int minute, int second);
double convertToJulianDate(int year, int month, int day);
string gpstime2Date(double gps_t);

string printTime(double time);
vector<double> readDateTime(string line);
string printDate(vector<double> epoch);
string replaceChars(string str, char ch1, char ch2);
void eraseSubStr(string & mainStr, const std::string & toErase);
string HHMMSS(double hours, double mins, double secs);
vector<double> dayOfYearToDateTime(int year, int day_of_year);

void fileSafeIn(string filename, ifstream &fin);
void fileSafeOut(string output_filename, ofstream &fout);
void firstCSVLine(ofstream &fout,int flag);
void writeCsvLine(ofstream &fout,string name,string datetime,double xk,double yk,double zk);
void writeCsvLine(ofstream &fout,string name,string datetime,double xk,double yk,double zk,double xk_dot,double yk_dot,double zk_dot,double dts);
void writeCsvLine(ofstream &fout,string name,string datetime,double xk,double yk,double zk,double clock);
double getRcm();



#endif /* UTILITY_H_ */
