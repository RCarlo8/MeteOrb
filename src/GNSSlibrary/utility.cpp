/*
* This class contains methods used to manage time and string
* input readed from files.
*/


#include "utility.h"
#include <ctime>



using namespace std;

/******************	FUNZIONI TIME ***********************************************************************************************/

string gpstime2Date(double gpsTime){

	time_t converted_time = gpsTime + GPS_ZERO - numLeapSeconds;
	struct tm* utc_tm = gmtime(&converted_time);
	
	string date = to_string(utc_tm->tm_year+1900)+"-"+to_string(utc_tm->tm_mon+1)+"-"+to_string(utc_tm->tm_mday);
	int hours = utc_tm->tm_hour;// (utc_timestamp / 3600) % 24; // 3600 seconds in an hour
	int minutes = utc_tm->tm_min;// (utc_timestamp / 60) % 60; // 60 seconds in a minute
	int seconds = utc_tm->tm_sec;// % 60; // Seconds component
	string time;
	// Display hours, minutes, and seconds with leading zeros
	//return to_string(utc_tm->tm_year+1900)+"-"+to_string(utc_tm->tm_mon+1)+"-"+to_string(utc_tm->tm_mday)+" "+to_string(utc_tm->tm_hour)+":"+to_string(utc_tm->tm_min)+":"+to_string(utc_tm->tm_sec);
	// Format Unix time as a date string
    if (hours < 10) time.append("0");
	time.append(to_string(hours)+":");
	
	if (minutes < 10) time.append("0");
	time.append(to_string(minutes)+":");
				
	if (seconds < 10) time.append("0");
	time.append(to_string(seconds));

	return date+" "+time;
}

// Transform date written from file into vector of numbers
vector<double> readDateTime(string line){
	vector<double> date;
	stringstream ss(line);
	string s;
	vector<string> words;

    istream_iterator<std::string> begin(ss);
    istream_iterator<std::string> end;
    vector<std::string> word(begin, end);
	//cout << line << endl;
    for (auto &s: word) {
        //std::cout << s << std::endl;
		date.push_back(stod(s));
	}
	// Return: Vector containing the Epoch Info
	return date;
}

// Function to convert date to Julian Date
double convertToJulianDate(int year, int month, int day) {
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    return day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045.5;
}


void fixTimeStamp(double* hr, double* min, double* sec){

	if ((*min >= 30) || ((*min == 29) && (*sec >= 30))) {
        // More than or equal to 30 minutes, round up to the next hour
        *hr += 1;
		
    }
	*min = 0;
	*sec = 0;

	return;
}

// Source: BOOK called GPS Theory Algorithm & Applications by Guochang Xu (Pg 18-20)
double gpsTime(std::vector<double> epochInfo) {

	// As precaution, check if we have required epoch info
	if (epochInfo.size() < 6) {
		// If so, we dont have enough epoch info to work with
		// We need to exit prematurely
		exit(1);
	}
	// Year Month Day of epoch
	double y = epochInfo.at(0); double m = epochInfo.at(1); double d = epochInfo.at(2);
	// Hour Minute Second of epoch
	double hr = epochInfo.at(3); double min = epochInfo.at(4); double sec = epochInfo.at(5);
	// UTC time in hours
  	// second
  	//fixTimeStamp(&hr, &min, &sec);
	//cout << hr << min << sec << endl;
  	double UTC = hr*3600 + min*60. + sec;

	// Taking care of month and year conditioning
	if (m <= 2) {

		y = y  - 1;
		m = m + 12;

	}
	// Julian Date
	double jDate = convertToJulianDate(y,m,d);

  	double WN = floor((jDate - JD_GPS_epoch) / DAYINWEEK);
  	double TOW = fmod((jDate - JD_GPS_epoch), DAYINWEEK) * SECPERDAY + UTC;
	// Return GPS Time in seconds
  	double gpsTime = (WN * SECPERWEEK) + TOW + numLeapSeconds;

	return gpsTime;
	
}

// A function used to replace a character ch1 in a string with another character ch2
string replaceChars(string str, char ch1, char ch2) {
	for (unsigned i = 0; i < str.length(); ++i) {
		if (str[i] == ch1)
			str[i] = ch2;
	}
	return str;
}

// A function to erase substring (first occurence of) from main string
void eraseSubStr(string & mainStr, const string & toErase) {
	// Search for the substring in string
	size_t pos = mainStr.find(toErase);
	if (pos != string::npos){
		// If found then erase it from string
		mainStr.erase(pos, toErase.length());
	}
}

// A function to define time in HH:MM:SS format
string HHMMSS(double hours, double mins, double secs) {
	stringstream ss;
	ss << std::setw(2) << std::setfill('0') << (int)hours << ":";
	ss << std::setw(2) << std::setfill('0') << (int)mins << ":";
	ss << std::setw(2) << std::setfill('0') << (int)secs;
	string hms = ss.str();
	return hms;
}

// Funzione per determinare se l'anno è bisestile
bool isLeapYear(int year) {
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

// Funzione per ottenere la data e l'ora da un numero del giorno dell'anno
std::vector<double> dayOfYearToDateTime(int day_of_year, int year) {
    int days_in_month[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    // Check if lap year
    if (isLeapYear(year)) {
        days_in_month[2] = 29; 
    }
	// Check on validity
    if (day_of_year < 1 || day_of_year > 365 + isLeapYear(year)) {
		perror("ERROR : DAY OF THE YEAR not ALLOWED!");
        //std::cerr << "Numero del giorno non valido." << std::endl;
        return {};
    }
    int month, day;
    // Trova il mese e il giorno del mese
    int cumulative_days = 0;
    for (month = 1; month <= 12; ++month) {
        cumulative_days += days_in_month[month];
        if (day_of_year <= cumulative_days) {
            day = day_of_year - (cumulative_days - days_in_month[month]);
            break;
        }
    }
    double hour = 0;
    double minute = 0;
    double second = 0;
    std::vector<double> datetime {static_cast<double>(year), static_cast<double>(month), static_cast<double>(day), hour, minute, second};
    return datetime;
}

/******************	FUNZIONI FILE ***********************************************************************************************/

// File Opener --> Initiates File pointer safely
void fileSafeIn(string filename, ifstream &fin)
{
	// Use input file stream, detect file
	ifstream inputfile(filename);
	// Check if file can be opened
	if (!inputfile.is_open())
	{
		perror("Error while opening file");
	}
	// Check is file can be read
	else if (inputfile.bad())
	{
		perror("Error while reading file");
	}
	else
	{
		// If no errors, we safely open the file...
		fin.open(filename);
	}
}

// A function to initialize the updated output datafile
void fileSafeOut(string output_filename, ofstream &fout)
{
	// Creating NEW file
	fout.open(output_filename);
	fout.close();
	fout.open(output_filename, ios::app);
	// Verifica se l'apertura del file ha avuto successo
    if (!fout.is_open()) {
        cerr << "Errore durante l'apertura del file di output: " << output_filename << endl;
        // Gestire l'errore appropriatamente, ad esempio lanciando un'eccezione
    }
}

void firstCSVLine(ofstream &fout,int flag){
	if(flag == 0)
		fout << "SV;DATETIME;X;Y;Z" << endl;
	if(flag == 1)
		fout << "SV;DATETIME;X;Y;Z;Xdot;Ydot;Zdot;dts" << endl;
	if(flag == 2)
		fout << "SV;DATETIME;X;Y;Z;CLK" << endl;
}

void writeCsvLine(ofstream &fout,string name,string datetime,double xk,double yk,double zk){
	fout << name << ";" << datetime << ";" << setprecision(12) << xk << ";" << yk << ";" << zk << endl; 
}
void writeCsvLine(ofstream &fout,string name,string datetime,double xk,double yk,double zk,double xk_dot,double yk_dot,double zk_dot,double dts){
	fout << name << ";" << datetime << ";" << setprecision(12) << xk << ";" << yk << ";" << zk << ";" << xk_dot << ";" << yk_dot << ";" << zk_dot << ";" << dts << endl;
}
void writeCsvLine(ofstream &fout,string name,string datetime,double xk,double yk,double zk,double clock){
	fout << name << ";" << datetime << ";" << setprecision(12) << xk << ";" << yk << ";" << zk << ";" << clock << endl; 
}

double getRcm(){
	return 6.3695e+06;
}
