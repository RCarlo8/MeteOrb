#include <iostream>
#include <stdio.h>
#include <armadillo>
#include "coresky.h"
//#include <filesystem>
#include <sys/stat.h>
#include "utility.h"

namespace fs = std::__fs::filesystem;
using namespace std;

/* 
* This script read all SP3 file in the data/ephemerides folder and compute interpolated coordinates 
* of the satellites each 30 seconds for the first day of each month from  February to December 2022
*/

int main(int argc, char const *argv[]) {

    //Select GNSS system
    // GALILEO e GPS
    vector<char> ss = {'G','E'};
    // Only GALILEO
    //vector<char> ss = {'E'};
    // Only GPS
    //vector<char> ss = {'G'};

    // Path to the directory
    string path
        = "data/ephemerides/";
 
    // directory
    struct stat sb;
    vector<string> files;
    // Looping until all the items of the directory are exhausted
    for (const auto& entry : fs::directory_iterator(path)) {
 
        // Converting the path to const char * in the subsequent lines
        std::__fs::filesystem::path outfilename = entry.path();
        std::string outfilename_str = outfilename.string();
        const char* path = outfilename_str.c_str();
        
        // Testing whether the path points to a non-directory or not If it does, displays path
        if (stat(path, &sb) == 0 && !(sb.st_mode & S_IFDIR))
            files.push_back(outfilename_str);
    } 
    // Sort file by name
    std::sort(files.begin(), files.end());

    // For each file coordinates are computed and stored in a csv file
    for (size_t i = 1; i < files.size(); i++) 
    {
        // Built coresky object using file
        CoreSky core(files[i],ss,ios::in);

        try
        {
            // Read operations
            core.readHeader();
            core.readCoord();
            // Compute polynomial coefficients 
            core.computePolyCoef();
            core.computeClockPolyCoef();
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        
        
        // Custom time interval
        int step = 1;
        vector<double> start = {2023.0, 5.0, 10.0, 6.0, 0.0, 0.0}; 
        vector<double> end = {2023,5,10,14.0,0.0,0.0}; 

        //vec time = regspace(core.getStartTime(),step,core.getEndTime());
        vec time = regspace(gpsTime(start),step,gpsTime(end));
        cout << "tempo da vettore : " << time(0) << endl;
        cout << "tempo da da input : " << gpsTime(start) << endl;
        vec sat = regspace(0,1,core.getNumOfsat()-1);
        
        vector<string> satIds = core.getSatIds();
        //return 0;
        cube X = core.coordInterpolate(time,sat);
        //cube clockCorr = core.getClock();

        string fileOut = "output/"+files[i].substr(17,34)+".csv";
        ofstream fout;
        // Creating NEW file for output
        fout.open(fileOut); fout.close();
        fout.open(fileOut,ios::app);
        vector<string> satId = core.getSatIds();
        
        firstCSVLine(fout,2);        
        for (int prn = 0; prn < sat.size(); prn++)
        {
            mat tmp = X.col_as_mat(sat[prn]);
        
            for(size_t i = 0; i < time.size(); i++)
            {   

                writeCsvLine(fout,satIds[sat[prn]],gpstime2Date(time(i)), tmp(i,0), tmp(i,1), tmp(i,2),tmp(i,3));

            }  
        }
        
        fout.close(); 
    }

    return 0;
}
