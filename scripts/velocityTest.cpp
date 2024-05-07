#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <armadillo>

#include "atmo.h"
#include "work.h"
#include <tuple>
#include "coresky.h"
#include "utility.h"

#include <filesystem>
#include <sys/stat.h>

using namespace std;
using namespace arma;

int main() {
    // Data
    time_t start,end;
    int scelta,row;
    // Timer
    wall_clock timer;
    // Atmo object
    Atmo atmosphere = Atmo();

    double lat,lon,h_ortho,rcm,lat_rad,lon_rad;
    mat el,az,el_rad,az_rad;
    vec ionoparams,sow,time;

    system("clear"); // Pulisci La Schermata

    
    // Create work object and load obs
    cout << " .::::::::::::::::::::::::::::.";
    cout << endl << " : Load receiver observations : ";
    cout << " :: Select an option to upload the data from receiver :: " << endl;
    cout << endl << " ::::::::::::::::::::::::::::::" << endl; 
    cout << endl << " :: 1 Automatic load         ::";
    cout << endl << " :: 2 Manual load            ::";
    cout << endl << " ::::::::::::::::::::::::::::::";
    cout << endl << " : Select an operation: > "; cin >> scelta;

    cout << endl << "Starting operations..." << endl;


    work rec = work((scelta == 1));

    lat = rec.get_lat();
    lon = rec.get_lon();
    lat_rad = rec.get_lat_rad();
    lon_rad = rec.get_lon_rad();
    el = rec.get_el();
    az = rec.get_az();
    sow = rec.get_sow();
    ionoparams = rec.get_iono();
    h_ortho = rec.get_h();

    int n_sat = el.n_rows;
    int t = el.n_cols;
    
    //timer.tic();

    mat delay(n_sat,t);
    
        /* code */
    
    timer.tic();
    for (int epoch = 0; epoch < rec.get_sow().n_elem; epoch++)
    {   
        delay.col(epoch) = atmosphere.klobuchar(lat,lon,az.col(epoch),el.col(epoch),sow(epoch),ionoparams);
        //delay.insert_cols(epoch, atmosphere.klobuchar(lat,lon,az.col(epoch),el.col(epoch),sow(epoch),ionoparams));        
    }
    
    cout << " : time taken for Klobuchar = " << timer.toc() << endl;//timer.toc() << endl;

    timer.tic();
    
    for (int epoch = 0; epoch < delay.n_cols; epoch++)
    {   
        delay.col(epoch) = atmosphere.saastamoinen(el.col(epoch),h_ortho);
                    //delay.insert_cols(epoch, atmosphere.saastamoinen(rec.get_el_rad().col(epoch),rec.get_h()));
    }
    cout << " : time taken for Saastamoinen = " << timer.toc()  << endl;

    // result matrix initialization
    mat gmfh(n_sat,t);
    mat gmfw(n_sat,t);
    vec gmf1,gmf2;

    timer.tic();

    for (int epoch = 0; epoch < rec.get_el().n_cols; epoch++) 
    {

        tie(gmf1,gmf2) = atmosphere.gmf(rec.get_time()(epoch),lat_rad,lon_rad,h_ortho,el.col(epoch));            
        gmfh.col(epoch) = gmf1;
        gmfw.col(epoch) = gmf2;
            
    }
    
    cout << " : time taken GMF = " << timer.toc() << endl;

    timer.tic();
    
    for (int epoch = 0; epoch < t; epoch++) 
    {
        delay.col(epoch) = atmosphere.getIonoMF(lat,h_ortho,el.col(epoch),rec.getRcm());    
    }

    cout << " : time for IonoMF  = " << timer.toc() << endl;

    /*          CORESKY TEST        */ 

    //Initialize coresky object
    // Select file path
    string file = "data/ephemerides/ESA0MGNFIN_20220320000_01D_05M_ORB.SP3";
    // Initialize core skyobject
    CoreSky core(file,{},ios::in);
    core.readHeader();
    core.readCoord();
    //
    timer.tic();
    core.computePolyCoef();
    cout << " : time for computing Pos Coeff  = " << timer.toc() << endl;

    //
    timer.tic();
    core.computeClockPolyCoef();
    cout << " : time for computing Clock Coeff  = " << timer.toc() << endl;

    
    vec gps_time = arma::regspace(core.getStartTime(),30,core.getEndTime()); 
    vec sat = arma::regspace(0,1,31);

    //cout << core.getEpoch() << endl;
    timer.tic();
    core.coordInterpolate(gps_time,sat);
    cout << " : time for interpolation  = " << timer.toc() << endl;

}   