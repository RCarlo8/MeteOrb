//
//

#ifndef ATMO_H
#define ATMO_H

#include<iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <armadillo>



using namespace std;
using namespace arma;

//double angle;
//#define sind(angle) (sin(fmod((x),360) * M_PI / 180));
//#define cosd(angle) (sin(fmod((x),360) * M_PI / 180));
const double NaN = datum::nan;
const double pi = datum::pi;


class Atmo
{

public:
    Atmo(/* args */);
    ~Atmo();


/**
 * Compute ionospheric error correction for multiple satellites in a single epoch of observations according 
 * to the Klobuchar model
 * 
 * @param
 * 
 *   [delay] = Atmosphere. klobuchar_model(lat, lon, az, el, sow, ionoparams)
 *
 * INPUT:
 *  
 *  lat = receiver latitude          [degrees] [nx1]
 *  lon = receiver longitude         [degrees] [nx1]
 *  az  = satellite azimuth          [degrees] [nx1]
 *  el = satellite elevation         [degrees] [nx1]
 *  sow = second of week                       [nx1]
 * 
 * OUTPUT:
 *
 *  delay = tropospheric error correction      [nx1]
 *
 * 
 * Algorithm from Leick, A. (2004) "GPS Satellite Surveying - 2nd Edition"
 * John Wiley & Sons, Inc., New York, pp. 301-303)*/
  
vec klobuchar(double lat, double lon, vec az, vec el, double sow, vec ionoparams);

/**
 * Compute the pseudorange correction due to troposheric refraction.
 * Saastamoinen algorithm using standard atmosphere accounting for height gradient for temperature pressure and humidity.
 *
 * 
 * @param
 * 
 *   [delay] = Atmosphere. saastamoinen(el,h)
 *
 * INPUT:
 *                
 *  h  = receiver orthometrical height [meters]
 *  el = satellite elevation           [rad] [nx1]
 *  
 * 
 * OUTPUT:
 *
 *  delay = ionospheric error correction      [nx1]
 *
 **/

vec saastamoinen(vec el, double h);

/**
 * Compute the Global Mapping Functions GMF
 * Reference: Boehm, J., A.E. Niell, P. Tregoning, H. Schuh (2006),
 * Global Mapping Functions (GMF): A new empirical mapping function based on numerical weather model data,
 * Geoph. Res. Letters, Vol. 33, L07304, doi:10.1029/2005GL025545.
 *
 * 
 * @param
 * 
 *   [gmfw,gmfh] = Atmosphere.gmf(el,h)
 *
 * INPUT:
 *
 *  timeMJD = modified Julian Date   [nx1]
 *  lat = receiver latitude          [rad] [nx1]
 *  lon = receiver longitude         [rad] [nx1]            
 *  htg  = receiver orthometrical height [meters]
 *  el = satellite elevation         [rad] [nx1]
 *  
 * 
 * OUTPUT:
 *
 *  delay = ionospheric error correction      [nx1]
 *
 **/

tuple<vec,vec> gmf(double timeMJD, double lat, double lon, double hgt, vec el);

  /** Get the iono mapping function: 
    * 
    * INPUT :
    *   
    *   lat_rad     ->      latitude of the receiver           [rad]
    *   h_ortho     ->      orthometric height of the receiver [m]
    *   el_rad      ->      elevation of the satellites        [rad]
    *   rcm         ->      meridian radius curvature          <optional>
    
    * OUTPUT :
    * 
    *   iono_mf     ->      iono mapping function
    *
    * SOURCES :
    *      [1] Handbook of Global Navigation Satellite System (2017)
    *
    * SYNTAX :
    *      [iono_mf] = getIonoMF(lat_rad, h_ortho, el_rad, rcm)
    **/

  vec getIonoMF(double lat, double h, vec el, double rcm);
  
/******+*/
  double pres,temp;
  void gpt(double mjd, double dlat, double dlon, double dhgt);
  vec saastamoinenModelGPT(double gps_time, double lat, double lon, double h, double undu, arma::vec el);
  vec saastamoinenModelPTH(double gps_time, double lat, double lon, double h, double undu, arma::vec el, double P, double T, double H) ;
  vec saast_wet(arma::vec T, double H, double h);
  vec saast_dry(arma::vec P, double h, double lat);
//////////////////////////////////////////////////////////////////////

private : 

    double getThinShellHeight();
    vec mfContinuedFractionForm(vec a,double b,vec c,vec el);
    vec mfContinuedFractionForm(vec a,double b,double c, vec el);
    vec mfContinuedFractionForm(double a,double b,double c,vec el);
    vec mfContinuedFractionForm();
    vec hydrostaticMFHeigthCorrection(double h_ell, vec el);


    /* data */
    double STD_TEMP = 291.15;
    double STD_PRES = 1013.25;
    double STD_HUMI = 50;

    bool cached = false;


    struct {
        vector<double> data;//       [], ...    % ionosphere single layer map [n_lat x _nlon x n_time]
        vector<double> first_lat; //,  [], ...    % first latitude
        vector<double> first_lon; // [], ...    % first longitude
        vector<double> d_lat;   //   [], ...    % lat spacing
        vector<double> d_lon;   //   [], ...    % lon_spacing
        vector<double> n_lat;   //   [], ...    % num lat
        vector<double> n_lon; //      [], ...    % num lon
        vector<double> first_time; //[], ...    % times [time] of the maps
        vector<double> first_time_double; // [], ...    % times [time] of the maps [seconds from GPS zero]
        vector<double> dt; //        [], ...    % time spacing
        vector<double> n_t;        //[], ...    % num of epocvhs
        vector<double> height; //     []  ...    % heigh of the layer
    } ionex;

    // variables used for gmf
    mat V,W,P;
    //
    double lat,lon,undu,ahm,aha,awm,awa,apm, apa, atm, ata;
    //
    mat delay;
};


#endif


