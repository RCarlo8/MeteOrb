#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <armadillo>

#include "atmo.h"
#include "work.h"
#include <tuple>
#include "atmo.h"
#include "utility.h"

using namespace std;
using namespace arma;


#define sind(x) (sin(fmod((x),360) * M_PI / 180));
#define cosd(x) (sin(fmod((x),360) * M_PI / 180));


#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class GPT {
private:
    double lon, lat;
    double apm, apa, atm, ata;
    mat P;
    double undu;

public:
    GPT() {}

    void gpt(double gps_time, double dlat, double dlon, double dhgt) {
        
        // Reference day is 28 January
        double doy = (gps_time / 86400 - 22) / 365.25;

        vec a_geoid[55] = { /* Values of a_geoid */ };
        vec b_geoid[55] = { /* Values of b_geoid */ };
        vec ap_mean[55] = { /* Values of ap_mean */ };
        vec bp_mean[55] = { /* Values of bp_mean */ };
        vec ap_amp[55] = { /* Values of ap_amp */ };
        vec bp_amp[55] = { /* Values of bp_amp */ };
        vec at_mean[55] = { /* Values of at_mean */ };
        vec bt_mean[55] = { /* Values of bt_mean */ };
        vec at_amp[55] = { /* Values of at_amp */ };
        vec bt_amp[55] = { /* Values of bt_amp */ };

        // Computing Legendre Polynomial
        double t = sin(dlat);

        int n = 9, m = 9;
        int nmax = 9;

        vec dfac(2 * n + 1);
        dfac(0) = 1;
        for (int i = 1; i <= 2 * n; i++)
            dfac(i) = dfac(i - 1) * i;

        
        // Determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62)
        mat P(n + 1, m + 1);
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= min(i, m); j++) {
                int ir = floor((i - j) / 2);
                double sum_t = 0;
                
                for (int k = 0; k <= ir; k++) {
                    sum_t += pow(-1, k) * dfac(2 * i - 2 * k + 1) / (dfac(k) * dfac(i - k) * dfac(i - j - 2 * k + 1)) *
                             pow(t, i - j - 2 * k);
                }
                // Legendre functions moved by 1
                P(i, j) = 1.0 / (pow(2, i) * sqrt(pow(1 - pow(t, 2), j)) * dfac(i - j)) * sum_t;
            }
        }

        // Spherical harmonics
        vec aP, bP;
        int idx = 0;
        for (int n = 0; n <= nmax; ++n) {
            for (int m = 0; m <= n; ++m) {
                idx++;
                aP.insert_rows(idx - 1, P.row(n) * cos(m * dlon));
                bP.insert_rows(idx - 1, P.row(n) * sin(m * dlon));
            }
        }

        if (undu == 0) {
            undu = sum(a_geoid) * aP + sum(b_geoid) * bP;
        }

        apm = sum(ap_mean) * aP + sum(bp_mean) * bP;
        apa = sum(ap_amp) * aP + sum(bp_amp) * bP;
        atm = sum(at_mean) * aP + sum(bt_mean) * bP;
        ata = sum(at_amp) * aP + sum(bt_amp) * bP;

        // Orthometric height
        double h_ort = dhgt - undu;

        // Height correction for pressure
        double pres0 = apm + apa * cos(doy * 2 * datum::pi);
        double pres = pres0 * pow(1 - 0.0000226 * h_ort, 5.225);

        // Height correction for temperature
        double temp0 = atm + ata * cos(doy * 2 * datum::pi);
        double temp = temp0 - 0.0065 * h_ort;

        cout << "Pressure: " << pres << " hPa" << endl;
        cout << "Temperature: " << temp << " Celsius" << endl;
        cout << "Undulation: " << undu << " m" << endl;
    }
};


/*

This subroutine determines Global Pressure and Temperature
based on Spherical Harmonics up to degree and order 9
            
input data
         ----------
        dmjd: modified julian date
        dlat: latitude in radians
        dlon: longitude in radians
        dhgt: ellipsoidal height in m
    
        output data
         -----------
        pres: pressure in hPa
        temp: temperature in Celsius
        undu: Geoid undulation in m (from a 9x9 EGM based model)
            
            
        
        Johannes Boehm, 2006 June 12
        rev 2006 June 16: geoid undulation is accounted for

        Reference:
        J. Boehm, R. Heinkelmann, and H. Schuh,
        Global Pressure and Temperature (GPT): A spherical harmonics expansion
        of annual pressure and temperature variations for geodetic applications,
        to be submitted to Journal of Geodesy, 2006
            
        reference day is 28 January
        this is taken from Niell (1996) to be consistent
*/



arma::vec saastamoinenModelGPT(this,double gps_time, double lat, double lon, double h, arma::vec undu, arma::vec el) {
    // Compute pressure and temperature using GPT model
    arma::vec pres, temp;
    this.gpt(pres, temp, gps_time, lat * arma::datum::pi / 180, lon * arma::datum::pi / 180, h, undu);

    // Compute tropospheric delays
    arma::vec t_h = h;
    t_h(undu > -300) -= undu(undu > -300);
    arma::vec ZHD_R = saast_dry(pres, t_h, lat);
    arma::vec ZWD_R = saast_wet(temp, STD_HUMI, t_h);

    // Apply GMF corrections
    arma::vec gmfh_R, gmfw_R;
    gmf(gmfh_R, gmfw_R, gps_time, lat * arma::datum::pi / 180, lon * arma::datum::pi / 180, h - undu, el * arma::datum::pi / 180);

    // Compute delay
    return gmfh_R % ZHD_R + gmfw_R % ZWD_R;
}

#include <armadillo>

arma::vec saastamoinenModelPTH(double gps_time, double lat, double lon, double h, arma::vec undu, arma::vec el, double P, double T, double H) {
    // Check if meteorological data is valid
    if (arma::vec({P, T, H}).has_nan()) {
        // Compute pressure and temperature using GPT model
        arma::vec pres, temp;
        gpt(pres, temp, gps_time, lat * arma::datum::pi / 180, lon * arma::datum::pi / 180, h, undu);

        // Set missing values
        if (isnan(P)) P = pres;
        if (isnan(T)) T = temp;
        if (isnan(H)) H = STD_HUMI;

        // Log warning
        //this->log.addWarning("No valid meteo data are present @ " + std::to_string(gps_time / 86400 + GPS_Time.GPS_ZERO) + "\nUsing standard GPT values\n - " + std::to_string(T) + " °C\n - " + std::to_string(P) + " hPa\n - humidity " + std::to_string(H), 100);
    }

    // Compute tropospheric delays
    arma::vec ZHD_R = saast_dry(P, h, lat);
    arma::vec ZWD_R = saast_wet(T, H, h);

    // Apply GMF corrections
    arma::vec gmfh_R, gmfw_R;
    gmf(gmfh_R, gmfw_R, gps_time, lat * arma::datum::pi / 180, lon * arma::datum::pi / 180, h - undu, el * arma::datum::pi / 180);

    // Compute delay
    return gmfh_R % ZHD_R + gmfw_R % ZWD_R;
}

/*  Saastamoinen wet
     
    DESCRIPTION: Zenith Wet Delay (ZWD) computation by Saastamoinen model.
    
    INPUT: 
    
        [ZWD] = saast_wet(T, H);
    
    INPUT:
        
        T = air temperature [in Celsius]
        H = humidity [valore himidity]
        h = heigth for correction

    OUTPUT:
       
        ZWD = Zenith Wet Delay
            
 */
arma::vec saast_wet(arma::vec T, arma::vec H, arma::vec h) {
    // Convert temperature from Celsius to Kelvin
    T += 273.15;
    
    // Perform height correction
    // (Comment out the following line if height correction is done before)
    // H %= arma::exp(-0.0006396 * h);
    
    // Convert humidity from percentage to fraction
    H /= 100.0;
    
    // Compute constant and exponent
    arma::vec c = -37.2465 + 0.213166 * T - 2.56908 * pow(10, -4) * arma::pow(T, 2);
    arma::vec e = H % arma::exp(c);

    // Compute Zenith Wet Delay (ZWD) using Saastamoinen model
    return 0.0022768 * (((1255.0 / T) + 0.05) % e);
}


/*  Saastamoinen dry
 *  
 *  DESCRIPTION:
 *  Zenith Hydrostatic Delay (ZHD) computation by Saastamoinen model.
 *   
 *  SYNTAX
 *       
 *      [ZHD] = saast_dry(P, h, lat);
 *
 *  INPUT:
 * 
 *      P = atmospheric pressure [hPa]
 *      h = orthometric height [m]
 *      lat = latitude [deg]
 *          
 *  OUTPUT:
 * 
 *      ZHD = Zenith Hydrostatic Delay
 * 
 */
arma::vec saast_dry(arma::vec P, arma::vec h, double lat) {
    // Compute Zenith Hydrostatic Delay (ZHD) using Saastamoinen model
    return 0.0022768 * P % (1 + 0.00266 * cosd((2 * lat)) + 0.00000028 * h);
    // Alternatively, for alternative formula:
    // return 0.0022767 * P / (1 - 0.00266 * arma::cosd(2 * lat) - 0.00000028 * h);
}


int main(int argc, char const *argv[])
{
    vec az = regspace(-180,0.5,180);
    vec el = regspace(0,0.5,90);
    
    int UTC;
    vec iono;
    vec time = regspace(0,2,22);
    cout << time << endl;
    double jDate = convertToJulianDate(2022,3,19);
    iono.load("data/satellite_atmo/ionoparams.csv");
    //cout << iono << endl;
    Atmo atmo = Atmo();

    double lat = (45 + (28 / 60) + (38.28 / 3600)); // degrees
    double lon =  (9 + (10 / 60) + (53.40 / 3600)); // degrees
    double h = 148.321; //ell

    cube F(el.n_elem,az.n_elem,time.n_elem);

    //matrix inizialization
    cube iono_effec2(el.n_elem,az.n_elem,time.n_elem);
    cube tropo_effect(el.n_elem,az.n_elem,time.n_elem);
    cube global_iono(el.n_elem,az.n_elem,time.n_elem);

    

    //time, elevation and azimuth cycle 
    for (size_t i = 0; i < time.size(); i++)
    {
        for (size_t j = 0; j < el.size(); j++)
        {
            for (size_t k = 0; k < az.size(); k++)
            {
                UTC = time(i)*3600;
                double TOW = fmod((jDate - JD_GPS_epoch), DAYINWEEK) * SECPERDAY + UTC;
                iono_effec2(j,k,i) = atmo.klobuchar(lat,lon,az.subvec(k,k),el.subvec(j,j),TOW,iono).at(0);
                tropo_effect(j,k,i) = atmo.saastamoinen(el.subvec(j,j),h).at(0);
                global_iono(j,k,i) = atmo.getIonoMF(lat,h,el.subvec(j,j),NULL).at(0);
                global_iono(j,k,i) = (datum::c_0) * global_iono(j,k,i) * F(j,k,i);
            }
            
        }
    }

     
    for (size_t i = 0; i < time.size(); i++)
    {
        iono_effec2.slice(i).save("/Users/carlocolli/Documents/MATLAB/Tesi/Test/Atmosphere/iono/iono_time_"+to_string(i)+".txt",csv_ascii);
        tropo_effect.slice(i).save("/Users/carlocolli/Documents/MATLAB/Tesi/Test/Atmosphere/tropo/tropo_time_"+to_string(i)+".txt",csv_ascii);
        global_iono.slice(i).save("/Users/carlocolli/Documents/MATLAB/Tesi/Test/Atmosphere/g_iono/g_iono_time_"+to_string(i)+".txt",csv_ascii);
    
    }
    iono_effec2.save("iono.data",arma_ascii);
    global_iono.save("ionoMF.data",arma_ascii);
    cout << lat << endl;
    cout << lon << endl;
     
    return 0;
}
