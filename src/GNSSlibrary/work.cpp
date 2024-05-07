#include <iostream>
#include <armadillo>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "work.h"
#include "utility"

using namespace std;
using namespace arma;

work::work(bool a) {

        if(a) {
            cout << " : Insert name of receiver : > "; cin >> rec_name;
            vec coord;
            coord.load("data/satellite_atmo/coord.csv",csv_ascii);
            lat = coord[0];
            lon = coord[1];
            h_ortho = coord[2];
            rcm = coord[3];
            el.load("data/satellite_atmo/el.csv",csv_ascii);
            az.load("data/satellite_atmo/az.csv",csv_ascii);
            ionoparams.load("data/satellite_atmo/ionoparams.csv",csv_ascii);
            sow.load("data/satellite_atmo/sow.csv",csv_ascii);
            time.load("data/satellite_atmo/time.csv",csv_ascii);
        
        } else {

            cout << " : Receiver position data : " << endl;
            cout << " : Insert values for latitude [deg] : > "; cin >> lat;
            cout << " : Insert values for longitude [deg] : > "; cin >> lon;
            cout << " : Insert values for orthometric heigth [deg] : > "; cin >> h_ortho;
            cout << " : Insert value for rcm [m] : > "; cin >>rcm;

            cout << " : Observations : " << endl;
            cout << " : Insert values from files that contains the receiver observations :" <<endl;
            cout << endl << " : Insert name of the file for azimuthal observations  >"; cin >> fnameAz;
            cout << endl << " : Insert name of the files for elevation observations >"; cin >> fnameEL;
            cout << endl << " : Insert name of the file for ionoparmas : >"; cin >> fnameIono;
            cout << endl << " : Insert name of the file for second of the week : >"; cin >> fnameSow;
            cout << endl << " : Insert name of file for time observations : > "; cin >> fnameTime; 

            // loading files
            el.load("data/satellite_atmo/"+fnameEL,csv_ascii);
            az.load("data/satellite_atmo/"+fnameAz,csv_ascii);
            ionoparams.load("data/satellite_atmo/"+fnameIono,csv_ascii);
            sow.load("data/satellite_atmo/"+fnameSow,csv_ascii);
            time.load("data/satellite_atmo/"+fnameTime,csv_ascii);

        }
        this->compute_rad();
       
    };

    string work::get_name(){
        return this->rec_name;
    }
    
    mat work::get_el(){
        return this->el;
    }
    mat work::get_az(){
        return this->az;
    }
    mat work::get_el_rad(){
        return this->el_rad;
    }
    mat work::get_az_rad(){
        return this->az_rad;
    }
    vec work::get_iono(){
        return this->ionoparams;
    }
    vec work::get_sow(){
        return this->sow;
    }
    vec work::get_time(){
        return this->time;
    }
    double work::get_lat(){
        return this->lat;
    }
    double work::get_lon(){
        return this->lon;
    }
    double work::get_lat_rad(){
        return this->lat_rad;
    }
    double work::get_lon_rad(){
        return this->lon_rad;
    }
    double work::get_h(){
        return this->h_ortho;
    }
    double work::getRcm(){
        return this->rcm;
    }


   
    
double work::compute_rad(){
        this->lat_rad = lat * datum::pi/180;
        this->lon_rad = lon * datum::pi/180;
        this->el_rad = el * datum::pi/180;
        this->az_rad = az * datum::pi/180;
    }