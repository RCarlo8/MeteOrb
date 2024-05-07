#include <iostream>
#include <armadillo>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;
using namespace arma;

/* Class work stores and manages observation of a specific receiver in order to store observation data and compute atmosphere delays */
class work
{
private:
    /* data */
    string rec_name,fnameAz,fnameEL,fnameIono,fnameSow,fnameTime;
    double lat,lon,h_ortho,rcm,lat_rad,lon_rad;
    mat el,az,el_rad,az_rad;
    vec ionoparams,sow,time;
    
public:
    work(bool a) {

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

    string get_name(){
        return this->rec_name;
    }
    
    mat get_el(){
        return this->el;
    }
    mat get_az(){
        return this->az;
    }
    mat get_el_rad(){
        return this->el_rad;
    }
    mat get_az_rad(){
        return this->az_rad;
    }
    vec get_iono(){
        return this->ionoparams;
    }
    vec get_sow(){
        return this->sow;
    }
    vec get_time(){
        return this->time;
    }
    double get_lat(){
        return this->lat;
    }
    double get_lon(){
        return this->lon;
    }
    double get_lat_rad(){
        return this->lat_rad;
    }
    double get_lon_rad(){
        return this->lon_rad;
    }
    double get_h(){
        return this->h_ortho;
    }
    double getRcm(){
        return this->rcm;
    }


    private :
    
    double compute_rad(){
        this->lat_rad = lat * datum::pi/180;
        this->lon_rad = lon * datum::pi/180;
        this->el_rad = el * datum::pi/180;
        this->az_rad = az * datum::pi/180;
    }

};

