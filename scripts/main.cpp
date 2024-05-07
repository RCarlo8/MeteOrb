#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <armadillo>

#include "atmo.h"
#include "work.h"
#include <tuple>

using namespace std;
using namespace arma;

int main() {
    // Data
    time_t start,end;
    int scelta,row;
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
    cout << endl << " ::::::::::::::::::::::::::::::";
    cout << endl << " :: 1 Automatic load         ::";
    cout << endl << " :: 2 Manual load            ::";
    cout << endl << " ::::::::::::::::::::::::::::::";
    cout << endl << " : Select an operation: > "; cin >> scelta;


    work rec = work((scelta == 1));

    lat = rec.get_lat();
    lon = rec.get_lon();
    el = rec.get_el();
    az = rec.get_az();
    sow = rec.get_sow();
    ionoparams = rec.get_iono();
    h_ortho = rec.get_h();



    do { // Inizio Ciclo D-WHILE
    
    //system("clear"); // Pulisci La Schermata

    cout << endl << " .::::::::::::::::::::::::::.";
    cout << endl << " :: Atmospheric  functions ::";
    cout << endl << " ::::::::::::::::::::::::::::";
    cout << endl << " :: 1 Klobuchar            ::";
    cout << endl << " :: 2 Saastamoinen         ::";
    cout << endl << " :: 3 gmf                  ::";
    cout << endl << " :: 4 gmfIono              ::";
    cout << "\n" << " ::::::::::::::::::::::::::::";
    cout << endl << " :: 5 End                  ::";
    cout << "\n" << " ::::::::::::::::::::::::::::";
   
    cout << endl << " : Select an operation: > "; cin >> scelta;

        if(scelta==1) { 

            cout << endl << " .::   Klobuchar Model     ::.";
            
            cout << endl << " : Select one option : ";
            cout << endl << " ::::::::::::::::::::::::::::::::::::::";
            cout << endl << " : 1. Compute single epoch delay      :";
            cout << endl << " : 2. Compute delay for all the epoch :";
            cout << endl << " ::::::::::::::::::::::::::::::::::::::";
            cout << endl << " : Select an operation: > "; cin >> scelta;
            
            if(scelta == 1 ){
                cout << endl << " : Select an epoch from 1 to " << rec.get_sow().n_elem << ": > "; cin >> row;
                row--;
                timer.tic();

                vec mono_delay = atmosphere.klobuchar(lat,lon,az.col(row),el.col(row),sow(row),ionoparams);
                cout << "time taken for sat operation = " << timer.toc() << endl;
                // save operation
                mono_delay.save("output/ATMO/"+rec.get_name()+"_Klobuchar.csv",csv_ascii);

            } else if(scelta == 2) {
                mat delay(32,2880);

                timer.tic();

                for (int epoch = 0; epoch < rec.get_sow().n_elem; epoch++)
                {   
                    delay.col(epoch) = atmosphere.klobuchar(lat,lon,az.col(epoch),el.col(epoch),sow(epoch),ionoparams);
                    //delay.insert_cols(epoch, atmosphere.klobuchar(lat,lon,az.col(epoch),el.col(epoch),sow(epoch),ionoparams));
                }
                //save operation
                cout << "time taken for sat operation = " << timer.toc() << endl;
                cout << "time taken for single operation = " << timer.toc()/ double(rec.get_el().n_cols) << endl;
                delay.save("output/ATMO/"+rec.get_name()+"_Klobuchar.csv",csv_ascii);

            } else {
                cout << " : Ops! something went wrong! : " << endl;
            }

            scelta = 1;
            system("pause>nul");   
}
      if(scelta==2) {  // Saastamoinen :: el ,h_orto

            cout << endl << " .::   Saastamoinen    ::.";
            cout << endl << " : Select one option : ";
            cout << endl << " ::::::::::::::::::::::::::::::::::::::";
            cout << endl << " : 1. Compute single epoch delay      :";
            cout << endl << " : 2. Compute delay for all the epoch :";
            cout << "\n" << " ::::::::::::::::::::::::::::::::::::::";
            cout << endl << " : Select an operation: > "; cin >> scelta;

            if(scelta == 1 ){
                cout << endl << " : Select an epoch from 1 to " << rec.get_sow().n_elem << ": > "; cin >> row;
                row--;
                vec mono_delay = atmosphere.saastamoinen(rec.get_el_rad().col(row),rec.get_h());
                // save operation
                mono_delay.save("output/ATMO/"+rec.get_name()+"_Saast.csv",csv_ascii);

            } else if(scelta == 2) {
                
                mat delay(32,2880);
                timer.tic();
                for (int epoch = 0; epoch < delay.n_cols; epoch++)
                {   
                    delay.col(epoch) = atmosphere.saastamoinen(el.col(epoch),h_ortho);
                    //delay.insert_cols(epoch, atmosphere.saastamoinen(rec.get_el_rad().col(epoch),rec.get_h()));
                }
                cout << "time taken = " << timer.toc()  << endl;
                delay.save("output/ATMO/"+rec.get_name()+"_Saast.csv",csv_ascii);

            } else {
                cout << "Ops! something went wrong!" << endl;
            }

            system("pause>nul");
        }

      if(scelta==3) {  
        
        cout << endl << " ::     GMF     ::" << endl;

        mat gmfh;
        mat gmfw;
        vec gmf1,gmf2;

        timer.tic();

        for (int epoch = 0; epoch < rec.get_el().n_cols; epoch++) {
            
            //Inserire metodo accesso diretto epoch ?
            tie(gmf1,gmf2) = atmosphere.gmf(rec.get_time()(epoch),rec.get_lat_rad(),rec.get_lon_rad(),rec.get_h(),rec.get_el_rad().col(epoch));
            
            gmfh.insert_cols(epoch,gmf1);
            gmfw.insert_cols(epoch,gmf2);
            
            }
        cout << "time taken = " << timer.toc() / double(rec.get_el().n_cols) << endl;
        cout << "time taken = " << timer.toc() << endl;
        gmfh.save("output/ATMO/"+rec.get_name()+"_gmfh.csv",csv_ascii);
        gmfw.save("output/ATMO/"+rec.get_name()+"_gmfw.csv",csv_ascii);
        system("pause>nul");
       
       }

       if(scelta==4) {

        cout << endl << " ::  Iono MF  ::" << endl;

        mat ionoMF(32,2880);

        timer.tic();
        for (int epoch = 0; epoch < rec.get_el().n_cols; epoch++) {
            timer.tic();
            atmosphere.getIonoMF(rec.get_lat_rad(),rec.get_h(),rec.get_el_rad().col(epoch),rec.getRcm());
            cout << "time  = " << timer.toc() << endl;
            }

        cout << "time  = " << timer.toc() << endl;

        
        for (int epoch = 0; epoch < rec.get_el().n_cols; epoch++) {

            ionoMF.insert_cols(epoch, atmosphere.getIonoMF(rec.get_lat_rad(),rec.get_h(),rec.get_el_rad().col(epoch),rec.getRcm()));
            
            }
        
        cout << "time taken = " << timer.toc() / double(rec.get_el().n_cols) << endl;
        cout << "time  = " << timer.toc() << endl;
        ionoMF.save("output/ATMO/"+rec.get_name()+"_ionoMf.csv",csv_ascii);
        
        system("pause>nul");

       }

    } while(scelta!=5);

    cout << endl << " :: Thank you!! See you later! :: " << endl << endl;;
    return 0;
}