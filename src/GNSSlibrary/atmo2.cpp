#include <iostream>
#include <armadillo>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include "atmo.h"
#include "utility.h"


using namespace std;
using namespace arma;


Atmo::Atmo(){};
Atmo::~Atmo(){};

vec Atmo::klobuchar(double lat, double lon, vec az, vec el, double sow, vec ionoparams) {
    // define number of satellite from observation
    int n_sat = max(size(az));
    // initialize result
    vec delay;
    delay.zeros(n_sat);

    // light speed constant
    double V_LIGHT = datum::c_0;

    // ionospheric parameters
    double a0 = ionoparams(0);
    double a1 = ionoparams(1);
    double a2 = ionoparams(2);
    double a3 = ionoparams(3);
    double b0 = ionoparams(4);
    double b1 = ionoparams(5);
    double b2 = ionoparams(6);
    double b3 = ionoparams(7);

    // elevation from 0 to 90 degrees
    el = abs(el);
    //conversion to semicircles
    lat = lat / 180;
    lon = lon / 180;
    az = az / 180;
    el = el / 180;

    vec f = 1 + 16 * (pow((0.53 - el),3));
    vec psi = (0.0137 /(el + 0.11)) - 0.022;
    vec phi = lat + (psi % cos(az*pi));

    for (int i = 0; i < n_sat; i++)
    {
        if(phi(i) > 0.416) phi(i) = 0.416;
        if(phi(i) < -0.416) phi(i) = -0.416;
    }

    vec lambda = lon + ((psi % sin(az*pi)) / cos(phi*pi));
    vec ro = phi + 0.064 * cos((lambda-1.617)*pi);
    vec t = lambda*43200 + sow;
    t = t - 86400 * floor(t/86400);
    vec a = a0 + a1*ro + a2*pow(ro,2) + a3*pow(ro,3);

    for (int i = 0; i < n_sat; i++)
    {
        if(a(i) < 0) a(i) = 0;
    }
    
    vec p = b0 + b1*ro + b2*pow(ro,2) + b3*pow(ro,3);
    
    for (int i = 0; i < n_sat; i++)
    {
        if(p(i) < 72000) p(i) = 72000;
    }

    vec x = (2*pi*(t - 50400)) / p;
    
    //ionospheric delay

    for (int i = 0; i < n_sat; ++i)
        {
            if(abs(x(i)) < 1.57) {
            delay(i) = V_LIGHT * f(i) * (5e-9 + a(i) * (1- (pow(x(i),2)/2 + pow(x(i),4)/24)));
            }

            if(abs(x(i)) >= 1.57) {
            delay(i) = V_LIGHT * f(i) * 5e-9;
            }
        }
    return delay;
}

vec Atmo::saastamoinen(vec el, double h){
    
    double STD_TEMP = 291.15;
    double STD_PRES = 1013.25;
    double STD_HUMI = 50;

    int n_sat = max(size(el));
    vec delay(n_sat);

    if (h < 5000)
    {
        el = abs(el);
        double P = STD_PRES * pow((1-0.0000226*h),5.225);
        double T = STD_TEMP - 0.0065*h;
        double H = STD_HUMI * exp(-0.0006396*h);
      
    
        vec h_a = {0, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000};
        vec B_a = {1.156, 1.079, 1.006, 0.938, 0.874, 0.813, 0.757, 0.654, 0.563};

        double t;
        double B;
        vec index(2);
        int a;
        int b;
      
        vec d = h_a - h;
        int j = index_min(abs(d)); 

        if(d(j) > 0) { // a prec e b succ
            a = j-1;
            b = j;
        } else {
            a = j;
            b = j+1;
        }

        t = (h - h_a(a)) / (h_a(b) - h_a(a));
        B = (1-t)*B_a(a) + t * B_a(b);

        double e = 0.01 * H * exp(-37.2465 + 0.213166*T - 0.000256908*pow(T,2));
        vec w_fun = (1 - pow(tan(el),2) / (pow(tan(el),2) + pow(tan(el+2),2)));

        for (int i = 0; i < n_sat; ++i){
            
            delay[i] = ((0.002277 / sin(el[i])) * (P - (B / (pow(tan(el[i]),2) +0.01*w_fun[i]))) + (0.002277 / sin(el[i])) * (1255 / T + 0.05) * e );
        }
        
    } else { 
        
        delay.zeros(max(size(el)));
    }
     
    return delay;
}

tuple<vec,vec> Atmo::gmf(double timeMJD, double lat, double lon, double hgt, vec el) {
    
    vec gmfh,gmfw;
    field<vec> result(2);

    double doy = timeMJD  - 44239 + 1;
    //pi = 3.14159265359;

    bool cached = (lon == this->lon) && (lat == this->lat) && (datum::nan == this->ahm) && (datum::nan==this->aha);

    if(cached) {

        mat W = this->W;
        mat V = this->V;

        double ahm = this->ahm;
        double aha = this->aha; 
        double awm = this->awm;
        double awa = this->awa;

    } else {

        vec ah_mean = {     125.170000000000,0.850300000000000,0.0693600000000000,-6.76000000000000,0.177100000000000,
                        0.0113000000000000,0.596300000000000,0.0180800000000000,0.00280100000000000,-0.00141400000000000,
                        -1.21200000000000,0.0930000000000000,0.00368300000000000,0.00109500000000000,4.67100000000000e-05,
                        0.395900000000000,-0.0386700000000000,0.00541300000000000,-0.000528900000000000,0.000322900000000000,
                        2.06700000000000e-05,0.300000000000000,0.0203100000000000,0.00590000000000000,0.000457300000000000,
                        -7.61900000000000e-05,2.32700000000000e-06,3.84500000000000e-06,0.118200000000000,0.0115800000000000,
                        0.00544500000000000,6.21900000000000e-05,4.20400000000000e-06,-2.09300000000000e-06,1.54000000000000e-07,
                        -4.28000000000000e-08,-0.475100000000000,-0.0349000000000000,0.00175800000000000,0.000401900000000000,
                        -2.79900000000000e-06,-1.28700000000000e-06,5.46800000000000e-07,7.58000000000000e-08,-6.30000000000000e-09,
                        -0.116000000000000,0.00830100000000000,0.000877100000000000,9.95500000000000e-05,-1.71800000000000e-06,
                        -2.01200000000000e-06,1.17000000000000e-08,1.79000000000000e-08,-1.30000000000000e-09,1.00000000000000e-10,
    }; 


    vec bh_mean = {     0,	0,	0.0324900000000000,	0,	0.0332400000000000,	
                        0.0185000000000000,	0,	-0.111500000000000,	0.0251900000000000,	0.00492300000000000,
                        0,	0.0273700000000000,	0.0159500000000000,	-0.000733200000000000,	0.000193300000000000,
                        0,	-0.0479600000000000,	0.00638100000000000,	-0.000159900000000000,	-0.000368500000000000,
                        1.81500000000000e-05,	0,	0.0703300000000000,	0.00242600000000000,	-0.00111100000000000,
                        -0.000135700000000000,	-7.82800000000000e-06,	2.54700000000000e-06,	0,	0.00577900000000000,
                        0.00313300000000000,	-0.000531200000000000,	-2.02800000000000e-05,	2.32300000000000e-07,	-9.10000000000000e-08,
                        -1.65000000000000e-08,	0,	0.0368800000000000,	-0.000863800000000000,	-8.51400000000000e-05,
                    	-2.82800000000000e-05,	5.40300000000000e-07,	4.39000000000000e-07,	1.35000000000000e-08,	1.80000000000000e-09,	
                        0,	-0.0273600000000000,	-0.000297700000000000,	8.11300000000000e-05,	2.32900000000000e-07,
                        8.45100000000000e-07,	4.49000000000000e-08,	-8.10000000000000e-09,	-1.50000000000000e-09,	2.00000000000000e-10
    };


    vec ah_amp = {      -0.273800000000000,	-2.83700000000000,	0.0129800000000000,	-0.358800000000000,	0.0241300000000000,
    	                0.0342700000000000,	-0.762400000000000,	0.0727200000000000,	0.0216000000000000,	-0.00338500000000000,
        	            0.442400000000000,	0.0372200000000000,	0.0219500000000000,	-0.00150300000000000,	0.000242600000000000,
            	        0.301300000000000,	0.0576200000000000,	0.0101900000000000,	-0.000447600000000000,	6.79000000000000e-05,
                	    3.22700000000000e-05,	0.312300000000000,	-0.0353500000000000,	0.00484000000000000,	3.02500000000000e-06,
                    	-4.36300000000000e-05,	2.85400000000000e-07,	-1.28600000000000e-06,	-0.672500000000000,	-0.0373000000000000,
                        0.000896400000000000,	0.000139900000000000,	-3.99000000000000e-06,	7.43100000000000e-06,	-2.79600000000000e-07,
                        -1.60100000000000e-07,	0.0406800000000000,	-0.0135200000000000,	0.000728200000000000,	9.59400000000000e-05,
                        2.07000000000000e-06,	-9.62000000000000e-08,	-2.74200000000000e-07,	-6.37000000000000e-08,	-6.30000000000000e-09,
                        0.0862500000000000,	-0.00597100000000000,	0.000470500000000000,	2.33500000000000e-05,	4.22600000000000e-06,
    	                2.47500000000000e-07,	-8.85000000000000e-08,	-3.60000000000000e-08,	-2.90000000000000e-09,	0
    };


    vec bh_amp = {      0,	0,	-0.113600000000000,	0,	-0.186800000000000,
    	                -0.0139900000000000,	0,	-0.104300000000000,	0.0117500000000000,	-0.00224000000000000,
        	            0,	-0.0322200000000000,	0.0133300000000000,	-0.00264700000000000,	-2.31600000000000e-05,
            	        0,	0.0533900000000000,	0.0110700000000000,	-0.00311600000000000,	-0.000107900000000000,
                	    -1.29900000000000e-05,	0,	0.00486100000000000,	0.00889100000000000,	-0.000644800000000000,
                    	-1.27900000000000e-05,	6.35800000000000e-06,	-1.41700000000000e-07,	0,	0.0304100000000000,
                        0.00115000000000000,	-0.000874300000000000,	-2.78100000000000e-05,	6.36700000000000e-07,	-1.14000000000000e-08,
                        -4.20000000000000e-08,	0,	-0.0298200000000000,	-0.00300000000000000,	1.39400000000000e-05,
                        -3.29000000000000e-05,	-1.70500000000000e-07,	7.44000000000000e-08,	2.72000000000000e-08,	-6.60000000000000e-09,
                        0,	0.0123600000000000,	-0.000998100000000000,	-3.79200000000000e-05,	-1.35500000000000e-05,
                        1.16200000000000e-06,	-1.78900000000000e-07,	1.47000000000000e-08,	-2.40000000000000e-09,	-4.00000000000000e-10
    };           
                
    
    vec aw_mean = {     56.4000000000000,	1.55500000000000,	-1.01100000000000,	-3.97500000000000,	0.0317100000000000,
                        0.106500000000000,	0.617500000000000,	0.137600000000000,	0.0422900000000000,	0.00302800000000000,
                        1.68800000000000,	-0.169200000000000,	0.0547800000000000,	0.0247300000000000,	0.000605900000000000,	
                        2.27800000000000,	0.00661400000000000, -0.000350500000000000,	-0.00669700000000000,	0.000840200000000000,	
                        0.000703300000000000,	-3.23600000000000,	0.218400000000000,	-0.0461100000000000,	-0.0161300000000000,
                        -0.00160400000000000,	5.42000000000000e-05,	7.92200000000000e-05,	-0.271100000000000,	-0.440600000000000,	
                        -0.0337600000000000,	-0.00280100000000000,	-0.000409000000000000,	-2.05600000000000e-05,	6.89400000000000e-06,	
                        2.31700000000000e-06,	1.94100000000000,	-0.256200000000000,	0.0159800000000000,	0.00544900000000000,	
                        0.000354400000000000,	1.14800000000000e-05,	7.50300000000000e-06,	-5.66700000000000e-07,	-3.66000000000000e-08,	
                        0.868300000000000,	-0.0593100000000000,	-0.00186400000000000,	-0.000127700000000000,	0.000202900000000000,	
                        1.26900000000000e-05,	1.62900000000000e-06,	9.66000000000000e-08,	-1.01500000000000e-07,	-5.00000000000000e-10
    };
                

    vec bw_mean = {     0,	0,	0.259200000000000,	0,	0.0297400000000000,
    	                -0.547100000000000,	0,	-0.592600000000000,	-0.103000000000000,	-0.0156700000000000,
        	            0,	0.171000000000000,	0.0902500000000000,	0.0268900000000000,	0.00224300000000000,
            	        0,	0.343900000000000,	0.0240200000000000,	0.00541000000000000,	0.00160100000000000,
                	    9.66900000000000e-05,	0,	0.0950200000000000,	-0.0306300000000000,	-0.00105500000000000,
                    	-0.000106700000000000,	-0.000113000000000000,	2.12400000000000e-05,	0,	-0.312900000000000,
                        0.00846300000000000,	0.000225300000000000,	7.41300000000000e-05,	-9.37600000000000e-05,	-1.60600000000000e-06,
                        2.06000000000000e-06,	0,	0.273900000000000,	0.00116700000000000,	-2.24600000000000e-05,
                        -0.000128700000000000,	-2.43800000000000e-05,	-7.56100000000000e-07,	1.15800000000000e-06,	4.95000000000000e-08,
                        0,	-0.134400000000000,	0.00534200000000000,	0.000377500000000000,	-6.75600000000000e-05,
                        -1.68600000000000e-06,	-1.18400000000000e-06,	2.76800000000000e-07,	2.73000000000000e-08,	5.70000000000000e-09

    };           


    vec aw_amp = {      0.102300000000000,	-2.69500000000000,	0.341700000000000,	-0.140500000000000,	0.317500000000000,
    	                0.211600000000000,	3.53600000000000,	-0.150500000000000,	-0.0166000000000000,	0.0296700000000000,	
                        0.381900000000000,	-0.169500000000000,	-0.0744400000000000,	0.00740900000000000,	-0.00626200000000000,	
                        -1.83600000000000,	-0.0175900000000000,	-0.0625600000000000,	-0.00237100000000000,	0.000794700000000000,	
                        0.000150100000000000,	-0.860300000000000,	-0.136000000000000,	-0.0362900000000000,	-0.00370600000000000,	
                        -0.000297600000000000,	1.85700000000000e-05,	3.02100000000000e-05,	2.24800000000000,	-0.117800000000000,	
                        0.0125500000000000,	0.00113400000000000,	-0.000216100000000000,	-5.81700000000000e-06,	8.83600000000000e-07,	
                        -1.76900000000000e-07,	0.731300000000000,	-0.118800000000000,	0.0114500000000000,	0.00101100000000000,	
                        0.000108300000000000,	2.57000000000000e-06,	-2.14000000000000e-06,	-5.71000000000000e-08,	2.00000000000000e-08,	
                        -1.63200000000000,	-0.00694800000000000,	-0.00389300000000000,	0.000859200000000000,	7.57700000000000e-05,	
                        4.53900000000000e-06,	-3.85200000000000e-07,	-2.21300000000000e-07,	-1.37000000000000e-08,	5.80000000000000e-09
    };           
                

    vec bw_amp = {      0,	0,	-0.113600000000000,	0,	-0.186800000000000,
    	                -0.0139900000000000,	0,	-0.104300000000000,	0.0117500000000000,	-0.00224000000000000,
                        0,	-0.0322200000000000, 0.0133300000000000,	-0.00264700000000000,	-2.31600000000000e-05,	
                        0,	0.0533900000000000,	0.0110700000000000,	-0.00311600000000000,	-0.000107900000000000,
            	        -1.29900000000000e-05,	0,	0.00486100000000000,	0.00889100000000000,	-0.000644800000000000,
                	    -1.27900000000000e-05,	6.35800000000000e-06,	-1.41700000000000e-07,	0,	0.0304100000000000,
                    	0.00115000000000000,	-0.000874300000000000,	-2.78100000000000e-05,	6.36700000000000e-07,	-1.14000000000000e-08,
                        -4.20000000000000e-08,	0,	-0.0298200000000000,	-0.00300000000000000,	1.39400000000000e-05,
                        -3.29000000000000e-05,	-1.70500000000000e-07,	7.44000000000000e-08,	2.72000000000000e-08,	-6.60000000000000e-09,
                        0,	0.0123600000000000,	-0.000998100000000000,	-3.79200000000000e-05,	-1.35500000000000e-05,
                        1.16200000000000e-06,	-1.78900000000000e-07,	1.47000000000000e-08,	-2.40000000000000e-09,	-4.00000000000000e-10
        
    };   

    //degree n and order m
    int nmax = 9;
                
                
    // unit vector
    double x = cos(lat)*cos(lon);
    double y = cos(lat)*sin(lon);
    double z = sin(lat);
                
    mat V = zeros(nmax + 1,nmax + 1);
    mat W = zeros(nmax + 1,nmax + 1);


    // Legendre polynomials
    V(0,0) = 1.0;
    W(0,0) = 0.0;
    V(1,0) = z * V(0,0);
    W(1,0) = 0.0;
    
    for (int n = 1; n < nmax; n++)
    {
        V(n+1,0) = ((2*(n+1)-1) * z * V(n,0) - (n) * V(n-1,0) ) / (n+1);
        W(n+1,0) = 0.0;
    }

    for (int m = 0; m < nmax ; m++)
    {
        V(m+1,m+1) = (2*(m +1)-1) * (x*V(m,m) - y*W(m,m));
        W(m+1,m+1) = (2*(m +1)-1) * (x*W(m,m) + y*V(m,m));

        if(m < (nmax - 1)) {
            
            V(m+2,m+1) = (2*(m+1)+1) * z * V(m+1,m+1);
            W(m+2,m+1) = (2*(m+1)+1) * z * W(m+1,m+1);
        }

        for (int n = m+1; n < nmax; n++)
        {
            V(n+1,m+1) = ((2*(n+1)-1)*z*V(n,m+1) - (n+m+1)*V(n-1,m+1)) / ((n+1)-(m+1));
            W(n+1,m+1) = ((2*(n+1)-1)*z*W(n,m+1) - (n+m+1)*W(n-1,m+1)) / ((n+1)-(m+1));
        }
        
    }


    this->V = V;
    this->W = W;
    this->lat = lat;
    this->lon = lon;
    


    // hysdrostatic

    double ahm = 0;
    double aha = 0;
    double awm = 0;
    double awa = 0;

    int i = 0;

    
    for (int n = 0; n < V.n_rows; n++) {
        
        for(int m = 0; m <= n; m++){
                
            ahm = ahm + (ah_mean(i)*V(n,m) + bh_mean(i)*W(n,m));
            aha = aha + (ah_amp(i) *V(n,m) + bh_amp(i) *W(n,m));
            awm = awm + (aw_mean(i)*V(n,m) + bw_mean(i)*W(n,m));
            awa = awa + (aw_amp(i) *V(n,m) + bw_amp(i) *W(n,m));
            i++;
        }
    }

        for (int n = 0; n < nmax + 1; n++)
        {   
            for (int m = 0; i <= n; m++)
            {
                ahm = ahm + (ah_mean(i)*V(n,m) + bh_mean(i)*W(n,m));
                aha = aha + (ah_amp(i) *V(n,m) + bh_amp(i) *W(n,m));
                awm = awm + (aw_mean(i)*V(n,m) + bw_mean(i)*W(n,m));
                awa = awa + (aw_amp(i) *V(n,m) + bw_amp(i) *W(n,m));

                i = i+1;
            }
            
        }
        
        this->ahm = ahm;
        this->aha = aha;
        this->awm = awm;
        this->awa = awa;
    }


    // Computation start

        double aw = (awm + awa*cos(doy*2*pi))*1e-5;
        double ah  = (ahm + aha*cos(doy*2*pi))*1e-5;

        double phi_h,c11_h,c10_h;
    
        if (lat<0) {      // southern hemisphere
                phi_h  = pi;
                c11_h = 0.007;
                c10_h = 0.002;
        
        } else {          // northern hemisphere
                phi_h  = 0;
                c11_h = 0.005;
                c10_h = 0.001;
        }

        double c0_h = 0.062;
        // hidrostatic b form Isobaric mapping function
        double bh = 0.002905;
        // c hydrostatic is taken from equation (7) in [1]
        double ch = c0_h + ((cos((doy - 28) / 365.25 * 2 * pi + phi_h) + 1) * c11_h / 2 + c10_h)*(1 - cos(lat));
        //wet b and c form Niell mapping function at 45ÿ lat tab 4 in [3]
        double bw = 0.00146;
        double cw = 0.04391;

        // Compute mapping functions
        // Ausiliary vectors
        vec ah_v(el.n_elem);
        ah_v.fill(ah);
        vec ch_v(el.n_elem);
        ch_v.fill(ch);
        vec aw_v(el.n_elem);
        aw_v.fill(aw);

        gmfh = mfContinuedFractionForm(ah_v,bh,ch_v,el);
        gmfw = mfContinuedFractionForm(aw_v,bw,cw,el);
        
        // correct hydrostatic for height
        vec ht_corr = hydrostaticMFHeigthCorrection(hgt,el);
    
        gmfh.replace(datum::nan,0);
        ht_corr.replace(datum::nan,0);

        gmfh = (gmfh + ht_corr);
        gmfh = (gmfh + ht_corr);
        
        //vec gmf = gmfh + gmfw;

        return make_tuple(gmfh,gmfw);
                
  }

  
vec Atmo::getIonoMF(double lat, double h, vec el, double rcm = NULL){
    if(rcm == NULL){

        rcm = getRcm();

    }

    double thin_shell_height = getThinShellHeight();

    el.transform([](double val) {if((val != datum::nan) && (val >= 0) ) {return (cos(val));}else {return val = datum::nan;}});
    vec iono_mf = (rcm + h) / (rcm + h + thin_shell_height) * el;
    iono_mf.for_each ( [] (double& val) {val = (pow(1 - (pow(val,2)),(-0.5))); });
    
    return iono_mf;
    
    }


    double Atmo::getThinShellHeight(){

        double thin_shell_height;

        if(this->ionex.height.empty()){
            thin_shell_height = 350 * 1e3; // if the ionex is not loaded use 350km
        } else {
            thin_shell_height = this->ionex.height[0] * 1e3;     // ionopshere thin shell height [km]
        }
        return thin_shell_height;
    }

    vec Atmo::mfContinuedFractionForm(vec a,double b,vec c,vec el)
    {   
        vec sine = sin(el);
        vec delay = (1 + (a / (1 + (b / (1 + c) )))) / (sine + (a / (sine + (b / (sine + c) ))));
        
        return delay;
    }
    vec Atmo::mfContinuedFractionForm(vec a,double b,double c, vec el)
    {   
        vec sine = sin(el);
        vec delay = (1 + (a / (1 + (b / (1 + c) )))) / (sine + (a / (sine + (b / (sine + c) ))));
        
        return delay;
    }
    vec Atmo::mfContinuedFractionForm(double a,double b,double c,vec el)
    {   
        vec sine = sin(el);
        vec delay = (1 + (a / (1 + (b / (1 + c) )))) / (sine + (a / (sine + (b / (sine + c) ))));
        
        return delay;
    }

    vec Atmo::hydrostaticMFHeigthCorrection(double h_ell, vec el)
    {
        // coorect the hysdrostatic mapping functions for the height
        // formulas and paramaters taken from :
        // Niell, A. E. "Global mapping functions for the atmosphere delay at radio wavelengths.
        // " Journal of Geophysical Research: Solid Earth 101.B2 (1996): 3227-3246.
        
        // height correction for the hydrostatic part

        // coefficent from tab 3
        double a_ht = 2.53e-5;
        double b_ht = 5.49e-3;
        double c_ht = 1.14e-3;
        double h_ell_km = h_ell/1000;   // convert height to km
        // eq (6)
        vec ht_corr_coef = 1 / sin(el) - mfContinuedFractionForm(a_ht,b_ht,c_ht,el);
        // eq (7)
        vec ht_corr = ht_corr_coef * h_ell_km;
    
        return ht_corr;

    }



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
void Atmo::gpt(double timeMJD, double dlat, double dlon, double dhgt) {
        
        vec pres,temp;
        double undu;

        double lat,lon,apm, apa, atm, ata;
        mat P;
        // Reference day is 28 January
        double doy = timeMJD  - 44239 + 1;
        //double doy = (gps_time / 86400 - 22) / 365.25;
        
        /* Values of a_geoid */
        vec a_geoid = { -5.6195e-001, -6.0794e-002, -2.0125e-001, -6.4180e-002, -3.6997e-002,
                        +1.0098e+001, +1.6436e+001, +1.4065e+001, +1.9881e+000, +6.4414e-001,
                        -4.7482e+000, -3.2290e+000, +5.0652e-001, +3.8279e-001, -2.6646e-002,
                        +1.7224e+000, -2.7970e-001, +6.8177e-001, -9.6658e-002, -1.5113e-002,
                        +2.9206e-003, -3.4621e+000, -3.8198e-001, +3.2306e-002, +6.9915e-003,
                        -2.3068e-003, -1.3548e-003, +4.7324e-006, +2.3527e+000, +1.2985e+000,
                        +2.1232e-001, +2.2571e-002, -3.7855e-003, +2.9449e-005, -1.6265e-004,
                        +1.1711e-007, +1.6732e+000, +1.9858e-001, +2.3975e-002, -9.0013e-004,
                        -2.2475e-003, -3.3095e-005, -1.2040e-005, +2.2010e-006, -1.0083e-006,
                        +8.6297e-001, +5.8231e-001, +2.0545e-002, -7.8110e-003, -1.4085e-004,
                        -8.8459e-006, +5.7256e-006, -1.5068e-006, +4.0095e-007, -2.4185e-008    };
        /* Values of b_geoid */
        vec b_geoid = { +0.0000e+000, +0.0000e+000, -6.5993e-002, +0.0000e+000, +6.5364e-002,
                        -5.8320e+000, +0.0000e+000, +1.6961e+000, -1.3557e+000, +1.2694e+000,
                        +0.0000e+000, -2.9310e+000, +9.4805e-001, -7.6243e-002, +4.1076e-002,
                        +0.0000e+000, -5.1808e-001, -3.4583e-001, -4.3632e-002, +2.2101e-003,
                        -1.0663e-002, +0.0000e+000, +1.0927e-001, -2.9463e-001, +1.4371e-003,
                        -1.1452e-002, -2.8156e-003, -3.5330e-004, +0.0000e+000, +4.4049e-001,
                        +5.5653e-002, -2.0396e-002, -1.7312e-003, +3.5805e-005, +7.2682e-005,
                        +2.2535e-006, +0.0000e+000, +1.9502e-002, +2.7919e-002, -8.1812e-003,
                        +4.4540e-004, +8.8663e-005, +5.5596e-005, +2.4826e-006, +1.0279e-006,
                        +0.0000e+000, +6.0529e-002, -3.5824e-002, -5.1367e-003, +3.0119e-005,
                        -2.9911e-005, +1.9844e-005, -1.2349e-006, -7.6756e-009, +5.0100e-008    };
        /* Values of ap_mean */ 
        vec ap_mean = { +1.0108e+003, +8.4886e+000, +1.4799e+000, -1.3897e+001, +3.7516e-003,
                        -1.4936e-001, +1.2232e+001, -7.6615e-001, -6.7699e-002, +8.1002e-003,
                        -1.5874e+001, +3.6614e-001, -6.7807e-002, -3.6309e-003, +5.9966e-004,
                        +4.8163e+000, -3.7363e-001, -7.2071e-002, +1.9998e-003, -6.2385e-004,
                        -3.7916e-004, +4.7609e+000, -3.9534e-001, +8.6667e-003, +1.1569e-002,
                        +1.1441e-003, -1.4193e-004, -8.5723e-005, +6.5008e-001, -5.0889e-001,
                        -1.5754e-002, -2.8305e-003, +5.7458e-004, +3.2577e-005, -9.6052e-006,
                        -2.7974e-006, +1.3530e+000, -2.7271e-001, -3.0276e-004, +3.6286e-003,
                        -2.0398e-004, +1.5846e-005, -7.7787e-006, +1.1210e-006, +9.9020e-008,
                        +5.5046e-001, -2.7312e-001, +3.2532e-003, -2.4277e-003, +1.1596e-004,
                        +2.6421e-007, -1.3263e-006, +2.7322e-007, +1.4058e-007, +4.9414e-009    };
        /* Values of bp_mean */
        vec bp_mean = { +0.0000e+000, +0.0000e+000, -1.2878e+000, +0.0000e+000, +7.0444e-001,
                        +3.3222e-001, +0.0000e+000, -2.9636e-001, +7.2248e-003, +7.9655e-003,
                        +0.0000e+000, +1.0854e+000, +1.1145e-002, -3.6513e-002, +3.1527e-003,
                        +0.0000e+000, -4.8434e-001, +5.2023e-002, -1.3091e-002, +1.8515e-003,
                        +1.5422e-004, +0.0000e+000, +6.8298e-001, +2.5261e-003, -9.9703e-004,
                        -1.0829e-003, +1.7688e-004, -3.1418e-005, +0.0000e+000, -3.7018e-001,
                        +4.3234e-002, +7.2559e-003, +3.1516e-004, +2.0024e-005, -8.0581e-006,
                        -2.3653e-006, +0.0000e+000, +1.0298e-001, -1.5086e-002, +5.6186e-003,
                        +3.2613e-005, +4.0567e-005, -1.3925e-006, -3.6219e-007, -2.0176e-008,
                        +0.0000e+000, -1.8364e-001, +1.8508e-002, +7.5016e-004, -9.6139e-005,
                        -3.1995e-006, +1.3868e-007, -1.9486e-007, +3.0165e-010, -6.4376e-010    };
        /* Values of ap_amp */
        vec ap_amp = {  -1.0444e-001, +1.6618e-001, -6.3974e-002, +1.0922e+000, +5.7472e-001,
                        -3.0277e-001, -3.5087e+000, +7.1264e-003, -1.4030e-001, +3.7050e-002,
                        +4.0208e-001, -3.0431e-001, -1.3292e-001, +4.6746e-003, -1.5902e-004,
                        +2.8624e+000, -3.9315e-001, -6.4371e-002, +1.6444e-002, -2.3403e-003,
                        +4.2127e-005, +1.9945e+000, -6.0907e-001, -3.5386e-002, -1.0910e-003,
                        -1.2799e-004, +4.0970e-005, +2.2131e-005, -5.3292e-001, -2.9765e-001,
                        -3.2877e-002, +1.7691e-003, +5.9692e-005, +3.1725e-005, +2.0741e-005,
                        -3.7622e-007, +2.6372e+000, -3.1165e-001, +1.6439e-002, +2.1633e-004,
                        +1.7485e-004, +2.1587e-005, +6.1064e-006, -1.3755e-008, -7.8748e-008,
                        -5.9152e-001, -1.7676e-001, +8.1807e-003, +1.0445e-003, +2.3432e-004,
                        +9.3421e-006, +2.8104e-006, -1.5788e-007, -3.0648e-008, +2.6421e-010    };
        
        /* Values of bp_amp */
        vec bp_amp = {  +0.0000e+000, +0.0000e+000, +9.3340e-001, +0.0000e+000, +8.2346e-001,
                        +2.2082e-001, +0.0000e+000, +9.6177e-001, -1.5650e-002, +1.2708e-003,
                        +0.0000e+000, -3.9913e-001, +2.8020e-002, +2.8334e-002, +8.5980e-004,
                        +0.0000e+000, +3.0545e-001, -2.1691e-002, +6.4067e-004, -3.6528e-005,
                        -1.1166e-004, +0.0000e+000, -7.6974e-002, -1.8986e-002, +5.6896e-003,
                        -2.4159e-004, -2.3033e-004, -9.6783e-006, +0.0000e+000, -1.0218e-001,
                        -1.3916e-002, -4.1025e-003, -5.1340e-005, -7.0114e-005, -3.3152e-007,
                        +1.6901e-006, +0.0000e+000, -1.2422e-002, +2.5072e-003, +1.1205e-003,
                        -1.3034e-004, -2.3971e-005, -2.6622e-006, +5.7852e-007, +4.5847e-008,
                        +0.0000e+000, +4.4777e-002, -3.0421e-003, +2.6062e-005, -7.2421e-005,
                        +1.9119e-006, +3.9236e-007, +2.2390e-007, +2.9765e-009, -4.6452e-009    };
        
        /* Values of at_mean */
        vec at_mean = { +1.6257e+001, +2.1224e+000, +9.2569e-001, -2.5974e+001, +1.4510e+000,
                        +9.2468e-002, -5.3192e-001, +2.1094e-001, -6.9210e-002, -3.4060e-002,
                        -4.6569e+000, +2.6385e-001, -3.6093e-002, +1.0198e-002, -1.8783e-003,
                        +7.4983e-001, +1.1741e-001, +3.9940e-002, +5.1348e-003, +5.9111e-003,
                        +8.6133e-006, +6.3057e-001, +1.5203e-001, +3.9702e-002, +4.6334e-003,
                        +2.4406e-004, +1.5189e-004, +1.9581e-007, +5.4414e-001, +3.5722e-001,
                        +5.2763e-002, +4.1147e-003, -2.7239e-004, -5.9957e-005, +1.6394e-006,
                        -7.3045e-007, -2.9394e+000, +5.5579e-002, +1.8852e-002, +3.4272e-003,
                        -2.3193e-005, -2.9349e-005, +3.6397e-007, +2.0490e-006, -6.4719e-008,
                        -5.2225e-001, +2.0799e-001, +1.3477e-003, +3.1613e-004, -2.2285e-004,
                        -1.8137e-005, -1.5177e-007, +6.1343e-007, +7.8566e-008, +1.0749e-009    };

        /* Values of bt_mean */
        vec bt_mean = { +0.0000e+000, +0.0000e+000, +1.0210e+000, +0.0000e+000, +6.0194e-001,
                        +1.2292e-001, +0.0000e+000, -4.2184e-001, +1.8230e-001, +4.2329e-002,
                        +0.0000e+000, +9.3312e-002, +9.5346e-002, -1.9724e-003, +5.8776e-003,
                        +0.0000e+000, -2.0940e-001, +3.4199e-002, -5.7672e-003, -2.1590e-003,
                        +5.6815e-004, +0.0000e+000, +2.2858e-001, +1.2283e-002, -9.3679e-003,
                        -1.4233e-003, -1.5962e-004, +4.0160e-005, +0.0000e+000, +3.6353e-002,
                        -9.4263e-004, -3.6762e-003, +5.8608e-005, -2.6391e-005, +3.2095e-006,
                        -1.1605e-006, +0.0000e+000, +1.6306e-001, +1.3293e-002, -1.1395e-003,
                        +5.1097e-005, +3.3977e-005, +7.6449e-006, -1.7602e-007, -7.6558e-008,
                        +0.0000e+000, -4.5415e-002, -1.8027e-002, +3.6561e-004, -1.1274e-004,
                        +1.3047e-005, +2.0001e-006, -1.5152e-007, -2.7807e-008, +7.7491e-009    };

        /* Values of at_amp */
        vec at_amp = {  -1.8654e+000, -9.0041e+000, -1.2974e-001, -3.6053e+000, +2.0284e-002,
                        +2.1872e-001, -1.3015e+000, +4.0355e-001, +2.2216e-001, -4.0605e-003,
                        +1.9623e+000, +4.2887e-001, +2.1437e-001, -1.0061e-002, -1.1368e-003,
                        -6.9235e-002, +5.6758e-001, +1.1917e-001, -7.0765e-003, +3.0017e-004,
                        +3.0601e-004, +1.6559e+000, +2.0722e-001, +6.0013e-002, +1.7023e-004,
                        -9.2424e-004, +1.1269e-005, -6.9911e-006, -2.0886e+000, -6.7879e-002,
                        -8.5922e-004, -1.6087e-003, -4.5549e-005, +3.3178e-005, -6.1715e-006,
                        -1.4446e-006, -3.7210e-001, +1.5775e-001, -1.7827e-003, -4.4396e-004,
                        +2.2844e-004, -1.1215e-005, -2.1120e-006, -9.6421e-007, -1.4170e-008,
                        +7.8720e-001, -4.4238e-002, -1.5120e-003, -9.4119e-004, +4.0645e-006,
                        -4.9253e-006, -1.8656e-006, -4.0736e-007, -4.9594e-008, +1.6134e-009    };

        /* Values of bt_amp */ 
        vec bt_amp = {  +0.0000e+000, +0.0000e+000, -8.9895e-001, +0.0000e+000, -1.0790e+000,
                        -1.2699e-001, +0.0000e+000, -5.9033e-001, +3.4865e-002, -3.2614e-002,
                        +0.0000e+000, -2.4310e-002, +1.5607e-002, -2.9833e-002, -5.9048e-003,
                        +0.0000e+000, +2.8383e-001, +4.0509e-002, -1.8834e-002, -1.2654e-003,
                        -1.3794e-004, +0.0000e+000, +1.3306e-001, +3.4960e-002, -3.6799e-003,
                        -3.5626e-004, +1.4814e-004, +3.7932e-006, +0.0000e+000, +2.0801e-001,
                        +6.5640e-003, -3.4893e-003, -2.7395e-004, +7.4296e-005, -7.9927e-006,
                        -1.0277e-006, +0.0000e+000, +3.6515e-002, -7.4319e-003, -6.2873e-004,
                        -8.2461e-005, +3.1095e-005, -5.3860e-007, -1.2055e-007, -1.1517e-007,
                        +0.0000e+000, +3.1404e-002, +1.5580e-002, -1.1428e-003, +3.3529e-005,
                        +1.0387e-005, -1.9378e-006, -2.7327e-007, +7.5833e-009, -9.2323e-009    };
        
        // Computing Legendre Polynomial
        double t = sin(dlat);

        int nmax = 9;

        vec dfac(2 * nmax + 1);
        dfac(0) = 1;
        for (int i = 1; i <= 2 * nmax; i++)
            dfac(i) = dfac(i - 1) * i;

        // Determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62)
        mat P(nmax + 1, nmax + 1);
        for (int i = 0; i <= nmax; i++) {
            for (int j = 0; j <= min(i, nmax); j++) {
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

        this->P = P;
        this->lat = dlat;
        this->lon = dlon;

        // Spherical harmonics
        vec aP, bP;
        // index
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

        double apm = sum(ap_mean) * aP + sum(bp_mean) * bP;
        double apa = sum(ap_amp) * aP + sum(bp_amp) * bP;
        double atm = sum(at_mean) * aP + sum(bt_mean) * bP;
        double ata = sum(at_amp) * aP + sum(bt_amp) * bP;

        this->apm = apm;
        this->apa = apa;
        this->atm = atm; 
        this->ata = ata;

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
};






vec Atmo::saastamoinenModelGPT(double mJD, double lat, double lon, double h, double undu, vec el) {
    // Compute pressure and temperature using GPT model
    arma::vec pres, temp;
    gpt(mJD, lat * arma::datum::pi / 180, lon * arma::datum::pi / 180, h);

    // Compute tropospheric delays
    double t_h = h;

    if (undu > -300) {
        t_h -= undu;
    }

    
    vec ZHD_R = saast_dry(pres, t_h, lat);
    vec ZWD_R = saast_wet(temp, STD_HUMI, t_h);

    // Apply GMF corrections
    vec gmfh_R, gmfw_R;
    tie(gmfh_R,gmfw_R) = gmf(mJD,lat * arma::datum::pi / 180, lon * arma::datum::pi / 180, h - undu, el * arma::datum::pi / 180);

    // Compute delay
    return gmfh_R * ZHD_R + gmfw_R * ZWD_R;
}


vec Atmo::saastamoinenModelPTH(double mJD, double lat, double lon, double h, arma::vec undu, arma::vec el, vec P, vec T, vec H) {
    // Check if meteorological data is valid
    if (any(vec({P, T, H}).has_nan())) {
        // Compute pressure and temperature using GPT model
        arma::vec pres, temp;
        gpt(mJD, lat * arma::datum::pi / 180, lon * arma::datum::pi / 180, h);

        // Set missing values
        if (arma::isnan(P)) P = pres;
        if (isnan(T)) T = temp;
        if (isnan(H)) H.fill(this->STD_HUMI);

        // Log warning
        //this->log.addWarning("No valid meteo data are present @ " + std::to_string(gps_time / 86400 + GPS_Time.GPS_ZERO) + "\nUsing standard GPT values\n - " + std::to_string(T) + " °C\n - " + std::to_string(P) + " hPa\n - humidity " + std::to_string(H), 100);
    }

    // Compute tropospheric delays
    arma::vec ZHD_R = saast_dry(P, h, lat);
    arma::vec ZWD_R = saast_wet(T, H, h);

    // Apply GMF corrections
    arma::vec gmfh_R, gmfw_R;
    tie(gmfh_R,gmfw_R) = this.gmf(mJD,lat * arma::datum::pi / 180, lon * arma::datum::pi / 180, h - undu, el * arma::datum::pi / 180);

    // Compute delay
    return gmfh_R * ZHD_R + gmfw_R * ZWD_R;
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
arma::vec saast_wet(vec T, vec H, double h) {
    // Convert temperature from Celsius to Kelvin
    T += 273.15;
    
    // Perform height correction
    // (Comment out the following line if height correction is done before)
    // H %= arma::exp(-0.0006396 * h);
    
    // Convert humidity from percentage to fraction
    H /= 100.0;
    
    // Compute constant and exponent
    arma::vec c = -37.2465 + 0.213166 * T - 2.56908 * pow(10, -4) * arma::pow(T, 2);
    arma::vec e = H * arma::exp(c);

    // Compute Zenith Wet Delay (ZWD) using Saastamoinen model
    return 0.0022768 * (((1255.0 / T) + 0.05) * e);
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
arma::vec saast_dry(vec P, double h, double lat) {
    double t = cosd((2 * lat));
    // Compute Zenith Hydrostatic Delay (ZHD) using Saastamoinen model
    return 0.0022768 * P * (1 + 0.00266 * t + 0.00000028 * h);
    // Alternatively, for alternative formula:
    // return 0.0022767 * P / (1 - 0.00266 * arma::cosd(2 * lat) - 0.00000028 * h);
}

