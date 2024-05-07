#include <sstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <armadillo>


#include "coresky.h"
#include "utility.h"

using namespace std;
using namespace arma;


// Constructor
CoreSky::CoreSky()
{
    pathFilename = "nofilename.out";
    fileMode = ios_base::out;
}

CoreSky::CoreSky(string inputFilePath, vector<char> inputSS ,ios_base::openmode mode)
{
    string record;
    pathFilename = inputFilePath;
    fileMode = mode;
    fileStream.open( pathFilename.c_str(), fileMode );
    if( ! fileStream ) {
      cerr << "Error: unable to open SP3 file in constructor: " << pathFilename << " using mode: " << fileMode << endl;
    }
    if(!inputSS.empty()){
      for (size_t i = 0; i < inputSS.size(); i++)
      { 
        if (find(defSatelliteSystems.begin(), defSatelliteSystems.end(), inputSS[i]) != defSatelliteSystems.end())
        {
          satelliteSystems.push_back(inputSS[i]);
        }
      }
      
    } else {
      satelliteSystems = defSatelliteSystems;
    }
}

// Destructor
CoreSky::~CoreSky()
{
}

// Setters

void CoreSky::setPathFilenameMode(string path, ios_base :: openmode mode)
{

    pathFilename = path;
    fileMode = mode;
    fileStream.open( pathFilename.c_str(), fileMode);
    if(fileStream.bad()){
        cerr << "Error: unable to open SP3 file in constructor: "
	     << pathFilename << " using mode: " << fileMode << endl;
//  Add error handler in opening files!!!
    }

}

double CoreSky::getStartTime()
{
  return this->SP3StartTime;
}

mat CoreSky::getClock_corr()
{
  return this->clock;
}

double CoreSky::getEndTime()
{
  return this->SP3EndTime;
}

int CoreSky::getNumOfsat()
{
  return this->satIds.size();
}

void CoreSky::printSatIds()
{

  vector<string> tmp = this->satIds;
  cout << "---------------------------------" << endl;
  cout << "|   Satellite ids and Indexes   |" << endl;
  cout << "|                               |" << endl;
  for (size_t i = 0; i < tmp.size(); i++)
  {
    cout << "| " << tmp[i] << " : " << to_string(i) << "\t\t\t|" << endl;
  }
  cout << "---------------------------------" << endl;
  return;
}

vector<string> CoreSky::getSatIds()
{
  return this->satIds;
}

vec CoreSky::getEpoch()
{
  return ref_epoch;
}

void CoreSky::setSS(vector<char> inputSS){

  satelliteSystems.clear();

  if(!inputSS.empty()){
      for (size_t i = 0; i < inputSS.size(); i++)
      { 
        if (find(defSatelliteSystems.begin(), defSatelliteSystems.end(), inputSS[i]) != defSatelliteSystems.end())
        {
          satelliteSystems.push_back(inputSS[i]);
        }
      }
      
  } else {
      cout << "Empty list of satellite systems, default list will be uploaded instead" << endl;
      satelliteSystems = defSatelliteSystems;
  }  
}

// Organizes Epoch Time Information into Vector
vector<double> CoreSky::SP3EpochTimeOrganizer(string line) {
	vector<double> epochRecord;
	// Splitting words in the line
	istringstream iss(line);
	vector<string> words{ istream_iterator<string>{iss}, istream_iterator<string>{} };
	for (string s : words) {
		epochRecord.push_back(stod(s));
	}
	// Return: Vector containing the Epoch Info
	return epochRecord;
}

void CoreSky::initHeader(){
   formatVersion = ' ';
   modeFlag = ' ';
   SP3StartTime = 0;
   SP3EndTime = 0;
   numberSP3epochs = 0;
   dataUsed = "";
   coordFrame = "";
   orbitType = "";
   sourceAgency = "";
   gpsWeek = 9999;
   secsOfWeek = 0.0;
   SP3interval = 0.0;
   SP3mjd = 999999;
   SP3fmjd = 0.0;
   numberSP3svs = 0;
}

void CoreSky::clearCoord(){
  coord.clear();
}

void CoreSky::cleaSetPoly(){
  coef_poly.clear();
}

void CoreSky::clearVariables(){
  this->initHeader();
  this->clearCoord();
  this->cleaSetPoly(); 
  this->isempty = 1;
}

int CoreSky::readHeader() {

  cout << "---------------------------------" << endl;
  cout << "     *** readHeader: ***" << endl << endl;
    
  string      recType; // usede to read first two digit of the row
  string      inputRec; // used to manage 
  string      temp;
  vector<double> date;
  char        nextFirstChar; // check next line first char
  char        inputRecC[256];

  //values variables
  int         indexSatId = 0;
  int         indexACC = 0;
  int         lineLength = 0;
  int         i,accur;
  
  
  // START READING HEADER
  if(isempty) {
    // Initialize header variables
    initHeader();
    // Change coresky status
    isempty = !isempty;
    // position the stream at the beginning of the file
    fileStream.seekg(0);
    // Start reading lines
    while( fileStream.getline(inputRecC, 255, '\n') ) {
      inputRec = inputRecC;
      recType = inputRec.substr( 0, 2 );
      // Check type of file
      if( recType[0] == '#' && recType[1] != '#' ) {
	      formatVersion = recType[1];
       
	      if( formatVersion != 'd' ) {
	        cerr << "This program reads only SP3 files with format-mode: aP "
	             << endl << "This file has format version: " << formatVersion << endl;
	        return( -1 );
        }

	    modeFlag = inputRec[2];
       
	    if( modeFlag != 'P') { 
        cerr << "This program reads only SP3 files with format-mode: aP "
	           << endl << "This file has mode flag: " << modeFlag << endl;
	      return( -1 );
      }

      // Save SP3 start time  
      temp = inputRec.substr(3, 30);
      date = SP3EpochTimeOrganizer(temp);
      SP3StartTime = gpsTime(date);
      cout << "Date: * " << temp << " *" << endl;

    
      // Save number of epoch in file
      temp = inputRec.substr(32, 7 );
      numberSP3epochs = (unsigned long) atol( temp.c_str() );
      // More info  
      dataUsed = inputRec.substr(40, 5 );
      coordFrame = inputRec.substr(46, 5 );
      orbitType = inputRec.substr(52, 3 );
      sourceAgency = inputRec.substr(56, 4 );
      
      }
      // Second line 
      else if( recType[0] == '#' && recType[1] == '#' ) {
        
        temp = inputRec.substr(3, 4 );
        gpsWeek = (unsigned long) atol( temp.c_str() );
        temp = inputRec.substr(8, 15 );
        secsOfWeek = atof( temp.c_str() );
        temp = inputRec.substr(24, 14 );
        SP3interval = atof( temp.c_str() );
        temp = inputRec.substr(39, 5 );
        SP3mjd = atol( temp.c_str() );
        temp = inputRec.substr(45, 15 );
        SP3fmjd = atof( temp.c_str() );
        
        if(SP3interval == 0){
          // default interval if not specified
          SP3interval = 10;
          numberSP3epochs = 86400/SP3interval;
        }
      }
      // Third line: list of sat ids in the file
      else if( recType[0] == '+' && recType[1] == ' ' )
      {
        // check if I already read first line
        if( numberSP3svs == 0 ) {
          // number of sat in file
          temp = inputRec.substr( 3, 3 );
          numberSP3svs = (unsigned short) atoi( temp.c_str() );
          // lenght of line
          lineLength = inputRec.size();
    
          for( i = 9; i < lineLength; i=i+3 ) {
            
            temp = inputRec.substr( i, 3 );
            if(find(satelliteSystems.begin(), satelliteSystems.end(), temp.front()) != satelliteSystems.end()){
              satIds.push_back(temp);
              int svn = atoi( temp.substr(1,2).c_str() );
              if( svn > 0 ) {
                sp3PRN.push_back(indexSatId);
                indexSatId++;
              }
            }
            
          }

        } else {
          //once first line has been read, read other not savings zeros 
          lineLength = inputRec.size();
          
          for( i = 9; i < lineLength; i=i+3 ) {
            temp = inputRec.substr( i, 3 );
            if(find(satelliteSystems.begin(), satelliteSystems.end(), temp.front()) != satelliteSystems.end()){
              int svn = atoi( temp.substr(1,2).c_str() );
            
              if( svn > 0 ) {
                satIds.push_back(temp);
                sp3PRN.push_back(indexSatId);
                indexSatId++;
              }
            }
          }
        }
      }

      // Check on accuracy data
      else if( recType[0] == '+' && recType[1] == '+' ) {
        
        lineLength = inputRec.size();
	      for( i = 9; i < lineLength; i=i+3 ) {
          temp = inputRec.substr( i, 3 );

          if(find(sp3PRN.begin(), sp3PRN.end(),indexACC ) != sp3PRN.end()){

	          accur = atoi( temp.c_str() );
	          if( accur > 0 ) {
              svAccu.push_back((unsigned short) accur);
	            if( accur <= 0 )
	            {
                cerr << "WARNING ! accuracy code ZERO for PRN: " <<
                satIds[indexACC] << endl;
                sp3PRN[indexACC] = (-1.0 * sp3PRN[indexACC]);

              }
	            indexACC++;
	          }
          }
        }
      }
      nextFirstChar = (char) fileStream.peek();
      if( nextFirstChar == '*' ) break;  // exit while loop
      
    } // end of while loop to read all SP3 header records
    
    numberSVparams = 4;
    numberGoodPRNs = (unsigned short) indexSatId;
    numberGoodACCURs = (unsigned short) indexACC;
  

    if( numberGoodPRNs <= 0  ||  numberGoodACCURs != numberGoodPRNs ) {

      cerr << "Error reading the PRNs from the header for SP3 file: "
      << endl << pathFilename << endl
      << "Number of SVs expected: " << setw(6) << numberSP3svs << endl
      << "Number of PRNs read in: " << setw(6) << numberGoodPRNs << endl;
      return(1);

    } else {
      
    // class saved elements
    cout << "---------------------------------" << endl;
    cout << "gpsWeek: " << gpsWeek << endl;
    cout << "secOfWeek: " << secsOfWeek << endl;
    cout << "SP3Interval: " << SP3interval << endl;
    cout << "SP3mjd: " << SP3mjd << endl;
    cout << "formatVersion: " << formatVersion << endl;
    cout << "Data used: " << dataUsed << endl;
    cout << "modeFlag: " << modeFlag << endl;
    cout << "sourceAgency: " << sourceAgency << endl;
    cout << "numberSP3svs: " << numberSP3svs << endl;
    cout << "numberSP3epoch: " << numberSP3epochs << endl << endl; 
    cout << "SP3 start time: " << SP3StartTime << endl;
    cout << "---------------------------------" << endl;
    cout << "End of reading header..." << endl;

    printSatIds();

    std::sort(satIds.begin(), satIds.end());

    printSatIds();

    
    return (0);
   }

  } else {

    cerr << "There is already a file loaded!" << endl;
    return(1);
    //metodo per controllare che header sia uguale ?
    // non più necessario sola lettura file singolo
  }
  
} // end of method SP3File::readHeader()


void CoreSky::readCoord(){

    cout << "---------------------------------" << endl;
    cout << "*** ReadCord ***" << endl;

    string      recType; // usede to read first two digit of the row
    string      inputRec; // used to manage readed line 
    string      temp;
    char        nextFirstChar; // check next line first char
    char        inputRecC[ 256 ];
    
    vector<double> date; //vector to handle date to conversion
    int epoch;
    int sat,cont,svn;
    long double num;
    
    coord = cube(numberSP3epochs,numberSP3svs,3);
    clock = mat(numberSP3epochs,numberSP3svs);
    ref_epoch = vec(numberSP3epochs);
    
 
    cont = 0;
    while( fileStream.getline(inputRecC, 255, '\n') )
    {
        
        inputRec = inputRecC;
        recType = inputRec.substr( 0, 2 );
        
        if (recType[0] == '*')
        {   // Reading time of observations

            temp = inputRec.substr(3,20);
            date = SP3EpochTimeOrganizer(temp);
            double epochsp3 = gpsTime(date);

            // dealing with week change in gps time
            if(epochsp3 == 0 && cont != 0){
              epochsp3 = SP3interval + epochs[cont-1];
            }

            epochs.push_back(epochsp3);
            epoch = (epochsp3 - SP3StartTime)/SP3interval;
            
            ref_epoch(cont) = epochsp3 - SP3StartTime;
            cont++;

            // Check if last time epoch
            if(cont == numberSP3epochs){
              // Save SP3 end time 
              SP3EndTime = epochsp3;
              cout << "EndTime: " << SP3EndTime << endl;
              
            }

            

        } else if(recType[0]== modeFlag){

            temp = inputRec.substr(1,3);
            sat = std::distance(satIds.begin(),find(satIds.begin(),satIds.end(),temp));

            if(sat >= numberSP3svs){
              cerr << "Error! There is an extra line in file not corresponding to the header";
              cout << endl << "Check following epoch * " << temp << " in the file " << endl;
              break;
            } 

            temp = inputRec.substr(5, 13);
            num = stold(temp.c_str())*1000;
            coord(epoch,sat,0) = num;

            temp = inputRec.substr(19, 13);
            num = stold(temp.c_str())*1000;
            coord(epoch,sat,1) = num;

            temp = inputRec.substr(33, 13);
            num = stold(temp.c_str())*1000;            
            coord(epoch,sat,2) = num;

            temp = inputRec.substr(47, 13);
            num = stold(temp.c_str());            
            clock(epoch,sat) = num;

            
        }

   }
   
  cout << "NumberSP3svs : " << numberSP3svs << endl;
  cout << "NumberSP3epochs : " << numberSP3epochs << endl;
  cout << "Time passed in seconds : " << (numberSP3epochs - 1)*300 << endl;
  
  cout << SP3StartTime << endl;
  cout << SP3EndTime << endl;
  cout << clock(0,0) << endl;
 
   return;
}


void CoreSky::computePolyCoef(){

  cout << "---------------------------------" << endl;
  cout << "*** compute polyCoef: ***" << endl;
  
  int obs = 10;
  int order = 11;
  int n_coef_set = numberSP3epochs - obs;  
  int n_sat = numberSP3svs;
  
  // set of coef 
  cube set_coef = cube(order,3,n_sat);
  coef_poly.reserve(n_coef_set);
  
  mat A = ones<mat>(order, order);
  
  vec x {-5,-4,-3,-2,-1,0,1,2,3,4,5};

  for (int i = 0; i < order; i++){
    A.col(i) = pow(x,i); 
  }

  mat A_inv = inv(A);
    
  for (int set = 0; set < n_coef_set; set++){
    
    for (int sat = 0; sat < n_sat; sat++){
      
      mat coord_set = coord.subcube(set, sat, 0, set + obs, sat, 2);
      // Compute the polynomial coefficients using the least squares method
      mat coef_set_sat = A_inv * coord_set;
      // Store the polynomial coefficients in the vector for the current satellite and epoch set
      set_coef.slice(sat) = coef_set_sat;           
    }
    
    coef_poly.push_back(set_coef);      
  } 
// Printing information about coordinates set
cout << "Number of sets : " << n_coef_set << endl;
cout << "Size of coef poly : " << coef_poly.size() << endl;
cout << "Size of mat in coef : " << coef_poly[0].size() << endl;

return;
}
void CoreSky::computeClockPolyCoef(){

  cout << "---------------------------------" << endl;
  cout << " *** compute Clock polyCoef: ***" << endl;
  
  int obs = 10;
  int order = 11;
  int n_coef_set = numberSP3epochs - obs;  
  int n_sat = numberSP3svs;
  
  // set of coef 
  mat set_coef = mat(order,n_sat);
  coef_poly_clock.reserve(n_coef_set);
  
  mat A = ones<mat>(order, order);
  
  vec x {-5,-4,-3,-2,-1,0,1,2,3,4,5};

  for (int i = 0; i < order; i++){
    A.col(i) = pow(x,i); 
  }

  mat A_inv = inv(A);
    
  for (int set = 0; set < n_coef_set; set++){
    
    for (int sat = 0; sat < n_sat; sat++){
      // Vector of clock parameters for a given set
      vec clock_set = clock.submat(set,sat,set+obs,sat);
      // Compute the polynomial coefficients using the least squares method
      vec coef_set_clock = A_inv * clock_set;
      // Store the polynomial coefficients in the vector for the current satellite and epoch set
      set_coef.col(sat) = coef_set_clock;          
    }
    
    coef_poly_clock.push_back(set_coef);      
  } 
// Printing file information about code
cout << "Number of sets : " << n_coef_set << endl;
cout << "Size of coef poly : " << coef_poly_clock.size() << endl;
cout << "Size of mat in coef : " << coef_poly_clock[0].size() << endl;
cout << coef_poly_clock[0].at(0) << endl;
return;
}


/*cube CoreSky::coordInterpolate(vec time, vec sat){

  cout << "---------------------------------" << endl;
  cout << "*** Coordinate Interpolation *** " << endl;
  cout << SP3EndTime -SP3StartTime << endl;

  //timer.tic();
  int n_sat = sat.size();
  int n_epoch = time.size();
  int n_border = 5;

  //n_border = ((size(this.coord, 1) - size(poly, 4)) / 2);

  //pid_ceil = ceil(t_diff / this.coord_rate) + 1 - n_border;
  // Ignore solution at the border of the polynomial
  // pid_ceil(pid_ceil < 1) = 1;
  // pid_ceil(pid_ceil > size(this.getPolyCoeff, 4)) = size(this.getPolyCoeff, 4);

  vec c_time = time - SP3StartTime;
  
  cube X(n_epoch,n_sat,3);
  mat clock_corr(n_epoch,n_sat);

  vec eval(11);
  vec t_fct_powers(11); // Cache powers of t_fct

  for (int i = 0; i < n_epoch; i++){ 
      
    int set_id = round((c_time(i)/SP3interval) + 1 - ((ref_epoch.size() - coef_poly.size()) / 2 ));
    int ceil_id = ceil((c_time(i)/SP3interval) + 1 - ((ref_epoch.size() - coef_poly.size()) / 2 ));
    
    if(set_id < 0) set_id = 0;
    if(set_id >= coef_poly.size()) set_id = coef_poly.size() - 1;
    
    if(ceil_id < 0) ceil_id = 0;
    if(ceil_id >= coef_poly.size()) ceil_id = coef_poly.size() - 1;

    double t_fct = (c_time(i) - ref_epoch(set_id + ((ref_epoch.size() - coef_poly.size())/ 2)))/SP3interval;
   
    t_fct_powers(0) = 1;
      
    for (int j = 1; j < 11; j++) {
          
        t_fct_powers(j) = t_fct_powers(j - 1) * t_fct;
    }

    eval = t_fct_powers;
      
    for (int j = 0; j < n_sat; j++) {

      mat coef_set_sat = coef_poly[set_id].slice(sat(j));
      vec coef_set_clock = coef_poly_clock[set_id].col(sat(j));
      X(i, sat(j), 0) = dot(eval, coef_set_sat.col(0));
      X(i, sat(j), 1) = dot(eval, coef_set_sat.col(1));
      X(i, sat(j), 2) = dot(eval, coef_set_sat.col(2));
      clock_corr(i,sat(j)) = dot(eval,coef_set_clock);
    }
  }
  cout << "clock " << endl << clock_corr(0,0) << endl;
  cout << clock(0,0) << endl;
  cout << "---------------------------------" << endl;
  cout << "***** End of coordinates computations *****" << endl;
  return X;


}*/

cube CoreSky::coordInterpolate(vec time, vec sat){

  cout << "---------------------------------" << endl;
  cout << "*** Coordinate Interpolation *** " << endl;
  
  int n_sat = sat.size();
  int n_epoch = time.size();
  int n_border = 5;

  vec c_time = time - SP3StartTime;
  
  cube X(n_epoch,n_sat,4);

  vec eval(11);
  vec t_fct_powers(11); // Cache powers of t_fct

  for (int i = 0; i < n_epoch; i++){ 
      
    int set_id = round((c_time(i)/SP3interval) + 1 - ((ref_epoch.size() - coef_poly.size()) / 2 ));
    int ceil_id = ceil((c_time(i)/SP3interval) + 1 - ((ref_epoch.size() - coef_poly.size()) / 2 ));
    
    if(set_id < 0) set_id = 0;
    if(set_id >= coef_poly.size()) set_id = coef_poly.size() - 1;
    
    if(ceil_id < 0) ceil_id = 0;
    if(ceil_id >= coef_poly.size()) ceil_id = coef_poly.size() - 1;

    double t_fct = (c_time(i) - ref_epoch(set_id + ((ref_epoch.size() - coef_poly.size())/ 2)))/SP3interval;
   
    t_fct_powers(0) = 1;
      
    for (int j = 1; j < 11; j++) {
          
        t_fct_powers(j) = t_fct_powers(j - 1) * t_fct;
    }

    eval = t_fct_powers;
      
    for (int j = 0; j < n_sat; j++) {

      mat coef_set_sat = coef_poly[set_id].slice(sat(j));
      vec coef_set_clock = coef_poly_clock[set_id].col(sat(j));
      
      X(i, sat(j), 0) = arma::dot(eval, coef_set_sat.col(0));
      X(i, sat(j), 1) = arma::dot(eval, coef_set_sat.col(1));
      X(i, sat(j), 2) = arma::dot(eval, coef_set_sat.col(2));
      X(i, sat(j), 3) = arma::dot(eval,coef_set_clock);

    }
  }
  cout << "---------------------------------" << endl;
  cout << "***** End of coordinates computations *****" << endl;
  return X;
}



   

