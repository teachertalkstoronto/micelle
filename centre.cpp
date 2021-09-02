#include <iostream>     //For output to shell
#include <string>
#include <iomanip>
#include "vecmat3.h"    //For vector/matrix operations
#include <cmath>        //For mathematical operations
#include <list>
#include "constants_1.h"	//Holds constant valued paramters
#include "classes.h"    //Holds different classes used                        
#include "functions.h"  //Functions
#include <fstream>		//Output binary data to binary file, and psf and xsf data


using namespace constants;		//Holds constant variables

int main() {
  //1. VARIABLES HEADLIST, CENTRE, AND SPEED.___________________________________________________
	std::list<Vector> headlist;  
    Vector centre(0.0,0.0,0.0);
	double speed;

  //2. READ XSF FILE WITH POSITIONS AND CREATE OUTPUT FILE FOR MICELLE SPEEDS.__________________
	std::string filename_Neutral = "data/neutral/" + molname + "/" 
	                               + Num_to_String(run) + "_" 
	                               + Num_to_String(lipTot) + "lipids"
							       + Num_to_String(headmonNum) + "headmonNum"
							       + Num_to_String(pHsemonNum) + "pHsemonNum"
								   + Num_to_String(tailmonNum) + "tailmonNum"
                                   + Num_to_String(branchtypemonTot(0)) + "branchheadmonTot"
							       + Num_to_String(branchtypemonTot(1)) + "branchpHsemonTot"
								   + Num_to_String(branchtypemonTot(2)) + "branchtailmonTot"
						           + Num_to_String(drugTot) + "drugs"
							       + Num_to_String(solvTot) + "solvent"
								   + Num_to_String(equilibstepTot) + "steps"; 									   
	
	std::string xsf_file_string = filename_Neutral + ".xsf";	
	const char* xsf_file = xsf_file_string.c_str();	
	
	std::string centre_bin_file_string = filename_Neutral + "_centre.bin";	
    const char* centre_bin_file = centre_bin_file_string.c_str();

	
	//2. SKIP FIRST LINES OF XSF FILE.__________________________________________________________
	std::string line;
	std::ifstream xsf(xsf_file);
	for (int i = 0; i < 12; i++) {
	getline (xsf,line);
	}
	
	//3. ATOM STEPS TO SKIP WHERE NO STABLE MICELLE IS FORMED. 
	for (int j = 1; j <= 649; j++) {
	  getline (xsf,line);
	  for (int i = 0; i < partTot; i++) {
	    getline (xsf,line);
	  }
	}
    
	std::ofstream centre_bin(centre_bin_file, std::ofstream::binary);
	//4. FIND THE SPEED OF THE MICELLE_______________________________________________________
	for (int j = 650; j <= 1001; j++) {
	  getline (xsf,line);
	  headlist.clear();
	  for (int i = 0; i < partTot; i++) {
	    std::string comp;
	    double xx,yy,zz;
	    Vector pos;
	    if (i == 0) {
	      getline(xsf, comp, ' ');
	  
	      getline(xsf, comp, ' ');
	      xx = atof(comp.c_str());  
	      getline(xsf, comp, ' ');
	      yy = atof(comp.c_str());
	      getline(xsf, comp, ' ');
	      zz = atof(comp.c_str());
		
	      pos = Vector(xx,yy,zz);
          //std::cout << pos << std::endl;
		  headlist.push_back(pos);
	    }
	    else {
	      getline(xsf, comp, ' ');
	      xx = atof(comp.c_str());  
	      getline(xsf, comp, ' ');
	      yy = atof(comp.c_str());
	      getline(xsf, comp, ' ');
	      zz = atof(comp.c_str()); 
		  if (i < monNum*lipTot + branchTot && i%(monNum + lipidbranchmonTot) == 0) {
	        pos = Vector(xx,yy,zz);
            //std::cout << i + 1 << " " << pos << std::endl;
			headlist.push_back(pos);
	      }
	    }
	  }
	  
	  Cluster_Centre(lipTot, length, speed, centre, headlist);
	  //std::cout << speed << std::endl;
	  
	  if (j != 650) {
	    centre_bin.write((char*) &speed, sizeof(double));
		std::cout << speed << std::endl;
	  }
    }
	
	centre_bin.close();
	
    return 0;
}
