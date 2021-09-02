//Constants

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <iostream>
#include <string>
#include <cmath>
#include <algorithm> 
#include "vecmat3.h"  


namespace constants {

  std::string molname = "control";
  int run = 1; 
  
  int intdataNum = 4;
  int doubledataNum = 22;

  int drug_bonded = 0; //Bool
  int drug_released = 0; //bool if drug bonded == 1 && protonated == 1
  int rev_lip = 0;     //bool
  
  //DIMENSIONLESS UNITS___________________________________________________________________
  double m = 1.0;                                                               //Mass
  double sigma = 1.0;	                                                        //Distance
  double epsilon = 1.0;	                                                            //Energy
  double t0 = sqrt(m*sigma*sigma*(1.0/epsilon));                                          //Time

  //TIMESTEP SIZES AND NUMBER OF TIMESTEPS________________________________________________
  double dt = 0.005*t0;                                                              //MD
  //double highdt = 0.0005*t0;
  int tau = 40;                                   
  //MPC
  int equilibstepTot = 2000000;
  int acidicstepTot = 40000;
  int xsfstep = 400;

  //MASS OF EACH TYPE OF ATOM/BEAD________________________________________________________
  double mh = m;
  double mp = m;  
  double mt = m;
  double md = m;
  double ms = 0.5*m;

  //TEMPERATURE___________________________________________________________________________
  double kT = 1.0*epsilon;
  
  //INTERACTION RANGES____________________________________________________________________
  double rc = pow(2.0, 1.0/6.0)*sigma;	            //Cut-off distance for LJ interaction.   
  double rinf = 1.5*sigma;                          //Maximum monomer bond length in lipid.
  double rclust = 2.6*sigma;  
  
  double whs = 1.65*sigma - rc;	
  double wtt = 2.6*sigma - rc;
  double wtd = 2.6*sigma - rc;
  double wdd = 2.6*sigma - rc;
  
  /*
  //Neutral pH_________________
  double wpp = 2.6*sigma - rc;   
  double wpt = 2.6*sigma - rc;   
  double wpd = 2.6*sigma - rc;
  //double wpd_bonded = 1.5*sigma - rc;   //NEWWWWWWWWWWWWWW BUT NEED TO CONTROL WHERE DRUG IS IN LATTICE IN BEG
  */
 
  //Acidic pH__________________    
  double wps = 1.65*sigma - rc;

 
  //LENNARD-JONES INTERACTION ENERGIES____________________________________________________
  double ehh = 0.5*epsilon;	                            //Interaction energy b/w h and h bead.
  double eht = epsilon; 
  double ehd = epsilon;
  double ehs = 0.05*epsilon;
  double ett = 0.5*epsilon;	                            //Interaction energy b/w t and t bead.
  double etd = 0.5*epsilon;
  double ets = 2.0*epsilon;	                            //Interaction energy b/w s and t bead.
  double edd = 0.5*epsilon; 
  double eds = 2.0*epsilon; 
  
  /*
  //Neutral pH:____________
  double ehp = epsilon;           
  double epp = 0.5*epsilon;          
  double ept = 0.5*epsilon;      
  double epd = 0.5*epsilon;            
  double eps = 2.0*epsilon;          
  */
   
  //Acidic conditions:____________
  double ehp = 0.5*epsilon; 
  double epp = 0.5*epsilon; 
  double ept = 1.0*epsilon;   
  double epd = 1.0*epsilon; 
  double eps = 0.05*epsilon;    
  
   
  //SPRING CONSTANTS______________________________________________________________________
  double kbond = 20.0*epsilon*(1.0/(sigma*sigma));
  double kbend = 2.5*epsilon*(1.0/(sigma*sigma));
  //double ktor = 0.3*epsilon;

  //SYSTEM DIMENSIONS_____________________________________________________________________ 
  double cellsize = 2.6*sigma;     //MD cell edge length
  int xind = 10;                   //x,y,z-direction cell indices go from 0 to 9.
  double mpc_cellsize = sigma;     //Choose this length such that xind*cellsize/mpc_cellsize gives an integer number
  int mpc_xind = 26;
  double length = xind*cellsize;   //or equivalently double length = mpc_xind*mpc_cellsize;
  
  double lpspacing = 0.65*sigma;
  int lpdim = 40;
  int lpsites = 4*pow(lpdim, 3);         //bl: 0.4596
  
  double lpspacing_lip = 1.3*sigma;
  int lpdim_lip = 20;
  int lpsites_lip = 4*pow(lpdim_lip, 3);      //bl: 0.919
  
  //SYSTEM SIZE___________________________________________________________________________
  int headmonNum = 1;
  int pHsemonNum = 1;
  int tailmonNum = 3;
  int monNum = 5;
  
  int branchNum[5] = {0,0,0,0,0};   //Set as 0,0,0,0,0, is branchTot is 0
  int branchmonNum[5][1] = {{0},{0},{0},{0},{0}}; //If branchTot = 0						  
  int prevbranchmonTot[5][1] = {{0},{0},{0},{0},{0}};
  vecmat3::Vector<int> branchtypemonTot(0,0,0);		//Must label as 0, 0 ,0 if banchTot 0;
  int lipidbranchmonTot = branchtypemonTot(0) + branchtypemonTot(1) + branchtypemonTot(2);
  
  /*
  int branchNum[5] = {0,0,0,0,0};   //Set as 0,0,0,0,0, is branchTot is 0
  int branchmonNum[5][1] = {{0},{0},{0},{0},{0}}; //If branchTot = 0						  
  int prevbranchmonTot[5][1] = {{0},{0},{0},{0},{0}};
  vecmat3::Vector<int> branchtypemonTot(0,0,0);		//Must label as 0, 0 ,0 if banchTot 0;
  int lipidbranchmonTot = branchtypemonTot(0) + branchtypemonTot(1) + branchtypemonTot(2);
  */
  
  int lipTot = 100;                     //Total number of h1, t1, t2, or t3 particles/beads.    //Total number of particles/beads of class Tail
  int branchTot = lipidbranchmonTot*lipTot;
  int drugTot = 100;
  //int drugTot = pHsemonNum*lipTot;
  //int acidTot = pHsemonNum*lipTot;
  int solvTot  = 68194 - (monNum*lipTot + branchTot + drugTot); 
  //int solvTot = 30000 - (monNum*lipTot + branchTot + drugTot + acidTot); 
  //int solvTot = 68194 - (monNum*lipTot + branchTot + drugTot + acidTot); 
  //int solvTot = lpsites - (monNum*lipTot + branchTot + drugTot + acidTot); 
 //int solvTot = 8;
  //68984 beads tot in 26.1 by 26.1 by 26.1 system
  //16384
  int partTot = (monNum*lipTot) + branchTot + drugTot + solvTot; //Total number of particles/beads.
  
  //PARAMTERS_______________________________________________________
  int mon1 = 1;
  int mon2 = monNum + lipidbranchmonTot;
  int mon3 = 1;
  int mon4 = monNum + lipidbranchmonTot;
 
  //INT DATA TO RECORD
  int step;
  int clustNum;
  int drugIn;
  int drugOut;  
  int MPCcellNum_with_no_nonsolv;  
      
  //DOUBLE DATA TO RECORD	  
  double TIME;
  double v_avg;
  double v_avg_solv;   
  double vmax;
  double vmax_solv;  
  double ekTot;  
  double uTot;
  double pressTot;  
  double eTot;  
  Vector lin_momTot;  
  double clustLipNum_avg;
  double clustLipNum_var;
  double clustRad_avg;
  double clustRad_var;
  double micRad_avg;
  double micRad_var;
  double clustDrugNum_avg; 
  double clustDrugNum_var;
  double drugInRad_div_clustRad_avg;
  double drugInRad_div_clustRad_var;
  double solvNum_avg_MPCcell_with_no_nonsolv;
  double solvNum_var_MPCcell_with_no_nonsolv;

            //Cutoff distance for particles to be considered of same cluster.
  
  //double solvNum_avg_MPCcell = solvTot/(mpc_xind*mpc_xind*mpc_xind);
  //int nonzerocells;
  //double solvNum_avg_MPC_cell_nonzero;

  


}

#endif

