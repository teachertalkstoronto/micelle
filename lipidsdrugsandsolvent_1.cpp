#include <iostream>     //For output to shell
#include <string>
#include <iomanip>
#include "vecmat3.h"    //For vector/matrix operations
#include <list>	        //To create particle list for each cell and neighbor list for each particle
#include <cmath>        //For mathematical operations
#include "forces.h"     //Force functions//#include <omp.h>	    //Just to record runtimes for various functions using OpenMP
#include "constants_1.h"	//Holds constant valued paramters
#include "classes.h"    //Holds different classes used 
#include "myrandom.h"   //To get random numbers                            
#include "functions.h"  //Functions
#include <fstream>		//Output binary data to binary file, and psf and xsf data
#include "psf.h"		//Holds function to create psf file for lipid and drug particles

using namespace constants;		//Holds constant variables

int main() {

  //1. FILENAME FOR PSF, XSF, BIN, AND TEXT FILES.________________________________________________________________________________
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
    std::string psf_file_string = filename_Neutral + ".psf";	
	std::string xsf_file_string = filename_Neutral + ".xsf";	
	std::string int_bin_file_string = filename_Neutral + "_int.bin";	
	std::string double_bin_file_string = filename_Neutral + "_double.bin";	
	std::string state_string = filename_Neutral + "_state.bin";	
	//std::string txt_file_string = filename_Neutral + ".txt";	
	const char* psf_file = psf_file_string.c_str();	
    const char* xsf_file = xsf_file_string.c_str();	
	const char* int_bin_file = int_bin_file_string.c_str();
	const char* double_bin_file = double_bin_file_string.c_str();
	const char* state_file = state_string.c_str();
	//const char* txt_file = txt_file_string.c_str();
		
  //2. CREATE BEADS AND ASSIGN TYPE, LIPID MONOMER INDICES.________________________________________________________________________
    Bead *particles = new Bead[partTot];
    
	//If drugs are not bonded
	if (drug_bonded == 0 ) {
	   //If order of lipid goes H - P - T
	   if (rev_lip == 0) 
         Assign_Type_Lipid_Monomer_Indices_1<Bead*>(particles,branchtypemonTot,lipidbranchmonTot,headmonNum,pHsemonNum,tailmonNum,monNum,lipTot,branchTot,drugTot,solvTot,partTot);
       //If order of lipid goes H - T - P
	   else
	     Assign_Type_Lipid_Monomer_Indices_2<Bead*>(particles,branchtypemonTot,lipidbranchmonTot,headmonNum,tailmonNum,pHsemonNum,monNum,lipTot,branchTot,drugTot,solvTot,partTot);
    }
	//If drugs are bonded
	else {
	  //If order of lipid goes H - P - T
	  if (rev_lip == 0)
        Assign_Type_Lipid_Monomer_Indices_3<Bead*>(particles,branchtypemonTot,lipidbranchmonTot,headmonNum,pHsemonNum,tailmonNum,monNum,lipTot,branchTot,drugTot,solvTot,partTot);
      //If order of lipid goes H - T - P
	  else
	    Assign_Type_Lipid_Monomer_Indices_4<Bead*>(particles,branchtypemonTot,lipidbranchmonTot,headmonNum,tailmonNum,pHsemonNum,monNum,lipTot,branchTot,drugTot,solvTot,partTot);
    }
	
    /*
    std::cout << "ptype" << std::setw(10) << "lipid" << std::setw(10) << "monomer" << std::setw(10) << "ind" << std::endl;
    for (int i = 0; i < partTot; i++) {
      std::cout << particles[i].ptype << std::setw(10) << particles[i].lipid << std::setw(10) << particles[i].monomer << std::setw(10) << particles[i].ind << std::endl;
    }
    std::cout << std::endl;
    */
 
  //3. FIX THE INITIAL VELOCITIES OF BEADS, AND FIND ekTot OF SYSTEM.____________________________________________________________________________________
    Initial_Velocities<Bead*>(particles, partTot, mh, mp, mt, md, ms, kT);
    Max_Speed_Avg_Speed<Bead*>(particles, solvTot, partTot, v_avg, v_avg_solv, vmax, vmax_solv);	
    KineticEnergy_LinearMomentum<Bead*>(particles, solvTot, partTot, ekTot, lin_momTot, mh, mp, mt, md, ms);

    /*
    std::cout << "Velocities" << std::endl;
    for (int i = 0; i < partTot; i++) {
      std::cout << particles[i].ptype << std::setw(15) << particles[i].v << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Initial kinetic energy " << std::setw(15) << ekTot << std::endl;
    std::cout << std::endl;
    */
  
  //4. FIX THE INITIAL POSITIONS OF BEADS, AND CLUSTER ANALYSIS OF DRUGS BASED ON DRUG POSITION.________________________________________________________________________
    LatticePos *latticepos_lip = new LatticePos[lpsites_lip];
    FCC_Lattice_Initial_Positions<LatticePos*>(latticepos_lip, lpspacing_lip, lpdim_lip, length);
    Lipid_Initial_Positions<Bead*, LatticePos*>(particles,latticepos_lip,drug_bonded,branchNum,branchmonNum,prevbranchmonTot,lipidbranchmonTot,monNum,lipTot,branchTot,drugTot,lpsites_lip,lpspacing_lip,length); 
    delete [] latticepos_lip;
    NonLipid_Initial_Positions_2<Bead*>(particles,drug_bonded,monNum,lipTot,branchTot,drugTot,solvTot,partTot,lpsites,lpspacing,length); 


    /*
    std::cout << "Positions" << std::endl;
    std::cout << "ptype" << std::setw(15) << "Positions" << std::endl;
    std::cout << particles[0].ptype << std::setw(15) << particles[0].r << std::endl;
    for (int i = 1; i < partTot; i++) {
      std::cout << particles[i].ptype << std::setw(15) << particles[i].r << std::setw(15) << dist(particles[i].r, particles[i - 1].r) << std::endl;
    }
    std::cout << std::endl;
    */

    Cluster_Drug_Analysis<Bead,Bead*,Cluster,std::list<Cluster>,std::list<Cluster>::iterator>(particles, 
	                                                                                          drug_bonded, drug_released,
																							  mon1,mon2,mon3,mon4,monNum,lipTot,lipidbranchmonTot,branchTot,drugTot,length,
																							  rclust,clustNum,clustLipNum_avg,clustLipNum_var,clustRad_avg,clustRad_var,micRad_avg,micRad_var, 
																							  drugIn,drugOut,clustDrugNum_avg,clustDrugNum_var,drugInRad_div_clustRad_avg,drugInRad_div_clustRad_var);
							
  //5. CREATE MPC CELLS WITHIN DOMAIN._____________________________________________________________________________________________________________________________________________
    MPC_Cell ***mpc_cell;
    alloc3Darray<MPC_Cell***,MPC_Cell**,MPC_Cell*,MPC_Cell>(mpc_cell,mpc_xind);
  
    /*
    for (int i = 0; i < mpc_xind; i++) {
      for (int j = 0; j < mpc_xind; j++) {
        for (int k = 0; k < mpc_xind; k++) {
	      if (mpc_cell[i][j][k].partlist.size() > 0) {
 	        std::cout << "For MPC cell " << i << " " << j << " " << k << " the initial particle list is:" << std::endl;
            for (std::list<Bead>::iterator it = mpc_cell[i][j][k].partlist.begin(); it != mpc_cell[i][j][k].partlist.end(); ++it) {
	          std::cout << (*it).r << std::setw(4) << (*it).ptype << std::setw(10) << (*it).mpc_cell << std::endl;
            }
		  } 
        }
      }
    }
    std::cout << std::endl;
    */
  
  //6. MAKE PARTICLE LIST FOR EACH CELL IN DOMAIN._______________________________________________________________________________________________________________________________
    Cell ***cell;
    alloc3Darray<Cell***, Cell**, Cell*, Cell>(cell, xind);
  
    Particle_Lists<Cell***, Bead*>(cell, particles, partTot, xind, cellsize);	
  
    /*
    for (int i = 0; i < xind; i++) {
      for (int j = 0; j < xind; j++) {
        for (int k = 0; k < xind; k++) {
	      if (cell[i][j][k].partlist.size() > 0) {
 	        std::cout << "For cell " << i << " " << j << " " << k << " the initial particle list is:" << std::endl;
            for (std::list<Bead>::iterator it = cell[i][j][k].partlist.begin(); it != cell[i][j][k].partlist.end(); ++it) {
	          std::cout << (*it).r << std::setw(4) << (*it).ind << std::setw(10) << (*it).cell << std::endl;
            }
		  } 
        }
      }
    }
    std::cout << std::endl;
    */ 

  //7. FIND ACCELERATION OF EACH BEAD BASED ON IC, and uTot, and pressTot._______________________________________________________________________________________
    if (drug_bonded == 0) {
      //Drug not attached
      Acceleration_uTot_FreeDrug_Neutral<Cell***,Bead*,std::list<Bead>::iterator,Bead>(cell, particles,
                                                                                       branchNum, branchmonNum, prevbranchmonTot, lipidbranchmonTot,
                                                                                       monNum, lipTot, branchTot, solvTot, partTot, xind, uTot, pressTot, length,
                                                                                       mh, mp, mt, md, ms,
                                                                                       sigma, rc, rinf, kbond, kbend, kT,
								                                                       ehh, ehp, eht, ehd, ehs, 
											                                           epp, ept, epd, eps, 																		 
																                       ett, etd, ets, 
																                       edd, eds,
								                                                       whs, 
																                       wpp, wpt, wpd, 																		 
																                       wtt, wtd, 
																                       wdd);	

    }
    else {
      //Drug attached
      Acceleration_uTot_BondedDrug_Neutral<Cell***, Bead*, std::list<Bead>::iterator, Bead>(cell, particles,
                                                                                            branchNum, branchmonNum, prevbranchmonTot, lipidbranchmonTot,
                                                                                            monNum, lipTot, branchTot, solvTot, partTot, xind, uTot, pressTot, length,
                                                                                            mh, mp, mt, md, ms,
                                                                                            sigma, rc, rinf, kbond, kbend, kT,
								                                                            ehh, ehp, eht, ehd, ehs, 
											                                                epp, ept, epd, eps, 																		 
																                            ett, etd, ets, 
																                            edd, eds,
								                                                            whs, 
																                            wpp, wpt, wpd, 																 
																                            wtt, wtd, 
																                            wdd);	
    }
 
    /*  
    for (int i = 0; i < partTot; i++) {
      std::cout << "For particle " << i << " of type " << particles[i].ptype << " the initial acceleration is " << particles[i].a << std::endl;
    }
    std::cout << std::endl;
    std::cout << uTot << std::endl;
 
   for (int i = 0; i < partTot; i++) {
     std::cout << particles[i].ptype << " " << i << std::endl;
   }
   */
 
  //8. TOTAL ENERGY OF THE SYSTEM.______________________________________________________________________________________________________________________________________________________________________											   
    step = 0;
	TIME = 0.0;
	eTot = ekTot + uTot;
    std::cout << "The initial total kinetic energy of the system: " << ekTot << std::endl;
    std::cout << "The initial total potential energy of the system: " << uTot << std::endl;
    std::cout << "Initial energy of the system: " << eTot << std::endl;
    std::cout << "Initial linear momentum of the system: " << lin_momTot << std::endl;std::cout << std::endl;
    std::cout << std::endl;
  
  //9. MAKES FILES.___________________________________________________________________________________________________________________________________________
    //Make psf file._________________________________________________________________
    Output_psf_File<Bead*>(psf_file,particles,monNum,lipTot,branchTot,drugTot,solvTot,branchNum,branchmonNum,prevbranchmonTot,lipidbranchmonTot,mh,mp,mt,md,ms,kT/epsilon);
  
    //Make xsf file._________________________________________________________________ 		
    std::ofstream xsf(xsf_file);
    xsf << "# " << "Lipids: "<< lipTot << std::endl
        << "# " << "Branched monomers: " << branchTot << std::endl
        << "# " << "Drugs: " << drugTot << std::endl
        << "# " << "Solvent: " << solvTot << std::endl
	    << "# " << "kT/eps: " << kT << std::endl
	    << "# " << "Timesteps: " << equilibstepTot << std::endl
	    << "# " << "xsf steps: " << xsfstep <<std::endl
        << "PRIMEVEC" << std::endl
        << length*Vector(1.0, 0.0, 0.0) << std::endl
        << length*Vector(0.0, 1.0, 0.0) << std::endl
	    << length*Vector(0.0, 0.0, 1.0) << std::endl
	    << "ANIMSTEPS " << 1 + floor(equilibstepTot/xsfstep) << std::endl
	    << "ATOMS " << 1 << std::endl;
    for (int i = 0; i < partTot; i++) {
      xsf << i + 1 << " " << std::setprecision(16) << particles[i].r << std::endl;
    }
    xsf.close();
  
    //Write int bin files._____________________________________
    std::ofstream int_bin(int_bin_file, std::ofstream::binary);
    int_bin.write((char*) &step, sizeof(int));
    int_bin.write((char*) &clustNum, sizeof(int));	
    int_bin.write((char*) &drugIn, sizeof(int));
    int_bin.write((char*) &drugOut, sizeof(int));
    int_bin.close();
	
	//Write double bin files._____________________________________
	std::ofstream double_bin(double_bin_file, std::ofstream::binary);
    double_bin.write((char*) &TIME, sizeof(double));
    double_bin.write((char*) &v_avg, sizeof(double));
    double_bin.write((char*) &v_avg_solv, sizeof(double));
    double_bin.write((char*) &vmax, sizeof(double));
    double_bin.write((char*) &vmax_solv, sizeof(double));
    double_bin.write((char*) &ekTot, sizeof(double));
    double_bin.write((char*) &uTot, sizeof(double));
    double_bin.write((char*) &pressTot, sizeof(double));
    double_bin.write((char*) &eTot, sizeof(double));
    double_bin.write((char*) &lin_momTot(0), sizeof(double));
	double_bin.write((char*) &lin_momTot(1), sizeof(double));
	double_bin.write((char*) &lin_momTot(2), sizeof(double));
    double_bin.write((char*) &clustLipNum_avg, sizeof(double));
	double_bin.write((char*) &clustLipNum_var, sizeof(double));
	double_bin.write((char*) &clustRad_avg, sizeof(double));
	double_bin.write((char*) &clustRad_var, sizeof(double));
    double_bin.write((char*) &micRad_avg, sizeof(double));
	double_bin.write((char*) &micRad_var, sizeof(double));
	double_bin.write((char*) &clustDrugNum_avg, sizeof(double));
	double_bin.write((char*) &clustDrugNum_var, sizeof(double));
	double_bin.write((char*) &drugInRad_div_clustRad_avg, sizeof(double));
	double_bin.write((char*) &drugInRad_div_clustRad_var, sizeof(double));
    double_bin.close();

  //Write IC to text file._________________________________________________________
  /*
  std::ofstream txtdata(txt_file);
  
  std::ifstream txtdatain("constants.h");
  std::string line;
  while (getline(txtdatain, line)) {
    txtdata << line;
  }
  txtdatain.close();
  
  txtdata << std::setw(15) << "step" 
          << std::setw(15) << "time"   
		  << std::setw(15) << "v_avg" 
	      << std::setw(15) << "v_avg_solv" 
		  << std::setw(15) << "vmax" 
	      << std::setw(15) << "vmax_solv" 
	      << std::setw(15) << "ekTot" 
	      << std::setw(15) << "uTot" 
		  << std::setw(15) << "pressTot" 
	      << std::setw(15) << "eTot" 
	      << std::setw(31) << "lin_momTot" 
		  << std::setw(15) << "drugIn"
	      << std::setw(15) << "drugOut"
		  << std::setw(15) << "drugIn_percent"
		  << std::setw(15) << "drugOut_percent"
		  << std::setw(15) << "clustNum"
		  << std::setw(15) << "clustSize_avg"
	      << std::setw(15) << "MPCcellNum_with_no_nonsolv"
		  << std::setw(15) << "solvNum_avg_MPCcell_with_no_nonsolv"
	      << std::endl;    
  txtdata << std::setw(15) << step 
	      << std::setw(15) << TIME 
		  << std::setw(15) << v_avg
	      << std::setw(15) << v_avg_solv 
		  << std::setw(15) << vmax 
	      << std::setw(15) << vmax_solv 
		  << std::setw(15) << ekTot 
		  << std::setw(15) << uTot 
		  << std::setw(15) << pressTot 
		  << std::setw(15) << eTot 
		  << std::setw(15) << lin_momTot 
		  << std::setw(15) << drugIn 
	      << std::setw(15) << drugOut 
		  << std::setw(15) << drugIn_percent 
		  << std::setw(15) << drugOut_percent
		  << std::setw(15) << clustNum
		  << std::setw(15) << clustSize_avg
	    //<< std::setw(15) << pressTot 
		  << std::endl;
  */
  
  //10. STORE DATA IN ARRAYS.____________________________________________________________________________________________________
    Vector *xsf_data = new Vector[(equilibstepTot/xsfstep)*partTot];
    int *int_data = new int[equilibstepTot*intdataNum];	
    double *double_data = new double[equilibstepTot*doubledataNum];
    //mpc________________________________________________________
    int *mpc_int_data = new int[equilibstepTot/tau];    
	double *mpc_double_data = new double[2*equilibstepTot/tau];
  
  //11. TIME EVOLUTION FROM IC TO EQUILIBRIUM._________________________________________________________
    for (step = 1; step <= equilibstepTot; step++) {
      //SRD___________________________________________________________________________________________________________________________________________________________
      if ((step%tau) == 0) {
	    Grid_Shift_Assign_vrand<Bead*>(particles, solvTot, partTot, mpc_xind, mpc_cellsize, length, kT, ms);
	    Calc_vmean<MPC_Cell***, Bead*>(mpc_cell, particles, solvTot, partTot, mpc_xind, MPCcellNum_with_no_nonsolv, solvNum_avg_MPCcell_with_no_nonsolv, solvNum_var_MPCcell_with_no_nonsolv, mpc_cellsize);	
	    Update_Solvent_Velocities<MPC_Cell***, Bead*>(mpc_cell, particles, solvTot, partTot);
	    Return_Positions<Bead*>(particles, solvTot, partTot);
	  }
  
	  //1. For full timestep update positions r(time) to r(time + dt) via TS expansion.
      int i;
	  //#pragma omp parallel default(none) shared(partTot, particles, dt, length,t0,m,epsilon,sigma) private(i) 
	  //{
	    //#pragma omp for
	    for (i = 0; i < partTot; i++) {
          particles[i].r += (particles[i].v)*dt  + 0.5*(particles[i].a)*dt*dt;  
          Boundary_Conditions<Bead*>(particles, i, length);
        }
      //}
	
	  //Cluster analysis after new positions.________________________________________																							
      Cluster_Drug_Analysis<Bead,Bead*,Cluster,std::list<Cluster>,std::list<Cluster>::iterator>(particles, 
	                                                                                            drug_bonded, drug_released,
																							    mon1,mon2,mon3,mon4,monNum,lipTot,lipidbranchmonTot,branchTot,drugTot,length,
																							    rclust,clustNum,clustLipNum_avg,clustLipNum_var,clustRad_avg,clustRad_var,micRad_avg,micRad_var, 
																								drugIn,drugOut,clustDrugNum_avg,clustDrugNum_var,drugInRad_div_clustRad_avg,drugInRad_div_clustRad_var);
																								
	  //2. For half timestep update velocities v(time) to v(time + 0.5*dt) via TS expansion.
	  //#pragma omp parallel default(none) shared(partTot, particles, dt) private(i) 
	  //{
	    //#pragma omp for
        for (i = 0; i < partTot; i++) {
          particles[i].v += 0.5*(particles[i].a)*dt;
	    }
      //}
	
	  //3. For full timestep update accelerations a(time) to a(time + dt) via a=f/m
	  Particle_Lists<Cell***, Bead*>(cell, particles, partTot, xind, cellsize);
	  
      if (drug_bonded == 0) {
        //Drug not attached
        Acceleration_uTot_FreeDrug_Neutral<Cell***,Bead*,std::list<Bead>::iterator,Bead>(cell, particles,
                                                                                         branchNum, branchmonNum, prevbranchmonTot, lipidbranchmonTot,
                                                                                         monNum, lipTot, branchTot, solvTot, partTot, xind, uTot, pressTot, length,
                                                                                         mh, mp, mt, md, ms,
                                                                                         sigma, rc, rinf, kbond, kbend, kT,
								                                                         ehh, ehp, eht, ehd, ehs, 
											                                             epp, ept, epd, eps, 																		 
																                         ett, etd, ets, 
																                         edd, eds,
								                                                         whs, 
																                         wpp, wpt, wpd, 																		 
																                         wtt, wtd, 
																                         wdd);	

      }
      else {
        //Drug attached
        Acceleration_uTot_BondedDrug_Neutral<Cell***, Bead*, std::list<Bead>::iterator, Bead>(cell, particles,
                                                                                              branchNum, branchmonNum, prevbranchmonTot, lipidbranchmonTot,
                                                                                              monNum, lipTot, branchTot, solvTot, partTot, xind, uTot, pressTot, length,
                                                                                              mh, mp, mt, md, ms,
                                                                                              sigma, rc, rinf, kbond, kbend, kT,
								                                                              ehh, ehp, eht, ehd, ehs, 
											                                                  epp, ept, epd, eps, 																		 
																                              ett, etd, ets, 
																                              edd, eds,
								                                                              whs, 
																                              wpp, wpt, wpd, 																		 
																                              wtt, wtd, 
																                              wdd);	
      }
	
	  //4. For Full timestep update velocities v(time) to v(time + dt) via TS expansion.
      //#pragma omp parallel default(none) shared(partTot, particles, dt) private(i) 
	  //{
	    //#pragma omp for
	    for (i = 0; i < partTot; i++) {
          particles[i].v += 0.5*(particles[i].a)*dt;       
        }
	  //}
	
	  //Calculate KE after new speeds.____________________________________________________________________________
	  //vmax = Max_Speed<Bead**>(particles, partTot);	//Max velocity of system
	  Max_Speed_Avg_Speed<Bead*>(particles, solvTot, partTot, v_avg, v_avg_solv, vmax, vmax_solv);	
      KineticEnergy_LinearMomentum<Bead*>(particles, solvTot, partTot, ekTot, lin_momTot, mh, mp, mt, md, ms);
	  
	  eTot = ekTot + uTot;
	  TIME += dt;  //Update time

	  //Put xsf data tp xsf array._________________________________
	  if ((step%xsfstep) == 0) {
	    int xsfstepcount = partTot*((step/xsfstep) - 1);
        for (int i = 0; i < partTot; i++) {
	      xsf_data[xsfstepcount + i] = particles[i].r;
        }  
	  }
	  //Put int data to array._____________________________________
	  int int_stepcount = intdataNum*(step-1);
	  int_data[int_stepcount] = step;
      int_data[int_stepcount+1] = clustNum;	  
      int_data[int_stepcount+2] = drugIn;
      int_data[int_stepcount+3] = drugOut;

	  //Put double data to array.__________________________________
	  int double_stepcount = doubledataNum*(step-1);
      double_data[double_stepcount] = TIME;
	  double_data[double_stepcount+1] = v_avg;
      double_data[double_stepcount+2] = v_avg_solv;
      double_data[double_stepcount+3] = vmax;
      double_data[double_stepcount+4] = vmax_solv;
      double_data[double_stepcount+5] = ekTot;
      double_data[double_stepcount+6] = uTot;
      double_data[double_stepcount+7] = pressTot;
      double_data[double_stepcount+8] = eTot;
      double_data[double_stepcount+9] =lin_momTot(0);
	  double_data[double_stepcount+10] =lin_momTot(1);
	  double_data[double_stepcount+11] =lin_momTot(2);
      double_data[double_stepcount+12] = clustLipNum_avg;
      double_data[double_stepcount+13] = clustLipNum_var;
	  double_data[double_stepcount+14] = clustRad_avg;
	  double_data[double_stepcount+15] = clustRad_var;
	  double_data[double_stepcount+16] = micRad_avg;
	  double_data[double_stepcount+17] = micRad_var;
	  double_data[double_stepcount+18] = clustDrugNum_avg;
	  double_data[double_stepcount+19] = clustDrugNum_var;	  
	  double_data[double_stepcount+20] = drugInRad_div_clustRad_avg;
	  double_data[double_stepcount+21] = drugInRad_div_clustRad_var;	 
	  //Put mpc int and double data to arrays._____________________
      if ((step%tau) == 0) {
	    mpc_int_data[(step/tau) - 1] = MPCcellNum_with_no_nonsolv;
	    mpc_double_data[2*(step/tau) - 2] = solvNum_avg_MPCcell_with_no_nonsolv;
		mpc_double_data[2*(step/tau) - 1] = solvNum_var_MPCcell_with_no_nonsolv;
	  }
	  
	  /*
	  //Put text data to text array._________________________________________________
      txtdata << std::setw(15) << step 
	          << std::setw(15) << TIME 
		      << std::setw(15) << v_avg
	          << std::setw(15) << v_avg_solv 
		      << std::setw(15) << vmax 
	          << std::setw(15) << vmax_solv 
		      << std::setw(15) << ekTot 
		      << std::setw(15) << uTot 
		      << std::setw(15) << pressTot 
		      << std::setw(15) << ekTot + uTot 
		      << std::setw(15) << lin_momTot 
		      << std::setw(15) << drugIn 
	          << std::setw(15) << drugOut 
		      << std::setw(15) << drugIn_percent 
		      << std::setw(15) << drugOut_percent
		      << std::setw(15) << clustNum
		      << std::setw(15) << clustSize_avg;
	  if ((step%tau) == 0) {
	    txtdata << std::setw(15) << MPCcellNum_with_no_nonsolv
		        << std::setw(15) << solvNum_avg_MPCcell_with_no_nonsolv;
	  }
	  txtdata << std::endl;
      */
	  
    }
  
  //12. WRITE FILES.__________________________________________________________________________________________________________________
    //Write xsf file_______________________________________________
    xsf.open(xsf_file, std::ofstream::app);
    for (int i = 1; i <= equilibstepTot/xsfstep; i++) {
      xsf << "ATOMS " << 1 + i << std::endl;
	  int xsfstepcount = partTot*(i - 1);
	  for (int j = 0; j < partTot; j++) {
        xsf << j + 1 << " " << xsf_data[xsfstepcount + j] << std::endl;
      }  
    }	
    xsf.close();
    delete [] xsf_data;
  
    /*
    std::string testfile = "data/textdata/tailpendant"+Num_to_String(lipTot)+"lipids"+Num_to_String(drugTot)+"drugs"+Num_to_String(acidTot)+"acid"+Num_to_String(solvTot)+"solvent"+Num_to_String(equilibstepTot)+"steps.txt";
    const char* testfile2 = testfile.c_str();
    std::ofstream testtxtdata(testfile2);
    for (int i = 0; i < partTot; i++) {
      testtxtdata << particles[i].v << std::endl;
    }
   testtxtdata.close();
 
    std::string postestfile = "data/textdata/postailpendant"+Num_to_String(lipTot)+"lipids"+Num_to_String(drugTot)+"drugs"+Num_to_String(acidTot)+"acid"+Num_to_String(solvTot)+"solvent"+Num_to_String(equilibstepTot)+"steps.txt";
    const char* postestfile2 = postestfile.c_str();
    std::ofstream postesttxtdata(postestfile2);
    for (int i = 0; i < partTot; i++) {
      postesttxtdata << particles[i].r << std::endl;
    }
   postesttxtdata.close();
   */
 
    //Write int binary data file_________________________________________
	int_bin.open(int_bin_file, std::ofstream::app|std::ofstream::binary);
    for (int i = 0; i < equilibstepTot*intdataNum; i++) { 
      int_bin.write((char*) &int_data[i], sizeof(int));
    }
    delete [] int_data;
    for (int i = 0; i < equilibstepTot/tau; i++) { 
      int_bin.write((char*) &mpc_int_data[i], sizeof(double));
    }
	int_bin.close();
    delete [] mpc_int_data;
	
	//Write double binary data file___________________________________________
	double_bin.open(double_bin_file, std::ofstream::app|std::ofstream::binary);
	for (int i = 0; i < equilibstepTot*doubledataNum; i++) {
      double_bin.write((char*) &double_data[i], sizeof(double));
    }
    delete [] double_data;
    for (int i = 0; i < 2*equilibstepTot/tau; i++) {
      double_bin.write((char*) &mpc_double_data[i], sizeof(double));
    }
    double_bin.close();
	delete [] mpc_double_data;
    
	//Write state of system for acidification._________________
	std::ofstream state(state_file);
    for (int j = 0; j < partTot; j++) {
      state.write((char*) &(particles[j].r(0)), sizeof(double));
	  state.write((char*) &(particles[j].r(1)), sizeof(double));
	  state.write((char*) &(particles[j].r(2)), sizeof(double));
    }
    for (int j = 0; j < partTot; j++) {
  	  state.write((char*) &(particles[j].v(0)), sizeof(double));
	  state.write((char*) &(particles[j].v(1)), sizeof(double));
	  state.write((char*) &(particles[j].v(2)), sizeof(double));
    }
	state.close();
	
    //txtdata.close();	//Close file holding configuration items of system. 

  //13. DEALLOCATE MEM FOR REMAINIGN ARRAYS._____________________________________________________________
    dealloc3Darray<MPC_Cell***>(mpc_cell, mpc_xind);
    dealloc3Darray<Cell***>(cell, xind);  //For 3D array cell of class Cell
    delete [] particles;
  
  return 0;
}
