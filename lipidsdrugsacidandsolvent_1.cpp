#include <iostream>     //For output to shell
#include <string>
#include <iomanip>
#include "vecmat3.h"    //For vector/matrix operations
#include <list>	        //To create particle list for each cell and neighbor list for each particle
#include <cmath>        //For mathematical operations
#include "forces.h"     //Force functions//#include <omp.h>	    //Just to record runtimes for various functions using OpenMP
#include "constants_acidic_1.h"	//Holds constant valued paramters
#include "classes.h"    //Holds different classes used 
#include "myrandom.h"   //To get random numbers                             
#include "functions.h"  //Functions
#include <fstream>		//Output binary data to binary file, and psf and xsf data
#include "psf.h"		//Holds function to create psf file for lipid and drug particles

using namespace constants;		//Holds constant variables

int main() {

  //1. FILENAME FOR PSF, XSF, BIN, AND TEXT FILES.________________________________________________________________________________
    std::string state_string = "data/neutral/" + molname + "/" 
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
							   + Num_to_String(equilibstepTot) + "steps"
							   + "_state.bin"; 
    const char* state_file = state_string.c_str();

    std::string filename_Acidic = "data/acidic/" + molname + "/" 
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
								  + Num_to_String(acidicstepTot) + "steps"; 
	std::string psf_file_string = filename_Acidic + ".psf";	
	std::string xsf_file_string = filename_Acidic + ".xsf";	
	std::string int_bin_file_string = filename_Acidic + "_int.bin";	
	std::string double_bin_file_string = filename_Acidic + "_double.bin";		
	const char* psf_file = psf_file_string.c_str();	
    const char* xsf_file = xsf_file_string.c_str();	
	const char* int_bin_file = int_bin_file_string.c_str();
	const char* double_bin_file = double_bin_file_string.c_str();

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
	
  //3. ASSIGN POSITIONS AND VELOCITIES BASED ON LAST STATE OF SYSTEM IN NEUTRAL CONDITIONS.
    std::ifstream state (state_file, std::ios::binary | std::ios::ate);	
    int state_file_size = state.tellg();
    std::cout << state_file_size << std::endl;
    char *memblock = new char[state_file_size];
    state.seekg(0, std::ios::beg);
    state.read (memblock, state_file_size);
    state.close();
    double* state_double_values = (double*)memblock;//reinterpret as doubles
	for (int i = 0; i < partTot; i++) {
      particles[i].r = Vector(state_double_values[3*i], state_double_values[3*i + 1], state_double_values[3*i + 2]);
    }
    for (int i = 0; i < partTot; i++) {
      particles[i].v = Vector(state_double_values[3*partTot + 3*i], state_double_values[3*partTot + 3*i + 1], state_double_values[3*partTot +3*i + 2]);
    }
	
    //Centre the micelle
	for (int i = 0; i < partTot; i++) {
        particles[i].r += Vector(-6.5,6.5,6.5);  
        Boundary_Conditions<Bead*>(particles, i, length);
    }

			  
    Max_Speed_Avg_Speed<Bead*>(particles, solvTot, partTot, v_avg, v_avg_solv, vmax, vmax_solv);	
    KineticEnergy_LinearMomentum<Bead*>(particles, solvTot, partTot, ekTot, lin_momTot, mh, mp, mt, md, ms);


	for (int i = 0; i < partTot; i++) {
      std::cout << particles[i].r << std::endl;
    }
	/*
	
	for (int i = 0; i < partTot; i++) {
      std::cout << particles[i].v << std::endl;
    }
	*/
  //4. CLUSTER ANALYSIS BASED ON POSITIONS
    Cluster_Drug_Analysis<Bead,Bead*,Cluster,std::list<Cluster>,std::list<Cluster>::iterator>(particles, 
	                                                                                          drug_bonded, drug_released,
																							  mon1,mon2,mon3,mon4,monNum,lipTot,lipidbranchmonTot,branchTot,drugTot,length,
																							  rclust,clustNum,clustLipNum_avg,clustLipNum_var,clustRad_avg,clustRad_var,micRad_avg,micRad_var, 
																							  drugIn,drugOut,clustDrugNum_avg,clustDrugNum_var,drugInRad_div_clustRad_avg,drugInRad_div_clustRad_var);
					
  //5. CREATE MPC CELLS WITHIN DOMAIN._____________________________________________________________________________________________________________________________________________
    MPC_Cell ***mpc_cell;
    alloc3Darray<MPC_Cell***,MPC_Cell**,MPC_Cell*,MPC_Cell>(mpc_cell,mpc_xind);	
	
  //6. MAKE PARTICLE LIST FOR EACH CELL IN DOMAIN._______________________________________________________________________________________________________________________________
    Cell ***cell;
    alloc3Darray<Cell***, Cell**, Cell*, Cell>(cell, xind);
  
    Particle_Lists<Cell***, Bead*>(cell, particles, partTot, xind, cellsize);			
  
  //7. FIND ACCELERATION OF EACH BEAD BASED ON IC, and uTot, and pressTot._______________________________________________________________________________________	
    if (drug_bonded == 0) {
      //Drug not attached
      Acceleration_uTot_FreeDrug_Acidic<Cell***, Bead*, std::list<Bead>::iterator, Bead>(cell, particles,
                                                                                         branchNum, branchmonNum, prevbranchmonTot, lipidbranchmonTot,
                                                                                         monNum, lipTot, branchTot, drugTot, solvTot, partTot, xind, uTot, pressTot, length,
                                                                                         mh, mp, mt, md, ms,
                                                                                         sigma, rc, rinf, kbond, kbend, kT,
	                                                                                     ehh, ehp, eht, ehd, ehs, 
							                                                             epp, ept, epd, eps,		     
							                                                             ett, etd, ets, 
							                                                             edd, eds,
							                                                             whs,
																					     wps, 
							                                                             wtt, wtd, 
							                                                             wdd);	
    }
    else {
      //Drug attached
      Acceleration_uTot_BondedDrug_Acidic<Cell***, Bead*, std::list<Bead>::iterator, Bead>(cell, particles,
                                                                                           branchNum, branchmonNum, prevbranchmonTot, lipidbranchmonTot,
                                                                                           monNum, lipTot, branchTot, drugTot, solvTot, partTot, xind, uTot, pressTot, length,
                                                                                           mh, mp, mt, md, ms,
                                                                                           sigma, rc, rinf, kbond, kbend, kT,
	                                                                                       ehh, ehp, eht, ehd, ehs, 
							                                                               epp, ept, epd, eps,		     
							                                                               ett, etd, ets, 
							                                                               edd, eds,
							                                                               whs,
																						   wps, 
							                                                               wtt, wtd, 
							                                                               wdd);	   
  }	
	
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
	    << "# " << "Timesteps: " << acidicstepTot << std::endl
	    << "# " << "xsf steps: " << xsfstep <<std::endl
        << "PRIMEVEC" << std::endl
        << length*Vector(1.0, 0.0, 0.0) << std::endl
        << length*Vector(0.0, 1.0, 0.0) << std::endl
	    << length*Vector(0.0, 0.0, 1.0) << std::endl
	    << "ANIMSTEPS " << 1 + floor(acidicstepTot/xsfstep) << std::endl
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
	
  //10. STORE DATA IN ARRAYS.____________________________________________________________________________________________________
    Vector *xsf_data = new Vector[(acidicstepTot/xsfstep)*partTot];
    int *int_data = new int[acidicstepTot*intdataNum];	
    double *double_data = new double[acidicstepTot*doubledataNum];
    //mpc________________________________________________________
    int *mpc_int_data = new int[acidicstepTot/tau];    
	double *mpc_double_data = new double[2*acidicstepTot/tau];
	
  //11. TIME EVOLUTION FROM IC TO EQUILIBRIUM._________________________________________________________
    for (step = 1; step <= acidicstepTot; step++) {
      //SRD___________________________________________________________________________________________________________________________________________________________
      if ((step%tau) == 0) {
	    Grid_Shift_Assign_vrand<Bead*>(particles, solvTot, partTot, mpc_xind, mpc_cellsize, length, kT, ms);
	    Calc_vmean<MPC_Cell***, Bead*>(mpc_cell, particles, solvTot, partTot, mpc_xind, MPCcellNum_with_no_nonsolv, solvNum_avg_MPCcell_with_no_nonsolv, solvNum_var_MPCcell_with_no_nonsolv, mpc_cellsize);	
	    Update_Solvent_Velocities<MPC_Cell***, Bead*>(mpc_cell, particles, solvTot, partTot);
	    Return_Positions<Bead*>(particles, solvTot, partTot);
	  }
	  std::cout << step <<std::endl;
  
	  //1. For full timestep update positions r(time) to r(time + dt) via TS expansion.
      int i;
	  #pragma omp parallel default(none) shared(partTot, particles, dt, length,t0,m,epsilon,sigma) private(i) 
	 {
	    #pragma omp for
	    for (i = 0; i < partTot; i++) {
          particles[i].r += (particles[i].v)*dt  + 0.5*(particles[i].a)*dt*dt;  
          Boundary_Conditions<Bead*>(particles, i, length);
        }
      }
	
	  //Cluster analysis after new positions.________________________________________																							
      Cluster_Drug_Analysis<Bead,Bead*,Cluster,std::list<Cluster>,std::list<Cluster>::iterator>(particles, 
	                                                                                            drug_bonded, drug_released,
																							    mon1,mon2,mon3,mon4,monNum,lipTot,lipidbranchmonTot,branchTot,drugTot,length,
																							    rclust,clustNum,clustLipNum_avg,clustLipNum_var,clustRad_avg,clustRad_var,micRad_avg,micRad_var, 
																								drugIn,drugOut,clustDrugNum_avg,clustDrugNum_var,drugInRad_div_clustRad_avg,drugInRad_div_clustRad_var);
																								
	  //2. For half timestep update velocities v(time) to v(time + 0.5*dt) via TS expansion.
	  #pragma omp parallel default(none) shared(partTot, particles, dt) private(i) 
	  {
	    #pragma omp for
        for (i = 0; i < partTot; i++) {
          particles[i].v += 0.5*(particles[i].a)*dt;
	    }
      }
	
	  //3. For full timestep update accelerations a(time) to a(time + dt) via a=f/m
	  Particle_Lists<Cell***, Bead*>(cell, particles, partTot, xind, cellsize);
	  
      if (drug_bonded == 0) {
        //Drug not attached
        Acceleration_uTot_FreeDrug_Acidic<Cell***, Bead*, std::list<Bead>::iterator, Bead>(cell, particles,
                                                                                           branchNum, branchmonNum, prevbranchmonTot, lipidbranchmonTot,
                                                                                           monNum, lipTot, branchTot, drugTot, solvTot, partTot, xind, uTot, pressTot, length,
                                                                                           mh, mp, mt, md, ms,
                                                                                           sigma, rc, rinf, kbond, kbend, kT,
	                                                                                       ehh, ehp, eht, ehd, ehs, 
							                                                               epp, ept, epd, eps,		     
							                                                               ett, etd, ets, 
							                                                               edd, eds,
							                                                               whs,
																					       wps, 
							                                                               wtt, wtd, 
							                                                               wdd);	
      }
      else {
        //Drug attached
        Acceleration_uTot_BondedDrug_Acidic<Cell***, Bead*, std::list<Bead>::iterator, Bead>(cell, particles,
                                                                                             branchNum, branchmonNum, prevbranchmonTot, lipidbranchmonTot,
                                                                                             monNum, lipTot, branchTot, drugTot, solvTot, partTot, xind, uTot, pressTot, length,
                                                                                             mh, mp, mt, md, ms,
                                                                                             sigma, rc, rinf, kbond, kbend, kT,
	                                                                                         ehh, ehp, eht, ehd, ehs, 
							                                                                 epp, ept, epd, eps,		     
							                                                                 ett, etd, ets, 
							                                                                 edd, eds,
							                                                                 whs,
																						     wps, 
							                                                                 wtt, wtd, 
							                                                                 wdd);	   
      }	
	
	  //4. For Full timestep update velocities v(time) to v(time + dt) via TS expansion.
      //#pragma omp parallel default(none) shared(partTot, particles, dt) private(i) 
	  {
	    #pragma omp for
	    for (i = 0; i < partTot; i++) {
          particles[i].v += 0.5*(particles[i].a)*dt;       
        }
	  }
	
	  //Calculate KE after new speeds.____________________________________________________________________________
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
	  
    }
  
  //12. WRITE FILES.__________________________________________________________________________________________________________________
    //Write xsf file_______________________________________________
    xsf.open(xsf_file, std::ofstream::app);
    for (int i = 1; i <= acidicstepTot/xsfstep; i++) {
      xsf << "ATOMS " << 1 + i << std::endl;
	  int xsfstepcount = partTot*(i - 1);
	  for (int j = 0; j < partTot; j++) {
        xsf << j + 1 << " " << xsf_data[xsfstepcount + j] << std::endl;
      }  
    }	
    xsf.close();
    delete [] xsf_data;
 
    //Write int binary data file_________________________________________
	int_bin.open(int_bin_file, std::ofstream::app|std::ofstream::binary);
    for (int i = 0; i < acidicstepTot*intdataNum; i++) { 
      int_bin.write((char*) &int_data[i], sizeof(int));
    }
    delete [] int_data;
    for (int i = 0; i < acidicstepTot/tau; i++) { 
      int_bin.write((char*) &mpc_int_data[i], sizeof(double));
    }
	int_bin.close();
    delete [] mpc_int_data;
	
	//Write double binary data file___________________________________________
	double_bin.open(double_bin_file, std::ofstream::app|std::ofstream::binary);
	for (int i = 0; i < acidicstepTot*doubledataNum; i++) {
      double_bin.write((char*) &double_data[i], sizeof(double));
    }
    delete [] double_data;
    for (int i = 0; i < 2*acidicstepTot/tau; i++) {
      double_bin.write((char*) &mpc_double_data[i], sizeof(double));
    }
    double_bin.close();
	delete [] mpc_double_data;

  //13. DEALLOCATE MEM FOR REMAINIGN ARRAYS.____________________________________________________________
    delete[] memblock;
    delete[] state_double_values;
    dealloc3Darray<MPC_Cell***>(mpc_cell, mpc_xind);
    dealloc3Darray<Cell***>(cell, xind);  //For 3D array cell of class Cell
    delete [] particles;
 
  return 0;	
	
}	
	
