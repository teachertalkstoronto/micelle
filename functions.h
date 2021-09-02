#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <omp.h>
#include "vecmat3.h"	//For vector operations
#include "myrandom.h"
#include <cmath>		//For math operations
#include <list>			//For PL and NL
#include <vector>
#include <sstream>
#include <string>



std::string Num_to_String(double num) {
  std::string stringnum;
  std::ostringstream convert;
  convert << num;
  stringnum = convert.str();
  return stringnum;
}

std::string Num_to_String(int num) {
  std::string stringnum;
  std::ostringstream convert;
  convert << num;
  stringnum = convert.str();
  return stringnum;
}


// For structure H - pHse - T
template <typename T1>
void Assign_Type_Lipid_Monomer_Indices_1(T1 &p,
                                      vecmat3::Vector<int> &btmonTot, int &lbmonTot,
                                      int &hmonNum, int &phsmonNum, int &tmonNum, int &mNum, int &ltot, int &btot, int &dtot, int &stot, int &ptot) {
    for (int i = 0; i < ptot; i++) {
	  p[i].ind = i; 
	}
	
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < hmonNum; j++) {
	    p[i*(mNum + lbmonTot) + j].ptype = 1;
		p[i*(mNum + lbmonTot) + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + j].monomer = j + 1;
	  }
	}
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < phsmonNum; j++) {
	    p[i*(mNum + lbmonTot) + hmonNum + j].ptype = 2;
		p[i*(mNum + lbmonTot) + hmonNum + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + hmonNum + j].monomer = hmonNum + j + 1;
	  }
	}
    for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < tmonNum; j++) {
	    p[i*(mNum + lbmonTot) + hmonNum + phsmonNum + j].ptype = 3;
	    p[i*(mNum + lbmonTot) + hmonNum + phsmonNum + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + hmonNum + phsmonNum + j].monomer = hmonNum + phsmonNum + j + 1;
	  }
	}
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < btmonTot(0); j++) {
	    p[mNum + i*(mNum + lbmonTot) + j].ptype = 1;
		p[mNum + i*(mNum + lbmonTot) + j].lipid = i + 1;
	  }
	}
    for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < btmonTot(1); j++) {
	    p[mNum + i*(mNum + lbmonTot) + btmonTot(0) + j].ptype = 2;
		p[mNum + i*(mNum + lbmonTot) + btmonTot(0) + j].lipid = i + 1;
	  }
	}
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < btmonTot(2); j++) {
	    p[mNum + i*(mNum + lbmonTot) + btmonTot(0) + btmonTot(1) + j].ptype = 3;
		p[mNum + i*(mNum + lbmonTot) + btmonTot(0) + btmonTot(1) + j].lipid = i + 1;
	  }
	}
	for (int i = mNum*ltot + btot; i < mNum*ltot + btot + dtot; i++) {
	  p[i].ptype = 4;
	  p[i].lipid = -p[i].ind;
	}
	for (int i = mNum*ltot + btot + dtot; i < ptot; i++) {
	  p[i].ptype = 5;
	  p[i].lipid = -p[i].ind;
	}
}


// For structure H - T - pHse
//Assign_Type_Lipid_Monomer_Indices_2<Bead*>(particles,branchtypemonTot, lipidbranchmonTot, headmonNum, tailmonNum, pHsemonNum, monNum, lipTot, branchTot, drugTot, acidTot, solvTot, partTot) {
template <typename T1>
void Assign_Type_Lipid_Monomer_Indices_2(T1 &p,
                                         vecmat3::Vector<int> &btmonTot, int &lbmonTot,
                                         int &hmonNum, int &tmonNum, int &phsmonNum, int &mNum, int &ltot, int &btot, int &dtot, int &stot, int &ptot) {
    for (int i = 0; i < ptot; i++) {
	  p[i].ind = i; 
	}
	
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < hmonNum; j++) {
	    p[i*(mNum + lbmonTot) + j].ptype = 1;
		p[i*(mNum + lbmonTot) + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + j].monomer = j + 1;
	  }
	}
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < tmonNum; j++) {
	    p[i*(mNum + lbmonTot) + hmonNum + j].ptype = 3;
		p[i*(mNum + lbmonTot) + hmonNum + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + hmonNum + j].monomer = hmonNum + j + 1;
	  }
	}
    for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < phsmonNum; j++) {
	    p[i*(mNum + lbmonTot) + hmonNum + tmonNum + j].ptype = 2;
	    p[i*(mNum + lbmonTot) + hmonNum + tmonNum + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + hmonNum + tmonNum + j].monomer = hmonNum + tmonNum + j + 1;                                     //Made them all protonated
	  }
	}
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < btmonTot(0); j++) {
	    p[i*(mNum + lbmonTot) + mNum + j].ptype = 1;
		p[i*(mNum + lbmonTot) + mNum + j].lipid = i + 1;
	  }
	}
    for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < btmonTot(1); j++) {
	    p[mNum + i*(mNum + lbmonTot) + btmonTot(0) + j].ptype = 3;
		p[mNum + i*(mNum + lbmonTot) + btmonTot(0) + j].lipid = i + 1;
	  }
	}
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < btmonTot(2); j++) {
	    p[i*(mNum + lbmonTot) + mNum + btmonTot(0) + btmonTot(1) + j].ptype = 2;
		p[i*(mNum + lbmonTot) + mNum + btmonTot(0) + btmonTot(1) + j].lipid = i + 1;
	  }
	}
	for (int i = mNum*ltot + btot; i < mNum*ltot + btot + dtot; i++) {
	  p[i].ptype = 4;
	  p[i].lipid = -p[i].ind;
	}
	for (int i = mNum*ltot + btot + dtot; i < ptot; i++) {
	  p[i].ptype = 5;
	  p[i].lipid = -p[i].ind;
	}
}









// For structure H - pHse - T
template <typename T1>
void Assign_Type_Lipid_Monomer_Indices_3(T1 &p,
                                      vecmat3::Vector<int> &btmonTot, int &lbmonTot,
                                      int &hmonNum, int &phsmonNum, int &tmonNum, int &mNum, int &ltot, int &btot, int &dtot, int &stot, int &ptot) {
    for (int i = 0; i < ptot; i++) {
	  p[i].ind = i; 
	}
	
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < hmonNum; j++) {
	    p[i*(mNum + lbmonTot) + j].ptype = 1;
		p[i*(mNum + lbmonTot) + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + j].monomer = j + 1;
	  }
	}
	int counter = 0;
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < phsmonNum; j++) {
	    p[i*(mNum + lbmonTot) + hmonNum + j].ptype = 2;
		p[i*(mNum + lbmonTot) + hmonNum + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + hmonNum + j].monomer = hmonNum + j + 1;
	    p[i*(mNum + lbmonTot) + hmonNum + j].drug = p[i*(mNum + lbmonTot) + hmonNum + j].ind;
		p[mNum*ltot + btot + counter].drug = p[i*(mNum + lbmonTot) + hmonNum + j].ind;
		counter++;
	  }
	}
    for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < tmonNum; j++) {
	    p[i*(mNum + lbmonTot) + hmonNum + phsmonNum + j].ptype = 3;
	    p[i*(mNum + lbmonTot) + hmonNum + phsmonNum + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + hmonNum + phsmonNum + j].monomer = hmonNum + phsmonNum + j + 1;
	  }
	}
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < btmonTot(0); j++) {
	    p[mNum + i*(mNum + lbmonTot) + j].ptype = 1;
		p[mNum + i*(mNum + lbmonTot) + j].lipid = i + 1;
	  }
	}
    for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < btmonTot(1); j++) {
	    p[i*(mNum + lbmonTot) + mNum + btmonTot(0) + j].ptype = 2;
		p[i*(mNum + lbmonTot) + mNum + btmonTot(0) + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + mNum + btmonTot(0) + j].drug = p[mNum + i*(mNum + lbmonTot) + btmonTot(0) + j].ind;         //part of no proton
		p[mNum*ltot + btot + counter].drug = p[mNum + i*(mNum + lbmonTot) + btmonTot(0) + j].ind;
		counter++;
	  }
	}
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < btmonTot(2); j++) {
	    p[i*(mNum + lbmonTot) + mNum + btmonTot(0) + btmonTot(1) + j].ptype = 3;
		p[i*(mNum + lbmonTot) + mNum + btmonTot(0) + btmonTot(1) + j].lipid = i + 1;
	  }
	}
	for (int i = mNum*ltot + btot; i < mNum*ltot + btot + dtot; i++) {
	  p[i].ptype = 4;
	  p[i].lipid = -p[i].ind;
	}
	for (int i = mNum*ltot + btot + dtot; i < ptot; i++) {
	  p[i].ptype = 5;
	  p[i].lipid = -p[i].ind;
	}
}




// For structure H - T - pHse
//Assign_Type_Lipid_Monomer_Indices_2<Bead*>(particles,branchtypemonTot, lipidbranchmonTot, headmonNum, tailmonNum, pHsemonNum, monNum, lipTot, branchTot, drugTot, acidTot, solvTot, partTot) {
template <typename T1>
void Assign_Type_Lipid_Monomer_Indices_4(T1 &p,
                                         vecmat3::Vector<int> &btmonTot, int &lbmonTot,
                                         int &hmonNum, int &tmonNum, int &phsmonNum, int &mNum, int &ltot, int &btot, int &dtot, int &stot, int &ptot) {
    for (int i = 0; i < ptot; i++) {
	  p[i].ind = i; 
	}
	
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < hmonNum; j++) {
	    p[i*(mNum + lbmonTot) + j].ptype = 1;
		p[i*(mNum + lbmonTot) + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + j].monomer = j + 1;
	  }
	}
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < tmonNum; j++) {
	    p[i*(mNum + lbmonTot) + hmonNum + j].ptype = 3;
		p[i*(mNum + lbmonTot) + hmonNum + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + hmonNum + j].monomer = hmonNum + j + 1;
	  }
	}
	int counter = 0;
    for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < phsmonNum; j++) {
	    p[i*(mNum + lbmonTot) + hmonNum + tmonNum + j].ptype = 2;
	    p[i*(mNum + lbmonTot) + hmonNum + tmonNum + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + hmonNum + tmonNum + j].monomer = hmonNum + tmonNum + j + 1;
		p[i*(mNum + lbmonTot) + hmonNum + tmonNum + j].drug = p[i*(mNum + lbmonTot) + hmonNum + tmonNum + j].ind;
		p[mNum*ltot + btot + counter].drug = p[i*(mNum + lbmonTot) + hmonNum + tmonNum + j].ind;
		counter++;
	  }
	}
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < btmonTot(0); j++) {
	    p[i*(mNum + lbmonTot) + mNum + j].ptype = 1;
		p[i*(mNum + lbmonTot) + mNum + j].lipid = i + 1;
	  }
	}
    for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < btmonTot(1); j++) {
	    p[mNum + i*(mNum + lbmonTot) + btmonTot(0) + j].ptype = 3;
		p[mNum + i*(mNum + lbmonTot) + btmonTot(0) + j].lipid = i + 1;
	  }
	}
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < btmonTot(2); j++) {
	    p[i*(mNum + lbmonTot) + mNum + btmonTot(0) + btmonTot(1) + j].ptype = 2;
		p[i*(mNum + lbmonTot) + mNum + btmonTot(0) + btmonTot(1) + j].lipid = i + 1;
		p[i*(mNum + lbmonTot) + mNum + btmonTot(0) + btmonTot(1) + j].drug = p[i*(mNum + lbmonTot) + mNum + btmonTot(0) + btmonTot(1) + j].ind;
		p[mNum*ltot + btot + counter].drug = p[i*(mNum + lbmonTot) + mNum + btmonTot(0) + btmonTot(1) + j].ind;	
		counter++;
	  }
	}
	for (int i = mNum*ltot + btot; i < mNum*ltot + btot + dtot; i++) {
	  p[i].ptype = 4;
	  p[i].lipid = -p[i].ind;
	}
	for (int i = mNum*ltot + btot + dtot; i < ptot; i++) {
	  p[i].ptype = 5;
	  p[i].lipid = -p[i].ind;
	}
}





template<typename T> 
void Initial_Velocities(T &p, int &ptot, double &m1, double &m2, double &m3, double &m4, double &m5, double &temp) {
  ENG eng;
  eng.seed(static_cast<unsigned int>(std::time(0)));
  
  N_DIST h_dist(0.0, sqrt(temp*(1.0/m1)));
  N_GEN h_gen(eng, h_dist);
  
  N_DIST phs_dist(0.0, sqrt(temp*(1.0/m2)));
  N_GEN phs_gen(eng, phs_dist);
  
  N_DIST t_dist(0.0, sqrt(temp*(1.0/m3)));
  N_GEN t_gen(eng, t_dist);
  
  N_DIST d_dist(0.0, sqrt(temp*(1.0/m4)));
  N_GEN d_gen(eng, d_dist);
  
  N_DIST s_dist(0.0, sqrt(temp*(1.0/m5)));
  N_GEN s_gen(eng, s_dist);
  
  for (int i = 0; i < ptot; i++) {
    //For head beads:
    if ( p[i].ptype == 1 )
	  (p[i]).v = Vector(h_gen(), h_gen(), h_gen());
	//For pHse beads:
	else if ( p[i].ptype == 2 ) 
	  p[i].v = Vector(phs_gen(), phs_gen(), phs_gen());    
    //For tail beads:	
	else if ( p[i].ptype == 3 ) 
	  p[i].v = Vector(t_gen(), t_gen(), t_gen());	
	//For drug beads:
    else if ( p[i].ptype == 4 ) 
	  p[i].v = Vector(d_gen(), d_gen(), d_gen());
	//For solvent beads:
	else if ( p[i].ptype == 5 ) 
	  p[i].v = Vector(s_gen(), s_gen(), s_gen()); 	  
  }
}



template<typename T>
void Max_Speed_Avg_Speed(T &p, int &stot, int &ptot, double &vavg, double &vavgs, double &vm, double &vms) {
  vavg = 0.0;
  vavgs = 0.0;
  vm = 0.0;
  vms = 0.0;
  for (int i = 0; i < ptot; i++) {
   vavg += p[i].v.nrm();
   if ( vm < p[i].v.nrm() ) 
     vm = p[i].v.nrm();
  }
  vavg *= (1.0/ptot);
  
  for (int i = ptot - stot; i < ptot; i++) {
   vavgs += p[i].v.nrm();
   if ( vms < p[i].v.nrm() ) 
     vms = p[i].v.nrm();
  }
  vavgs *= (1.0/stot);
  
}

template<typename T>
void KineticEnergy_LinearMomentum(T &p, int &stot, int &ptot, double &ektot, Vector &lin_momtot, double &m1, double &m2, double &m3, double &m4, double &m5) {
  int i;
  ektot = 0.0;
  lin_momtot = Vector(0.0, 0.0, 0.0);
  //#pragma omp parallel default(none) shared(p, ptot, m1, m2, m3, m4) private(i) reduction(+:ek,linmom) 
  //{
    //#pragma omp for
    for (i = 0; i < ptot - stot; i++) {
	  if ( p[i].ptype == 1 ) {
	    ektot += 0.5*m1*(p[i].v.nrm2());
	    lin_momtot += m1*p[i].v;
	  }
	  else if ( p[i].ptype == 2 ) {
	    ektot += 0.5*m2*(p[i].v.nrm2());
	    lin_momtot += m2*p[i].v;
	  }
      else if ( p[i].ptype == 3 ) {
	    ektot += 0.5*m3*(p[i].v.nrm2());
	    lin_momtot += m3*p[i].v;
	  }
      else if ( p[i].ptype == 4 ) {
	    ektot += 0.5*m4*(p[i].v.nrm2());
	    lin_momtot += m4*p[i].v;
	  }
    }
	
	for (int i = ptot - stot; i < ptot; i++) {
	  ektot += 0.5*m5*(p[i].v.nrm2());
	  lin_momtot += m5*p[i].v;
	}
	
  //}
}


//Boundary Conditions in 1D
void BC_1D(Vector &a, int jj, double &len) {
  while ( a(jj) >= 0.5*len ) 
    a(jj) -= len;
  while ( a(jj) < -0.5*len )
    a(jj) += len;
}

//Boundary conditions in 3D
void BC(Vector &a, double &len) {
  BC_1D(a, 0, len);
  BC_1D(a, 1, len);
  BC_1D(a, 2, len);
}


//Boundary Conditions in 1D
template <typename T>
void Boundary_Conditions_1D(T &p, int ii, int jj, double &len) {
  while ( p[ii].r(jj) >= 0.5*len ) 
    p[ii].r(jj) -= len;
  while ( p[ii].r(jj) < -0.5*len )
    p[ii].r(jj) += len;
}

//Boundary conditions in 3D
template <typename T1>
void Boundary_Conditions(T1 &p, int iii, double &len) {

  Boundary_Conditions_1D<T1>(p, iii, 0, len);
  Boundary_Conditions_1D<T1>(p, iii, 1, len);
  Boundary_Conditions_1D<T1>(p, iii, 2, len);

  
  /*
  while ( p[iii].r(0) >= 0.5*len ) 
    p[iii].r(0) -= len;
  while ( p[iii].r(0) < -0.5*len )
    p[iii].r(0) += len;
	
  while ( p[iii].r(1) >= 0.5*len ) 
    p[iii].r(1) -= len;
  while ( p[iii].r(1) < -0.5*len )
    p[iii].r(1) += len;
	
  while ( p[iii].r(2) >= 0.5*len ) 
    p[iii].r(2) -= len;
  while ( p[iii].r(2) < -0.5*len )
    p[iii].r(2) += len;
	*/
}



//Give initial positions to lipid monomer beads that are at the end of the tail along xy plane.
template <typename T>
void FCC_Lattice_Initial_Positions(T &lp, double &lps, int &lpd, double &len) {
  int counter = 0;
  //top and bottom of fcc unit
  for (int i = 0; i < lpd; i++) {
    for (int j = 0; j < lpd; j++) {
      for (int k = 0; k < lpd; k++) {
	    (lp[counter]).r = Vector(-0.5*len + 0.75*lps + i*lps, -0.5*len + 0.75*lps + j*lps, -0.5*len + 0.25*lps + k*lps);
        counter++;
	  }
	}
  }
  for (int i = 0; i < lpd; i++) {
    for (int j = 0; j < lpd; j++) {
      for (int k = 0; k < lpd; k++) {
	    (lp[counter]).r = Vector(-0.5*len + 0.25*lps + i*lps, -0.5*len + 0.25*lps + j*lps, -0.5*len + 0.25*lps + k*lps);
        counter++;
	  }
	}
  }
  //middle layer of fcc unit
  for (int i = 0; i < lpd; i++) {
    for (int j = 0; j < lpd; j++) {
      for (int k = 0; k < lpd; k++) {
	    (lp[counter]).r = Vector(-0.5*len + 0.25*lps + i*lps, -0.5*len + 0.75*lps + j*lps, -0.5*len + 0.75*lps + k*lps);
        counter++;
	  }
	}
  }
  for (int i = 0; i < lpd; i++) {
    for (int j = 0; j < lpd; j++) {
      for (int k = 0; k < lpd; k++) {
	    (lp[counter]).r = Vector(-0.5*len + 0.75*lps + i*lps, -0.5*len + 0.25*lps + j*lps, -0.5*len + 0.75*lps + k*lps);
        counter++;
	  }
	}
  }
  
}



//For drug attached
template <typename T1, typename T2, size_t size_x, size_t size_y>
void Lipid_Initial_Positions(T1 &p, T2 &lp, int &drug_bonded,
							int (&bNum)[size_x], int (&bmonNum)[size_x][size_y], int (&pbmonTot)[size_x][size_y], int &lbmonTot, 
                            int &mNum, int &ltot, int &btot, int &dtot, int &lpsit, double &lps, double &len) {
  ENG eng;
  eng.seed(static_cast<unsigned int>(std::time(0)));
  
  U_INT_DIST dist1(0, lpsit - 1);
  U_INT_GEN gen1(eng, dist1);
  
  U_INT_DIST dist2(0, 1);
  U_INT_GEN gen2(eng, dist2);
  
  U_INT_DIST dist3(1, 3);
  U_INT_GEN gen3(eng, dist3);
  
  //Assign position for lipid beads one lipid at a time._________________________________________________________________________
  for (int i = 0; i < ltot; i++) {
    //Start positioning from first main chain monomer (and branched units) to last (and its branched units).
    for (int j = 0; j < mNum; j++) {
	  //If position first main chain monomer, you can pick any random lattice site.
	  if (j == 0) {
        int x = gen1();
        while ( lp[x].taken == 1 ) {
          x = gen1();
        }
		p[i*(mNum + lbmonTot)].r = (lp[x]).r;
		(lp[x]).taken = 1;
	  }
	  //If positioning following main chain monomers, random position depends on previous monomer's position due to bond.
	  else {
	    int neighbor = 0;
	    while (neighbor == 0) {
          int xgen = gen2();
	      if (xgen == 0) xgen = -1;
          int ygen = gen2();
	      if (ygen == 0) ygen = -1;
		  Vector shift(0.0, 0.0, 0.0);
          int plane = gen3();
		  if (plane == 1) {
			shift = Vector(0.0, xgen*0.5*lps, ygen*0.5*lps);
			p[i*(mNum + lbmonTot) + j].r = p[i*(mNum + lbmonTot) + j - 1].r + shift;
		    }
	      else if (plane == 2) { 
			shift = Vector(xgen*0.5*lps, 0.0, ygen*0.5*lps);
			p[i*(mNum + lbmonTot) + j].r = p[i*(mNum + lbmonTot) + j - 1].r + shift; 
		  }
	 	  else if (plane == 3) {
			shift = Vector(xgen*0.5*lps, ygen*0.5*lps, 0.0);
			p[i*(mNum + lbmonTot) + j].r = p[i*(mNum + lbmonTot) + j - 1].r + shift;
		  }
          Boundary_Conditions<T1>(p, i*(mNum + lbmonTot) + j, len);
          for (int q = 0; q < lpsit; q++) {
            if ( dist((lp[q]).r, p[i*(mNum + lbmonTot) + j].r) < 0.00005 ) {
		      if ( (lp[q]).taken == 0 ) {
                (lp[q]).taken = 1;
				//std:: cout << lp[q].r << std::endl;
			    neighbor = 1;
			    break;
			  }
	        }
	      }
	    }
	  }  
	  //After positioning main chain monomer, position the branched chain monomers.
	  for (int k = 0; k < bNum[j]; k++) {
	    //One branch at a time ladies.
	    for (int m = 0; m < bmonNum[j][k]; m++) {
		  int neighbor = 0;
	      while (neighbor == 0) {
            int xgen = gen2();
	        if (xgen == 0) xgen = -1;
            int ygen = gen2();
	        if (ygen == 0) ygen = -1;
			Vector shift(0.0, 0.0, 0.0);
            int plane = gen3();
	        if (m == 0) {
		      if (plane == 1) {
				shift = Vector(0.0, xgen*0.5*lps, ygen*0.5*lps);
				p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k]].r = p[i*(mNum + lbmonTot) + j].r + shift;
	          }
	          else if (plane == 2) {
				shift = Vector(xgen*0.5*lps, 0.0, ygen*0.5*lps);
				p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k]].r = p[i*(mNum + lbmonTot) + j].r + shift;   
	          }
	 	      else if (plane == 3) {
				shift = Vector(xgen*0.5*lps, ygen*0.5*lps, 0.0);
				p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k]].r = p[i*(mNum + lbmonTot) + j].r + shift;
	          }
              Boundary_Conditions<T1>(p, i*(mNum + lbmonTot) + mNum + pbmonTot[j][k], len);
              for (int q = 0; q < lpsit; q++) {
                if (dist((lp[q]).r, p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].r) < 0.00005) {
		          if ((lp[q]).taken == 0) {
				    (lp[q]).taken = 1;
			        neighbor = 1;
					break;
				  }
	            }
              }
	        }
		    else {
			  if (plane == 1) {
               
				shift = Vector(0.0, xgen*0.5*lps, ygen*0.5*lps);
				p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].r = p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m - 1].r + shift;
	          }
	          else if (plane == 2) {
				shift = Vector(xgen*0.5*lps, 0.0, ygen*0.5*lps);
				p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].r = p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m - 1].r + shift;   
	          }
	 	      else if (plane == 3) {
				shift = Vector(xgen*0.5*lps, ygen*0.5*lps, 0.0);
				p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].r = p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m - 1].r + shift;
	          }
              Boundary_Conditions<T1>(p, i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m, len);
              for (int q = 0; q < lpsit; q++) {
                if (dist((lp[q]).r, p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].r) < 0.00005) {
		          if ((lp[q]).taken == 0) {
					(lp[q]).taken = 1;
			        neighbor = 1;
					break;
				  }
                }
              }	
			}
		  }
	    }
	  }
    }
  }
  
  //If the drug is bonded then assign its positions too
  if (drug_bonded == 1) {
    for (int i = mNum*ltot + btot; i <  mNum*ltot + btot + dtot; i++) {
      int neighbor = 0;
	  while (neighbor == 0) {
	    int xgen = gen2();
	    if (xgen == 0) xgen = -1;
        int ygen = gen2();
	    if (ygen == 0) ygen = -1;
	    Vector shift(0.0, 0.0, 0.0);
        int plane = gen3();
	    if (plane == 1) {
	      shift = Vector(0.0, xgen*0.5*lps, ygen*0.5*lps);
		  p[i].r = p[p[i].drug].r + shift;
	    }
	    else if (plane == 2) {
          shift = Vector(xgen*0.5*lps, 0.0, ygen*0.5*lps); 		
		  p[i].r = p[p[i].drug].r + shift;  
        }		
        else if (plane == 3) {
		  shift = Vector(xgen*0.5*lps, ygen*0.5*lps, 0.0);
		  p[i].r = p[p[i].drug].r + shift;
	    }
        Boundary_Conditions<T1>(p, i, len);
        for (int q = 0; q < lpsit; q++) {
          if (dist((lp[q]).r, p[i].r) < 0.00005) {
		    if ((lp[q]).taken == 0) {
		      (lp[q]).taken = 1;
			  neighbor = 1;
			  break;
		    }
	      }
        }
	  }
    }
  }
 
  
}

//For drug attached
template <typename T1>
void NonLipid_Initial_Positions_2(T1 &p,int &drug_bonded,int &mNum,int &ltot,int &btot,int &dtot,int &stot,int &ptot,int &lpsit,double &lps,double &len) {
								  
  ENG eng;
  eng.seed(static_cast<unsigned int>(std::time(0)));
  
  U_INT_DIST dist1(0, lpsit - 1);
  U_INT_GEN gen1(eng, dist1);
  
  U_INT_DIST dist2(0, 1);
  U_INT_GEN gen2(eng, dist2);
  
  U_INT_DIST dist3(1, 3);
  U_INT_GEN gen3(eng, dist3);		

  U_DIST dist4(-0.5*len, 0.5*len);
  U_GEN gen4(eng, dist4);  
  
  //If the drug isn't bonded to the lipid. Make sure everything else is at least one lattice point away from a lipid lattice point.
  if (drug_bonded == 0) {
	//Position nonbonded rug beads at least one unit away from lipid beads (randomly).
    for (int i = mNum*ltot + btot; i < mNum*ltot + btot + dtot; i++) {
	  int stop = 0;
	  while (stop == 0) {
	    p[i].r = Vector(gen4(), gen4(), gen4());
	    int tooclose = 0;
	    for (int j = 0; j < i; j++) {
	      Vector pj_r = p[j].r;
	      if ( abs(p[i].r(0) - p[j].r(0)) > 0.5*len ) {
	         if ( p[j].r(0) < 0 ) pj_r(0) += len;
	         else if ( p[j].r(0) > 0 ) pj_r(0) -= len;
          }	
	      if ( abs(p[i].r(1) - p[j].r(1)) > 0.5*len ) {
	        if ( p[j].r(1) < 0 ) pj_r(1)  += len;
	        else if ( p[j].r(1) > 0 ) pj_r(1) -= len;
          }	
	      if ( abs(p[i].r(2) - p[j].r(2)) > 0.5*len ) {
	        if ( p[j].r(2) < 0 ) pj_r(2) += len;
	        else if ( p[j].r(2) > 0 ) pj_r(2) -= len;
          }
		  if (dist(p[i].r, pj_r) < 1.5*lps) {
		    tooclose = 1;
		    break;
		  }
	    }
		if (tooclose == 0) stop = 1;
	  }	  
	}
  }

  //Position solvent beads randomly.
  for (int i = ptot - stot; i < ptot; i++) {
	int stop = 0;
	while (stop == 0) {
	  p[i].r = Vector(gen4(), gen4(), gen4());
	  int tooclose = 0;
	  for (int j = 0; j < mNum*ltot + btot + dtot; j++) {
	    Vector pj_r = p[j].r;
	    if ( abs(p[i].r(0) - p[j].r(0)) > 0.5*len ) {
	       if ( p[j].r(0) < 0 ) pj_r(0) += len;
	       else if ( p[j].r(0) > 0 ) pj_r(0) -= len;
        }	
        if ( abs(p[i].r(1) - p[j].r(1)) > 0.5*len ) {
	      if ( p[j].r(1) < 0 ) pj_r(1)  += len;
	      else if ( p[j].r(1) > 0 ) pj_r(1) -= len;
        }	
	    if ( abs(p[i].r(2) - p[j].r(2)) > 0.5*len ) {
	      if ( p[j].r(2) < 0 ) pj_r(2) += len;
	      else if ( p[j].r(2) > 0 ) pj_r(2) -= len;
        }
		if (dist(p[i].r, pj_r) < 1.5*lps) {
		  tooclose = 1;
		  break;
		}
	  }
      if (tooclose == 0) stop = 1;
	} 
  }
  
}

//______________________________________________________________________________________________________________________________________________________________________

void gcpstns(Vector &b_r, Vector &a, Vector &b, double len) {	
    if ( std::abs(a(0) - b(0)) > 0.5*len ) {
	  if ( b(0) < 0 )
	    b_r(0) += len;
	  else if ( b(0) > 0 )
	    b_r(0) -= len;
    }	
	if ( std::abs(a(1) - b(1)) > 0.5*len ) {
	  if ( b(1) < 0 )
	    b_r(1) += len;
	  else if ( b(1) > 0 )
	    b_r(1) -= len;
    }	
	if ( std::abs(a(2) - b(2)) > 0.5*len ) {
	  if ( b(2) < 0 )
	    b_r(2) += len;
	  else if ( b(2) > 0 )
		b_r(2) -= len;
    }
}

template<typename T1, typename T2, typename T5, typename T6, typename T7>
void Cluster_Drug_Analysis(T2 &p, 
                           int &drug_bonded, int &drug_released, 
                           int &mon1, int &mon2, int &mon3, int &mon4, int &mNum, int &ltot, int lbmonTot, int &btot, int &dtot, double &len,
						   double &rclust, int &cNum, double &cLipNum_avg, double &cLipNum_var, double &cRad_avg, double &cRad_var, double &mRad_avg, double &mRad_var, 
						   int &dIn, int &dOut, double &cDrugNum_avg, double &cDrugNum_var, double &dInRad_div_cRad_avg, double &dInRad_div_cRad_var) 
						  {
						   
  //List of the lipid clusters in the domain__________________________________________________________
  T6 clustlist;
  clustlist.clear();
  
  //First say that all lipids and drugs are not part of a cluster____________________________________
  for (int i = 0; i < mNum*ltot + btot + dtot; i++) {
    p[i].cluster = -1;
	p[i].counted = 0;
  }

 //Assign lipids to cluster________________________________________________________________
  int cluster = 0;
  //For each lipid
  for (int i = 0; i < (ltot - 1); i++) {
    //If the lipid has not been checked to belong to a cluster
    if (p[i*(mNum + lbmonTot)].cluster == -1) {
	  //All monomers in lipid labelled to belong to "cluster" and the size of cluster is 1.
	  int clustLipNum = 1;
	  cluster++;
	  for (int j = 0; j < mNum + lbmonTot; j++) {
	    p[i*(mNum + lbmonTot) + j].cluster = cluster;
		p[i*(mNum + lbmonTot) + j].counted = 1;
      }
	  //Potential list of first monomers start with this lipid's first monomer.
	  std::list<Vector> headlist;
	  headlist.push_back(p[i*(mNum + lbmonTot)].r);
	  //std::cout << "Cluster: " << cluster << std::endl;
	  //std::cout << p[i*(mNum + lbmonTot)].r << "                 " << p[i*(mNum + lbmonTot)].ind << std::endl;
	  //For all lipids
	  for (int j = 0; j < ltot; j++) {
	    //If the lipid has not been checked to belong to a cluster
	    if (p[j*(mNum + lbmonTot)].cluster == -1) {
		  int mstop = 0;
		  for (int m = mon1; m < mon2 && mstop == 0; m++) {
	        int nstop = 0;
		    for (int n = mon1; n < mon2 && nstop == 0; n++) {
		  	  Vector pj_r = p[j*(mNum + lbmonTot) + n].r;
              gcpstns(pj_r, p[i*(mNum + lbmonTot) + m].r, p[j*(mNum + lbmonTot) + n].r, len);
		      if (dist(p[i*(mNum + lbmonTot) + m].r, pj_r) < 2.6) {
			    headlist.push_back(p[j*(mNum + lbmonTot)].r);
				//std::cout << p[j*(mNum + lbmonTot)].r << "                 " << p[j*(mNum + lbmonTot)].ind << std::endl;
			    for (int k = 0; k < mNum + lbmonTot; k++) {
	              p[j*(mNum + lbmonTot) + k].cluster = cluster;
                }
		        clustLipNum++;
			    nstop = 1;
			    mstop = 1;
			  } 
		    }	
	      }
	    }
	  }
	  
	  //Now so, if that oen lipid found another lipid._______________________________________________________________
	  if (clustLipNum == 1) {
	    for (int j = 0; j < mNum + lbmonTot; j++) {
	      p[i*(mNum + lbmonTot) + j].cluster = -2;
        }
	    cluster--;
      }
      else if (clustLipNum > 1 ) {
	    int stop = 0;
		while (stop == 0) {
		  int count = 0;
		  //For all lipids
  	      for (int j = 0; j < ltot; j++) {
		    //If the cluster is the same, but it has not been counted 
	        if (p[j*(mNum + lbmonTot)].cluster == cluster && p[j*(mNum + lbmonTot)].counted == 0) {
			  for (int k = 0; j < mNum + lbmonTot; j++) {
		        p[j*(mNum + lbmonTot) + k].counted = 1;
              }
		      for (int k = 0; k < ltot; k++) {
		        if (p[k*(mNum + lbmonTot)].cluster == -1) {
		          int mstop = 0;
		          for (int m = mon1; m < mon2 && mstop == 0; m++) {
	                int nstop = 0;
		            for (int n = mon1; n < mon2 && nstop == 0; n++) {
		  	          Vector pk_r = p[k*(mNum + lbmonTot) + n].r;
                      gcpstns(pk_r, p[j*(mNum + lbmonTot) + m].r, p[k*(mNum + lbmonTot) + n].r, len);
		              if (dist(p[j*(mNum + lbmonTot) + m].r, pk_r) < 2.6) { 
					  	headlist.push_back(p[k*(mNum + lbmonTot)].r);
						//std::cout << p[k*(mNum + lbmonTot)].r << "                 " << p[k*(mNum + lbmonTot)].ind << std::endl;
			            for (int q = 0; q < mNum + lbmonTot; q++) {
	                      p[k*(mNum + lbmonTot) + q].cluster = cluster;
                        }
					    count++;
		                clustLipNum++;
				        nstop = 1;
				        mstop = 1;
				      } 
					}
			      }
			    }
		      }	
	        }
		  }
		  if (count == 0) {
		    T5 clu;
	        clu.clustLipNum = clustLipNum;
	        clu.cluster = cluster;
			
/*
			std::cout << "Cluster: " << cluster << std::endl;
	        for (std::list<Vector>::iterator its =  headlist.begin(); its != headlist.end(); ++its) {
             std::cout << *its << std::endl;
	        }
			*/
			
			std::list<Vector> newheadlist;
			Vector reference = headlist.front();
			//std::cout << reference << std::endl;
			newheadlist.push_back(reference);
			
			
            int countheads = 0;
            for (std::list<Vector>::iterator its = headlist.begin(); its != headlist.end(); ++its) {
              if (countheads != 0) {
			    Vector pn_r = *its;
		        gcpstns(pn_r, reference, *its, len); 
				newheadlist.push_back(pn_r);
				//std::cout << pn_r << std::endl;
			  }
			  countheads++;
			}
			clu.firstmons = newheadlist;
			
			/*
			std::cout << "Cluster: " << cluster << std::endl;
	        for (std::list<Vector>::iterator its =  newheadlist.begin(); its != newheadlist.end(); ++its) {
             std::cout << *its << std::endl;
	        }
            */
			
            clustlist.push_back(clu); 
		    stop = 1;
		  }
		}
	  }
	}
  }
		      
  //Calc number of clusters____________________________________________________________________________________________
  cNum = clustlist.size();

  //Drug in cluster, drug outside of cluster______________________________________________________________
  
  //1. No clusters_____________________________________________________________________________________________________
  if (cNum == 0) {
    cLipNum_avg = 0.0;
    cLipNum_var = 0.0;
    cRad_avg = 0.0;
    cRad_var = 0.0;
    mRad_avg = 0.0;
    mRad_var = 0.0;

	//If drug is not bonded to lipid.
    if (drug_bonded == 0) {
   	  dIn = 0;
      dOut = dtot;
    }
	//Else drug is bonded to lipid.
    else {
	  //If drug is not released (neutral conditions).
	   if (drug_released == 0) {
	   	 dIn = dtot;
         dOut = 0;
	   }
	   //Else drug is released (acidic conditions). 
	   else {
	     dIn = 0;
		 dOut = dtot;
	   }  
    }

    cDrugNum_avg = 0.0; 
    cDrugNum_var = 0.0;
    dInRad_div_cRad_avg = 0.0;
    dInRad_div_cRad_var = 0.0;
  }
  //2. One cluster aka the micelle_____________________________________________________________________________________
  else if (cNum == 1) {
    for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
      cLipNum_avg = (*it).clustLipNum;
	}
    cLipNum_var = 0.0;
	
	//Find clustRad_avg
	cRad_avg = cRad_var = 0.0;
	for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
      (*it).centre = Vector(0.0,0.0,0.0);
	  for (std::list<Vector>::iterator its = (*it).firstmons.begin(); its != (*it).firstmons.end(); ++its) {
          (*it).centre += *its;
	  }
	  (*it).centre *= 1.0/(*it).firstmons.size();  
	  (*it).radius = 0.0;
	  for (std::list<Vector>::iterator its = (*it).firstmons.begin(); its != (*it).firstmons.end(); ++its) {
        (*it).radius += dist((*it).centre, *its);
	  }
	  (*it).radius *= (1.0/(*it).firstmons.size()); 
	  
	  cRad_avg += (*it).radius;	  
	}
	
    //Find micRad_avg
	mRad_avg = cRad_avg;
	//Find micRad_var
    mRad_var = 0.0; 
	for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
	  for (std::list<Vector>::iterator its = (*it).firstmons.begin(); its != (*it).firstmons.end(); ++its) {
        mRad_var += (dist((*it).centre, *its) - mRad_avg)*(dist((*it).centre, *its) - mRad_avg);
	  }
	  mRad_var += (*it).clustLipNum;
	}
	
    //NEWWWWWWWWWWWWWWWWWWWW
	dIn = dOut = 0;
    if (dtot > 0) {
	  //If drug is not bonded to lipid.
	  if (drug_bonded == 0) {
	    for (int i = mNum*ltot + btot; i < mNum*ltot + btot + dtot; i++) {
	    int foundclust = 0;
	    for (int j = 0; j < ltot && foundclust == 0; j++) {
	      if (p[j*(mNum + lbmonTot)].cluster > 0) {
		    int kkstop = 0;
	  	    for (int k = mon3; k < mon4 && kkstop == 0; k++) {
		      if (p[j*(mNum + lbmonTot) + k].monomer != 0) {
		        Vector pj_r = p[j*(mNum + lbmonTot) + k].r;
		        gcpstns(pj_r, p[i].r, p[j*(mNum + lbmonTot) + k].r, len);
	            if (dist(p[i].r, pj_r) < 2.6) {
			      int mmstop = 0;
			      for (int m = 0; m < ltot && mmstop == 0; m++) {
			        if (p[m*(mNum + lbmonTot)].ind != p[j*(mNum + lbmonTot)].ind && p[m*(mNum + lbmonTot)].cluster == p[j*(mNum + lbmonTot)].cluster) {
			  	      int nnstop = 0;
				      for (int n = mon3; n < mon4 && nnstop == 0; n++) {
					    if (p[m*(mNum + lbmonTot) + n].monomer != 0) {
				          Vector pn_r = p[m*(mNum + lbmonTot) + n].r;
		                  gcpstns(pn_r, p[i].r, p[m*(mNum + lbmonTot) + n].r, len);
			              if (dist(p[i].r, pn_r) < 2.6) {
				            p[i].cluster = p[j*(mNum + lbmonTot)].cluster;
						    for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
					          if ((*it).cluster == p[i].cluster) {
	                            (*it).drugs.push_back(p[i].r); 
							  }
                            }
					        nnstop = 1;
					        mmstop = 1;
					        kkstop = 1;
				            foundclust = 1;
						  }
			            }
			          }
			        }
		          }
		        }
		      }
            }
	      }
	    }
	    if (p[i].cluster == -1)
	      dOut++;
	    else if (p[i].cluster > 0) 
	      dIn++;
		}
	  }
	  //Else drug is bonded to lipid.
	  else {
	    //If the drug is not released (neutral conditions).
	    if (drug_released == 0) {
		  dIn = dtot;
		  dOut = 0; 
		  for (int i = mNum*ltot + btot; i < mNum*ltot + btot + dtot; i++) {
	        int foundclust = 0;
	        for (int j = 0; j < ltot && foundclust == 0; j++) {
	          if (p[j*(mNum + lbmonTot)].cluster > 0) {
		        int kkstop = 0;
	  	        for (int k = mon3; k < mon4 && kkstop == 0; k++) {
		          if (p[j*(mNum + lbmonTot) + k].monomer != 0) {
		            Vector pj_r = p[j*(mNum + lbmonTot) + k].r;
		            gcpstns(pj_r, p[i].r, p[j*(mNum + lbmonTot) + k].r, len);
	                if (dist(p[i].r, pj_r) < 2.6) {
			          int mmstop = 0;
			          for (int m = 0; m < ltot && mmstop == 0; m++) {
			            if (p[m*(mNum + lbmonTot)].ind != p[j*(mNum + lbmonTot)].ind && p[m*(mNum + lbmonTot)].cluster == p[j*(mNum + lbmonTot)].cluster) {
			  	          int nnstop = 0;
				          for (int n = mon3; n < mon4 && nnstop == 0; n++) {
					        if (p[m*(mNum + lbmonTot) + n].monomer != 0) {
				              Vector pn_r = p[m*(mNum + lbmonTot) + n].r;
		                      gcpstns(pn_r, p[i].r, p[m*(mNum + lbmonTot) + n].r, len);
			                  if (dist(p[i].r, pn_r) < 2.6) {
				                p[i].cluster = p[j*(mNum + lbmonTot)].cluster;
						        for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
					              if ((*it).cluster == p[i].cluster) {
	                                (*it).drugs.push_back(p[i].r); 
							      }
                                }
					            nnstop = 1;
					            mmstop = 1;
					            kkstop = 1;
				                foundclust = 1;
						      }
			                }
			              }
			            }
		              }
		            }
		          }
                }
	          }
	        }
		  }
		}
		//Else the drug is released (acidic conditions).
		else {
		  for (int i = mNum*ltot + btot; i < mNum*ltot + btot + dtot; i++) {
	        int foundclust = 0;
	        for (int j = 0; j < ltot && foundclust == 0; j++) {
	          if (p[j*(mNum + lbmonTot)].cluster > 0) {
		        int kkstop = 0;
	  	        for (int k = mon3; k < mon4 && kkstop == 0; k++) {
		          if (p[j*(mNum + lbmonTot) + k].monomer != 0) {
		            Vector pj_r = p[j*(mNum + lbmonTot) + k].r;
		            gcpstns(pj_r, p[i].r, p[j*(mNum + lbmonTot) + k].r, len);
	                if (dist(p[i].r, pj_r) < 2.6) {
			          int mmstop = 0;
			          for (int m = 0; m < ltot && mmstop == 0; m++) {
			            if (p[m*(mNum + lbmonTot)].ind != p[j*(mNum + lbmonTot)].ind && p[m*(mNum + lbmonTot)].cluster == p[j*(mNum + lbmonTot)].cluster) {
			  	          int nnstop = 0;
				          for (int n = mon3; n < mon4 && nnstop == 0; n++) {
					        if (p[m*(mNum + lbmonTot) + n].monomer != 0) {
				              Vector pn_r = p[m*(mNum + lbmonTot) + n].r;
		                      gcpstns(pn_r, p[i].r, p[m*(mNum + lbmonTot) + n].r, len);
			                  if (dist(p[i].r, pn_r) < 2.6) {
				                p[i].cluster = p[j*(mNum + lbmonTot)].cluster;
						        for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
					              if ((*it).cluster == p[i].cluster) {
	                                (*it).drugs.push_back(p[i].r); 
							      }
                                }
					            nnstop = 1;
					            mmstop = 1;
					            kkstop = 1;
				                foundclust = 1;
						      }
			                }
			              }
			            }
		              }
		            }
		          }
                }
	          }
	        }
	        if (p[i].cluster == -1)
	          dOut++;
	        else if (p[i].cluster > 0) 
	          dIn++;
		  }
	    }
	  }
    }
	
	//What is the real amount of enclosed drugs in clusters. 
	int dIn_real = 0;
    for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
	    dIn_real  = (*it).drugs.size();	  
	}
	
    cDrugNum_avg = cDrugNum_var = 0.0;
    dInRad_div_cRad_avg = dInRad_div_cRad_var = 0.0;
    if (dIn_real > 0) { 
	  cDrugNum_avg = dIn_real;
      cDrugNum_var = 0.0;
	  
	  //Find average
      for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
	    std::list<Vector> newdruglist;
	    Vector reference = (*it).centre;
        for (std::list<Vector>::iterator its = (*it).drugs.begin(); its != (*it).drugs.end(); ++its) {
	      Vector pn_r = *its;
	      gcpstns(pn_r, reference, *its, len); 
	      newdruglist.push_back(pn_r);
	    }
	    (*it).drugs = newdruglist;
	
	    for (std::list<Vector>::iterator its = (*it).drugs.begin(); its != (*it).drugs.end(); ++its) {
	      dInRad_div_cRad_avg += dist(*its, (*it).centre)*(1.0/(*it).radius);
	    }
      }
      dInRad_div_cRad_avg *= 1.0/dIn_real;
	  //Variance
      for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
	    for (std::list<Vector>::iterator its = (*it).drugs.begin(); its != (*it).drugs.end(); ++its) {
	      dInRad_div_cRad_var += (dist(*its, (*it).centre)*(1.0/(*it).radius) - dInRad_div_cRad_avg)*(dist(*its, (*it).centre)*(1.0/(*it).radius) - dInRad_div_cRad_avg);
	    }
      }  
      dInRad_div_cRad_var *= 1.0/dIn_real;
    }
	
  }
  //Multiple clusters forming______________________________________________________________________________________________________
  else {
	cLipNum_avg = cLipNum_var = 0.0;
	for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
	  cLipNum_avg += (*it).clustLipNum;	  
	}
	cLipNum_avg *= 1.0/(cNum);	  
    for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
	  cLipNum_var += ((*it).clustLipNum - cLipNum_avg)*((*it).clustLipNum - cLipNum_avg);	  
	}
	cLipNum_var *= 1.0/(cNum);	  
	

	//Find clustRad_avg
	cRad_avg = cRad_var = 0.0;
	for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
      (*it).centre = Vector(0.0,0.0,0.0);
	  //std::cout << "Cluster: " << (*it).cluster << std::endl;
	  for (std::list<Vector>::iterator its = (*it).firstmons.begin(); its != (*it).firstmons.end(); ++its) {
	  	  //std::cout << *its<< std::endl;
          (*it).centre += *its;
	  }
	  (*it).centre *= 1.0/(*it).firstmons.size();  
	  (*it).radius = 0.0;
	  for (std::list<Vector>::iterator its = (*it).firstmons.begin(); its != (*it).firstmons.end(); ++its) {
        (*it).radius += dist((*it).centre, *its);
	  }
	  (*it).radius *= 1.0/(*it).firstmons.size();  
	  //std::cout << (*it).radius << std::endl;
	  cRad_avg += (*it).radius;	  
	}
	cRad_avg *= 1.0/cNum;	 
	for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
	  cRad_var += ((*it).radius - cRad_avg)*((*it).radius - cRad_avg);	  
	}
	cRad_var *= 1.0/(cNum);	  
	
    //Find micRad_avg
	mRad_avg = mRad_var = 0.0; 
	
	//NEWWWWWWWWWWWWWWWWWWWW
	dIn = dOut = 0;
    if (dtot > 0) {
	  //If drug is not bonded to lipid.
	  if (drug_bonded == 0) {
	    for (int i = mNum*ltot + btot; i < mNum*ltot + btot + dtot; i++) {
	    int foundclust = 0;
	    for (int j = 0; j < ltot && foundclust == 0; j++) {
	      if (p[j*(mNum + lbmonTot)].cluster > 0) {
		    int kkstop = 0;
	  	    for (int k = mon3; k < mon4 && kkstop == 0; k++) {
		      if (p[j*(mNum + lbmonTot) + k].monomer != 0) {
		        Vector pj_r = p[j*(mNum + lbmonTot) + k].r;
		        gcpstns(pj_r, p[i].r, p[j*(mNum + lbmonTot) + k].r, len);
	            if (dist(p[i].r, pj_r) < 2.6) {
			      int mmstop = 0;
			      for (int m = 0; m < ltot && mmstop == 0; m++) {
			        if (p[m*(mNum + lbmonTot)].ind != p[j*(mNum + lbmonTot)].ind && p[m*(mNum + lbmonTot)].cluster == p[j*(mNum + lbmonTot)].cluster) {
			  	      int nnstop = 0;
				      for (int n = mon3; n < mon4 && nnstop == 0; n++) {
					    if (p[m*(mNum + lbmonTot) + n].monomer != 0) {
				          Vector pn_r = p[m*(mNum + lbmonTot) + n].r;
		                  gcpstns(pn_r, p[i].r, p[m*(mNum + lbmonTot) + n].r, len);
			              if (dist(p[i].r, pn_r) < 2.6) {
				            p[i].cluster = p[j*(mNum + lbmonTot)].cluster;
						    for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
					          if ((*it).cluster == p[i].cluster) {
	                            (*it).drugs.push_back(p[i].r); 
							  }
                            }
					        nnstop = 1;
					        mmstop = 1;
					        kkstop = 1;
				            foundclust = 1;
						  }
			            }
			          }
			        }
		          }
		        }
		      }
            }
	      }
	    }
	    if (p[i].cluster == -1)
	      dOut++;
	    else if (p[i].cluster > 0) 
	      dIn++;
		}
	  }
	  //Else drug is bonded to lipid.
	  else {
	    //If the drug is not released (neutral conditions).
	    if (drug_released == 0) {
		  dIn = dtot;
		  dOut = 0; 
		  for (int i = mNum*ltot + btot; i < mNum*ltot + btot + dtot; i++) {
	        int foundclust = 0;
	        for (int j = 0; j < ltot && foundclust == 0; j++) {
	          if (p[j*(mNum + lbmonTot)].cluster > 0) {
		        int kkstop = 0;
	  	        for (int k = mon3; k < mon4 && kkstop == 0; k++) {
		          if (p[j*(mNum + lbmonTot) + k].monomer != 0) {
		            Vector pj_r = p[j*(mNum + lbmonTot) + k].r;
		            gcpstns(pj_r, p[i].r, p[j*(mNum + lbmonTot) + k].r, len);
	                if (dist(p[i].r, pj_r) < 2.6) {
			          int mmstop = 0;
			          for (int m = 0; m < ltot && mmstop == 0; m++) {
			            if (p[m*(mNum + lbmonTot)].ind != p[j*(mNum + lbmonTot)].ind && p[m*(mNum + lbmonTot)].cluster == p[j*(mNum + lbmonTot)].cluster) {
			  	          int nnstop = 0;
				          for (int n = mon3; n < mon4 && nnstop == 0; n++) {
					        if (p[m*(mNum + lbmonTot) + n].monomer != 0) {
				              Vector pn_r = p[m*(mNum + lbmonTot) + n].r;
		                      gcpstns(pn_r, p[i].r, p[m*(mNum + lbmonTot) + n].r, len);
			                  if (dist(p[i].r, pn_r) < 2.6) {
				                p[i].cluster = p[j*(mNum + lbmonTot)].cluster;
						        for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
					              if ((*it).cluster == p[i].cluster) {
	                                (*it).drugs.push_back(p[i].r); 
							      }
                                }
					            nnstop = 1;
					            mmstop = 1;
					            kkstop = 1;
				                foundclust = 1;
						      }
			                }
			              }
			            }
		              }
		            }
		          }
                }
	          }
	        }
		  }
		}
		//Else the drug is released (acidic conditions).
		else {
		  for (int i = mNum*ltot + btot; i < mNum*ltot + btot + dtot; i++) {
	        int foundclust = 0;
	        for (int j = 0; j < ltot && foundclust == 0; j++) {
	          if (p[j*(mNum + lbmonTot)].cluster > 0) {
		        int kkstop = 0;
	  	        for (int k = mon3; k < mon4 && kkstop == 0; k++) {
		          if (p[j*(mNum + lbmonTot) + k].monomer != 0) {
		            Vector pj_r = p[j*(mNum + lbmonTot) + k].r;
		            gcpstns(pj_r, p[i].r, p[j*(mNum + lbmonTot) + k].r, len);
	                if (dist(p[i].r, pj_r) < 2.6) {
			          int mmstop = 0;
			          for (int m = 0; m < ltot && mmstop == 0; m++) {
			            if (p[m*(mNum + lbmonTot)].ind != p[j*(mNum + lbmonTot)].ind && p[m*(mNum + lbmonTot)].cluster == p[j*(mNum + lbmonTot)].cluster) {
			  	          int nnstop = 0;
				          for (int n = mon3; n < mon4 && nnstop == 0; n++) {
					        if (p[m*(mNum + lbmonTot) + n].monomer != 0) {
				              Vector pn_r = p[m*(mNum + lbmonTot) + n].r;
		                      gcpstns(pn_r, p[i].r, p[m*(mNum + lbmonTot) + n].r, len);
			                  if (dist(p[i].r, pn_r) < 2.6) {
				                p[i].cluster = p[j*(mNum + lbmonTot)].cluster;
						        for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
					              if ((*it).cluster == p[i].cluster) {
	                                (*it).drugs.push_back(p[i].r); 
							      }
                                }
					            nnstop = 1;
					            mmstop = 1;
					            kkstop = 1;
				                foundclust = 1;
						      }
			                }
			              }
			            }
		              }
		            }
		          }
                }
	          }
	        }
	        if (p[i].cluster == -1)
	          dOut++;
	        else if (p[i].cluster > 0) 
	          dIn++;
		  }
	    }
	  }
    }
	
	//What is the real amount of enclosed drugs in clusters. 
	int dIn_real = 0;
    for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
	    dIn_real  += (*it).drugs.size();	  
	}
	
	cDrugNum_avg = cDrugNum_var = 0.0;
    dInRad_div_cRad_avg = dInRad_div_cRad_var = 0.0;
    if (dIn_real > 0) { 
	  //Find average 
	  cDrugNum_avg = dIn_real*(1.0/cNum); 
	  
      for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
	    //std::cout << "Cluster: " << (*it).cluster  << " with centre " << (*it).centre << " and radius: " << (*it).radius <<  std::endl;
	    std::list<Vector> newdruglist;
	    Vector reference = (*it).centre;
        for (std::list<Vector>::iterator its = (*it).drugs.begin(); its != (*it).drugs.end(); ++its) {
	      Vector pn_r = *its;
	      gcpstns(pn_r, reference, *its, len); 
	      newdruglist.push_back(pn_r);
	    }
	    (*it).drugs = newdruglist;
	
	    for (std::list<Vector>::iterator its = (*it).drugs.begin(); its != (*it).drugs.end(); ++its) {
		  //td::cout << dist(*its, (*it).centre) << std::endl;
	      dInRad_div_cRad_avg += dist(*its, (*it).centre)*(1.0/(*it).radius);
	    }
      }
      dInRad_div_cRad_avg *= 1.0/dIn_real;
	  //Variance
      for (T7 it = clustlist.begin(); it != clustlist.end(); ++it) {
	    cDrugNum_var += ((*it).drugs.size() - cDrugNum_avg)*((*it).drugs.size() - cDrugNum_avg);	
	    for (std::list<Vector>::iterator its = (*it).drugs.begin(); its != (*it).drugs.end(); ++its) {
	      dInRad_div_cRad_var += (dist(*its, (*it).centre)*(1.0/(*it).radius) - dInRad_div_cRad_avg)*(dist(*its, (*it).centre)*(1.0/(*it).radius) - dInRad_div_cRad_avg);
	    }
      }
	  cDrugNum_var *= 1.0/cNum;	  
      dInRad_div_cRad_var *= 1.0/dIn_real;
    }
	
  }

}





void Cluster_Centre(int &ltot, double &len, double &spee, Vector &c, std::list<Vector> &headlist) {			
  std::list<Vector> newheadlist;
  Vector reference = headlist.front();
  newheadlist.push_back(reference);

  int countheads = 0;
  for (std::list<Vector>::iterator its = headlist.begin(); its != headlist.end(); ++its) {
    if (countheads != 0) {
	  Vector pn_r = *its;
	  gcpstns(pn_r, reference, *its, len); 
	  newheadlist.push_back(pn_r);
	}
    countheads++;
  }

  Vector oldc = c;
  c = Vector(0.0,0.0,0.0);
  for (std::list<Vector>::iterator its = newheadlist.begin(); its != newheadlist.end(); ++its) {
    c += *its;
  }
  c *= 1.0/ltot;  
  
  Vector oldc_r = oldc;
  gcpstns(oldc_r, c, oldc, len); 
  
  spee = 0.1*dist(c, oldc_r);

}






//Allocate 3d cell array_____________________________________________________________________________________________________________________________________________
template <typename T1, typename T2, typename T3, typename T4>
void alloc3Darray(T1 &c, int x) {
  //Allocate space for 3d dynamically allocated array for the indexed cells
  c = new T2 [x];
  
  c[0] =  new T3[x*x];
  
  for (int i = 1; i < x; i++)
    c[i] = &c[0][i*x];

  for (int i = 0; i < x; i++)
    c[i][0] = new T4 [x*x];
  
  for (int i = 0; i < x; i++) {
    for (int j = 1; j < x; j++) {
	  c[i][j] = &c[i][0][j*x];
	}
  }
}

//Deallocate 3d cell array
template <typename T>
void dealloc3Darray(T &c, int x) {
  for (int i = 0; i < x; i++) {
    delete [] c[i][0];
  }
  delete [] c[0];
  delete [] c;
}



//Create particle list based on domain dimension and particle location___________________________________________________________________________________________________
template <typename T1, typename T2>
void Particle_Lists(T1 &c, T2 &p, int &ptot, int &x, double &cell_len) {
  int i, j, k;
  //#pragma omp parallel default(none) shared(c,x) private(i,j,k)
  //{
    //#pragma omp for collapse(3)
    for (i = 0; i < x; i++) {
      for (j = 0; j < x; j++) {
	    for (k = 0; k < x; k++) {
	      c[i][j][k].partlist.clear();
	    }
	  }
    }
  //}
  
  double max, min;
  //#pragma omp parallel default(none) shared(c,p,ptot,x,cell_len) private(i,j,min,max)
  //{
    //#pragma omp for
    for (i = 0; i < ptot; i++) {
      p[i].cell = -1, -1, -1;
	  
	  for (j = 0; j < x && ( p[i].cell(0) < 0 || p[i].cell(1) < 0 || p[i].cell(2) < 0); j++) {
	    min = (-0.5*x + j)*cell_len;
		max = (-0.5*x + j + 1)*cell_len;
	    if (p[i].cell(0) < 0) {
          if ( p[i].r(0) >= min && p[i].r(0) < max ) {
		    p[i].cell(0) = j;
		  }
		}
		if (p[i].cell(1) < 0) {
		  if ( p[i].r(1) >= min && p[i].r(1) < max ) {
	        p[i].cell(1) = j;		
          }	
	    }
		if (p[i].cell(2) < 0) {
		  if ( p[i].r(2) >= min && p[i].r(2) < max ) {
	        p[i].cell(2) = j;
		  }
		}
	  }
	  if (p[i].cell(0) == -1 || p[i].cell(1) == -1 || p[i].cell(2) == -1) {
        std::cerr << "Particle " << i << " is not within the domain of the box." << std::endl;
	  }
      c[ p[i].cell(0) ][ p[i].cell(1) ][ p[i].cell(2) ].partlist.push_back(p[i]);
    }	
  //}	
  
}	

//_______________________________________________________________________________________________________________________________________________________________________
template <typename T1, typename T2>
void MPC_Particle_Lists(T1 &mpc_c, T2 &p, int &stot, int &ptot, int &mpc_x, double &mpc_cell_len) {
  int i, j, k;
  
  //#pragma omp parallel default(none) shared(mpc_c,mpc_x) private(i,j,k)
  //{
    //#pragma omp for collapse(3)
    for (i = 0; i < mpc_x; i++) {
      for (j = 0; j < mpc_x; j++) {
	    for (k = 0; k < mpc_x; k++) {
	      mpc_c[i][j][k].partlist.clear();
	    }
	  }
    }
  //}
  
  double min, max;
  //#pragma omp parallel default(none) shared(mpc_c,p,stot,mpc_x,mpc_cell_len) private(i,j,min,max)
  //{
    //#pragma omp for
    for (i = ptot - stot; i < ptot; i++) {
      p[i].mpc_cell = -1, -1, -1;    
	  
	  for (j = 0; j < mpc_x && ( p[i].mpc_cell(0) < 0 || p[i].mpc_cell(1) < 0 || p[i].mpc_cell(2) < 0 ); j++) {
	    min = (-0.5*mpc_x + j)*mpc_cell_len;
		max = (-0.5*mpc_x + j + 1)*mpc_cell_len;
		if (p[i].mpc_cell(0) < 0) {
          if ( p[i].r(0) >= min && p[i].r(0) < max ) {
		    p[i].mpc_cell(0) = j;
		  }
		}
		if (p[i].mpc_cell(1) < 0) {
		  if ( p[i].r(1) >= min && p[i].r(1) < max ) {
	        p[i].mpc_cell(1) = j;
		  }
		}
		if (p[i].mpc_cell(2) < 0) {
	      if ( p[i].r(2) >= min && p[i].r(2) < max ) {
	        p[i].mpc_cell(2) = j;
		  }
		}	
	  }
   
      mpc_c[ p[i].mpc_cell(0) ][ p[i].mpc_cell(1) ][ p[i].mpc_cell(2) ].partlist.push_back(p[i]);
    }	
  //}
  
}	





//Get rotation operator and mean cell velocity
template<typename T>
void Calc_rotangle(T &mpc_c, int &mpc_x) { 
  ENG eng;
  eng.seed(static_cast<unsigned int>(std::time(0)));
  U_DIST angle_dist(0, M_PI);
  U_GEN angle_gen(eng, angle_dist);
  
  for (int i = 0; i < mpc_x; i++) {
    for (int j = 0; j < mpc_x; j++) {
      for (int k = 0; k < mpc_x; k++) {
	    mpc_c[i][j][k].rotangle = angle_gen();	            //Get angle for rotation operator for cell[i][j][k].
	  }
	}
  }
}


/*
//______________________________________________________________________________________________________________________________________________________________________
template<typename T11>
void gcpstns(Vector &b_r, T11 &a, T11 &b, double len) {	
    if ( abs(a.r(0) - b.r(0)) > 0.5*len ) {
	  if ( b.r(0) < 0 )
	    b_r(0) += len;
	  else if ( b.r(0) > 0 )
	    b_r(0) -= len;
    }	
	if ( abs(a.r(1) - b.r(1)) > 0.5*len ) {
	  if ( b.r(1) < 0 )
	    b_r(1) += len;
	  else if ( b.r(1) > 0 )
	    b_r(1) -= len;
    }	
	if ( abs(a.r(2) - b.r(2)) > 0.5*len ) {
	  if ( b.r(2) < 0 )
	    b_r(2) += len;
	  else if ( b.r(2) > 0 )
		b_r(2) -= len;
    }
}

//______________________________________________________________________________________________________________________________________________________________________
template<typename T11>
void Neighbor_List_Cutoff_1(T11 &a, T11 &b, double &nlco, double &cutoffrad, 
                          double &wha, double &whs, 
						  double &wpp, double &wpt, double &wpd,						  
						  double &wtt, double &wtd,  
						  double &wdd,
						  double &waa, double &was) {
  //If a is a head particle
  if ( a.ptype == 1 ) {
    if (b.ptype == 5 ) 
      nlco = wha + cutoffrad;
    else if ( b.ptype == 6 ) 
      nlco = whs + cutoffrad;
  }
  //If a is a pHse particle
  else if ( a.ptype == 2 ) {
    if ( b.ptype == 2 )  
	  nlco = wpp + cutoffrad;
    else if ( b.ptype == 3 )
      nlco = wpt + cutoffrad;
    else if ( b.ptype == 4 ) 
      nlco = wpd + cutoffrad;
  }
  //If a is a tail particle
  else if ( a.ptype == 3 ) {
    if ( b.ptype == 2 )  
	  nlco = wpt + cutoffrad;
    else if ( b.ptype == 3 ) 
      nlco = wtt + cutoffrad;
    else if ( b.ptype == 4 )
      nlco = wtd + cutoffrad;
  }
  //If a is a drug particle
  else if ( a.ptype == 4 ) {
    if ( b.ptype == 4 )
      nlco = wdd + cutoffrad;
  }
  //If a is an acid particle
  else if ( a.ptype == 5 ) {
    if ( b.ptype == 5 )
	  nlco = waa + cutoffrad;
    else if ( b.ptype == 6 ) 
      nlco = was + cutoffrad;
  }
 
}  

//void Neighbor_Lists_Third(T1 &c, T2 &p, int &mNum, int &ltot, int &dtot, int &x, double &len, double &cell_len) {
//
template <typename T1, typename T2, typename T3, typename T4>
void Neighbor_Lists_1(T1 &c, T2 &p, int &stot, int &ptot, int &x, double &len, double &RC, 
                      double &WHA, double &WHS, 
				      double &WPP, double &WPT, double &WPD, 						  
					  double &WTT, double &WTD, 
					  double &WDD,
					  double &WAA, double &WAS) {
  //For each particle make a list of particles in neighbouring cells 
  int i, x_step, y_step, z_step, A, B, C;
  T3 plel;
  Vector plel_r;
  
  for (i = 0; i < ptot - stot; i++) { 
    p[i].neighlist.clear();  	  
    for (x_step  = -1; x_step <= 1; x_step++) { 
	  for (y_step  = -1; y_step <= 1; y_step++) {
		for (z_step = -1; z_step <= 1; z_step++) {
		  A = p[i].cell(0) + x_step;
		  B = p[i].cell(1) + y_step;
		  C = p[i].cell(2) + z_step;
		  if (A == -1) A += x;
		  else if (A == x) A -= x;
	      if (B == -1) B += x;
		  else if (B == x) B -= x;
		  if (C == -1) C += x;
		  else if (C == x) C -= x;
		  //Check adjacent cells to cell[i][j][k]
		  for (plel = c[A][B][C].partlist.begin(); plel != c[A][B][C].partlist.end(); ++plel) {
		    if ( p[i].lipid != (*plel).lipid && p[i].ind < (*plel).ind ) {
              double nlcutoff = RC;
			 Neighbor_List_Cutoff_1<T4>(p[i], *plel, nlcutoff, RC, WHA, WHS, WPP, WPT, WPD, WTT, WTD, WDD, WAA, WAS);	
			  double distance = dist(p[i].r, (*plel).r);
			  if (distance < nlcutoff ) {
			    p[i].neighlist.push_back(*plel);	
			  }
              else if (distance > 0.5*len) {
			    plel_r = (*plel).r;
				gcpstns<T4>(plel_r, p[i], *plel, len);
				if (dist(p[i].r, plel_r) < nlcutoff )
				  p[i].neighlist.push_back(*plel);
              }								
		    }   
	      }
	    }
	  }
	}
  }
	
}
*/


//______________________________________________________________________________________________________________________________________________________________________
template<typename T11>
void Neighbor_List_Cutoff_2(T11 &a, T11 &b, double &nlco, double &cutoffrad, 
  				                 double &wha, double &whs, 
					             double &wpa, double &wps, 								 
					             double &wtt, double &wtd, 
					             double &wdd,
					             double &was) {
  //If a is a head particle
  if ( a.ptype == 1 ) {
    if ( b.ptype == 5 ) 
      nlco = wha + cutoffrad;
    else if ( b.ptype == 6 ) 
      nlco = whs + cutoffrad;
  }
  //If a is a pHse particle
  else if ( a.ptype == 2 ) {
    if( b.ptype == 5 )
      nlco = wpa + cutoffrad;
    else if( b.ptype == 6 ) 
      nlco = wps + cutoffrad;
  }
  //If a is a tail particle
  else if ( a.ptype == 3) {
    if ( b.ptype == 3 ) 
      nlco = wtt + cutoffrad;
    else if ( b.ptype == 4 )
      nlco = wtd + cutoffrad;
  }
  //If a is a drug particle
  else if( a.ptype == 4 && b.ptype == 4 ) {
    nlco = wdd + cutoffrad;
  }
  //If a is an acid particle
  else if ( a.ptype == 5 && b.ptype == 6 ) {
    nlco = was + cutoffrad;
  }
  
}  

/*
template <typename T1, typename T2, typename T3, typename T4>
void Neighbor_Lists_2(T1 &c, T2 &p, int &stot, int &ptot, int &x, double &len, double &RC,
                           double &WHA, double &WHS, 
						   double &WPA, double &WPS, 						   
						   double &WTT, double &WTD, 
						   double &WDD,
						   double &WAS) {
  //For each particle make a list of particles in neighbouring cells 
  
  int i, x_step, y_step, z_step, A, B, C;
  T3 plel;
  Vector plel_r;
  
  //#pragma omp parallel default(none) shared(c,p,mNum,ltot,dtot,x,cell_len) private (i,x_step,y_step,z_step,A,B,C,plel,plel_r)
  //{
  //#pragma omp for
  for (i = 0; i < ptot - stot; i++) { 
    p[i].neighlist.clear();  	  
    for (x_step  = -1; x_step <= 1; x_step++) { 
	  for (y_step  = -1; y_step <= 1; y_step++) {
		for (z_step = -1; z_step <= 1; z_step++) {
		  A = p[i].cell(0) + x_step;
		  B = p[i].cell(1) + y_step;
		  C = p[i].cell(2) + z_step;
		  if (A == -1) 
		    A += x;
		  else if (A == x)
		    A -= x;
	      if (B == -1) 
		    B += x;
		  else if (B == x)
		    B -= x;
		  if (C == -1) 
		    C += x;
		  else if (C == x)
		    C -= x;
		  //Check adjacent cells to cell[i][j][k]
		  for (plel = c[A][B][C].partlist.begin(); plel != c[A][B][C].partlist.end(); ++plel) {
	        if ( p[i].lipid != (*plel).lipid && p[i].ind < (*plel).ind ) {
             // double nlcutoff = RC;
			  //Neighbor_List_Cutoff_2<T4>(p[i], *plel, nlcutoff, RC, WHA, WHS, WPA, WPS, WTT, WTD, WDD, WAS);	
		      double distance = dist(p[i].r, (*plel).r);
			  if (distance < 2.6) {
			    p[i].neighlist.push_back(*plel);	
			  }
              else if (distance > 0.5*len) {
			    plel_r = (*plel).r;
				gcpstns<T4>(plel_r, p[i], *plel, len);
				if (dist(p[i].r, plel_r) < 2.6)
				  p[i].neighlist.push_back(*plel);
              }								
		    }   
	      }
	    }
	  }
	}
  }
  
}
*/


//______________________________________________________________________________________________________________________________________________________________________



//MPC: SRD_____________________________________________________________________________________________________________________________________________________________

//Grid shift
template<typename T>
void Grid_Shift_Assign_vrand(T &p, int &stot, int &ptot, int &mpc_x, double &mpc_cell_len, double &len, double &temp, double &ms) {
  ENG eng; 
  eng.seed(static_cast<unsigned int>(std::time(0)));
  N_DIST s_dist(0.0, sqrt(temp*(1.0/ms)));
  N_GEN s_gen(eng, s_dist);
  U_DIST gs_dist(-0.5*mpc_cell_len, 0.5*mpc_cell_len);
  U_GEN gs_gen(eng, gs_dist);
  
  Vector shift(gs_gen(), gs_gen(), gs_gen());
   
  //#pragma omp parallel default(none) shared(p,stot,len,shift) private(i)
  //{
    //#pragma omp for
    for (int i = ptot - stot; i < ptot; i++) {
	  p[i].vrand = s_gen(), s_gen(), s_gen();    
      p[i].gsr = p[i].r; 
	  p[i].r += shift;
	  Boundary_Conditions<T>(p, i, len);
	}
  //}
  
}




template <typename T1, typename T2>
void Calc_vmean(T1 &mpc_c, T2 &p, int &stot, int &ptot, int &mpc_x, int &MPCcellNum_with_no_nonsolv, double &solvNum_avg_MPCcell_with_no_nonsolv, double &solvNum_var_MPCcell_with_no_nonsolv, double &mpc_cell_len) {
  int i, j, k;
  
  //#pragma omp parallel default(none) shared(mpc_c,mpc_x) private(i,j,k)
  //{
    //#pragma omp for collapse(3)
    for (i = 0; i < mpc_x; i++) {
      for (j = 0; j < mpc_x; j++) {
	    for (k = 0; k < mpc_x; k++) {
	      mpc_c[i][j][k].size = 0;
		  mpc_c[i][j][k].vmean = Vector(0.0, 0.0, 0.0);
		  mpc_c[i][j][k].vrandsum = Vector(0.0, 0.0, 0.0);
		  mpc_c[i][j][k].has_nonsolv = 0;
	    }
	  }
    }
  //}

    double min, max;
  //See how many mpc cells have a nonsolvent
  for (i = 0; i < ptot - stot; i++) {
    p[i].mpc_cell = -1, -1, -1;    
	  
	for (j = 0; j < mpc_x && ( p[i].mpc_cell(0) < 0 || p[i].mpc_cell(1) < 0 || p[i].mpc_cell(2) < 0 ); j++) {
	  min = (-0.5*mpc_x + j)*mpc_cell_len;
      max = (-0.5*mpc_x + j + 1)*mpc_cell_len;
	  if (p[i].mpc_cell(0) < 0) {
        if ( p[i].r(0) >= min && p[i].r(0) < max ) 
		  p[i].mpc_cell(0) = j;
	  }
	  if (p[i].mpc_cell(1) < 0) {
		if ( p[i].r(1) >= min && p[i].r(1) < max ) 
	      p[i].mpc_cell(1) = j;
	  }
	  if (p[i].mpc_cell(2) < 0) {
	    if ( p[i].r(2) >= min && p[i].r(2) < max ) 
	        p[i].mpc_cell(2) = j;
	   }	
	}
    mpc_c[ p[i].mpc_cell(0) ][ p[i].mpc_cell(1) ][ p[i].mpc_cell(2) ].has_nonsolv = 1;
  }	
  
  //#pragma omp parallel default(none) shared(mpc_c,p,stot,mpc_x,mpc_cell_len) private(i,j,min,max)
  //{
    //#pragma omp for
    for (i = ptot - stot; i < ptot; i++) {
      p[i].mpc_cell = -1, -1, -1;    
	  
	  for (j = 0; j < mpc_x && ( p[i].mpc_cell(0) < 0 || p[i].mpc_cell(1) < 0 || p[i].mpc_cell(2) < 0 ); j++) {
	    min = (-0.5*mpc_x + j)*mpc_cell_len;
		max = (-0.5*mpc_x + j + 1)*mpc_cell_len;
		if (p[i].mpc_cell(0) < 0) {
          if ( p[i].r(0) >= min && p[i].r(0) < max ) {
		    p[i].mpc_cell(0) = j;
		  }
		}
		if (p[i].mpc_cell(1) < 0) {
		  if ( p[i].r(1) >= min && p[i].r(1) < max ) {
	        p[i].mpc_cell(1) = j;
		  }
		}
		if (p[i].mpc_cell(2) < 0) {
	      if ( p[i].r(2) >= min && p[i].r(2) < max ) {
	        p[i].mpc_cell(2) = j;
		  }
		}	
	  }
   
      mpc_c[ p[i].mpc_cell(0) ][ p[i].mpc_cell(1) ][ p[i].mpc_cell(2) ].size++;
	  mpc_c[ p[i].mpc_cell(0) ][ p[i].mpc_cell(1) ][ p[i].mpc_cell(2) ].vmean += p[i].v;
	  mpc_c[ p[i].mpc_cell(0) ][ p[i].mpc_cell(1) ][ p[i].mpc_cell(2) ].vrandsum += p[i].vrand;
    }	
  //}
  
  //ENG eng;
  //eng.seed(static_cast<unsigned int>(std::time(0)));
  //U_SPH_DIST unitv_dist(3);
  //U_SPH_GEN unitv_gen(eng, unitv_dist);
  
  MPCcellNum_with_no_nonsolv = 0;
  solvNum_avg_MPCcell_with_no_nonsolv = 0.0;
  solvNum_var_MPCcell_with_no_nonsolv = 0.0;
  //nonzerocells = 0;
  //solvNum_avg_MPC_cell_nonzero = 0.0;
  for (i = 0; i < mpc_x; i++) {
    for (j = 0; j < mpc_x; j++) {
	  for (k = 0; k < mpc_x; k++) {
	    if ( mpc_c[i][j][k].has_nonsolv == 0 ) {
		  solvNum_avg_MPCcell_with_no_nonsolv += mpc_c[i][j][k].size;
		  MPCcellNum_with_no_nonsolv++;
		}
		if ( mpc_c[i][j][k].size > 0 ) {
		  mpc_c[i][j][k].vmean *= (1.0/mpc_c[i][j][k].size);		  	
		}
		//std::vector<double> uv = unitv_gen();	        //Get unit vector for rotation operator for cell[i][j][k].
	    //mpc_c[i][j][k].rotunitvec = Vector(uv.at(0), uv.at(1), uv.at(2));
	  }
	}
  }
  //solvNum_avg_MPC_cell_nonzero *= (1.0/nonzerocells);
  solvNum_avg_MPCcell_with_no_nonsolv *= (1.0/MPCcellNum_with_no_nonsolv);
  
  for (i = 0; i < mpc_x; i++) {
    for (j = 0; j < mpc_x; j++) {
	  for (k = 0; k < mpc_x; k++) {
	    if ( mpc_c[i][j][k].has_nonsolv == 0 ) {
		  solvNum_var_MPCcell_with_no_nonsolv += (mpc_c[i][j][k].size - solvNum_avg_MPCcell_with_no_nonsolv)*(mpc_c[i][j][k].size - solvNum_avg_MPCcell_with_no_nonsolv);
		}
	  }
	}
  }
  solvNum_var_MPCcell_with_no_nonsolv *= (1.0/MPCcellNum_with_no_nonsolv);
  
}

/*
void Average_solvNum_MPC_Cell_
    int nonzerocells = 0;
    int solvNum_avg_MPC_cell_nonzero = 0;
    int solvNum_avg_MPC_cell = 0;
    for (i = 0; i < mpc_x; i++) {
      for (j = 0; j < mpc_x; j++) {
	    for (k = 0; k < mpc_x; k++) {
		  if (mpc_c[i][j][k].size() > 0) {
		    solvNum_avg_MPC_cell_nonzero += mpc_c[i][j][k].size();
			nonzerocells++;
		  }
	      solvNum_avg_MPC_cell += mpc_c[i][j][k].size();
	    }
	  }
    }
	solvNum_avg_MPC_cell_nonzero *= (1.0/nonzerocells)
	solvNum_avg_MPC_cell *= (1.0/(mpc_xind*mpc_xind*mpc_xind));






	




//Get rotation operator and mean cell velocity
template<typename T1, typename T2>
void Calc_vmean_rotunitvec(T1 &mpc_c, int &mpc_x) { 
  
  ENG eng;
  eng.seed(static_cast<unsigned int>(std::time(0)));
  U_SPH_DIST unitv_dist(3);
  U_SPH_GEN unitv_gen(eng, unitv_dist);
  
  //Assign rotation operator to cell and calc vmean of each cell
  for (int i = 0; i < mpc_x; i++) {
    for (int j = 0; j < mpc_x; j++) {
      for (int k = 0; k < mpc_x; k++) {
		std::vector<double> uv = unitv_gen();	        //Get unit vector for rotation operator for cell[i][j][k].
	    mpc_c[i][j][k].rotunitvec = Vector(uv.at(0), uv.at(1), uv.at(2));
	    mpc_c[i][j][k].vmean = Vector(0.0, 0.0, 0.0);	//Get mean velocity for cell[i][j][k].
        for (T2 it = mpc_c[i][j][k].partlist.begin(); it != mpc_c[i][j][k].partlist.end(); ++it) {
		     mpc_c[i][j][k].vmean += ((*it).v);
        }
        mpc_c[i][j][k].vmean *= 1.0/mpc_c[i][j][k].partlist.size();
	  }
	}
  }
}

*/


//Update solvent velocities
template<typename T1, typename T2>
void Update_Solvent_Velocities(T1 &mpc_c, T2 &p, int &stot, int &ptot) {
  //Update solvent particle velocities
  int i, size;
  Vector vmean, vrandsum; 
  //double rotangle;
  
  //#pragma omp parallel default(none) shared(mpc_c,s,stot) private(i)
  //{
   // #pragma omp for
    for (i = ptot - stot; i < ptot; i++) { 
	  int X = p[i].mpc_cell(0);
	  int Y = p[i].mpc_cell(1);
	  int Z = p[i].mpc_cell(2);
	  vmean = mpc_c[X][Y][Z].vmean;
	  vrandsum = mpc_c[X][Y][Z].vrandsum;
	  size = mpc_c[X][Y][Z].size;
	  //rotunitvec = mpc_c[X][Y][Z].rotunitvec;
	  //rotangle = mpc_c[X][Y][Z].rotangle;
	  
	  //p[i].v = vmean 
	         //  + dotProduct(rotunitvec, p[i].v - vmean)*rotunitvec 
	         //  + cos(rotangle)*(p[i].v - vmean) 
			 //  - cos(rotangle)*dotProduct(rotunitvec, p[i].v - vmean)*rotunitvec 
			 //  - sin(rotangle)*crossProduct(rotunitvec, p[i].v - vmean);
			   
		p[i].v = vmean + p[i].vrand - vrandsum*(1.0/size);
			   
    }
  //}
  
}



//Return solvent positions 
template<typename T>
void Return_Positions(T &p, int &stot, int &ptot) {
  int i;
  //#pragma omp parallel default(none) shared(p,stot,ptot) private(i)
  //{
    //#pragma omp for
    for (i = ptot - stot; i < ptot; i++) {
      p[i].r = p[i].gsr;
	}
  //}
  
}








/*
core_pos += p[i(mNum + lbmonTot)].r; = 0.0;
core_pos = 0.0;
for (int i = 0; i < ltot; i++) {
  core_pos += p[i(mNum + lbmonTot)].r;
  micellerad_avg += p[i(mNum + lbmonTot)].r.nrm();
}
core_pos *= (1.0/ltot);
core_pos *= (1.0/ltot);

int counter = 0;
drug_coredist_avg = 0.0;
for (int i = mNum*ltot + btot; i < mNum*ltot + btot + dtot; i++) {
  if (p[i].drug_released == 0) {
    drug_coredist_avg += dist(p[i].r, core_pos);
	counter++;
  }
}
drug_coredist_avg *= (1.0/counter);




//clustNum clustsizeAvg drugsin drugsout %drugsin %drugsout
	  








//Grid shift
template<typename T>
void Grid_Shift(T &s, int &stot, int &mpc_x, double &mpc_cell_len, double &len) {
  ENG eng;
  eng.seed(static_cast<unsigned int>(std::time(0)));
  U_DIST gs_dist(-0.5*mpc_cell_len, 0.5*mpc_cell_len);
  U_GEN gs_gen(eng, gs_dist);
  
  int i;
  Vector shift(gs_gen(), gs_gen(), gs_gen());
   
  #pragma omp parallel default(none) shared(s,stot,len,shift) private(i)
  {
    #pragma omp for
    for (i = 0; i < stot; i++) {
      s[i].gsr = s[i].r;
	  s[i].r += shift;
	  Boundary_Conditions_1Part<T>(s, i, len);
	}
  }
  
}



//Get rotation operator and mean cell velocity
template<typename T1, typename T2>
void Calc_vmean_rotunitvec(T1 &mpc_c, int &mpc_x) { 
  
  ENG eng;
  eng.seed(static_cast<unsigned int>(std::time(0)));
  U_SPH_DIST unitv_dist(3);
  U_SPH_GEN unitv_gen(eng, unitv_dist);
  
  //Assign rotation operator to cell and calc vmean of each cell
  for (int i = 0; i < mpc_x; i++) {
    for (int j = 0; j < mpc_x; j++) {
      for (int k = 0; k < mpc_x; k++) {
		std::vector<double> uv = unitv_gen();	        //Get unit vector for rotation operator for cell[i][j][k].
	    mpc_c[i][j][k].rotunitvec = Vector(uv.at(0), uv.at(1), uv.at(2));
		
	    mpc_c[i][j][k].vmean = Vector(0.0, 0.0, 0.0);	//Get mean velocity for cell[i][j][k].
        for (T2 it = mpc_c[i][j][k].partlist.begin(); it != mpc_c[i][j][k].partlist.end(); ++it) {
		     mpc_c[i][j][k].vmean += ((*it).v);
        }
        mpc_c[i][j][k].vmean *= 1.0/mpc_c[i][j][k].partlist.size();
	  }
	}
  }
}




//Update solvent velocities
template<typename T1, typename T2>
void Update_Solvent_Velocities(T1 &mpc_c, T2 &s, int &stot) {
  //Update solvent particle velocities
  int i;
  #pragma omp parallel default(none) shared(mpc_c,s,stot) private(i)
  {
    #pragma omp for
    for (i = 0; i < stot; i++) {
	  s[i].v = mpc_c[ s[i].mpc_cell(0) ][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].vmean 
	           + dotProduct(mpc_c[ s[i].mpc_cell(0) ][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].rotunitvec, s[i].v - mpc_c[ s[i].mpc_cell(0) ][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].vmean)*mpc_c[ s[i].mpc_cell(0) ][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].rotunitvec 
	           + cos(mpc_c[ s[i].mpc_cell(0) ][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].rotangle)*(s[i].v - mpc_c[ s[i].mpc_cell(0)][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].vmean) 
			   - cos(mpc_c[ s[i].mpc_cell(0) ][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].rotangle)*dotProduct(mpc_c[ s[i].mpc_cell(0) ][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].rotunitvec, s[i].v - mpc_c[ s[i].mpc_cell(0) ][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].vmean)*mpc_c[ s[i].mpc_cell(0) ][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].rotunitvec 
			   - sin(mpc_c[ s[i].mpc_cell(0) ][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].rotangle)*crossProduct(mpc_c[ s[i].mpc_cell(0) ][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].rotunitvec, s[i].v - mpc_c[ s[i].mpc_cell(0) ][ s[i].mpc_cell(1) ][ s[i].mpc_cell(2) ].vmean);
    }
  }
  
}



//Return solvent positions 
template<typename T>
void Return_Positions(T &s, int &stot) {
  int i;
  #pragma omp parallel default(none) shared(s,stot) private(i)
  {
    #pragma omp for
    for (i = 0; i < stot; i++) {
      s[i].r = s[i].gsr;
	}
  }
  
}


template<typename T1>
void Calc_clustNum_clustSize_avg(T1 &t, int &tmonNum, int &ltot, int &x, int & clNum, double &cell_len, double &clSize_avg, double &rcl) {
  //Define all lipids as not in a cluster.
  for (int i = (tmonNum - 1)*ltot; i < tmonNum*ltot; i++) {
    t[i].cluster = 1;  
  }
  
  //Build clustlist
  std::list<int> clustlist;
  for (int i = (tmonNum - 1)*ltot; i < (tmonNum*ltot - 1); i++) {
    if (t[i].cluster == 1) {
      int counter = 1;
      for (int j = i + 1; j < tmonNum*ltot; j++) {
        if (dist(t[i].r, t[j].r) < rcl) {
	      t[j].cluster = 0;
	      counter++;
	    }
		else if (dist(t[i].r, t[j].r) > 0.5*x*cell_len) {
		  Vector two_r = t[j].r;
          if (abs(t[i].r(0) - t[j].r(0)) > 0.5*x*cell_len) {
	        if (t[j].r(0) < 0)
	          two_r(0) += x*cell_len;
	        else if (t[j].r(0) > 0)
	          two_r(0) -= x*cell_len;
          }	
	      if (abs(t[i].r(1) - t[j].r(1)) > 0.5*x*cell_len) {
	        if (t[j].r(1) < 0)
	          two_r(1) += x*cell_len;
	        else if (t[j].r(1) > 0)
	          two_r(1) -= x*cell_len;
          }	
	      if (abs(t[i].r(2) - t[j].r(2)) > 0.5*x*cell_len) {
	        if (t[j].r(2) < 0)
	          two_r(2) += x*cell_len;
	        else if (t[j].r(2) > 0)
		      two_r(2) -= x*cell_len;
          }
		  if (dist(t[i].r, two_r) < rcl) {
	        t[j].cluster = 0;
	        counter++;
	      }
		}
      }	
	  
      if (counter > 1) {
        clustlist.push_back(counter);
      }
	  
    }
  }

  //Get the number of clusters.
  clNum = clustlist.size();

  //Get the average cluster size.
  clSize_avg = 0.0;  
  if (clNum > 1) {
    for (std::list<int>::iterator it = clustlist.begin(); it != clustlist.end(); ++it) {
      clSize_avg += *it;
    }
    clSize_avg /= clustlist.size();
  }
  
  clustlist.clear();

}

template<typename T1>
void Calc_clustNum_clustSize_avg(T1 &t, int &tmonNum, int &ltot, int &x, int & clNum, double &cell_len, double &clSize_avg, double &rcl) {
  //Define all lipids as not in a cluster.
  for (int i = (tmonNum - 1)*ltot; i < tmonNum*ltot; i++) {
    t[i].cluster = 1;  
  }
  
  //Build clustlist
  std::list<int> clustlist;
  for (int i = 0; i < (tmonNum*ltot - 1); i++) {
      for (int j = i + 1; j < tmonNum*ltot; j++) {
        if (dist(t[i].r, t[j].r) < rcl) {
	      t[j].cluster = 0;
	      counter++;
	    }
		else if (dist(t[i].r, t[j].r) > 0.5*x*cell_len) {
		  Vector two_r = t[j].r;
          if (abs(t[i].r(0) - t[j].r(0)) > 0.5*x*cell_len) {
	        if (t[j].r(0) < 0)
	          two_r(0) += x*cell_len;
	        else if (t[j].r(0) > 0)
	          two_r(0) -= x*cell_len;
          }	
	      if (abs(t[i].r(1) - t[j].r(1)) > 0.5*x*cell_len) {
	        if (t[j].r(1) < 0)
	          two_r(1) += x*cell_len;
	        else if (t[j].r(1) > 0)
	          two_r(1) -= x*cell_len;
          }	
	      if (abs(t[i].r(2) - t[j].r(2)) > 0.5*x*cell_len) {
	        if (t[j].r(2) < 0)
	          two_r(2) += x*cell_len;
	        else if (t[j].r(2) > 0)
		      two_r(2) -= x*cell_len;
          }
		  if (dist(t[i].r, two_r) < rcl) {
	        t[j].cluster = 0;
	        counter++;
	      }
		}
      }	
	  
      if (counter > 1) {
        clustlist.push_back(counter);
      }
	  
    }
  }
  
  
  

  //Get the number of clusters.
  clNum = clustlist.size();

  //Get the average cluster size.
  clSize_avg = 0.0;  
  if (clNum > 1) {
    for (std::list<int>::iterator it = clustlist.begin(); it != clustlist.end(); ++it) {
      clSize_avg += *it;
    }
    clSize_avg /= clustlist.size();
  }
  
  clustlist.clear();

}



//MPC: Anderson Thermostat and conserved angular momentum______________________________________________________________________________________
template<typename T1, typename T2>
void Assign_cmv(T1 &c, int &x) { 
  ENG eng;
  eng.seed(static_cast<unsigned int>(std::time(0)));
  U_SPH_DIST unitv_dist(3);
  U_SPH_GEN unitv_gen(eng, unitv_dist);
  
  //Assign rotation operator to cell and calc vmean of each cell
  for (int i = 0; i < x; i++) {
    for (int j = 0; j < x; j++) {
      for (int k = 0; k < x; k++) {
	    c[i][j][k].cmv = Vector(0.0, 0.0, 0.0);	//Get mean velocity for cell[i][j][k].
		double counter = 0;
        for (T2 it = c[i][j][k].partlist.begin(); it != c[i][j][k].partlist.end(); ++it) {
		   if ((**it).ptype() == 2) {
		     c[i][j][k].cmv += ((**it).v);
			 counter += 1.0;
		   }
        }
        //c[i][j][k].cmv /= c[i][j][k].partlist.size();
		c[i][j][k].cmv *= 1.0/counter;
	  }
	}
  }
}



template<typename T>
void Assign_vrand(T &solv, int &stot, double &m1, double &temp)
  ENG eng;
  eng.seed(static_cast<unsigned int>(std::time(0)));
  N_DIST dist(0, sqrt(temp/m1));
  N_GEN gen(eng, dist);
  
  for (int i = 0; i < stot; i++) {
    solv[i].vrand = Vector(gen(), gen(), gen());
  }
}


//Moment of inertia tensor
template<typename T>
void Assign_Iij(T &c, int &x, double &m1) {
  for (int i = 0; i < x; i++) {
    for (int j = 0; j < x; j++) {
      for (int k = 0; k < x; k++) {
         double IXX = 0.0;
		 double IYY = 0.0;
		 double IZZ = 0.0;
		 double IXY = 0.0;
		 double IYZ = 0.0
		 double IXZ = 0.0;
        for (T2 it = c[i][j][k].partlist.begin(); it != c[i][j][k].partlist.end(); ++it) {
          if ((**it).ptype() == 2) {
		    IXX += ((**it).r(1)*(**it).r(1) + (**it).r(2)*(**it).r(2))*m1;
			IYY += ((**it).r(0)*(**it).r(0) + (**it).r(2)*(**it).r(2))*m1;
			IZZ += ((**it).r(0)*(**it).r(0) + (**it).r(1)*(**it).r(1))*m1;
			IXY += -(**it).r(0)*(**it).r(1)*m1;
			IYZ += -(**it).r(1)*(**it).r(2)*m1;
			IXZ += -(**it).r(0)*(**it).r(2)*m1;
		  }
		}
		c[i][j][k].Iij = IXX, IXY, IXZ,
		                 IXY, IYY, IYZ,
						 IXZ, IYZ, IZZ;
	  }
	}
  }
}

temaplate<typename T>
void Assign_cmr(T1 &c, int &x) { 
  for (int i = 0; i < x; i++) {
    for (int j = 0; j < x; j++) {
      for (int k = 0; k < x; k++) {
	    c[i][j][k].cmr = Vector(0.0, 0.0, 0.0);
		int counter = 0.0;
	    for (T2 it = c[i][j][k].partlist.begin(); it != c[i][j][k].partlist.end(); ++it) {
          if ((**it).ptype() == 2) {
             c[i][j][k].cmr += m1*(**it).r;
			 counter += 1.0;
		  }
		}
		c[i][j][k].cmr *= 1.0/(counter*m1);
      }
	}
  }
}

temaplate<typename T>
void Assign_cancelvrandP(T1 &c, int &x) { 
  for (int i = 0; i < x; i++) {
    for (int j = 0; j < x; j++) {
      for (int k = 0; k < x; k++) {
	    c[i][j][k].cancelvrandP = Vector(0.0, 0.0, 0.0);
		int counter = 0.0;
	    for (T2 it = c[i][j][k].partlist.begin(); it != c[i][j][k].partlist.end(); ++it) {
          if ((**it).ptype() == 2) {
            c[i][j][k].cancelvrandP += (**it).vrand;
			counter += 1.0;
		  }
		}
		c[i][j][k].cancelvrandP *= 1.0/counter;
	  }
	}
  }
}


temaplate<typename T>
void Assign_cancelvrandL(T1 &c, int &x) { 
  for (int i = 0; i < x; i++) {
    for (int j = 0; j < x; j++) {
      for (int k = 0; k < x; k++) {
	    c[i][j][k].cancelvrandL = Vector(0.0, 0.0, 0.0);
	    for (T2 it = c[i][j][k].partlist.begin(); it != c[i][j][k].partlist.end(); ++it) {
          if ((**it).ptype() == 2) {
            c[i][j][k].cancelvrandL += m1*Inverse(c[i][j][k].Iij)*crossProduct( crossProduct( (**it).r - c[i][j][k].cmr, (**it).v - (**it).vrand ), (**it).r - c[i][j][k].cmr );
		  }
		}
	  }
	}
  }
}


//Update solvent velocities
template<typename T1, typename T2>
void Update_Solvent_Velocities(T1 &c, T2 &solv, int &sTot) {
  int n, i, j, k;
    for (n = 0; n < sTot; n++) {
      i = solv[n].Cell(0);
	  j = solv[n].Cell(1);
	  k = solv[n].Cell(2);
	
	  solv[n].v = c[i][j][k].cmv + c[i][j][k].vrand - c[i][j][k].cancelvrandP - c[i][j][k].cancelvrandL;
    }
}





//Clusters

//if the distance between two or more diff t3 beads are within a certain distance of each other then they are froming a cluster
//index this cluster
*/

#endif





