//forces.h

#ifndef FORCES_H
#define FORCES_H

#include <iostream>
#include <cmath>
#include "vecmat3.h"
#include <omp.h>

//______________________________________________________________________________________________________________________________________________________________________
template<typename T11>
void gcpstns(Vector &b_r, T11 &a, T11 &b, double len) {	
    if ( std::abs(a.r(0) - b.r(0)) > 0.5*len ) {
	  if ( b.r(0) < 0 )
	    b_r(0) += len;
	  else if ( b.r(0) > 0 )
	    b_r(0) -= len;
    }	
	if ( std::abs(a.r(1) - b.r(1)) > 0.5*len ) {
	  if ( b.r(1) < 0 )
	    b_r(1) += len;
	  else if ( b.r(1) > 0 )
	    b_r(1) -= len;
    }	
	if ( std::abs(a.r(2) - b.r(2)) > 0.5*len ) {
	  if ( b.r(2) < 0 )
	    b_r(2) += len;
	  else if ( b.r(2) > 0 )
		b_r(2) -= len;
    }
}




//Just your standard LJ repulsion force for r < rc________________________________________________________________________________
Vector RepForce(double &sigma, double &strength, Vector &disty, double &scald) {
  if (scald < 0.000000001)
    std::cout << "particles too close together" << std::endl;
  return (48.0*strength*(1.0/(sigma*sigma)))*(pow(sigma*(1.0/scald), 14.0) - 0.5*pow(sigma*(1.0/scald), 8.0))*disty;
}
	double RepInt_1(double &sigma, double &strength, double &scald) {
	  return 4.0*strength*(pow(sigma*(1.0/scald), 12.0) - pow(sigma*(1.0/scald), 6.0) + 0.25);
	}
	double RepInt_2(double &sigma, double &strength, double &scald) {
      return 4.0*strength*(pow(sigma*(1.0/scald), 12.0) - pow(sigma*(1.0/scald), 6.0));
	}


//Attractive force between solvent and head particles, and tail particles of different lipids_______________________________________
Vector AttForce(double &strength, double &range, double &cutoff, Vector &disty, double &scald) {
  return ( -0.5*M_PI*strength*(1.0/(range*scald)) )*sin( M_PI*(scald - cutoff)*(1.0/range) )*disty;
}
	double AttInt(double &strength, double &range, double &cutoff, double &scald) {
      return -strength*pow( cos(0.5*M_PI*(scald - cutoff)*(1.0/range)), 2.0 );
	}


//Bonding force between beads in lipid that are next to each other for r < rc < rinf_______________________________________________
Vector BondForce(double &strength, double &maxlen, Vector &disty, double &scald) {
    if (scald > 1.5)
    std::cout << "Lipid broke" << std::endl;
  return ( -strength*(1.0/(1.0 - pow(scald*(1.0/maxlen), 2.0))) )*disty;
}
	double BondInt(double &strength, double &maxlen, double scald) {
	  return -0.5*strength*maxlen*maxlen*log(1.0 - pow(scald*(1.0/maxlen), 2.0));
	}


//Bending force for next nearest neighbor beads in lipid.___________________________________________________________________________
Vector BendForce(double &strength, double &range, Vector &disty, double &scald) {
  return -strength*( 1.0 - (4.0*range*(1.0/scald)) )*disty;
}
	double BendInt(double &strength, double &range, double scald) {
      return 0.5*strength*pow(scald - 4.0*range, 2.0);
    }



/*

//Bonding force between pHse beads and acid beads_______________________________________________
Vector AttForce_H(double &strength, double &maxlen, double &cutoff, Vector &disty, double &scald) {
  return (-strength*(scald - cutoff)/(1.0 - pow((scald - cutoff)/maxlen, 2.0)))*disty*(1.0/scald);
}

double AttInt_H(double &strength, double &epsy, double &maxlen, double &cutoff, double scald) {
 return -0.5*strength*maxlen*maxlen*log(1.0 - pow((scald - cutoff)/maxlen, 2.0)) - epsy;
}

void vector_unsigned_zero(Vector &a, double &meps) {
  if ( std::abs(a(0)) < meps ) {
	a(0) = 0.0;
  }
  if ( std::abs(a(1)) < meps ) {
	a(1) = 0.0;
  }
  if ( std::abs(a(2)) < meps ) {
	a(2) = 0.0;
  }
}

void double_unsigned_zero(double &a, double &meps) {
  if ( std::abs(a) < meps ) {
	a = 0.0;
  }
}



template<typename T1, typename T2>
void TorsForce_TorsInt(int ind1, int ind2, T2 &p, double &utot, double &length, double &strength) {
  double meps = pow(10.0, -15.0);              //Machine epsilon.
  //___________________________________________________________________________________________________
  Vector j = p[ind2 - 2].r;
  Vector k = p[ind2 - 1].r;
  Vector l = p[ind2].r;
  gcpstns<T1>(j, p[ind1], p[ind2 - 2], length);
  gcpstns<T1>(k, p[ind1], p[ind2 - 1], length);
  gcpstns<T1>(l, p[ind1], p[ind2], length);
  
  
  std::cout << "ij: " << dist(j,a.r) <<  std::endl
            << "kj: " << dist(k,j)  <<  std::endl
			<< "kl: " << dist(l,k) <<  std::endl;

  std::cout << "i pos: " << a.r <<  std::endl
            << "j pos: " << j <<  std::endl
			<< "k pos: " << k <<  std::endl
			<< "l pos: " << l <<  std::endl;
  
  std::cout << "tor force between " << a.ind << " " << p[ind2 - 2].ind  << " " <<p[ind2 - 1].ind << " " << ind2 <<  std::endl;
  
  //_______________________________________________________________________________________________________
  Vector ij = p[ind1].r - j;
  Vector kj = k - j;
  Vector kl = k - l;
  Vector m = crossProduct(ij, kj);
  Vector n = crossProduct(kj, kl);
  
  //Calc phi and sign______________________________________________________________________________________
  double sign_phi = dotProduct(ij, n);
  //std::cout << "sign of phi: " << sign_phi <<  std::endl;
  double phi = 0.0;
  if ( std::abs(sign_phi) > meps ) { 
    if ( sign_phi > 0.0 ) {
      phi  = acos(dotProduct(m*(1.0/m.nrm()), n*(1.0/n.nrm())));
    }
    else if ( sign_phi < 0.0 ) {
      phi  = -acos(dotProduct(m*(1.0/m.nrm()), n*(1.0/n.nrm())));
	}
  }   
  //std::cout << "phi: " << phi <<  std::endl;
  //std::cout << "cos(3.0*phi): " << cos(3.0*phi) <<  std::endl;
  
  //________________________________________________________________________________________________________
  if (std::abs(cos(3.0*phi)) < meps) { 
    utot += strength;
	//std::cout << "Add to uTot: " << strength <<  std::endl;
  }
  else {
    double ui = strength*( 1 + cos(3.0*phi) );
	double_unsigned_zero(ui, meps);
    utot += ui;	
    //std::cout << "Add to uTot: " << ui <<  std::endl;
  }
  
  //_________________________________________________________________________________________________________
  if (std::abs(sin(3.0*phi)) > meps) { 
    Vector fi = 3.0*strength*sin(3.0*phi)*kj.nrm()*(1.0/m.nrm2())*m;
    Vector fl = -3.0*strength*sin(3.0*phi)*kj.nrm()*(1.0/n.nrm2())*n;
    Vector fj = -fi + (1.0/kj.nrm2())*dotProduct(ij,kj)*fi - (1.0/kj.nrm2())*dotProduct(kl,kj)*fl;
    Vector fk = -fl - (1.0/kj.nrm2())*dotProduct(ij,kj)*fi + (1.0/kj.nrm2())*dotProduct(kl,kj)*fl;
    
	vector_unsigned_zero(fi, meps);
	vector_unsigned_zero(fl, meps);
	vector_unsigned_zero(fj, meps);
	vector_unsigned_zero(fk, meps);
	
    //a.a += fi;
    //d.a += fl;
    p[ind1].a += fi;
	p[ind2].a += fl;
	p[ind2 - 2].a += fj;
	p[ind2 - 1].a += fk;
	
	std::cout << "fi: " << fi << std::endl
	          << "fj: " << fj << std::endl
			  << "fk: " << fk << std::endl
			  << "fl: " << fl << std::endl; 
	//std::cout << "tor force between " << a.ind << " " << p[ind2 - 2].ind  << " " <<p[ind2 - 1].ind << " " << ind2 <<  std::endl;	
	
  }
  
  //std::cout << "tor force between " << a.ind << " " << p[ind2 - 2].ind  << " " <<p[ind2 - 1].ind << " " << ind2 <<  std::endl;	
  
  //std::cout << std::endl;
}

*/


//template<typename T1>
//Vector ExVol_Force(T1 &one, T1 &two, double &utot, double &pressTot, double &length, double &sigma, double &rc, double &strength) {
Vector ExVol_Force(Vector &one, Vector &two, double &utot, double &pressTot, double &length, double &sigma, double &rc, double &strength) {
  double u = 0.0;
  Vector f(0.0, 0.0, 0.0);
  
  //Vector two_r = two.r;
  //gcpstns_1<T1>(two_r, one, two, length);
  //Vector distance = one.r - two_r;
  Vector distance = one - two;
  double scalard = distance.nrm();
  
  if (scalard < rc) {
    f += RepForce(sigma, strength, distance, scalard);
	u += RepInt_1(sigma, strength, scalard);
  }
  
  utot += u;
  pressTot += (1.0/3.0)*dotProduct(distance, f);
  return f;
}


//template<typename T1>
Vector Att_Force(Vector &one, Vector &two, double &utot, double &pressTot, double &length, double &sigma, double &rc, double &strength, double &range) {
  Vector f(0.0, 0.0, 0.0);
  double u = 0.0;
  
  //Vector two_r = two.r;
  //gcpstns_1<T1>(two_r, one, two, length);
  //Vector distance = one.r - two_r;
  Vector distance = one - two;
  double scalard = distance.nrm();
  
  if (scalard < rc) {
      f += RepForce(sigma, strength, distance, scalard);
	  u += RepInt_2(sigma, strength, scalard);
  }
  else if (scalard <= (rc + range)) {
      f += AttForce(strength, range, rc, distance, scalard);
	  u += AttInt(strength, range, rc, scalard);
  }
  
  utot += u;
  pressTot += (1.0/3.0)*dotProduct(distance, f);
  return f;
  
}


template<typename T1>
Vector Bond_Force(T1 &one, T1 &two, double &utot, double &pressTot, double &length, double &rinf, double &kbond) {
  double u = 0.0;
  Vector f(0.0, 0.0, 0.0);
  
  Vector two_r = two.r;
  gcpstns<T1>(two_r, one, two, length);
  Vector distance = one.r - two_r;
  double scalard = distance.nrm();
  
  f += BondForce(kbond, rinf, distance, scalard);
  u += BondInt(kbond, rinf, scalard);
  
  //std::cout << scalard << " " << "bond" << std::endl; 
  
  utot += u;
  pressTot += (1.0/3.0)*dotProduct(distance, f);
  return f;
}


template<typename T1>
Vector Bend_Force(T1 &one, T1 &two, double &utot, double &pressTot, double &length, double &sigma, double &kbend) {
  double u = 0.0;
  Vector f(0.0, 0.0, 0.0);
  
  Vector two_r = two.r;
  gcpstns<T1>(two_r, one, two, length);
  Vector distance = one.r - two_r;
  double scalard = distance.nrm();

  f += BendForce(kbend, sigma, distance, scalard);
  u += BendInt(kbend, sigma, scalard);
  
  utot += u;
 pressTot += (1.0/3.0)*dotProduct(distance, f);
  return f;
}


//Forces____________________________________________________________________________________________________


Vector Drug_Bond_Force(Vector &one, Vector &two, double &utot, double &pressTot, double &length, double &sigma, double &rc, double &rinf, double &kbond, double &strength) {
  double u = 0.0;
  Vector f(0.0, 0.0, 0.0);
  
  Vector distance = one - two;
  double scalard = distance.nrm();
  //std::cout << scalard << std::endl;
  
  if (scalard < rc) {
    f += RepForce(sigma, strength, distance, scalard);
	u += RepInt_1(sigma, strength, scalard);
  }
 
  f += BondForce(kbond, rinf, distance, scalard);
  u += BondInt(kbond, rinf, scalard);
  
  utot += u;
  pressTot += (1.0/3.0)*dotProduct(distance, f);
  return f;
}


//Forces____________________________________________________________________________________________________

//1_____________________________________________For unprot lipid ____________________________________________________________
template<typename T1>
Vector ExVol_FORCE(T1 &a, T1 &b, double &utot, double &pressTot, double &length, double &sigma, double &rc, 
                   double &ehh, double &ehp, double &eht, double &epp, double &ept, double &ett) {
  Vector b_r = b.r;
  gcpstns<T1>(b_r, a, b, length);	

  int type1 = a.ptype;
  int type2 = b.ptype;
				  
  if ( type1 == 1 ) {
    if ( type2 == 1 )
      return ExVol_Force(a.r, b_r, utot, pressTot, length, sigma, rc, ehh);
    else if ( type2 == 2 )
      return ExVol_Force(a.r, b_r, utot, pressTot, length, sigma, rc, ehp);
    else if ( type2 == 3 )
      return ExVol_Force(a.r, b_r, utot, pressTot, length, sigma, rc, eht);
	else
	  return Vector(0.0, 0.0, 0.0);
  }
  else if ( type1 == 2 ) {
    if ( type2 == 1 )
      return ExVol_Force(a.r, b_r, utot, pressTot, length, sigma, rc, ehp);
    else if ( type2 == 2 )
      return ExVol_Force(a.r, b_r, utot, pressTot, length, sigma, rc, epp);
    else if ( type2 == 3 )
      return ExVol_Force(a.r, b_r, utot, pressTot, length, sigma, rc, ept);
	else
	  return Vector(0.0, 0.0, 0.0);
  }
  else if ( type1 == 3 ) {
    if ( type2 == 1 )
      return ExVol_Force(a.r, b_r, utot, pressTot, length, sigma, rc, eht);
    else if ( type2 == 2 )
      return ExVol_Force(a.r, b_r, utot, pressTot, length, sigma, rc, ept);
    else if ( type2 == 3 )
      return ExVol_Force(a.r, b_r, utot, pressTot, length, sigma, rc, ett);
	else
	  return Vector(0.0, 0.0, 0.0);
  }
  else {
    return Vector(0.0, 0.0, 0.0);
  }
}


//Function directs to nonbonded interaction between two beads___________________________________________________________
Vector Nonbond_FORCE_FreeDrug_Neutral(Vector &a, Vector &b, int &type1, int &type2, double &utot, double &pressTot, double &length, double &sigma, double &rc, 
				                      double &ehh, double &ehp, double &eht, double &ehd, double &ehs, 
					                  double &epp, double &ept, double &epd, double &eps, 					 
					                  double &ett, double &etd, double &ets, 
					                  double &edd, double &eds,
									  double &whs, 
					                  double &wpp, double &wpt, double &wpd, 
					                  double &wtt, double &wtd, 
					                  double &wdd) {
					  
	//std::cout << "Nonbond force between: " << a.ptype << " " << b.ptype << std::endl;				  
  //If a is a head particle
  if ( type1 == 1 ) {
    if ( type2 == 1 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehh);
    else if ( type2 == 2 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehp);
    else if ( type2 == 3 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eht);
    else if ( type2 == 4 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehd);
    else if ( type2 == 5 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, ehs, whs);
	else
	  return Vector(0.0, 0.0, 0.0);
  }
  //If a is a pHse particle
  else if ( type1 == 2 ) {
    if ( type2 == 1 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehp);
    else if ( type2 == 2 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, epp, wpp);	
    else if ( type2 == 3 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, ept, wpt);
    else if ( type2 == 4 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, epd, wpd);
    else if ( type2 == 5 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eps);
	else
	  return Vector(0.0, 0.0, 0.0);
  }  
  //If a is a tail particle
  else if ( type1 == 3 ) {
    if ( type2 == 1 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eht);
    else if ( type2 == 2 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, epp, wpt);	 
    else if ( type2 == 3 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, ett, wtt);
    else if ( type2 == 4 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, etd, wtd);
    else if ( type2 == 5 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ets);
	else
	  return Vector(0.0, 0.0, 0.0);
  }
  //If a is a drug particle
  else if ( type1 == 4 ) {
    if ( type2 == 4 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, edd, wdd);
    else if ( type2 == 5 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eds);
	else
	  return Vector(0.0, 0.0, 0.0);
  }	
    else {
      return Vector(0.0, 0.0, 0.0);
  }
  
}




//For neutral solution and drug attached to lipid
//Function directs to nonbonded interaction between two beads___________________________________________________________
template<typename T1>
Vector Nonbond_FORCE_BondedDrug_Neutral(T1 p, int &ind1, int &ind2, Vector &a, Vector &b, int &type1, int &type2, double &utot, double &pressTot, double &length, double &sigma, double &rc, double &rinf, double &kbond,
				                        double &ehh, double &ehp, double &eht, double &ehd, double &ehs, 
					                    double &epp, double &ept, double &epd, double &eps, 					 
					                    double &ett, double &etd, double &ets, 
					                    double &edd, double &eds,
				                        double &whs, 
					                    double &wpp, double &wpt, double &wpd, 
					                    double &wtt, double &wtd, 
					                    double &wdd) {
					  
	//std::cout << "Nonbond force between: " << a.ptype << " " << b.ptype << std::endl;				  
  //If a is a head particle
  if ( type1 == 1 ) {
    if ( type2 == 1 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehh);
    else if ( type2 == 2 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehp);
    else if ( type2 == 3 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eht);
    else if ( type2 == 4 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehd);
    else if ( type2 == 5 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, ehs, whs);
	else
	  return Vector(0.0, 0.0, 0.0);
  }
  //If a is a pHse particle
  else if ( type1 == 2 ) {
    if ( type2 == 1 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehp);
    else if ( type2 == 2 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, epp, wpp);	
    else if ( type2 == 3 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, ept, wpt);
    else if ( type2 == 4 ) {
	  if ( p[ind1].drug == p[ind2].drug ) 
	    return Drug_Bond_Force(a, b, utot, pressTot, length, sigma, rc, rinf, kbond, epd);
		//return Att_Force(a, b, utot, pressTot, length, sigma, rc, epd, wpd_bonded);  
	  else 
	    return Att_Force(a, b, utot, pressTot, length, sigma, rc, epd, wpd);  
	}
    else if ( type2 == 5 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eps);
	else
	  return Vector(0.0, 0.0, 0.0);
  }  
  //If a is a tail particle
  else if ( type1 == 3 ) {
    if ( type2 == 1 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eht);
    else if ( type2 == 2 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, epp, wpt);	 
    else if ( type2 == 3 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, ett, wtt);
    else if ( type2 == 4 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, etd, wtd);
    else if ( type2 == 5 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ets);
	else
	  return Vector(0.0, 0.0, 0.0);
  }
  //If a is a drug particle
  else if ( type1 == 4 ) {
    if ( type2 == 4 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, edd, wdd);
    else if ( type2 == 5 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eds);
	else
	  return Vector(0.0, 0.0, 0.0);
  }	
  else {
      return Vector(0.0, 0.0, 0.0);
  }
  
}




//Function directs to nonbonded interaction between two beads___________________________________________________________
template<typename T1>
Vector Nonbond_FORCE_FreeDrug_Acidic(T1 p, int &ind1, int &ind2, Vector &a, Vector &b, int &type1, int &type2, double &utot, double &pressTot, double &length, double &sigma, double &rc, double &rinf, double &kbond,
					                 double &ehh, double &ehp, double &eht, double &ehd, double &ehs, 
					                 double &epp, double &ept, double &epd, double &eps,  							    
				                     double &ett, double &etd, double &ets, 
					                 double &edd, double &eds,
					                 double &whs,
					                 double &wps, 								
					                 double &wtt, double &wtd, 
					                 double &wdd) {   	   
  //a is head particle_____________________________________________________________________
  if ( type1 == 1 ) {
    if ( type2 == 1 ) 
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehh);
    else if ( type2 == 2 ) 
	  return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehp);
    else if ( type2 == 3 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eht);
    else if ( type2 == 4 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehd);
    else if ( type2 == 5 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, ehs, whs);
	else  
	  return Vector(0.0, 0.0, 0.0);
  }
  
  //a is pHse particle_____________________________________________________________
  else if ( type1 == 2 ) {
    if ( type2 == 1 ) 
	  return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehp);
    else if ( type2 == 2 ) 
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, epp);
    else if ( type2 == 3 ) 
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ept);
    else if ( type2 == 4 ) 
        return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, epd);
    else if ( type2 == 5 ) 
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, eps, wps);
	else  
	  return Vector(0.0, 0.0, 0.0);
  }
  
  //a is tail particle________________________________________________________________
  else if ( type1 == 3 ) {
    if ( type2 == 1 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eht);
    else if ( type2 == 2 ) 
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ept);
    else if ( type2 == 3 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, ett, wtt);
    else if ( type2 == 4 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, etd, wtd);
    else if ( type2 == 5)
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ets);
	else  
	  return Vector(0.0, 0.0, 0.0);
  }
  
  //a is drug particle___________________________________________________
  else if ( type1 == 4 ) {
    if ( type2 == 4 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, edd, wdd);
    else if ( type2 == 5 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eds);
    else  
	  return Vector(0.0, 0.0, 0.0);
  }
  
  else  {
	return Vector(0.0, 0.0, 0.0);
  }  
}



//For acidic solution with drug attached to lipid or released from lipid
//Function directs to nonbonded interaction between two beads___________________________________________________________
template<typename T1>
Vector Nonbond_FORCE_BondedDrug_Acidic(T1 p, int &ind1, int &ind2, Vector &a, Vector &b, int &type1, int &type2, double &utot, double &pressTot, double &length, double &sigma, double &rc, double &rinf, double &kbond,
					                   double &ehh, double &ehp, double &eht, double &ehd, double &ehs, 
					                   double &epp, double &ept, double &epd, double &eps,  							    
				                       double &ett, double &etd, double &ets, 
					                   double &edd, double &eds,
					                   double &whs,
									   double &wps, 								
					                   double &wtt, double &wtd, 
					                   double &wdd) {   	   
  //a is head particle_____________________________________________________________________
  if ( type1 == 1 ) {
    if ( type2 == 1 ) 
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehh);
    else if ( type2 == 2 ) 
	  return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehp);
    else if ( type2 == 3 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eht);
    else if ( type2 == 4 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehd);
    else if ( type2 == 5 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, ehs, whs);
	else
	  return Vector(0.0, 0.0, 0.0);
  }
  
  //a is pHse particle_____________________________________________________________
  else if ( type1 == 2 ) {
    if ( type2 == 1 ) 
	  return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ehp);
    else if ( type2 == 2 ) 
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, epp);
    else if ( type2 == 3 ) 
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ept);
    else if ( type2 == 4) 
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, epd);
    else if ( type2 == 5 ) 
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, eps, wps);
	else
	  return Vector(0.0, 0.0, 0.0);
  }
  
  //a is tail particle________________________________________________________________
  else if ( type1 == 3 ) {
    if ( type2 == 1 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eht);
    else if ( type2 == 2 ) 
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ept);
    else if ( type2 == 3 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, ett, wtt);
    else if ( type2 == 4 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, etd, wtd);
    else if ( type2 == 5)
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, ets);
    else
	  return Vector(0.0, 0.0, 0.0);
  }
  
  //a is drug particle___________________________________________________
  else if ( type1 == 4 ) {
    if ( type2 == 4 )
      return Att_Force(a, b, utot, pressTot, length, sigma, rc, edd, wdd);
    else if ( type2 == 5 )
      return ExVol_Force(a, b, utot, pressTot, length, sigma, rc, eds);
    else
	  return Vector(0.0, 0.0, 0.0);
  }
  
  else {
      return Vector(0.0, 0.0, 0.0);
  } 

}



//void Neighbor_Lists_Third(T1 &c, T2 &p, int &mNum, int &ltot, int &dtot, int &x, double &len, double &cell_len) {
//
template <typename T1, typename T2, typename T3, typename T4, size_t size_x, size_t size_y>
void Acceleration_uTot_FreeDrug_Neutral(T1 &c, T2 &p,
                                        int (&bNum)[size_x], int (&bmonNum)[size_x][size_y], int (&pbmonTot)[size_x][size_y], int &lbmonTot,
                                        int &mNum, int &ltot, int &btot, int &stot, int &ptot, int &x, double &utot, double &pressTot, double &length,
                                        double &mh, double &mp, double &mt, double &md, double &ms,
                                        double &sigma, double &rc, double &rinf, double &kbond, double &kbend, double &kbT,
				                        double &ehh, double &ehp, double &eht, double &ehd, double &ehs, 
					                    double &epp, double &ept, double &epd, double &eps, 								
					                    double &ett, double &etd, double &ets, 
					                    double &edd, double &eds,
				                        double &whs, 
					                    double &wpp, double &wpt, double &wpd, 			            
					                    double &wtt, double &wtd, 
					                    double &wdd) {					  
	//double meps = pow(10.0, -15.0); 				 			  
  for (int i = 0; i < ptot; i++) {
    p[i].a = Vector(0.0, 0.0, 0.0);
  }
  utot = 0.0;
  pressTot = ptot*kbT;
  
  Vector f;
  //Bonded___________________________________________________________________________________________________________________________________________________________________
  //LJ rep between lipid beads
  for (int i = 0; i < ltot; i++) {
    for (int j = 0; j < mNum + lbmonTot - 1; j++) {
      for (int k = j + 1; k < mNum + lbmonTot; k++) {
	    f = ExVol_FORCE<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + k], utot, pressTot, length, sigma, rc, ehh, ehp, eht, epp, ept, ett);
	    p[i*(mNum + lbmonTot) + j].a += f;
        p[i*(mNum + lbmonTot) + k].a -= f; 
	  }
	}	
  }
  //Bond, bend, tor between main chain beads
  for (int i = 0; i < ltot; i++) {
    for (int j = 0; j < mNum - 1; j++) {  
	  f = Bond_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + j + 1], utot, pressTot, length, rinf, kbond);	
	  p[i*(mNum + lbmonTot) + j].a += f;
      p[i*(mNum + lbmonTot) + j + 1].a -= f; 
	  if ( j + 2 < mNum ) {
	    f = Bend_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + j + 2], utot, pressTot, length, sigma, kbend);
	    p[i*(mNum + lbmonTot) + j].a += f;
        p[i*(mNum + lbmonTot) + j + 2].a -= f; 
	  }
	}
  }
  if (btot > 0 ) {
    //Bond bend between main chain monomers nd side chain monomers 
    for (int i = 0; i < ltot; i++) {
      for (int j = 0; j < mNum; j++) {
	    for (int k = 0; k < bNum[j]; k++) {
	      f = Bond_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k]], utot, pressTot, length, rinf, kbond);	
	      p[i*(mNum + lbmonTot) + j].a += f;
          p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k]].a -= f; 
	      if ( 1 < bmonNum[j][k] ) {
	        f = Bend_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + 1], utot, pressTot, length, sigma, kbend);
	        p[i*(mNum + lbmonTot) + j].a += f;
            p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + 1].a -= f;   
	      }
	    }
      }
    }
    //Side chains bond bend tor shit
    for (int i = 0; i < ltot; i++) {
      for (int j = 0; j < mNum; j++) {
	    for (int k = 0; k < bNum[j]; k++) {
	      for (int m = 0; m < bmonNum[j][k] - 1; m++) {
	        f = Bond_Force<T4>(p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 1], utot, pressTot, length, rinf, kbond);	
	        p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].a += f;
            p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 1].a -= f; 
		    //std::cout << "Side chain bond force between " << p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].ind << " " << p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 1].ind  << " " << f << std::endl;
	        if ( m + 2 < bmonNum[j][k] ) {
	          f = Bend_Force<T4>(p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 2], utot, pressTot, length, sigma, kbend);
	          p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].a += f;
              p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 2].a -= f;
	          // std::cout << "Side chain bend force between " << p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].ind << " " << p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 2].ind  << " " << f << std::endl;			
		      //if ( m + 3 < bmonNum[j][k]) {
		        //Tors_Force<T4, T2>(p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 3], p, utot, pressTot, length, ktor);	
		        //std::cout << "Side chain tor force between " << p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].ind << " " << p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 3].ind  << std::endl;	
	          //}     
	        }
	      }
        }
      }
    }	
  }
  
  //Nonbonded_________________________________________________________________________________________________________
  for (int i = 0; i < ptot - stot; i++) {  	  
    for (int x_step  = -1; x_step <= 1; x_step++) { 
	  for (int y_step  = -1; y_step <= 1; y_step++) {
		for (int z_step = -1; z_step <= 1; z_step++) {
		  int A = p[i].cell(0) + x_step;
		  int B = p[i].cell(1) + y_step;
		  int C = p[i].cell(2) + z_step;
		  if (A == -1) A += x;
		  else if (A == x) A -= x;
	      if (B == -1) B += x;
		  else if (B == x) B -= x;
		  if (C == -1) C += x;
		  else if (C == x) C -= x;
		  //Check adjacent cells to cell[i][j][k]
		  for (T3 plel = c[A][B][C].partlist.begin(); plel != c[A][B][C].partlist.end(); ++plel) {
		    if ( p[i].lipid != (*plel).lipid && p[i].ind < (*plel).ind ) {
              //double nlcutoff = rc;
			  //Neighbor_List_Cutoff_1(p[i].ptype, (*plel).ptype, nlcutoff, rc, wha, whs, wpp, wpt, wpd, wtt, wtd, wdd, waa, was);	
			  double distance = dist(p[i].r, (*plel).r);
			  
			  if (distance < 2.6 ) {
			    f = Nonbond_FORCE_FreeDrug_Neutral(p[i].r, (*plel).r, p[i].ptype, (*plel).ptype, utot, pressTot, length, sigma, rc, 
	                                               ehh, ehp, eht, ehd, ehs, 
							                       epp, ept, epd, eps, 						
						                           ett, etd, ets, 
							                       edd, eds,
							                       whs, 
							                       wpp, wpt, wpd, 							
							                       wtt, wtd, 
							                       wdd);							
	              p[i].a += f;
                  p[(*plel).ind].a -= f; 
			  }
              else if (distance > 0.5*length) {
			    Vector plel_r = (*plel).r;
				gcpstns<T4>(plel_r, p[i], *plel, length);
				if (dist(p[i].r, plel_r) < 2.6 ) {
			      f = Nonbond_FORCE_FreeDrug_Neutral(p[i].r, plel_r, p[i].ptype, (*plel).ptype, utot, pressTot, length, sigma, rc, 
	                                                 ehh, ehp, eht, ehd, ehs, 
							                         epp, ept, epd, eps, 						
						                             ett, etd, ets, 
							                         edd, eds,
													 whs, 
							                         wpp, wpt, wpd, 							
							                         wtt, wtd, 
							                         wdd);							
	              p[i].a += f;
                  p[(*plel).ind].a -= f; 
				}
              }		
			  
		    }   
	      }
	    }
	  }
	}
  }
  
  for (int i = 0; i < ptot - stot; i++) {
    if ( p[i].ptype == 1 )
      p[i].a *= 1.0/mh;
    else if ( p[i].ptype == 2 )
      p[i].a *= 1.0/mp;
	else if ( p[i].ptype == 3 )
      p[i].a *= 1.0/mt;
	else if ( p[i].ptype == 4 )
      p[i].a *= 1.0/md;
  }
  for (int i = ptot - stot; i < ptot; i++) {
    p[i].a *= 1.0/ms;
  }
  
  pressTot *= (1.0/(length*length*length));
	
}



template <typename T1, typename T2, typename T3, typename T4, size_t size_x, size_t size_y>
void Acceleration_uTot_BondedDrug_Neutral(T1 &c, T2 &p,
                         int (&bNum)[size_x], int (&bmonNum)[size_x][size_y], int (&pbmonTot)[size_x][size_y], int &lbmonTot,
                         int &mNum, int &ltot, int &btot, int &stot, int &ptot, int &x, double &utot, double &pressTot, double &length,
                         double &mh, double &mp, double &mt, double &md, double &ms,
                         double &sigma, double &rc, double &rinf, double &kbond, double &kbend, double &kbT,
				         double &ehh, double &ehp, double &eht, double &ehd, double &ehs, 
					     double &epp, double &ept, double &epd, double &eps, 								
					     double &ett, double &etd, double &ets, 
					     double &edd, double &eds,
				         double &whs, 
					     double &wpp, double &wpt, double &wpd, 		            
					     double &wtt, double &wtd, 
					     double &wdd) {					  
	//double meps = pow(10.0, -15.0); 				 			  
  for (int i = 0; i < ptot; i++) {
    p[i].a = Vector(0.0, 0.0, 0.0);
  }
  utot = 0.0;
  pressTot = ptot*kbT;
  
  Vector f;
  //Bonded___________________________________________________________________________________________________________________________________________________________________
  //LJ rep between lipid beads
  for (int i = 0; i < ltot; i++) {
    for (int j = 0; j < mNum + lbmonTot - 1; j++) {
      for (int k = j + 1; k < mNum + lbmonTot; k++) {
	    f = ExVol_FORCE<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + k], utot, pressTot, length, sigma, rc, ehh, ehp, eht, epp, ept, ett);
	    p[i*(mNum + lbmonTot) + j].a += f;
        p[i*(mNum + lbmonTot) + k].a -= f; 
	  }
	}	
  }
  //Bond, bend, tor between main chain beads
  for (int i = 0; i < ltot; i++) {
    for (int j = 0; j < mNum - 1; j++) {  
	  f = Bond_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + j + 1], utot, pressTot, length, rinf, kbond);	
	  p[i*(mNum + lbmonTot) + j].a += f;
      p[i*(mNum + lbmonTot) + j + 1].a -= f; 
	  if ( j + 2 < mNum ) {
	    f = Bend_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + j + 2], utot, pressTot, length, sigma, kbend);
	    p[i*(mNum + lbmonTot) + j].a += f;
        p[i*(mNum + lbmonTot) + j + 2].a -= f; 
	  }
	}
  }
  if (btot > 0 ) {
    //Bond bend between main chain monomers nd side chain monomers 
    for (int i = 0; i < ltot; i++) {
      for (int j = 0; j < mNum; j++) {
	    for (int k = 0; k < bNum[j]; k++) {
	      f = Bond_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k]], utot, pressTot, length, rinf, kbond);	
	      p[i*(mNum + lbmonTot) + j].a += f;
          p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k]].a -= f; 
	      if ( 1 < bmonNum[j][k] ) {
	        f = Bend_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + 1], utot, pressTot, length, sigma, kbend);
	        p[i*(mNum + lbmonTot) + j].a += f;
            p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + 1].a -= f;    
	      }
	    }
      }
    }
    //Side chains bond bend tor shit
    for (int i = 0; i < ltot; i++) {
      for (int j = 0; j < mNum; j++) {
	    for (int k = 0; k < bNum[j]; k++) {
	      for (int m = 0; m < bmonNum[j][k] - 1; m++) {
	        f = Bond_Force<T4>(p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 1], utot, pressTot, length, rinf, kbond);	
	        p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].a += f;
            p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 1].a -= f; 
	        if ( m + 2 < bmonNum[j][k] ) {
	          f = Bend_Force<T4>(p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 2], utot, pressTot, length, sigma, kbend);
	          p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].a += f;
              p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 2].a -= f;   
	        }
	      }
        }
      }
    }	
  }
  
  //Nonbonded_________________________________________________________________________________________________________
  for (int i = 0; i < ptot - stot; i++) {  	  
    for (int x_step  = -1; x_step <= 1; x_step++) { 
	  for (int y_step  = -1; y_step <= 1; y_step++) {
		for (int z_step = -1; z_step <= 1; z_step++) {
		  int A = p[i].cell(0) + x_step;
		  int B = p[i].cell(1) + y_step;
		  int C = p[i].cell(2) + z_step;
		  if (A == -1) A += x;
		  else if (A == x) A -= x;
	      if (B == -1) B += x;
		  else if (B == x) B -= x;
		  if (C == -1) C += x;
		  else if (C == x) C -= x;
		  //Check adjacent cells to cell[i][j][k]
		  for (T3 plel = c[A][B][C].partlist.begin(); plel != c[A][B][C].partlist.end(); ++plel) {
		    if ( p[i].lipid != (*plel).lipid && p[i].ind < (*plel).ind ) {
			  double distance = dist(p[i].r, (*plel).r);
			  if (distance < 2.6 ) {
			    f = Nonbond_FORCE_BondedDrug_Neutral<T2>(p, p[i].ind, (*plel).ind, p[i].r, (*plel).r, p[i].ptype, (*plel).ptype, utot, pressTot, length, sigma, rc, rinf, kbond,
	                                ehh, ehp, eht, ehd, ehs, 
							        epp, ept, epd, eps, 						
						            ett, etd, ets, 
							        edd, eds,
									whs, 
							        wpp, wpt, wpd, 						
							        wtt, wtd, 
							        wdd);							
	              p[i].a += f;
                  p[(*plel).ind].a -= f; 
			  }
              else if (distance > 0.5*length) {
			    Vector plel_r = (*plel).r;
				gcpstns<T4>(plel_r, p[i], *plel, length);
				if (dist(p[i].r, plel_r) < 2.6 ) {
			      f = Nonbond_FORCE_BondedDrug_Neutral<T2>(p, p[i].ind, (*plel).ind, p[i].r, plel_r, p[i].ptype, (*plel).ptype, utot, pressTot, length, sigma, rc, rinf, kbond,
	                                  ehh, ehp, eht, ehd, ehs, 
							          epp, ept, epd, eps, 						
						              ett, etd, ets, 
							          edd, eds,
							          whs, 
							          wpp, wpt, wpd, 					
							          wtt, wtd, 
							          wdd);							
	              p[i].a += f;
                  p[(*plel).ind].a -= f; 
				}
              }		
			  
		    }   
	      }
	    }
	  }
	}
  }
  
  for (int i = 0; i < ptot - stot; i++) {
    if ( p[i].ptype == 1 )
      p[i].a *= 1.0/mh;
    else if ( p[i].ptype == 2 )
      p[i].a *= 1.0/mp;
	else if ( p[i].ptype == 3 )
      p[i].a *= 1.0/mt;
	else if ( p[i].ptype == 4 )
      p[i].a *= 1.0/md;
  }
  for (int i = ptot - stot; i < ptot; i++) {
    p[i].a *= 1.0/ms;
  }
  
  pressTot *= (1.0/(length*length*length));
	
}






































//void Neighbor_Lists_Third(T1 &c, T2 &p, int &mNum, int &ltot, int &dtot, int &x, double &len, double &cell_len) {
//
template <typename T1, typename T2, typename T3, typename T4, size_t size_x, size_t size_y>
void Acceleration_uTot_FreeDrug_Acidic(T1 &c, T2 &p,
                         int (&bNum)[size_x], int (&bmonNum)[size_x][size_y], int (&pbmonTot)[size_x][size_y], int &lbmonTot,
                         int &mNum, int &ltot, int &btot, int &dtot, int &stot, int &ptot, int &x, double &utot, double &pressTot, double &length,
                         double &mh, double &mp, double &mt, double &md, double &ms,
                         double &sigma, double &rc, double &rinf, double &kbond, double &kbend, double &kbT,
					     double &ehh, double &ehp, double &eht, double &ehd, double &ehs, 
					     double &epp, double &ept, double &epd, double &eps,  							    
				         double &ett, double &etd, double &ets, 
					     double &edd, double &eds,
					     double &whs,
					     double &wps, 								
					     double &wtt, double &wtd, 
					     double &wdd) {  			
  //double meps = pow(10.0, -15.0); 				 			  
  for (int i = 0; i < ptot; i++) {
    p[i].a = Vector(0.0, 0.0, 0.0);
  }
  utot = 0;
  pressTot = ptot*kbT;
  
  Vector f;
 
  //Bonded___________________________________________________________________________________________________________________________________________________________________
  //LJ rep between lipid beads
  for (int i = 0; i < ltot; i++) {
    for (int j = 0; j < mNum + lbmonTot - 1; j++) {
      for (int k = j + 1; k < mNum + lbmonTot; k++) {
	    f = ExVol_FORCE<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + k], utot, pressTot, length, sigma, rc, ehh, ehp, eht, epp, ept, ett);
	    p[i*(mNum + lbmonTot) + j].a += f;
        p[i*(mNum + lbmonTot) + k].a -= f; 
	    //std::cout << "Rep force between " << p[i*(mNum + lbmonTot) + j].ind << " " << p[i*(mNum + lbmonTot) + k].ind  << " " << f << std::endl;
	  }
	}	
  }
  //Bond, bend, tor between main chain beads
  for (int i = 0; i < ltot; i++) {
    for (int j = 0; j < mNum - 1; j++) {  
	  f = Bond_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + j + 1], utot, pressTot, length, rinf, kbond);	
	  p[i*(mNum + lbmonTot) + j].a += f;
      p[i*(mNum + lbmonTot) + j + 1].a -= f; 
	  //std::cout << "Main chain bond force between " << p[i*(mNum + lbmonTot) + j].ind  << " " << p[i*(mNum + lbmonTot) + j + 1].ind  << " " << f << std::endl;
	  if ( j + 2 < mNum ) {
	    f = Bend_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + j + 2], utot, pressTot, length, sigma, kbend);
	    p[i*(mNum + lbmonTot) + j].a += f;
        p[i*(mNum + lbmonTot) + j + 2].a -= f; 
		//std::cout << "uTot before tor: " << utot << std::endl;
		//std::cout << "Main chain bend force between " << p[i*(mNum + lbmonTot) + j].ind  << " " << p[i*(mNum + lbmonTot) + j + 2].ind  << " " << f << std::endl;
	  }
	}
  }
  if (btot > 0 ) {
    //Bond bend between main chain monomers nd side chain monomers 
    for (int i = 0; i < ltot; i++) {
      for (int j = 0; j < mNum; j++) {
	    for (int k = 0; k < bNum[j]; k++) {
	      f = Bond_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k]], utot, pressTot, length, rinf, kbond);	
	      p[i*(mNum + lbmonTot) + j].a += f;
          p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k]].a -= f; 
	      if ( 1 < bmonNum[j][k] ) {
	        f = Bend_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + 1], utot, pressTot, length, sigma, kbend);
	        p[i*(mNum + lbmonTot) + j].a += f;
            p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + 1].a -= f; 
	      }
	    }
      }
    }
    //Side chains bond bend tor shit
    for (int i = 0; i < ltot; i++) {
      for (int j = 0; j < mNum; j++) {
	    for (int k = 0; k < bNum[j]; k++) {
	      for (int m = 0; m < bmonNum[j][k] - 1; m++) {
	        f = Bond_Force<T4>(p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 1], utot, pressTot, length, rinf, kbond);	
	        p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].a += f;
            p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 1].a -= f; 
	        if ( m + 2 < bmonNum[j][k] ) {
	          f = Bend_Force<T4>(p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 2], utot, pressTot, length, sigma, kbend);
	          p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].a += f;
              p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 2].a -= f;
	        }
	      }
        }
      }
    }	
  }
  
  //Nonbonded_________________________________________________________________________________________________________
  for (int i = 0; i < ptot - stot; i++) {  	  
    for (int x_step  = -1; x_step <= 1; x_step++) { 
	  for (int y_step  = -1; y_step <= 1; y_step++) {
		for (int z_step = -1; z_step <= 1; z_step++) {
		  int A = p[i].cell(0) + x_step;
		  int B = p[i].cell(1) + y_step;
		  int C = p[i].cell(2) + z_step;
		  if (A == -1) A += x;
		  else if (A == x) A -= x;
	      if (B == -1) B += x;
		  else if (B == x) B -= x;
		  if (C == -1) C += x;
		  else if (C == x) C -= x;
		  for (T3 plel = c[A][B][C].partlist.begin(); plel != c[A][B][C].partlist.end(); ++plel) {
		    if ( p[i].lipid != (*plel).lipid && p[i].ind < (*plel).ind) {
			  double distance = dist(p[i].r, (*plel).r);
			  if (distance < 2.6 ) {
			    f = Nonbond_FORCE_FreeDrug_Acidic<T2>(p, p[i].ind, (*plel).ind, p[i].r, (*plel).r, p[i].ptype, (*plel).ptype, utot, pressTot, length, sigma, rc, rinf, kbond,
					                                  ehh, ehp, eht, ehd, ehs, 
					                                  epp, ept, epd, eps,  							    
				                                      ett, etd, ets, 
					                                  edd, eds,
					                                  whs,
					                                  wps, 								
					                                  wtt, wtd, 
					                                  wdd);				
	              p[i].a += f;
                  p[(*plel).ind].a -= f; 
			  }
              else if (distance > 0.5*length) {
			    Vector plel_r = (*plel).r;
				gcpstns<T4>(plel_r, p[i], *plel, length);
				if (dist(p[i].r, plel_r) < 2.6 ) {
			      f = Nonbond_FORCE_FreeDrug_Acidic<T2>(p, p[i].ind, (*plel).ind, p[i].r, plel_r, p[i].ptype, (*plel).ptype, utot, pressTot, length, sigma, rc, rinf, kbond,
					                                    ehh, ehp, eht, ehd, ehs, 
					                                    epp, ept, epd, eps,  							    
				                                        ett, etd, ets, 
					                                    edd, eds,
					                                    whs,
														wps, 								
					                                    wtt, wtd, 
					                                    wdd);						
	              p[i].a += f;
                  p[(*plel).ind].a -= f;  
				}
              }		
		    }   
	      }
	    }
	  }
	}
  }
  
  for (int i = 0; i < ptot - stot; i++) {
    if ( p[i].ptype == 1 )
      p[i].a *= 1.0/mh;
    else if ( p[i].ptype == 2 )
      p[i].a *= 1.0/mp;
	else if ( p[i].ptype == 3 )
      p[i].a *= 1.0/mt;
	else if ( p[i].ptype == 4 )
      p[i].a *= 1.0/md;
  }
  for (int i = ptot - stot; i < ptot; i++) {
    p[i].a *= 1.0/ms;
  }
  
  pressTot *= (1.0/(length*length*length));
	
}



//void Neighbor_Lists_Third(T1 &c, T2 &p, int &mNum, int &ltot, int &dtot, int &x, double &len, double &cell_len) {
//
template <typename T1, typename T2, typename T3, typename T4, size_t size_x, size_t size_y>
void Acceleration_uTot_BondedDrug_Acidic(T1 &c, T2 &p,
                                         int (&bNum)[size_x], int (&bmonNum)[size_x][size_y], int (&pbmonTot)[size_x][size_y], int &lbmonTot,
                                         int &mNum, int &ltot, int &btot, int &dtot, int &stot, int &ptot, int &x, double &utot, double &pressTot, double &length,
                                         double &mh, double &mp, double &mt, double &md, double &ms,
                                         double &sigma, double &rc, double &rinf, double &kbond, double &kbend, double &kbT,
					                     double &ehh, double &ehp, double &eht, double &ehd, double &ehs, 
					                     double &epp, double &ept, double &epd, double &eps,  							    
				                         double &ett, double &etd, double &ets, 
					                     double &edd, double &eds,
					                     double &whs,
					                     double &wps, 								
					                     double &wtt, double &wtd, 
					                     double &wdd) {  			
  //double meps = pow(10.0, -15.0); 				 			  
  for (int i = 0; i < ptot; i++) {
    p[i].a = Vector(0.0, 0.0, 0.0);
  }
  utot = 0;
  pressTot = ptot*kbT;
  
  Vector f;
  
  //Bonded___________________________________________________________________________________________________________________________________________________________________
  //LJ rep between lipid beads
  for (int i = 0; i < ltot; i++) {
    for (int j = 0; j < mNum + lbmonTot - 1; j++) {
      for (int k = j + 1; k < mNum + lbmonTot; k++) {
	    f = ExVol_FORCE<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + k], utot, pressTot, length, sigma, rc, ehh, ehp, eht, epp, ept, ett);
	    p[i*(mNum + lbmonTot) + j].a += f;
        p[i*(mNum + lbmonTot) + k].a -= f; 
	    //std::cout << "Rep force between " << p[i*(mNum + lbmonTot) + j].ind << " " << p[i*(mNum + lbmonTot) + k].ind  << " " << f << std::endl;
	  }
	}	
  }
  //Bond, bend, tor between main chain beads
  for (int i = 0; i < ltot; i++) {
    for (int j = 0; j < mNum - 1; j++) {  
	  f = Bond_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + j + 1], utot, pressTot, length, rinf, kbond);	
	  p[i*(mNum + lbmonTot) + j].a += f;
      p[i*(mNum + lbmonTot) + j + 1].a -= f; 
	  //std::cout << "Main chain bond force between " << p[i*(mNum + lbmonTot) + j].ind  << " " << p[i*(mNum + lbmonTot) + j + 1].ind  << " " << f << std::endl;
	  if ( j + 2 < mNum ) {
	    f = Bend_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + j + 2], utot, pressTot, length, sigma, kbend);
	    p[i*(mNum + lbmonTot) + j].a += f;
        p[i*(mNum + lbmonTot) + j + 2].a -= f; 
		//std::cout << "uTot before tor: " << utot << std::endl;
		//std::cout << "Main chain bend force between " << p[i*(mNum + lbmonTot) + j].ind  << " " << p[i*(mNum + lbmonTot) + j + 2].ind  << " " << f << std::endl;
	  }
	}
  }
  if (btot > 0 ) {
    //Bond bend between main chain monomers nd side chain monomers 
    for (int i = 0; i < ltot; i++) {
      for (int j = 0; j < mNum; j++) {
	    for (int k = 0; k < bNum[j]; k++) {
	      f = Bond_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k]], utot, pressTot, length, rinf, kbond);	
	      p[i*(mNum + lbmonTot) + j].a += f;
          p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k]].a -= f; 
		  //std::cout << "Main to side chain bond force between " << p[i*(mNum + lbmonTot) + j].ind << " " << p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k]].ind  << " " << f << std::endl;
	      if ( 1 < bmonNum[j][k] ) {
	        f = Bend_Force<T4>(p[i*(mNum + lbmonTot) + j], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + 1], utot, pressTot, length, sigma, kbend);
	        p[i*(mNum + lbmonTot) + j].a += f;
            p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + 1].a -= f; 
		    //std::cout << "Main to side chain bend force between " << p[i*(mNum + lbmonTot) + j].ind << " " << p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + 1].ind  << " " << f << std::endl;
	      }
	    }
      }
    }
    //Side chains bond bend tor shit
    for (int i = 0; i < ltot; i++) {
      for (int j = 0; j < mNum; j++) {
	    for (int k = 0; k < bNum[j]; k++) {
	      for (int m = 0; m < bmonNum[j][k] - 1; m++) {
	        f = Bond_Force<T4>(p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 1], utot, pressTot, length, rinf, kbond);	
	        p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].a += f;
            p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 1].a -= f; 
		    //std::cout << "Side chain bond force between " << p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].ind << " " << p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 1].ind  << " " << f << std::endl;
	        if ( m + 2 < bmonNum[j][k] ) {
	          f = Bend_Force<T4>(p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m], p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 2], utot, pressTot, length, sigma, kbend);
	          p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].a += f;
              p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 2].a -= f;
	          // std::cout << "Side chain bend force between " << p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m].ind << " " << p[i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 2].ind  << " " << f << std::endl;			    
	        }
	      }
        }
      }
    }	
  }
  
  //Nonbonded_________________________________________________________________________________________________________
  for (int i = 0; i < ptot - stot; i++) {   
    for (int x_step  = -1; x_step <= 1; x_step++) { 
	  for (int y_step  = -1; y_step <= 1; y_step++) {
		for (int z_step = -1; z_step <= 1; z_step++) {
		  int A = p[i].cell(0) + x_step;
		  int B = p[i].cell(1) + y_step;
		  int C = p[i].cell(2) + z_step;
		  if (A == -1) A += x;
		  else if (A == x) A -= x;
	      if (B == -1) B += x;
		  else if (B == x) B -= x;
		  if (C == -1) C += x;
		  else if (C == x) C -= x;
		  for (T3 plel = c[A][B][C].partlist.begin(); plel != c[A][B][C].partlist.end(); ++plel) {
		    if ( p[i].lipid != (*plel).lipid && p[i].ind < (*plel).ind ) {
			  double distance = dist(p[i].r, (*plel).r);
			  if (distance < 2.6 ) {
			    f = Nonbond_FORCE_BondedDrug_Acidic<T2>(p, p[i].ind, (*plel).ind, p[i].r, (*plel).r, p[i].ptype, (*plel).ptype, utot, pressTot, length, sigma, rc, rinf, kbond,
					                                     ehh, ehp, eht, ehd, ehs, 
					                                     epp, ept, epd, eps,  							    
				                                         ett, etd, ets, 
					                                     edd, eds,
					                                     whs,
														 wps, 								
					                                     wtt, wtd, 
					                                     wdd);				
	              p[i].a += f;
                  p[(*plel).ind].a -= f; 
			  }
              else if (distance > 0.5*length) {
			    Vector plel_r = (*plel).r;
				gcpstns<T4>(plel_r, p[i], *plel, length);
				if (dist(p[i].r, plel_r) < 2.6 ) {
			      f = Nonbond_FORCE_BondedDrug_Acidic<T2>(p, p[i].ind, (*plel).ind, p[i].r, plel_r, p[i].ptype, (*plel).ptype, utot, pressTot, length, sigma, rc, rinf, kbond,
					                                      ehh, ehp, eht, ehd, ehs, 
					                                      epp, ept, epd, eps,  							    
				                                          ett, etd, ets, 
					                                      edd, eds,
					                                      whs,
					                                      wps, 								
					                                      wtt, wtd, 
					                                      wdd);						
	              p[i].a += f;
                  p[(*plel).ind].a -= f;  
				}
              }		
			  
		    }   
	      }
	    }
	  }
	}
  }
  
  for (int i = 0; i < ptot - stot; i++) {
    if ( p[i].ptype == 1 )
      p[i].a *= 1.0/mh;
    else if ( p[i].ptype == 2 )
      p[i].a *= 1.0/mp;
	else if ( p[i].ptype == 3 )
      p[i].a *= 1.0/mt;
	else if ( p[i].ptype == 4 )
      p[i].a *= 1.0/md;
  }
  for (int i = ptot - stot; i < ptot; i++) {
    p[i].a *= 1.0/ms;
  }
  
  pressTot *= (1.0/(length*length*length));
	
}

#endif