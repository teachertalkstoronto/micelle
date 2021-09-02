//classes

#ifndef CLASSES_H
#define CLASSES_H

#include "vecmat3.h"

//Beads/particles have positions, velocities, accelerations, and neighbor lists
class Bead {
  public:
    int ptype, monomer, lipid, ind, cluster, drug, counted;
    vecmat3::Vector<int> mpc_cell;
    Vector gsr, vrand;
    Vector r, v, a;                   //Each bead has a position, velocity, acceleration, and is in a specific cell
	vecmat3::Vector<int> cell;
	std::list<Bead> neighlist;             //Neighbor list is list of Bead pointers that point to specific type of particle (head/tail/solv)
    //virtual int ptype() {
	  //return -1;                            //Particle type? 0 for head particle, 1 for tail particle, and 2 for solvent particle
	//}
}; 

/*

//Head, Tail, and Solv classes are derived from base class Bead and each have their own particle type designated by 0,1, or 2.

class Head : public Bead {
    public: 
    int ptype() {
      return 1;
    }
};

class pHSe : public Bead {
  public: 
    int ptype() {
      return 2;
    }
	pHSe();
};
pHSe::pHSe() {
    protonated = 0; //true 1 false 0
    proton = -1;    // which proton is it bonded to
};

class Tail : public Bead {
  public: 
    int ptype() {
      return 3;
    }
	int cluster;
};


class Drug : public Bead {
  public: 
    int ptype() {
      return 4;
    }
};

class Acid : public Bead {
  public: 
    int ptype() {
      return 5;
    }
	Acid();
};
Acid::Acid() {
  proton = -1;    // which proton is it bonded to
};

class Solv : public Bead {
  public: 
    int ptype() {
      return 6;
    }
};

class Lipid {
  public:
    int cluster;
	int lipid;
	
	Bead** c
};
*/

class LatticePos {
  public:
    int taken;
	Vector r;
	LatticePos();
};

LatticePos::LatticePos() {
  taken = 0;
};









class Cell {
  public:
    std::list<Bead> partlist;
};

class MPC_Cell {
  public:
    //std::list<Bead> partlist;
	Vector vmean, vrandsum;
	//Vector rotunitvec;
	//double rotangle;
	int size;
	int size_solv;
	int has_nonsolv;
};

class Cluster {
 public:
   int cluster;
   int clustLipNum;  //lipids
   std::list<Vector> drugs;
   std::list<Vector> firstmons;
   Vector centre;
   double radius;
};





  
#endif  