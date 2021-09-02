#include "vecmat3.h"
#include <fstream>
#include <cmath>

/*
template <typename T>
void Output_xsf_File_1(const char* xsf_file, T &p, const int &mNum, int &ltot, int &btot, int &dtot, int &atot, int &stot, int steTot, int &xsfste, int &x, double &cells, double temp) {
  std::ofstream xsf(xsf_file);
  xsf << "# " << "Lipids: "<< ltot << std::endl
      << "# " << "Branched monomers: " << btot << std::endl
      << "# " << "Drugs: " << dtot << std::endl
	  << "# " << "H: " << atot << std::endl
      << "# " << "Solvent: " << stot << std::endl
	  << "# " << "kT/eps: " << temp << std::endl
	  << "# " << "Timesteps: " << steTot << std::endl
	  << "# " << "xsf steps: " << xsfste <<std::endl
       << "PRIMEVEC" << std::endl
       << cells*x*Vector(1.0, 0.0, 0.0) << std::endl
       << cells*x*Vector(0.0, 1.0, 0.0) << std::endl
	   << cells*x*Vector(0.0, 0.0, 1.0) << std::endl
	   << "ANIMSTEPS " << 1 + floor(steTot/xsfste) << std::endl
	   << "ATOMS " << 1 << std::endl;
  for (int i = 0; i < mNum*ltot + btot + dtot + atot + stot; i++) {
    xsf << i + 1 << " " << std::setprecision(16) << p[i].r << std::endl;
  }
  xsf.close();
}
  
  //___________________________
  
template <typename T> 
void Output_xsf_File_2(const char* xsf_file, T &p, const int &mNum, int &ltot, int &btot, int &dtot, int &atot, int &stot, int &ste, int &xsfste) {
  std::ofstream xsf(xsf_file, std::ofstream::app);
  xsf << "ATOMS " << 1 + (ste/xsfste) << std::endl;
  for (int i = 0; i < mNum*ltot + btot + dtot + atot + stot; i++) {
    xsf << i + 1 << " " << p[i].r << std::endl;
  }
  xsf.close();
}
*/

template <typename T>
void Output_xsf_File_1(std::ofstream xsf, T &p, const int &mNum, int &ltot, int &btot, int &dtot, int &atot, int &stot, int steTot, int &xsfste, int &x, double &cells, double temp) {
  xsf << "# " << "Lipids: "<< ltot << std::endl
      << "# " << "Branched monomers: " << btot << std::endl
      << "# " << "Drugs: " << dtot << std::endl
	  << "# " << "H: " << atot << std::endl
      << "# " << "Solvent: " << stot << std::endl
	  << "# " << "kT/eps: " << temp << std::endl
	  << "# " << "Timesteps: " << steTot << std::endl
	  << "# " << "xsf steps: " << xsfste <<std::endl
       << "PRIMEVEC" << std::endl
       << cells*x*Vector(1.0, 0.0, 0.0) << std::endl
       << cells*x*Vector(0.0, 1.0, 0.0) << std::endl
	   << cells*x*Vector(0.0, 0.0, 1.0) << std::endl
	   << "ANIMSTEPS " << 1 + floor(steTot/xsfste) << std::endl
	   << "ATOMS " << 1 << std::endl;
  for (int i = 0; i < mNum*ltot + btot + dtot + atot + stot; i++) {
    xsf << i + 1 << " " << std::setprecision(16) << p[i].r << std::endl;
  }
}
  
  //___________________________
  
template <typename T> 
void Output_xsf_File_2(std::ofstream xsf, T &p, const int &mNum, int &ltot, int &btot, int &dtot, int &atot, int &stot, int &ste, int &xsfste) {
  xsf << "ATOMS " << 1 + (ste/xsfste) << std::endl;
  for (int i = 0; i < mNum*ltot + btot + dtot + atot + stot; i++) {
    xsf << i + 1 << " " << p[i].r << std::endl;
  }
}