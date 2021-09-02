#include <fstream>
#include <iomanip>

template<typename T, size_t size_x, size_t size_y>
void Output_psf_File(const char* psf_file, T &p, 
                    int &mNum, int &ltot, int &btot, int &dtot, int &stot, 
					int (&bNum)[size_x], int (&bmonNum)[size_x][size_y], int (&pbmonTot)[size_x][size_y], int &lbmonTot,
					double &m1, double &m2, double &m3, double &m4, double &m5, double temp) {
  std::ofstream psf(psf_file);
  psf << "PSF" << std::endl
       << std::endl
	   << std::setw(8) << 1 << " " << "!NTITLE" << std::endl
	   << " REMARKS Bond specification belonging to xsf file." << std::endl
	   << " REMARKS " << "Lipids: " << ltot << std::endl
	   << " REMARKS " << "Drugs: " << dtot << std::endl
	   << " REMARKS " << "Solvent: " << stot << std::endl
	   << " REMARKS " << "kT/eps: " << temp << std::endl
	   << std::endl
	   << std::setw(8) << mNum*ltot + btot + dtot + stot << " " << "!NATOM" << std::endl;
	   
  for (int i = 0; i < mNum*ltot + btot + dtot + stot; i++) {
    if ( p[i].ptype == 1 ) {
      psf  << std::right << std::setw(8) << i + 1 << std::right << std::setw(5) << "MIC " << std::left << std::setw(5) << p[i].lipid  
	      << std::right << std::setw(3) << "LIP" << std::right << std::setw(6) << "H" << std::right << std::setw(7) << "H" << std::right << std::setw(5) 
		  << 0 << std::right << std::setw(16) << m1 << std::right << std::setw(11) << 0 << std::endl;
	}
	else if ( p[i].ptype == 2 ) {
      psf  << std::right << std::setw(8) << i + 1 << std::right << std::setw(5) << "MIC " << std::left << std::setw(5) << p[i].lipid  
	      << std::right << std::setw(3) << "LIP" << std::right << std::setw(6) << "O" << std::right << std::setw(7) << "P" << std::right << std::setw(5) 
		  << 0 << std::right << std::setw(16) << m2 << std::right << std::setw(11) << 0 << std::endl;
	}
	else if ( p[i].ptype == 3 ) {
      psf  << std::right << std::setw(8) << i + 1 << std::right << std::setw(5) << "MIC " << std::left << std::setw(5) << p[i].lipid  
	      << std::right << std::setw(3) << "LIP" << std::right << std::setw(6) << "N" << std::right << std::setw(7) << "T" << std::right << std::setw(5) 
		  << 0 << std::right << std::setw(16) << m3 << std::right << std::setw(11) << 0 << std::endl;
	}
    else if ( p[i].ptype == 4 ) {
      psf  << std::right << std::setw(8) << i + 1 << std::right << std::setw(5) << "MIC " << std::left << std::setw(5) << p[i].lipid 
	      << std::right << std::setw(3) << "DRG" << std::right << std::setw(6) << "C" << std::right << std::setw(7) << "D" << std::right << std::setw(5) 
		  << 0 << std::right << std::setw(16) << m4 << std::right << std::setw(11) << 0 << std::endl;
	}
	else if ( p[i].ptype == 5 ) {
      psf  << std::right << std::setw(8) << i + 1 << std::right << std::setw(5) << "MIC " << std::left << std::setw(5) << p[i].lipid 
	      << std::right << std::setw(3) << "SOL" << std::right << std::setw(6) << "S" << std::right << std::setw(7) << "S" << std::right << std::setw(5) 
		  << 0 << std::right << std::setw(16) << m5 << std::right << std::setw(11) << 0 << std::endl;
	}
  }
	psf << std::endl
	    << std::endl;
		
    int bondNum = (mNum - 1)*ltot;
	for (int i = 0; i < mNum; i++) {
	  for (int j = 0; j < bNum[i]; j++) {
	    bondNum += (1 + (bmonNum[i][j] - 1))*ltot;
	  }
	}    
	//pbc set {26 26 26} -all
		  
    psf	<< std::setw(8) << bondNum << " " << "!NBOND" << std::endl;
	int counter = 0;
	for (int i = 0; i < ltot; i++) {
	  for (int j = 0; j < mNum - 1; j++) {
	    psf << std::setw(8) << i*(mNum + lbmonTot) + j + 1 << std::setw(8) << i*(mNum + lbmonTot) + j + 2;
        counter += 2;	  
	    if (counter%8 == 0)
		  psf << std::endl;
	  }
	  for (int j = 0; j < mNum; j++) {
	    for (int k = 0; k < bNum[j]; k++) {
	      psf << std::setw(8) << i*(mNum + lbmonTot) + j + 1 << std::setw(8) << i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + 1;
          counter += 2;	  
	      if (counter%8 == 0)
	  	    psf << std::endl;
		}
	  }
	  for (int j = 0; j < mNum; j++) {
	    for (int k = 0; k < bNum[j]; k++) {
		  for (int m = 0; m < bmonNum[j][k] - 1; j++) {
	        psf << std::setw(8) << i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 1 << std::setw(8) << i*(mNum + lbmonTot) + mNum + pbmonTot[j][k] + m + 2;
            counter += 2;	  
	        if (counter%8 == 0)
		      psf << std::endl;
	      }
		}
	  }
	}
	
	
	
	
	
	
	
	
	
	/*


	  psf << std::setw(8) << i + 2 + ltot << std::setw(8) << i + 3;
	  counter += 2;
	  if (counter%8 == 0)
	    psf << std::endl;
	  psf << std::setw(8) << i + 3 << std::setw(8) << i + 1 + 3*ltot;
	  counter += 2;
	  if (counter%8 == 0)
	    psf << std::endl;
	}
	psf << std::endl
	    << std::endl;
	
	
	
	psf << std::setw(8) << 2*ltot << " " << "!NTHETA" << std::endl;
    counter = 0;
	for (int i = 0; i < ltot; i++) {
	  psf << std::setw(8) << i + 1 << std::setw(8) << i + 1 + ltot << std::setw(8) << i + 1 + 2*ltot;    
      counter += 3;	  
	  if (counter%9 == 0)
	    psf << std::endl;
	  psf << std::setw(8) << i + 1 + ltot << std::setw(8) << i + 1 + 2*ltot << std::setw(8) << i + 1 + 3*ltot;
	  counter += 3;
	  if (counter%9 == 0)
	    psf << std::endl;
	}
	psf << std::endl
	    << std::endl;
	
    psf << std::setw(8) << ltot << " " << "!NPHI" << std::endl;
    counter = 0;
	for (int i = 0; i < ltot; i++) {
	  psf << std::setw(8) << i + 1 << std::setw(8) << i + 1 + ltot << std::setw(8) << i + 1 + 2*ltot << std::setw(8) << i + 1 + 3*ltot; 
      counter += 4;	  
	  if (counter%8 == 0)
	    psf << std::endl;
	}
	*/
	
	psf.close();
}