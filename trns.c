/************************************************************************** 
*                        TwoRns - Two Fluid RNS
* 
**************************************************************************/
//
#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h> 

#include "consts.h"
#include "struct.h"

#include "nrutil.h"
#include "equil.h"
#include "equil_util.h"
#include "findmodel.h"
#include "surface.h"
#include "stableorbit.h"
#include "interpol.h"


/* Main; where it all starts and ends */

int main(int argc, char **argv)     /* Number of command line arguments, Command line arguments */
{ NeutronStar star;
  EOS eos;
  NeutronStar starDM;
  EOSDM eosDM;
  int i, ierr;
  double
    e_min,e_minDM,  e_max, e_max_mass, ratio_r = 1.0,
    e_center=1e15,  e_centerDM=1e15,                   /* central en. density */
    B,                            /* Quark Bag Constant */
    K=3.0,                        /* Second parameter in "quark" eos */
    spin_freq=100,                  /* Spin Frequency */
    Gamma_P=0.0,                      /* Gamma for polytropic EOS */
    m_chi = 0.0,                  /* DM particle mass*/
    y_chi = 0.0;                  /* DM self-interaction strength */
                
  int j;

  int a = 0, numseq=2;
  int spin_lim = 0;
  float e_c[4], M_0[4];
  float M0, Mtot,Mstat, Rstat, Radius, freq,freqK,energy_value, temp_energy, ej,Volp;
  float maxmass, maxradius;   // Mass and radius of the maximum mass neutron star for an EOS
  float T, W;
  long int angmom;
  //double Kfreq, Kfreq_j;

  //andreas
  double lumi;
  double u_lumi;
  double poten[4];
  int call;
  double ratio_ch=0.005;
  double energy_min=0.2;//10^{15}$ g/$cm^3
  double e_ch=0.05;//10^{15}$ g/$cm^3
  FILE *out;

  //andreas

  FILE *fpointer;

  char eos_file[80] = "no EOS file specified";   /* EOS file name */
  char eos_type[80] = "tab";                     /* EOS type (poly or tab) */
    char eos_typeDM[80] = "f";                     /* EOS type (b for bosonic, f for fermionic, EOS file name for tabulated) */
  char data_dir[80] = "junk";                    /* Data output directory */
  char filename[100] = "Many Fermion";


  //andreas
    
  /* READ IN THE COMMAND LINE OPTIONS */
  for(i=1;i<argc;i++) 
    if(argv[i][0]=='-'){
      switch(argv[i][1]){

      case 'b':
	/* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
	   NAME OF THE FILE */
	sscanf(argv[i+1],"%s",eos_file);
	break;
	
      case 'd':
              /* CHOOSE WHETHER THE DM IS BOSONIC OR
                 FERMIONIC OR TABULATED: 'f' for fermionic, 'b' for bosonic,
               IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
                  NAME OF THE FILE */
	sscanf(argv[i+1],"%s",eos_typeDM);
	break;
              
      case 'm':
    /* CHOOSE THE DM PARTICLE MASS (IN MeV) */
    sscanf(argv[i+1],"%lf",&m_chi);
    m_chi *= 1.60218e-6*G/(sqrt(KAPPA)*C*C*C*C);
    break;
              
      case 'y':
    /* CHOOSE THE DM SELF-INTERACTION STRENGTH */
    sscanf(argv[i+1],"%lf",&y_chi);
    break;
	
      case 'e':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN MeV/fm^3) */
	sscanf(argv[i+1],"%lf",&e_min);
	if(strcmp(eos_type,"poly")!=0)
	  e_min *= 1.60218e-6*KSCALE/(1.0e-13*1.0e-13*1.0e-13);
              printf("e_Bc_cgs = %.15e ", e_min*C*C/(KAPPA*G));
	break;
	
      case 'c':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   DARK MATTER (IN MeV/fm^3) */
	sscanf(argv[i+1],"%lf",&e_minDM);
	if(strcmp(eos_typeDM,"poly")!=0)
	  e_minDM *= 1.60218e-6*KSCALE/(1.0e-13*1.0e-13*1.0e-13);
	break;
	
      case 'r':
	/* r_ratio */
	sscanf(argv[i+1],"%lf",&ratio_r);
	break;
      }
    }
    
  //strncat(filename, eos_file, 65);
  strncat(filename, ".txt", 65);
  //printf("%s\n", filename);
  fpointer = fopen(filename, "a");


  //strncat(filename, eos_file, 65);
  //strncat(filename, "_table.txt", 65);
  //printf("%s\n", filename);
  out = fopen(filename, "a"); 
  /* PRINT THE HEADER */
  if(strcmp(eos_type,"tab")==0)
    printf("%s,  MDIVxSDIV=%dx%d\n",eos_file,MDIV,SDIV);
  
  if(strcmp(eos_type,"quark")==0)
    printf("Quark star with B=%f, MDIVxSDIV=%dx%d\n",B/1.602e33/KSCALE,MDIV,SDIV);
    

  /* SetUpStar loads in the eos and sets up the grid */
  /* Source code for SetUpStar can be found in findmodel.c */


  //printf("The star infrastructure has been set up! \n");

  e_center = e_min;
  e_centerDM = e_minDM;
  //temp_energy = e_center;

    
    //Computing the non-rotating spherical neutron star
      ierr = SetUpStar(eos_file, eos_type, data_dir, Gamma_P, B, K,&eos, &star, eos_typeDM, m_chi, y_chi ,&eosDM, &starDM);
    
      ierr = MakeSphere(&eos, &star, e_center,&eosDM, &starDM, e_centerDM);


      // The "1.0" means, compute a star with r_ratio = 1.0
      rns(1.0, e_center, &eos, &star, e_centerDM, &eosDM, &starDM);
   
  
      // Compute a spinning star
      if (ratio_r < 1.0)
	rns(ratio_r, e_center, &eos, &star, e_centerDM, &eosDM, &starDM);  
 
   //if(ratio_r==1.00) printf(" %f \n",star.MassDM/MSUN/(star.Mass/MSUN+star.MassDM/MSUN));
      /*  if(!((isnan((star.Mass))&&isnan((star.MassDM)))||(isnan((star.OmegaDM))&&isnan((star.Omega))))){

      printf("%s %g %g %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.3f %.3f %.5f\n",
       eos_fileDM, star.e_center, star.e_centerDM, star.Mass/MSUN+star.MassDM/MSUN, star.Mass/MSUN, star.MassDM/MSUN, star.R_e*1e-5, star.R_eDM*1e-5, star.Ratio_sch,ratio_r, star.r_ratioDM, star.Omega/(2.0*PI),star.OmegaDM/(2.0*PI), star.Omega_K/(2.0*PI));
           
      fprintf(fpointer,"%s %g %g %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.3f %.3f %.5f\n",
       eos_fileDM, star.e_center, star.e_centerDM, star.Mass/MSUN+star.MassDM/MSUN, star.Mass/MSUN, star.MassDM/MSUN, star.R_e*1e-5, star.R_eDM*1e-5, star.Ratio_sch,ratio_r, star.r_ratioDM, star.Omega/(2.0*PI),star.OmegaDM/(2.0*PI), star.Omega_K/(2.0*PI));
     
  //if(star.MassDM/(star.Mass+star.MassDM)<0.005) exit(0);

  }*/
   

 


  fclose(out);//andreas
  fclose(fpointer);
  return 0;
}
