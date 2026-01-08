/*****************************************************************************
*	equil.c
*
*		The code in this file is a set of procedures written by
*	Nikolaos Stergioulas. These are the procedures used to integrate
*	the field equations for a rapidly rotating neutron star.
*
* 	The most important procedures are:
*	
*	make_grid:	Create the MDIV x SDIV grid. 
*			MDIV = number of divisions of variable mu=cos theta
*			SDIV = number of divisions of radial variable s
*	load_eos:	Load the equation of state file
*	make_center:	Calculate the central pressure and enthalpy
*	sphere:		Compute the metric of a spherical star 
*	TOV:		Integrates the Tolman-Oppenheimer-Volkoff
*			equations for spherically symmetric star
*	spin:		Integrates the equations for a rapidly rotating
*			neutron star with oblateness = r_ratio = 
*				radius of pole/radius of equator
*	mass_radius:	Calculates the gravitational mass and equatorial
*			radius of the rotating star, along with other
*			equilibrium quantities. 
*
******************************************************************************/



#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <stdlib.h>
#include "equil_util.h"
#include "consts.h"
#include "nrutil.h"
#include "equil.h"
#include <stdbool.h>

//#define C 2.9979e10                  /* speed of light in vacuum */
//#define G 6.6732e-8                  /* gravitational constant */ 
//#define KAPPA 1.346790806509621e+13  /* square of length scale = 1e-15*C*C/G */
//#define KSCALE 1.112668301525780e-36 /* KAPPA*G/(C*C*C*C) */  
//#define MSUN 1.987e33                /* Mass of Sun */
//#define PI 3.1415926535  
//#define Out_cond 4
/*******************************************************************/
/* Create computational grid.                                      */
/* Points in the mu-direction are stored in the array mu[i].       */
/* Points in the s-direction are stored in the array s_gp[j].      */
/*******************************************************************/

double max(double v1,double v2){
     if(v1>=v2){
     return v1;
     }
     else{ 
     return v2;
     }
}

double min(double v1,double v2){
     if(v1>v2){
     return v2;
     }
     else{
     return v1;
     }
}


void make_grid(double s_gp[SDIV+1], 
               double mu[MDIV+1])                        
{ 
  int m, s;                         /* counters */
    
      for(s=1;s<=SDIV;s++) 
         s_gp[s] = SMAX*(s-1.0)/(SDIV-1.0);

	/* s_gp[1] = 0.0     corresponds to the center of the star
	   s_gp[SDIV] = SMAX corresponds to infinity */

	/* SMAX is defined in the file consts.h */

      for(m=1;m<=MDIV;m++) 
         mu[m] = (m-1.0)/(MDIV-1.0);

	/* mu[1] = 0.0    corresponds to the plane of the equator 
	   mu[MDIV] = 1.0 corresponds to the axis of symmetry */

	/* s_gp[0] and mu[0] are not used by the program */

}



/*************************************************************************/
/* Load EOS file.                                                        */ 
/*************************************************************************/
void load_eos( char eos_file[], 
               double log_e_tab[2001], 
               double log_p_tab[2001], 
               double log_h_tab[2001],
               double log_n0_tab[2001], 
               int *n_tab)
{
 int i;                    /* counter */

 double p,                 /* pressure */
        rho,               /* density */
        h,                 /* enthalpy */
        n0;                /* number density */    
        //g;                 /* Gamma */

 FILE *f_eos;              /* pointer to eos_file */
  

    /* OPEN FILE TO READ */

    if((f_eos=fopen(eos_file,"r")) == NULL ) {    
       printf("cannot open file:  %s\n",eos_file); 
       exit(0);
    }

 
    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_eos,"%d\n",n_tab);


    /* READ EOS, H, N0 AND MAKE THEM DIMENSIONLESS */
 
    for(i=1;i<=(*n_tab);i++) {  
      /*fscanf(f_eos,"%lf %lf %lf %lf %lf\n",&rho,&p,&h,&n0,&g) ; */
       fscanf(f_eos,"%lf %lf %lf %lf\n",&rho,&p,&h,&n0) ;
       log_e_tab[i]=log10(rho*C*C*KSCALE);     /* multiply by C^2 to get */ 
       log_p_tab[i]=log10(p*KSCALE);           /* energy density. */
       log_h_tab[i]=log10(h/(C*C));        
       log_n0_tab[i]=log10(n0);
       /*Gamma_tab[i]=g;*/
    }
}

/*************************************************************************/
/* DM EOS.                                                        */
/*************************************************************************/
double epsilon_D_EOS(char eos_typeDM[1],
                     double m_chi,
                     double y_chi,
                     double x_D)
{
    return (pow(m_chi, 4.0)/(pow(HBAR, 3.0)*pow(PI, 2.0)))*((1.0/8.0)*((2.0*pow(x_D, 3.0) + x_D)*sqrt(1.0 + pow(x_D, 2.0)) - asinh(x_D)) + pow(y_chi, 2.0)*pow(x_D, 6.0)/(9.0*pow(PI, 2.0)));
}

double P_D_EOS(char eos_typeDM[1],
               double m_chi,
               double y_chi,
               double x_D)
{
    return (pow(m_chi, 4.0)/(3.0*pow(HBAR, 3.0)*pow(PI, 2.0)))*((1.0/8.0)*((2.0*pow(x_D, 3.0) - 3.0*x_D)*sqrt(1.0 + pow(x_D, 2.0)) + 3.0*asinh(x_D)) + pow(y_chi, 2.0)*pow(x_D, 6.0)/(3.0*pow(PI, 2.0)));
}

double h_D_EOS(char eos_typeDM[1],
               double m_chi,
               double y_chi,
               double x_D)
{
    return log(sqrt(1.0 + pow(x_D, 2.0)) + 2.0*pow(x_D, 3.0)*pow(y_chi, 2.0)/(3.0*pow(PI, 2.0)));
}




/*******************************************************************/
double e_of_rho0(double rho0, double Gamma_P)
{
 return(pow(rho0,Gamma_P)/(Gamma_P-1.0)+rho0);
}
   

/*C*/
/*******************************************************************/
double e_at_p(double pp, 
              double pp_surface,
              double log_e_tab[2001], 
              double log_p_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt,
              char eos_type[],
              double Gamma_P)
{
 if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"DM")==0)){
   if(pp<pp_surface){
     return 0;
   }else{
     return pow(10.0,interp(log_p_tab,log_e_tab,n_tab,log10(pp), n_nearest_pt));
   }
 }
}

/*C*/
/*******************************************************************/
double e_at_p_DM(double pp,
                 double pp_surface,
                 char eos_typeDM[1],
                 double m_chi,
                 double y_chi)
{
    if(pp<pp_surface){
        return 0;
    }else{
        double int_end = X_D_MAX;
        double int_st = X_D_MIN;
        double root = int_st;
        
        double mid_pt = (int_end + int_st) / 2;
        double mid_pt_old = mid_pt;
          
        while (fabs (root - mid_pt_old) / mid_pt_old >= 1.0e-15)
        {
            mid_pt_old = mid_pt;
            // chcck initial num * mid_pt is less than 0
            if ( (P_D_EOS(eos_typeDM, m_chi, y_chi, int_st) - pp) * (P_D_EOS(eos_typeDM, m_chi, y_chi, mid_pt) - pp) < 0)
            {
                int_end = mid_pt; // assign the mid_pt to int_end
            }
            else
            {
                int_st = mid_pt; // else it assign the mid_pt to int_st
            }
            
            mid_pt = (int_end + int_st) / 2;
            root = mid_pt;
        }
        
        return epsilon_D_EOS(eos_typeDM, m_chi, y_chi, root);
    }
}

/*C*/
/*******************************************************************/
double e_at_h(double hh,
              double hh_surface,
              double log_e_tab[2001],
              double log_h_tab[2001],
              int    n_tab,
              int    *n_nearest_pt)
{

 if(hh<hh_surface){
  return 0;
 }else{
 return pow(10.0,interp(log_h_tab,log_e_tab,n_tab,log10(hh), n_nearest_pt));
 }
}

/*C*/
/*******************************************************************/
double e_at_h_DM(double hh,
                 double hh_surface,
                 char eos_typeDM[1],
                 double m_chi,
                 double y_chi)
{
 if(hh<hh_surface){
  return 0;
 }else{
     double int_end = X_D_MAX;
     double int_st = X_D_MIN;
     double root = int_st;
     
     double mid_pt = (int_end + int_st) / 2;
     double mid_pt_old = mid_pt;
       
     while (fabs (root - mid_pt_old) / mid_pt_old >= 1.0e-15)
     {
         mid_pt_old = mid_pt;
         // chcck initial num * mid_pt is less than 0
         if ( (h_D_EOS(eos_typeDM, m_chi, y_chi, int_st) - hh) * (h_D_EOS(eos_typeDM, m_chi, y_chi, mid_pt) - hh) < 0)
         {
             int_end = mid_pt; // assign the mid_pt to int_end
         }
         else
         {
             int_st = mid_pt; // else it assign the mid_pt to int_st
         }
           
         mid_pt = (int_end + int_st) / 2;
         root = mid_pt;
     }
     
     return epsilon_D_EOS(eos_typeDM, m_chi, y_chi, root);
 }
}

/*C*/
/*******************************************************************/
double p_at_e(double ee,
              double ee_surface, 
              double log_p_tab[2001], 
              double log_e_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt)
{ 
 if(ee<ee_surface){
  return 0;
 }else{
 return pow(10.0,interp(log_e_tab,log_p_tab,n_tab,log10(ee), n_nearest_pt));
 }
} 

/*C*/
/*******************************************************************/
double p_at_e_DM(double ee,
                 double ee_surface,
                 char eos_typeDM[1],
                 double m_chi,
                 double y_chi)
{
 if(ee<ee_surface){
  return 0;
 }else{
     double int_end = X_D_MAX;
     double int_st = X_D_MIN;
     double root = int_st;
     
     double mid_pt = (int_end + int_st) / 2;
     double mid_pt_old = mid_pt;
       
     while (fabs (root - mid_pt_old) / mid_pt_old >= 1.0e-15)
     {
         mid_pt_old = mid_pt;
         // chcck initial num * mid_pt is less than 0
         if ( (epsilon_D_EOS(eos_typeDM, m_chi, y_chi, int_st) - ee) * (epsilon_D_EOS(eos_typeDM, m_chi, y_chi, mid_pt) - ee) < 0)
         {
             int_end = mid_pt; // assign the mid_pt to int_end
         }
         else
         {
             int_st = mid_pt; // else it assign the mid_pt to int_st
         }
           
         mid_pt = (int_end + int_st) / 2;
         root = mid_pt;
     }
     
     return P_D_EOS(eos_typeDM, m_chi, y_chi, root);
 }
}

/*C*/
/*******************************************************************/
double p_at_h(double hh,
              double hh_surface, 
              double log_p_tab[2001], 
              double log_h_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt)
{

 if(hh<hh_surface){
  return 0;
 }else{
 return pow(10.0,interp(log_h_tab,log_p_tab,n_tab,log10(hh), n_nearest_pt));
 }
}

/*C*/
/*******************************************************************/
double p_at_h_DM(double hh,
                 double hh_surface,
                 char eos_typeDM[1],
                 double m_chi,
                 double y_chi)
{
 if(hh<hh_surface){
  return 0;
 }else{
     double int_end = X_D_MAX;
     double int_st = X_D_MIN;
     double root = int_st;
     
     double mid_pt = (int_end + int_st) / 2;
     double mid_pt_old = mid_pt;
       
     while (fabs (root - mid_pt_old) / mid_pt_old >= 1.0e-15)
     {
         mid_pt_old = mid_pt;
         // chcck initial num * mid_pt is less than 0
         if ( (h_D_EOS(eos_typeDM, m_chi, y_chi, int_st) - hh) * (h_D_EOS(eos_typeDM, m_chi, y_chi, mid_pt) - hh) < 0)
         {
             int_end = mid_pt; // assign the mid_pt to int_end
         }
         else
         {
             int_st = mid_pt; // else it assign the mid_pt to int_st
         }
           
         mid_pt = (int_end + int_st) / 2;
         root = mid_pt;
     }
     
     return P_D_EOS(eos_typeDM, m_chi, y_chi, root);
 }
}

/*C*/
/*******************************************************************/
double h_at_p(double pp, 
              double pp_surface,
              double log_h_tab[2001], 
              double log_p_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt)
{

 if(pp<pp_surface){
  return 0;
 }else{
 return pow(10.0,interp(log_p_tab,log_h_tab,n_tab,log10(pp), n_nearest_pt));
 }
}

/*C*/
/*******************************************************************/
double h_at_p_DM(double pp,
                 double pp_surface,
                 char eos_typeDM[1],
                 double m_chi,
                 double y_chi)
{
 if(pp<pp_surface){
  return 0;
 }else{
     double int_end = X_D_MAX;
     double int_st = X_D_MIN;
     double root = int_st;
     
     double mid_pt = (int_end + int_st) / 2;
     double mid_pt_old = mid_pt;
       
     while (fabs (root - mid_pt_old) / mid_pt_old >= 1.0e-15)
     {
         mid_pt_old = mid_pt;
           
         // chcck initial num * mid_pt is less than 0
         if ( (P_D_EOS(eos_typeDM, m_chi, y_chi, int_st) - pp) * (P_D_EOS(eos_typeDM, m_chi, y_chi, mid_pt) - pp) < 0)
         {
             int_end = mid_pt; // assign the mid_pt to int_end
         }
         else
         {
             int_st = mid_pt; // else it assign the mid_pt to int_st
         }
           
         mid_pt = (int_end + int_st) / 2;
         root = mid_pt;
     }

     return h_D_EOS(eos_typeDM, m_chi, y_chi, root);
 }
}

/*C*/
/*******************************************************************/
double n0_at_e(double ee, 
	       double ee_surface,
               double log_n0_tab[2001], 
               double log_e_tab[2001],
               int    n_tab, 
               int    *n_nearest_pt)
{
 if(ee<ee_surface){
  return 0;
 }else{
  return pow(10.0,interp(log_e_tab,log_n0_tab,n_tab,log10(ee), n_nearest_pt));
 }
}
 
/*C*/
/***************************************************************/
void make_center(
	       char eos_file[], 
               double log_e_tab[2001], 
               double log_p_tab[2001], 
               double log_h_tab[2001],
               double log_n0_tab[2001], 
               int n_tab,                 
	       char eos_type[],
	       double Gamma_P, 
	       double e_center,
	       double *p_center, 
	       double *h_center,
	       double e_surface,
	       double p_surface)

{

 int n_nearest;

 double rho0_center;

 n_nearest=n_tab/2; 


   (*p_center) = p_at_e( e_center, e_surface, log_p_tab, log_e_tab, n_tab, &n_nearest);
   (*h_center) = h_at_p( (*p_center), p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest);
    
}

/*C*/
/***************************************************************/
void make_centerDM(
           char eos_typeDM[1],
           double m_chi,
           double y_chi,
           double e_centerDM,
           double *p_centerDM,
           double *h_centerDM,
           double e_surfaceDM,
           double p_surfaceDM)

{
    
   (*p_centerDM) = p_at_e_DM( e_centerDM, e_surfaceDM, eos_typeDM, m_chi, y_chi);
   (*h_centerDM) = h_at_p_DM( (*p_centerDM), p_surfaceDM, eos_typeDM, m_chi, y_chi);

}

/*C*/
/***********************************************************************/
/* Computes the gravitational mass, equatorial radius, angular momentum
 *	of the star
 * 	and the velocity of co- and counter-rotating particles      
 *	with respect to a ZAMO                                         */
/***********************************************************************/
void mass_radius(
		 double s_gp[SDIV+1],
		 double mu[MDIV+1],
		 double log_e_tab[2001], 
		 double log_p_tab[2001], 
		 double log_h_tab[2001],
		 double log_n0_tab[2001], 
		 int n_tab,                 
		 char eos_type[],
		 double Gamma_P,
		 char eos_typeDM[1],
         double m_chi,
         double y_chi,
		 double **rho,
		 double **gama,
		 double **alpha,
		 double **omega, //20
		 double **energy,
		 double **pressure,
		 double **enthalpy,
		 double **velocity_sq,
		 double **energyDM,
		 double **pressureDM,
		 double **enthalpyDM,
		 double **velocity_sqDM,
                 double **Omega_hDM,
                 double r_ratio,
                 double r_ratioDM, //30
                 double *Ratio_sch,
                 double e_center,
                 double e_centerDM,
		 double e_surface,
		 double e_surfaceDM,
                 double r_e,
                 double r_eDM,
                 double Omega,
                 double OmegaDM,
                 double *Mass,
		 double *Mass_0,
		 double *ang_mom,
                 double *R_e, //40
                 double *MassDM, 
		 double *Mass_0DM,
		 double *ang_momDM,
                 double *R_eDM,
		 double *v_plus,
		 double *v_minus,
		 double *Omega_K,
		 double *Vp,
		 double *Mp) //49

{


 int s,
     m,
     n_nearest,
     n_nearestDM;

int index;

 
 double   
   **rho_0, /*rest mass density*/
   **velocity,
   gama_equator,              /* gama at equator */
   rho_equator,
    gama_pole,
    rho_pole,
   gama_equatorDM,              /* gama at equator */
   rho_equatorDM,               /* rho at equator */
    gama_poleDM,
    rho_poleDM,
   omega_equator,             /* omega at equator */
   s1,
   s_1,r_eDM_old,
   d_gama_s,
   d_rho_s,
   d_omega_s,
   sqrt_v,
   D_m[SDIV+1],               /* int. quantity for M */
    D_mRB[SDIV+1],               /* int. quantity for M */
   D_m_0[SDIV+1],             /* int. quantity for M_0 */
   D_J[SDIV+1],               /* int. quantity for J */
   D_mDM[SDIV+1],               /* int. quantity for M */
    D_mDMRB[SDIV+1],               /* int. quantity for M */
   D_m_0DM[SDIV+1],             /* int. quantity for M_0 */
   D_JDM[SDIV+1],               /* int. quantity for J */
    D_s[MDIV+1],               /* int. quantity for M */
    D_sDM[MDIV+1],               /* int. quantity for M */
   s_e,
   s_eDM,                 
   d_o_e[SDIV+1],
   d_g_e[SDIV+1],
   d_r_e[SDIV+1],
   d_v_e[SDIV+1],
   doe,
   dge, 
   dre,
   dve,
   vek, 
   ratio_old,ratio_oldDM,    
   gama_mu_0[SDIV+1],                   
   rho_mu_0[SDIV+1], 
   gama_mu_1[SDIV+1],                   
   rho_mu_1[SDIV+1],                      
   omega_mu_0[SDIV+1],
   J,
   JDM,
   r_p,
    r_pDM,
   s_p,
   s_pDM,
   D_vp[SDIV+1],
   D_mp[SDIV+1],
   Rv,r_eq,r_out;        
   r_eDM_old=r_eDM;
   /* Circumferential radius */
 
   r_eq=r_e;
   switch(Out_cond){
     case 1:
        r_out= ((r_eq>=r_eDM)?r_eq:r_eDM);
        break;
     case 2:
        r_out= ((r_eq<=r_eDM)?r_eq:r_eDM);
        break; 
     case 3:
        r_out= r_eq;//((r_eq>=r_eDM)?r_eq:r_eDM);
        break;
     case 4:
        r_out= r_eDM;//((r_eq<=r_eDM)?r_eq:r_eDM);
        break;   
     case 5:
        r_out= (r_eq+r_eDM)/2.;
        break;    
   }
   
   
   
   r_p= r_ratio*r_eq;                              /* radius at pole */
   s_p=r_ratio*r_eq/(r_ratio*r_eq+r_out);            
   s_e=r_eq/(r_eq+r_out);    


//   s_pDM=r_ratioDM*r_eDM/(r_ratioDM*r_eDM+r_out);            
   s_eDM=r_eDM/(r_eDM+r_out);            

       
   rho_0 = dmatrix(1,SDIV,1,MDIV);
   velocity = dmatrix(1,SDIV,1,MDIV);
    
    int s_pDM_found = 0;

   for(s=1;s<=SDIV;s++) {               
      gama_mu_0[s]=gama[s][1];                   
      rho_mu_0[s]=rho[s][1];                                                    
      gama_mu_1[s]=gama[s][MDIV];                   
      rho_mu_1[s]=rho[s][MDIV];
       
       if(enthalpyDM[s][MDIV] < 1.0/(C*C) && s_pDM_found == 0)
       {
           s_pDM = s_gp[s];
           s_pDM_found += 1;
       }
   }
    if(r_eDM != 0.0)
    {
        r_pDM = r_out*s_pDM/(1.0 - s_pDM);
    } else
    {
        r_pDM = 0.0;
    }
    

   n_nearest= SDIV/2;
   gama_equator=interp(s_gp,gama_mu_0,SDIV,s_e, &n_nearest);
   rho_equator=interp(s_gp,rho_mu_0,SDIV,s_e, &n_nearest);
    gama_pole = interp(s_gp, gama_mu_1, SDIV, s_p, &n_nearest);
    rho_pole = interp(s_gp, rho_mu_1, SDIV, s_p, &n_nearest);

   n_nearestDM= SDIV/2;
   gama_equatorDM=interp(s_gp,gama_mu_0,SDIV,s_eDM, &n_nearestDM);  
   rho_equatorDM=interp(s_gp,rho_mu_0,SDIV,s_eDM, &n_nearestDM);
    gama_poleDM = interp(s_gp, gama_mu_1, SDIV, s_pDM, &n_nearestDM);
    rho_poleDM = interp(s_gp, rho_mu_1, SDIV, s_pDM, &n_nearestDM);

/* Circumferential radius */

   (*R_e) = sqrt(KAPPA)*r_eq*exp((gama_equator-rho_equator)/2.0);
   (*R_eDM) =  sqrt(KAPPA)*r_eDM*exp((gama_equatorDM-rho_equatorDM)/2.0);
    double R_p = sqrt(KAPPA)*r_p*exp((gama_pole - rho_pole)/2.0);
    double R_pDM = sqrt(KAPPA)*r_pDM*exp((gama_poleDM - rho_poleDM)/2.0);
   (*Ratio_sch)=R_p/(*R_e);


 /* Masses and angular momentum */
    
    double Mass_RB;
    double MassDM_RB;
    
    Mass_RB = 0.0;              /* initialize */
    MassDM_RB = 0.0;              /* initialize */
 
   (*Mass) = 0.0;              /* initialize */
   (*Mass_0) = 0.0;
   (*MassDM) = 0.0;              /* initialize */
   (*Mass_0DM) = 0.0;
   (*Vp) = 0.0;
   (*Mp) = 0.0;
   J=0.0;
   JDM=0.0;
   Rv=0.0;
   /* CALCULATE THE REST MASS DENSITY */
 if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"DM")==0)) {
   n_nearest=n_tab/2;
   for(s=1;s<=SDIV;s++)
      for(m=1;m<=MDIV;m++) {
           if(energy[s][m]>e_surface)
             rho_0[s][m]=n0_at_e(energy[s][m],e_surface, log_n0_tab, log_e_tab, n_tab,
                                             &n_nearest)*MB*KSCALE*SQ(C);
           else
             rho_0[s][m]=0.0;
      }  
 }

   for(s=1;s<=SDIV;s++) {
    D_m[s]=0.0;           /* initialize */
    D_mDM[s]=0.0;           /* initialize */
       D_mRB[s]=0.0;           /* initialize */
       D_mDMRB[s]=0.0;           /* initialize */
    D_m_0[s]=0.0;
    D_J[s]=0.0;
    D_JDM[s]=0.0;
    D_vp[s]=0.0;
    D_mp[s]=0.0;
    
    for(m=1;m<=MDIV-2;m+=2) {

     D_m[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+gama[s][m])*
              (((energy[s][m]+pressure[s][m])/(1.0-velocity_sq[s][m]))*
              (1.0+velocity_sq[s][m]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_out*omega[s][m]*
              exp(-rho[s][m])) + 2.0*pressure[s][m])

            + 4.0*exp(2.0*alpha[s][m+1]+gama[s][m+1])*
              (((energy[s][m+1]+pressure[s][m+1])/(1.0-velocity_sq[s][m+1]))*
              (1.0+velocity_sq[s][m+1]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+1])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_out*omega[s][m+1]*
              exp(-rho[s][m+1])) + 2.0*pressure[s][m+1]) 

            + exp(2.0*alpha[s][m+2]+gama[s][m+2])*
              (((energy[s][m+2]+pressure[s][m+2])/(1.0-velocity_sq[s][m+2]))*
              (1.0+velocity_sq[s][m+2]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+2])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_out*omega[s][m+2]*
              exp(-rho[s][m+2])) + 2.0*pressure[s][m+2]));


    D_mDM[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+gama[s][m])*
              (((energyDM[s][m]+pressureDM[s][m])/(1.0-velocity_sqDM[s][m]))*
              (1.0+velocity_sqDM[s][m]+(2.0*s_gp[s]*sqrt(velocity_sqDM[s][m])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_out*omega[s][m]*
              exp(-rho[s][m])) + 2.0*pressureDM[s][m])

            + 4.0*exp(2.0*alpha[s][m+1]+gama[s][m+1])*
              (((energyDM[s][m+1]+pressureDM[s][m+1])/(1.0-velocity_sqDM[s][m+1]))*
              (1.0+velocity_sqDM[s][m+1]+(2.0*s_gp[s]*sqrt(velocity_sqDM[s][m+1])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_out*omega[s][m+1]*
              exp(-rho[s][m+1])) + 2.0*pressureDM[s][m+1]) 

            + exp(2.0*alpha[s][m+2]+gama[s][m+2])*
              (((energyDM[s][m+2]+pressureDM[s][m+2])/(1.0-velocity_sqDM[s][m+2]))*
              (1.0+velocity_sqDM[s][m+2]+(2.0*s_gp[s]*sqrt(velocity_sqDM[s][m+2])/
              (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_out*omega[s][m+2]*
              exp(-rho[s][m+2])) + 2.0*pressureDM[s][m+2]));
              //printf("%d %d %f %f %f\n",s,m, energyDM[s][m+2]+pressureDM[s][m+2], energyDM[s][m+1]+pressureDM[s][m+1],energyDM[s][m]+pressureDM[s][m]);
        
        if(enthalpy[s][m] > 0.0)
        {
            D_mRB[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+gama[s][m])*
                     (((energy[s][m]+pressure[s][m])/(1.0-velocity_sq[s][m]))*
                     (1.0+velocity_sq[s][m]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m])/
                     (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_out*omega[s][m]*
                     exp(-rho[s][m])) + 2.0*pressure[s][m])

                   + 4.0*exp(2.0*alpha[s][m+1]+gama[s][m+1])*
                     (((energy[s][m+1]+pressure[s][m+1])/(1.0-velocity_sq[s][m+1]))*
                     (1.0+velocity_sq[s][m+1]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+1])/
                     (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_out*omega[s][m+1]*
                     exp(-rho[s][m+1])) + 2.0*pressure[s][m+1])

                   + exp(2.0*alpha[s][m+2]+gama[s][m+2])*
                     (((energy[s][m+2]+pressure[s][m+2])/(1.0-velocity_sq[s][m+2]))*
                     (1.0+velocity_sq[s][m+2]+(2.0*s_gp[s]*sqrt(velocity_sq[s][m+2])/
                     (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_out*omega[s][m+2]*
                     exp(-rho[s][m+2])) + 2.0*pressure[s][m+2]));
            
            D_mDMRB[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+gama[s][m])*
                      (((energyDM[s][m]+pressureDM[s][m])/(1.0-velocity_sqDM[s][m]))*
                      (1.0+velocity_sqDM[s][m]+(2.0*s_gp[s]*sqrt(velocity_sqDM[s][m])/
                      (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_out*omega[s][m]*
                      exp(-rho[s][m])) + 2.0*pressureDM[s][m])

                    + 4.0*exp(2.0*alpha[s][m+1]+gama[s][m+1])*
                      (((energyDM[s][m+1]+pressureDM[s][m+1])/(1.0-velocity_sqDM[s][m+1]))*
                      (1.0+velocity_sqDM[s][m+1]+(2.0*s_gp[s]*sqrt(velocity_sqDM[s][m+1])/
                      (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_out*omega[s][m+1]*
                      exp(-rho[s][m+1])) + 2.0*pressureDM[s][m+1])

                    + exp(2.0*alpha[s][m+2]+gama[s][m+2])*
                      (((energyDM[s][m+2]+pressureDM[s][m+2])/(1.0-velocity_sqDM[s][m+2]))*
                      (1.0+velocity_sqDM[s][m+2]+(2.0*s_gp[s]*sqrt(velocity_sqDM[s][m+2])/
                      (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_out*omega[s][m+2]*
                      exp(-rho[s][m+2])) + 2.0*pressureDM[s][m+2]));
            
        }
 


     D_m_0[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+(gama[s][m]
              -rho[s][m])/2.0)*rho_0[s][m]/sqrt(1.0-velocity_sq[s][m])

             + 4.0* exp(2.0*alpha[s][m+1]+(gama[s][m+1]
             -rho[s][m+1])/2.0)*rho_0[s][m+1]/sqrt(1.0-velocity_sq[s][m+1])
         
             + exp(2.0*alpha[s][m+2]+(gama[s][m+2]
             -rho[s][m+2])/2.0)*rho_0[s][m+2]/sqrt(1.0-velocity_sq[s][m+2]));

     D_J[s] += (1.0/(3.0*(MDIV-1)))*( sqrt(1.0-mu[m]*mu[m])*
              exp(2.0*alpha[s][m]+gama[s][m]-rho[s][m])*(energy[s][m]
              +pressure[s][m])*sqrt(velocity_sq[s][m])/(1.0-velocity_sq[s][m])
  
              +4.0*sqrt(1.0-mu[m+1]*mu[m+1])*
              exp(2.0*alpha[s][m+1]+gama[s][m+1]-rho[s][m+1])*(energy[s][m+1]
              +pressure[s][m+1])*sqrt(velocity_sq[s][m+1])/
              (1.0-velocity_sq[s][m+1])

              + sqrt(1.0-mu[m+2]*mu[m+2])*
              exp(2.0*alpha[s][m+2]+gama[s][m+2]-rho[s][m+2])*(energy[s][m+2]
              +pressure[s][m+2])*sqrt(velocity_sq[s][m+2])/
              (1.0-velocity_sq[s][m+2]));


     D_JDM[s] += (1.0/(3.0*(MDIV-1)))*( sqrt(1.0-mu[m]*mu[m])*
              exp(2.0*alpha[s][m]+gama[s][m]-rho[s][m])*(energyDM[s][m]
              +pressureDM[s][m])*sqrt(velocity_sqDM[s][m])/(1.0-velocity_sqDM[s][m])
  
              +4.0*sqrt(1.0-mu[m+1]*mu[m+1])*
              exp(2.0*alpha[s][m+1]+gama[s][m+1]-rho[s][m+1])*(energyDM[s][m+1]
              +pressureDM[s][m+1])*sqrt(velocity_sqDM[s][m+1])/
              (1.0-velocity_sqDM[s][m+1])

              + sqrt(1.0-mu[m+2]*mu[m+2])*
              exp(2.0*alpha[s][m+2]+gama[s][m+2]-rho[s][m+2])*(energyDM[s][m+2]
              +pressureDM[s][m+2])*sqrt(velocity_sqDM[s][m+2])/
              (1.0-velocity_sqDM[s][m+2]));

     D_mp[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+(gama[s][m]
              -rho[s][m])/2.0)*((energy[s][m]+pressure[s][m])/sqrt(1.0-velocity_sq[s][m]))

             + 4.0* exp(2.0*alpha[s][m+1]+(gama[s][m+1]
             -rho[s][m+1])/2.0)*((energy[s][m+1]+pressure[s][m+1])/sqrt(1.0-velocity_sq[s][m+1]))
         
             + exp(2.0*alpha[s][m+2]+(gama[s][m+2]
             -rho[s][m+2])/2.0)*((energy[s][m+2]+pressure[s][m+2])/sqrt(1.0-velocity_sq[s][m+2])));

     if(energy[s][m+1]+pressure[s][m+1]==0){
      D_vp[s]+=0.0;
     }else{
 
     D_vp[s] += (1.0/(3.0*(MDIV-1)))*( exp(2.0*alpha[s][m]+(gama[s][m]
              -rho[s][m])/2.0)*(1.0/sqrt(1.0-velocity_sq[s][m]))

             + 4.0* exp(2.0*alpha[s][m+1]+(gama[s][m+1]
             -rho[s][m+1])/2.0)*(1.0/sqrt(1.0-velocity_sq[s][m+1]))
         
             + exp(2.0*alpha[s][m+2]+(gama[s][m+2]
             -rho[s][m+2])/2.0)*(1.0/sqrt(1.0-velocity_sq[s][m+2])));
             Rv= s_gp[s]*r_out*sqrt(KAPPA)/((1.-s_gp[s])*100000.);
    }    

    }
   }
   index = 1;
    for(s=1;s<=SDIV-2;s+=2) { 
        //printf("%f %f %f ",log(velocity_sq[s][0]),energy[s][0],omega[s][0]);
     (*Mass) += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_m[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m[s+2]);

     (*MassDM) += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_mDM[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_mDM[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_mDM[s+2]);
        
        Mass_RB += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
             D_mRB[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_mRB[s+1]
             +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_mRB[s+2]);

        MassDM_RB += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
             D_mDMRB[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_mDMRB[s+1]
             +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_mDMRB[s+2]);
          
     //printf("%d \t %g %g %g\n", index, (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
      //    D_mDM[s],4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_mDM[s+1],pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_mDM[s+2])); 

     //printf("%d \t %g\n", index, (*Mass));      /***************************/

     (*Mass_0) += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_m_0[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m_0[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m_0[s+2]);
 
     J += (SMAX/(3.0*(SDIV-1)))*((pow(s_gp[s],3.0)/pow(1.0-s_gp[s],5.0))*
          D_J[s]+ 4.0*(pow(s_gp[s+1],3.0)/pow(1.0-s_gp[s+1],5.0))*
          D_J[s+1] + (pow(s_gp[s+2],3.0)/pow(1.0-s_gp[s+2],5.0))*
          D_J[s+2]);

     JDM += (SMAX/(3.0*(SDIV-1)))*((pow(s_gp[s],3.0)/pow(1.0-s_gp[s],5.0))*
          D_JDM[s]+ 4.0*(pow(s_gp[s+1],3.0)/pow(1.0-s_gp[s+1],5.0))*
          D_JDM[s+1] + (pow(s_gp[s+2],3.0)/pow(1.0-s_gp[s+2],5.0))*
          D_JDM[s+2]);

     (*Mp) += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_mp[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_mp[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_mp[s+2]);

     (*Vp) += (SMAX/(3.0*(SDIV-1)))*(pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
          D_vp[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_vp[s+1]
          +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_vp[s+2]);
  
     index++;
    }
    
    double M_cloud;
   
    if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"DM")==0)) {
      (*Mass) *= 4.*PI*sqrt(KAPPA)*C*C*pow(r_out,3.0)/G;
      (*MassDM) *= 4.*PI*sqrt(KAPPA)*C*C*pow(r_out,3.0)/G;
        Mass_RB *= 4.*PI*sqrt(KAPPA)*C*C*pow(r_out,3.0)/G;
        MassDM_RB *= 4.*PI*sqrt(KAPPA)*C*C*pow(r_out,3.0)/G;
      (*Mass_0) *= 4.*PI*sqrt(KAPPA)*C*C*pow(r_out,3.0)/G;
      (*Vp) *= 4.*PI*sqrt(KAPPA)*sqrt(KAPPA)*sqrt(KAPPA)/3.0;
      (*Mp) *= 4.*PI*sqrt(KAPPA)*C*C*pow(r_out,3.0)/G;

      printf("mass_radius: MB = %lf  MD = %lf \n", *Mass/MSUN, *MassDM/MSUN);
      printf("M_tot = %lf\n", (*Mass+*MassDM)/MSUN);
        printf("mass_radius: MB(RB) = %lf  MD(RB) = %lf \n", Mass_RB/MSUN, MassDM_RB/MSUN);
        printf("M_tot(RB) = %lf\n", (Mass_RB+MassDM_RB)/MSUN);
        
        M_cloud = *Mass + *MassDM - Mass_RB - MassDM_RB;
        
        printf("M_cloud = %lf\n", M_cloud/MSUN);
        
        printf("RBe_Sch = %lf\n", (*R_e)*1e-5);
        printf("RBp_Sch = %lf\n", R_p*1e-5);
        printf("RDe_Sch = %lf\n", (*R_eDM)*1e-5);
        printf("RDp_Sch = %lf\n", R_pDM*1e-5);
        
        printf("Omega/(2pi) = %lf\n", Omega/(2.0*PI));
        printf("OmegaDM/(2pi) = %lf\n", OmegaDM/(2.0*PI));
  
    }
     if(isnan((*Mass))&&isnan((*MassDM))){
       printf("The system is unstable\n");
       exit(0);
       return;
     }
    if(r_ratio==1.0) 
         J=0.0; 
    else {    
          if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"poly")==0))
              J *= 4.0*PI*KAPPA*C*C*C*pow(r_out,4.0)/G;
    }

    (*ang_mom) = J;

    if(r_ratioDM==1.0) 
         JDM=0.0; 
    else {    
          if((strcmp(eos_type,"b")==0) || (strcmp(eos_type,"f")==0))
              JDM *= 4.0*PI*KAPPA*C*C*C*pow(r_out,4.0)/G;
    }

    (*ang_momDM) = JDM;

    //printf(" J = %g \n", J);


  /* Compute the velocities of co-rotating and counter-rotating particles
	with respect to a ZAMO 	*/

  for(s=1+(SDIV-1)/2;s<=SDIV;s++) {
    s1= s_gp[s]*(1.0-s_gp[s]);
    s_1=1.0-s_gp[s];
        
    d_gama_s=deriv_s(gama,s,1);
    d_rho_s=deriv_s(rho,s,1);
    d_omega_s=deriv_s(omega,s,1);

    sqrt_v= exp(-2.0*rho[s][1])*r_out*r_out*pow(s_gp[s],4.0)*pow(d_omega_s,2.0) 
            + 2.0*s1*(d_gama_s+d_rho_s)+s1*s1*(d_gama_s*d_gama_s-d_rho_s*d_rho_s);

    if(sqrt_v>0.0) sqrt_v= sqrt(sqrt_v);
     else {
      sqrt_v=0.0;
     }

    v_plus[s]=(exp(-rho[s][1])*r_out*s_gp[s]*s_gp[s]*d_omega_s + sqrt_v)/
              (2.0+s1*(d_gama_s-d_rho_s));

    v_minus[s]=(exp(-rho[s][1])*r_out*s_gp[s]*s_gp[s]*d_omega_s - sqrt_v)/
               (2.0+s1*(d_gama_s-d_rho_s));
  }


/* Kepler angular velocity */

   for(s=1;s<=SDIV;s++) { 
     d_o_e[s]=deriv_s(omega,s,1);
     d_g_e[s]=deriv_s(gama,s,1);
     d_r_e[s]=deriv_s(rho,s,1);
     d_v_e[s]=deriv_s(velocity,s,1);
     /* Value of omega on the equatorial plane*/
     omega_mu_0[s] = omega[s][1];
   }

   n_nearest=SDIV/2;
   doe=interp(s_gp,d_o_e,SDIV,0.5, &n_nearest);
   dge=interp(s_gp,d_g_e,SDIV,0.5, &n_nearest);
   dre=interp(s_gp,d_r_e,SDIV,0.5, &n_nearest);
   dve=interp(s_gp,d_v_e,SDIV,0.5, &n_nearest);

  vek=(doe/(8.0+dge-dre))*r_out*exp(-rho_equator) + sqrt(((dge+dre)/(8.0+dge
        -dre)) + pow((doe/(8.0+dge-dre))*r_out*exp(-rho_equator),2.0));


  //if (r_ratio ==1.0)
  //  omega_equator = 0.0;
  //else
    omega_equator = interp(s_gp,omega_mu_0,SDIV,0.5, &n_nearest);




   (*Omega_K) = (C/sqrt(KAPPA))*(omega_equator+vek*exp(rho_equator)/r_out);
    
    
    
    
    
    
    
//    double enthalpy_m[SDIV+1];
//    double enthalpyDM_m[SDIV+1];
//    double gama_m[SDIV+1];
//    double rho_m[SDIV+1];
//    double omega_m[SDIV+1];
    
    double s_RB[MDIV+1];
    double gama_RB[MDIV+1];
    double rho_RB[MDIV+1];
    double omega_RB[MDIV+1];
    
    double s_RD[MDIV+1];
    double gama_RD[MDIV+1];
    double rho_RD[MDIV+1];
    double omega_RD[MDIV+1];
    
    int RB_found;
    int RD_found;
    
    for(m=1;m<=MDIV;m++) {
        
        RB_found = 0;
        RD_found = 0;
        
        for(s=1;s<=SDIV;s++) {
            
            if(enthalpy[s][m] < 1.0/(C*C) && RB_found == 0)
            {
                s_RB[m] = s_gp[s];
                gama_RB[m] = gama[s][m];
                rho_RB[m] = rho[s][m];
                omega_RB[m] = omega[s][m];
                
                RB_found += 1;
            }
            
            if(enthalpyDM[s][m] < 1.0/(C*C) && RD_found == 0)
            {
                s_RD[m] = s_gp[s];
                gama_RD[m] = gama[s][m];
                rho_RD[m] = rho[s][m];
                omega_RD[m] = omega[s][m];
                
                RD_found += 1;
            }
            
            if(RB_found != 0 && RD_found != 0)
            {
                break;
            }
        }
        
//        printf("s_RB[%d] = %.12e, %d \n", m, s_RB[m]);
    }
    

   free_dmatrix(velocity,1,SDIV,1,MDIV);
   free_dmatrix(rho_0,1,SDIV,1,MDIV);

              //printf("%.5f\n",Rv);


}

/*C*/
/**************************************************************************/
double dm_dr_is(double r_is, 
                double r, 
                double m, 
                double p, 
                double e_center, 
                double p_surface,
                double log_e_tab[SDIV+1],
                double log_p_tab[SDIV+1],
                int    n_tab,
                int    *n_nearest_pt,
                char eos_type[],
                double Gamma_P)
{
 double dmdr,
        e_d;

 if(p<p_surface) 
    e_d=0.0;
 else  
    e_d = e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, n_nearest_pt, eos_type, 
                                                                  Gamma_P);
 
 if(r_is<RMIN)
    dmdr=4.0*PI*e_center*r*r*(1.0+4.0*PI*e_center*r*r/3.0);
 else
    dmdr=4.0*PI*e_d*r*r*r*sqrt(1.0-2.0*m/r)/r_is;
 
return dmdr;
}

/*C*/
/**************************************************************************/
double dm_dr(double r_is,
                double r,
                double m,
                double p,
                double e_center,
                double p_surface,
                double log_e_tab[SDIV+1],
                double log_p_tab[SDIV+1],
                int    n_tab,
                int    *n_nearest_pt,
                char eos_type[],
                double Gamma_P)
{
 double dmdr,
        e_d;

 if(p<p_surface)
    e_d=0.0;
 else
    e_d = e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, n_nearest_pt, eos_type,
                                                                  Gamma_P);
 
 if(r<RMIN)
    dmdr=4.0*PI*e_center*r*r;
 else
    dmdr=4.0*PI*e_d*r*r;
 
return dmdr;
}

/*C*/
/**************************************************************************/
double dm_dr_is_DM(double r_is,
                   double r,
                   double m,
                   double p,
                   double e_center,
                   double p_surface,
                   char eos_typeDM[1],
                   double m_chi,
                   double y_chi)
{
 double dmdr,
        e_d;

 if(p<p_surface)
    e_d=0.0;
 else
    e_d = e_at_p_DM(p,p_surface, eos_typeDM, m_chi, y_chi);
 
 if(r_is<RMIN)
    dmdr=4.0*PI*e_center*r*r*(1.0+4.0*PI*e_center*r*r/3.0);
 else
    dmdr=4.0*PI*e_d*r*r*r*sqrt(1.0-2.0*m/r)/r_is;
 
return dmdr;
}

/*C*/
/**************************************************************************/
double dm_dr_DM(double r_is,
                   double r,
                   double m,
                   double p,
                   double e_center,
                   double p_surface,
                   char eos_typeDM[1],
                   double m_chi,
                   double y_chi)
{
 double dmdr,
        e_d;

 if(p<p_surface)
    e_d=0.0;
 else
    e_d = e_at_p_DM(p,p_surface, eos_typeDM, m_chi, y_chi);
 
 if(r<RMIN)
    dmdr=4.0*PI*e_center*r*r;
 else
    dmdr=4.0*PI*e_d*r*r;
 
return dmdr;
}
 
/*C*/
/**************************************************************************/
double dp_dr_is(double r_is, 
                double r, 
                double m, 
                double p,
                double p_other, 
                double e_center,
                double e_center_other,
                double p_surface,
                double log_e_tab[SDIV+1],
                double log_p_tab[SDIV+1],
                int    n_tab,
                int    *n_nearest_pt,
                char eos_type[],
                double Gamma_P)
{ double dpdr,
         e_d; 

  if(p<p_surface){
   e_d=0.0;
   p=0.0; 
  }
  else        
   e_d=e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, n_nearest_pt, eos_type, 
                                                                  Gamma_P);
  
  if(r_is<RMIN) dpdr = -4.0*PI*(e_center+p)*(e_center + e_center_other +3.0*(p+p_other))*r*(1.0
                     +4.0*(e_center + e_center_other)*r*r/3.0)/3.0;

  else 
   dpdr = -(e_d+p)*(m+4.0*PI*r*r*r*(p+p_other))/(r*r_is*sqrt(1.0-2.0*m/r));

 return dpdr;
}

/*C*/
/**************************************************************************/
double dp_dr(double r_is,
                double r,
                double m,
                double p,
                double p_other,
                double e_center,
                double e_center_other,
                double p_surface,
                double log_e_tab[SDIV+1],
                double log_p_tab[SDIV+1],
                int    n_tab,
                int    *n_nearest_pt,
                char eos_type[],
                double Gamma_P)
{ double dpdr,
         e_d;

  if(p<p_surface){
   e_d=0.0;
   p=0.0;
  }
  else
   e_d=e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, n_nearest_pt, eos_type,
                                                                  Gamma_P);
  
  if(r<RMIN) dpdr = -4.0*PI*(e_center+p)*((e_center + e_center_other)/3.0 + p+p_other)*r/(1.0 - 8.0*PI*(e_center + e_center_other)*r*r/3.0);

  else
   dpdr = -(e_d+p)*(m+4.0*PI*r*r*r*(p+p_other))/(r*(r-2.0*m));

 return dpdr;
}

/*C*/
/**************************************************************************/
double dp_dr_is_DM(double r_is,
                   double r,
                   double m,
                   double p,
                   double p_other,
                   double e_center,
                   double e_center_other,
                   double p_surface,
                   char eos_typeDM[1],
                   double m_chi,
                   double y_chi)
{ double dpdr,
         e_d;

  if(p<p_surface){
   e_d=0.0;
   p=0.0;
  }
  else
   e_d=e_at_p_DM(p,p_surface, eos_typeDM, m_chi, y_chi);
  
  if(r_is<RMIN) dpdr = -4.0*PI*(e_center+p)*(e_center + e_center_other +3.0*(p+p_other))*r*(1.0
                     +4.0*(e_center + e_center_other)*r*r/3.0)/3.0;

  else
   dpdr = -(e_d+p)*(m+4.0*PI*r*r*r*(p+p_other))/(r*r_is*sqrt(1.0-2.0*m/r));

 return dpdr;
}

/*C*/
/**************************************************************************/
double dp_dr_DM(double r_is,
                   double r,
                   double m,
                   double p,
                   double p_other,
                   double e_center,
                   double e_center_other,
                   double p_surface,
                   char eos_typeDM[1],
                   double m_chi,
                   double y_chi)
{ double dpdr,
         e_d;

  if(p<p_surface){
   e_d=0.0;
   p=0.0;
  }
  else
   e_d=e_at_p_DM(p,p_surface, eos_typeDM, m_chi, y_chi);
  
  if(r<RMIN) dpdr = -4.0*PI*(e_center+p)*((e_center + e_center_other)/3.0 + p+p_other)*r/(1.0 - 8.0*PI*(e_center + e_center_other)*r*r/3.0);

  else
   dpdr = -(e_d+p)*(m+4.0*PI*r*r*r*(p+p_other))/(r*(r-2.0*m));

 return dpdr;
}

/**************************************************************************/
double dr_dr_is(double r_is, double r, double m)
{
 double drdris;

 if(r<RMIN) drdris=1.0;
//    if(r<RMIN) drdris=(r/r_is)*sqrt(1.0-8.0*PI*(e_center + e_centerDM)*r*r/3.0);
  else
   drdris=(r/r_is)*sqrt(1.0-2.0*m/r);

 return drdris;
}

/**************************************************************************/
double dr_is_dr(double r_is, double r, double m)
{
 return 1.0/dr_dr_is(r_is, r, m);
}

/*C*/
/**************************************************************************/
double dr2_dh(double r2, double m, double p)
{ double dr2dh;
    
    dr2dh = -2.0*r2*(sqrt(r2) - 2.0*m)/(m + 4.0*PI*pow(r2, 1.5)*p);

 return dr2dh;
}

/*C*/
/**************************************************************************/
double dm_dh(double r2, double m, double p, double e)
{ double dmdh;
    
    dmdh = 2.0*PI*sqrt(r2)*e*dr2_dh(r2, m, p);

 return dmdh;
}

/*C*/
/************************************************************************/
void TOV_enthalpy(
           int    i_check,
               char   eos_type[],
               double e_center,
               double p_center,
               double enthalpy_center,
               double p_surface,
               double e_surface,
               double Gamma_P,
               double log_e_tab[2001],
               double log_p_tab[2001],
               double log_n0_tab[2001],
               double log_h_tab[2001],
               int    n_tab,
//               double r_is_gp[RDIV+1],
//               double lambda_gp[RDIV+1],
                  double r_is_gp[],
                  double lambda_gp[],
               char   eos_typeDM[1],
               double m_chi,
               double y_chi,
               double e_centerDM,
               double p_centerDM,
               double enthalpy_centerDM,
               double p_surfaceDM,
               double e_surfaceDM,
               double nu_gp[RDIV+1],
//           double enthalpy_gp[RDIV+1],
//           double enthalpy_gpDM[RDIV+1],
                  double enthalpy_gp[RDIV*10],
                  double enthalpy_gpDM[RDIV*10],
               double *r_is_final,
               double *r_final,
               double *m_final,
               double *rDM_is_final,
               double *rDM_final,
               double *mDM_final)
{
  int i=2,
      n_nearest,n_nearestDM;
  int b_index, d_index;

  double r2,                           /* radius */
         r_is,                        /* isotropic radial coordinate */
         r_is_est,                    /* estimate on final isotr. radius */
         r_is_check,                  /*                      */
         dr_is_save,                  /* r_is saving interval */
         rho0,
         e_d,                         /* density */
         p,                           /* pressure */
         hh,
         h,                           /* stepsize during integration */
         m,                           /* mass   */
         rho0DM,
         e_dDM,                         /* density */
         pDM,                           /* pressure */
         hhDM,
         mDM,                           /* mass   */
         nu_s,
         a1,a2,a3,a4,b1,b2,b3,b4,     /* coeff. in Runge-Kutta equations */
         c1,c2,c3,c4,
         a1DM,a2DM,a3DM,a4DM,b1DM,b2DM,b3DM,b4DM,     /* coeff. in Runge-Kutta equations */
         c1DM,c2DM,c3DM,c4DM,
         k_rescale,
//         r_gp[RDIV+1],
//         m_gp[RDIV+1],
//         e_d_gp[RDIV+1],
//         p_d_gp[RDIV+1],
//         m_gpDM[RDIV+1],
//         e_d_gpDM[RDIV+1],
//         p_d_gpDM[RDIV+1];
    r_gp[RDIV*10],
    m_gp[RDIV*10],
    e_d_gp[RDIV*10],
    p_d_gp[RDIV*10],
    m_gpDM[RDIV*10],
    e_d_gpDM[RDIV*10],
    p_d_gpDM[RDIV*10];

  double enthalpy_min;
  enthalpy_min = 1.0/(C*C);//1e-10;
//    enthalpy_min = 0.0;

  // SMM: Add enthalpy vectors
  // double enthalpy_gp[RDIV+1],
  //enthalpy_gpDM[RDIV+1];


//    if(i_check==1) {
//      if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"DM")==0))
//        r_is_est=1.5e6/sqrt(KAPPA);
//      h=r_is_est/(10000.0);
//    }
//    else {
//          r_is_est= max(*rDM_is_final,*r_is_final);
//          h=r_is_est/100000.0;
//            dr_is_save = max(*rDM_is_final,*r_is_final)/RDIV;
//          r_is_check = dr_is_save;
//     }
    
    h = -min(enthalpy_center, enthalpy_centerDM)/1.0e4;

    // Initialize variables at centre of star
    r_is=0.0;                            /* initial isotropic radius */
    r2=0.0;                               /* initial radius */
    m=0.0;                               /* initial mass */
    p=p_center;                          /* initial pressure */
    hh=enthalpy_center;
    mDM=0.0;
    pDM=p_centerDM;
    hhDM=enthalpy_centerDM;

    
    double RBM,RDM,RBMis,RDMis=0.0;
    n_nearest = n_tab/2;
//    r_is_gp[1]=0.0;
    r_gp[1]=0.0;
    m_gp[1]=0.0;
    m_gpDM[1]=0.0;
    lambda_gp[1]=0.0;
    e_d_gp[1] = e_at_p(p_center,p_surface, log_e_tab, log_p_tab, n_tab, &n_nearest, eos_type, Gamma_P);
    e_d_gpDM[1] = e_at_p_DM(pDM,p_surfaceDM, eos_typeDM, m_chi, y_chi);
    
    p_d_gp[1] = p;
    p_d_gpDM[1] = pDM;
    
    enthalpy_gp[1] = h_at_p(p_center,p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest);
    enthalpy_gpDM[1] = h_at_p_DM(p_centerDM,p_surfaceDM, eos_typeDM, m_chi, y_chi);

    printf("TOV_enthalpy: central enthalpy = %lf  dark central enthalpy = %lf \n", enthalpy_gp[1], enthalpy_gpDM[1]);
//    if(i_check==1)
//    {
//        printf("TOV: central enthalpy = %lf  dark central enthalpy = %lf \n", enthalpy_gp[1], enthalpy_gpDM[1]);
//    }
    

    while ( (hh>=enthalpy_min) || (hhDM>=enthalpy_min) ) {
 
      e_d = e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, &n_nearest, eos_type,
                                                                    Gamma_P);
      e_dDM = e_at_p_DM(pDM,p_surfaceDM, eos_typeDM, m_chi, y_chi);
//     if((i_check<=3) && (r_is>r_is_check) && (i<=RDIV)) {
      r_is_gp[i]=r_is;
      r_gp[i]=sqrt(r2);
      if((hh>=enthalpy_min)){
        m_gp[i]=m;
        e_d_gp[i]=e_d;
        p_d_gp[i]=p;
      }else{
        m_gp[i]=0.0;
        e_d_gp[i]=0.0;
        p_d_gp[i]=0.0;
    //enthalpy_gp[i] = 0.0;
      }
      if((pDM>=p_surfaceDM)){
        m_gpDM[i]=mDM;
        e_d_gpDM[i]=e_dDM;
        p_d_gpDM[i]=pDM;
      }else{
        m_gpDM[i]=0.0;
        e_d_gpDM[i]=0.0;
        p_d_gpDM[i]=0.0;
    //enthalpy_gpDM[i] = 0.0;
      }

      enthalpy_gp[i] = h_at_p(p,p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest);
      enthalpy_gpDM[i] = h_at_p_DM(pDM,p_surfaceDM, eos_typeDM, m_chi, y_chi);
      
//       if (i < RDIV)
//      printf("i=%d r_is=%lf enth_B=%lf enth_D=%lf \n", i, r_is_gp[i],enthalpy_gp[i],enthalpy_gpDM[i]);
      
      i++;
//      r_is_check += dr_is_save;
//     }

     
     

     if((hh>=enthalpy_min)){
       RBM=sqrt(r2);
       RBMis=r_is;
       (*r_is_final)=r_is;
       (*r_final)=sqrt(r2);
       b_index = i-1;
     }
     if((hhDM>=enthalpy_min)){
      RDM=sqrt(r2);
      RDMis=r_is;
      (*rDM_is_final)=r_is;
      (*rDM_final)=sqrt(r2);
      d_index = i-1;
     }
     (*m_final)=m;
     (*mDM_final)=mDM;

 
        if(r2 >= pow(0.4/(sqrt(KAPPA)*1e-5), 2.0))
        {
            
            
            a1=dr_is_dr(r_is,sqrt(r2),m+mDM)*dr2_dh(r2, m+mDM, p+pDM)/2.0/sqrt(r2);
            b1=dr2_dh(r2, m+mDM, p+pDM);
            
            if(hh>enthalpy_min){
                c1=dm_dh(r2, m, p, e_d);
            }
            
            if(hhDM>enthalpy_min){
                c1DM=dm_dh(r2, mDM, pDM, e_dDM);
            }
            
            
            a2=dr_is_dr(r_is+h*a1/2.0,sqrt(r2+h*b1/2.0),m+h*c1/2.0+mDM+h*c1DM/2.0)*dr2_dh(r2+h*b1/2.0, m+h*c1/2.0+mDM+h*c1DM/2.0, p_at_h(hh + h/2.0,pow(10.0,log_h_tab[1]), log_p_tab, log_h_tab, n_tab, &n_nearest) + p_at_h_DM(hhDM + h/2.0,enthalpy_min, eos_typeDM, m_chi, y_chi))/2.0/sqrt(r2+h*b1/2.0);
            b2=dr2_dh(r2+h*b1/2.0, m+h*c1/2.0+mDM+h*c1DM/2.0, p_at_h(hh + h/2.0,pow(10.0,log_h_tab[1]), log_p_tab, log_h_tab, n_tab, &n_nearest) + p_at_h_DM(hhDM + h/2.0,enthalpy_min, eos_typeDM, m_chi, y_chi));
            
            if(hh+h/2.0>enthalpy_min){
                c2=dm_dh(r2+h*b1/2.0, m+h*c1/2.0, p_at_h(hh + h/2.0,pow(10.0,log_h_tab[1]), log_p_tab, log_h_tab, n_tab, &n_nearest), e_at_h(hh + h/2.0,pow(10.0,log_h_tab[1]), log_e_tab, log_h_tab, n_tab, &n_nearest));
            }
            if(hhDM+h/2.0>enthalpy_min){
                c2DM=dm_dh(r2+h*b1/2.0, mDM+h*c1DM/2.0, p_at_h_DM(hhDM + h/2.0,enthalpy_min, eos_typeDM, m_chi, y_chi), e_at_h_DM(hhDM + h/2.0,enthalpy_min, eos_typeDM, m_chi, y_chi));
            }
            
            
            a3=dr_is_dr(r_is+h*a2/2.0,sqrt(r2+h*b2/2.0),m+h*c2/2.0+mDM+h*c2DM/2.0)*dr2_dh(r2+h*b2/2.0, m+h*c2/2.0+mDM+h*c2DM/2.0, p_at_h(hh + h/2.0,pow(10.0,log_h_tab[1]), log_p_tab, log_h_tab, n_tab, &n_nearest) + p_at_h_DM(hhDM + h/2.0,enthalpy_min, eos_typeDM, m_chi, y_chi))/2.0/sqrt(r2+h*b2/2.0);
            b3=dr2_dh(r2+h*b2/2.0, m+h*c2/2.0+mDM+h*c2DM/2.0, p_at_h(hh + h/2.0,pow(10.0,log_h_tab[1]), log_p_tab, log_h_tab, n_tab, &n_nearest) + p_at_h_DM(hhDM + h/2.0,enthalpy_min, eos_typeDM, m_chi, y_chi));
            
            if(hh+h/2.0>enthalpy_min){
                c3=dm_dh(r2+h*b2/2.0, m+h*c2/2.0, p_at_h(hh + h/2.0,pow(10.0,log_h_tab[1]), log_p_tab, log_h_tab, n_tab, &n_nearest), e_at_h(hh + h/2.0,pow(10.0,log_h_tab[1]), log_e_tab, log_h_tab, n_tab, &n_nearest));
            }
            if(hhDM+h/2.0>enthalpy_min){
                c3DM=dm_dh(r2+h*b2/2.0, mDM+h*c2DM/2.0, p_at_h_DM(hhDM + h/2.0,enthalpy_min, eos_typeDM, m_chi, y_chi), e_at_h_DM(hhDM + h/2.0,enthalpy_min, eos_typeDM, m_chi, y_chi));
            }
            
            a4=dr_is_dr(r_is+h*a3, sqrt(r2+h*b3), m + h*c3 + mDM + h*c3DM)*dr2_dh(r2 + h*b3, m + h*c3 + mDM + h*c3DM, p_at_h(hh + h,pow(10.0,log_h_tab[1]), log_p_tab, log_h_tab, n_tab, &n_nearest) + p_at_h_DM(hhDM + h,enthalpy_min, eos_typeDM, m_chi, y_chi))/2.0/sqrt(r2+h*b3);
            b4=dr2_dh(r2+h*b3, m+h*c3+mDM+h*c3DM, p_at_h(hh + h,pow(10.0,log_h_tab[1]), log_p_tab, log_h_tab, n_tab, &n_nearest) + p_at_h_DM(hhDM + h,enthalpy_min, eos_typeDM, m_chi, y_chi));
            
            if(hh+h>enthalpy_min){
                c4=dm_dh(r2+h*b3, m+h*c3, p_at_h(hh + h,pow(10.0,log_h_tab[1]), log_p_tab, log_h_tab, n_tab, &n_nearest), e_at_h(hh + h,pow(10.0,log_h_tab[1]), log_e_tab, log_h_tab, n_tab, &n_nearest));
            }
            if(hhDM+h>enthalpy_min){
                c4DM=dm_dh(r2+h*b3, mDM+h*c3DM, p_at_h_DM(hhDM + h,enthalpy_min, eos_typeDM, m_chi, y_chi), e_at_h_DM(hhDM + h,enthalpy_min, eos_typeDM, m_chi, y_chi));
            }
            
            
            
            r_is += (h/6.0)*(a1+2.*a2+2.*a3+a4);
            r2 += (h/6.0)*(b1+2.*b2+2.*b3+b4);
            if(hh+h>enthalpy_min){
                m += (h/6.0)*(c1+2.*c2+2.*c3+c4);
            }else{
                m += 0.0;
            }
            
            if(hhDM+h>enthalpy_min){
                mDM += (h/6.0)*(c1DM+2.*c2DM+2.*c3DM+c4DM);
            }else{
                mDM += 0.0;
            }
            
        } else
        {
            r2 = 3.0*(enthalpy_center + enthalpy_centerDM - hh - hhDM)/(2.0*PI*(3.0*(p_center + p_centerDM) + e_center + e_centerDM));
            r_is = sqrt(r2);
            m = 4.0*PI*e_center*pow(r2, 3.0/2.0)/3.0;
            mDM = 4.0*PI*e_centerDM*pow(r2, 3.0/2.0)/3.0;
        }

     
     hh += h;
     hhDM += h;
        
    p = p_at_h(hh,pow(10.0,log_h_tab[1]), log_p_tab, log_h_tab, n_tab, &n_nearest);
    pDM = p_at_h_DM(hhDM,enthalpy_min, eos_typeDM, m_chi, y_chi);
        
    printf("%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf \n", hh, hhDM, sqrt(r2)*sqrt(KAPPA)*1e-5, r_is*sqrt(KAPPA)*1e-5, m, mDM, (m + mDM)/(pow(r2, 3.0/2.0)), 1.0 - 2.0*(m + mDM)/sqrt(r2));

    }
   
//      printf("TOV_enthalpy: r_B_is = %lf \n", RBMis);
//      printf("TOV_enthalpy: r_B = %lf \n", RBM);
    
    // Interpolate to find the surface
//      enthalpy_min = 1.0/(C*C);
      n_nearest=b_index-1;
//      RBMis = interp(enthalpy_gp,r_is_gp,b_index,enthalpy_min,&n_nearest);
//      RBM = interp(enthalpy_gp,r_gp,b_index,enthalpy_min,&n_nearest);
////       printf("TOV: icheck=%d r_is_interp = %lf \n", i_check, RBMis);
////      printf("TOV: icheck=%d r_interp = %lf \n", i_check, RBM);
////
////       printf("TOV: icheck=%d r_B_isDM = %lf \n", i_check, RDMis);
//      n_nearest=d_index-1;
//    if(RDMis != 0.0)
//    {
//        RDMis = interp(enthalpy_gpDM,r_is_gp,d_index,enthalpy_min,&n_nearest);
//        RDM = interp(enthalpy_gpDM,r_gp,d_index,enthalpy_min,&n_nearest);
//    }
////      printf("TOV: icheck=%d r_is_interpDM = %lf \n", i_check, RDMis);
    
      
    
//    (*r_is_final)=RBMis;
//    (*r_final)=RBM;
//    (*rDM_is_final)=RDMis;
//    (*rDM_final)=RDM;
//    r_is_gp[RDIV]=max(RBMis,RDMis);
//    r_gp[RDIV]=max(RBM,RDM);
//    m_gp[RDIV]=(*m_final);
//    m_gpDM[RDIV]=(*mDM_final);


  


    
/* Rescale r_is and compute lambda */

//    if(i_check==3) {
      k_rescale=0.5*(max(RBM,RDM)/max(RBMis,RDMis))*(1.0-(*m_final+*mDM_final)/max(RBM,RDM)+
                sqrt(1.0-2.0*(*m_final+*mDM_final)/max(RBM,RDM)));
 
//      (*r_is_final) *= k_rescale;
//      (*rDM_is_final) *= k_rescale;
       
      nu_s = log((1.0-(*m_final+*mDM_final)/(2.0*max(RBMis,RDMis)*k_rescale))/(1.0+(*m_final+*mDM_final)/
                                        (2.0*max(RBMis,RDMis)*k_rescale)));
      for(i=1;i<=RDIV;i++) {
    r_is_gp[i] *= k_rescale;
 
    if(i==1){
      lambda_gp[1]= log(1.0/k_rescale);
         }
    else {lambda_gp[i]=log(r_gp[i]/r_is_gp[i]); }

    if(((e_d_gp[i]+e_d_gpDM[i])<(e_surfaceDM+e_surface))){
      hh=0.0;
    }
    else{
      p=p_at_e_DM(e_d_gpDM[i], e_surfaceDM, eos_typeDM, m_chi, y_chi);
      hh=h_at_p_DM(p, p_surfaceDM, eos_typeDM, m_chi, y_chi);
    }
    nu_gp[i]=nu_s-hh;
      }
      nu_gp[RDIV]=nu_s;
//    }
//    if(i_check==3) {
      printf(" \n" );
      printf("TOV_enthalpy: RK4 for non-rotating NS \n" );


      printf("R_iso = %lf ,       RD_iso = %lf\n",
         *r_is_final*sqrt(KAPPA)*1e-5,  *rDM_is_final*sqrt(KAPPA)*1e-5 );
      printf("R_sh  = %lf,        RD_sh  = %lf \n", *r_final*sqrt(KAPPA)*1e-5, *rDM_final*sqrt(KAPPA)*1e-5);
      printf("M     = %lf,        MD     = %lf \n",
         *m_final*sqrt(KAPPA)*C*C/(G*MSUN), *mDM_final*sqrt(KAPPA)*C*C/(G*MSUN));
      printf("M_tot = %lf \n",
         (*m_final + *mDM_final)*sqrt(KAPPA)*C*C/(G*MSUN));

      
      printf("GM/Rc^2 = %lf (At RD) \n", (*m_final+*mDM_final)/(*rDM_final));
      
      //printf("GM/R_isoc^2 = %lf (At RD) \n", (*m_final+*mDM_final)/(*rDM_is_final));
        // printf("radius at RD = %lf \n", *rDM_is_final*sqrt(KAPPA)*1e-5 *
        //pow( 1 + 0.5 * (*m_final+*mDM_final)/(*rDM_is_final), 2));


    
      printf(" \n" );
//    }
}

/*C*/
/************************************************************************/
void TOV_r(
           int    i_check,
               char   eos_type[],
               double e_center,
               double p_center,
               double p_surface,
               double e_surface,
               double Gamma_P,
               double log_e_tab[2001],
               double log_p_tab[2001],
               double log_n0_tab[2001],
               double log_h_tab[2001],
               int    n_tab,
               double r_is_gp[RDIV+1],
               double lambda_gp[RDIV+1],
               char   eos_typeDM[1],
               double m_chi,
               double y_chi,
               double e_centerDM,
               double p_centerDM,
               double p_surfaceDM,
               double e_surfaceDM,
               double nu_gp[RDIV+1],
           double enthalpy_gp[RDIV+1],
           double enthalpy_gpDM[RDIV+1],
               double *r_is_final,
               double *r_final,
               double *m_final,
               double *rDM_is_final,
               double *rDM_final,
               double *mDM_final)
{
    int i=2,
    n_nearest,n_nearestDM;
    int b_index, d_index;
    
    double r,                           /* radius */
    r_is,                        /* isotropic radial coordinate */
    r_est,                    /* estimate on final isotr. radius */
    r_check,                  /*                      */
    dr_save,                  /* r_is saving interval */
    rho0,
    e_d,                         /* density */
    p,                           /* pressure */
    h,                           /* stepsize during integration */
    m,                           /* mass   */
    rho0DM,
    e_dDM,                         /* density */
    pDM,                           /* pressure */
    mDM,                           /* mass   */
    nu_s,
    hh,
    a1,a2,a3,a4,b1,b2,b3,b4,     /* coeff. in Runge-Kutta equations */
    c1,c2,c3,c4,
    a1DM,a2DM,a3DM,a4DM,b1DM,b2DM,b3DM,b4DM,     /* coeff. in Runge-Kutta equations */
    c1DM,c2DM,c3DM,c4DM,
    k_rescale,
    r_gp[RDIV+1],
    m_gp[RDIV+1],
    e_d_gp[RDIV+1],
    p_d_gp[RDIV+1],
    m_gpDM[RDIV+1],
    e_d_gpDM[RDIV+1],
    p_d_gpDM[RDIV+1];
    
    double enthalpy_min;
    enthalpy_min = 1.0/(C*C);
    
    // SMM: Add enthalpy vectors
    // double enthalpy_gp[RDIV+1],
    //enthalpy_gpDM[RDIV+1];
    
    
    if(i_check==1) {
        if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"DM")==0))
            r_est=1.5e6/sqrt(KAPPA);
        h=r_est/10000.0;
    }
    else {
        r_est= max(*rDM_final,*r_final);
        printf("r_est = %.12e \n", r_est);
        h=r_est/100000.0;
        printf("h = %.12e \n", h);
        dr_save = max(*rDM_final,*r_final)/RDIV;
        r_check = dr_save;
    }
    
    // Initialize variables at centre of star
    r_is=0.0;                            /* initial isotropic radius */
    r=0.0;                               /* initial radius */
    m=0.0;                               /* initial mass */
    p=p_center;                          /* initial pressure */
    mDM=0.0;
    pDM=p_centerDM;
    
    
    double RBM,RDM,RBMis,RDMis=0.0;
    n_nearest = n_tab/2;
    r_is_gp[1]=0.0;
    r_gp[1]=0.0;
    m_gp[1]=0.0;
    m_gpDM[1]=0.0;
    lambda_gp[1]=0.0;
    e_d_gp[1] = e_at_p(p_center,p_surface, log_e_tab, log_p_tab, n_tab, &n_nearest, eos_type, Gamma_P);
    e_d_gpDM[1] = e_at_p_DM(p_centerDM,p_surfaceDM, eos_typeDM, m_chi, y_chi);
    p_d_gp[1] =p;
    p_d_gpDM[1] = pDM;
    enthalpy_gp[1] = h_at_p(p_center,p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest);
    enthalpy_gpDM[1] = h_at_p_DM(p_centerDM,p_surfaceDM, eos_typeDM, m_chi, y_chi);
    
    //    printf("TOV: icheck=%d central enthalpy = %lf  dark central enthalpy = %lf \n", i_check, enthalpy_gp[1], enthalpy_gpDM[1]);
    
    
    printf("TOV_r: P_B_c_cgs_found = %.12e  P_D_c_cgs_found = %.12e \n", p_center*C*C*C*C/(G*KAPPA), p_centerDM*C*C*C*C/(G*KAPPA));
    printf("TOV_r: Correct epsilon_B_c_cgs = %.12e  Correct epsilon_D_c_cgs = %.12e \n", e_center*C*C/(G*KAPPA), e_centerDM*C*C/(G*KAPPA));
    printf("TOV_r: epsilon_B_c_cgs found = %.12e  epsilon_D_c_cgs found = %.12e \n", e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, &n_nearest, eos_type, Gamma_P)*C*C/(G*KAPPA), e_at_p_DM(pDM,p_surfaceDM, eos_typeDM, m_chi, y_chi)*C*C/(G*KAPPA));
    printf("TOV_r: epsilon_B_c relative error = %.12e  epsilon_D_c relative error = %.12e \n", fabs(e_center - e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, &n_nearest, eos_type, Gamma_P))/e_center, fabs(e_centerDM - e_at_p_DM(pDM,p_surfaceDM, eos_typeDM, m_chi, y_chi))/e_centerDM);
    printf("TOV_r: h_B_c_cgs = %.12e  h_D_c_cgs = %.12e \n", enthalpy_gp[1]*C*C, enthalpy_gpDM[1]*C*C);
    printf("TOV_r: p_surf_cgs = %lf  p_surfDM = %lf \n", p_surface*C*C*C*C/(KAPPA*G), p_surfaceDM*C*C*C*C/(KAPPA*G));
    
    while ( (p>=p_surface) || (pDM>=p_surfaceDM) ) {
        
        e_d = e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, &n_nearest, eos_type,
                     Gamma_P);
        e_dDM = e_at_p_DM(pDM,p_surfaceDM, eos_typeDM, m_chi, y_chi);
        if((i_check<=3) && (r>r_check) && (i<=RDIV)) {
            r_is_gp[i]=r_is;
            r_gp[i]=r;
            if((p>=p_surface)){
                m_gp[i]=m;
                e_d_gp[i]=e_d;
                p_d_gp[i]=p;
            }else{
                m_gp[i]=0.0;
                e_d_gp[i]=0.0;
                p_d_gp[i]=0.0;
                //enthalpy_gp[i] = 0.0;
            }
            if((pDM>=p_surfaceDM)){
                m_gpDM[i]=mDM;
                e_d_gpDM[i]=e_dDM;
                p_d_gpDM[i]=pDM;
            }else{
                m_gpDM[i]=0.0;
                e_d_gpDM[i]=0.0;
                p_d_gpDM[i]=0.0;
                //enthalpy_gpDM[i] = 0.0;
            }
            
            enthalpy_gp[i] = h_at_p(p,p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest);
            enthalpy_gpDM[i] = h_at_p_DM(pDM,p_surfaceDM, eos_typeDM, m_chi, y_chi);
        
//            
//            i++;
//            r_check += dr_save;
        }
            
            
            
            
            if((p>=p_surface)){
                RBM=r;
                RBMis=r_is;
                (*r_is_final)=r_is;
                (*r_final)=r;
                b_index = i-1;
            }
            if((pDM>=p_surfaceDM)){
                RDM=r;
                RDMis=r_is;
                (*rDM_is_final)=r_is;
                (*rDM_final)=r;
                d_index = i-1;
            }
            (*m_final)=m;
            (*mDM_final)=mDM;
            
            
            a1=dr_is_dr(r_is,r,m+mDM);
            
            if(p>=p_surface){
                b1=dm_dr(r_is,r,m+mDM,p, e_center, p_surface, log_e_tab, log_p_tab, n_tab,
                         &n_nearest, eos_type, Gamma_P);
                c1=dp_dr(r_is,r,m+mDM,p,pDM, e_center, e_centerDM, p_surface, log_e_tab, log_p_tab, n_tab,
                         &n_nearest, eos_type, Gamma_P);
            }
            
            if(pDM>=p_surfaceDM){
                b1DM=dm_dr_DM(r_is,r,m+mDM,pDM, e_centerDM, p_surfaceDM, eos_typeDM, m_chi, y_chi);
                c1DM=dp_dr_DM(r_is,r,m+mDM,pDM,p, e_centerDM, e_center, p_surfaceDM, eos_typeDM, m_chi, y_chi);
            }
            
            
            a2=dr_is_dr(r_is+h*a1/2.0, r+h/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0);
            
            if(p+h*c1/2.0>=p_surface){
                b2=dm_dr(r_is+h*a1/2.0, r+h/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0, p+h*c1/2.0, e_center,
                         p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest,
                         eos_type, Gamma_P);
                
                c2=dp_dr(r_is+h*a1/2.0, r+h/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0, p+h*c1/2.0, pDM+h*c1DM/2.0, e_center, e_centerDM,
                         p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest,
                         eos_type, Gamma_P);
            } else
            {
                b2=m;
                c2=0.0;
            }
            if(pDM+h*c1DM/2.0>=p_surfaceDM){
                b2DM=dm_dr_DM(r_is+h*a1/2.0, r+h/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0, pDM+h*c1DM/2.0, e_centerDM,
                              p_surfaceDM, eos_typeDM, m_chi, y_chi);
                
                c2DM=dp_dr_DM(r_is+h*a1/2.0, r+h/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0,
                              pDM+h*c1DM/2.0, p+h*c1/2.0, e_centerDM, e_center, p_surfaceDM, eos_typeDM, m_chi, y_chi);
            } else
            {
                b2DM=mDM;
                c2DM=0.0;
            }
            
            
            a3=dr_is_dr(r_is+h*a2/2.0, r+h/2.0, m+h*b2/2.0+mDM+h*b2DM/2.0);
            
            if(p+h*c2/2.0>=p_surface){
                b3=dm_dr(r_is+h*a2/2.0, r+h/2.0, m+h*b2/2.0+mDM+h*b2DM/2.0, p+h*c2/2.0, e_center,
                         p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest,
                         eos_type, Gamma_P);
                
                c3=dp_dr(r_is+h*a2/2.0, r+h/2.0, m+h*b2/2.0+mDM+h*b2DM/2.0, p+h*c2/2.0, pDM+h*c2DM/2.0, e_center, e_centerDM,
                         p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest,
                         eos_type, Gamma_P);
            } else
            {
                b3=m;
                c3=0.0;
            }
            
            
            if(pDM+h*c2DM/2.0>=p_surfaceDM){
                b3DM=dm_dr_DM(r_is+h*a2/2.0, r+h/2.0, m+h*b2/2.0+mDM+h*b2DM/2.0, pDM+h*c2DM/2.0, e_centerDM,
                              p_surfaceDM, eos_typeDM, m_chi, y_chi);
                
                c3DM=dp_dr_DM(r_is+h*a2/2.0, r+h/2.0, m+h*b2/2.0+mDM+h*b2DM/2.0, pDM+h*c2DM/2.0, p+h*c2/2.0, e_centerDM, e_center,
                              p_surfaceDM, eos_typeDM, m_chi, y_chi);
            } else
            {
                b3DM=mDM;
                c3DM=0.0;
            }
            
            a4=dr_is_dr(r_is+h*a3, r+h, m+h*b3+mDM+h*b3DM);
            
            if(p+h*c3>=p_surface){
                b4=dm_dr(r_is+h*a3, r+h, m+h*b3+mDM+h*b3DM, p+h*c3, e_center, p_surface,
                         log_e_tab, log_p_tab, n_tab,&n_nearest,
                         eos_type, Gamma_P);
                
                c4=dp_dr(r_is+h*a3, r+h, m+h*b3+mDM+h*b3DM, p+h*c3,pDM+h*c3DM, e_center, e_centerDM, p_surface,
                         log_e_tab, log_p_tab, n_tab,&n_nearest,
                         eos_type, Gamma_P);
            } else
            {
                b4=m;
                c4=0.0;
            }
            if(pDM+h*c3DM>=p_surfaceDM){
                b4DM=dm_dr_DM(r_is+h*a3, r+h, m+h*b3+mDM+h*b3DM, pDM+h*c3DM, e_centerDM, p_surfaceDM,
                              eos_typeDM, m_chi, y_chi);
                
                c4DM=dp_dr_DM(r_is+h*a3, r+h, m+h*b3+mDM+h*b3DM, pDM+h*c3DM,p+h*c3, e_centerDM, e_center, p_surfaceDM,
                              eos_typeDM, m_chi, y_chi);
            } else
            {
                b4DM=mDM;
                c4DM=0.0;
            }
            
            
            
            r_is += (h/6.0)*(a1+2.*a2+2.*a3+a4);
            if(p+h*c1/2.0<p_surface || p+h*c2/2.0<p_surface || p+h*c3<p_surface || (p+(h/6.0)*(c1+2.*c2+2.*c3+c4))<p_surface){
                p=0.0;
                m=m;
            }else{
                p += (h/6.0)*(c1+2.*c2+2.*c3+c4);
                m += (h/6.0)*(b1+2.*b2+2.*b3+b4);
            }
            
            if(pDM+h*c1DM/2.0<p_surfaceDM || pDM+h*c2DM/2.0<p_surfaceDM || pDM+h*c3DM<p_surfaceDM || (pDM+(h/6.0)*(c1DM+2.*c2DM+2.*c3DM+c4DM))<p_surfaceDM){
                pDM=0.0;
                mDM=mDM;
            }else{
                pDM += (h/6.0)*(c1DM+2.*c2DM+2.*c3DM+c4DM);
                mDM += (h/6.0)*(b1DM+2.*b2DM+2.*b3DM+b4DM);
            }
            
            
            r += h;
            
        }
        
        
        
        // Interpolate to find the surface
        enthalpy_min = 1.0/(C*C);
        n_nearest=b_index-1;
        //      RBMis = interp(enthalpy_gp,r_is_gp,b_index,enthalpy_min,&n_nearest);
        //      RBM = interp(enthalpy_gp,r_gp,b_index,enthalpy_min,&n_nearest);
        //       printf("TOV: icheck=%d r_is_interp = %lf \n", i_check, RBMis);
        //      printf("TOV: icheck=%d r_interp = %lf \n", i_check, RBM);
        //
        //       printf("TOV: icheck=%d r_B_isDM = %lf \n", i_check, RDMis);
        n_nearest=d_index-1;
        
        //    printf("Before Interp: RDM = %lf \n", RDM*sqrt(KAPPA)*1e-5);
        //    if(RDMis != 0.0)
        //    {
        //        RDMis = interp(enthalpy_gpDM,r_is_gp,d_index,enthalpy_min,&n_nearest);
        //        RDM = interp(enthalpy_gpDM,r_gp,d_index,enthalpy_min,&n_nearest);
        //    }
        //    printf("After Interp: RDM = %lf \n", RDM*sqrt(KAPPA)*1e-5);
        //      printf("TOV: icheck=%d r_is_interpDM = %lf \n", i_check, RDMis);
        
        
        
        (*r_is_final)=RBMis;
        (*r_final)=RBM;
        (*rDM_is_final)=RDMis;
        (*rDM_final)=RDM;
        r_is_gp[RDIV]=max(RBMis,RDMis);
        r_gp[RDIV]=max(RBM,RDM);
        m_gp[RDIV]=(*m_final);
        m_gpDM[RDIV]=(*mDM_final);
        
        
        
        
        
        
        /* Rescale r_is and compute lambda */
        
        if(i_check==3) {
            k_rescale=0.5*(r_gp[RDIV]/r_is_gp[RDIV])*(1.0-(*m_final+*mDM_final)/r_gp[RDIV]+
                                                      sqrt(1.0-2.0*(*m_final+*mDM_final)/r_gp[RDIV]));
            
            (*r_is_final) *= k_rescale;
            (*rDM_is_final) *= k_rescale;
            
            nu_s = log((1.0-(*m_final+*mDM_final)/(2.0*r_is_gp[RDIV]*k_rescale))/(1.0+(*m_final+*mDM_final)/
                                                                                  (2.0*r_is_gp[RDIV]*k_rescale)));
            for(i=1;i<=RDIV;i++) {
                r_is_gp[i] *= k_rescale;
                
                if(i==1){
                    lambda_gp[1]= log(1.0/k_rescale);
                }
                else {lambda_gp[i]=log(r_gp[i]/r_is_gp[i]); }
                
                if(((e_d_gp[i]+e_d_gpDM[i])<(e_surfaceDM+e_surface))){
                    hh=0.0;
                }
                else{
                    p=p_at_e_DM(e_d_gpDM[i], e_surfaceDM, eos_typeDM, m_chi, y_chi);
                    hh=h_at_p_DM(p, p_surfaceDM, eos_typeDM, m_chi, y_chi);
                }
                nu_gp[i]=nu_s-hh;
            }
            nu_gp[RDIV]=nu_s;
        }
        //    if(i_check==3) {
        printf(" \n" );
        printf("TOV_r: RK4 for non-rotating NS \n" );
        
        
        printf("R_iso = %lf ,       RD_iso = %lf\n",
               *r_is_final*sqrt(KAPPA)*1e-5,  *rDM_is_final*sqrt(KAPPA)*1e-5 );
        printf("R_sh  = %lf,        RD_sh  = %lf \n", *r_final*sqrt(KAPPA)*1e-5, *rDM_final*sqrt(KAPPA)*1e-5);
        printf("M     = %lf,        MD     = %lf \n",
               *m_final*sqrt(KAPPA)*C*C/(G*MSUN), *mDM_final*sqrt(KAPPA)*C*C/(G*MSUN));
        printf("M_tot = %lf \n",
               (*m_final + *mDM_final)*sqrt(KAPPA)*C*C/(G*MSUN));
        
        
        printf("GM/Rc^2 = %lf (At RD) \n", (*m_final+*mDM_final)/(*rDM_final));
        
        //printf("GM/R_isoc^2 = %lf (At RD) \n", (*m_final+*mDM_final)/(*rDM_is_final));
        // printf("radius at RD = %lf \n", *rDM_is_final*sqrt(KAPPA)*1e-5 *
        //pow( 1 + 0.5 * (*m_final+*mDM_final)/(*rDM_is_final), 2));
        
        
        
        printf(" \n" );
        //    }
    }


/*C*/
/************************************************************************/
void TOV(
	       int    i_check, 
               char   eos_type[],
               double e_center,
               double p_center,
               double p_surface,
               double e_surface,
               double Gamma_P, 
               double log_e_tab[2001],
               double log_p_tab[2001],
               double log_n0_tab[2001],
               double log_h_tab[2001],
               int    n_tab,
               double r_is_gp[RDIV+1], 
               double lambda_gp[RDIV+1], 
               char   eos_typeDM[1],
               double m_chi,
               double y_chi,
               double e_centerDM,
               double p_centerDM,
               double p_surfaceDM,
               double e_surfaceDM,
               double nu_gp[RDIV+1],
	       double enthalpy_gp[RDIV+1],
	       double enthalpy_gpDM[RDIV+1],
               double *r_is_final, 
               double *r_final, 
               double *m_final,
               double *rDM_is_final, 
               double *rDM_final, 
               double *mDM_final,
         double *mDM_RB_final)
{
  int i=2,
      n_nearest,n_nearestDM;
  int b_index, d_index;

  double r,                           /* radius */
         r_is,                        /* isotropic radial coordinate */
         r_is_est,                    /* estimate on final isotr. radius */ 
         r_is_check,                  /*                      */    
         dr_is_save,                  /* r_is saving interval */  
         rho0,
         e_d,                         /* density */
         p,                           /* pressure */
         h,                           /* stepsize during integration */
         m,                           /* mass   */
         rho0DM,
         e_dDM,                         /* density */
         pDM,                           /* pressure */
         mDM,                           /* mass   */
         nu_s,
         hh,
         a1,a2,a3,a4,b1,b2,b3,b4,     /* coeff. in Runge-Kutta equations */
         c1,c2,c3,c4,
         a1DM,a2DM,a3DM,a4DM,b1DM,b2DM,b3DM,b4DM,     /* coeff. in Runge-Kutta equations */
         c1DM,c2DM,c3DM,c4DM,
         k_rescale, 
         r_gp[RDIV+1],
         m_gp[RDIV+1],
         e_d_gp[RDIV+1],
         p_d_gp[RDIV+1],
         m_gpDM[RDIV+1],
         e_d_gpDM[RDIV+1],
         p_d_gpDM[RDIV+1];
    
  double enthalpy_min;
  enthalpy_min = 1.0/(C*C);

  // SMM: Add enthalpy vectors
  // double enthalpy_gp[RDIV+1],
  //enthalpy_gpDM[RDIV+1];


    if(i_check==1) {
      if((strcmp(eos_type,"tab")==0) || (strcmp(eos_type,"DM")==0))
        r_is_est=1.5e6/sqrt(KAPPA);
      h=r_is_est/10000.0;
    }
    else {
          r_is_est= max(*rDM_is_final,*r_is_final);
          h=r_is_est/100000.0;
      	  dr_is_save = max(*rDM_is_final,*r_is_final)/RDIV;
    	  r_is_check = dr_is_save;
	 }

    // Initialize variables at centre of star
    r_is=0.0;                            /* initial isotropic radius */
    r=0.0;                               /* initial radius */
    m=0.0;                               /* initial mass */
    p=p_center;                          /* initial pressure */ 
    mDM=0.0;
    pDM=p_centerDM;

    
    double RBM,RDM,RBMis,RDMis=0.0;
    double m_RB_final = 0.0;
    *mDM_RB_final = 0.0;
    n_nearest = n_tab/2;
    r_is_gp[1]=0.0;
    r_gp[1]=0.0;
    m_gp[1]=0.0;
    m_gpDM[1]=0.0;
    lambda_gp[1]=0.0;
    e_d_gp[1] = e_at_p(p_center,p_surface, log_e_tab, log_p_tab, n_tab, &n_nearest, eos_type, Gamma_P);
    e_d_gpDM[1] = e_at_p_DM(pDM,p_surfaceDM, eos_typeDM, m_chi, y_chi);
    
  // test for files not existing.
//  if (out_file == NULL)
//    {
//      printf("Error! Could not open file\n");
//      exit(-1); // must include stdlib.h
//    }
    
    p_d_gp[1] =p;
    p_d_gpDM[1] = pDM;
    
    enthalpy_gp[1] = h_at_p(p_center,p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest);
    enthalpy_gpDM[1] = h_at_p_DM(pDM,p_surfaceDM, eos_typeDM, m_chi, y_chi);

//    printf("TOV: icheck=%d central enthalpy = %lf  dark central enthalpy = %lf \n", i_check, enthalpy_gp[1], enthalpy_gpDM[1]);
    if(i_check==1)
    {
        printf("TOV: central enthalpy = %lf  dark central enthalpy = %lf \n", enthalpy_gp[1], enthalpy_gpDM[1]);
        printf("TOV: central pressure = %.10e  dark central pressure = %.12e \n", p_center*C*C*C*C/(G*KAPPA), p_centerDM*C*C*C*C/(G*KAPPA));
    }
    printf("TOV: p_surf_cgs = %lf  p_surfDM = %lf \n", p_surface*C*C*C*C/(KAPPA*G), p_surfaceDM*C*C*C*C/(KAPPA*G));

    while ( (p>=p_surface) || (pDM>=p_surfaceDM) ) {
 
      e_d = e_at_p(p,p_surface, log_e_tab, log_p_tab, n_tab, &n_nearest, eos_type, 
                                                                    Gamma_P);
      e_dDM = e_at_p_DM(pDM,p_surfaceDM, eos_typeDM, m_chi, y_chi);
     if((i_check<=3) && (r_is>r_is_check) && (i<=RDIV)) {
      r_is_gp[i]=r_is;
      r_gp[i]=r;
      if((p>=p_surface)){
        m_gp[i]=m;
        e_d_gp[i]=e_d;
        p_d_gp[i]=p;
      }else{
        m_gp[i]=0.0;
        e_d_gp[i]=0.0;
        p_d_gp[i]=0.0;
	//enthalpy_gp[i] = 0.0;
      } 
      if((pDM>=p_surfaceDM)){
        m_gpDM[i]=mDM;
        e_d_gpDM[i]=e_dDM; 
        p_d_gpDM[i]=pDM;
      }else{
        m_gpDM[i]=0.0;
        e_d_gpDM[i]=0.0;       
        p_d_gpDM[i]=0.0;
	//enthalpy_gpDM[i] = 0.0;
      }

      enthalpy_gp[i] = h_at_p(p,p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest);
      enthalpy_gpDM[i] = h_at_p_DM(pDM,p_surfaceDM, eos_typeDM, m_chi, y_chi);
      
      i++;   
      r_is_check += dr_is_save;     
     }    

     
     

     if((p>=p_surface)){
       RBM=r;
       RBMis=r_is;
       (*r_is_final)=r_is;
       (*r_final)=r;
         m_RB_final = m;
         *mDM_RB_final = mDM;
       b_index = i-1;
     }
     if((pDM>=p_surfaceDM)){
      RDM=r;
      RDMis=r_is;
      (*rDM_is_final)=r_is;
      (*rDM_final)=r;
      d_index = i-1;
     }
     (*m_final)=m;
     (*mDM_final)=mDM;

 
     a1=dr_dr_is(r_is,r,m+mDM);
     
     if(p>=p_surface){
       b1=dm_dr_is(r_is,r,m+mDM,p, e_center, p_surface, log_e_tab, log_p_tab, n_tab,
                                                &n_nearest, eos_type, Gamma_P);
       c1=dp_dr_is(r_is,r,m+mDM,p,pDM, e_center, e_centerDM, p_surface, log_e_tab, log_p_tab, n_tab,
                                                &n_nearest, eos_type, Gamma_P);
     }

     if(pDM>=p_surfaceDM){
       b1DM=dm_dr_is_DM(r_is,r,m+mDM,pDM, e_centerDM, p_surfaceDM, eos_typeDM, m_chi, y_chi);
       c1DM=dp_dr_is_DM(r_is,r,m+mDM,pDM,p, e_centerDM, e_center, p_surfaceDM, eos_typeDM, m_chi, y_chi);
     }

     
     a2=dr_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0);

     if(p+h*c1/2.0>=p_surface){
       b2=dm_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0, p+h*c1/2.0, e_center,
                          p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest, 
                          eos_type, Gamma_P);

       c2=dp_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0, p+h*c1/2.0, pDM+h*c1DM/2.0, e_center, e_centerDM,
                          p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest,
                          eos_type, Gamma_P);
     } else
     {
         b2=m;
         c2=0.0;
     }
     if(pDM+h*c1DM/2.0>=p_surfaceDM){
       b2DM=dm_dr_is_DM(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0, pDM+h*c1DM/2.0, e_centerDM,
                          p_surfaceDM, eos_typeDM, m_chi, y_chi);

       c2DM=dp_dr_is_DM(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0+mDM+h*b1DM/2.0,
		     pDM+h*c1DM/2.0, p+h*c1/2.0, e_centerDM, e_center, p_surfaceDM, eos_typeDM, m_chi, y_chi);
     } else
     {
         b2DM=mDM;
         c2DM=0.0;
     }
                          

     a3=dr_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0+ mDM+h*b2DM/2.0);
     
     if(p+h*c2/2.0>=p_surface){
       b3=dm_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0+ mDM+h*b2DM/2.0, p+h*c2/2.0, e_center,
                          p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest, 
                          eos_type, Gamma_P);

       c3=dp_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0+ mDM+h*b2DM/2.0, p+h*c2/2.0, pDM+h*c2DM/2.0, e_center, e_centerDM,
                          p_surface, log_e_tab, log_p_tab, n_tab,&n_nearest,
                          eos_type, Gamma_P);
     } else
     {
         b3=m;
         c3=0.0;
     }
     
     
     if(pDM+h*c2DM/2.0>=p_surfaceDM){
       b3DM=dm_dr_is_DM(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0+ mDM+h*b2DM/2.0, pDM+h*c2DM/2.0, e_centerDM,
                          p_surfaceDM, eos_typeDM, m_chi, y_chi);

       c3DM=dp_dr_is_DM(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0+ mDM+h*b2DM/2.0, pDM+h*c2DM/2.0, p+h*c2/2.0, e_centerDM, e_center,
                          p_surfaceDM, eos_typeDM, m_chi, y_chi);
     } else
     {
         b3DM=mDM;
         c3DM=0.0;
     }

     a4=dr_dr_is(r_is+h, r+h*a3, m+h*b3+mDM+h*b3DM);

     if(p+h*c3>=p_surface){
       b4=dm_dr_is(r_is+h, r+h*a3, m+h*b3+mDM+h*b3DM, p+h*c3, e_center, p_surface,
                                log_e_tab, log_p_tab, n_tab,&n_nearest, 
                                eos_type, Gamma_P);

       c4=dp_dr_is(r_is+h, r+h*a3, m+h*b3+mDM+h*b3DM, p+h*c3,pDM+h*c3DM, e_center, e_centerDM, p_surface,
                                log_e_tab, log_p_tab, n_tab,&n_nearest,
                                eos_type, Gamma_P);
     } else
     {
         b4=m;
         c4=0.0;
     }
     if(pDM+h*c3DM>=p_surfaceDM){
       b4DM=dm_dr_is_DM(r_is+h, r+h*a3, m+h*b3+mDM+h*b3DM, pDM+h*c3DM, e_centerDM, p_surfaceDM,
                        eos_typeDM, m_chi, y_chi);

       c4DM=dp_dr_is_DM(r_is+h, r+h*a3, m+h*b3+mDM+h*b3DM, pDM+h*c3DM,p+h*c3, e_centerDM, e_center, p_surfaceDM,
                        eos_typeDM, m_chi, y_chi);
     } else
     {
         b4DM=mDM;
         c4DM=0.0;
     }

     
          
     r += (h/6.0)*(a1+2.*a2+2.*a3+a4);
     if(p+h*c1/2.0<p_surface || p+h*c2/2.0<p_surface || p+h*c3<p_surface || (p+(h/6.0)*(c1+2.*c2+2.*c3+c4))<p_surface){
         p=0.0;
           m=m;
       }else{
           p += (h/6.0)*(c1+2.*c2+2.*c3+c4);
           m += (h/6.0)*(b1+2.*b2+2.*b3+b4);
       }
        
        if(pDM+h*c1DM/2.0<p_surfaceDM || pDM+h*c2DM/2.0<p_surfaceDM || pDM+h*c3DM<p_surfaceDM || (pDM+(h/6.0)*(c1DM+2.*c2DM+2.*c3DM+c4DM))<p_surfaceDM){
            pDM=0.0;
              mDM=mDM;
        }else{
            pDM += (h/6.0)*(c1DM+2.*c2DM+2.*c3DM+c4DM);
            mDM += (h/6.0)*(b1DM+2.*b2DM+2.*b3DM+b4DM);
        }

     
     r_is += h;

    }
    
    // Interpolate to find the surface
//    enthalpy_min = 1.0/(C*C);
      n_nearest=b_index-1;
//      RBMis = interp(enthalpy_gp,r_is_gp,b_index,enthalpy_min,&n_nearest);
//      RBM = interp(enthalpy_gp,r_gp,b_index,enthalpy_min,&n_nearest);
//       printf("TOV: icheck=%d r_is_interp = %lf \n", i_check, RBMis);
//      printf("TOV: icheck=%d r_interp = %lf \n", i_check, RBM);
//
//       printf("TOV: icheck=%d r_B_isDM = %lf \n", i_check, RDMis);
      n_nearest=d_index-1;
    
//    printf("Before Interp: RDM = %lf \n", RDM*sqrt(KAPPA)*1e-5);
//    if(RDMis != 0.0)
//    {
//        RDMis = interp(enthalpy_gpDM,r_is_gp,d_index,enthalpy_min,&n_nearest);
//        RDM = interp(enthalpy_gpDM,r_gp,d_index,enthalpy_min,&n_nearest);
//    }
//    printf("After Interp: RDM = %lf \n", RDM*sqrt(KAPPA)*1e-5);
//      printf("TOV: icheck=%d r_is_interpDM = %lf \n", i_check, RDMis);
    
      
    
    (*r_is_final)=RBMis;
    (*r_final)=RBM;
    (*rDM_is_final)=RDMis;
    (*rDM_final)=RDM;
    r_is_gp[RDIV]=max(RBMis,RDMis);
    r_gp[RDIV]=max(RBM,RDM);
    m_gp[RDIV]=(*m_final);
    m_gpDM[RDIV]=(*mDM_final);


  


    
/* Rescale r_is and compute lambda */

    if(i_check==3) {
      k_rescale=0.5*(r_gp[RDIV]/r_is_gp[RDIV])*(1.0-(*m_final+*mDM_final)/r_gp[RDIV]+
                sqrt(1.0-2.0*(*m_final+*mDM_final)/r_gp[RDIV]));
 
      (*r_is_final) *= k_rescale;
      (*rDM_is_final) *= k_rescale;
       
      nu_s = log((1.0-(*m_final+*mDM_final)/(2.0*r_is_gp[RDIV]*k_rescale))/(1.0+(*m_final+*mDM_final)/
									    (2.0*r_is_gp[RDIV]*k_rescale)));      
      for(i=1;i<=RDIV;i++) {
	r_is_gp[i] *= k_rescale;
 
	if(i==1){ 
	  lambda_gp[1]= log(1.0/k_rescale);
         }
	else {lambda_gp[i]=log(r_gp[i]/r_is_gp[i]); }

	if(((e_d_gp[i]+e_d_gpDM[i])<(e_surfaceDM+e_surface))){
	  hh=0.0;          
	}
	else{ 
        p = p_at_e(e_d_gp[i], e_surface, log_p_tab, log_e_tab, n_tab, &n_nearest);
        pDM = p_at_e_DM(e_d_gpDM[i], e_surfaceDM, eos_typeDM, m_chi, y_chi);
        hh = h_at_p( p, p_surface, log_h_tab, log_p_tab, n_tab, &n_nearest) + h_at_p_DM(pDM, p_surfaceDM, eos_typeDM, m_chi, y_chi);
//        hh = h_at_p_DM(pDM, p_surfaceDM, eos_typeDM, m_chi, y_chi);
	}
	nu_gp[i]=nu_s-hh;
      }
      nu_gp[RDIV]=nu_s;
    }
//    if(i_check==3) {
      printf(" \n" );
      printf("TOV: RK4 for non-rotating NS \n" );


      printf("R_iso = %lf ,       RD_iso = %lf\n",
	     *r_is_final*sqrt(KAPPA)*1e-5,  *rDM_is_final*sqrt(KAPPA)*1e-5 );
      printf("R_sh  = %lf,        RD_sh  = %lf \n", *r_final*sqrt(KAPPA)*1e-5, *rDM_final*sqrt(KAPPA)*1e-5);
      printf("M     = %lf,        MD     = %lf \n",
	     *m_final*sqrt(KAPPA)*C*C/(G*MSUN), *mDM_final*sqrt(KAPPA)*C*C/(G*MSUN));
      printf("M_tot = %lf \n",
	     (*m_final + *mDM_final)*sqrt(KAPPA)*C*C/(G*MSUN));
    printf("M(RB)     = %lf,        MD(RB)     = %lf \n",
       m_RB_final*sqrt(KAPPA)*C*C/(G*MSUN), *mDM_RB_final*sqrt(KAPPA)*C*C/(G*MSUN));
    printf("M_tot(RB) = %lf \n",
       (m_RB_final + *mDM_RB_final)*sqrt(KAPPA)*C*C/(G*MSUN));

      
      printf("GM/Rc^2 = %lf (At RD) \n", (*m_final+*mDM_final)/max(*r_final,*rDM_final));
      
      //printf("GM/R_isoc^2 = %lf (At RD) \n", (*m_final+*mDM_final)/(*rDM_is_final));
	    // printf("radius at RD = %lf \n", *rDM_is_final*sqrt(KAPPA)*1e-5 *
	    //pow( 1 + 0.5 * (*m_final+*mDM_final)/(*rDM_is_final), 2));


    
      printf(" \n" );
//    }
}

/*C*/
/*************************************************************************/
void sphere(double s_gp[SDIV+1], 
	    double log_e_tab[2001], 
	    double log_p_tab[2001], 
	    double log_h_tab[2001],
	    double log_n0_tab[2001], 
	    int n_tab,                 
	    char eos_type[],
	    double Gamma_P, 
	    double e_center,
	    double p_center, 
	    double h_center,
	    double p_surface,
	    double e_surface,
	    char eos_typeDM[1],
        double m_chi,
        double y_chi,
	    double e_centerDM,
	    double p_centerDM, 
	    double h_centerDM,
	    double p_surfaceDM,
	    double e_surfaceDM,
	    double **rho,
	    double **gama,
	    double **alpha,
	    double **omega,
	    double **enthalpy,
	    double **enthalpyDM,
	    double *r_e,
	    double *r_eDM)

{
 int s,
     m,
     n_nearest;

 double r_is_s,
   r_is_final,
   r_final, 
   m_final,
   rDM_is_final,
   rDM_final, 
   mDM_final,
    mDM_RB_final,
   lambda_s,
   nu_s,
   enthalpy_s,
   enthalpy_sDM,
   r_is_gp[RDIV+1],
   lambda_gp[RDIV+1],
   nu_gp[RDIV+1],
   r_is_sp[SDIV+1],
   gama_mu_0[SDIV+1],
   rho_mu_0[SDIV+1],
   enthalpy_gp[RDIV+1],
   enthalpy_gpDM[RDIV+1],
   enthalpy_sp[SDIV+1],
   enthalpy_spDM[SDIV+1],
   gama_eq,
   rho_eq,
   R_is,
   s_e=0.5;
    
 FILE *output;
 

 // printf("Entered sphere: h_c = %lf  h_cDM = %lf \n", h_center, h_centerDM);

 
 /* The function TOV integrates the TOV equations. The function
	can be found in the file equil.c */
    
//    TOV_r(1, eos_type, e_center, p_center, p_surface, e_surface, Gamma_P,
//        log_e_tab, log_p_tab,log_n0_tab, log_h_tab, n_tab, r_is_gp, lambda_gp,
//        eos_typeDM, m_chi, y_chi, e_centerDM, p_centerDM, p_surfaceDM, e_surfaceDM,
//        nu_gp, enthalpy_gp, enthalpy_gpDM,
//        &r_is_final, &r_final, &m_final,&rDM_is_final, &rDM_final, &mDM_final);
//
//    TOV_r(2, eos_type, e_center, p_center, p_surface, e_surface, Gamma_P,
//        log_e_tab, log_p_tab,log_n0_tab, log_h_tab, n_tab, r_is_gp, lambda_gp,
//        eos_typeDM, m_chi, y_chi, e_centerDM, p_centerDM, p_surfaceDM, e_surfaceDM,
//        nu_gp,  enthalpy_gp, enthalpy_gpDM,
//        &r_is_final, &r_final, &m_final,&rDM_is_final, &rDM_final , &mDM_final);
//                 
//    TOV_r(3, eos_type, e_center, p_center, p_surface, e_surface, Gamma_P,
//        log_e_tab, log_p_tab,log_n0_tab, log_h_tab, n_tab, r_is_gp, lambda_gp,
//        eos_typeDM, m_chi, y_chi, e_centerDM, p_centerDM, p_surfaceDM, e_surfaceDM,
//        nu_gp, enthalpy_gp, enthalpy_gpDM,
//        &r_is_final, &r_final, &m_final,&rDM_is_final, &rDM_final , &mDM_final);

 TOV(1, eos_type, e_center, p_center, p_surface, e_surface, Gamma_P,
     log_e_tab, log_p_tab,log_n0_tab, log_h_tab, n_tab, r_is_gp, lambda_gp, 
     eos_typeDM, m_chi, y_chi, e_centerDM, p_centerDM, p_surfaceDM, e_surfaceDM,
     nu_gp, enthalpy_gp, enthalpy_gpDM,
     &r_is_final, &r_final, &m_final,&rDM_is_final, &rDM_final, &mDM_final, &mDM_RB_final);

 TOV(2, eos_type, e_center, p_center, p_surface, e_surface, Gamma_P,
     log_e_tab, log_p_tab,log_n0_tab, log_h_tab, n_tab, r_is_gp, lambda_gp, 
     eos_typeDM, m_chi, y_chi, e_centerDM, p_centerDM, p_surfaceDM, e_surfaceDM,
     nu_gp,  enthalpy_gp, enthalpy_gpDM,
     &r_is_final, &r_final, &m_final,&rDM_is_final, &rDM_final , &mDM_final, &mDM_RB_final);
              
 TOV(3, eos_type, e_center, p_center, p_surface, e_surface, Gamma_P,
     log_e_tab, log_p_tab,log_n0_tab, log_h_tab, n_tab, r_is_gp, lambda_gp, 
     eos_typeDM, m_chi, y_chi, e_centerDM, p_centerDM, p_surfaceDM, e_surfaceDM,
     nu_gp, enthalpy_gp, enthalpy_gpDM,
     &r_is_final, &r_final, &m_final,&rDM_is_final, &rDM_final , &mDM_final, &mDM_RB_final);
    
// TOV_enthalpy(1, eos_type, e_center, p_center, h_center, p_surface, e_surface, Gamma_P,
//        log_e_tab, log_p_tab,log_n0_tab, log_h_tab, n_tab, r_is_gp, lambda_gp,
//        eos_typeDM, m_chi, y_chi, e_centerDM, p_centerDM, h_centerDM, p_surfaceDM, e_surfaceDM,
//        nu_gp, enthalpy_gp, enthalpy_gpDM,
//        &r_is_final, &r_final, &m_final,&rDM_is_final, &rDM_final, &mDM_final);

 // printf("Sphere: TOV3: R_B = %lf  R_D = %lf \n", r_final*sqrt(KAPPA)*1e-5, rDM_final*sqrt(KAPPA)*1e-5);
 // printf("Sphere: TOV3: Riso_B = %lf  Riso_D = %lf \n", r_is_final*sqrt(KAPPA)*1e-5, rDM_is_final*sqrt(KAPPA)*1e-5);


 R_is=r_is_gp[RDIV];
 n_nearest=RDIV/2;
 double r_out,r_is_out;
 switch(Out_cond){
     case 1:
        r_out= ((r_final>=rDM_final)?r_final:rDM_final);
        r_is_out= ((r_is_final>=rDM_is_final)?r_is_final:rDM_is_final);
        break;
     case 2:
        r_out= ((r_final<=rDM_final)?r_final:rDM_final);
        r_is_out= ((r_is_final<=rDM_is_final)?r_is_final:rDM_is_final);
        break; 
     case 3:
        r_out= r_final;//((r_final>=rDM_final)?r_final:rDM_final);
        r_is_out= r_is_final;
        break;
     case 4:
        r_out= rDM_final;//((r_final<=rDM_final)?r_final:rDM_final);
        r_is_out= rDM_is_final;
        break;   
     case 5:
        r_out= (r_final+rDM_final)/2.;
        r_is_out= (r_is_final+rDM_is_final)/2.;
        break;    
   }   

 //output = fopen("enthalpy-sphere.txt","w");
 // fprintf(output,"#s r_is hB1 hB2 hD1 hD2 0.5*(gama+rho) error\n");

   
 for(s=1;s<=SDIV;s++) {
    r_is_s=r_is_out*(s_gp[s]/(1.0-s_gp[s]));
    r_is_sp[s] = r_is_s;

    if((r_is_s<R_is) ||(r_is_s<r_is_final)||(r_is_s<rDM_is_final)){
      lambda_s=interp(r_is_gp,lambda_gp,RDIV,r_is_s,&n_nearest);
      nu_s=interp(r_is_gp,nu_gp, RDIV,r_is_s,&n_nearest);
      enthalpy_s = interp(r_is_gp,enthalpy_gp,RDIV,r_is_s,&n_nearest);
      if (enthalpy_s < 0)
	enthalpy_s = 0.0;
    }
    else {
      lambda_s=2.0*log(1.0+(m_final+mDM_final)/(2.0*r_is_s));
      nu_s=log((1.0-(m_final+mDM_final)/(2.0*r_is_s))/(1.0+(m_final+mDM_final)/(2.0*r_is_s)));
      enthalpy_s = 0.0;
    }

    enthalpy_sDM = interp(r_is_gp,enthalpy_gpDM,RDIV,r_is_s,&n_nearest);
 
    gama[s][1]=nu_s+lambda_s;
    rho[s][1]=nu_s-lambda_s;

    //Compute enthalpy
    if (s==1){
      enthalpy[s][1] = h_center;
      enthalpyDM[s][1] = h_centerDM;
    }
    else{
      enthalpy[s][1] = h_center + 0.5*(gama[1][1]+rho[1][1] -gama[s][1] - rho[s][1]);
      enthalpyDM[s][1] = h_centerDM + 0.5*(gama[1][1]+rho[1][1] -gama[s][1] - rho[s][1]);


      enthalpy_sp[s] = enthalpy[s][1];
      //enthalpy_sp[s] = enthalpy_s;
      enthalpy_spDM[s] = enthalpyDM[s][1];
      
    }

    //   fprintf(output,"%d %lf %lf %lf %lf %lf %lf %6.5e\n", s, r_is_s,
    //	    enthalpy_s, enthalpy[s][1], enthalpy_sDM, enthalpyDM[s][1],
    //	    0.5*(gama[s][1]+rho[s][1]), 0.5*(gama[s][1]+rho[s][1]) - 0.5*(gama[1][1]+rho[1][1]) +  enthalpy[s][1] - h_center  );
    
    for(m=1;m<=MDIV;m++) {
        gama[s][m]=gama[s][1];        
        rho[s][m]=rho[s][1];
        alpha[s][m]=(gama[s][1]-rho[s][1])/2.0;
        omega[s][m]=0.0;

	enthalpy[s][m] = enthalpy[s][1];
	enthalpyDM[s][m] = enthalpyDM[s][1];	
    }
 
    gama_mu_0[s]=gama[s][1];                   /* gama at \mu=0 */
    rho_mu_0[s]=rho[s][1];                     /* rho at \mu=0 */

 }

   n_nearest=SDIV/2;
   int n_nearestDM=SDIV/2;

   s_e=r_is_final/(r_is_final+r_is_out);
   double s_eDM=rDM_is_final/(rDM_is_final+r_is_out);
   
   gama_eq = interp(s_gp,gama_mu_0,SDIV,s_e,&n_nearest); /* gama at equator */
   rho_eq = interp(s_gp,rho_mu_0,SDIV,s_e,&n_nearest);   /* rho at equator */
   double h_eq = interp(s_gp,enthalpy_sp,SDIV,s_e,&n_nearest); /* enthalpy at equator */

   // printf("Sphere: Version 1: Error in HSE = %6.5e    h_eq=%6.5e \n",
   //	  0.5*(gama_eq+rho_eq) - 0.5*(gama[1][1]+rho[1][1]) + h_eq - h_center, h_eq);

   
   double gama_eqDM = interp(s_gp,gama_mu_0,SDIV,s_eDM,&n_nearestDM); /* gama at equator */
   double rho_eqDM = interp(s_gp,rho_mu_0,SDIV,s_eDM,&n_nearestDM);   /* rho at equator */

     (*r_e)= r_final*exp(0.5*(rho_eq-gama_eq)); 
     (*r_eDM)= rDM_final*exp(0.5*(rho_eqDM-gama_eqDM));


        printf("Sphere: Isotropic rB_e = %lf \n", *r_e*sqrt(KAPPA)*1e-5);
        printf("Sphere: Isotropic rD_e = %lf \n", *r_eDM*sqrt(KAPPA)*1e-5);

     //printf("Sphere: Error in Isotropic Baryon Radius = %6.5e \n", (r_is_final - *r_e)/r_is_final);
     //printf("Sphere: Error in Isotropic Dark M Radius = %6.5e \n", (rDM_is_final - *r_eDM)/rDM_is_final);


     // Interpolation based on enthalpy
//     double r_test = interp(enthalpy_sp,r_is_sp,SDIV,0.0,&n_nearest);
//     printf("Sphere: r_isoB = %lf \n", r_test*sqrt(KAPPA)*1e-5);

     // Interpolation based on enthalpy
//     double r_testDM = interp(enthalpy_spDM,r_is_sp,SDIV,0.0,&n_nearest);
//     printf("Sphere: r_isoD = %lf \n", r_testDM*sqrt(KAPPA)*1e-5);

//     (*r_e) = r_test;
//     s_e=r_test/(r_test+r_is_out);
     gama_eq = interp(r_is_sp,gama_mu_0,SDIV,*r_e,&n_nearest); /* gama at equator */
     rho_eq = interp(r_is_sp,rho_mu_0,SDIV,*r_e,&n_nearest);   /* rho at equator */
     h_eq = interp(r_is_sp,enthalpy_sp,SDIV,*r_e,&n_nearest); /* enthalpy at equator */

       printf("Sphere: Error in HSE = %6.5e    h_eqB=%6.5e \n",
	  0.5*(gama_eq+rho_eq) - 0.5*(gama[1][1]+rho[1][1]) + h_eq - h_center, h_eq);
   
       //   fclose(output);
    
}


/*C*/
/*************************************************************************/
/* Main iteration cycle for computation of the rotating star's metric    */
/*************************************************************************/
void spin(double s_gp[SDIV+1],
	  double mu[MDIV+1],
	  double log_e_tab[2001], 
	  double log_p_tab[2001], 
	  double log_h_tab[2001],
	  double log_n0_tab[2001], 
	  int n_tab,                 
	  char eos_type[],
	  double Gamma_P, 
	  double h_center,
	  double enthalpy_min,
	  char eos_typeDM[1],
      double m_chi,
      double y_chi,
	  double h_centerDM,
	  double enthalpy_minDM,
	  double **rho,
	  double **gama,
	  double **alpha,
	  double **omega,
	  double **energy,
	  double **pressure,
	  double **enthalpy,
	  double **velocity_sq,
	  double **energyDM,
	  double **pressureDM,
	  double **enthalpyDM,
	  double **velocity_sqDM,
          double **Omega_hDM,
	  int    a_check,
	  double accuracy,
	  double cf,
	  double r_ratio,
	  double r_ratioDM,
	  double *r_e_new,
	  double *rDM_e_new,
	  double *Omega,
  	  double *OmegaDM)

 {
 int m,                      /* counter */
     s,                      /* counter */
     n,                      /* counter */
     k,                      /* counter */
     n_of_it=0,              /* number of iterations */
     n_nearest,
     n_nearestDM,
     print_dif = 0,
     i,
     j;

 double dark_h_RBP;
 
double **D2_rho,
  **D2_gama,
  **D2_omega;

float  ***f_rho,
       ***f_gama;

 
double   sum_rho=0.0,         /* intermediate sum in eqn for rho */
	 sum_gama=0.0,        /* intermediate sum in eqn for gama */
	 sum_omega=0.0,       /* intermediate sum in eqn for omega */
         r_e_old,             /* equatorial radius in previus cycle */
   	 dif=1.0,dif2=1.0,             /* difference | r_e_old - r_e | */
         d_gama_s,            /* derivative of gama w.r.t. s */
         d_gama_m,            /* derivative of gama w.r.t. m */
         d_rho_s,             /* derivative of rho w.r.t. s */
         d_rho_m,             /* derivative of rho w.r.t. m */
         d_omega_s,           /* derivative of omega w.r.t. s */
         d_omega_m,           /* derivative of omega w.r.t. m */
         d_gama_ss,           /* 2nd derivative of gama w.r.t. s */
         d_gama_mm,           /* 2nd derivative of gama w.r.t. m */
         d_gama_sm,           /* derivative of gama w.r.t. m and s */
         temp1,                /* temporary term in da_dm */ 
         temp2, 
         temp3,
         temp4,
         temp5,
         temp6,
         temp7,
         temp8,
         m1,                  
         s1,
         s2,
         ea,eaDM,
         rsm,
         gsm,
         omsm,
         esm,
         psm,
         esmDM,
         psmDM,
         v2sm,
         v2smDM,
    OmsmDM,
         mum,
         sgp,
         s_1,
         e_gsm,e_gsmDM,
         e_rsm,e_rsmDM, 
         rho0sm,
         term_in_Omega_h,
         term_in_Omega_hDM,
         r_p,
         s_p,
         gama_pole_h,                  /* gama^hat at pole */  
         gama_pole_hDM,                  /* gama^hat at pole */  
         gama_center_h,                /* gama^hat at center */
         gama_equator_h,               /* gama^hat at equator */
         gama_equator_hDM,               /* gama^hat at equator */
         rho_pole_h,                   /* rho^hat at pole */ 
         rho_pole_hDM,                   /* rho^hat at pole */ 
         rho_center_h,                 /* rho^hat at center */
         rho_equator_h,                /* rho^hat at equator */ 
         rho_equator_hDM,                /* rho^hat at equator */ 
         omega_equator_h,              /* omega^hat at equator */         
         omega_equator_hDM,              /* omega^hat at equator */         
         gama_mu_1[SDIV+1],            /* gama at \mu=1 */
         gama_mu_0[SDIV+1],            /* gama at \mu=0 */
         rho_mu_1[SDIV+1],             /* rho at \mu=1 */
         rho_mu_0[SDIV+1],             /* rho at \mu=0 */
         omega_mu_0[SDIV+1],           /* omega at \mu=0 */
    Omega_hDM_mu_0[SDIV+1],           /* omega at \mu=0 */
  enthalpy_0[SDIV+1],  /* enthalpy at mu=0 equator */
  enthalpy_1[SDIV+1],  /* enthalpy at mu=1 pole */
    enthalpyDM_0[SDIV+1],  /* enthalpy at mu=0 equator */
    enthalpyDM_1[SDIV+1],  /* enthalpy at mu=1 pole */
  h_eq, h_pole,
       **da_dm,
       **dgds,
       **dgdm,
       **D1_rho,
       **D1_gama,
       **D1_omega,
       **S_gama,
       **S_rho,
       **S_omega,
       **f2n,
       **P_2n,   
       **P1_2n_1,
         Omega_h,
//         Omega_hDM,
    Omega_h_eDM,
         sin_theta[MDIV+1],
         theta[MDIV+1],
         sk,
         sj,
         sk1,
         sj1,
         r_e,
         dum_e,
         k0,k1,k2,k3,
         kDM0,kDM1,kDM2,kDM3,
         l0,l1,l2,l3,
         lDM0,lDM1,lDM2,lDM3,r_eq;

         double s_e,
         s_eDM_new,
         s_pDM_new,
         s_e_new,
         r_eDM0,
         r_eDM, r_eDM_old,r_out,
         ratio_old, ratio_oldDM;

	 FILE *output;
	 char filename[80]="enthalpy.txt";


	 printf("\nSpin: \n");


    f2n = dmatrix(1,LMAX+1,1,SDIV);
    f_rho = f3tensor(1,SDIV,1,LMAX+1,1,SDIV);
    f_gama = f3tensor(1,SDIV,1,LMAX+1,1,SDIV);
 
    P_2n = dmatrix(1,MDIV,1,LMAX+1);   
    P1_2n_1 = dmatrix(1,MDIV,1,LMAX+1);


    // Create the Greens Functions (f_rho, f_gama)
    
    for(n=0;n<=LMAX;n++) 
       for(i=2;i<=SDIV;i++) f2n[n+1][i] = pow((1.0-s_gp[i])/s_gp[i],2.0*n);

    if(SMAX!=1.0) {

     for(j=2;j<=SDIV;j++)
        for(n=1;n<=LMAX;n++)
           for(k=2;k<=SDIV;k++) {
                 sk=s_gp[k];
                 sj=s_gp[j];
                 sk1=1.0-sk;
                 sj1=1.0-sj;

                 if(k<j) {   
                          f_rho[j][n+1][k] = f2n[n+1][j]*sj1/(sj*
                                  f2n[n+1][k]*sk1*sk1);
                          f_gama[j][n+1][k] = f2n[n+1][j]/(f2n[n+1][k]*sk*sk1);
	                 }else {     
                          f_rho[j][n+1][k] = f2n[n+1][k]/(f2n[n+1][j]*sk*sk1);
                          f_gama[j][n+1][k] = f2n[n+1][k]*sj1*sj1*sk/(sj*sj
                                            *f2n[n+1][j]*sk1*sk1*sk1);
                 }
	    }
     j=1;
 
       n=0; 
       for(k=2;k<=SDIV;k++) {
          sk=s_gp[k];
          f_rho[j][n+1][k]=1.0/(sk*(1.0-sk));
       }

       n=1;
       for(k=2;k<=SDIV;k++) {
          sk=s_gp[k];
          sk1=1.0-sk;         
          f_rho[j][n+1][k]=0.0;
          f_gama[j][n+1][k]=1.0/(sk*sk1);
       }

       for(n=2;n<=LMAX;n++)
          for(k=1;k<=SDIV;k++) {
             f_rho[j][n+1][k]=0.0;
             f_gama[j][n+1][k]=0.0;
          }


     k=1;

       n=0;
       for(j=1;j<=SDIV;j++)
          f_rho[j][n+1][k]=0.0;

       for(j=1;j<=SDIV;j++)
          for(n=1;n<=LMAX;n++) {
             f_rho[j][n+1][k]=0.0;
             f_gama[j][n+1][k]=0.0;
          }


     n=0;
     for(j=2;j<=SDIV;j++)
        for(k=2;k<=SDIV;k++) {
               sk=s_gp[k];
               sj=s_gp[j];
               sk1=1.0-sk;
               sj1=1.0-sj;

               if(k<j) 
                 f_rho[j][n+1][k] = sj1/(sj*sk1*sk1);
               else     
                 f_rho[j][n+1][k] = 1.0/(sk*sk1);
      }

   }
   else{      
        for(j=2;j<=SDIV-1;j++)
           for(n=1;n<=LMAX;n++)
              for(k=2;k<=SDIV-1;k++) {
                 sk=s_gp[k];
                 sj=s_gp[j];
                 sk1=1.0-sk;
                 sj1=1.0-sj;

                 if(k<j) {   
                          f_rho[j][n+1][k] = f2n[n+1][j]*sj1/(sj*
                                           f2n[n+1][k]*sk1*sk1);
                          f_gama[j][n+1][k] = f2n[n+1][j]/(f2n[n+1][k]*sk*sk1);
                 }else {     
                          f_rho[j][n+1][k] = f2n[n+1][k]/(f2n[n+1][j]*sk*sk1);

                          f_gama[j][n+1][k] = f2n[n+1][k]*sj1*sj1*sk/(sj*sj
                                            *f2n[n+1][j]*sk1*sk1*sk1);
                 }
	      }
   
        j=1;
 
          n=0; 
          for(k=2;k<=SDIV-1;k++) {
             sk=s_gp[k];
             f_rho[j][n+1][k]=1.0/(sk*(1.0-sk));
          }

          n=1;
          for(k=2;k<=SDIV-1;k++) {
             sk=s_gp[k];
             sk1=1.0-sk;         
             f_rho[j][n+1][k]=0.0;
             f_gama[j][n+1][k]=1.0/(sk*sk1);
          }

          for(n=2;n<=LMAX;n++)
             for(k=1;k<=SDIV-1;k++) {
                f_rho[j][n+1][k]=0.0;
                f_gama[j][n+1][k]=0.0;
             }

        k=1;
 
          n=0;
          for(j=1;j<=SDIV-1;j++)
             f_rho[j][n+1][k]=0.0;

          for(j=1;j<=SDIV-1;j++)
             for(n=1;n<=LMAX;n++) {
                f_rho[j][n+1][k]=0.0;
                f_gama[j][n+1][k]=0.0;
             }
 
 
        n=0;
          for(j=2;j<=SDIV-1;j++)
             for(k=2;k<=SDIV-1;k++) {
                sk=s_gp[k];
                sj=s_gp[j];
                sk1=1.0-sk;
                sj1=1.0-sj;

                if(k<j) 
                  f_rho[j][n+1][k] = sj1/(sj*sk1*sk1);
                else     
                  f_rho[j][n+1][k] = 1.0/(sk*sk1);
             }
 
        j=SDIV;
          for(n=1;n<=LMAX;n++)
             for(k=1;k<=SDIV;k++) {
                f_rho[j][n+1][k] = 0.0;
                f_gama[j][n+1][k] = 0.0;
             }

        k=SDIV;
          for(j=1;j<=SDIV;j++)
              for(n=1;n<=LMAX;n++) {
                 f_rho[j][n+1][k] = 0.0;
                 f_gama[j][n+1][k] = 0.0;
              }
   }

  n=0;
   for(i=1;i<=MDIV;i++)
      P_2n[i][n+1]=legendre(2*n,mu[i]);

   for(i=1;i<=MDIV;i++)
     for(n=1;n<=LMAX;n++) {
      P_2n[i][n+1]=legendre(2*n,mu[i]);
      P1_2n_1[i][n+1] = plgndr(2*n-1 ,1,mu[i]);
    }

  free_dmatrix(f2n,1,LMAX+1,1,SDIV);


  for(m=1;m<=MDIV;m++) { 
     sin_theta[m] = sqrt(1.0-mu[m]*mu[m]);  
     theta[m] = asin(sin_theta[m]);
  }

  
      r_eq =(*r_e_new);
      r_eDM=(*rDM_e_new);


   switch(Out_cond){
//    switch(4){
     case 1:
        r_out= (((*r_e_new)>=(*rDM_e_new))?(*r_e_new):(*rDM_e_new));
        break;
     case 2:
        r_out= (((*r_e_new)<=(*rDM_e_new))?(*r_e_new):(*rDM_e_new));
        break;
     case 3:
        r_out= (*r_e_new);//(((*r_e_new)>=(*rDM_e_new))?(*r_e_new):(*rDM_e_new));
        break;
     case 4:
        r_out= (*rDM_e_new);//(((*r_e_new)<=(*rDM_e_new))?(*r_e_new):(*rDM_e_new));
        break;   
     case 5:
        r_out= ((*r_e_new)+(*rDM_e_new))/2.;
        break;    
   }   

   // accuracy = 1e-4;
   
   // printf("%g %g %g\n",r_eq, r_eDM, r_out);
   //printf("spin: dif = %lf accuracy = %lf dif2=%lf \n", dif, accuracy, dif2);
    print_dif = 1;

    //    while(((dif> accuracy)||(dif2> accuracy)) || n_of_it<2) {
   while( (dif > accuracy) || n_of_it<2) {

      if(print_dif!=0)
	printf("dif=%4.3e  dif2=%4.3e r_out = %lf \n",dif,dif2, r_out*sqrt(KAPPA)*1e-5);
//       printf("r_eq = %f \n", r_eq*sqrt(KAPPA)*1e-5);

      /*sprintf(filename,"enthalpy_%d",n_of_it);
	 output = fopen(filename,"w");
	 printf(output,"#s hB hD 0.5*(rho+gama)\n");
      */

      
      /* Rescale potentials and construct arrays with the potentials along
       | the equatorial and polar directions.
      */        

      for(s=1;s<=SDIV;s++) {
		if (s <= 1)
        {
            printf("Spin C: s=%d enthalpy[s][1]=%f h_center=%f enthalpyDM[s][1]=%f h_centerDM=%f s=%d nu[s]=%lf\n",
                   s,enthalpy[s][1],h_center,enthalpyDM[s][1],h_centerDM,s, 0.5*(rho[s][1]+gama[s][1]));
//            printf("rho[s][1] = %f    gama[s][1] = %f \n", rho[s][1], gama[s][1]);
        }

//	fprintf(output,"%d %lf %lf %lf \n", s, enthalpy[s][1], enthalpyDM[s][1],0.5*(rho[s][1]+gama[s][1]));
	
	
         for(m=1;m<=MDIV;m++) {
            rho[s][m] /= SQ(r_out);
            gama[s][m] /= SQ(r_out);
            alpha[s][m] /= SQ(r_out);
            omega[s][m] *= r_out;
            Omega_hDM[s][m] = omega[s][m];
         }
         rho_mu_0[s]=rho[s][1];     
         gama_mu_0[s]=gama[s][1];   
         omega_mu_0[s]=omega[s][1]; 
         rho_mu_1[s]=rho[s][MDIV];  
         gama_mu_1[s]=gama[s][MDIV];
	 enthalpy_0[s] = enthalpy[s][1];
     enthalpy_1[s] = enthalpy[s][MDIV];
          enthalpyDM_0[s] = enthalpyDM[s][1];
          enthalpyDM_1[s] = enthalpyDM[s][MDIV];
          Omega_hDM_mu_0[s]=Omega_hDM[s][1];
      }
 
      /* Compute new r_e. */ 

      gama_center_h=gama[1][1]; 
      rho_center_h=rho[1][1]; 
     
      r_e_old=r_eq;
      s_p = r_ratio*r_eq/(r_ratio*r_eq+r_out);
       printf("r_ratio = %f \n", r_ratio);
       printf("s_p = %f \n", s_p);
       s_e=r_eq/(r_eq+r_out);
       printf("s_e = %f \n", s_e);
       
      n_nearest= SDIV/2;      
      gama_pole_h=interp(s_gp,gama_mu_1,SDIV,s_p,&n_nearest);
      rho_pole_h=interp(s_gp,rho_mu_1,SDIV,s_p,&n_nearest);
      h_pole=interp(s_gp,enthalpy_1,SDIV,s_p,&n_nearest);

      gama_equator_h=interp(s_gp,gama_mu_0,SDIV,s_e,&n_nearest);
      rho_equator_h=interp(s_gp,rho_mu_0,SDIV,s_e,&n_nearest); 
      h_eq=interp(s_gp,enthalpy_0,SDIV,s_e,&n_nearest);

      /* printf("Spin: Error in HSE equator = %6.5e    h_eq=%6.5e\n",
	     SQ(r_out)*(0.5*(gama_equator_h+rho_equator_h) - 0.5*(gama[1][1]+rho[1][1])) + h_eq - enthalpy[1][1],
	     h_eq);
      printf("Spin: Error in HSE pole    = %6.5e    h_p =%6.5e\n",
	     SQ(r_out)*(0.5*(gama_pole_h+rho_pole_h) - 0.5*(gama[1][1]+rho[1][1])) + h_pole - enthalpy[1][1],
	     h_pole);
      */
      
       r_eDM_old=r_eDM;
       n_nearestDM= SDIV/2;
       if(h_centerDM != 0.0)
       {
           //      s_pDM_new=r_ratioDM*r_eDM/(r_ratioDM*r_eDM+r_out); // This method assumes an r_ratio for the dark fluid
           
           // The following method DOES NOT assume any r_ratio for the dark fluid. The dark fluid differentially rotates
           // with zero angular momentum and angular velocity equal to the frame-dragging velocity of spacetime caused by the rigidly rotating
           // baryonic fluid.
           
           // This method is the fastest but has interpolation errors, resulting in unequal polar radii for dark fluid for spherical stars
//          s_pDM_new = interp(enthalpyDM_1, s_gp, SDIV, enthalpy_minDM, &n_nearest);
           
           // The following method is slower but more accurate
           for(s=1;s<=SDIV;s++)
           {
               if(enthalpyDM[s][MDIV] < enthalpy_minDM)
               {
                   s_pDM_new = s_gp[s];
                   break;
               }
           }
           
           
           s_eDM_new = r_eDM/(r_eDM + r_out);
       } else
       {
           s_pDM_new = 0.0;
           s_eDM_new = 0.0;
       }
       printf("s_pDM_new = %f \n", s_pDM_new);
       printf("s_eDM_new = %f \n", s_eDM_new);
//       printf("r_ratioDM = %f \n", r_);

      gama_pole_hDM=interp(s_gp,gama_mu_1,SDIV,s_pDM_new,&n_nearestDM);                  
      gama_equator_hDM=interp(s_gp,gama_mu_0,SDIV,s_eDM_new,&n_nearestDM);      


      rho_pole_hDM=interp(s_gp,rho_mu_1,SDIV,s_pDM_new,&n_nearestDM);   
      rho_equator_hDM=interp(s_gp,rho_mu_0,SDIV,s_eDM_new,&n_nearestDM);


      // Compute the value of the Dark enthalpy at the north pole of the baryonic surface
      dark_h_RBP = 0.5*( gama_pole_hDM+rho_pole_hDM - (gama_pole_h+rho_pole_h) )*SQ(r_out);
      // printf("Dark enthalpy at the NP of the Baryonic surface = %lf \n",dark_h_RBP);
       
       
       // Since the euatorial baryonic fluid radius is the parameter used for convergence, Equation A28 Cook, Shapiro, T
       // is used for calcuting the equatorial baryonic radius. The equatorial dark radius is set where dark enthalpy vanishes.
       
       // This method is the fastest but has interpolation errors, resulting in unequal polar radii for dark fluid for spherical stars
//      double s_eDM_new = interp(enthalpyDM_0, s_gp, SDIV, enthalpy_minDM, &n_nearest);
       
       // The following method is slower but more accurate
       if(h_centerDM != 0.0)
       {
           for(s=1;s<=SDIV;s++)
           {
               if(enthalpyDM[s][1] < enthalpy_minDM)
               {
                   s_eDM_new = s_gp[s];
                   break;
               }
           }
       } else
       {
           s_eDM_new = 0.0;
       }
       
       r_eDM = r_out*(s_eDM_new/(1.0 - s_eDM_new));
//       double r_pDM = r_out*(s_pDM_new/(1.0 - s_pDM_new));
       
      r_eq=sqrt(SQ(r_e_old)*2.0*(  h_center )/((gama_pole_h+rho_pole_h-gama_center_h-rho_center_h)*SQ(r_out))); // Equation A28 Cook, Shapiro, T
//       printf("r_e_old = %f km   h_center = %f  gama_pole_h = %f    rho_pole_h = %f gama_center_h = %f  rho_center_h = %f   r_out = %f km \n", r_e_old*sqrt(KAPPA)*1e-5, h_center, gama_pole_h, rho_pole_h, gama_center_h, rho_center_h, r_out*sqrt(KAPPA)*1e-5);
//       printf("r_eq = %f km \n", r_eq*sqrt(KAPPA)*1e-5);
      

   double r_out_old=r_out;
   
   switch(Out_cond){
     case 1:
        r_out= ((r_eq>=r_eDM)?r_eq:r_eDM);
        break;
     case 2:
        r_out= ((r_eq<=r_eDM)?r_eq:r_eDM);
        break;   
     case 3:
        r_out= r_eq;//((r_eq>=r_eDM)?r_eq:r_eDM);
        break;
     case 4:
        r_out= r_eDM;//((r_eq<=r_eDM)?r_eq:r_eDM);
        break;   
     case 5:
        r_out= (r_eq+r_eDM)/2.;
        break;    
   }       

   //printf("spin: D rout=%lf \n",  r_out*sqrt(KAPPA)*1e-5);

   
      s_e=r_eq/(r_eq+r_out);    
      s_eDM_new=r_eDM/(r_eDM+r_out);                        
    
      
      /* Compute angular velocity Omega. */
      
      
      if(r_ratio==1.0) {
          r_ratioDM=1.0;
        Omega_h=0.0;
        omega_equator_h=0.0;
      } 
      else {
            omega_equator_h=interp(s_gp,omega_mu_0,SDIV,s_e, &n_nearest);
            term_in_Omega_h=1.0-exp(SQ(r_out)*(gama_pole_h+rho_pole_h-gama_equator_h-rho_equator_h));
          printf("spin: gama_pole_h+rho_pole_h-gama_equator_h-rho_equator_h=%lf \n",  gama_pole_h+rho_pole_h-gama_equator_h-rho_equator_h);
          printf("spin: term_in_Omega_h=%lf \n",  term_in_Omega_h);
            if(term_in_Omega_h>=0.0)
               Omega_h = omega_equator_h + (1.0-s_e)/s_e*exp(SQ(r_out)*rho_equator_h)*sqrt(term_in_Omega_h);
            else {
                Omega_h=0.0;
	    }
      }
       printf("spin: Omega_h=%lf         omega_equator_h =%lf \n",  Omega_h, omega_equator_h);

       if(r_ratioDM==1.0) {
//        Omega_hDM=0.0;
           Omega_h_eDM=0.0;
        omega_equator_hDM=0.0;
      } 
      else {
            omega_equator_hDM=interp(s_gp,omega_mu_0,SDIV,s_eDM_new, &n_nearestDM);
//            term_in_Omega_hDM=1.0-exp(SQ(r_out)*(gama_pole_hDM+rho_pole_hDM-gama_equator_hDM-rho_equator_hDM));
//          term_in_Omega_hDM=1.0-exp(SQ(r_eDM)*(gama_pole_hDM+rho_pole_hDM-gama_equator_hDM-rho_equator_hDM));
          printf("spin: gama_pole_hDM+rho_pole_hDM-gama_equator_hDM-rho_equator_hDM=%lf \n",  gama_pole_hDM+rho_pole_hDM-gama_equator_hDM-rho_equator_hDM);
//          printf("spin: term_in_Omega_hDM=%lf \n",  term_in_Omega_hDM);
//            if(term_in_Omega_hDM>=0.0) {
//               Omega_hDM = omega_equator_hDM + (1.0-s_eDM_new)/s_eDM_new*exp(SQ(r_out)*rho_equator_hDM) *sqrt(term_in_Omega_hDM);
//            }
//            else {
//                Omega_hDM=0.0;
//            }
//          Omega_hDM = omega_equator_hDM + (1.0-s_eDM_new)/s_eDM_new*exp(SQ(r_out)*rho_equator_hDM) *sqrt(term_in_Omega_hDM);
//          Omega_hDM = omega_equator_hDM;
          Omega_h_eDM = interp(s_gp,Omega_hDM_mu_0,SDIV,s_eDM_new, &n_nearestDM);
      }
       printf("spin: Omega_h_eDM=%lf         omega_equator_hDM =%lf \n",  Omega_h_eDM, omega_equator_hDM);


      /* Compute velocity, energy density and pressure. */
      
      n_nearest=n_tab/2; 
      
      for(s=1;s<=SDIV;s++) {
         sgp=s_gp[s];

         for(m=1;m<=MDIV;m++) {
            rsm=rho[s][m];
            
            if((r_ratio==1.0)&&(r_ratioDM==1.0)) 
                velocity_sq[s][m]=0.0;
            else 
	      velocity_sq[s][m]=SQ((Omega_h-omega[s][m])*(sgp/(1.0-sgp))
                                  *sin_theta[m]*exp(-rsm*SQ(r_out)));

            if(velocity_sq[s][m]>=1.0) 
              velocity_sq[s][m]=0.0;

            if((r_ratio==1.0)&&(r_ratioDM==1.0)) 
                velocity_sqDM[s][m]=0.0;
            else               
//	      velocity_sqDM[s][m]=SQ((Omega_hDM-omega[s][m])*(sgp/(1.0-sgp))
//                                  *sin_theta[m]*exp(-rsm*SQ(r_out)));
             velocity_sqDM[s][m]=0.0;

            if(velocity_sqDM[s][m]>=1.0) 
              velocity_sqDM[s][m]=0.0;
              
                  
	    // Enthalpy definition based on values at NP
//	    enthalpy[s][m]=enthalpy_min + 0.5*(SQ(r_out)*(gama_pole_h+rho_pole_h
//                           -gama[s][m]-rsm)-log(1.0-velocity_sq[s][m]));

	    //	if (s <= 10 && m==1)
	    //printf("NP: enthalpy[s][1]=%f h_centre=%f eDM=%f hDM=%f s=%d\n",enthalpy[s][1],h_center,enthalpyDM[s][1],h_centerDM,s);
	    

	    // Enthalpy definition based on central value
	    enthalpy[s][m]=h_center + 0.5*(SQ(r_out)*(gama_center_h+rho_center_h
						      -gama[s][m]-rsm)-log(1.0-velocity_sq[s][m]));

	    //if (s <= 10 && m==1)
	    //printf("C : enthalpy[s][1]=%f h_centre=%f eDM=%f hDM=%f s=%d\n",enthalpy[s][1],h_center,enthalpyDM[s][1],h_centerDM,s);


            if((enthalpy[s][m]<enthalpy_min) || (sgp>s_e)) {
                  pressure[s][m]=0.0;
                  energy[s][m]=0.0; 
	    }
            else {

                     pressure[s][m]=p_at_h(enthalpy[s][m],pow(10.0,log_h_tab[1]), log_p_tab, 
                                           log_h_tab, n_tab, &n_nearest);
                     energy[s][m]=e_at_p(pressure[s][m],pow(10.0,log_p_tab[1]), log_e_tab, 
                                      log_p_tab, n_tab, &n_nearest, eos_type,
                                       Gamma_P);

	           }
	           
	    // enthalpy definition based on NP
//	    enthalpyDM[s][m]=enthalpy_minDM + 0.5*(SQ(r_out)*(gama_pole_hDM+rho_pole_hDM
//							      -gama[s][m]-rsm)-log(1.0-velocity_sqDM[s][m]));


	    // Enthalpy definition based on central value
	    enthalpyDM[s][m]=h_centerDM + 0.5*(SQ(r_out)*(gama_center_h+rho_center_h
						      -gama[s][m]-rsm)-log(1.0-velocity_sqDM[s][m]));
  
            if((enthalpyDM[s][m]<enthalpy_minDM) || (sgp>s_eDM_new)) {

                  pressureDM[s][m]=0.0;
                  energyDM[s][m]=0.0; 
	    }
	    
            else {
                     pressureDM[s][m]=p_at_h_DM(enthalpyDM[s][m],enthalpy_minDM, eos_typeDM, m_chi, y_chi);
                     energyDM[s][m]=e_at_p_DM(pressureDM[s][m],enthalpy_minDM, eos_typeDM, m_chi, y_chi);
	           }

           


            rho[s][m] *= SQ(r_out);
            gama[s][m] *= SQ(r_out);
            alpha[s][m] *= SQ(r_out);
	 }
      
      }
      // fclose(output);


         



      /* Compute metric potentials */
     // r_eq=sqrt(2*( h_center-enthalpy_min)/(gama_pole_h+rho_pole_h-gama_center_h-rho_center_h));
      S_gama = dmatrix(1,SDIV,1,MDIV);
      S_rho = dmatrix(1,SDIV,1,MDIV);
      S_omega = dmatrix(1,SDIV,1,MDIV);


      for(s=1;s<=SDIV;s++){
         for(m=1;m<=MDIV;m++) {
            rsm=rho[s][m];
            gsm=gama[s][m];
            omsm=omega[s][m];

            e_gsm=exp(0.5*gsm);
            e_rsm=exp(-rsm);            
            
            
            esm=energy[s][m];
            psm=pressure[s][m];
            v2sm=velocity_sq[s][m];            
            
            esmDM=energyDM[s][m];
            psmDM=pressureDM[s][m];
            v2smDM=velocity_sqDM[s][m];
            OmsmDM=Omega_hDM[s][m];
            
            


            mum=mu[m];            
            m1=1.0-SQ(mum);
            sgp=s_gp[s];
            s_1=1.0-sgp;
            s1=sgp*s_1;
            s2=SQ(sgp/s_1);  

            ea=16.0*PI*exp(2.0*alpha[s][m])*SQ(r_out);

            if(s==1) {
              d_gama_s=0.0;
              d_gama_m=0.0;
              d_rho_s=0.0;
              d_rho_m=0.0;
              d_omega_s=0.0;
              d_omega_m=0.0;
            }else{
                 d_gama_s=deriv_s(gama,s,m);
                 d_gama_m=deriv_m(gama,s,m);
                 d_rho_s=deriv_s(rho,s,m);
                 d_rho_m=deriv_m(rho,s,m);
                 d_omega_s=deriv_s(omega,s,m);
                 d_omega_m=deriv_m(omega,s,m);
	     }
             
//             if (s <= 1)
//             {
//                 printf("e_gsm = %f ea = %f esm = %f psm = %f s2 = %f v2sm = %f esmDM = %f psmDM = %f v2smDM = %f \n", e_gsm, ea, esm, psm, s2, v2sm, esmDM, psmDM, v2smDM);
//             }

            S_rho[s][m] = e_gsm*(0.5*ea*(esm + psm)*s2*(1.0+v2sm)/(1.0-v2sm) + 0.5*ea*(esmDM + psmDM)*s2*(1.0+v2smDM)/(1.0-v2smDM)
  
                          + s2*m1*SQ(e_rsm)*(SQ(s1*d_omega_s) 
                       
                          + m1*SQ(d_omega_m))
                         
                          + s1*d_gama_s - mum*d_gama_m + 0.5*rsm*(ea*(psm+psmDM)*s2  
 
                          - s1*d_gama_s*(0.5*s1*d_gama_s+1.0) 
 
                          - d_gama_m*(0.5*m1*d_gama_m-mum)));
             /*if((r_ratio==0.9)&&(S_rho[s][m]>100)) {
             
             printf("%g %g %g %g %g %g %g %d %d\n",S_rho[s][m],gsm,rsm,omsm,esm,psm,v2sm,s,m);
             }*/

            S_gama[s][m] = e_gsm*(ea*(psm+psmDM)*s2 + 0.5*gsm*(ea*(psm+psmDM)*s2 - 0.5*SQ(s1

                           *d_gama_s) - 0.5*m1*SQ(d_gama_m)));

            S_omega[s][m]=e_gsm*e_rsm*( -ea*(OmsmDM-omsm)*((esmDM+psmDM))

                          *s2/(1.0-v2smDM)  -ea*(Omega_h-omsm)*(esm+psm) *s2/(1.0-v2sm) 
                          
                          + omsm*( -0.5*ea*(((1.0+v2sm)*(esm) 
                           
                          + 2.0*v2sm*(psm))/(1.0-v2sm))*s2 
                          
                          -0.5*ea*(((1.0+v2smDM)*(esmDM) 
                           
                          + 2.0*v2smDM*(psmDM))/(1.0-v2smDM))*s2 

                          - s1*(2.0*d_rho_s+0.5*d_gama_s)

                          + mum*(2.0*d_rho_m+0.5*d_gama_m) + 0.25*SQ(s1)*(4.0

                          *SQ(d_rho_s)-SQ(d_gama_s)) + 0.25*m1*(4.0*SQ(d_rho_m)

                          - SQ(d_gama_m)) - m1*SQ(e_rsm)*(SQ(SQ(sgp)*d_omega_s)

                          + s2*m1*SQ(d_omega_m))));
          



	 }

     }

      /* ANGULAR INTEGRATION */
   
      D1_rho = dmatrix(1,LMAX+1,1,SDIV);
      D1_gama = dmatrix(1,LMAX+1,1,SDIV);
      D1_omega = dmatrix(1,LMAX+1,1,SDIV);

      n=0;
      for(k=1;k<=SDIV;k++) {      

         for(m=1;m<=MDIV-2;m+=2) {
               sum_rho += (DM/3.0)*(P_2n[m][n+1]*S_rho[k][m]
                          + 4.0*P_2n[m+1][n+1]*S_rho[k][m+1] 
                          + P_2n[m+2][n+1]*S_rho[k][m+2]);

	 }

         D1_rho[n+1][k]=sum_rho;
         D1_gama[n+1][k]=0.0;
         D1_omega[n+1][k]=0.0;
         sum_rho=0.0;

      }
  
      for(n=1;n<=LMAX;n++)
         for(k=1;k<=SDIV;k++) {      
            for(m=1;m<=MDIV-2;m+=2) {
                
//                if (m <= 1)
//                {
//                    printf("P_2n[m][n+1] = %f    S_rho[k][m] = %f \n", P_2n[m][n+1], S_rho[k][m]);
//                }

               sum_rho += (DM/3.0)*(P_2n[m][n+1]*S_rho[k][m]
                          + 4.0*P_2n[m+1][n+1]*S_rho[k][m+1] 
                          + P_2n[m+2][n+1]*S_rho[k][m+2]);
                       
               sum_gama += (DM/3.0)*(sin((2.0*n-1.0)*theta[m])*S_gama[k][m]
                           +4.0*sin((2.0*n-1.0)*theta[m+1])*S_gama[k][m+1]
                           +sin((2.0*n-1.0)*theta[m+2])*S_gama[k][m+2]);
  
               sum_omega += (DM/3.0)*(sin_theta[m]*P1_2n_1[m][n+1]*S_omega[k][m]
                            +4.0*sin_theta[m+1]*P1_2n_1[m+1][n+1]*S_omega[k][m+1]
                            +sin_theta[m+2]*P1_2n_1[m+2][n+1]*S_omega[k][m+2]);
                
//                if (m <= 1)
//                {
//                    printf("sum_rho = %f    sum_gama = %f \n", sum_rho, sum_gama);
//                }
	    }
            D1_rho[n+1][k]=sum_rho;
            D1_gama[n+1][k]=sum_gama;
            D1_omega[n+1][k]=sum_omega;
            sum_rho=0.0;
            sum_gama=0.0;
            sum_omega=0.0;
	}


      free_dmatrix(S_gama,1,SDIV,1,MDIV);
      free_dmatrix(S_rho,1,SDIV,1,MDIV);
      free_dmatrix(S_omega,1,SDIV,1,MDIV);



      /* RADIAL INTEGRATION */

      D2_rho = dmatrix(1,SDIV,1,LMAX+1);
      D2_gama = dmatrix(1,SDIV,1,LMAX+1);
      D2_omega = dmatrix(1,SDIV,1,LMAX+1);



      n=0;
      for(s=1;s<=SDIV;s++) {
            for(k=1;k<=SDIV-2;k+=2) { 
               sum_rho += (DS/3.0)*( f_rho[s][n+1][k]*D1_rho[n+1][k] 
                          + 4.0*f_rho[s][n+1][k+1]*D1_rho[n+1][k+1]
                          + f_rho[s][n+1][k+2]*D1_rho[n+1][k+2]);
                          
 	    }

	    D2_rho[s][n+1]=sum_rho;
	    D2_gama[s][n+1]=0.0;
	    D2_omega[s][n+1]=0.0;
            sum_rho=0.0;
	 }

 
      for(s=1;s<=SDIV;s++)
         for(n=1;n<=LMAX;n++) {
            for(k=1;k<=SDIV-2;k+=2) { 
               sum_rho += (DS/3.0)*( f_rho[s][n+1][k]*D1_rho[n+1][k] 
                          + 4.0*f_rho[s][n+1][k+1]*D1_rho[n+1][k+1]
                          + f_rho[s][n+1][k+2]*D1_rho[n+1][k+2]);

             

               sum_gama += (DS/3.0)*( f_gama[s][n+1][k]*D1_gama[n+1][k] 
                           + 4.0*f_gama[s][n+1][k+1]*D1_gama[n+1][k+1]
                           + f_gama[s][n+1][k+2]*D1_gama[n+1][k+2]);
     
               if(k<s && k+2<=s) 
                 sum_omega += (DS/3.0)*( f_rho[s][n+1][k]*D1_omega[n+1][k] 
                              + 4.0*f_rho[s][n+1][k+1]*D1_omega[n+1][k+1]
                              + f_rho[s][n+1][k+2]*D1_omega[n+1][k+2]);
               else {
                 if(k>=s) 
                   sum_omega += (DS/3.0)*( f_gama[s][n+1][k]*D1_omega[n+1][k] 
                                + 4.0*f_gama[s][n+1][k+1]*D1_omega[n+1][k+1]
                                + f_gama[s][n+1][k+2]*D1_omega[n+1][k+2]);
                 else
                   sum_omega += (DS/3.0)*( f_rho[s][n+1][k]*D1_omega[n+1][k] 
                                + 4.0*f_rho[s][n+1][k+1]*D1_omega[n+1][k+1]
                                + f_gama[s][n+1][k+2]*D1_omega[n+1][k+2]);
               }
	    }
	    D2_rho[s][n+1]=sum_rho;
	    D2_gama[s][n+1]=sum_gama;
	    D2_omega[s][n+1]=sum_omega;
            sum_rho=0.0;
            sum_gama=0.0;
            sum_omega=0.0;

	 }

      free_dmatrix(D1_rho,1,LMAX+1,1,SDIV);
      free_dmatrix(D1_gama,1,LMAX+1,1,SDIV);
      free_dmatrix(D1_omega,1,LMAX+1,1,SDIV);


      /* SUMMATION OF COEFFICIENTS */

      for(s=1;s<=SDIV;s++) 
         for(m=1;m<=MDIV;m++) {

            gsm=gama[s][m];
            rsm=rho[s][m];
            omsm=omega[s][m];             
            e_gsm=exp(-0.5*gsm);
            e_rsm=exp(rsm);
            temp1=sin_theta[m];
             
//             if (s <= 1)
//             {
//                 printf("e_gsm = %f    P_2n[m][0+1] = %f    D2_rho[s][0+1] = %f \n", e_gsm, P_2n[m][0+1], D2_rho[s][0+1]);
//             }

            sum_rho += -e_gsm*P_2n[m][0+1]*D2_rho[s][0+1]; 

            for(n=1;n<=LMAX;n++) {

               sum_rho += -e_gsm*P_2n[m][n+1]*D2_rho[s][n+1]; 

               if(m==MDIV) {             
                 sum_omega += 0.5*e_rsm*e_gsm*D2_omega[s][n+1]; 
                 sum_gama += -(2.0/PI)*e_gsm*D2_gama[s][n+1];   
	       }
               else { 
                     sum_omega += -e_rsm*e_gsm*(P1_2n_1[m][n+1]/(2.0*n
                                  *(2.0*n-1.0)*temp1))*D2_omega[s][n+1];
  
                     sum_gama += -(2.0/PI)*e_gsm*(sin((2.0*n-1.0)*theta[m])
                                 /((2.0*n-1.0)*temp1))*D2_gama[s][n+1];   
	       }
	    }

            rho[s][m]=rsm + cf*(sum_rho-rsm);
            gama[s][m]=gsm + cf*(sum_gama-gsm);
            omega[s][m]=omsm + cf*(sum_omega-omsm);

            sum_omega=0.0;
            sum_rho=0.0;
            sum_gama=0.0; 
	  }

      free_dmatrix(D2_rho,1,SDIV,1,LMAX+1);
      free_dmatrix(D2_gama,1,SDIV,1,LMAX+1);
      free_dmatrix(D2_omega,1,SDIV,1,LMAX+1);


      /* CHECK FOR DIVERGENCE */

      if(fabs(omega[2][1])>100.0 || fabs(rho[2][1])>100.0 
         || fabs(gama[2][1])>300.0) {
         a_check=200; 
         break;
      }


      /* TREAT SPHERICAL CASE */
      
      if((r_ratio==1.0)&&(r_ratioDM==1.0)) {
        for(s=1;s<=SDIV;s++)
           for(m=1;m<=MDIV;m++) {
              rho[s][m]=rho[s][1];
              gama[s][m]=gama[s][1];
              omega[s][m]=0.0;          
	   }
      }


      /* TREAT INFINITY WHEN SMAX=1.0 */

      if(SMAX==1.0) {
         for(m=1;m<=MDIV;m++) {
            rho[SDIV][m]=0.0;
            gama[SDIV][m]=0.0;
            omega[SDIV][m]=0.0;
	 }
      } 


      /* COMPUTE FIRST ORDER DERIVATIVES OF GAMA */ 


      da_dm = dmatrix(1,SDIV,1,MDIV);
      dgds = dmatrix(1,SDIV,1,MDIV);
      dgdm = dmatrix(1,SDIV,1,MDIV); 
 
      for(s=1;s<=SDIV;s++)
         for(m=1;m<=MDIV;m++) {
            dgds[s][m]=deriv_s(gama,s,m);
            dgdm[s][m]=deriv_m(gama,s,m);
	 }



      /* ALPHA */
 
      if((r_ratio==1.0)&&(r_ratioDM==1.0)) {
        for(s=1;s<=SDIV;s++)
           for(m=1;m<=MDIV;m++)
              da_dm[s][m]=0.0; 
      } 
      else {
            for(s=2;s<=SDIV;s++)
               for(m=1;m<=MDIV;m++) {

                  da_dm[1][m]=0.0; 
       
                  sgp=s_gp[s];
                  s1=sgp*(1.0-sgp);
                  mum=mu[m]; 
                  m1=1.0-SQ(mum);
          
                  d_gama_s=dgds[s][m];
                  d_gama_m=dgdm[s][m];
                  d_rho_s=deriv_s(rho,s,m);
                  d_rho_m=deriv_m(rho,s,m);
                  d_omega_s=deriv_s(omega,s,m);
                  d_omega_m=deriv_m(omega,s,m);
                  d_gama_ss=s1*deriv_s(dgds,s,m)+(1.0-2.0*sgp)
                                               *d_gama_s;
                  d_gama_mm=m1*deriv_m(dgdm,s,m)-2.0*mum*d_gama_m;  
                  d_gama_sm=deriv_sm(gama,s,m);

           temp1=2.0*SQ(sgp)*(sgp/(1.0-sgp))*m1*d_omega_s*d_omega_m

                *(1.0+s1*d_gama_s) - (SQ(SQ(sgp)*d_omega_s) - 
 
                SQ(sgp*d_omega_m/(1.0-sgp))*m1)*(-mum+m1*d_gama_m); 
  
           temp2=1.0/(m1 *SQ(1.0+s1*d_gama_s) + SQ(-mum+m1*d_gama_m));

           temp3=s1*d_gama_ss + SQ(s1*d_gama_s);
  
           temp4=d_gama_m*(-mum+m1*d_gama_m);
   
           temp5=(SQ(s1*(d_rho_s+d_gama_s)) - m1*SQ(d_rho_m+d_gama_m))

                 *(-mum+m1*d_gama_m);

           temp6=s1*m1*(0.5*(d_rho_s+d_gama_s)* (d_rho_m+d_gama_m) 
  
                + d_gama_sm + d_gama_s*d_gama_m)*(1.0 + s1*d_gama_s); 

           temp7=s1*mum*d_gama_s*(1.0+s1*d_gama_s);

           temp8=m1*exp(-2.0*rho[s][m]);
 
          da_dm[s][m] = -0.5*(d_rho_m+d_gama_m) - temp2*(0.5*(temp3 - 

            d_gama_mm - temp4)*(-mum+m1*d_gama_m) + 0.25*temp5 

            - temp6 +temp7 + 0.25*temp8*temp1);	 
       }
   }

      for(s=1;s<=SDIV;s++) {
         alpha[s][1]=0.0;
         for(m=1;m<=MDIV-1;m++) 
            alpha[s][m+1]=alpha[s][m]+0.5*DM*(da_dm[s][m+1]+
                          da_dm[s][m]);
      } 
 

   free_dmatrix(da_dm,1,SDIV,1,MDIV);
   free_dmatrix(dgds,1,SDIV,1,MDIV);
   free_dmatrix(dgdm,1,SDIV,1,MDIV);


      for(s=1;s<=SDIV;s++) {        
         for(m=1;m<=MDIV;m++) {   

            alpha[s][m] += -alpha[s][MDIV]+0.5*(gama[s][MDIV]-rho[s][MDIV]);


            if(alpha[s][m]>=300.0) {
              a_check=200; 
              break;
            }
            omega[s][m] /= r_out;
         } 


              	     }
      if(SMAX==1.0) {
         for(m=1;m<=MDIV;m++)      
            alpha[SDIV][m] = 0.0;
      }

      if(a_check==200)
        break;



        dif=fabs(r_e_old-r_eq)/r_eq;
       
       if(h_centerDM != 0.0)
       {
           dif2=fabs(r_eDM_old-r_eDM)/r_eDM;
       }



      n_of_it++;

 }   /* end while */
//printf("N of it= %d\n", n_of_it);
//printf("%g %g %g %g %g\n",velocity_sq[1][1],alpha[1][1], gama[1][1], rho[1][1], omega[1][1]);

     /* COMPUTE OMEGA */  
 
    /* UPDATE r_e_new */

    (*r_e_new) = r_eq;
    (*rDM_e_new) = r_eDM;
    (*Omega) = Omega_h*C/(r_out*sqrt(KAPPA));
    (*OmegaDM) = Omega_h_eDM*C/(r_out*sqrt(KAPPA));
    
    printf("Omega/(2pi) = %lf\n", *Omega/(2.0*PI));
    printf("OmegaDM/(2pi) = %lf\n", *OmegaDM/(2.0*PI));


//printf("%g %g %g %g\n",log(Omega_h), log(Omega_hDM), log(r_ratio),  log(r_ratioDM));
    free_f3tensor(f_rho, 1,SDIV,1,LMAX+1,1,SDIV);
    free_f3tensor(f_gama,1,SDIV,1,LMAX+1,1,SDIV);
    free_dmatrix(P_2n,   1,MDIV,1,LMAX+1);   
    free_dmatrix(P1_2n_1,1,MDIV,1,LMAX+1);  
}
