double f(double ptot,double etot);

double h_integ_tab(double p1[RDIV+1],double e1[RDIV+1], double p2[RDIV+1],double e2[RDIV+1],int Num);

double max(double v1,double v2);

double min(double v1,double v2);

double dN_dr_is(double r_is, double r,double m, double p, double pDM, double p_surface, double p_surfaceDM,double e_center);

void make_grid( double s_gp[SDIV+1], 
                double mu[MDIV+1]);                        

void load_eos( char eos_file[], 
               double log_e_tab[2001], 
               double log_p_tab[2001], 
               double log_h_tab[2001],
               double log_n0_tab[2001], 
               int *n_tab);

double epsilon_D_EOS(char eos_typeDM[1],
                     double m_chi,
                     double y_chi,
                     double x_D);

double P_D_EOS(char eos_typeDM[1],
               double m_chi,
               double y_chi,
               double x_D);

double h_D_EOS(char eos_typeDM[1],
               double m_chi,
               double y_chi,
               double x_D);

double e_of_rho0(double rho0, double Gamma_P);

double e_at_p(double pp, 
              double pp_surface,
              double log_e_tab[2001], 
              double log_p_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt,
              char eos_type[],
              double Gamma_P);

double e_at_p_DM(double pp,
                 double pp_surface,
                 char eos_typeDM[1],
                 double m_chi,
                 double y_chi);

double e_at_h(double hh,
              double hh_surface,
              double log_e_tab[2001],
              double log_h_tab[2001],
              int    n_tab,
              int    *n_nearest_pt);

double e_at_h_DM(double hh,
                 double hh_surface,
                 char eos_typeDM[1],
                 double m_chi,
                 double y_chi);

double p_at_e(double ee,
              double ee_surface,
              double log_p_tab[2001], 
              double log_e_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt);

double p_at_e_DM(double ee,
                 double ee_surface,
                 char eos_typeDM[1],
                 double m_chi,
                 double y_chi);

double p_at_h(double hh, 
              double hh_surface,
              double log_p_tab[2001], 
              double log_h_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt);

double p_at_h_DM(double hh,
                 double hh_surface,
                 char eos_typeDM[1],
                 double m_chi,
                 double y_chi);

double h_at_p(double pp, 
              double pp_surface,
              double log_h_tab[2001], 
              double log_p_tab[2001],
              int    n_tab, 
              int    *n_nearest_pt);

double h_at_p_DM(double pp,
                 double pp_surface,
                 char eos_typeDM[1],
                 double m_chi,
                 double y_chi);

double n0_at_e(double ee, 
               double ee_surface,
               double log_n0_tab[2001], 
               double log_e_tab[2001],
               int    n_tab, 
               int    *n_nearest_pt);


double dpds(double p,
            double p_surface,
	    double e,
	    double rho_s,
	    double gamma_s,
            double s,
            double v,
            double Omega,
            double omega,
            double omega_s);

double dpdm(double p,
            double p_surface,
	    double e,
	    double rho_m,
	    double gamma_m,
            double m,
            double v,
            double Omega,
            double omega,
            double omega_m);

double y_at_x(double x,
	      double x1,
	      double x2,
	      double y1,
	      double y2);





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
	       double p_surface);





void make_centerDM(
           char eos_typeDM[1],
           double m_chi,
           double y_chi,
           double e_centerDM,
           double *p_centerDM,
           double *h_centerDM,
           double e_surfaceDM,
           double p_surfaceDM);




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
		 double *Mp); //49

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
                char   eos_type[],
                double Gamma_P);

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
                char   eos_type[],
                double Gamma_P);

double dm_dr_is_DM(double r_is,
                   double r,
                   double m,
                   double p,
                   double e_center,
                   double p_surface,
                   char eos_typeDM[1],
                   double m_chi,
                   double y_chi);

double dm_dr_DM(double r_is,
                   double r,
                   double m,
                   double p,
                   double e_center,
                   double p_surface,
                   char eos_typeDM[1],
                   double m_chi,
                   double y_chi);

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
                double Gamma_P);

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
                double Gamma_P);

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
                   double y_chi);

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
                   double y_chi);

double dr_dr_is(double r_is, double r, double m);

double dr_is_dr(double r_is, double r, double m);

double dr2_dh(double r2, double m, double p);

double dm_dh(double r2, double m, double p, double e);

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
               double r_is_gp[RDIV+1],
               double lambda_gp[RDIV+1],
               char   eos_typeDM[1],
               double m_chi,
               double y_chi,
               double e_centerDM,
               double p_centerDM,
               double enthalpy_centerDM,
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
                  double *mDM_final);

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
               double *mDM_final);

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
         double *mDM_RB_final);

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
	    double *r_eDM);

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
	  double *r_eDM_new,
	  double *Omega,
	  double *OmegaDM) ;
