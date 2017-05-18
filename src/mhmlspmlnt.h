/*
 *  mhmlspmlnt.h
 *  MuLoNeTests
 *
 *  Created by sonsins on Wed Mar 12 2003.
 *
 */


struct var /*structure to introduce in the coalescent and ML programs. All data and results ml. ONLY ONE VECTOR.*/
{
    int n_loci;

	/*MC coalescent simulations*/
    long unsigned n_iter;
    long int seed1;
    double tlimit;
	
    double time_spec;
	int no_rec_males;
    
	int n_events;
	int neutrality;
	int subdivision;
	int changesizepop;
	int split_pop;
	int pselection;

	/*ML theta estimation*/
    long unsigned n_iter2;
    long int seed2;
	int steps;
	double thetamax_bp;
	double thetamin_bp;

	/*ml estimation*/
	double theta1_ml[1];
	double ml1;
	double p1;
	double LR12;
	
	double theta2_ml[2];
	int nloci2m[2];
	double ml2[3];
	double p2;
	double LR23;
	
	double theta3_ml[3];
	int nloci3m[3];
	double ml3[4];
	double p3;
	double LR34;
	
	double theta4_ml[4];
	int nloci4m[4];
	double ml4[5];
	double p4;
	double LR4l;

	double mli_theta;
	int bestmodel;
};

struct var2b/*structure for working in ms2.c and ML. Data for one locus. NLOCI vectors*/
{
    int nsam;
    double r;
    unsigned long nsites;
    unsigned long nmhits;
	
    double theta;
    int nloci;
    double ratio_sv;
	int S;
	double factor_chrn;

	/*each locus has its own information about its evolutive process*/
    int definedmodel;
	
	int neutrality;

    int subdivision;/*but the samples are coming from a SINGLE population*/
    int npop;
    double mig_rate;

    double *factor_pop;
    double *freq;
   
    int split_pop;
    double time_split;
    double time_scoal;
    double factor_anc;
    
    int changesizepop;
	int iflogistic;
    int nintn;
    double *nrec;
    double *npast;
    double *tpast;
	double ts;       

    int pselection;
    int npsloci;	/*number of loci under selection*/
    long int pop_size;	/*effective population size N*/
    double pop_sel;	/*4Ns for each locus*/
    long int sel_nt;
    /*nt position of the selected mutation. 0 is the left value of the sequence. negative is allowed!!*/
    double sinit;/*time in 4N generations the selection process ended. negative value means sel is not finished..*/

	/*ml estimation*/
	double thetaml[6]; /*0 nothing, 1 theta1, 2 theta2, ... , 5 thetaallloci*/
	double mltheta;
};

struct var2/*structure for working in ms2.c. Data for one locus. 1 vector*/
{
    long unsigned howmany;
    
    int nsam;
    double r;
    unsigned long nsites;
    double theta;
	double factor_chrn;

    int nloci;
    int tloci;
    int mhits;
    double ratio_sv;
    double T_out;
     
    long int npop;
    double migrate;
    int *config;
    double *factor_pop;
   
    int nintn;
	int iflogistic;
    double *nrec;
    double *npast;
    double *tpast;
	double ts;    

    int split_pop;
    double time_split;
    double time_scoal;
    double factor_anc;
    double *freq;
    
    int ifselection;
    long int pop_size;
    double pop_sel;
    double sel_nt;
    double sinit;

    double tlimit;    
    int Sout;
};

struct statistisim /*structure made for keeping all results from coalescent simulations. THERE ARE NLOCI x NITER VECTORS*/
{
    long int nsites;
	double factor_chrn;
    int nsamples;
    int biallsites;
    int fixed;
    int nhapl;
    double nhaplsam;
    double hapldiv;
	double ewtest;
    double za;
    int b;
    int q;
	double Rvpi;
	int Rm;
	
	int mhsites;
	
    double theta_wat;
    double theta_taj;
    double theta_fw;
    double theta_fuli;
    double theta_fulin;
    double theta_L;
    double ndivergence;

    double theta_wat_nut;
    double theta_taj_nut;
    double theta_fw_nut;
    double theta_fuli_nut;
    double theta_fulin_nut;
    double theta_L_nut;
    double ndivergence_nut;

    double tajimaD;
    double fuliD;
    double fuliDn;
    double fuliF;
    double fuliFn;
    double fuFs;
    double faywuH;
    double faywuHo;
    double zengE;
    double rZA;
    double wB;
    double wQ;
    double R2;
    double hka;
    double hka_theta;
};

struct statistisimmuloc /*structure that keeps all multilocus results. Sum and Sum squared. THERE ARE NITER VECTORS*/
{
    long int Snsites;
    long int Sbiallsites;
    long int Sfixed;
    double Snhapl;
    double Snhaplsam;
    double Shapldiv;
	double Sewtest;
    double Sza;
    long int Sb;
    long int Sq;
	double SRvpi;
	long int SRm;
    double Stheta_wat;
    double Stheta_taj;
    double Stheta_fw;
    double Stheta_fuli;
    double Stheta_fulin;
    double Stheta_L;
    double Sndivergence;
    double Sndivergencejc;
    double Stheta_wat_nut;
    double Stheta_taj_nut;
    double Stheta_fw_nut;
    double Stheta_fuli_nut;
    double Stheta_fulin_nut;
    double Stheta_L_nut;
    double Sndivergence_nut;
    double Sndivergencejc_nut;
    
    long int S2biallsites;
    long int S2fixed;
    double S2nhapl;
    double S2nhaplsam;
    double S2hapldiv;
	double S2ewtest;
    double S2za;
    long int S2b;
    long int S2q;
	double S2Rvpi;
	long int S2Rm;
    double S2theta_wat;
    double S2theta_taj;
    double S2theta_fw;
    double S2theta_fuli;
    double S2theta_fulin;
    double S2theta_L;
    double S2ndivergence;
    double S2ndivergencejc;
    double S2theta_wat_nut;
    double S2theta_taj_nut;
    double S2theta_fw_nut;
    double S2theta_fuli_nut;
    double S2theta_fulin_nut;
    double S2theta_L_nut;
    double S2ndivergence_nut;
    double S2ndivergencejc_nut;

    double StajimaD;
    double SfuliD;
    double SfuliDn;
    double SfuliF;
    double SfuliFn;
    double SfuFs;
    double SfaywuH;
    double SfaywuHo;
    double SzengE;
    double SrZA;
    double SwB;
    double SwQ;
    double SR2;
    double S2tajimaD;
    double S2fuliD;
    double S2fuliDn;
    double S2fuliF;
    double S2fuliFn;
    double S2fuFs;
    double S2faywuH;
    double S2faywuHo;
    double S2zengE;
    double S2rZA;
    double S2wB;
    double S2wQ;
    double S2R2;

    double Shka;
    double hka_T;

    int nltajD;
    int nlflD;
    int nlflDn;
    int nlflF;
    int nlflFn;
    int nlFs;
    int nlH;
    int nlHo;
    int nlE;
    int nlZ;
    int nlB;
    int nlQ;
    int nldivjc;
    int nlR2;
};

struct horizontalstatsml /*structure that keeps the average for each loci from each statistic. THERE ARE NLOCI VECTORs*/
{
    double biallsites;
    double fixed;
    double nhapl;
    double nhaplsam;
    double hapldiv;
	double ewtest;
    double za;
    double b;
    double q;
	double Rvpi;
	double Rm;
    
	double theta_wat;
    double theta_taj;
    double theta_fw;
    double theta_fuli;
    double theta_fulin;
    double theta_L;
    double ndivergence;
	
    double theta_wat_nut;
    double theta_taj_nut;
    double theta_fw_nut;
    double theta_fuli_nut;
    double theta_fulin_nut;
    double theta_L_nut;
    double ndivergence_nut;
	
    double tajimaD;
    double fuliD;
    double fuliDn;
    double fuliF;
    double fuliFn;
    double fuFs;
    double faywuH;
    double faywuHo;
    double zengE;
    double rZA;
    double wB;
    double wQ;
    double R2;
};

/*Structure that keeps the LRT distributions. There are n_iter2 vectors*/
struct LRTdist
{
	double LRT12;
	double LRT23;
	double LRT34;
	double LRT4l;
};

/*Structure that keeps the simulated theta distributions for each model and loci. There are n_iter2 * NLOCI vectors*/
/*
struct MLthetasim
{
	double thetaml1;
	double thetaml2;
	double thetaml3;
	double thetaml4;
};
*/

