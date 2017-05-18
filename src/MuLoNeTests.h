/*
 *  MuLoNeTests.h
 *  MuLoNeTests
 *
 *  Created by sebas on Tue Feb 25 2003.
 *
 */

#define COMMAND_LINE 1
#define UNIX_CL 1
#define DOS_CL 0
#define GFF_ACTIVE 1
#define MANVa "MANVa: Software for Multilocus Analyses of Nucleotide Variation within Population. v0.9892beta (20120515)\n       Ramos-Onsins, S. E., Windsor, A. and Mitchell-Olds, T.\n"
#define MLtheta_open 0
#define CALCREC 0

#define _CRT_SECURE_NO_DEPRECATE

struct globvar 
{
    int observed_data;
    int n_loci;
    int outgroup;
    char name_outgroup[50];
    char name_excluded[50];
    char name_ingroups[50];
    /*int maxloc;*/

    int montecarlo_sim; /*once mc simulations are done*/
    int nlocimat;/*=nlocimatt. nloci in mc simulations and in thetaml*/
	long int nitermat;/*niter in mc simulations. once is finished*/
    int dataforsim;/*indicator mc simulated parameters are included. ready to simulations */
    int onlymulo;/*in case only multilocus test can be performed*/

    int dataindata;/*data included in data mc =1 and in ML =2 simulations, (both=3). Not necessary a model*/	
    int dataequalsim;/*indicate same values in data than in mc simulations or thetaml estimations*/

	int ml0done;/*once estimation of parameters is done*/
    int nlocimatt;/*=nlocimat. nloci in thetaml and in mc sims. */
	long int nitermatt;/*niter in estimation theta parameters section*/
	int dataforsimt;/*indicator theta-ml simulated parameters are included. ready to start */
	int mltheta;/*0 not done, 1 2 3 4 5(all) means the best model for theta (no recombination)*/
	
	int gfffiles; /*0 when no defined, -1 if not gff, 1 if gff*/
	char subset_positions[50];/*if gffactive, define the subset of positions that are analyzed. Include in headers*/
	int ifgencode; /*if genetic code must be calculated 1, else 0*/
	char code_name[50];/*name of the genetic code*/
	char genetic_code[64];/*aa for each triplet*/
};

struct statistics /*structure that contains all the statistics from observed data. ONE VECTOR FOR EACH LOCUS.*/
{
    char gene[21];
    int nsamples;
    int noutgroups;
    double nsites;
	double factor_chrn;

	int nmhits;
    
	int transitions;
	int transversions;
	
    int biallsites;
    int biallsitesn;
    int shared;
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
    double theta_wat;
    double theta_watn;
    double theta_taj;
    double theta_tajn;
    double theta_fw;
    double theta_fuli;
    double theta_fulin;
	double theta_L;
	double theta_mlS;
    double ndivergence;
    
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
	double Sobshka;
	double Dobshka;
	double Sexphka;
	double Dexphka;
    double hka_theta;

};

struct statmulo 
/*structure that contains sum and sum square for all statistics (calculate avg and variance for all loci).*/
{
    int nloci;
    double Snsites;
	
	long int Snmhits;
	    
    long int Sbiallsites;
    long int Sbiallsitesn;
    long int Sshared;
    long int Sfixed;
    long int S2biallsites;
    long int S2biallsitesn;
    long int S2shared;
    long int S2fixed;
	
	long int Stransitions;
	long int Stransversions;

    long int Snhapl;
    double Snhaplsam;
    double Shapldiv;
	double Sewtest;
    double Sza;
    long int Sb;
    long int Sq;   
    double SRvpi;
    long int SRm;
    long int S2nhapl;
    double S2nhaplsam;
    double S2hapldiv;
	double S2ewtest;
    double S2za;
    long int S2b;
    long int S2q;
	double S2Rvpi;
	long int S2Rm;
	
    double Stheta_wat;
    double Stheta_watn;
    double Stheta_taj;
    double Stheta_tajn;
    double Stheta_fw;
    double Stheta_fuli;
    double Stheta_fulin;
    double Stheta_L;
    double Sndivergence;
    double Sndivergencejc;
    double S2theta_wat;
    double S2theta_watn;
    double S2theta_taj;
    double S2theta_tajn;
    double S2theta_fw;
    double S2theta_fuli;
    double S2theta_fulin;
    double S2theta_L;
    double S2ndivergence;
    double S2ndivergencejc;

    double Stheta_wat_nut;
    double Stheta_watn_nut;
    double Stheta_taj_nut;
    double Stheta_tajn_nut;
    double Stheta_fw_nut;
    double Stheta_fuli_nut;
    double Stheta_fulin_nut;
    double Stheta_L_nut;
    double Sndivergence_nut;
    double Sndivergencejc_nut;
    double S2theta_wat_nut;
    double S2theta_watn_nut;
    double S2theta_taj_nut;
    double S2theta_tajn_nut;
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
    int nlR2;
    int nldivjc;
    int nlhka;
    
    double Shka;
    double signif_hka;
    double hka_T;
	
};

struct lianinput 
{
	int nsam;
	int nloci;
	char **namesam;
	int **samhap;
};

