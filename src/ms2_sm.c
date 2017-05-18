/*HUDSON'S MS PROGRAM (modified)*/
/*also Wall's Rm function*/

#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SITESINC 100

unsigned long maxsites = 1000;	/* la llargada de la regió */
struct node {
    int abv;
    int ndes;
    double time;
};
struct segl {
    long int beg;
    struct node *ptree;
    long int next;
};

struct dnapar {
    double k;
    unsigned long S;
    int nhapl;
    int *fhapl;
    int B1;
    int Q1;
    int *freq;
    int *unic;
	double C;
	int Rm;
	double thetaL;
	int mhsites;
};

double **coef = NULL;
struct dnapar *neutpar = NULL;    
unsigned long *posit;	/*important! externes perque al reallocar memoria no localitza la nova posicio. */
                        /*Tambe es fa posant la posicio de memoria*/
char **list;

int ms(struct var2 *inputp,struct statistisim **matrix_test,struct statistisimmuloc **matrixmlsim,struct horizontalstatsml **avgstatloci,int onlymulo,FILE *file_output)
{
    unsigned long count;
    unsigned long segsites;
    int i;

    char **cmatrix(int,unsigned long);
    unsigned long gensam(long int,int,int *,unsigned long,double/*,unsigned long */,double,double,double,double,int,
        /*unsigned long,*/double * /*,double **/,int,double,double,double,long int,double,int *,int,double *,double *,double *,
        int, double, double, double, double *, double,int,double,double,double);
    int mod_mhits(unsigned long, struct var2 **);
    double ran1(void);
    
    void init_coef(double *,int);
    int calc_neutpar(double,int,unsigned long,struct var2 *,struct dnapar *);
    double tajima_d(double,int, double *);
    double Fs(int,double,int);
    double fl_d(/*int,*/int,int,double *);
    double fl_f(/*int,*/int,int,double,double *);
    double fl_d2(int,double,int,double *);
    double fl_f2(int,double,int,double,double *);
    double fay_wu(int,int *,double);
    double fay_wu_original(int,int *,double);
    double ZnA_(int,unsigned long,struct var2 *);
    double koutgJC(int,int,int *,unsigned long);
    double koutg(int,int,int * /*,unsigned long*/);
    double R2(int *,double,int,int);
	double testHap(int, int *);
	double E_zeng(int,double,double,double,double *); /*outgroup*/
	double EWtest(int, int *);
	  
    double stf;
    /*double lengtht;*/
    
    static double counterp,restp;
    static long int p;
	int x;

    if((*inputp).nloci == 0) {
		counterp = (double)(((*inputp).howmany) * (*inputp).tloci)/(double)50;
		p = 1;
		restp = (double)0;
    }
    
    if((*inputp).theta >= 0.0) {        
        if((list = cmatrix((*inputp).nsam,maxsites+1)) == 0) return 0;
        if(!(posit = (unsigned long *)malloc((unsigned long)(maxsites*sizeof(unsigned long))))) {
            puts("ms error.1");
            if(file_output) fputs("Error in ms",file_output);
            return 0;
        }
    }
    if((neutpar = (struct dnapar *)calloc(1,sizeof(struct dnapar))) == NULL) {
        perror("calloc error ms.1d1");
        exit(1);
    }
    if((coef = (double **)calloc(1,sizeof(double *))) == NULL) {
        puts("calloc error ms.1d1");
        if(file_output) fputs("Error in ms",file_output);
        return 0;
    }
    for(i=0;i<1;i++) {
        if((coef[i] = (double *)calloc(12,sizeof(double))) == NULL) {
            puts("calloc error ms.1d1");
            if(file_output) fputs("Error in ms",file_output);
            return 0;
        }
    }
    init_coef(coef[0],(*inputp).nsam);

    if((neutpar[0].freq = (int *)calloc((*inputp).nsam,sizeof(int))) == NULL) {
        puts("calloc error ms.1d3");
        if(file_output) fputs("Error in ms",file_output);
        return 0;
    }
    if((neutpar[0].fhapl = (int *)calloc((*inputp).nsam,sizeof(int))) == NULL) {
        puts("calloc error ms.1d3");
        if(file_output) fputs("Error in ms",file_output);
        return 0;
    }
    if((neutpar[0].unic = (int *)calloc((*inputp).nsam,sizeof(int))) == NULL) {
        puts("calloc error ms.1d3");
        if(file_output) fputs("Error in ms",file_output);
        return 0;
    }
    
    /************************************************  ITERATIONS ****************************************************/
    count = 0; 
    while(((*inputp).howmany) - count++) {    
        if((segsites = gensam((long int)(*inputp).npop,(int)(*inputp).nsam,(int *)(*inputp).config,(unsigned long)(*inputp).nsites,
            (double)(*inputp).theta/*,(unsigned long)0*/,(double)(*inputp).r,(double)0.,(double)0.,
            (double)(*inputp).migrate,(int)(*inputp).mhits/*,(unsigned long)count*/,(double *)(*inputp).factor_pop, 
            /*(double *)&lengtht,*/(int)(*inputp).ifselection, 
            (double)(*inputp).pop_sel,(double)(*inputp).sinit,(double)(*inputp).pop_size,(long int)(*inputp).sel_nt,(double)(*inputp).T_out,(int *)&(*inputp).Sout,
            (int)(*inputp).nintn,(double *)(*inputp).nrec,(double *)(*inputp).npast,(double *)(*inputp).tpast,
            (int)(*inputp).split_pop,(double)(*inputp).time_split,(double)(*inputp).time_scoal,(double)(*inputp).factor_anc,(double *)(*inputp).freq,
            (double)(*inputp).tlimit,
			(int)(*inputp).iflogistic,(double)(*inputp).ts,(double)(*inputp).factor_chrn,(double)(*inputp).ratio_sv)) == 100000) return 0;                    
        /******** mhits ******************/
        if((*inputp).mhits) if(mod_mhits(segsites,&inputp) == 0) return 0;
        /************ calculate statistics and neutrality test and include in matrix_test *****************/
        if(calc_neutpar((double)(*inputp).r,0,segsites,inputp,neutpar+0) == 0) return 0;
        if(onlymulo == 0) {
            matrix_test[(*inputp).nloci][count-1].nsites = (long int)inputp[0].nsites - (long int)neutpar[0].mhsites;
            matrix_test[(*inputp).nloci][count-1].mhsites = (long int)neutpar[0].mhsites;
            matrix_test[(*inputp).nloci][count-1].factor_chrn = (double)inputp[0].factor_chrn;
            matrix_test[(*inputp).nloci][count-1].nsamples = (int)inputp[0].nsam;
            matrix_test[(*inputp).nloci][count-1].biallsites = (int)neutpar[0].S;
            if((*inputp).T_out > 0.) matrix_test[(*inputp).nloci][count-1].fixed = (int)(*inputp).Sout;
            else matrix_test[(*inputp).nloci][count-1].fixed = (int)-10000;
            matrix_test[(*inputp).nloci][count-1].nhapl = (int)neutpar[0].nhapl;
            matrix_test[(*inputp).nloci][count-1].nhaplsam = (double)neutpar[0].nhapl/(double)(*inputp).nsam;
            matrix_test[(*inputp).nloci][count-1].hapldiv = (double)testHap((int)(*inputp).nsam,(int *)neutpar[0].fhapl);
            matrix_test[(*inputp).nloci][count-1].ewtest = (double)EWtest((int)(*inputp).nsam,(int *)neutpar[0].fhapl);
            matrix_test[(*inputp).nloci][count-1].za = (double)ZnA_(0,segsites,inputp);
            matrix_test[(*inputp).nloci][count-1].b = (int)neutpar[0].B1;
            matrix_test[(*inputp).nloci][count-1].q = (int)neutpar[0].Q1;
            matrix_test[(*inputp).nloci][count-1].Rvpi = neutpar[0].C;
            matrix_test[(*inputp).nloci][count-1].Rm = (int)neutpar[0].Rm;
            
			if(neutpar[0].S != -10000) 
				matrix_test[(*inputp).nloci][count-1].theta_wat = (double)neutpar[0].S/coef[0][0];
			else 
				matrix_test[(*inputp).nloci][count-1].theta_wat = -10000.;
            matrix_test[(*inputp).nloci][count-1].theta_taj = (double)neutpar[0].k;
            if(neutpar[0].k != -10000) 
				matrix_test[(*inputp).nloci][count-1].theta_fw = (double)neutpar[0].k - (double)fay_wu_original((*inputp).nsam,neutpar[0].freq,neutpar[0].k);
			else
				matrix_test[(*inputp).nloci][count-1].theta_fw = -10000.;
            if(neutpar[0].S != -10000) 
				matrix_test[(*inputp).nloci][count-1].theta_fuli = (double)(neutpar[0].freq[1]);
			else
				matrix_test[(*inputp).nloci][count-1].theta_fuli = -10000.;
            if(neutpar[0].S != -10000) 
				matrix_test[(*inputp).nloci][count-1].theta_fulin = (double)(neutpar[0].freq[1]+neutpar[0].freq[(*inputp).nsam-1])*((double)(*inputp).nsam-(double)1.0)/(double)(*inputp).nsam;
            else
				matrix_test[(*inputp).nloci][count-1].theta_fulin = -10000.;
			matrix_test[(*inputp).nloci][count-1].theta_L = (double)(neutpar[0].thetaL);
            if((*inputp).T_out > 0.) 
                matrix_test[(*inputp).nloci][count-1].ndivergence = (double)koutg((*inputp).nsam,(*inputp).Sout,neutpar[0].freq/*,inputp[0].nsites*/);
            else matrix_test[(*inputp).nloci][count-1].ndivergence = (double)-10000.;
 
			if(neutpar[0].S != -10000) {
				matrix_test[(*inputp).nloci][count-1].theta_wat_nut = (double)matrix_test[(*inputp).nloci][count-1].theta_wat/(double)matrix_test[(*inputp).nloci][count-1].nsites;
				matrix_test[(*inputp).nloci][count-1].theta_taj_nut = (double)matrix_test[(*inputp).nloci][count-1].theta_taj/(double)matrix_test[(*inputp).nloci][count-1].nsites;
				matrix_test[(*inputp).nloci][count-1].theta_fw_nut  = (double)matrix_test[(*inputp).nloci][count-1].theta_fw /(double)matrix_test[(*inputp).nloci][count-1].nsites;
				matrix_test[(*inputp).nloci][count-1].theta_fuli_nut = (double)matrix_test[(*inputp).nloci][count-1].theta_fuli/(double)matrix_test[(*inputp).nloci][count-1].nsites;
				matrix_test[(*inputp).nloci][count-1].theta_fulin_nut = (double)matrix_test[(*inputp).nloci][count-1].theta_fulin/(double)matrix_test[(*inputp).nloci][count-1].nsites;
				matrix_test[(*inputp).nloci][count-1].theta_L_nut = (double)matrix_test[(*inputp).nloci][count-1].theta_L/(double)matrix_test[(*inputp).nloci][count-1].nsites;
			}
			else {
				matrix_test[(*inputp).nloci][count-1].theta_wat_nut = -10000.;
				matrix_test[(*inputp).nloci][count-1].theta_taj_nut = -10000.;
				matrix_test[(*inputp).nloci][count-1].theta_fw_nut  = -10000.;
				matrix_test[(*inputp).nloci][count-1].theta_fuli_nut = -10000.;
				matrix_test[(*inputp).nloci][count-1].theta_fulin_nut = -10000.;
				matrix_test[(*inputp).nloci][count-1].theta_L_nut = -10000.;
			}
            if((*inputp).mhits == 1) 
                matrix_test[(*inputp).nloci][count-1].ndivergence_nut = (double)matrix_test[(*inputp).nloci][count-1].ndivergence/(double)matrix_test[(*inputp).nloci][count-1].nsites;
            else matrix_test[(*inputp).nloci][count-1].ndivergence_nut = (double)-10000.;
			
			matrix_test[(*inputp).nloci][count-1].tajimaD = (double)tajima_d(neutpar[0].k,neutpar[0].S,coef[0]);
            matrix_test[(*inputp).nloci][count-1].fuliD = (double)fl_d(/*(*inputp).nsam,*/neutpar[0].freq[1],neutpar[0].S,coef[0]);
            matrix_test[(*inputp).nloci][count-1].fuliF = (double)fl_f(/*(*inputp).nsam,*/neutpar[0].freq[1],neutpar[0].S,neutpar[0].k,coef[0]);
            matrix_test[(*inputp).nloci][count-1].fuliDn = (double)fl_d2((*inputp).nsam,(double)neutpar[0].freq[1]+neutpar[0].freq[(*inputp).nsam-1],neutpar[0].S,coef[0]);
            matrix_test[(*inputp).nloci][count-1].fuliFn = (double)fl_f2((*inputp).nsam,(double)neutpar[0].freq[1]+neutpar[0].freq[(*inputp).nsam-1],neutpar[0].S, neutpar[0].k,coef[0]);
            matrix_test[(*inputp).nloci][count-1].faywuH = (double)fay_wu((*inputp).nsam,neutpar[0].freq,neutpar[0].k);
            if(neutpar[0].k != -10000) matrix_test[(*inputp).nloci][count-1].faywuHo = (double)neutpar[0].k - (double)matrix_test[(*inputp).nloci][count-1].theta_fw;
			else matrix_test[(*inputp).nloci][count-1].faywuHo = -10000.;
            matrix_test[(*inputp).nloci][count-1].zengE = (double)E_zeng((*inputp).nsam,(double)neutpar[0].thetaL,(double)neutpar[0].S/coef[0][0],(double)neutpar[0].S,coef[0]);
            matrix_test[(*inputp).nloci][count-1].fuFs = (double)Fs((*inputp).nsam,(double)neutpar[0].k,neutpar[0].nhapl);
            matrix_test[(*inputp).nloci][count-1].R2 = (double)R2(neutpar[0].unic,neutpar[0].k,(*inputp).nsam,(int)neutpar[0].S);
            if(neutpar[0].S > 1) {
                matrix_test[(*inputp).nloci][count-1].wB = (double)neutpar[0].B1/((double)neutpar[0].S - (double)1.0);
                matrix_test[(*inputp).nloci][count-1].wQ = (double)neutpar[0].Q1/((double)neutpar[0].S);
                matrix_test[(*inputp).nloci][count-1].rZA = (double)ZnA_(0,segsites,inputp)/((double)neutpar[0].S - (double)1.0);
            }
            else {
                matrix_test[(*inputp).nloci][count-1].wB = (double)-10000;
                matrix_test[(*inputp).nloci][count-1].wQ = (double)-10000;
                matrix_test[(*inputp).nloci][count-1].rZA = (double)-10000;
            }
        }
        /*multilocus matrix*/
        matrixmlsim[0][count-1].Snsites += (long int)inputp[0].nsites - (long int)neutpar[0].mhsites;
        if(neutpar[0].S != -10000) {
			matrixmlsim[0][count-1].Sbiallsites += (long int)neutpar[0].S;
			matrixmlsim[0][count-1].S2biallsites += (long int)neutpar[0].S*(long int)neutpar[0].S;
			matrixmlsim[0][count-1].Sfixed += (long int)(*inputp).Sout;
			matrixmlsim[0][count-1].S2fixed += (long int)(*inputp).Sout*(long int)(*inputp).Sout;
			matrixmlsim[0][count-1].Snhapl += (long int)neutpar[0].nhapl;
			matrixmlsim[0][count-1].S2nhapl += (long int)neutpar[0].nhapl*(long int)neutpar[0].nhapl;
			matrixmlsim[0][count-1].Snhaplsam += (double)neutpar[0].nhapl/(double)inputp[0].nsam;
			matrixmlsim[0][count-1].S2nhaplsam += (double)neutpar[0].nhapl/(double)inputp[0].nsam*(double)neutpar[0].nhapl/(double)inputp[0].nsam;
			matrixmlsim[0][count-1].Shapldiv += (double)testHap((int)inputp[0].nsam,(int *)neutpar[0].fhapl);
			matrixmlsim[0][count-1].S2hapldiv += (double)testHap((int)inputp[0].nsam,(int *)neutpar[0].fhapl)*(double)testHap((int)inputp[0].nsam,(int *)neutpar[0].fhapl);
			matrixmlsim[0][count-1].Sewtest += (double)EWtest((int)(*inputp).nsam,(int *)neutpar[0].fhapl);
			matrixmlsim[0][count-1].S2ewtest += (double)EWtest((int)(*inputp).nsam,(int *)neutpar[0].fhapl) * (double)EWtest((int)(*inputp).nsam,(int *)neutpar[0].fhapl);
			matrixmlsim[0][count-1].Sza += (double)ZnA_(0,segsites,inputp);
			matrixmlsim[0][count-1].S2za += (double)ZnA_(0,segsites,inputp)*(double)ZnA_(0,segsites,inputp);
			matrixmlsim[0][count-1].Sb += (long int)neutpar[0].B1;
			matrixmlsim[0][count-1].S2b += (long int)neutpar[0].B1*(long int)neutpar[0].B1;
			matrixmlsim[0][count-1].Sq += (long int)neutpar[0].Q1;
			matrixmlsim[0][count-1].S2q += (long int)neutpar[0].Q1*(long int)neutpar[0].Q1;
			matrixmlsim[0][count-1].SRvpi += neutpar[0].C;
			matrixmlsim[0][count-1].S2Rvpi += neutpar[0].C*neutpar[0].C;
			matrixmlsim[0][count-1].SRm += (long int)neutpar[0].Rm;
			matrixmlsim[0][count-1].S2Rm += (long int)neutpar[0].Rm*(long int)neutpar[0].Rm;

			if((double)neutpar[0].S != (double)-10000) {
				stf = (double)neutpar[0].S/coef[0][0] * ((double)1/(double)(*inputp).factor_chrn);
				matrixmlsim[0][count-1].Stheta_wat += stf;
				matrixmlsim[0][count-1].S2theta_wat += stf*stf;
				matrixmlsim[0][count-1].Stheta_wat_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
				matrixmlsim[0][count-1].S2theta_wat_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites) * stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
			}
			if((double)neutpar[0].k != (double)-10000) {
				stf = (double)neutpar[0].k * ((double)1/(double)(*inputp).factor_chrn);
				matrixmlsim[0][count-1].Stheta_taj += stf;
				matrixmlsim[0][count-1].S2theta_taj += stf*stf;
				matrixmlsim[0][count-1].Stheta_taj_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
				matrixmlsim[0][count-1].S2theta_taj_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites) * stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
			}
			if((double)neutpar[0].k != (double)-10000) {
				stf = ((double)neutpar[0].k - (double)fay_wu((*inputp).nsam,neutpar[0].freq,neutpar[0].k) * ((double)1/(double)(*inputp).factor_chrn));
				matrixmlsim[0][count-1].Stheta_fw += stf;
				matrixmlsim[0][count-1].S2theta_fw += stf*stf;
				matrixmlsim[0][count-1].Stheta_fw_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
				matrixmlsim[0][count-1].S2theta_fw_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites) * stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
			}
			if((double)(neutpar[0].freq[1]) != (double)-10000) {
				stf = (double)(neutpar[0].freq[1])* ((double)1/(double)(*inputp).factor_chrn);
				matrixmlsim[0][count-1].Stheta_fuli += stf;
				matrixmlsim[0][count-1].S2theta_fuli += stf*stf;
				matrixmlsim[0][count-1].Stheta_fuli_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
				matrixmlsim[0][count-1].S2theta_fuli_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites) * stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
			}
			if((double)neutpar[0].freq[1] != (double)-10000) {
				stf = (double)(neutpar[0].freq[1]+neutpar[0].freq[(*inputp).nsam-1])*((double)(*inputp).nsam-(double)1.0)/(double)(*inputp).nsam;
				matrixmlsim[0][count-1].Stheta_fulin += stf;
				matrixmlsim[0][count-1].S2theta_fulin += stf*stf;
				matrixmlsim[0][count-1].Stheta_fulin_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
				matrixmlsim[0][count-1].S2theta_fulin_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites) * stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
			}
			if(((double)neutpar[0].thetaL) != (double)-10000) {
				stf = (double)neutpar[0].thetaL * ((double)1/(double)(*inputp).factor_chrn);
				matrixmlsim[0][count-1].Stheta_L += stf;
				matrixmlsim[0][count-1].S2theta_L += stf*stf;
				matrixmlsim[0][count-1].Stheta_L_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
				matrixmlsim[0][count-1].S2theta_L_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites) * stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
			}
		}
        if((*inputp).T_out > 0.) {
            if((stf = (double)koutg((*inputp).nsam,(*inputp).Sout,neutpar[0].freq/*,inputp[0].nsites*/)) != (double)-10000) {
                matrixmlsim[0][count-1].Sndivergence += stf;
                matrixmlsim[0][count-1].S2ndivergence += stf*stf;
                matrixmlsim[0][count-1].Sndivergence_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
                matrixmlsim[0][count-1].S2ndivergence_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites) * stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
            }
            if((stf = (double)koutgJC((*inputp).nsam,(*inputp).Sout,neutpar[0].freq,inputp[0].nsites)) != (double)-10000) {
                matrixmlsim[0][count-1].Sndivergencejc += stf;
                matrixmlsim[0][count-1].S2ndivergencejc += stf*stf;
                matrixmlsim[0][count-1].Sndivergencejc_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
                matrixmlsim[0][count-1].S2ndivergencejc_nut += stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites) * stf/(double)((long int)inputp[0].nsites - (long int)neutpar[0].mhsites);
                matrixmlsim[0][count-1].nldivjc += 1;
            }
        }
        if((stf = (double)tajima_d(neutpar[0].k,neutpar[0].S,coef[0])) != (double)-10000) {
            matrixmlsim[0][count-1].StajimaD += stf;
            matrixmlsim[0][count-1].S2tajimaD += stf*stf;
            matrixmlsim[0][count-1].nltajD += 1;
        }
        if((stf = (double)fl_d2((*inputp).nsam,(double)(neutpar[0].freq[1]+neutpar[0].freq[(*inputp).nsam-1]),neutpar[0].S,coef[0])) != (double)-10000) {
            matrixmlsim[0][count-1].SfuliDn += stf;
            matrixmlsim[0][count-1].S2fuliDn += stf*stf;
            matrixmlsim[0][count-1].nlflDn += 1;
        }
        if((stf = (double)fl_d(/*(*inputp).nsam,*/neutpar[0].freq[1],neutpar[0].S,coef[0])) != (double)-10000) {
            matrixmlsim[0][count-1].SfuliD += stf;
            matrixmlsim[0][count-1].S2fuliD += stf*stf;
            matrixmlsim[0][count-1].nlflD += 1;
        }
        if((stf = (double)fl_f2((*inputp).nsam,(double)(neutpar[0].freq[1]+neutpar[0].freq[(*inputp).nsam-1]),neutpar[0].S,neutpar[0].k,coef[0])) != (double)-10000) {
            matrixmlsim[0][count-1].SfuliFn += stf;
            matrixmlsim[0][count-1].S2fuliFn += stf*stf;
            matrixmlsim[0][count-1].nlflFn += 1;
        }
        if((stf = (double)fl_f(/*(*inputp).nsam,*/neutpar[0].freq[1],neutpar[0].S, neutpar[0].k,coef[0])) != (double)-10000) {
            matrixmlsim[0][count-1].SfuliF += stf;
            matrixmlsim[0][count-1].S2fuliF += stf*stf;
            matrixmlsim[0][count-1].nlflF += 1;
        }
        if((stf = (double)Fs((*inputp).nsam,neutpar[0].k,neutpar[0].nhapl)) != (double)-10000) {
            matrixmlsim[0][count-1].SfuFs += stf;
            matrixmlsim[0][count-1].S2fuFs += stf*stf;
            matrixmlsim[0][count-1].nlFs += 1;
        }
        if(neutpar[0].k != -10000) {
			if((stf = (double)fay_wu((*inputp).nsam,neutpar[0].freq,neutpar[0].k)) != (double)-10000) {
				matrixmlsim[0][count-1].SfaywuH += stf;
				matrixmlsim[0][count-1].S2faywuH += stf*stf;
				matrixmlsim[0][count-1].nlH += 1;
			}
			if((stf = (double)neutpar[0].k - (double)matrix_test[(*inputp).nloci][count-1].theta_fw) != (double)-10000) {
				matrixmlsim[0][count-1].SfaywuHo += stf;
				matrixmlsim[0][count-1].S2faywuHo += stf*stf;
				matrixmlsim[0][count-1].nlHo += 1;
			}
		}
        if((stf = (double)E_zeng((*inputp).nsam,(double)neutpar[0].thetaL,(double)neutpar[0].S/coef[0][0],(double)neutpar[0].S,coef[0])) != (double)-10000) {
            matrixmlsim[0][count-1].SzengE += stf;
            matrixmlsim[0][count-1].S2zengE += stf*stf;
            matrixmlsim[0][count-1].nlE += 1;
        }
        if(neutpar[0].S > 1) {
            if((double)ZnA_(0,segsites,inputp) != (double)-10000) {
				stf = (double)ZnA_(0,segsites,inputp)/((double)neutpar[0].S - (double)1.0);
                matrixmlsim[0][count-1].SrZA += stf;
                matrixmlsim[0][count-1].S2rZA += stf*stf;
                matrixmlsim[0][count-1].nlZ += 1;
            }
            if(((double)neutpar[0].B1) != (double)-10000) {
                stf = (double)neutpar[0].B1/((double)neutpar[0].S - (double)1.0);
				matrixmlsim[0][count-1].SwB += stf;
                matrixmlsim[0][count-1].S2wB += stf*stf;
                matrixmlsim[0][count-1].nlB += 1;
            }
            if((double)neutpar[0].Q1 != (double)-10000) {
				stf = (double)neutpar[0].Q1/((double)neutpar[0].S);
                matrixmlsim[0][count-1].SwQ += stf;
                matrixmlsim[0][count-1].S2wQ += stf*stf;
                matrixmlsim[0][count-1].nlQ += 1;
            }
        }
        if((stf = (double)R2(neutpar[0].unic,neutpar[0].k,(*inputp).nsam,(int)neutpar[0].S)) != (double)-10000) {
            matrixmlsim[0][count-1].SR2 += stf;
            matrixmlsim[0][count-1].S2R2 += stf*stf;
            matrixmlsim[0][count-1].nlR2 += 1;
        }

        if(onlymulo == 1) {
            /*calculate the average in case onlymulo = 0*/
            if(neutpar[0].S != -10000) {
				avgstatloci[0][(*inputp).nloci].biallsites += (double)neutpar[0].S;
				avgstatloci[0][(*inputp).nloci].fixed += (double)(*inputp).Sout;
				avgstatloci[0][(*inputp).nloci].nhapl += (double)neutpar[0].nhapl;
				avgstatloci[0][(*inputp).nloci].nhaplsam += (double)neutpar[0].nhapl/(double)(*inputp).nsam;
				avgstatloci[0][(*inputp).nloci].hapldiv += (double)testHap((int)(*inputp).nsam,(int *)neutpar[0].fhapl);
				avgstatloci[0][(*inputp).nloci].ewtest += (double)EWtest((int)(*inputp).nsam,(int *)neutpar[0].fhapl);
				avgstatloci[0][(*inputp).nloci].za += (double)ZnA_(0,segsites,inputp);
				avgstatloci[0][(*inputp).nloci].b += (double)neutpar[0].B1;
				avgstatloci[0][(*inputp).nloci].q += (double)neutpar[0].Q1;
				avgstatloci[0][(*inputp).nloci].Rvpi += (double)neutpar[0].C;
				avgstatloci[0][(*inputp).nloci].Rm += (double)neutpar[0].Rm;
				avgstatloci[0][(*inputp).nloci].theta_wat += (double)neutpar[0].S/coef[0][0] * ((double)1/(double)(*inputp).factor_chrn);
				avgstatloci[0][(*inputp).nloci].theta_taj += (double)neutpar[0].k * ((double)1/(double)(*inputp).factor_chrn);
				avgstatloci[0][(*inputp).nloci].theta_fw += ((double)neutpar[0].k - (double)fay_wu((*inputp).nsam,neutpar[0].freq,neutpar[0].k)) 
															* ((double)1/(double)(*inputp).factor_chrn);
				avgstatloci[0][(*inputp).nloci].theta_fuli += (double)(neutpar[0].freq[1]) * ((double)1/(double)(*inputp).factor_chrn);
				avgstatloci[0][(*inputp).nloci].theta_fulin += (double)(neutpar[0].freq[1]+neutpar[0].freq[(*inputp).nsam-1])*((double)(*inputp).nsam-(double)1.0)/(double)(*inputp).nsam
															* ((double)1/(double)(*inputp).factor_chrn);
				avgstatloci[0][(*inputp).nloci].theta_L += (double)(neutpar[0].thetaL);
            }
            if((*inputp).T_out > 0.) {
                avgstatloci[0][(*inputp).nloci].ndivergence += (double)koutg((*inputp).nsam,(*inputp).Sout,neutpar[0].freq/*,inputp[0].nsites*/);
            }
            
            if((stf = (double)tajima_d(neutpar[0].k,neutpar[0].S,coef[0])) != (double)-10000) {
                avgstatloci[0][(*inputp).nloci].tajimaD += stf;
            }
            if((stf = (double)fl_d2((*inputp).nsam,(double)(neutpar[0].freq[1]+neutpar[0].freq[(*inputp).nsam-1]),neutpar[0].S,coef[0])) != (double)-10000) {
                avgstatloci[0][(*inputp).nloci].fuliDn += stf;
            }
            if((stf = (double)fl_d(/*(*inputp).nsam,*/neutpar[0].freq[1],neutpar[0].S,coef[0])) != (double)-10000) {
                avgstatloci[0][(*inputp).nloci].fuliD += stf;
            }
            if((stf = (double)fl_f2((*inputp).nsam,(double)(neutpar[0].freq[1]+neutpar[0].freq[(*inputp).nsam-1]),neutpar[0].S,neutpar[0].k,coef[0])) != (double)-10000) {
                avgstatloci[0][(*inputp).nloci].fuliFn += stf;
            }
            if((stf = (double)fl_f(/*(*inputp).nsam,*/neutpar[0].freq[1],neutpar[0].S, neutpar[0].k,coef[0])) != (double)-10000) {
                avgstatloci[0][(*inputp).nloci].fuliF += stf;
            }
            if((stf = (double)Fs((*inputp).nsam,neutpar[0].k,neutpar[0].nhapl)) != (double)-10000) {
                avgstatloci[0][(*inputp).nloci].fuFs += stf;
            }
            if((stf = (double)fay_wu((*inputp).nsam,neutpar[0].freq,neutpar[0].k)) != (double)-10000) {
                avgstatloci[0][(*inputp).nloci].faywuH += stf;
            }
            if(neutpar[0].S != -10000) {
				if(((double)neutpar[0].k) != (double)-10000) {
					stf = (double)neutpar[0].k - (double)matrix_test[(*inputp).nloci][count-1].theta_fw;
					avgstatloci[0][(*inputp).nloci].faywuHo += stf;
				}
			}
            if((stf = (double)E_zeng((*inputp).nsam,(double)neutpar[0].thetaL,(double)neutpar[0].S/coef[0][0],(double)neutpar[0].S,coef[0])) != (double)-10000) {
                avgstatloci[0][(*inputp).nloci].zengE += stf;
            }
			if(neutpar[0].S > 1) {
			   if(((double)ZnA_(0,segsites,inputp)) != (double)-10000) {
					stf = (double)ZnA_(0,segsites,inputp)/((double)neutpar[0].S - (double)1.0);
					avgstatloci[0][(*inputp).nloci].rZA += stf;
				}
				if(((double)neutpar[0].B1) != (double)-10000) {
					stf = (double)neutpar[0].B1/((double)neutpar[0].S - (double)1.0);
					avgstatloci[0][(*inputp).nloci].wB += stf;
				}
				if(((double)neutpar[0].Q1) != (double)-10000) {
					stf = (double)neutpar[0].Q1/((double)neutpar[0].S);
					avgstatloci[0][(*inputp).nloci].wQ += stf;
				}
			}
			if((stf = (double)R2(neutpar[0].unic,neutpar[0].k,(*inputp).nsam,(int)neutpar[0].S)) != (double)-10000) {
				avgstatloci[0][(*inputp).nloci].R2 += stf;
			}
		}
        /*printing a point every 2% of the total iterations*/
		if((double)p+restp >= counterp) {
			restp += (double)p - counterp; 
			if((double)restp/(double)counterp > (double)1) {
				for(x=0;x<(int)floor(restp/counterp);x++) printf(".");
				restp -= (double)floor(restp/counterp) * counterp;
			}
			else printf(".");
			fflush(stdout);
			p = 1;
		}
		else p += 1;
    }
    /* alliberar les matrius i vectors */
    free(posit);
    for(i=0;i<(*inputp).nsam;i++)
        free(list[i]);
    free(list);
    
    free(neutpar[0].freq);
    free(neutpar[0].fhapl);
    free(neutpar[0].unic);
    free(neutpar);
    
    free(coef[0]);
    free(coef);
    
    return 0;
}

int mod_mhits(unsigned long segsites, struct var2 **inputp)
{
    long int x,y,z,nsit;
    int r,h,i,j,k;
    char a[1];
    int Sout;
    double rr,ratio,r_transv,r_transc;
    double ran1(void);
    int *mhsout;
	double poisso(double);
    long int segsitesin = 0; /*here is fixed to 0. That was a variable in case use S instead theta*/
    
	ratio = (*inputp)->ratio_sv;
    r_transc = ratio/(ratio + 1.);
    r_transv = (/*ratio + */0.5)/(ratio + 1.); /*suma de transc + 1nt a transversio, quan sumen els 2nt es total = 1*/
    Sout = (*inputp)->Sout;
    nsit = (*inputp)->nsites;
    x = 0;
        
    if(segsites > 0) {
		while(x<(int)segsites-1) {
			k = 0;
			while(posit[x] == posit[x+1]) {	/* buscar els mhits */
				k++;
				x++;
				if(x == (long int)segsites-1) break;
			}
			if(k) {				/* ordenar mhits de mes antic a mes recent */
				for(y=(x-k);y<x;y++) {
					j = 0;
					for(i=0;i<(*inputp)->nsam;i++) if(list[i][y] != '0') j++;
					for(z=y+1;z<(x+1);z++) {
						h = 0;
						for(i=0;i<(*inputp)->nsam;i++) if(list[i][z] != '0') h++;
						if(j < h) {
							for(i=0;i<(*inputp)->nsam;i++) {
								*a = list[i][y];
								list[i][y] = list[i][z];
								list[i][z] = *a;
							}
							j = h;
						}
					}
				}
				if(segsitesin == 0) {
					for(y=(x-k+1);y<x+1;y++) {	/* afegir les mutacions a la posicio mes antiga */
						r = -1;
						for(i=0;i<(*inputp)->nsam;i++) {
							if(list[i][y] != '0') {
								if(r == -1) {
									rr = (double)ran1();	/* inclou ratio trans/transv */
									if(rr < r_transc) r = 0;	/* transicio, la resta transversions */
									else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = 1;
										else r = 2;
									if(list[i][x-k] == '0') {
										if(r==0) list[i][x-k] = *a = '1';
										else if(r==1) list[i][x-k] = *a = '2';
											else if(r==2) list[i][x-k] = *a = '3';
									}
									else
										if(list[i][x-k] == '1') {
											if(r==0) list[i][x-k] = *a = '0';
											else if(r==1) list[i][x-k] = *a = '2';
												else if(r==2) list[i][x-k] = *a = '3';
										}
										else
											if(list[i][x-k] == '2') {
												if(r==0) list[i][x-k] = *a = '3';
												else if(r==1) list[i][x-k] = *a = '0';
													else if(r==2) list[i][x-k] = *a = '1';
											}
											else
												if(list[i][x-k] == '3') {
													if(r==0) list[i][x-k] = *a = '2';
													else if(r==1) list[i][x-k] = *a = '1';
														else if(r==2) list[i][x-k] = *a = '0';
												}
								}
								else list[i][x-k] = *a;
							}
						}
					}
				}
				else {
					for(y=(x+1-k),r=2;y<x+1;y++,r++) {	/* mutacions fix, totes s'han de veure (fins a tres mutacions) */
						for(i=0;i<(*inputp)->nsam;i++) 
							if(list[i][y] != '0') 
								list[i][x-k] = r + '0';
					}
					if(r>4) {
						perror("Error: more than 3 mutations in one position");
						exit(1);
					}
				}
			}
			/* incorporacio per mhits en especiacio.
			 Els mhits al outgroup no s'observen, nomes canviem el nt ancestral, pero no estan indicades!*/
			/*nombre de mutacions de Sout a la pos x-k*/
			/*
			 if((p = (double)poisso((double)Sout/(double)nsit)) > (double)0) {
				Sout -= (int)p;
				nsit--;
				
				rr = (double)ran1();
				if(rr < r_transc) r = 1;
				else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = 2;
					else r = 3;
				for(i=0;i<(*inputp)->nsam;i++) {
					if(list[i][x-k] == '0') list[i][x-k] = (char)r + '0';
					else if (list[i][x-k] == (char)r + '0') list[i][x-k] = '0';
				}
			}
			 */
			/*in case no outgroup but mhits: assign '0' randomly to the mhits positions (no '0' is observed but in putative outgroup)*/
			if((*inputp)->T_out == 0.) {
				r = 0;
				for(i=0;i<(*inputp)->nsam;i++) if(list[i][x-k] == '0') r++;
				if(r==0) {
					r = list[0][x-k];
					for(i=0;i<(*inputp)->nsam;i++) if(list[i][x-k] == r) list[i][x-k] = '0';
				}
			}
			x++;
		}
	}
    /* Treure els mhits de la branca de l'outgroup (no es mostra l'outgroup i per tant es veuen menys). */
    /* Per Sout mutacions en (nsites - segsites) */
    if(Sout) {
        /*
		nsit = (*inputp)->nsites-segsites;
        if(nsit <= 0) nsit = 0;
		if(nsit) {
			if(!(mhsout = (int *)calloc((unsigned long)(nsit),sizeof(int)))) {
				puts("calloc error in mhits.0b");
				return 0;
			}
			for(x=0;x<Sout;x++) {
				y = (long int)((double)floor((double)ran1()*(double)nsit));
				rr = (double)ran1();
				if(rr < r_transc) r = 1;
				else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = 2;
					else r = 3 ;
				if(mhsout[y] == 0) mhsout[y] = r;
				else if(mhsout[y] == 1) {
						if(r==1) mhsout[y] = 0;
						else mhsout[y] = r;
					}
					else if(mhsout[y] == 2) {
							if(r==1) mhsout[y] = 3;
							else mhsout[y] = r-2;
						}else if(mhsout[y] == 3) {
								if(r==1) mhsout[y] = 2;
								else mhsout[y] = r-2;
							}
			}
			Sout = 0;
			for(x=0;x<nsit;x++) if(mhsout[x] > 0) Sout++;
			free(mhsout);
		}
		*/
        if(!(mhsout = (int *)calloc((long int)((*inputp)->nsites),sizeof(int)))) {
            perror("calloc error in mhits.0b");
            exit(1);
        }
        for(x=0;x<Sout;x++) {
			y = (long int)floor((double)ran1()*(double)((*inputp)->nsites));
            
			rr = (double)ran1();
            if(rr < r_transc) r = 1;	/* transicio, la resta transversions */
            else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = 2;
			else r = 3 ;
			
			if(mhsout[y] == 0) 
				mhsout[y] = r;
            else {
				if(mhsout[y] == 1) {
					if(r==1) mhsout[y] = 0;
					else mhsout[y] = r;
				}
				else {
					if(mhsout[y] == 2) {
						if(r==1) mhsout[y] = 3;
						else mhsout[y] = r-2;
					}
					else {
						if(mhsout[y] == 3) {
							if(r==1) mhsout[y] = 2;
							else mhsout[y] = r-2;
						}
					}
				}
			}
			/*mutation outgroup-ingroup: Els mhits al outgroup no s'observen, nomes canviem el nt ancestral, pero no estan indicades!*/
			for(k=0;k<(int)segsites;k++) if(posit[k] == y) break;
			if(k<(int)segsites) {/*nombre de mutacions de Sout a la pos k*/
				mhsout[y] = 0;
				
				for(i=0;i<(*inputp)->nsam;i++) {
					if(list[i][k] == '0') list[i][k] = (char)r + '0';
					else if (list[i][k] == (char)r + '0') list[i][k] = '0';
				}
			}
        }
        Sout = 0;
        for(x=0;x<(long int)(*inputp)->nsites;x++) if(mhsout[x] > 0) Sout++;
        free(mhsout);
    }
    (*inputp)->Sout = Sout;

	return 1;
}
char **cmatrix(int nsam,unsigned long len)	/* defineix l'espai per col.locar els polimorfismes */
{
    int i;
    char **m;
    
    if(!(m=(char **)malloc((unsigned)nsam*sizeof(char *)))) {
        puts("alloc error in cmatrix");
        return 0;
    }
    for(i=0;i<nsam;i++)
        if(!(m[i] = (char *)malloc((unsigned long)len*sizeof(char)))) {
            puts("alloc error in cmatrix.2");
            return 0;
        }
    return(m);
}

unsigned long gensam(long int npop,int nsam,int inconfig[],unsigned long nsites,double theta/*,unsigned long segsites*/,
	double r,double f,double track_len,double mig_rate,int mhits/*,unsigned long iteration*/, double *factor/*, double *lengtht*/,
	int ifselection, double pop_sel, double sinit,double pop_size,long int sel_nt,double T_out, int *Sout,int nintn,double *nrec,
	double *npast,double *tpast,int split_pop, double time_split, double time_scoal, double factor_anc, double *freq, 
	double tlimit,int iflogistic,double ts, double factor_chrn, double rsv)
{
	int i,ii;
    unsigned long nsegs,seg,ns,start,end,len,segsit,k; 
    struct segl *seglst;
    double /*nsinv,*/tseg,tt;
    double tout;
	double r_transc,r_transv;
    struct segl *segtre_mig(long int ,int,int *,unsigned long,double,double,double,
        double ,unsigned long *,/*unsigned long,*/
        double *,int,double,double,double,long int,int *,int *,int,double *,double *,double *,
        int, double, double, double, double *,double,int,double,double);/* used to be: [MAXSEG]; */
    double ttime(struct node *, int);
    double poisso(double), ran1(void);
    int biggerlist(int);
    void make_gametes(int,struct node *,double,unsigned long,unsigned long,int,double,double);
	void locate(unsigned long,/*double*/unsigned long,/*double*/unsigned long,/*double **/unsigned long *,int /*nou*/);
    /*partial selection*/
    int all_sel,*selnsam; /*number of lines under selection and the vector with the lines (the first all_sel lines)*/
    unsigned long segsit_sel=0;/*parameter for partial selection*/    
    /*void make_gametes_psel(int,unsigned long,int,int *);*/
    void locate_psel(unsigned long,unsigned long,unsigned long,unsigned long *,int,unsigned long);
    /*nsinv = 1./nsites;*/
    if(!(selnsam = (int *)malloc((nsam)*sizeof(int)))) {
		puts("malloc error. gensam.");
		return 100000;
	}
	
    all_sel = 0;
    seglst = segtre_mig(npop,nsam,inconfig,nsites,r,f,track_len,mig_rate,&nsegs,
        /*iteration,*/factor,ifselection,pop_sel,sinit,pop_size,sel_nt,&all_sel,selnsam,
        nintn, nrec, npast, tpast,split_pop,time_split,time_scoal,factor_anc,freq,tlimit,iflogistic,ts,factor_chrn);
    
    r_transc = rsv/(rsv + 1.);
    r_transv = ((double)0.5)/(rsv + 1.); /*suma de transc + 1nt a transversio, quan sumen els 2nt es total = 1*/

	ns = 0;
	*Sout = 0;
    if(theta >= 0.0) {
        for(seg=0,k=0;k<nsegs;seg=seglst[seg].next,k++) {
            end = (k<nsegs-1 ? (unsigned long)seglst[seglst[seg].next].beg -1 : nsites-1); /*next és l'index, beg és el punt físic */
            start = seglst[seg].beg;
            len = end - start + 1;
            tseg = len*(theta/nsites); 		/* part de la theta que li correspon al segment */
            tt = (double)ttime(seglst[seg].ptree,nsam); /* Ttot pel segment, en funció de 4No respecte la pob 0*/
            if(mhits) {/*T_out es el temps a la pop ancestral*/
				if(T_out <= 0.) 
					*Sout = 0;
				else {
					tout = 2.0 * T_out;/*2 branches. T_out is considered as a fixed value, not a parameter.*/
					tout += -1.0*(double)log((double)1-ran1());/*time after divergence, assuming equal No=1*/
					tout -= (seglst[seg].ptree + 2*nsam-2)->time; /*substract the distance of the sample*/
					if(tout < (double)0) tout = (double)0;	/*Outgroup t can not accumulate negative mutations*/
					*Sout += (int)poisso((double)(tseg*tout));	/* Sout needed to calculate hka and mhits */
				}
            }
            segsit = (unsigned long)poisso((double)(tseg*tt));		/* nombre de mutacions al segment */
            if(segsit == (unsigned long)0 && all_sel > (int)0 && all_sel < nsam && sel_nt >= (long int)start && sel_nt <= (long int)end)
				segsit = 1;/*we force the selective mut*/
            if(segsit > len && mhits == 0) segsit = len;	/* mutacions discretes ...: afegit.. */
            if((segsit + ns) >= maxsites) {	/* refem la matriu dels polimorfismes */
                maxsites = segsit + ns + SITESINC;
                posit = (/*double **/unsigned long *)realloc(posit,(unsigned long)maxsites*sizeof(/*double*/unsigned long)); 
                /* canvia mida del vector dels nombres dels polimorfismes */
                if(posit==NULL) {
					puts("realloc error. gensam.1");
					return 100000;
				}
                if(biggerlist(nsam) == 0) {
					puts("realloc error. gensam.1");
					return 100000;
				}	/* refem la llista dels polimorfismes */
            }
            /*partial selection*//*not well debugged yet*/ 
            if(all_sel > (int)0 && all_sel < (int)nsam && sel_nt >= (long int)start && sel_nt <= (long int)end) {
                /*make_gametes_psel(nsam,ns,all_sel,selnsam);*//*locate the sel mut in the selnsam lines*/
                segsit_sel = (unsigned long)1;
            }
            make_gametes(nsam,seglst[seg].ptree,tt,(unsigned long int)(segsit-segsit_sel),(unsigned long int)(ns+segsit_sel),mhits,r_transc,r_transv);	/* posa les mutacions  a list*/
            free(seglst[seg].ptree);
            /*partial selection*/
            if(segsit_sel == (unsigned long)1) {
                locate_psel(segsit,start,len,posit+ns,mhits,sel_nt-start);
				ii = ns;
				while((long int)posit[ii] != sel_nt) ii++; 
				for(i=0;i<nsam;i++) {
					list[i][ns] = list[i][ii];
					if(all_sel>i) list[i][ii] = '1';
					else list[i][ii] = '0';
				}
                segsit_sel = (unsigned long)0;
            }
            else locate(segsit,start,len,posit+ns,mhits);/* posa el nombre de les mutacions a la matriu */
            ns += segsit;
        }
    }
    free(selnsam);
    return(ns);
}
int biggerlist(int nsam)	/* fa més gran la matriu dels polimorfismes */
{
    int i;
    for(i=0;i<nsam;i++) {
        list[i] = (char *)realloc(list[i],maxsites*sizeof(char));
        if(list[i] == NULL) {
            puts("realloc error. biggerlist");
            return 0;
        }
    }
    return 1;
}
double ttime(struct node *ptree, int nsam)	/* la Ttot de l'arbre */
{
    double t;
    int i;
    
    t = (ptree + 2*nsam-2)->time;
    for(i=nsam;i<2*nsam-1;i++) t += (ptree + i)->time;
    return(t);
}
void make_gametes_psel(int nsam,unsigned long ns,int all_sel,int *selnsam)
{	/* posa la mutacio sel (partial) a la matriu list */
    int tip;
    for(tip=0;tip<nsam;tip++)
        list[tip][ns] = '0';
    for(tip=0;tip<all_sel;tip++) 
        list[selnsam[tip]][ns] = '1'; 
}
void make_gametes(int nsam, struct node *ptree, double tt,unsigned long newsites,unsigned long ns, int mhits, double r_transc,double r_transv)
{					/* posa les mutacions a la matriu list */
    unsigned long j;
    int tip,node;
    int pickb(int, struct node *,double);
    int tdesn(struct node *,int,int);
       
	double rr;
	double ran1();
	char r;
	
    for(j=ns;j<ns+newsites;j++) {
        node = pickb(nsam,ptree,tt);	/* busca una branca al'atzar en funcio de la mida de t*/
		if(mhits) {
			rr = (double)ran1();
			if(rr < r_transc) r = '1';/* transicio, la resta transversions */
			else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = '2';
				else r = '3';
		}
		else r='1';
        for(tip=0;tip<nsam;tip++) {
            if(tdesn(ptree,tip,node))	/* posa mutació si la mostra té relació amb la branca */
                list[tip][j] = r; 	
            else
                list[tip][j] = '0';
        }
    }
}
int pickb(int nsam, struct node *ptree,double tt) /* agafa la branca a on ha caigut la mutació */
{
    double x,y,z;
    int i;
    double ran1(void);
    
    x = (double)ran1()*(double)tt;
    for(i=0,y=0;i<2*nsam-2;i++) {
        z = (ptree + (ptree+i)->abv)->time - (ptree+i)->time;
        if(z < 0.) z = 0.; /* en cas d'extincio, per si hi han problemes de precissio */
        y += z;
        if(y >= x) return(i);
    }
    return(i);
}
int tdesn(struct node *ptree, int tip, int node) /* mira si la mostra que mirem esta relacionada amb la branca */
{
    int k;
    
    for(k=tip;k<node;k = (ptree+k)->abv);
    
    if(k==node) return(1);
    else return(0);
}
void mnmial(unsigned long n,unsigned long nclass,double *p,unsigned long *rv)
{   
    double x,s;
    long int i,j;
    double ran1(void);
    
    for(i=0;i<(long int)nclass;i++) rv[i]=0;	/* inicialitzar */
    for(i=0;i<(long int)n;i++) {	/* posa les n mutacions */
        x = (double)ran1();	/* agafa un nombre [0,1) */
        j = 0;
        s = p[0];	/* p és el temps total de cada segment relatiu a 1 */
        while((x>s) && (j<(long int)(nclass-1))) s += p[++j]; /* busca el segment fins que ran sigui major que s, posa a j+1 */
        rv[j]++; 	/* afegeix la mutació al segment j */
    }
}
void mnmial2(unsigned long n,unsigned long nclass,double *p,unsigned long *rv,unsigned long *len)
{
    double x,s;
    long int i,j;
    double ran1(void);
    
    for(i=0;i<(long int)nclass;i++) rv[i]=0;	/* inicialitzar */
    i = 0;
    while(i<(long int)n) {	/* posa les n mutacions */
        x = (double)ran1();	/* agafa un nombre [0,1) */
        j = 0;
        s = p[0];	/* p és el temps total de cada segment relatiu a 1 */
        while((x>s) && (j<(long int)(nclass-1))) s += p[++j]; /* busca el segment fins que ran sigui major que s, posa a j+1 */
        if(len[j] > 0) {
            len[j]--;	
            rv[j]++; 	/* afegeix la mutació al segment j */
            i++;
        }
    }
}
void locate_psel(unsigned long n,unsigned long beg,unsigned long len,unsigned long *ptr,int mhits,unsigned long sel_nt)
/* localitza les mutacions en un fragment */
{
     unsigned long i;
     void ordran_psel(unsigned long,unsigned long *,unsigned long,int,unsigned long);
     
     ordran_psel(n,ptr,len,mhits,sel_nt);	/* mutacions en [0,len) ordenades de major a menor */
     for(i=0;i<n;i++) ptr[i] = beg + ptr[i];	/* les escala entre inici i final del fragment */
}

void ordran_psel(unsigned long n,unsigned long *pbuf,unsigned long len,int mhits,unsigned long sel_nt)
{
    void ranvec_psel(unsigned long,unsigned long *,unsigned long,int,unsigned long);
    void order(unsigned long,unsigned long *);

    ranvec_psel(n,pbuf,len,mhits,sel_nt);
    order(n,pbuf);
}
void ranvec_psel(unsigned long n,unsigned long *pbuf,unsigned long len,int mhits,unsigned long sel_nt)
 /* posa un nombre entre [0,len) */
{
    long int i,x;
    double ran1(void);
    
    pbuf[0] = sel_nt;
    for(i=1;i<(long int)n;i++) {
        while((pbuf[i] = (unsigned long)((double)ran1()*(double)len)) == pbuf[0]);   
        if(!mhits) {	/* bucle per no mhits (mhits=0) */
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    pbuf[i] = (unsigned long)((double)ran1()*(double)len);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}

void locate(unsigned long n,/*double*/unsigned long beg,/*double*/unsigned long len,/*double*/unsigned long *ptr,int mhits/*nou*/)
/* localitza les mutacions en un fragment */
{
     unsigned long i;
     void ordran(unsigned long,/*double **/unsigned long *,/*nou->*/unsigned long,int /*nou*/);
     
     ordran(n,ptr,len,mhits);	/* mutacions en [0,len) ordenades de major a menor */
     for(i=0;i<n;i++) ptr[i] = beg + ptr[i]/**len*/;/* les escala entre inici i final del fragment */
}

void ordran(unsigned long n, /*double */unsigned long *pbuf,/*nou->*/unsigned long len,int mhits/*nou*/)
{
    void ranvec(unsigned long,/*double **/unsigned long *,/*nou->*/unsigned long,int/*nou*/);
    /* posa un nombre entre [0,len) */
    void order(unsigned long,/*double **/unsigned long *);/* ordena els valors */

    ranvec(n,pbuf,len,mhits);
    order(n,pbuf);
}

void ranvec(unsigned long n, /*double */unsigned long *pbuf,/*nou->*/unsigned long len,int mhits /*nou*/) /* posa un nombre entre [0,len) */
{
    long int i,x;
    double ran1(void);
    
    
    for(i=0;i<(long int)n;i++) {
        pbuf[i] = (unsigned long)((double)ran1()*(double)len);        
        if(!mhits) {	/* bucle per no mhits (mhits=0) */
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    pbuf[i] = (unsigned long)((double)ran1()*(double)len);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}

void order(unsigned long n, /*double */unsigned long *pbuf)/* ordena els valors */
{
    long int gap,i,j;
    /*double*/unsigned long temp;
    
    for(gap= n/2; gap>0;gap /= 2)
        for(i=gap;i<(long int)n;i++)
            for(j=i-gap;j>=0 && pbuf[j]>pbuf[j+gap];j -= gap) {
                temp = pbuf[j];
                pbuf[j] = pbuf[j+gap];
                pbuf[j+gap] = temp;
            }
}

double rest_branch(int nsam,struct node *ptree, double tt)/*agafa una branca, en cas de l'arrel agafa les dues branques, perque son complementaries*/
{
    double x,y,z;    
    int i,j;
    double ran1(void);
    
    x= (double)ran1()*tt;
    for(i=0,y=0;i<2*nsam-2;i++) {
        z = (ptree + (ptree+i)->abv)->time - (ptree+i)->time;
        if(z < 0.) z = 0.; /* en cas d'extincio, si hi ha problemes de precissio */
        y += z;
        if(y >= x) { /* ja ha arribat a la branca */
            if((ptree+i)->abv == 2*nsam-1) { /* si toca l'ultim node agafem les dues branques*/
                for(j=0;j<2*nsam-2;j++) {
                    if((ptree+j)->abv == 2*nsam-1 && j!=i) {
                        z = ((ptree + (ptree+i)->abv)->time - (ptree+i)->time) + ((ptree + (ptree+j)->abv)->time - (ptree+j)->time);
                        if(z < 0.) z = 0.;
                        return(z);
                    }
                }
            }
            else {
                z = (ptree + (ptree+i)->abv)->time - (ptree+i)->time;
                if(z < 0.) z = 0.;
                return(z);
            }
        }
    }
    z = (ptree + (ptree+i)->abv)->time - (ptree+i)->time;
    if(z < 0.) z = 0.;
    return(z);
}

int calc_neutpar(double valuer,int valuep,unsigned long segsit,struct var2 *inputp, struct dnapar *ntpar)
{
    unsigned long pi;
    double k_;
    unsigned long S;
    int nhapl;
    int B;
    int A;
    int *haplotype = 0;
    long int segsitesm1;
    
    int inits,nsam;
    int val10,val20,val21;
    unsigned long j,k;
    int a,b,c,d,h,i,nmh,x;
	unsigned long int comb;
    char *hapl;
    int ispolnomhit(unsigned long,int,int,int);
	int **veca;

	double k2_,hij,Shij,Shj2;
	double C=0.;
	/*double x1,x2,xh2,xacc;*/
	int zbracC(double *,double *,double,int);
	double zriddrC(double,double,double,double,int);
	int Min_rec(int,int,int,int,int);
    
    if(valuep == 0) {
        inits = 0;
        nsam  = (*inputp).nsam;
    }
    else {
        if(valuep == 1) {
            inits = 0;
            nsam  = (*inputp).config[0];
        }
        else {
            inits = (*inputp).config[0];
            nsam  = (*inputp).config[1];
        } 
    }
    comb = (unsigned long int)(((double)nsam*((double)nsam-(double)1.0))/(double)2.0);
    
   if(segsit == 0 || nsam < 2) {
		if(segsit == 0) {
			(ntpar)->B1 = 0;
			(ntpar)->Q1 = 0;
			(ntpar)->k = (double)0;
			(ntpar)->S = 0;
			(ntpar)->nhapl = 1;
			(ntpar)->C = (double)0;
			for(a=0;a<2;a++) {
				(ntpar)->freq[a] = 0;
				(ntpar)->unic[a] = 0;
				(ntpar)->fhapl[a] = 0;
			}
			(ntpar)->fhapl[0] = nsam;
			(ntpar)->Rm = 0;
			(ntpar)->thetaL = (double)0;
			(ntpar)->mhsites = 0;
		}
	   if(nsam < 2) {
		   (ntpar)->B1 = -10000;
		   (ntpar)->Q1 = -10000;
		   (ntpar)->k = (double)-10000;
		   (ntpar)->S = -10000;
		   (ntpar)->nhapl = -10000;
		   (ntpar)->C = (double)0;
		   for(a=0;a<2;a++) {
			   (ntpar)->freq[a] = -10000;
			   (ntpar)->unic[a] = -10000;
			   (ntpar)->fhapl[a] = -10000;
		   }
		   (ntpar)->fhapl[0] = nsam;
		   (ntpar)->Rm = 0;
		   (ntpar)->thetaL = (double)-10000;
		   (ntpar)->mhsites = 0;
	   }
    }
    else {
		if((veca = (int **)calloc(segsit,sizeof(int *))) == 0) {
			puts("calloc error veca.1");
			return -1;
		}
		for(d=0;d<(long int)segsit;d++) {
			if((veca[d] = (int *)calloc((*inputp).nsam,sizeof(int))) == 0) {
				puts("calloc error veca.2");
				return -1;
			}
		}
    
        /* calcul B' i Q, pero sense dividir per S (Wall) */
        B = 0;
        A = 0;

        a = 0;
        if((segsitesm1 = segsit-1) < 0) segsitesm1 = 0;
        for(j=0;j<(unsigned long)segsitesm1;) {
            k = j;
            while(k+1 < segsit) { /*calcular k*/
                if(ispolnomhit(k,inits,nsam,(*inputp).nsam) > 0) break;
                else {
					k++;
				}
            }
            j = k+1;
            while(j < segsit) { /*calcular j*/
                if(ispolnomhit(j,inits,nsam,(*inputp).nsam) > 0) break;
                else j++;
            }
            if(j < segsit) {                
                val20 = val21 = -1;
                b = 0;
                for(i=inits;i<inits+nsam;i++) {
                    val10 = (list[i][k] - 48)*4 + (list[i][j]- 48);
                    if(val20 == -1) val20 = val10;
                    else{
                        if(val20 != val10) {
                            if(val21 == -1) val21 = val10;
                            else  if(val21 != val10) {
                                b = 1;/*incongruent*/
                                break;
                            }
                        }
                    }
                }
				if(!b) {/*congruent*/
					B += 1;
					if(!A) {
						for(i=inits;i<inits+nsam;i++) {
							if(list[i][j] > '0') x = '1';
							else x = '0';
							veca[A][i] = x;
						}
						A += 1;
					}
					else {
						a = 0;
						for(c=0;c<A;c++) {
							d = 0;
							for(i=inits;i<inits+nsam;i++) {
								if(list[i][j] > '0') x = '1';
								else x = '0';
								if(veca[c][i] == x) {
									d += 1;
								}
							}
							if(d == nsam || d == 0) {
								a = 1;
								break;
							}
						}
						if(!a) {
							for(i=inits;i<inits+nsam;i++) {
								if(list[i][j] > '0') x = '1';
								else x = '0';
								veca[A][i] = x;
							}
							A += 1;
						}
					}
				}
			}
		}
		
        (ntpar)->Q1 = B+A;
        (ntpar)->B1 = B;
        
        /* calcul de S, pi, fr, nhapl i unics/seq (excluding mhits)*/
        if((hapl = (char *)calloc(nsam*segsit,sizeof(char))) == NULL) {
            puts("calloc error calc_neutpar.0");
            return 0;
        }
        k_ = (double)0;
		k2_= (double)0;
        for(a=0;a<nsam;a++) {
            (ntpar)->freq[a] = 0;
            (ntpar)->unic[a] = 0;
        }
        
		hij = Shij = (double)0;
		Shj2 = (double)0;            
		
		S = 0;
        nhapl = 0;
		nmh = 0;
        for(j=0;j<segsit;j++) {
            pi = 0;
            while(j < segsit) {
                if((h=ispolnomhit(j,inits,nsam,(*inputp).nsam)) > 0) break; /*h is the frequency of the new mutation*/
                else {
					if(h == -2) nmh++;
					j++;
				}
            }            
            if(j<segsit) {
                (ntpar)->freq[h] += 1;
                for(a=inits;a<inits+nsam;a++) {
                    hapl[(a-inits)*segsit+S] = list[a][j];
                    /*new two lines: for singleton mutations (no outgroup)*/
                    if(h == 1) if(list[a][j] != '0') (ntpar)->unic[a-inits] += 1;
                    if(h == nsam-1) if(list[a][j] == '0') (ntpar)->unic[a-inits] += 1;
                }
                S++;
                for(a=inits;a<inits+nsam-1;a++)
                    for(b=a+1;b<inits+nsam;b++)
                        if(list[a][j] != list[b][j]) pi++;
                k_ += (double)pi;
#if CALCREC == 1
				hij = (double)1 - ((double)h*(double)h + (double)(nsam-h)*(double)(nsam-h))/((double)nsam*(double)nsam);
				Shij += hij; 
				Shj2 += (hij * hij);
#endif
            }
        }
        (ntpar)->k = k_/(double)comb;
        (ntpar)->S = S;
		(ntpar)->mhsites = nmh;

        /*calcular thetaL*/
		(ntpar)->thetaL = 0;
		for(a=0;a<nsam;a++) {
            (ntpar)->thetaL += (double)a * (double)ntpar->freq[a];
		}
		(ntpar)->thetaL *= (double)1/(double)(nsam-1);

#if CALCREC == 1
		/*recombination*/
		k_ = (double)2*k_/((double)nsam*(double)nsam);/*km*/
		for(a=inits;a<inits+nsam-1;a++) {
			for(b=a+1;b<inits+nsam;b++) {
				pi = (unsigned long)0;
				for(j=0;j<segsit;j++) {
					while(j < segsit) {
						if((h=ispolnomhit(j,inits,nsam,(*inputp).nsam)) > 0) break; 
						else j++;
					}            
					if(j<segsit) if(list[a][j] != list[b][j]) pi++;
				}
				k2_+= ((double)pi - k_) * ((double)pi - k_);
			}
        }
		Sk2 = ((double)2*k2_ + (double)nsam*k_*k_)/((double)nsam * (double)nsam);
		thetam = Shij * ((double)nsam / ((double)nsam - (double)1));
		gcn = (Sk2 - Shij + Shj2)/(thetam*thetam);	
		/*calculate the range*/
		x1 = (double)0.00001;
		x2 = xh2 = (double)segsit;
		while(zbracC(&x1,&x2,gcn,nsam) == 0) {
			if(x2 > 1e10) {
				C = 1e10;
				return 1;
			}
			x2 = xh2 * (double)10;
		}
		/*if(zbracC(&x1,&x2,gcn,nsam) == 0) {*//*if 0, that means that funtionC is not crossing 0*/
		/*	C = (double)-10000;*//*na*/
		/*}*/
		/*estimate the value of C*/
		if(x2 < 1e10) {
			xacc = (double)1e-6*x2; /* accuracy of five numbers*/
			C = (double)zriddrC(x1,x2,xacc,gcn,nsam);
			if(C == -1.11e30) C = (double)-10000; /*na*/
		}
#else
		C = -10000.;
#endif
		(ntpar)->C = C;
	
		/*haplotypes*/
        if((haplotype = (int *)malloc(nsam*sizeof(int))) == NULL) {
            puts("calloc error calc_neutpar.1");
            return 0;
        }
        (ntpar)->nhapl = 0;
        for(a=0;a<nsam;a++) haplotype[a] = 1;            
        for(a=0;a<nsam-1;a++) {
            if(haplotype[a]) {
                (ntpar)->nhapl += 1;
                for(b=a+1;b<nsam;b++) {
                    if(haplotype[b]) {
                        if(memcmp(hapl+a*segsit,hapl+b*segsit,S) == 0) { 
                            haplotype[a] += 1;
                            haplotype[b] = 0;
                        }
                    }
                }
            }
        }
        if(haplotype[a]) (ntpar)->nhapl += 1;
        
        for(a=0;a<nsam;a++) (ntpar)->fhapl[a] = haplotype[a]; /*calcular freq. haplotips*/
/*
        b = 0;
        for(a=0;a<nsam;a++) if(haplotype[b] < haplotype[a]) b=a;
        printf("%d\r",haplotype[b]);
*/
		if(valuer || (*inputp).mhits) (ntpar)->Rm = Min_rec(0,segsit,nsam,inits,(*inputp).nsam);
		else (ntpar)->Rm = (int)0;

        free(hapl);
        free(haplotype);

		for(d=0;d<(long int)segsit;d++) free(veca[d]);
		free(veca);
    }
    return 1;
}

double ZnA_(int valuep,unsigned long segsites,struct var2 *inputp)
{
    double ZnA;
    long k,j,comb;
    int i,inits,nsam;
    int vala,valb,val00,a,b,na,nb;
    double A,B,C;
    int ispolnomhit(unsigned long,int,int,int);
    
    if(valuep == 0) {
        inits = 0;
        nsam  = (*inputp).nsam;
    }
    else {
        if(valuep == 1) {
            inits = 0;
            nsam  = (*inputp).config[0];
        }
        else {
            inits = (*inputp).config[0];
            nsam  = (*inputp).config[1];
        } 
    }
    
    ZnA = (double)0;
    j = 0;
    comb = 0;
    while(j+1 < (long int)segsites) {
        while(j < (long int)segsites) {
            if(ispolnomhit(j,inits,nsam,(*inputp).nsam) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < (long int)segsites) {
            if(ispolnomhit(k,inits,nsam,(*inputp).nsam) > 0) break;
            else k++;
        }
        if(k < (long int)segsites) {
            /*calcular freqs p1,q1*/
            vala = valb = 1;
            a = list[inits][j];
            b = list[inits][k];
            for(i=1+inits;i<inits+nsam;i++) {
                if(list[i][j] == a) vala++;
                else na = list[i][j];
                if(list[i][k] == b) valb++;                                
                else nb = list[i][k];
            }
            if(nsam - vala > vala) {
                a = na;
                vala = nsam - vala;
            }
            if(nsam - valb > valb) {
                b = nb;
                valb = nsam - valb;
            }
            /*calcular p1q1*/
            val00 = 0;
            for(i=0+inits;i<inits+nsam;i++) if(list[i][j] == a && list[i][k] == b) val00++;
            
            /*calcular R2 = (p1q1 - (p1*q1))"2 / (p1*(1-p1) * (q1*(1-q1))) */
            A = (double)vala/(double)nsam;
            B = (double)valb/(double)nsam;
            C = (double)val00/(double)nsam;
            /*if(vala != ((double)nsam-1.0) && valb != ((double)nsam-1.0)) if only informative */
            ZnA += ((C - (A*B)) * (C - (A*B))) / (A*((double)1.0 - A)*B*((double)1.0 - B));
            comb++;
        }
        j++;
    }
    /*
    if(comb) ZnA = ZnA/(double)comb;
    else return((double) -10000);
    */
    return(ZnA);
}

int ispolnomhit(unsigned long j,int init, int nsam,int totalsam)
{
    int i,h,g;
    int a0,a1;
    
    if(j) if(posit[j-1] == posit[j]) return(-1);/*estem a la segona mutació o més*/

    h = g = a1 = 0;
	a0 = '0';
    for(i=0;i<totalsam;i++) {
        if((int)list[i][j] != a0) {
            if(!a1) a1 = (int)list[i][j];
            if((int)list[i][j] != a1) return(-2);/*mhit, including outgroup*/
        }
    }
 	
    h = g = a1 = 0;
	a0 = '0';
    for(i=init;i<init+nsam;i++) {
        if((int)list[i][j] == a0) h++;
        else {
            if(!a1) a1 = (int)list[i][j];
            if((int)list[i][j] == a1) g++;
            else return(-2);/*en el cas de la primera mutacio hagin almenys tres variants*/
        }
    }
    if(h==nsam || h == 0) return(0);
    return(nsam-h);
}

double koutgJC(int nsam, int Sout, int *freq, unsigned long nsites)
{
    /**/int x;/**/
    double div;
    
    if(nsam) {
        div = (double)Sout;
        /**/
		for(x=1;x<nsam;x++) 
            div += (double)freq[x]*(double)x/(double)nsam;
        /**/
		div/= (double)nsites;
        if(div >= (double).75) return (double) -10000;
        else {
            div = (double)-.75*(double)log((double)1-(double)4/(double)3*div);/*Jukes and Cantor correction...*/
            return div*(double)nsites;
        }
    }
    else return (double) -10000;
}
double koutg(int nsam, int Sout, int *freq/*, unsigned long nsites*/)
{
    /**/int x;/*d*/
    double div;
    
    if(nsam) {
        div = (double)Sout;
        /**/
		for(x=1;x<nsam;x++) 
            div += (double)freq[x]*(double)x/(double)nsam;
        /**/
		return div;
    }
    else return (double) -10000;
}

/*Wall's program for calculating minimum recombination events*/
int Min_rec(int x, int segsit, int nsam, int inits,int totalsam)
{  /* Calculate min # rec. events */
  int a, b, c, e, gtest, flag = 0;
  int h;
  int t11,t12,t21,t22;
  int ispolnomhit(unsigned long,int,int,int);
  
	if (segsit<2 || x >= (segsit-1)) return (0);
	
	for (a=x+1; a<segsit; ++a) {
		while(a < segsit) {
			if((h=ispolnomhit((unsigned long)a,inits,nsam,totalsam)) > 0) break; /*h is the frequency of the new mutation*/
			else a++;
		}            
		if(a < segsit) {
			for (b=x; b<a; ++b) {
				while(b < a) {
					if((h=ispolnomhit((unsigned long)b,inits,nsam,totalsam)) > 0) break; /*h is the frequency of the new mutation*/
					else b++;
				}
				if(b < a) {
					t21 = list[0][b];
					t11 = list[0][a];
					for (e=inits+1; e<inits+nsam; ++e) {
						if(list[e][b] != t21) t22 = list[e][b];
						if(list[e][a] != t11) t12 = list[e][a];
					}
					
					gtest = 0;
					for (e=inits; e<inits+nsam; ++e) {
						if (list[e][b] == t21 && list[e][a] == t11) {
							++gtest;
						break;
						}
					}
					for (e=inits; e<inits+nsam; ++e) {
						if (list[e][b] == t21 && list[e][a] == t12) {
							++gtest;
							break;
						}
					}
					for (e=inits; e<inits+nsam; ++e) {
						if (list[e][b] == t22 && list[e][a] == t11) {
							++gtest;
							break;
						}
					}
					for (e=inits; e<inits+nsam; ++e) {
						if (list[e][b] == t22 && list[e][a] == t12) {
							++gtest;
							break;
						}
					}
					if (gtest == 4) {
						flag = 1;
						break;
					}
				}
			}
		}
		if (flag == 1) break;
	}
	if (a >= segsit) return (0);
	else {
		c = Min_rec(a,segsit,nsam,inits,totalsam);
		return (1+c);
	}
	return 0;
}


