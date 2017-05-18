/*
 *  calcobsneuttests.c
 *  MuLoNeTests
 *
 *  Created by sebas on Wed Feb 26 2003.
 *
 */

#include "MuLoNeTests.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int calc_neutests(struct statistics **matrix,int *unic,int *outgroup, int *n_loci/*, FILE *file_output*/)
{
    double p[12];
    void init_coef(double *,int);
    double tajima_d(double, int, double *);/*not shared*/
    double Fs(int, double, int);/*not shared*/
    double fl_d(/*int,*/int,int, double *); /* outgroup */
    double fl_f(/*int, */int, int, double, double *); /* outgroup */
    double fl_d2(int,double,int, double *); /* not outgroup */
    double fl_f2(int,double, int, double, double *); /* not outgroup */
	double fay_wu_obs(int,double,double,double,double *,double);/*outgroup*/
    double R2(int *,double,int,int);/*not outgroup*/
	double E_zeng(int,double,double,double,double *); /*outgroup*/
	  
    /*be careful with divisions by zero: sample sizes and biallsites. The result is in that case -10000 */
    
    matrix[0][*n_loci].tajimaD = (double)0.;
    matrix[0][*n_loci].fuliD = (double)0.;
    matrix[0][*n_loci].fuliDn = (double)0.;
    matrix[0][*n_loci].fuliF = (double)0.;
    matrix[0][*n_loci].fuliFn = (double)0.;
    matrix[0][*n_loci].fuFs = (double)0.;
    matrix[0][*n_loci].faywuH = (double)0.;
    matrix[0][*n_loci].faywuHo = (double)0.;
    matrix[0][*n_loci].R2 = (double)0.;
    matrix[0][*n_loci].rZA = (double)0.;
    matrix[0][*n_loci].wB = (double)0.;
    matrix[0][*n_loci].wQ = (double)0.;
	matrix[0][*n_loci].zengE = (double)0.;
    matrix[0][*n_loci].hka = (double)-1.0;

    init_coef(p,matrix[0][*n_loci].nsamples);
    
    matrix[0][*n_loci].tajimaD = tajima_d(matrix[0][*n_loci].theta_taj,matrix[0][*n_loci].biallsites,p);/*not shared*/
    matrix[0][*n_loci].fuFs = (double)Fs(matrix[0][*n_loci].nsamples,
        (double)matrix[0][*n_loci].theta_taj,matrix[0][*n_loci].nhapl);/*not shared*/
    matrix[0][*n_loci].R2 =
        R2(unic,matrix[0][*n_loci].theta_taj,matrix[0][*n_loci].nsamples,matrix[0][*n_loci].biallsites);/*not shared*/
    
    if(*outgroup) {
        if(matrix[0][*n_loci].biallsites > 1) {
            matrix[0][*n_loci].rZA = matrix[0][*n_loci].za/(matrix[0][*n_loci].biallsites - (double)1.);/*not shared*/
            matrix[0][*n_loci].wB  = (double)matrix[0][*n_loci].b/(matrix[0][*n_loci].biallsites - (double)1.);/*not shared*/
            matrix[0][*n_loci].wQ  = (double)matrix[0][*n_loci].q/(double)matrix[0][*n_loci].biallsites;/*not shared*/
        }
        else {
            matrix[0][*n_loci].rZA = (double)-10000;
            matrix[0][*n_loci].wB  = (double)-10000;
            matrix[0][*n_loci].wQ  = (double)-10000;
        }
        matrix[0][*n_loci].fuliD  = 
            fl_d(/*matrix[0][*n_loci].nsamples,*/(int)matrix[0][*n_loci].theta_fuli,matrix[0][*n_loci].biallsites,p);
        matrix[0][*n_loci].fuliF  = 
            fl_f(/*matrix[0][*n_loci].nsamples,*/(int)matrix[0][*n_loci].theta_fuli,
                matrix[0][*n_loci].biallsites, matrix[0][*n_loci].theta_taj,p);
		matrix[0][*n_loci].faywuH = (double)fay_wu_obs((int)matrix[0][*n_loci].nsamples,(double)matrix[0][*n_loci].theta_L,(double)matrix[0][*n_loci].theta_wat,(double)matrix[0][*n_loci].biallsites,p,(double)matrix[0][*n_loci].theta_taj);
									/*double fay_wu_obs(int n,double thetaL,double thetaw,double S,double *coef,double pi)*/
									/*matrix[0][*n_loci].theta_taj - matrix[0][*n_loci].theta_fw;*/
		matrix[0][*n_loci].faywuHo = matrix[0][*n_loci].theta_taj - matrix[0][*n_loci].theta_fw;
		matrix[0][*n_loci].zengE = (double)E_zeng((int)matrix[0][*n_loci].nsamples,(double)matrix[0][*n_loci].theta_L,(double)matrix[0][*n_loci].theta_wat,(double)matrix[0][*n_loci].biallsites,p);
        matrix[0][*n_loci].fuliDn  = fl_d2(matrix[0][*n_loci].nsamples,
            ((matrix[0][*n_loci].theta_fulin)*((double)matrix[0][*n_loci].nsamples)/
            (matrix[0][*n_loci].nsamples - (double)1.)),matrix[0][*n_loci].biallsitesn,p);/*not shared*/
        matrix[0][*n_loci].fuliFn  = fl_f2(matrix[0][*n_loci].nsamples,
            ((matrix[0][*n_loci].theta_fulin)*((double)matrix[0][*n_loci].nsamples)/
            (matrix[0][*n_loci].nsamples - (double)1.)),matrix[0][*n_loci].biallsitesn, 
            matrix[0][*n_loci].theta_tajn,p);/*not shared*/
    }
    else {
        matrix[0][*n_loci].fuliD = (double)-10000;
        matrix[0][*n_loci].fuliF = (double)-10000;
        matrix[0][*n_loci].faywuH = (double)-10000;
        matrix[0][*n_loci].faywuHo = (double)-10000;
        matrix[0][*n_loci].zengE = (double)-10000;
        
        if(matrix[0][*n_loci].biallsitesn > 1) {
            matrix[0][*n_loci].rZA = matrix[0][*n_loci].za/(matrix[0][*n_loci].biallsitesn - (double)1.);
            matrix[0][*n_loci].wB  = (double)matrix[0][*n_loci].b/(matrix[0][*n_loci].biallsitesn - (double)1.);
            matrix[0][*n_loci].wQ  = (double)matrix[0][*n_loci].q/(double)matrix[0][*n_loci].biallsitesn;
        }
        else {
            matrix[0][*n_loci].rZA = (double)-10000;
            matrix[0][*n_loci].wB  = (double)-10000;
            matrix[0][*n_loci].wQ  = (double)-10000;
        }
        matrix[0][*n_loci].fuliDn  = fl_d2(matrix[0][*n_loci].nsamples,
            (matrix[0][*n_loci].theta_fulin)*((double)matrix[0][*n_loci].nsamples)/
            (matrix[0][*n_loci].nsamples - (double)1.),matrix[0][*n_loci].biallsitesn,p);
        matrix[0][*n_loci].fuliFn  = fl_f2(matrix[0][*n_loci].nsamples,
            (matrix[0][*n_loci].theta_fulin)*((double)matrix[0][*n_loci].nsamples)/
            (matrix[0][*n_loci].nsamples - (double)1.),matrix[0][*n_loci].biallsitesn, matrix[0][*n_loci].theta_tajn,p);
    }
    return 1;
}


double tajima_d(double k_, int S_, double *coef_taj)
{
	double an,ut,vt;
	double S_D;
        
        if(S_ == 0 || *(coef_taj+0) < (double)1.51) return((double) -10000); 
        
	an = *(coef_taj+0);
	ut = *(coef_taj+3);
	vt = *(coef_taj+2);

	S_D = (k_ - ((double)S_/an)) / ((double)sqrt((ut*(double)S_) + (vt*(double)S_*((double)S_))));
	
	if (fabs(S_D) < (double)1.0E-15)
		S_D = (double)0.0;

	return S_D;
}

double fl_d(/*int sample_size,*/int fr1,int S, double *coef) /* fu and Li D amb outgroup */
{
	double an;
	double ud,vd;
	int re;
	double D;
		
        if(S == 0 || *(coef+0) < (double)1.5) return((double) -10000);
                
	re = fr1;	
	an = *(coef+0);
	
	vd = *(coef+8);
	ud = *(coef+9);
	D  = ((double)S - an*(double)re) / 
	     (double)sqrt(ud*(double)S + vd*(double)S*(double)S);

	return D;
}

double fl_d2(int sample_size,double fr1w,int S, double *coef) /* NO outgroup */
{
	double an;
	int n;
	double ud2,vd2;
	double rs;
	double D2;
        
        if(S == 0 || *(coef+0) < (double)1.51) return((double) -10000);
	
	rs = fr1w;

	n = sample_size;
	an = *(coef+0);
	
	vd2 = *(coef+4);
	ud2 = *(coef+5);
	D2  = ((double)S/an - (double)rs*(((double)n-(double)1.0)/(double)n)) /
	      (double)sqrt(ud2*(double)S + vd2*(double)S*(double)S);

	return D2;
}

double fl_f(/*int sample_size,*/int fr1, int S, double pi, double *coef) /* Fu and Li F amb outgroup */
{
	double uf,vf;
	int re;
	double F;
	
	if(S == 0 || *(coef+0) < (double)1.5) return((double) -10000);

	re = fr1;		
	vf = *(coef+10);
	uf = *(coef+11);

	F  = (pi - (double)re) / (double)sqrt(uf*(double)S + vf*(double)S*(double)S);

	return F;
}

double fl_f2(int sample_size,double fr1w, int S, double pi, double *coef) /* NO outgroup */
{
	int n;
	double uf2,vf2;
	double rs;
	double F2;
        
	if(S == 0 || *(coef+0) < (double)1.51) return((double) -10000);
	
	rs = fr1w;
	
	n   = sample_size;
	vf2 = *(coef+6);
	uf2 = *(coef+7);
	
	F2  = (pi - ((((double)n-(double)1.0)/(double)n)*(double)rs)) / 
	        (double)sqrt(uf2*(double)S + vf2*(double)S*(double)S);

	return F2;
}

double fay_wu_original(int sample_size,int *fr,double pi) /* Fay and Wu H nomes outgroup */
{
    int i;
    double Th,H;
    
    if(pi == (double)0.0 || sample_size < 2) return(0);/*EXCEPTION TO HAVE A VALUE IN theta. It should be (double) -10000...*/
    
    Th = (double)0.;
    for(i=1;i<sample_size;i++) Th += ((double)*(fr+i))*((double)i*(double)i);
    Th *= (double)2.0/((double)sample_size*(sample_size-(double)1));
    
    H = pi - Th;

    return H;
}

double fay_wu(int n,int *fr,double pi) /* Fay and Wu H nomes outgroup NORMALIZED (Zeng et al. Genetics 2006 174: 1431-9)*/
{
    int i;
    double TL,H,varpiTL,thetaw,an,bn,S;
    
    if(pi == (double)0.0 || n < 2) return(-10000);
    
    TL = thetaw = an = bn = (double)0;
	for(i=1;i<n;i++) {
		TL += ((double)*(fr+i))*((double)i);
		thetaw += (double)*(fr+i); 
		an += (double)1/(double)i;
		bn += (double)1.0/((double)i*(double)i);
	}
    TL *= (double)1.0/((double)(n-(double)1));
    S = thetaw;
	thetaw = thetaw/an;
	varpiTL = thetaw * ((double)(n-(double)2))/((double)6*((double)(n-(double)1))) + 
	          S*(S-(double)1)/(an*an+bn) * 
			  ((double)18*(double)n*(double)n*((double)3*(double)n+(double)2)*(bn+(double)1.0/((double)n*(double)n)) - 
			  ((double)88*(double)n*(double)n*(double)n + (double)9*(double)n*(double)n - 13*(double)n + (double)6)) /
			  ((double)9*((double)n*(n-(double)1)*(n-(double)1)));
	
	H = (pi - TL)/(double)sqrt(varpiTL);

    return H;
}

double fay_wu_obs(int n,double thetaL,double thetaw,double S,double *coef,double pi) /* Fay and Wu H nomes outgroup NORMALIZED (Zeng et al. Genetics 2006 174: 1431-9)*/
{
    double H,varpiTL,an,bn;
    
    if(pi == (double)0.0 || n < 2) return(-10000);
    an = coef[0];
	bn = coef[1];
	varpiTL = thetaw * ((double)(n-(double)2))/((double)6*((double)(n-(double)1))) + 
	          S*(S-(double)1)/(an*an+bn) * 
			  ((double)18*(double)n*(double)n*((double)3*(double)n+(double)2)*(bn+(double)1.0/((double)n*(double)n)) - 
			  ((double)88*(double)n*(double)n*(double)n + (double)9*(double)n*(double)n - 13*(double)n + (double)6)) /
			  ((double)9*((double)n*(n-(double)1)*(n-(double)1)));
	
	H = (pi - thetaL)/(double)sqrt(varpiTL);

    return H;
}

double E_zeng(int n,double thetaL,double thetaw,double S,double *coef) /* (Zeng et al. Genetics 2006 174: 1431-9)*/
{
    double E,varLW,an,bn;
    
    if(thetaw == (double)0.0 || n < 2) return(-10000);
    an = coef[0];
	bn = coef[1];
	varLW = thetaw * ((double)n/((double)2*(double)(n-1)) - (double)1/an) +
			S*(S-(double)1)/(an*an+bn) * 
			(bn/(an*an) + (double)2*bn*((double)n/(double)(n-1))*((double)n/(double)(n-1)) - 
			 (double)2*((double)n*bn-(double)n+(double)1)/((double)(n-1)*an) - 
			 ((double)3*(double)n+(double)1)/((double)(n-1)));
	
	E = (thetaL - thetaw)/(double)sqrt(varLW);

    return E;
}

/*R2 Ramos-Onsins&Rozas:"*unic" is the number of singletons in each sequence (in comparison to the sample studied)*/
double R2(int *unic,double pi,int sample_size,int S)
{
    double sm2 = 0.0;
    int i;
    
    if(S == 0 || sample_size == 0) return((double) -10000);
    for (i=0;i<sample_size;i++)
            sm2 += ((double)unic[i] - pi/2.0)*((double)unic[i] - pi/2.0);
    
    sm2 = sqrt(sm2/((double)sample_size))/(double)S;
            
    if (sm2 < 1.0E-15)
            sm2 = 0.0;

    return (double)sm2;
}

double Fs(int Nsample,double pi,int NumAlelos)
{
    /* Rozas program */
	
    double SumaP;
    double RestaP;
    int AleloI;
    unsigned long i;       
    double ValorFs;
    double *qew;
    double est_var;
    double FunEq23Ewens(int, int, double, double *);

    if(pi == 0.0 || Nsample < 2) return((double) -10000);	
    est_var = pi;

    if((qew  = (double *)malloc((unsigned long)Nsample*(unsigned long)Nsample*sizeof(double))) == 0)
        return (double) -10000;
    
    for(i=0;i<(unsigned long)Nsample*(unsigned long)Nsample;i++)
    	qew[i] = -1.0;

    SumaP=RestaP= (double)0;

    for (AleloI=1;AleloI<NumAlelos;AleloI++) {
        /* calculo q(n,aleloI)   ecuacion 21 (recurrente con eq. 19 y 20) */
        SumaP += FunEq23Ewens(Nsample, AleloI, est_var, qew);
    }

    if(SumaP > 1.-1E-37) {
    	for (AleloI = NumAlelos;AleloI <= Nsample; AleloI++)
            RestaP += FunEq23Ewens(Nsample, AleloI, est_var, qew);
            	 	
        if(RestaP < 1E-37) 
            return (double) -10000;   		
   	ValorFs = (double)log(RestaP) - (double)log(1.0-RestaP);
    }
    else {
	if(SumaP < 1E-37)
            return +10000;
        else
            ValorFs = (double)log(1.0-SumaP) - (double)log(SumaP);
    }
	
    if (fabs(ValorFs) < 1E-15)
        ValorFs = 0.0;
	    
    free(qew);	
    return ValorFs;
}

double FunEq23Ewens(int N,int i,double theta, double *qew_)
{                  
    unsigned long acceso; 
    int jj;
    double ValorN;  /* log del numerador */
    double ValorD;  /* log del denominador */

    acceso= (unsigned long)(N-1) * (unsigned long)N + (unsigned long)i - (unsigned long)1;    
    ValorN=0.0;
    ValorD=0.0;        
    if (qew_[acceso] < 0.0)  {   
        if (i==1) {
            /* calculo de qj,1   (i = 1)   Antigua equacion 19  */
            if (N > 2)  {             
                for (jj=2;jj<N;jj++)
                    ValorN = ValorN + (double)log((double)jj);  
            }
            ValorN = ValorN + (double)log(theta);
            for (jj=0;jj<N;jj++)
                ValorD  = ValorD + (double)log(theta + (double)jj);      
            qew_[acceso] = exp(ValorN - ValorD); 
        }    
        if (i==N) {          
            /* calculo de qj,j   (n = i)   antigua equacion 20 */
            ValorN = (double)log(theta) * (double)N;
            for (jj=0;jj<N;jj++)     
                ValorD  = ValorD + (double)log(theta + (double)jj);
            qew_[acceso] = exp(ValorN - ValorD);
        }
        if (i>1 && i<N) {    
            /*  recursividad  */
            qew_[acceso] = FunEq23Ewens( N - 1, i,theta, qew_) * ((double)(N - 1) / (theta + (double)N - 1.0))
                            + FunEq23Ewens( N - 1, i - 1, theta, qew_) * (theta / (theta + (double)N - 1.0));
        }    
    }  
    return(qew_[acceso]);
}

void init_coef(double *p,int sample_size)
{
/*
 Tajima and Fu coefficients.
*/
	int x,n;
	double an,bn,cn;
		
	an = bn = (double)0.0;
	n = sample_size;

	for(x=1;x<n;x++)
	{
		an += (double)1.0/(double)x;
		bn += (double)1.0/((double)x*(double)x);
	}
	
	p[0] = an;
	p[1] = bn;

	/* vt */
	p[2] = ((double)2.0*((double)n*(double)n + (double)n + (double)3.0)/((double)9.0*(double)n*((double)n-(double)1.0))
	        -((double)n+(double)2.0)/(an * (double)n) + bn/(an * an)) / (an*an + bn);
	        
	/* ut */
	p[3] = (( ((double)n+(double)1.0)/((double)3.0*((double)n-(double)1.0)) - (double)1.0/an)/an) - p[2];
	
	/* vd* */
	p[4] = (bn/(an*an) - (double)2.0/(double)n *((double)1.0 + (double)1.0/an - an + an/(double)n) - (double)1.0/((double)n*(double)n))
	       / (an*an + bn);
	
	/* ud* */
	p[5] = ((((double)n-(double)1.0)/(double)n - (double)1.0/an) / an) - p[4];
	
	/* vf* */
	p[6] = ((((double)2.0*(double)n*(double)n*(double)n + (double)110.0 * (double)n*(double)n - (double)255.0 * (double)n + (double)153.0)
	        / ((double)9.0 * (double)n * (double)n * ((double)n-(double)1.0)) + ((double)2.0*((double)n-(double)1.0) *an)/ ((double)n*(double)n) 
	        - ((double)8.0 * bn)/(double)n)) / (an*an + bn);
	
	/* uf* */
	p[7] = ((((double)4.0*(double)n*(double)n + (double)19.0*(double)n + (double)3.0 - (double)12.0 * ((double)n +(double)1.0) * (an + (double)1.0/(double)n))
	      / ((double)3.0*(double)n *((double)n-(double)1.0))) / an) - p[6];

	/* cn */
	cn = (double)2.0 * ((double)n*an-(double)2.0*((double)n-(double)1.0)) / (((double)n-(double)1.0)*((double)n-(double)2.0));
	
	/* vd */
	p[8] = (double)1.0 + (an*an/(bn+an*an)) * (cn - (((double)n+(double)1.0)/((double)n-(double)1.0)));
	
	/* ud* */
	p[9] = an -(double)1.0 - p[8];
	
	/* vf */
	p[10] = (cn + ((double)2.0*((double)n*(double)n+(double)n+(double)3.0))/((double)9.0*(double)n*((double)n-(double)1.0)) - (double)2.0
	       /((double)n-(double)1.0)) / (an*an + bn);
	
	/* uf */
	p[11] = ((double)1.0 + ((double)n+(double)1.0)/((double)3.0*((double)n-(double)1.0)) - (double)4*((double)n+(double)1.0)
	       /(((double)n-(double)1.0)*((double)n-(double)1.0)) * (an+((double)1.0/(double)n) - (double)2.0*(double)n/((double)n+(double)1.0)))
	       / an  -  p[10];
}
