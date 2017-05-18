/*
 *  hka.c
 *  MuLoNeTests
 *
 *  Created by sonsins on Wed Mar 19 2003. 
 *  Thanks to J. Schumacher (BGC-MPI Jena, Germany) to explain me how to do the algorithm.
 *  Biometry functions included
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define JCCorrection 1

int dohka_obs(struct statistics **matrix,struct statmulo **matrixml,int *outgroup, int *n_loci)
{
    int i,j,k,nlocihka;
    double *S;
    double *D,*a,*theta,T,b,*f;
	double coef[12];
    int hka_par(double *,double *, double *,double *,int, double **, double *);
	void init_coef(double *,int);

    if(*outgroup == 1 && *n_loci > 1) {
        matrixml[0][0].Shka = (double)0.0;
        matrixml[0][0].hka_T = (double)-10000;
    
        if((S = (double *) calloc(*n_loci,sizeof(double))) == 0) {
            return 0;
        }
        if((D = (double *) calloc(*n_loci,sizeof(double))) == 0) {
            free(S);
            return 0;
        }
        if((a = (double *) calloc(*n_loci,sizeof(double))) == 0) {
            free(S);
            free(D);
            return 0;
        }
        if((theta = (double *) calloc(*n_loci,sizeof(double))) == 0) {
            free(S);
            free(D);
            free(a);
            return 0;
        }

        if((f = (double *) calloc(*n_loci,sizeof(double))) == 0) {
            free(S);
            free(D);
            free(a);
            free(theta);
            return 0;
        }
        for(i=0,k=0;i<*n_loci;i++) {
			init_coef(coef,matrix[0][i].nsamples);
            a[k] = (double)0.;
			for(j=1;j<matrix[0][i].nsamples;j++) a[k] += (double)1./(double)j;
            f[k] = matrix[0][i].factor_chrn;
            /*with no correction... debugging.*/
			#if JCCorrection == 0
            S[k] = (double)matrix[0][i].biallsitesn;
            D[k] = matrix[0][i].ndivergence;
			#else
            /*or with correction of Jukes and Cantor.*/
            if(matrix[0][i].ndivergence/(double)matrix[0][i].nsites < (double)0.75) {
				if(matrix[0][i].ndivergence == (double)0) D[k] = (double)0;
				else D[k] = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * matrix[0][i].ndivergence/(double)matrix[0][i].nsites) * ((double)matrix[0][i].nsites + (double)matrix[0][i].nmhits);
                
				if(matrix[0][i].biallsitesn == 0) S[k] = (double)0;
				else { 
					/*S[k] = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * ((double)matrix[0][i].biallsitesn/coef[0])/(double)matrix[0][i].nsites) * (double)matrix[0][i].nsites * coef[0];*/
					S[k] = D[k]/matrix[0][i].ndivergence * (double)matrix[0][i].biallsitesn;
                }
            }
			#endif
            /*in case S and D data are zero, not counting*/
            if((matrix[0][i].biallsitesn > 0 || matrix[0][i].ndivergence > (double)0) /**/&&
               (matrix[0][i].ndivergence/(double)matrix[0][i].nsites < (double)0.75)/**/){
                k++;
            }
        }        
        nlocihka = k;/*substracting loci with no variation*/  
        
        if(nlocihka) {
            if(!(hka_par(S,D,a,f,nlocihka/**n_loci*/,&theta,&T))) {
				free(S);
				free(D);
				free(a);
				free(theta);
				free(f);
				return 0;
			}
        }else {
			free(S);
			free(D);
			free(a);
			free(theta);
			free(f);
			return 0;
        }
        for(i=0,k=0;i<*n_loci;i++) {
            if((matrix[0][i].biallsitesn > 0 || matrix[0][i].ndivergence > (double)0) /**/&&
               (matrix[0][i].ndivergence/(double)matrix[0][i].nsites < (double)0.75)/**/){
                b = (double)0.;
                for(j=1;j<matrix[0][i].nsamples;j++) b += (double)1./((double)j*(double)j);
				matrix[0][i].hka = 
                    (S[k] - theta[k]*f[k]*a[k]) * (S[k] - theta[k]*f[k]*a[k]) / (theta[k]*f[k]*a[k] + theta[k]*f[k]*theta[k]*f[k]*b) 
                    +
                    (D[k] - theta[k]*(T+(double)1*f[k])) * (D[k] - theta[k]*(T+(double)1*f[k])) / (theta[k]*(T+(double)1*f[k]) + theta[k]*theta[k]*((double)1*f[k])*((double)1*f[k]));
				/*keep resuts and data*/
				matrix[0][i].Sexphka = theta[k]*f[k]*a[k];
				matrix[0][i].Dexphka = theta[k]*(T+(double)1*f[k]);
				matrix[0][i].Sobshka = S[k];
				matrix[0][i].Dobshka = D[k];
				matrixml[0][0].Shka += matrix[0][i].hka;
				matrix[0][i].hka_theta = theta[k];
                k++;
            }
            /*in case S[i]+D[i] is zero, we decide to include in the analysis, and Chi-square at this locus is zero.*/
            else {
				matrix[0][i].hka = (double)0.;
				matrix[0][i].hka_theta = (double)0;
			}            
        }
        matrixml[0][0].hka_T = T;
                
		free(S);
        free(D);
        free(a);
        free(theta);
        free(f);
        
        return 1;
    }
    else return 0;
    
    return 0;
}

int dohka_sim(struct statistisim ***matrix,struct statistisimmuloc **matrixml,double *outgroup, int *n_loci,long int it)
{
    int i,j,k,nlocihka;
    static double *S = 0;
    static int nlocibefore;
    static double *D,*a,*theta,T,b,*f;
	double coef[12];
	void init_coef(double *,int);
    int hka_par(double *,double *, double *,double *,int, double **, double *);

    if(*outgroup >= (double)0. && *n_loci > 1) {
        matrixml[0][it].Shka = (double)0.0;
        matrixml[0][it].hka_T = (double)-10000;
        
        if(S == 0) {
            if((S = (double *) calloc(*n_loci,sizeof(double))) == 0) {
                return 0;
            }
            if((D = (double *) calloc(*n_loci,sizeof(double))) == 0) {
                free(S);
                S = 0;
                return 0;
            }
            if((a = (double *) calloc(*n_loci,sizeof(double))) == 0) {
                free(S);
                free(D);
                S = 0;
                return 0;
            }
            if((theta = (double *) calloc(*n_loci,sizeof(double))) == 0) {
                free(S);
                free(D);
                free(a);
                S = 0;
                return 0;
            }
            if((f = (double *) calloc(*n_loci,sizeof(double))) == 0) {
                free(S);
                free(D);
                free(a);
                free(theta);
                S = 0;
                return 0;
            }
            nlocibefore = *n_loci;
        }
        else if(nlocibefore != *n_loci) {
            if((S = (double *) realloc(S,(*n_loci)*sizeof(double))) == 0) {
                return 0;
            }
            if((D = (double *) realloc(D,(*n_loci)*sizeof(double))) == 0) {
                free(S);
                S = 0;
                return 0;
            }
            if((a = (double *) realloc(a,(*n_loci)*sizeof(double))) == 0) {
                free(S);
                free(D);
                S = 0;
                return 0;
            }
            if((theta = (double *) realloc(theta,(*n_loci)*sizeof(double))) == 0) {
                free(S);
                free(D);
                free(a);
                S = 0;
                return 0;
            }
            if((f = (double *) realloc(f,(*n_loci)*sizeof(double))) == 0) {
                free(S);
                free(D);
                free(a);
                free(theta);
                S = 0;
                return 0;
            }
        }

        for(i=0,k=0;i<*n_loci;i++) {
            a[k] = (double)0.;
			init_coef(coef,matrix[0][i][it].nsamples);
            for(j=1;j<matrix[0][i][it].nsamples;j++) a[k] += (double)1./(double)j;
            f[k] = matrix[0][i][it].factor_chrn;
			#if JCCorrection == 0
            /*with no correction for multiple hits*/
            S[k] = (double)matrix[0][i][it].biallsites;
            D[k] = matrix[0][i][it].ndivergence;
			#else
            /*or with correction of Jukes and Cantor.*/
            if(matrix[0][i][it].ndivergence/(double)matrix[0][i][it].nsites < (double)0.75) {
				if(matrix[0][i][it].ndivergence == (double)0) D[k] = (double)0;
				else D[k] = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * matrix[0][i][it].ndivergence/(double)matrix[0][i][it].nsites) * ((double)matrix[0][i][it].nsites + (double)matrix[0][i][it].mhsites);
                if(matrix[0][i][it].biallsites == 0) S[k] = (double)0;
				else {
					/*S[k] = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * ((double)matrix[0][i][it].biallsites/coef[0])/(double)matrix[0][i][it].nsites) * (double)matrix[0][i][it].nsites * coef[0];*/
					S[k] = D[k]/matrix[0][i][it].ndivergence * (double)matrix[0][i][it].biallsites;
				}
            }
			#endif
			/*in case S and D data are zero ot divergence too high, not counting*/
            if((matrix[0][i][it].biallsites > 0 || matrix[0][i][it].ndivergence > (double)0) /**/&&
               (matrix[0][i][it].ndivergence/(double)matrix[0][i][it].nsites < (double)0.75)/**/){
                k++;
            }
			/*printf("\nS[%d] = %f\t D[%d] = %f",k-1,S[k-1],k-1,D[k-1]);*/
        }
        nlocihka = k;/*substracting loci with no variation*/  
        
        if(nlocihka) {
			if(!(hka_par(S,D,a,f,nlocihka/**n_loci*/,&theta,&T))) 
				return 0;
        } else return 0;
        
        for(i=0,k=0;i<*n_loci;i++) {
            if((matrix[0][i][it].biallsites > 0 || matrix[0][i][it].ndivergence > (double)0) /**/&&
               (matrix[0][i][it].ndivergence/(double)matrix[0][i][it].nsites < (double)0.75)/**/){
                b = (double)0.;
                for(j=1;j<matrix[0][i][it].nsamples;j++) b += (double)1./((double)j*(double)j);
                matrix[0][i][it].hka = 
                    (S[k] - theta[k]*f[k]*a[k]) * (S[k] - theta[k]*f[k]*a[k]) / (theta[k]*f[k]*a[k] + theta[k]*f[k]*theta[k]*f[k]*b) 
                    +
                    (D[k] - theta[k]*(T+(double)1*f[k])) * (D[k] - theta[k]*(T+(double)1*f[k])) / (theta[k]*(T+(double)1*f[k]) + theta[k]*theta[k]*((double)1*f[k])*((double)1*f[k]));
                matrixml[0][it].Shka += matrix[0][i][it].hka;
				matrix[0][i][it].hka_theta = theta[k];
                k++;
            }
            /*in case S[i]+D[i] is zero, we decide to do the hka, but the Chi-square at this locus is zero.*/
            else {
				matrix[0][i][it].hka = (double)0.;
				matrix[0][i][it].hka_theta = (double)0;
			}
        }
        matrixml[0][it].hka_T = T;
        
        return 1;
    }
    else return 0;
    
    return 0;
}

int hka_par(double *S,double *D, double *a,double *f,int nloci,double **theta,double *T)
{
    int i;
    double x1,x2;
    double xacc;
    double functionT(double,double *,double *,double *, double *,int);
    int zbrac(double *,double *,double *,double *,double *, double *,int);
    double zriddr(double,double,double,double *,double *,double *, double *,int);
    
    if(functionT(-0.5,S,D,f,a,nloci) < (double)0) { /*check if negative at T=-0.5*/
        /*calculate the rangs*/
        x1 = (double)0.;
        x2 = 10.;
        if(zbrac(&x1,&x2,S,D,f,a,nloci) == 1) {
            /*estimate the value of T*/
            xacc = (double)1e-6*x2; /* accuracy of five decimals*/
            *T = (double)zriddr(x1,x2,xacc,S,D,a,f,nloci);
            if(*T == -1.11e30) return 0;
            /*estimate all thetas*/
            for(i=0;i<nloci;i++)
                theta[0][i] = (S[i] + D[i]) / ((*T+(double)1.0*f[i])+a[i]*f[i]);
        }
        else return 0;/*if not, that means that T is less than 0.5*/
    }
    else return 0;

    return 1;
}

double functionT(double T,double *S,double *D, double *a,double *f,int nloci) /*solution for T is obtained when the result is 0*/
{
	/*solution for T is obtained when the result is 0*/
	/*isolate theta_L for sum(Si) equation*/
	/*isolate theta_L for sum/Di) equation*/
	/*isolate theta_i for Si+Di equation*/
	/*equality between eq1 and eq2 and left only T by including eq3 in both*/
    int i;
    double fdT;
    double s1,s2,s3,s4;
    
    s1 = s2 = s3 = s4 = (double)0.0;
    for(i=0;i<nloci-1;i++) {
        s1 += S[i]/(a[nloci-1]*f[nloci-1]);
        s2 += (S[i] + D[i])/(a[i]*f[i]+(T+(double)1*f[i])) * a[i]*f[i];
        s3 += D[i]/(T+(double)1*f[i]);
		s4 += (S[i] + D[i])/(a[i]*f[i]+(T+(double)1*f[i]));
    }
    s1 += S[nloci-1]/(a[nloci-1]*f[nloci-1]);
	s2  = s2/(a[nloci-1]*f[nloci-1]);
	s3 += D[nloci-1]/(T+(double)1*f[nloci-1]);
	s3  = (s3 - s4);
	
    fdT = s1 - s2 - s3;
	
    return fdT;
}

int zbrac(double *x1,double *x2,double *S,double *D, double *a,double *f,int nloci)
{
	/* Based on Numerical Recipes in C. Press et al. 1992. 
    We need *x1 and *x2 be the range where a root is within them. We expand geometrically the range until finding
	(one a positive and one a negative value). If not, return 0.
	*/

    double f1,f2;
    int k=60;
    double functionT(double,double *,double *,double *, double *,int);
    
    if(*x1 == *x2) return 0;

    f1 = functionT(*x1,S,D,a,f,nloci);
    f2 = functionT(*x2,S,D,a,f,nloci);
	
	if(f1*f2 < (double)0) return 1;

    while(k--) {
        if(fabs(f1) < fabs(f2)) {
            *x1 += (double)1.5 * (*x1 - *x2);
            f1 = functionT(*x1,S,D,a,f,nloci);
        }
        else {
            *x2 += (double)1.5 * (*x2 - *x1);
            f2 = functionT(*x2,S,D,a,f,nloci);
        }
        if(f1*f2 < (double)0) return 1;
    }
    return 0;
}

double zriddr(double xlow,double xhigh,double xacc,double *S,double *D, double *a,double *f,int nloci)
{
	/* Based on Numerical Recipes in C. Press et al. 1992., p. 358 an on
	Ridders, 1979, IEEE Transactions on Circuits and systems, Vol. Cas-26, No. 11, pp. 979-980.
	*/
    int k=60;
	double flow,fhigh;
	double f1,f2,f3,f4;
	double x1,x2,x3,x4;
	double den,num,nsign;
    double functionT(double,double *,double *,double *, double *,int);
	

    flow  = functionT(xlow,S,D,a,f,nloci);
    fhigh = functionT(xhigh,S,D,a,f,nloci);

	if(flow  == (double)0) return xlow;
	if(fhigh == (double)0) return xhigh;
	if(flow*fhigh > (double)0) 
		return (double)-1e32;
	
	x1 = xlow;
	x2 = xhigh;
	f1 = flow;
	f2 = fhigh;
		
	while(k--) {
		x3 = (x1+x2)/(double)2;
		f3 = functionT(x3,S,D,a,f,nloci);
		if(f1 - f2 < (double)0) nsign = (double)-1;
		else nsign = (double)1;
		num = (x3-x1) * f3 * nsign;
		den = (double)sqrt((double)f3*(double)f3 - (double)f1*(double)f2);
		if(den <= xacc && -den <= xacc) return x3;
		x4 = x3 + num/den;
		f4 = functionT(x4,S,D,a,f,nloci);
		if(f4 <= xacc && -f4 <= xacc) return x4;
		if(f3*f4<(double)0) {
			x1 = x3;
			f1 = f3;
			x2 = x4;
			f2 = f4;
		}
		else {
			if(f1*f4<(double)0) {
				x2 = x4;
				f2 = f4;
			}
			else {
				if(f2*f4<(double)0) {
					x1 = x4;
					f1 = f4;
				}
			}
		}
		if(fabs(x1-x2) <= xacc) return x1;
	}	
	return (double)-1e32;
}

void print_prob_obshka(struct statmulo *matrixml, int *n_loci, FILE *file_output)
{

    double beta;
/*
    double z;
    double table_chisquare(int,double);
    double table_standardnormal(double);
    double chitonormal(int,double);
    double qgaus(double,double);
*/
    double probQ_chisquare(int,double);
    
    beta = probQ_chisquare(*n_loci-1,matrixml[0].Shka);
    
    if(beta > (double)0.05) {
        if(!file_output) printf("Significance of Chi-square: Non-Significant, P(dgf=%d) = %.3f\n",*n_loci-1,beta);
        if(file_output) 
            fprintf(file_output,"Significance of Chi-square: Non-Significant, P(dgf=%d) = %.3f\n",*n_loci-1,beta);
    }
    else {
        if(!file_output) printf("Significance of Chi-square: Significant, P(dgf=%d) = %.3f\n",*n_loci-1,beta);
        if(file_output)
            fprintf(file_output,"Significance of Chi-square: Significant, P(dgf=%d) = %.3f\n",*n_loci-1,beta);
    }

/*    
    if(!(file_output)) {
        //probability.
        if(*n_loci < 31) {//chi-square
            beta = table_chisquare(*n_loci-1,0.001);
            if(beta < matrixml[0].Shka) 
                printf("Significance of Chi-square: Significant, Chi(dgf=%d,P=%.3f) = %.3f\n",*n_loci-1,0.001,beta);
            else {
                beta = table_chisquare(*n_loci-1,0.01);
                if(beta < matrixml[0].Shka)
                    printf("Significance of Chi-square: Significant, Chi(dgf=%d,P=%.2f) = %.3f\n",*n_loci-1,0.01,beta);
                else {
                    beta = table_chisquare(*n_loci-1,0.05);
                    if(beta < matrixml[0].Shka)
                        printf("Significance of Chi-square: Significant, Chi(dgf=%d,P=%.2f) = %.3f\n",*n_loci-1,0.05,beta);
                    else {
                        beta = table_chisquare(*n_loci-1,0.01);
                        if(beta < matrixml[0].Shka)
                            printf("Significance of Chi-square: Non-Significant, Chi(dgf=%d,P=%.2f) = %.3f\n",*n_loci-1,0.05,beta);
                        else
                            printf("Significance of Chi-square: Non-Significant, Chi(dgf=%d,P=%.1f) = %.3f\n",*n_loci-1,0.10,beta);
                    }
                }
            }
        }
        else {//gaussian approach for dgf > 30.
            z = chitonormal(*n_loci-1,matrixml[0].Shka);
*/
    /*
            beta = (double)1.0 - qgaus(z,-z);
            if(beta > (double)0.05)
                printf("Significance of Chi-square: Non-Significant, P(dgf=%d) = %.3f\n",*n_loci-1,beta);
            else
                printf("Significance of Chi-square: Significant, P(dgf=%d) = %.3f\n",*n_loci-1,beta);
    */
    /*
            beta = table_standardnormal(0.001);
            if(beta < matrixml[0].Shka) 
                printf("Significance of HKA test: Significant at %.3f percent",0.001);
            else {
                beta = table_standardnormal(0.01);
                if(beta < matrixml[0].Shka)
                    printf("Significance of HKA test: Significant at %.2f percent",0.01);
                else {
                    beta = table_standardnormal(0.05);
                    if(beta < matrixml[0].Shka)
                        printf("Significance of HKA test: Significant at %.2f percent",0.05);
                    else {
                        beta = table_standardnormal(0.1);
                        if(beta < matrixml[0].Shka)
                            printf("Significance of Chi-square: Non-Significant, higher than %.3f\n",0.05);
                        else
                            printf("Significance of Chi-square: Non-Significant, higher than %.3f\n",0.10);
                    }
                }
            }
        */
/*
        }
    }
    else {
        //probability.
        if(*n_loci < 31) {//chi-square
            beta = table_chisquare(*n_loci-1,0.001);
            if(beta < matrixml[0].Shka) 
                fprintf(file_output,"Significance of Chi-square: Significant, Chi(dgf=%d,P=%.3f) = %.3f\n",*n_loci-1,0.001,beta);
            else {
                beta = table_chisquare(*n_loci-1,0.01);
                if(beta < matrixml[0].Shka)
                    fprintf(file_output,"Significance of Chi-square: Significant, Chi(dgf=%d,P=%.2f) = %.3f\n",*n_loci-1,0.01,beta);
                else {
                    beta = table_chisquare(*n_loci-1,0.05);
                    if(beta < matrixml[0].Shka)
                        fprintf(file_output,"Significance of Chi-square: Significant, Chi(dgf=%d,P=%.2f) = %.3f\n",*n_loci-1,0.05,beta);
                    else {
                        beta = table_chisquare(*n_loci-1,0.01);
                        if(beta < matrixml[0].Shka)
                            fprintf(file_output,"Significance of Chi-square: Non-Significant, Chi(dgf=%d,P=%.2f) = %.3f\n",*n_loci-1,0.05,beta);
                        else
                            fprintf(file_output,"Significance of Chi-square: Non-Significant, Chi(dgf=%d,P=%.1f) = %.3f\n",*n_loci-1,0.10,beta);
                    }
                }
            }
        }
        else {//gaussian approach for dgf > 30.
            z = chitonormal(*n_loci-1,matrixml[0].Shka);
*/
    /*
            beta = (double)1.0 - qgaus(z,-z);
            if(beta > (double)0.05)
                fprintf(file_output,"Significance of Chi-square: Non-Significant, P(dgf=%d) = %.3f\n",*n_loci-1,beta);
            else
                fprintf(file_output,"Significance of Chi-square: Significant, P(dgf=%d) = %.3f\n",*n_loci-1,beta);
    */
    /*
            beta = table_standardnormal(0.001);
            if(beta < matrixml[0].Shka) 
                fprintf(file_output,"Significance of HKA test: Significant at %.3f percent",0.001);
            else {
                beta = table_standardnormal(0.01);
                if(beta < matrixml[0].Shka)
                    fprintf(file_output,"Significance of HKA test: Significant at %.2f percent",0.01);
                else {
                    beta = table_standardnormal(0.05);
                    if(beta < matrixml[0].Shka)
                        fprintf(file_output,"Significance of HKA test: Significant at %.2f percent",0.05);
                    else {
                        beta = table_standardnormal(0.1);
                        if(beta < matrixml[0].Shka)
                            fprintf(file_output,"Significance of Chi-square: Non-Significant, higher than %.3f\n",0.05);
                        else
                            fprintf(file_output,"Significance of Chi-square: Non-Significant, higher than %.3f\n",0.10);
                    }
                }
            }
        */
/*
        }
    }
*/
    return;
}

