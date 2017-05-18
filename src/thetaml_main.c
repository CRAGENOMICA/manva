/*thetaml_main.c created by Sebastian E. Ramos-Onsins. April 7th 2004*/

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#define NLOCISEQ 20

/* 1:Tavare eq.  0:Hudson recursive eq.*/
#define EQTAVARE 1
#define RNDSIM 500
#define MHSIM 500
#define INITSEEDM 269764

struct parameters_thetaml {
    int nloci;
    long int *numberloci;
    int *nsamA;
    double *factor_chrn;
    int *SA;
    double *length;
    double thetamin_user;
    double thetamax_user;
    int steps_user;
    long int ncoalsim;
    long int seed;
};

struct Prob {
    double probc12;
    double probc13;
    double probc14;
    double probc1l;
    double probc23;
    double probc24;
    double probc2l;
    double probc34;
    double probc3l;
    double probc4l;
};

int thetamlS(struct globvar **global,struct var **datams,struct var2b **inputms,struct statistics *matrix,struct LRTdist **lrdist/*,struct MLthetasim ***mthetasim*/,FILE *file_output)
{
    long int i;
    int j,k,l,m,n;
    double thetaA,thetaM1,mlM1;
    
    struct parameters_thetaml *data;

    double **ml;
    double *thetalocusml;
    double *locusml;
    double *sumlocus;
    int *sortedlociml;
    int *stillloci;
    double LR1l=0./*,probLRT1l*/;
    int  dgf1l;
    
    int *vector2ml1=0;
    int *vector2ml2=0;
    int nloci2m1=0,nloci2m2=0;
    double *thetaM2=0;
    double *mlM2=0;
    double LR12=0./*,probLRT12*/;
    int  dgf12;
    double LR2l=0./*,probLRT2l*/;
    int  dgf2l;
    double sortthetamlmin;
    int sortlocimin;
    int findml2sorted(int,int *,int *,int,int,double *,double *,double **);
    int compare(const void *,const void *);
    
    int *vector3ml1=0;
    int *vector3ml2=0;
    int *vector3ml3=0;
    double *thetaM3=0;
    double *mlM3=0;
    int nloci3m1=0,nloci3m2=0,nloci3m3=0;
    double LR23=0.,LR3l=0.;
    /*double probLRT23,probLRT3l;*/
    int dgf23,dgf3l;
    int findml3sorted(int,int,int *,int *,int *,int,int,double *,double *,double **);
    int findml3sortedmh(int,int,int *,int *,int *,int *,int *,int,int,double *,
        double *,double *,double *,double **,double);
    double *relthetaM3=0;
    double *relmlM3=0;
    int relnloci3m1=0,relnloci3m2=0;
	double relM3totc=0.;
    int vnloci3m1=0,vnloci3m2=0;
    int relvnloci3m1=0,relvnloci3m2=0;
    
    int *vector4ml1=0;
    int *vector4ml2=0;
    int *vector4ml3=0;
    int *vector4ml4=0;
    double *thetaM4=0;
    double *mlM4=0;
    int nloci4m1=0,nloci4m2=0,nloci4m3=0,nloci4m4=0;
    double LR34=0.,LR4l=0.;
    /*double probLRT34,probLRT4l;*/
    int dgf34,dgf4l;
    int findml4sorted(int,int,int,int *,int *,int *,int *,int,int,double *,double *,double **);
    int findml4sortedmh(int,int,int,int *,int *,int *,int *,int *,int *,int *,int,int,double *,
        double *,double *,double *,double **,double);
    double *relthetaM4=0;
    double *relmlM4=0;
    int relnloci4m1=0,relnloci4m2=0,relnloci4m3=0;
	double relM4totc=0.;
    int vnloci4ml1=0,vnloci4ml2=0,vnloci4ml3=0;
    int relvnloci4ml1=0,relvnloci4ml2=0,relvnloci4ml3=0;
    
    double LR13=0.,LR14=0.,LR24=0.;
    /*double probLRT13,probLRT14,probLRT24;*/
    int dgf13,dgf14,dgf24;
	
    int suml1,suml2,suml3;
    double ran,logranloci;
    double ran1();
    void init_seed1(long int);
    int printdot;    
    
    struct Prob Pcoal;
    int probcoalthhaetasml(struct parameters_thetaml *,double,double *,double *,double *,int *,int *,int *,int *,int *,
                            int *,int *,int *,int *,int,int,int,int,int,int,int,int,int,double,double,double,double,double,
                            double,double,double,double,double,struct Prob *,struct LRTdist ** /*,struct MLthetasim ****/);

    double probQ_chisquare(int,double);
    int input_datathetasml(struct var **, struct var2b **,struct parameters_thetaml **);
    #if EQTAVARE == 0
    double **Q,**P;
    int maxS,maxn;
    double dfhudson(int,double,int,double ***,double ***);
    #else
    double dftavare(int,double,int);
    #endif
	double valuedf;
	
	double AIC1,AIC2,AIC3,AIC4,AICn,AICmin,sumwAIC;
	int bestAIC;
	
	if(file_output) {
		fprintf(file_output,"\n\n     MENU:\n     2. Estimation of parameters by Maximum likelihood:");
		fprintf(file_output,"\n     2.0. Estimation of levels of variation under neutral panmictic model with null recombination:\n\n");
		fprintf(file_output,"\n     2.0.1 Run program to estimate parameters.\n");
	}

	/*initialize*/	
    Pcoal.probc12 = (double)1;
    Pcoal.probc13 = (double)1;
    Pcoal.probc14 = (double)1;
    Pcoal.probc1l = (double)1;
    Pcoal.probc23 = (double)1;
    Pcoal.probc24 = (double)1;
    Pcoal.probc2l = (double)1;
    Pcoal.probc34 = (double)1;
    Pcoal.probc3l = (double)1;
    Pcoal.probc4l = (double)1;
    /*introduce data and define files*/
    if(!(data = (struct parameters_thetaml *)malloc(sizeof(struct parameters_thetaml)))) {
        puts("calloc error.0");
        return(1);
    }

    if((i=input_datathetasml(datams,inputms,&data)) != 0)
        return(1);

    /*define matrix and vectors*/
    init_seed1(data[0].seed);
	printdot = (int)ceil((double)data[0].nloci/(double)10);
    if(!(ml = (double **)malloc((data[0].nloci+1)*sizeof(double *)))) {
        puts("calloc error.0a");
        return(1);
    }
    for(j=0;j<=data[0].nloci;j++) {
        if(!(ml[j] = (double *)malloc((data[0].steps_user)*sizeof(double)))) {
            puts("calloc error.0b");
            return(1);
        }
    }
    if(!(thetalocusml = (double *)malloc((data[0].nloci+1)*sizeof(double)))) {
        puts("calloc error.0c");
        return(1);
    }
    if(!(locusml = (double *)malloc((data[0].nloci+1)*sizeof(double)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(sumlocus = (double *)malloc((data[0].nloci)*sizeof(double)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(sortedlociml = (int *)malloc((data[0].nloci)*sizeof(int)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(stillloci = (int *)malloc((data[0].nloci)*sizeof(int)))) {
        puts("calloc error.0d");
        return(1);
    }
    for(j=0;j<data[0].nloci;j++) stillloci[j] = j;

    /*vectors for M2 model*/
    if(!(thetaM2 = (double *)calloc(2,sizeof(double)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(mlM2 = (double *)calloc(2,sizeof(double)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(vector2ml1 = (int *)malloc((data[0].nloci + 1)*sizeof(int)))) {
        puts("calloc error.0d");
        return(1);
    }    
    if(!(vector2ml2 = (int *)malloc((data[0].nloci + 1)*sizeof(int)))) {
        puts("calloc error.0d");
        return(1);
    }
    
    /*vectors for M3 model*/
    if(!(thetaM3 = (double *)calloc(3,sizeof(double)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(mlM3 = (double *)calloc(3,sizeof(double)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(vector3ml1 = (int *)malloc((data[0].nloci)*sizeof(int)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(vector3ml2 = (int *)malloc((data[0].nloci)*sizeof(int)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(vector3ml3 = (int *)malloc((data[0].nloci)*sizeof(int)))) {
        puts("calloc error.0d");
        return(1);
    }

    /*vectors for M4 model*/
    if(!(thetaM4 = (double *)calloc(4,sizeof(double)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(mlM4 = (double *)calloc(4,sizeof(double)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(vector4ml1 = (int *)malloc((data[0].nloci)*sizeof(int)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(vector4ml2 = (int *)malloc((data[0].nloci)*sizeof(int)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(vector4ml3 = (int *)malloc((data[0].nloci)*sizeof(int)))) {
        puts("calloc error.0d");
        return(1);
    }
    if(!(vector4ml4 = (int *)malloc((data[0].nloci)*sizeof(int)))) {
        puts("calloc error.0d");
        return(1);
    }
    
    /*MHMCMC*/
    if(data[0].nloci > NLOCISEQ) {
        if(!(relthetaM3 = (double *)malloc(3*sizeof(double)))) {
            puts("calloc error.0d");
            return(1);
        }
        if(!(relmlM3 = (double *)malloc(3*sizeof(double)))) {
            puts("calloc error.0d");
            return(1);
        }
        if(!(relthetaM4 = (double *)malloc(4*sizeof(double)))) {
            puts("calloc error.0d");
            return(1);
        }
        if(!(relmlM4 = (double *)malloc(4*sizeof(double)))) {
            puts("calloc error.0d");
            return(1);
        }
    }

    #if EQTAVARE == 0
    /*Define max bound for S and n (useful for defining Q and P)*/
    maxS = maxn = 0;
    for(j=0;j<data[0].nloci;j++) {
        if(data[0].nsamA[j] > maxn) maxn = data[0].nsamA[j];
        if(data[0].SA[j] > maxS) maxS = data[0].SA[j];
    }
    
    /*Define the matrix Q[n][S] and P[n][S]*/
    if(!(Q = (double **)malloc((maxn+1)*sizeof(double *)))) {
        puts("calloc error.10");
        return(1);
    }
    if(!(P = (double **)malloc((maxn+1)*sizeof(double *)))) {
        puts("calloc error.10");
        return(1);
    }
    for(j=0;j<=maxn;j++) {
        if(!(Q[j] = (double *)malloc((maxS+1)*sizeof(double)))) {
            puts("calloc error.11");
            return(1);
        }
        if(!(P[j] = (double *)malloc((maxS+1)*sizeof(double)))) {
            puts("calloc error.11");
            return(1);
        }
    }
    #endif
    
    /*PRINT DATA ON THE SCREEN and define bounds of theta and file name*/
    printf("\nESTIMATION BY MAXIMUM LIKELIHOOD OF THETA GIVEN THE NUMBER OF SEGREGATING SITES.\n");
    printf("\nnumber of loci: %d",data[0].nloci);
    puts("\nnid_locus\tlength\tnsam\tS\n");
    for(j=0;j<data[0].nloci;j++) printf("%ld\t%.2f\t%d\t%d\n",data[0].numberloci[j],data[0].length[j],data[0].nsamA[j],
        data[0].SA[j]);
    printf("\ntheta/nt(min): %g\ntheta/nt(max): %g\n\n",data[0].thetamin_user,data[0].thetamax_user);
    
    /******************************************** MAIN PROGRAM *************************************************/
    /*Calculate the probabilities for the user bounds*/
    printf("\nCalculate maximum likelihood ");
    for(i=0;i<data[0].steps_user;i++) 
        ml[data[0].nloci][i] = (double)0.;/*init*/
    for(j=0;j<data[0].nloci;j++) {
        for(i=0;i<data[0].steps_user;i++) { 
            thetaA = ((double)i/((double)data[0].steps_user-(double)1.))*
                     (data[0].thetamax_user*(double)data[0].length[j] - data[0].thetamin_user*(double)data[0].length[j]) + 
                      data[0].thetamin_user*(double)data[0].length[j];/*uniform*/
            #if EQTAVARE == 0
            valuedf =  dfhudson(data[0].SA[j],thetaA*data[0].factor_chrn[j],data[0].nsamA[j],&Q,&P);
			if(valuedf <= (double)1e-323)
				ml[j][i] = (double)-1000;
			else 
				ml[j][i] = (double)log(valuedf);/*logprobability*/
            #else
            valuedf =  dftavare(data[0].SA[j],thetaA*data[0].factor_chrn[j],data[0].nsamA[j]);
			if(valuedf <= (double)1e-323)
				ml[j][i] = (double)-1000;
			else 
				ml[j][i] = (double)log(valuedf);/*logprobability*/
            #endif
            ml[data[0].nloci][i] += ml[j][i];/*sum of likelihood for each iteration of all loci*/
            if(i==0 || locusml[j] < ml[j][i]) {
                locusml[j] = ml[j][i];
                thetalocusml[j] = (double)i;
            }
        }
        if((double)j/(double)printdot == (int)j/(int)printdot) {
            printf(".");
            fflush(stdout);
        }
    }
    /*total maximum likelihood. M1 and Mnloci models*/
    locusml[data[0].nloci] = (double)0.;
    for(j=0;j<data[0].nloci;j++) {
        thetalocusml[j] = (thetalocusml[j]/((double)data[0].steps_user-(double)1.))*
                          (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user;
        locusml[data[0].nloci] += locusml[j];
		inputms[0][j].thetaml[5] = thetalocusml[j];
    }
	datams[0][0].mli_theta = locusml[data[0].nloci];

    mlM1 = thetaM1 = (double)0.;
    for(i=0;i<data[0].steps_user;i++) {
        if(i==0 || mlM1 < ml[data[0].nloci][i]) {
            mlM1 = ml[data[0].nloci][i];
            thetaM1 = ((double)i/((double)data[0].steps_user-(double)1.))*
                      (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user;
        }
    }
	datams[0][0].theta1_ml[0] = thetaM1;
	datams[0][0].ml1 = mlM1;
	for(j=0;j<data[0].nloci;j++) {
		inputms[0][j].thetaml[1] = datams[0][0].theta1_ml[0];
		inputms[0][j].mltheta = locusml[j];
	}
	
    if(data[0].nloci > 2) {
        /*loci are sorted by thetaM1 value*/
        for(j=0;j<data[0].nloci;j++) {
            sortthetamlmin = (double)0.;
            for(k=0;k<data[0].nloci;k++) {
                if(stillloci[k] > -1) {
                    if(sortthetamlmin == (double)0. || sortthetamlmin > thetalocusml[k]) {
                        sortthetamlmin = thetalocusml[k];
                        sortlocimin = k;
                    }
                }
            }
            stillloci[sortlocimin] = -1;
            sortedlociml[j] = sortlocimin;
        }
        /*M2 model*/
        /*do all lineal combinations and calculate ml*/
        printf("\nCalculation given all lineal combinations of maximum likelihood for two independent thetas ");
        nloci2m1 = -1;
        for(j=0;j<data[0].nloci-1;j++) {
            if(findml2sorted(j+1,sortedlociml,&nloci2m1,data[0].nloci,data[0].steps_user,thetaM2,mlM2,ml))
                return(1);
            if((double)j/(double)printdot == (int)j/(int)printdot) {
                printf(".");
                fflush(stdout);
            }
        }
        thetaM2[0] = ((double)thetaM2[0]/((double)data[0].steps_user-(double)1.))*
                    (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user;
        thetaM2[1] = ((double)thetaM2[1]/((double)data[0].steps_user-(double)1.))*
                    (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user;
        /*assign the vector2ml1 and 2*/
        for(j=0;j<nloci2m1;j++) 
            vector2ml1[j] = sortedlociml[j];
        vector2ml1[j] = -1;
        qsort(vector2ml1,j,sizeof(int),compare);
        for(j=nloci2m1,k=0;j<data[0].nloci;j++,k++) 
            vector2ml2[k] = sortedlociml[j];
        vector2ml2[k] = -1;
        qsort(vector2ml2,k,sizeof(int),compare);
        nloci2m2 = k;

		datams[0][0].theta2_ml[0] = thetaM2[0];
		datams[0][0].theta2_ml[1] = thetaM2[1];
		datams[0][0].ml2[0] = mlM2[0];
		datams[0][0].ml2[1] = mlM2[1];
		datams[0][0].ml2[2] = mlM2[0] + mlM2[1];
		datams[0][0].nloci2m[0] = nloci2m1;
		datams[0][0].nloci2m[1] = nloci2m2;

        j = 0;
        while(vector2ml1[j] > -1) {
			inputms[0][vector2ml1[j]].thetaml[2] = datams[0][0].theta2_ml[0];
			j++;
		}
        j = 0;
        while(vector2ml2[j] > -1) {
			inputms[0][vector2ml2[j]].thetaml[2] = datams[0][0].theta2_ml[1];
			j++;
		}
    }
    if(data[0].nloci > 3) {
        /*M3 model*/
        if(data[0].nloci <= NLOCISEQ) {
            /* all lineal combinations*/
            printf("\nCalculation given all lineal combinations of maximum likelihood for three independent thetas ");
            vnloci3m1 = vnloci3m2 = -1;
            for(j=0;j<data[0].nloci-2;j++) {
                for(k=j+1,l=0;k<data[0].nloci-1;k++,l++) {
                    if(findml3sorted(j+1,l+1,sortedlociml,&vnloci3m1,&vnloci3m2,data[0].nloci,
                        data[0].steps_user,thetaM3,mlM3,ml))
                        return(1);
                }
                if((double)j/(double)printdot == (int)j/(int)printdot) {
                    printf(".");
                    fflush(stdout);
                }
            }
        }
        else {
			/*Find randomly. First do lineal and then do randomly*/
			/*lineal given the groups of M2 model*/
			printf("\nHeuristic calculation of maximum likelihood for three independent thetas ");
			vnloci3m1 = vnloci3m2 = -1;
			/*random using iterations (RNDSIM/MHSIM)*/
			/*random values*/
			do {
				for(n=1;n<RNDSIM;n++) {
					nloci3m1 = (int)(ran1()*(double)(data[0].nloci-2) + (double)1);
					nloci3m2 = (int)(ran1()*(double)(data[0].nloci-1-nloci3m1) + (double)1);
					if(findml3sorted(nloci3m1,nloci3m2,sortedlociml,&vnloci3m1,&vnloci3m2,data[0].nloci,
						data[0].steps_user,thetaM3,mlM3,ml))
						return(1);
					if((double)RNDSIM/(double)n == (int)RNDSIM/n) {
						printf(".");
						fflush(stdout);
					}
				}
				/*calculate nloci3ml for the best ml*/
				nloci3m1 = vnloci3m1;
				nloci3m2 = vnloci3m2;
				/*assign best values to relative values (not the max but the current) (to do mhmcmc)*/
				relvnloci3m1 = vnloci3m1;
				relvnloci3m2 = vnloci3m2;
				relthetaM3[0] = thetaM3[0];
				relthetaM3[1] = thetaM3[1];
				relthetaM3[2] = thetaM3[2];
				relmlM3[0] = mlM3[0];
				relmlM3[1] = mlM3[1];
				relmlM3[2] = mlM3[2];
				relnloci3m1 = nloci3m1;
				relnloci3m2 = nloci3m2;
				suml1 = nloci3m1;
				suml2 = suml1 + nloci3m2;
				/*Now similar to 'MHMCMC'. values close to the best ml*/
				for(n=1;n<MHSIM;n++) {
					if((ran=ran1()) < (double)0.333) {
						if(suml1 == 1) suml1 = 2;
						else suml1 -= 1;
						if((ran=ran1()) < (double)0.333) {
							if(suml2 <= suml1 + 1) suml2 = suml1 + 1;
							else if(suml2 > suml1 + 1) suml2 -= 1;
						}
						else {
							if(ran > (double)0.666)
								if(suml2 < data[0].nloci - 2) suml2 += 1;
						}
					}
					else {
						if(ran > (double)0.666) {
							if(suml1 >= data[0].nloci - 3) {
								suml1 = data[0].nloci - 2;
								suml2 = suml1 + 1;
							}
							else {
								suml1 += 1;
								if((ran=ran1()) < (double)0.333) {
									if(suml2 <= suml1 + 1) suml2 = suml1 + 1;
									else suml2 -= 1;
								}
								if(ran > (double)0.666)
									if(suml2 < data[0].nloci - 2) suml2 += 1;
							}
						}
						else {
							if((ran=ran1()) < (double)0.333) {
								if(suml2 <= suml1 + 1) suml2 = suml1 + 1;
								else suml2 -= 1;
							}
							if(ran > (double)0.666)
								if(suml2 < data[0].nloci - 2) suml2 += 1;
						}
					}
					relnloci3m1 = suml1;
					relnloci3m2 = suml2 - suml1;

					logranloci = (double)log(ran1()/(double)data[0].nloci);
					relM3totc = relmlM3[0] + relmlM3[1] + relmlM3[2];
					if(findml3sortedmh(relnloci3m1,relnloci3m2,sortedlociml,&relvnloci3m1,&relvnloci3m2,
						&vnloci3m1,&vnloci3m2,data[0].nloci,data[0].steps_user,relthetaM3,relmlM3,thetaM3,
						mlM3,ml,logranloci))
						return(1);
					if(relM3totc == relmlM3[0] + relmlM3[1] + relmlM3[2]) {
						/*rejected*/
						relnloci3m1 = nloci3m1;
						relnloci3m2 = nloci3m2;
					}
					else {
						/*accepted*/
						nloci3m1 = relnloci3m1;
						nloci3m2 = relnloci3m2;
					}
					if((double)MHSIM/(double)n == (int)MHSIM/n) {
						printf(".");
						fflush(stdout);
					}
				}
			}while(datams[0][0].ml2[2] > mlM3[0]+mlM3[1]+mlM3[2]);
        }
        thetaM3[0] = ((double)thetaM3[0]/((double)data[0].steps_user-(double)1.))*
                    (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user;
        thetaM3[1] = ((double)thetaM3[1]/((double)data[0].steps_user-(double)1.))*
                    (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user;
        thetaM3[2] = ((double)thetaM3[2]/((double)data[0].steps_user-(double)1.))*
                    (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user;
        /*assign the vector3ml1,2 and 3*/
        for(j=0;j<vnloci3m1;j++) 
            vector3ml1[j] = sortedlociml[j];
        vector3ml1[j] = -1;
        qsort(vector3ml1,j,sizeof(int),compare);
        nloci3m1 = vnloci3m1;
        for(j=vnloci3m1,k=0;j<vnloci3m1+vnloci3m2;j++,k++) 
            vector3ml2[k] = sortedlociml[j];
        vector3ml2[k] = -1;
        qsort(vector3ml2,k,sizeof(int),compare);
        nloci3m2 = k;
        for(j=vnloci3m1+vnloci3m2,k=0;j<data[0].nloci;j++,k++) 
            vector3ml3[k] = sortedlociml[j];
        vector3ml3[k] = -1;
        qsort(vector3ml3,k,sizeof(int),compare);
        nloci3m3 = k;

		datams[0][0].theta3_ml[0] = thetaM3[0];
		datams[0][0].theta3_ml[1] = thetaM3[1];
		datams[0][0].theta3_ml[2] = thetaM3[2];
		datams[0][0].ml3[0] = mlM3[0]; 
		datams[0][0].ml3[1] = mlM3[1];
		datams[0][0].ml3[2] = mlM3[2];
		datams[0][0].ml3[3] = mlM3[1] + mlM3[1] + mlM3[2];
		datams[0][0].nloci3m[0] = nloci3m1;
		datams[0][0].nloci3m[1] = nloci3m2;
		datams[0][0].nloci3m[2] = nloci3m3;

        j = 0;
        while(vector3ml1[j] > -1) {
			inputms[0][vector3ml1[j]].thetaml[3] = datams[0][0].theta3_ml[0];
			j++;
		}
        j = 0;
        while(vector3ml2[j] > -1) {
			inputms[0][vector3ml2[j]].thetaml[3] = datams[0][0].theta3_ml[1];
			j++;
		}
        j = 0;
        while(vector3ml3[j] > -1) {
			inputms[0][vector3ml3[j]].thetaml[3] = datams[0][0].theta3_ml[2];
			j++;
		}
    }
    if(data[0].nloci > 4) {
        /*M4 model*/
        if(data[0].nloci <= NLOCISEQ) {
            /*heuristic lineal all*/
            printf("\nCalculation given all lineal combinations of maximum likelihood for four independent thetas ");
            vnloci4ml1 = vnloci4ml2 = vnloci4ml3 = -1;
            for(j=0;j<data[0].nloci-3;j++) {
                l=0;
                for(k=j+1,l=0;k<data[0].nloci-2;k++,l++) {
                    for(n=k+1,m=0;n<data[0].nloci-1;n++,m++) {
                        if(findml4sorted(j+1,l+1,m+1,sortedlociml,
                            &vnloci4ml1,&vnloci4ml2,&vnloci4ml3,data[0].nloci,data[0].steps_user,thetaM4,mlM4,ml))
                            return(1);
                    }
                }
                printf(".");
                fflush(stdout);
            }
        }
        else {
            /*Find randomly. First do lineal and then do randomly*/
            /*lineal given the groups of M3 model*/
            printf("\nHeuristic calculation of maximum likelihood for four independent thetas ");
            vnloci4ml1 = vnloci4ml2 = vnloci4ml3 = -1;
            do {
				/*random using iterations (RNDSIM/MHSIM)*/
				for(n=1;n<RNDSIM;n++) {/*random values*/
					nloci4m1 = (int)(ran1()*(double)(data[0].nloci-3)) + 1;
					nloci4m2 = (int)(ran1()*(double)(data[0].nloci-2-nloci4m1)) + 1;
					nloci4m3 = (int)(ran1()*(double)(data[0].nloci-1-nloci4m1-nloci4m2)) + 1;
					if(findml4sorted(nloci4m1,nloci4m2,nloci4m3,sortedlociml,
						&vnloci4ml1,&vnloci4ml2,&vnloci4ml3,data[0].nloci,data[0].steps_user,thetaM4,mlM4,ml))
						return(1);
					if((double)RNDSIM/(double)n == (int)RNDSIM/n) {
						printf(".");
						fflush(stdout);
					}
				}
				/*calculate nloci4ml for the best ml*/
				nloci4m1 = vnloci4ml1;
				nloci4m2 = vnloci4ml2; 
				nloci4m3 = vnloci4ml3; 
				/*assign best values to relative values (not the max but the current) (to do mhmcmc)*/
				relvnloci4ml1 = vnloci4ml1;
				relvnloci4ml2 = vnloci4ml2;
				relvnloci4ml3 = vnloci4ml3;
				relthetaM4[0] = thetaM4[0];
				relthetaM4[1] = thetaM4[1];
				relthetaM4[2] = thetaM4[2];
				relthetaM4[3] = thetaM4[3];
				relmlM4[0] = mlM4[0];
				relmlM4[1] = mlM4[1];
				relmlM4[2] = mlM4[2];
				relmlM4[3] = mlM4[3];
				relnloci4m1 = nloci4m1;
				relnloci4m2 = nloci4m2;
				relnloci4m3 = nloci4m3;
				suml1 = nloci4m1;
				suml2 = suml1 + nloci4m2;
				suml3 = suml2 + nloci4m3;
				/*Now similar to MHMCMC. values close to the best ml*/
				for(n=1;n<MHSIM;n++) {
					if((ran=ran1()) < (double)0.333) {
						if(suml1 > 1) suml1 -= 1;
						if((ran=ran1()) < (double)0.333) {
							if(suml2 <= suml1 + 1) suml2 = suml1 + 1;
							else suml2 -= 1;
							if((ran=ran1()) < (double)0.333) {
								if(suml3 <= suml2 + 1) suml3 = suml2 + 1;
								else suml3 -= 1;
							}
							else {
								if(ran > (double)0.666)
									if(suml3 < data[0].nloci - 2) suml3 += 1;
							}
						}
						else {
							if(ran > (double)0.666) {
								if(suml2 < data[0].nloci - 3) suml2 += 1;
								if((ran=ran1()) < (double)0.333) {
									if(suml3 <= suml2 + 1) suml3 = suml2 + 1;
									else suml3 -= 1;
								}
								else {
									if(ran > (double)0.666)
										if(suml3 < data[0].nloci - 2) suml3 += 1;
								}
							}
						}
					}
					else {
						if(ran > (double)0.666) {
							if(suml1 >= data[0].nloci - 4) {
								suml1 = data[0].nloci - 3;
								suml2 = suml1 + 1;
							}
							else {
								suml1 += 1;
								if((ran=ran1()) < (double)0.333) {
									if(suml2 <= suml1 + 1) suml2 = suml1 + 1;
									else suml2 -= 1;
									if((ran=ran1()) < (double)0.333) {
										if(suml3 <= suml2 + 1) suml3 = suml2 + 1;
										else suml3 -= 1;
									}
									else {
										if(ran > (double)0.666)
											if(suml3 < data[0].nloci - 2) suml3 += 1;
									}
								}
								if(ran > (double)0.666) {
									if(suml2 < data[0].nloci - 3) suml2 += 1;
									if((ran=ran1()) < (double)0.333) {
										if(suml3 <= suml2 + 1) suml3 = suml2 + 1;
										else suml3 -= 1;
									}
									else {
										if(ran > (double)0.666)
											if(suml3 < data[0].nloci - 2) suml3 += 1;
									}
								}
							}
						}
						else {
							if((ran=ran1()) < (double)0.333) {
								if(suml2 <= suml1 + 1) suml2 = suml1 + 1;
								else suml2 -= 1;
								if((ran=ran1()) < (double)0.333) {
									if(suml3 <= suml2 + 1) suml3 = suml2 + 1;
									else suml3 -= 1;
								}
								else {
									if(ran > (double)0.666)
										if(suml3 < data[0].nloci - 2) suml3 += 1;
								}
							 }
							if(ran > (double)0.666) {
								if(suml2 < data[0].nloci - 3) suml2 += 1;
								if((ran=ran1()) < (double)0.333) {
									if(suml3 <= suml2 + 1) suml3 = suml2 + 1;
									else suml3 -= 1;
								}
								else {
									if(ran > (double)0.666)
										if(suml3 < data[0].nloci - 2) suml3 += 1;
								}
							}
						}
					}
					relnloci4m1 = suml1;
					relnloci4m2 = suml2 - suml1;
					relnloci4m3 = suml3 - suml2;

					logranloci = (double)log(ran1()/(double)data[0].nloci);
					relM4totc = relmlM4[0] + relmlM4[1] + relmlM4[2] + relmlM4[3];
					if(findml4sortedmh(relnloci4m1,relnloci4m2,relnloci4m3,
						sortedlociml,&relvnloci4ml1,&relvnloci4ml2,&relvnloci4ml3,&vnloci4ml1,&vnloci4ml2,
						&vnloci4ml3,data[0].nloci,data[0].steps_user,relthetaM4,relmlM4,thetaM4,mlM4,
						ml,logranloci))
						return(1);
					if(relM4totc == relmlM4[0] + relmlM4[1] + relmlM4[2] + relmlM4[3]) {
						/*rejected*/
						relnloci4m1 = nloci4m1;
						relnloci4m2 = nloci4m2;
						relnloci4m3 = nloci4m3;
					}
					else {
						/*accepted*/
						nloci4m1 = relnloci4m1;
						nloci4m2 = relnloci4m2;
						nloci4m3 = relnloci4m3;
					}
					if((double)MHSIM/(double)n == (int)MHSIM/n) {
						printf(".");
						fflush(stdout);
					}
				}
			}while(datams[0][0].ml3[3] > mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]);
        }
        thetaM4[0] = ((double)thetaM4[0]/((double)data[0].steps_user-(double)1.))*
                    (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user;
        thetaM4[1] = ((double)thetaM4[1]/((double)data[0].steps_user-(double)1.))*
                    (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user;
        thetaM4[2] = ((double)thetaM4[2]/((double)data[0].steps_user-(double)1.))*
                    (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user;
        thetaM4[3] = ((double)thetaM4[3]/((double)data[0].steps_user-(double)1.))*
                    (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user;
        /*assign the vector4ml1,2,3 and 4*/
        for(j=0;j<vnloci4ml1;j++) 
            vector4ml1[j] = sortedlociml[j];
        vector4ml1[j] = -1;
        qsort(vector4ml1,j,sizeof(int),compare);
        nloci4m1 = vnloci4ml1;
        for(j=vnloci4ml1,k=0;j<vnloci4ml1+vnloci4ml2;j++,k++) 
            vector4ml2[k] = sortedlociml[j];
        vector4ml2[k] = -1;
        qsort(vector4ml2,k,sizeof(int),compare);
        nloci4m2 = k;
        for(j=vnloci4ml1+vnloci4ml2,k=0;j<vnloci4ml1+vnloci4ml2+vnloci4ml3;j++,k++) 
            vector4ml3[k] = sortedlociml[j];
        vector4ml3[k] = -1;
        qsort(vector4ml3,k,sizeof(int),compare);
        nloci4m3 = k;
        for(j=vnloci4ml1+vnloci4ml2+vnloci4ml3,k=0;j<data[0].nloci;j++,k++) 
            vector4ml4[k] = sortedlociml[j];
        vector4ml4[k] = -1;
        qsort(vector4ml4,k,sizeof(int),compare);
        nloci4m4 = k;

		datams[0][0].theta4_ml[0] = thetaM4[0];
		datams[0][0].theta4_ml[1] = thetaM4[1];
		datams[0][0].theta4_ml[2] = thetaM4[2];
		datams[0][0].theta4_ml[3] = thetaM4[3];
		datams[0][0].ml4[0] = mlM4[0];
		datams[0][0].ml4[1] = mlM4[1];
		datams[0][0].ml4[2] = mlM4[2];
		datams[0][0].ml4[3] = mlM4[3];
		datams[0][0].ml4[4] = mlM4[0] + mlM4[1] + mlM4[2] + mlM4[3];
		datams[0][0].nloci4m[0] = nloci4m1;
		datams[0][0].nloci4m[1] = nloci4m2;
		datams[0][0].nloci4m[2] = nloci4m3;
		datams[0][0].nloci4m[3] = nloci4m4;

        j = 0;
        while(vector4ml1[j] > -1) {
			inputms[0][vector4ml1[j]].thetaml[4] = datams[0][0].theta4_ml[0];
			j++;
		}
        j = 0;
        while(vector4ml2[j] > -1) {
			inputms[0][vector4ml2[j]].thetaml[4] = datams[0][0].theta4_ml[1];
			j++;
		}
        j = 0;
        while(vector4ml3[j] > -1) {
			inputms[0][vector4ml3[j]].thetaml[4] = datams[0][0].theta4_ml[2];
			j++;
		}
        j = 0;
        while(vector4ml4[j] > -1) {
			inputms[0][vector4ml4[j]].thetaml[4] = datams[0][0].theta4_ml[3];
			j++;
		}
    }
    
    /*likelihood ratio tests*/
    LR1l = (double)2. * (-mlM1 + locusml[data[0].nloci]);
    dgf1l = 2 * (data[0].nloci - 1);
    /*probLRT1l = probQ_chisquare(dgf1l,LR1l);*/

    if(data[0].nloci > 1) {
        LR12 = (double)2. * (-mlM1 + (mlM2[0]+mlM2[1]));
        dgf12 = 2 * (3 - 1);
        /*probLRT12 = probQ_chisquare(dgf12,LR12);*/
        
        if(data[0].nloci > 2) {
			LR2l = (double)2. * (-(mlM2[0]+mlM2[1]) + locusml[data[0].nloci]);
			if(data[0].nloci > 3) {
				dgf2l = 2 * (data[0].nloci - 3);
				/*probLRT2l = probQ_chisquare(dgf2l,LR2l);*/
			}
		}
    }
    if(data[0].nloci > 2) {
        LR13 = (double)2. * (-mlM1 + (mlM3[0]+mlM3[1]+mlM3[2]));
        dgf13 = 2 * (5 - 1);
        /*probLRT13 = probQ_chisquare(dgf13,LR13);*/

        LR23 = (double)2. * (-(mlM2[0]+mlM2[1]) + (mlM3[0]+mlM3[1]+mlM3[2]));
        dgf23 = 2 * (5 - 3);
        /*probLRT23 = probQ_chisquare(dgf23,LR23);*/
        
		if(data[0].nloci > 3) {
			LR3l = (double)2. * (-(mlM3[0]+mlM3[1]+mlM3[2]) + locusml[data[0].nloci]);
			if(data[0].nloci > 5) {
				dgf3l = 2 * (data[0].nloci - 5);
				/*probLRT3l = probQ_chisquare(dgf3l,LR3l);*/
			}
		}
    }
    if(data[0].nloci > 3) {
        LR14 = (double)2. * (-mlM1 + (mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]));
        dgf14 = 2 * (7 - 1);
        /*probLRT14 = probQ_chisquare(dgf14,LR14);*/
        
        LR24 = (double)2. * (-(mlM2[0]+mlM2[1]) + (mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]));
        dgf24 = 2 * (7 - 3);
        /*probLRT24 = probQ_chisquare(dgf24,LR24);*/
        
        LR34 = (double)2. * (-(mlM3[0]+mlM3[1]+mlM3[2]) + (mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]));
        dgf34 = 2 * (7 - 5);
        /*probLRT34 = probQ_chisquare(dgf34,LR34);*/
        
		if(data[0].nloci > 4) {
			LR4l = (double)2. * (-(mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]) + locusml[data[0].nloci]);
			if(data[0].nloci > 7) {
				dgf4l = 2 * (data[0].nloci - 7);
				/*probLRT4l = probQ_chisquare(dgf4l,LR4l);*/
			}
		}
    }
    /****************************************** END MAIN PROGRAM (EXCEPT COAL. SIMULATIONS) ***********************/

    /*PRINT THE DATA AND RESULTS*/
    /*ON THE SCREEN*/
    printf("\n\nMaximum likelihood estimate considering one theta.\n");
    puts("\nnid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM1/region\n");
    for(j=0;j<data[0].nloci;j++) {
        if(global[0][0].dataequalsim)
			printf("%ld:%s\t%.2f\t%.2f\t%d\t%d\t%g\n",data[0].numberloci[j],matrix[j].gene,data[0].factor_chrn[j],data[0].length[j],data[0].nsamA[j],data[0].SA[j],thetaM1*(double)data[0].length[j]);
		else 
			printf("%ld\t%.2f\t%.2f\t%d\t%d\t%g\n",data[0].numberloci[j],data[0].factor_chrn[j],data[0].length[j],data[0].nsamA[j],data[0].SA[j],thetaM1*(double)data[0].length[j]);
	}
    printf("\ntheta/nt\n %g",thetaM1);
    printf("\nML\n %g\n",mlM1);
    if(data[0].nloci > 2) {
        puts("\nMaximum likelihood estimate considering two groups of independent loci");
        j=0;
        while(vector2ml1[j] > -1) j++;
        printf("\nlocus_group1: (%d)\n",j);
        puts("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM2/region\n");
        j = 0;
        while(vector2ml1[j] > -1) {
            if(global[0][0].dataequalsim) 
				printf("%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector2ml1[j]],matrix[vector2ml1[j]].gene,data[0].factor_chrn[vector2ml1[j]],data[0].length[vector2ml1[j]],
					data[0].nsamA[vector2ml1[j]],data[0].SA[vector2ml1[j]],thetaM2[0]*(double)data[0].length[vector2ml1[j]]);
			else 
				printf("%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector2ml1[j]],data[0].factor_chrn[vector2ml1[j]],data[0].length[vector2ml1[j]],
					data[0].nsamA[vector2ml1[j]],data[0].SA[vector2ml1[j]],thetaM2[0]*(double)data[0].length[vector2ml1[j]]);
            j++;
        }
        printf("\nlocus_group2: (%d)\n",data[0].nloci-j);
        puts("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM2/region\n");
        k = 0;
        while(vector2ml2[k] > -1) {
			if(global[0][0].dataequalsim)
				printf("%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector2ml2[k]],matrix[vector2ml2[k]].gene,data[0].factor_chrn[vector2ml2[k]],data[0].length[vector2ml2[k]],
					data[0].nsamA[vector2ml2[k]],data[0].SA[vector2ml2[k]],thetaM2[1]*(double)data[0].length[vector2ml2[k]]);
            else 
				printf("%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector2ml2[k]],data[0].factor_chrn[vector2ml2[k]],data[0].length[vector2ml2[k]],
					data[0].nsamA[vector2ml2[k]],data[0].SA[vector2ml2[k]],thetaM2[1]*(double)data[0].length[vector2ml2[k]]);
            k++;
        }
        printf("\ntheta1ml/nt\ttheta2ml/nt\n%g\t%g\n",thetaM2[0],thetaM2[1]);
        printf("\nmlgroup1\tmlgroup2\tML\n%g\t%g\t%g\n",mlM2[0],mlM2[1],mlM2[0]+mlM2[1]);
    }
    if(data[0].nloci > 3) {
        puts("\nMaximum likelihood estimate considering three groups of independent loci");
        j=0;
        while(vector3ml1[j] > -1) j++;
        printf("\nlocus_group1: (%d)\n",j);
        puts("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM3/region\n");
        j = 0;
        while(vector3ml1[j] > -1) {
            if(global[0][0].dataequalsim)
				printf("%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector3ml1[j]],matrix[vector3ml1[j]].gene,data[0].factor_chrn[vector3ml1[j]],data[0].length[vector3ml1[j]],
					data[0].nsamA[vector3ml1[j]],data[0].SA[vector3ml1[j]],thetaM3[0]*(double)data[0].length[vector3ml1[j]]);
			else
				printf("%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector3ml1[j]],data[0].factor_chrn[vector3ml1[j]],data[0].length[vector3ml1[j]],
					data[0].nsamA[vector3ml1[j]],data[0].SA[vector3ml1[j]],thetaM3[0]*(double)data[0].length[vector3ml1[j]]);
            j++;
        }
        k=0;
        while(vector3ml2[k] > -1) k++;
        printf("\nlocus_group2: (%d)\n",k);
        puts("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM3/region\n");
        k = 0;
        while(vector3ml2[k] > -1) {
            if(global[0][0].dataequalsim) 
				printf("%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector3ml2[k]],matrix[vector3ml2[k]].gene,data[0].factor_chrn[vector3ml2[k]],data[0].length[vector3ml2[k]],
					data[0].nsamA[vector3ml2[k]],data[0].SA[vector3ml2[k]],thetaM3[1]*(double)data[0].length[vector3ml2[k]]);
			else
				printf("%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector3ml2[k]],data[0].factor_chrn[vector3ml2[k]],data[0].length[vector3ml2[k]],
					data[0].nsamA[vector3ml2[k]],data[0].SA[vector3ml2[k]],thetaM3[1]*(double)data[0].length[vector3ml2[k]]);
            k++;
        }
        printf("\nlocus_group3: (%d)\n",data[0].nloci-j-k);
        puts("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM3/region\n");
        l=0;
        while(vector3ml3[l] > -1) {
            if(global[0][0].dataequalsim) 
				printf("%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector3ml3[l]],matrix[vector3ml3[l]].gene,data[0].factor_chrn[vector3ml3[l]],data[0].length[vector3ml3[l]],
					data[0].nsamA[vector3ml3[l]],data[0].SA[vector3ml3[l]],thetaM3[2]*(double)data[0].length[vector3ml3[l]]);
			else
				printf("%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector3ml3[l]],data[0].factor_chrn[vector3ml3[l]],data[0].length[vector3ml3[l]],
					data[0].nsamA[vector3ml3[l]],data[0].SA[vector3ml3[l]],thetaM3[2]*(double)data[0].length[vector3ml3[l]]);
            l++;
        }
        printf("\ntheta1ml/nt\ttheta2ml/nt\ttheta3ml/nt\n%g\t%g\t%g\n",thetaM3[0],thetaM3[1],thetaM3[2]);
        printf("\nmlgroup1\tmlgroup2\tmlgroup3\tML\n%g\t%g\t%g\t%g\n",mlM3[0],mlM3[1],mlM3[2],mlM3[0]+mlM3[1]+mlM3[2]);
    }
    if(data[0].nloci > 4) {    
        puts("\nMaximum likelihood estimate considering four groups of independent loci");
        j=0;
        while(vector4ml1[j] > -1) j++;
        printf("\nlocus_group1: (%d)\n",j);
        puts("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM4/region\n");
        j = 0;
        while(vector4ml1[j] > -1) {
            if(global[0][0].dataequalsim)
				printf("%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml1[j]],matrix[vector4ml1[j]].gene,data[0].factor_chrn[vector4ml1[j]],data[0].length[vector4ml1[j]],
					data[0].nsamA[vector4ml1[j]],data[0].SA[vector4ml1[j]],thetaM4[0]*(double)data[0].length[vector4ml1[j]]);
			else
				printf("%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml1[j]],data[0].factor_chrn[vector4ml1[j]],data[0].length[vector4ml1[j]],
					data[0].nsamA[vector4ml1[j]],data[0].SA[vector4ml1[j]],thetaM4[0]*(double)data[0].length[vector4ml1[j]]);
            j++;
        }
        k=0;
        while(vector4ml2[k] > -1) k++;
        printf("\nlocus_group2: (%d)\n",k);
        puts("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM4/region\n");
        k = 0;
        while(vector4ml2[k] > -1) {
            if(global[0][0].dataequalsim)
				printf("%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml2[k]],matrix[vector4ml2[k]].gene,data[0].factor_chrn[vector4ml2[k]],data[0].length[vector4ml2[k]],
					data[0].nsamA[vector4ml2[k]],data[0].SA[vector4ml2[k]],thetaM4[1]*(double)data[0].length[vector4ml2[k]]);
			else
				printf("%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml2[k]],data[0].factor_chrn[vector4ml2[k]],data[0].length[vector4ml2[k]],
					data[0].nsamA[vector4ml2[k]],data[0].SA[vector4ml2[k]],thetaM4[1]*(double)data[0].length[vector4ml2[k]]);
            k++;
        }
        l=0;
        while(vector4ml3[l] > -1) l++;
        printf("\nlocus_group3: (%d)\n",l);
        puts("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM4/region\n");
        l = 0;
        while(vector4ml3[l] > -1) {
            if(global[0][0].dataequalsim)
				printf("%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml3[l]],matrix[vector4ml3[l]].gene,data[0].factor_chrn[vector4ml3[l]],data[0].length[vector4ml3[l]],
					data[0].nsamA[vector4ml3[l]],data[0].SA[vector4ml3[l]],thetaM4[2]*(double)data[0].length[vector4ml3[l]]);
			else
				printf("%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml3[l]],data[0].factor_chrn[vector4ml3[l]],data[0].length[vector4ml3[l]],
					data[0].nsamA[vector4ml3[l]],data[0].SA[vector4ml3[l]],thetaM4[2]*(double)data[0].length[vector4ml3[l]]);
            l++;
        }
        printf("\nlocus_group4: (%d)\n",data[0].nloci-j-k-l);
        puts("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM4/region\n");
        m = 0;
        while(vector4ml4[m] > -1) {
			if(global[0][0].dataequalsim)
				printf("%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml4[m]],matrix[vector4ml4[m]].gene,data[0].factor_chrn[vector4ml4[m]],data[0].length[vector4ml4[m]],
					data[0].nsamA[vector4ml4[m]],data[0].SA[vector4ml4[m]], thetaM4[3]*(double)data[0].length[vector4ml4[m]]);
            else
				printf("%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml4[m]],data[0].factor_chrn[vector4ml4[m]],data[0].length[vector4ml4[m]],
					data[0].nsamA[vector4ml4[m]],data[0].SA[vector4ml4[m]], thetaM4[3]*(double)data[0].length[vector4ml4[m]]);
            m++;
        }
        printf("\ntheta1ml/nt\ttheta2ml/nt\ttheta3ml/nt\ttheta4ml/nt\n%g\t%g\t%g\t%g\n",thetaM4[0],thetaM4[1],thetaM4[2],thetaM4[3]);
        printf("\nmlgroup1\tmlgroup2\tmlgroup3\tmlgroup4\tML\n%g\t%g\t%g\t%g\t%g\n",mlM4[0],mlM4[1],mlM4[2],mlM4[3],mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]);
    }

    puts("\nMaximum likelihood estimate considering each locus independently");
    puts("\nnid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaMl/nt\tthetaMl/region\tml\n");
    for(j=0;j<data[0].nloci;j++) {
        if(global[0][0].dataequalsim)
			printf("%ld:%s\t%.2f\t%g\t%d\t%d\t%g\t%g\t%g\n",data[0].numberloci[j],matrix[j].gene,data[0].factor_chrn[j],data[0].length[j],data[0].nsamA[j],data[0].SA[j],
				thetalocusml[j],thetalocusml[j]*(double)data[0].length[j],locusml[j]);
		else
			printf("%ld\t%.2f\t%g\t%d\t%d\t%g\t%g\t%g\n",data[0].numberloci[j],data[0].factor_chrn[j],data[0].length[j],data[0].nsamA[j],data[0].SA[j],
				thetalocusml[j],thetalocusml[j]*(double)data[0].length[j],locusml[j]);
	}
    printf("ML:\t\t\t\t\t\t\t%g",locusml[data[0].nloci]);
    puts("\n");

	puts("\nSUMMARY\n");
    puts("\nMaximum likelihood theta values estimated using different models:\n");
    printf("\nM1: %g (%d loci)",thetaM1,data[0].nloci);
    if(data[0].nloci > 2) 
        printf("\nM2: %g (%d loci), %g (%d loci)",thetaM2[0],nloci2m1,thetaM2[1],nloci2m2);
    if(data[0].nloci > 3) 
        printf("\nM3: %g (%d loci), %g (%d loci), %g (%d loci)",thetaM3[0],nloci3m1,thetaM3[1],nloci3m2,thetaM3[2],nloci3m3);
    if(data[0].nloci > 4) {
        printf("\nM4: %g (%d loci), %g (%d loci), %g (%d loci), %g (%d loci)",
        thetaM4[0],nloci4m1,thetaM4[1],nloci4m2,thetaM4[2],nloci4m3,thetaM4[3],nloci4m4);
    }
    puts("\n");

    puts("\nModels:\n");
    printf(" M1: One theta value for all loci.\tML = %g\n",mlM1);
    if(data[0].nloci > 2) {
        printf(" M2: Two independent theta values.\tML = %g\n",mlM2[0]+mlM2[1]);
        if(data[0].nloci > 3) {
            printf(" M3: Three independent theta values.\tML = %g\n",mlM3[0]+mlM3[1]+mlM3[2]);
            if(data[0].nloci > 4) {
                printf(" M4: Four independent theta values.\tML = %g\n",mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]);
            }
        }
    }
    if(data[0].nloci > 1) printf(" M%d: %d independent theta values.\tML = %g\n",data[0].nloci,data[0].nloci,locusml[data[0].nloci]);
/*    
    if(data[0].nloci > 1) {
        puts("\nLikelihood ratio test (LRT) assuming a Chi-square distribution:\n");
        if(data[0].nloci > 2)
            printf("\nLR(M1 vs M2) = %g with 2 degree of freedom. Probability = %g\n",LR12,probLRT12);
        if(data[0].nloci > 3)
            printf("LR(M1 vs M3) = %g with 4 degrees of freedom. Probability = %g\n",LR13,probLRT13);
        if(data[0].nloci > 4)
            printf("LR(M1 vs M4) = %g with 6 degrees of freedom. Probability = %g\n",LR14,probLRT14);
        printf("LR(M1 vs M%d) = %g with %d degrees of freedom. Probability = %g\n",data[0].nloci,LR1l,dgf1l,probLRT1l);
    }
    if(data[0].nloci > 2) {
        if(data[0].nloci > 3)
            printf("\nLR(M2 vs M3) = %g with 2 degree of freedom. Probability = %g\n",LR23,probLRT23);
        if(data[0].nloci > 4)
            printf("LR(M2 vs M4) = %g with 4 degrees of freedom. Probability = %g\n",LR24,probLRT24);
        if(data[0].nloci > 3) 
			printf("LR(M2 vs M%d) = %g with %d degrees of freedom. Probability = %g\n",data[0].nloci,LR2l,dgf2l,probLRT2l);
    }
    if(data[0].nloci > 3) {
        if(data[0].nloci > 4)
            printf("\nLR(M3 vs M4) = %g with 2 degree of freedom. Probability = %g\n",LR34,probLRT34);
        if(data[0].nloci > 5) 
			printf("LR(M3 vs M%d) = %g with %d degrees of freedom. Probability = %g\n",data[0].nloci,LR3l,dgf3l,probLRT3l);
    }
    if(data[0].nloci > 4)
        if(data[0].nloci > 7)
			printf("\nLR(M4 vs M%d) = %g with %d degrees of freedom. Probability = %g\n",data[0].nloci,LR4l,dgf4l,probLRT4l);
*/
	/******** AKAIKE INFORMATION CRITERION (AIC) **********/
    /*calculate AIC: the number of free parameters are not clear for me... I put the highest number...*/
	AIC1 = -(double)2. * (mlM1 + (double)1.);
	if(data[0].nloci > 2) AIC2 = -(double)2. * (mlM2[0] + mlM2[1] + (double)3.);
	if(data[0].nloci > 3) AIC3 = -(double)2. * (mlM3[0] + mlM3[1] + mlM3[2] + (double)5.);
	if(data[0].nloci > 4) AIC4 = -(double)2. * (mlM4[0] + mlM4[1] + mlM4[2] + mlM4[3] + (double)7.);
	if(data[0].nloci > 1) AICn = -(double)2. * (locusml[data[0].nloci] + (double)data[0].nloci);
	/*calculate minimum AIC value*/
	AICmin = AIC1; bestAIC = 1;
	if(data[0].nloci > 2) if(AIC2 < AICmin) {AICmin = AIC2; bestAIC = 2;}
	if(data[0].nloci > 3) if(AIC3 < AICmin) {AICmin = AIC3; bestAIC = 3;}
	if(data[0].nloci > 4) if(AIC4 < AICmin) {AICmin = AIC4; bestAIC = 4;}
	if(data[0].nloci > 1) if(AICn < AICmin) {AICmin = AICn; bestAIC = data[0].nloci;}
	/*calculate sumwr*/
	sumwAIC = (double)exp(-0.5*(AIC1-AICmin));
	if(data[0].nloci > 2) sumwAIC += (double)exp(-0.5*(AIC2-AICmin));
	if(data[0].nloci > 3) sumwAIC += (double)exp(-0.5*(AIC3-AICmin));
	if(data[0].nloci > 4) sumwAIC += (double)exp(-0.5*(AIC4-AICmin));
	if(data[0].nloci > 1) sumwAIC += (double)exp(-0.5*(AICn-AICmin));	
	/*print*/
	/*
	printf("\nAKAIKE INFORMATION CRITERION (AIC):\n");
	printf("Warning: These results are orientative. Statistical inference can obtained by Monte Carlo simulations\n");
	printf("\nModel\tlnL\tK\tAIC\tdelta\tweight\n");
	printf("M1\t%g\t%d\t%g\t%g\t%g\n",mlM1,1,AIC1,AIC1-AICmin,(double)exp(-0.5*(AIC1-AICmin))/sumwAIC);
	if(data[0].nloci > 2) printf("M2\t%g\t%d\t%g\t%g\t%g\n",mlM2[0]+mlM2[1],3,AIC2,AIC2-AICmin,(double)exp(-0.5*(AIC2-AICmin))/sumwAIC);
	if(data[0].nloci > 3) printf("M3\t%g\t%d\t%g\t%g\t%g\n",mlM3[0]+mlM3[1]+mlM3[2],5,AIC3,AIC3-AICmin,(double)exp(-0.5*(AIC3-AICmin))/sumwAIC);
	if(data[0].nloci > 4) printf("M4\t%g\t%d\t%g\t%g\t%g\n",mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3],7,AIC4,AIC4-AICmin,(double)exp(-0.5*(AIC4-AICmin))/sumwAIC);
	if(data[0].nloci > 1) printf("M%d\t%g\t%d\t%g\t%g\t%g\n",data[0].nloci,locusml[data[0].nloci],data[0].nloci,AICn,AICn-AICmin,(double)exp(-0.5*(AICn-AICmin))/sumwAIC);
	
	if(bestAIC == 1) printf("\nThe selected model is M1\n");
	else if(bestAIC == 2) printf("\nThe selected model is M2\n");
		else if(bestAIC == 3) printf("\nThe selected model is M3\n");
			else if(bestAIC == 4) printf("\nThe selected model is M4\n");
				else if(bestAIC == data[0].nloci) printf("\nThe selected model is M%d\n",data[0].nloci);
	printf("Supported models should be those with delta values lower than 2.\n");
	printf("Not supported models should be those with delta values larger than 10.\n");
	*/
	/***********************************************/
    fflush(stdout);
    
	/*IN THE OUTPUT FILE*/
	
	if(file_output) {
		fputs("\nESTIMATION BY MAXIMUM LIKELIHOOD OF THETA GIVEN THE NUMBER OF SEGREGATING SITES.\n",file_output);
		#if EQTAVARE == 0
		fputs("\nCalculation using the recursive equation n. 12 from Hudson 1990. Oxford surveys in evol. biol. 1-44\n",
			file_output);
		fputs("\n     (Precission until 5e-12 for each locus).\n",file_output);
		#else
		fputs("\nCalculation using the equation 9.5 from Tavaré Theor. Pop. Biol. 26: 119-164\n",file_output);
		#endif
		fprintf(file_output,"\ntheta/nt(min): %g\ntheta/nt(max): %g\n\n",data[0].thetamin_user,data[0].thetamax_user);
		fputs("\nMaximum likelihood estimate considering one theta:\n",file_output);
		fprintf(file_output,"\nnumber of loci: %d",data[0].nloci);
		fputs("\nnid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM1/region\n",file_output);
		for(j=0;j<data[0].nloci;j++) {
			if(global[0][0].dataequalsim) 
				fprintf(file_output,"%ld:%s\t%.2f\t%.2f\t%d\t%d\t%g\n",data[0].numberloci[j],matrix[j].gene,data[0].factor_chrn[j],data[0].length[j],data[0].nsamA[j],data[0].SA[j],thetaM1*(double)data[0].length[j]);
			else 
				fprintf(file_output,"%ld\t%.2f\t%.2f\t%d\t%d\t%g\n",data[0].numberloci[j],data[0].factor_chrn[j],data[0].length[j],data[0].nsamA[j],data[0].SA[j],thetaM1*(double)data[0].length[j]);
		}
		fprintf(file_output,"\ntheta/nt\n %g",thetaM1);
		fprintf(file_output,"\nML\n %g\n",mlM1);
		if(data[0].nloci > 2) {
			fputs("\nMaximum likelihood estimate considering two groups of independent loci",file_output);
			j=0;
			while(vector2ml1[j] > -1) j++;
			fprintf(file_output,"\nlocus_group1: (%d)\n",j);
			
			fputs("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM2/region\n",file_output);
			j = 0;
			while(vector2ml1[j] > -1) {
				if(global[0][0].dataequalsim) 
					fprintf(file_output,"%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector2ml1[j]],matrix[vector2ml1[j]].gene,data[0].factor_chrn[vector2ml1[j]],data[0].length[vector2ml1[j]],
						data[0].nsamA[vector2ml1[j]],data[0].SA[vector2ml1[j]],thetaM2[0]*(double)data[0].length[vector2ml1[j]]);
				else 
					fprintf(file_output,"%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector2ml1[j]],data[0].factor_chrn[vector2ml1[j]],data[0].length[vector2ml1[j]],
						data[0].nsamA[vector2ml1[j]],data[0].SA[vector2ml1[j]],thetaM2[0]*(double)data[0].length[vector2ml1[j]]);
				j++;
			}
			fprintf(file_output,"\nlocus_group2: (%d)\n",data[0].nloci-j);
			fputs("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM2/region\n",file_output);
			k = 0;
			while(vector2ml2[k] > -1) {
				if(global[0][0].dataequalsim)
					fprintf(file_output,"%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector2ml2[k]],matrix[vector2ml2[k]].gene,data[0].factor_chrn[vector2ml2[k]],data[0].length[vector2ml2[k]],
						data[0].nsamA[vector2ml2[k]],data[0].SA[vector2ml2[k]],thetaM2[1]*(double)data[0].length[vector2ml2[k]]);
				else 
					fprintf(file_output,"%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector2ml2[k]],data[0].factor_chrn[vector2ml2[k]],data[0].length[vector2ml2[k]],
						data[0].nsamA[vector2ml2[k]],data[0].SA[vector2ml2[k]],thetaM2[1]*(double)data[0].length[vector2ml2[k]]);
				k++;
			}
			fprintf(file_output,"\ntheta1ml/nt\ttheta2ml/nt\n%g\t%g\n",thetaM2[0],thetaM2[1]);
			fprintf(file_output,"\nmlgroup1\tmlgroup2\tML\n%g\t%g\t%g\n",mlM2[0],mlM2[1],mlM2[0]+mlM2[1]);
		}
		if(data[0].nloci > 3) {
			fputs("\nMaximum likelihood estimate considering three groups of independent loci",file_output);
			j=0;
			while(vector3ml1[j] > -1) j++;
			fprintf(file_output,"\nlocus_group1: (%d)\n",j);
			fputs("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM3/region\n",file_output);
			j = 0;
			while(vector3ml1[j] > -1) {
				if(global[0][0].dataequalsim)
					fprintf(file_output,"%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector3ml1[j]],matrix[vector3ml1[j]].gene,data[0].factor_chrn[vector3ml1[j]],data[0].length[vector3ml1[j]],
						data[0].nsamA[vector3ml1[j]],data[0].SA[vector3ml1[j]],thetaM3[0]*(double)data[0].length[vector3ml1[j]]);
				else
					fprintf(file_output,"%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector3ml1[j]],data[0].factor_chrn[vector3ml1[j]],data[0].length[vector3ml1[j]],
						data[0].nsamA[vector3ml1[j]],data[0].SA[vector3ml1[j]],thetaM3[0]*(double)data[0].length[vector3ml1[j]]);
				j++;
			}
			k=0;
			while(vector3ml2[k] > -1) k++;
			fprintf(file_output,"\nlocus_group2: (%d)\n",k);
			fputs("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM3/region\n",file_output);
			k = 0;
			while(vector3ml2[k] > -1) {
				if(global[0][0].dataequalsim) 
					fprintf(file_output,"%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector3ml2[k]],matrix[vector3ml2[k]].gene,data[0].factor_chrn[vector3ml2[k]],data[0].length[vector3ml2[k]],
						data[0].nsamA[vector3ml2[k]],data[0].SA[vector3ml2[k]],thetaM3[1]*(double)data[0].length[vector3ml2[k]]);
				else
					fprintf(file_output,"%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector3ml2[k]],data[0].factor_chrn[vector3ml2[k]],data[0].length[vector3ml2[k]],
						data[0].nsamA[vector3ml2[k]],data[0].SA[vector3ml2[k]],thetaM3[1]*(double)data[0].length[vector3ml2[k]]);
				k++;
			}
			fprintf(file_output,"\nlocus_group3: (%d)\n",data[0].nloci-j-k);
			fputs("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM3/region\n",file_output);
			l=0;
			while(vector3ml3[l] > -1) {
				if(global[0][0].dataequalsim) 
					fprintf(file_output,"%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector3ml3[l]],matrix[vector3ml3[l]].gene,data[0].factor_chrn[vector3ml3[l]],data[0].length[vector3ml3[l]],
						data[0].nsamA[vector3ml3[l]],data[0].SA[vector3ml3[l]],thetaM3[2]*(double)data[0].length[vector3ml3[l]]);
				else
					fprintf(file_output,"%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector3ml3[l]],data[0].factor_chrn[vector3ml3[l]],data[0].length[vector3ml3[l]],
						data[0].nsamA[vector3ml3[l]],data[0].SA[vector3ml3[l]],thetaM3[2]*(double)data[0].length[vector3ml3[l]]);
				l++;
			}
			fprintf(file_output,"\ntheta1ml/nt\ttheta2ml/nt\ttheta3ml/nt\n%g\t%g\t%g\n",thetaM3[0],thetaM3[1],thetaM3[2]);
			fprintf(file_output,"\nmlgroup1\tmlgroup2\tmlgroup3\tML\n%g\t%g\t%g\t%g\n",mlM3[0],mlM3[1],mlM3[2],
																						mlM3[0]+mlM3[1]+mlM3[2]);
		}
		if(data[0].nloci > 4) {
			fputs("\nMaximum likelihood estimate considering four groups of independent loci",file_output);
			j=0;
			while(vector4ml1[j] > -1) j++;
			fprintf(file_output,"\nlocus_group1: (%d)\n",j);
			fputs("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM4/region\n",file_output);
			j = 0;
			while(vector4ml1[j] > -1) {
				if(global[0][0].dataequalsim)
					fprintf(file_output,"%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml1[j]],matrix[vector4ml1[j]].gene,data[0].factor_chrn[vector4ml1[j]],data[0].length[vector4ml1[j]],
						data[0].nsamA[vector4ml1[j]],data[0].SA[vector4ml1[j]],thetaM4[0]*(double)data[0].length[vector4ml1[j]]);
				else
					fprintf(file_output,"%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml1[j]],data[0].factor_chrn[vector4ml1[j]],data[0].length[vector4ml1[j]],
						data[0].nsamA[vector4ml1[j]],data[0].SA[vector4ml1[j]],thetaM4[0]*(double)data[0].length[vector4ml1[j]]);
				j++;
			}
			k=0;
			while(vector4ml2[k] > -1) k++;
			fprintf(file_output,"\nlocus_group2: (%d)\n",k);
			fputs("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM4/region\n",file_output);
			k = 0;
			while(vector4ml2[k] > -1) {
				if(global[0][0].dataequalsim)
					fprintf(file_output,"%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml2[k]],matrix[vector4ml2[k]].gene,data[0].factor_chrn[vector4ml2[k]],data[0].length[vector4ml2[k]],
						data[0].nsamA[vector4ml2[k]],data[0].SA[vector4ml2[k]],thetaM4[1]*(double)data[0].length[vector4ml2[k]]);
				else
					fprintf(file_output,"%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml2[k]],data[0].factor_chrn[vector4ml2[k]],data[0].length[vector4ml2[k]],
						data[0].nsamA[vector4ml2[k]],data[0].SA[vector4ml2[k]],thetaM4[1]*(double)data[0].length[vector4ml2[k]]);
				k++;
			}
			l=0;
			while(vector4ml3[l] > -1) l++;
			fprintf(file_output,"\nlocus_group3: (%d)\n",l);
			fputs("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM4/region\n",file_output);
			l = 0;
			while(vector4ml3[l] > -1) {
				if(global[0][0].dataequalsim)
					fprintf(file_output,"%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml3[l]],matrix[vector4ml3[l]].gene,data[0].factor_chrn[vector4ml3[l]],data[0].length[vector4ml3[l]],
						data[0].nsamA[vector4ml3[l]],data[0].SA[vector4ml3[l]],thetaM4[2]*(double)data[0].length[vector4ml3[l]]);
				else
					fprintf(file_output,"%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml3[l]],data[0].factor_chrn[vector4ml3[l]],data[0].length[vector4ml3[l]],
						data[0].nsamA[vector4ml3[l]],data[0].SA[vector4ml3[l]],thetaM4[2]*(double)data[0].length[vector4ml3[l]]);
				l++;
			}
			fprintf(file_output,"\nlocus_group4: (%d)\n",data[0].nloci-j-k-l);
			fputs("nid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaM4/region\n",file_output);
			m = 0;
			while(vector4ml4[m] > -1) {
				if(global[0][0].dataequalsim)
					fprintf(file_output,"%ld:%s\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml4[m]],matrix[vector4ml4[m]].gene,data[0].factor_chrn[vector4ml4[m]],data[0].length[vector4ml4[m]],
						data[0].nsamA[vector4ml4[m]],data[0].SA[vector4ml4[m]], thetaM4[3]*(double)data[0].length[vector4ml4[m]]);
				else
					fprintf(file_output,"%ld\t%.2f\t%g\t%d\t%d\t%g\n",data[0].numberloci[vector4ml4[m]],data[0].factor_chrn[vector4ml4[m]],data[0].length[vector4ml4[m]],
						data[0].nsamA[vector4ml4[m]],data[0].SA[vector4ml4[m]], thetaM4[3]*(double)data[0].length[vector4ml4[m]]);
				m++;
			}
			fprintf(file_output,"\ntheta1ml/nt\ttheta2ml/nt\ttheta3ml/nt\ttheta4ml/nt\n%g\t%g\t%g\t%g\n",thetaM4[0],
																				thetaM4[1],thetaM4[2],thetaM4[3]);
			fprintf(file_output,"\nmlgroup1\tmlgroup2\tmlgroup3\tmlgroup4\tML\n%g\t%g\t%g\t%g\t%g\n",
																					mlM4[0],mlM4[1],mlM4[2],mlM4[3],
																					mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]);
		}
		fputs("\nMaximum likelihood estimate considering each locus independently",file_output);
		fputs("\nnid_locus\tfactor_chrn\tlength\tnsam\tS\tthetaMl/nt\tthetaMl/region\tml\n",file_output);
		for(j=0;j<data[0].nloci;j++) {
			if(global[0][0].dataequalsim)
				fprintf(file_output,"%ld:%s\t%.2f\t%g\t%d\t%d\t%g\t%g\t%g\n",data[0].numberloci[j],matrix[j].gene,data[0].factor_chrn[j],data[0].length[j],data[0].nsamA[j],data[0].SA[j],
					thetalocusml[j],thetalocusml[j]*(double)data[0].length[j],locusml[j]);
			else
				fprintf(file_output,"%ld\t%.2f\t%g\t%d\t%d\t%g\t%g\t%g\n",data[0].numberloci[j],data[0].factor_chrn[j],data[0].length[j],data[0].nsamA[j],data[0].SA[j],
					thetalocusml[j],thetalocusml[j]*(double)data[0].length[j],locusml[j]);
		}
		fprintf(file_output,"ML:\t\t\t\t\t\t%g",locusml[data[0].nloci]);
		fputs("\n",file_output);

		
		fputs("\nAll iterations for the theta/nt and the likelihood considering all loci\n",file_output);
		fputs("theta/nt\t",file_output);
		fputs("lik[total]\t",file_output);
		for(j=0;j<data[0].nloci;j++) fprintf(file_output,"lik[%d]\t",j);
		fputs("\n",file_output);
		for(i=0;i<data[0].steps_user;i++) {
			fprintf(file_output,"%g\t",((double)i/((double)data[0].steps_user-(double)1))*
								  (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user);
			fprintf(file_output,"%g\t",ml[j][i]);
			for(j=0;j<data[0].nloci;j++) fprintf(file_output,"%g\t",ml[j][i]);
			fputs("\n",file_output);
		}
		fputs("\n",file_output);

		fputs("\nSUMMARY\n",file_output);
		fputs("\nMaximum likelihood theta values estimated using different models:\n",file_output);
		fprintf(file_output,"\nM1: %g (%d loci)",thetaM1,data[0].nloci);
		if(data[0].nloci > 2) 
			fprintf(file_output,"\nM2: %g (%d loci), %g (%d loci)",thetaM2[0],nloci2m1,thetaM2[1],nloci2m2);
		if(data[0].nloci > 3) {
			fprintf(file_output,"\nM3: %g (%d loci), %g (%d loci), %g (%d loci)",
			thetaM3[0],nloci3m1,thetaM3[1],nloci3m2,thetaM3[2],nloci3m3);
		}
		if(data[0].nloci > 4) {
			fprintf(file_output,"\nM4: %g (%d loci), %g (%d loci), %g (%d loci), %g (%d loci)",
			thetaM4[0],nloci4m1,thetaM4[1],nloci4m2,thetaM4[2],nloci4m3,thetaM4[3],nloci4m4);
		}
		fputs("\n",file_output);

		fputs("\nModels:\n\n",file_output);
		fprintf(file_output," M1: One theta value for all loci.\tML = %g\n",mlM1);
		if(data[0].nloci > 2) {
			fprintf(file_output," M2: Two independent theta values.\tML = %g\n",mlM2[0]+mlM2[1]);
			if(data[0].nloci > 3) {
				fprintf(file_output," M3: Three independent theta values.\tML = %g\n",mlM3[0]+mlM3[1]+mlM3[2]);
				if(data[0].nloci > 4) {
					fprintf(file_output," M4: Four independent theta values.\tML = %g\n",mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]);
				}
			}
		}
		if(data[0].nloci > 1) fprintf(file_output," M%d: %d independent theta values.\tML = %g\n",data[0].nloci,data[0].nloci,
				locusml[data[0].nloci]);
		/*
		if(data[0].nloci > 1) {
			fputs("\nLikelihood ratio test (LRT) assuming a Chi-square distribution:\n",file_output);
			if(data[0].nloci > 2)
				fprintf(file_output,"\nLR(M1 vs M2) = %g with 2 degree of freedom. Probability = %g\n",LR12,probLRT12);
			if(data[0].nloci > 3)
				fprintf(file_output,"LR(M1 vs M3) = %g with 4 degrees of freedom. Probability = %g\n",LR13,probLRT13);
			if(data[0].nloci > 4)
				fprintf(file_output,"LR(M1 vs M4) = %g with 6 degrees of freedom. Probability = %g\n",LR14,probLRT14);
			fprintf(file_output,"LR(M1 vs M%d) = %g with %d degrees of freedom. Probability = %g\n",data[0].nloci,
					LR1l,dgf1l,probLRT1l);
		}
		if(data[0].nloci > 2) {
			if(data[0].nloci > 3)
				fprintf(file_output,"\nLR(M2 vs M3) = %g with 2 degree of freedom. Probability = %g\n",LR23,probLRT23);
			if(data[0].nloci > 4)
				fprintf(file_output,"LR(M2 vs M4) = %g with 4 degrees of freedom. Probability = %g\n",LR24,probLRT24);
			if(data[0].nloci > 3)
				fprintf(file_output,"LR(M2 vs M%d) = %g with %d degrees of freedom. Probability = %g\n",data[0].nloci,
				LR2l,dgf2l,probLRT2l);
		}
		if(data[0].nloci > 3) {
			if(data[0].nloci > 4)
				fprintf(file_output,"\nLR(M3 vs M4) = %g with 2 degree of freedom. Probability = %g\n",LR34,probLRT34);
			if(data[0].nloci > 5)
				fprintf(file_output,"LR(M3 vs M%d) = %g with %d degrees of freedom. Probability = %g\n",data[0].nloci,
				LR3l,dgf3l,probLRT3l);
		}
		if(data[0].nloci > 4)
			if(data[0].nloci > 7)	
				fprintf(file_output,"\nLR(M4 vs M%d) = %g with %d degrees of freedom. Probability = %g\n",data[0].nloci,
				LR4l,dgf4l,probLRT4l);
		*/
		    
		/******** AKAIKE INFORMATION CRITERION (AIC) **********/
		/*calculate AIC: the number of free parameters are not clear for me... I put the highest number...*/
		AIC1 = -(double)2. * (mlM1 + (double)1.);
		if(data[0].nloci > 2) AIC2 = -(double)2. * (mlM2[0] + mlM2[1] + (double)3.);
		if(data[0].nloci > 3) AIC3 = -(double)2. * (mlM3[0] + mlM3[1] + mlM3[2] + (double)5.);
		if(data[0].nloci > 4) AIC4 = -(double)2. * (mlM4[0] + mlM4[1] + mlM4[2] + mlM4[3] + (double)7.);
		if(data[0].nloci > 1) AICn = -(double)2. * (locusml[data[0].nloci] + (double)data[0].nloci);
		/*calculate minimum AIC value*/
		AICmin = AIC1; bestAIC = 1;
		if(data[0].nloci > 2) if(AIC2 < AICmin) {AICmin = AIC2; bestAIC = 2;}
		if(data[0].nloci > 3) if(AIC3 < AICmin) {AICmin = AIC3; bestAIC = 3;}
		if(data[0].nloci > 4) if(AIC4 < AICmin) {AICmin = AIC4; bestAIC = 4;}
		if(data[0].nloci > 1) if(AICn < AICmin) {AICmin = AICn; bestAIC = data[0].nloci;}
		/*calculate sumwr*/
		sumwAIC = (double)exp(-0.5*(AIC1-AICmin));
		if(data[0].nloci > 2) sumwAIC += (double)exp(-0.5*(AIC2-AICmin));
		if(data[0].nloci > 3) sumwAIC += (double)exp(-0.5*(AIC3-AICmin));
		if(data[0].nloci > 4) sumwAIC += (double)exp(-0.5*(AIC4-AICmin));
		if(data[0].nloci > 1) sumwAIC += (double)exp(-0.5*(AICn-AICmin));	
		/*print*/
		/*
		fprintf(file_output,"\nAKAIKE INFORMATION CRITERION (AIC):\n");
		fprintf(file_output,"Warning: These results are orientative. Statistical inference is obtained by Monte Carlo simulations\n");
		fprintf(file_output,"\nModel\tlnL\tK\tAIC\tdelta\tweight\n");
		fprintf(file_output,"M1\t%g\t%d\t%g\t%g\t%g\n",mlM1,1,AIC1,AIC1-AICmin,(double)exp(-0.5*(AIC1-AICmin))/sumwAIC);
		if(data[0].nloci > 2) fprintf(file_output,"M2\t%g\t%d\t%g\t%g\t%g\n",mlM2[0]+mlM2[1],3,AIC2,AIC2-AICmin,(double)exp(-0.5*(AIC2-AICmin))/sumwAIC);
		if(data[0].nloci > 3) fprintf(file_output,"M3\t%g\t%d\t%g\t%g\t%g\n",mlM3[0]+mlM3[1]+mlM3[2],5,AIC3,AIC3-AICmin,(double)exp(-0.5*(AIC3-AICmin))/sumwAIC);
		if(data[0].nloci > 4) fprintf(file_output,"M4\t%g\t%d\t%g\t%g\t%g\n",mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3],7,AIC4,AIC4-AICmin,(double)exp(-0.5*(AIC4-AICmin))/sumwAIC);
		if(data[0].nloci > 1) fprintf(file_output,"M%d\t%g\t%d\t%g\t%g\t%g\n",data[0].nloci,locusml[data[0].nloci],data[0].nloci,AICn,AICn-AICmin,(double)exp(-0.5*(AICn-AICmin))/sumwAIC);
		
		if(bestAIC == 1) fprintf(file_output,"\nThe selected model is M1\n");
		else if(bestAIC == 2) fprintf(file_output,"\nThe selected model is M2\n");
			else if(bestAIC == 3) fprintf(file_output,"\nThe selected model is M3\n");
				else if(bestAIC == 4) fprintf(file_output,"\nThe selected model is M4\n");
					else if(bestAIC == data[0].nloci) fprintf(file_output,"\nThe selected model is M%d\n",data[0].nloci);
		fprintf(file_output,"Supported models should be those with delta values lower than 2.\n");
		fprintf(file_output,"Not supported models should be those with delta values larger than 10.\n");
		*/
		/***********************************************/
		
		fflush(file_output);
		/*END PRINTING DATA RESULTS WITHOUT COALESCENT SIMULATIONS*/
    }
    /*CALCULATE PROBABILIY WITH COALESCENT SIMULATIONS*/
    if(data[0].nloci > 1) {
		printf("\nSTATISTICAL INFERENCE OBTAINED BY MONTE CARLO COALESCENT SIMULATIONS");
		if(file_output) fprintf(file_output,"\nSTATISTICAL INFERENCE OBTAINED BY MONTE CARLO COALESCENT  SIMULATIONS");
		
		init_seed1(data[0].seed);
		if(probcoalthhaetasml(data,thetaM1,thetaM2,thetaM3,thetaM4,vector2ml1,vector2ml2,vector3ml1,vector3ml2,vector3ml3,
								vector4ml1,vector4ml2,vector4ml3,vector4ml4,nloci2m1,nloci2m2,nloci3m1,nloci3m2,nloci3m3,
								nloci4m1,nloci4m2,nloci4m3,nloci4m4,LR12,LR13,LR14,LR1l,LR23,LR24,LR2l,LR34,LR3l,LR4l,
								&Pcoal,lrdist/*,mthetasim*/)) return(1);    
    }
    /*Assign results to datams*/
	if(data[0].nloci > 2) {
		datams[0][0].LR12 = LR12;
		datams[0][0].p1 = Pcoal.probc12;
	}
	if(data[0].nloci == 2) {
		datams[0][0].LR12 = LR1l;
		datams[0][0].p1 = Pcoal.probc1l;
	}
	if(data[0].nloci > 3) {
		datams[0][0].LR23 = LR23;
		datams[0][0].p2 = Pcoal.probc23;
	}
	if(data[0].nloci == 3) {
		datams[0][0].LR23 = LR2l;
		datams[0][0].p2 = Pcoal.probc2l;
	}
	if(data[0].nloci > 4) {
		datams[0][0].LR34 = LR34;
		datams[0][0].p3 = Pcoal.probc34;
	}
	if(data[0].nloci == 4) {
		datams[0][0].LR34 = LR3l;
		datams[0][0].p3 = Pcoal.probc3l;
	}
	datams[0][0].LR4l = LR4l;
	datams[0][0].p4 = Pcoal.probc4l;
	
	if(data[0].nloci == 1) datams[0][0].bestmodel = 1;
	if(data[0].nloci == 2) {
		if(Pcoal.probc1l > (double)0.05) datams[0][0].bestmodel = 1;
		else datams[0][0].bestmodel = 2;
	}
	if(data[0].nloci > 2) {
		if(Pcoal.probc12 > (double)0.05) datams[0][0].bestmodel = 1;
		else {
			if(data[0].nloci == 3) {
				if(Pcoal.probc2l > (double)0.05) datams[0][0].bestmodel = 2;
				else datams[0][0].bestmodel = 3;
			}
			if(data[0].nloci > 3) {
				if(Pcoal.probc23 > (double)0.05) datams[0][0].bestmodel = 2;
				else {
					if(data[0].nloci == 4) {
						if(Pcoal.probc3l > (double)0.05) datams[0][0].bestmodel = 3;
						else datams[0][0].bestmodel = 4;
					}
					if(data[0].nloci > 4) {
						if(Pcoal.probc34 > (double)0.05) datams[0][0].bestmodel = 3;
						else {
							if(Pcoal.probc4l > (double)0.05) datams[0][0].bestmodel = 4;
							else datams[0][0].bestmodel = 5;
						}
					}
				}
			}
		}
	}
	global[0][0].mltheta = datams[0][0].bestmodel;
	
	/*PRINT OUTPUT COALESCENT SIMULATIONS*/
	if(data[0].nloci > 1) {
        printf("\nLikelihood ratio test (LRT) probabilities obtained with %ld coalescent simulations:\n",data[0].ncoalsim);
        if(data[0].nloci > 2) {
            if(Pcoal.probc12 <= (double)1.)
                printf("\nLR(M1 vs M2) = %g. Probability = %g\n",LR12,Pcoal.probc12);
            else
                printf("\nLR(M1 vs M2) = %g. Probability not calculated\n",LR12);
		}
		/*
		if(data[0].nloci > 3)
            if(Pcoal.probc13 <= (double)1.)
                printf("LR(M1 vs M3) = %g. Probability = %g\n",LR13,Pcoal.probc13);
            else
                printf("LR(M1 vs M3) = %g. Probability not calculated\n",LR13);
        if(data[0].nloci > 4)
            if(Pcoal.probc14 <= (double)1.)
                printf("LR(M1 vs M4) = %g. Probability = %g\n",LR14,Pcoal.probc14);
            else
                printf("LR(M1 vs M4) = %g. Probability not calculated\n",LR14);
        if(Pcoal.probc1l <= (double)1.)
            printf("LR(M1 vs M%d) = %g. Probability = %g\n",data[0].nloci,LR1l,Pcoal.probc1l);
        else
            printf("LR(M1 vs M%d) = %g. Probability not calculated\n",data[0].nloci,LR1l);
		*/
        if(data[0].nloci == 2) {
            if(Pcoal.probc1l <= (double)1.)
                printf("\nLR(M1 vs M2) = %g. Probability = %g\n",LR1l,Pcoal.probc1l);
            else
                printf("\nLR(M1 vs M2) = %g. Probability not calculated\n",LR1l);
		}
    }
    if(data[0].nloci > 2) {
        if(data[0].nloci > 3) {
            if(Pcoal.probc23 <= (double)1.)
                printf("LR(M2 vs M3) = %g. Probability = %g\n",LR23,Pcoal.probc23);
            else
                printf("LR(M2 vs M3) = %g. Probability not calculated\n",LR23);
        }
		/*
        if(data[0].nloci > 4) {
            if(Pcoal.probc24 <= (double)1.)
                printf("LR(M2 vs M4) = %g. Probability = %g\n",LR24,Pcoal.probc24);
            else
                printf("LR(M2 vs M4) = %g. Probability not calculated\n",LR24);
        }
        if(Pcoal.probc2l <= (double)1.)
            printf("LR(M2 vs M%d) = %g. Probability = %g\n",data[0].nloci,
                LR2l,Pcoal.probc2l);
        else
            printf("LR(M2 vs M%d) = %g. Probability not calculated\n",data[0].nloci,LR2l);
		*/
        if(data[0].nloci == 3) {
            if(Pcoal.probc2l <= (double)1.)
                printf("LR(M2 vs M3) = %g. Probability = %g\n",LR2l,Pcoal.probc2l);
            else
                printf("LR(M2 vs M3) = %g. Probability not calculated\n",LR2l);
        }
    }
    if(data[0].nloci > 3) {
        if(data[0].nloci > 4) {
            if(Pcoal.probc34 <= (double)1.)
                printf("LR(M3 vs M4) = %g. Probability = %g\n",LR34,Pcoal.probc34);
            else
                printf("LR(M3 vs M4) = %g. Probability not calculated\n",LR34);
        }
		/*
        if(Pcoal.probc3l <= (double)1.)
            printf("LR(M3 vs M%d) = %g. Probability = %g\n",data[0].nloci,
                LR3l,Pcoal.probc3l);
        else
            printf("LR(M3 vs M%d) = %g. Probability not calculated\n",data[0].nloci,LR3l);
		*/
        if(data[0].nloci == 4) {
            if(Pcoal.probc3l <= (double)1.)
                printf("LR(M3 vs M4) = %g. Probability = %g\n",LR3l,Pcoal.probc3l);
            else
                printf("LR(M3 vs M4) = %g. Probability not calculated\n",LR3l);
        }
    }
    if(data[0].nloci > 4) {
        if(Pcoal.probc4l <= (double)1.)
            printf("LR(M4 vs M%d) = %g. Probability = %g\n",data[0].nloci,
                LR4l,Pcoal.probc4l);
        else
            printf("LR(M4 vs M%d) = %g. Probability not calculated\n",data[0].nloci,LR4l);
    }

	/*THE BEST MODEL*/
	if(datams[0][0].bestmodel == 1) printf("\nTHE BEST MODEL is M%d. (A single theta/nt per all loci).\n\n",datams[0][0].bestmodel);
	else {
		if(datams[0][0].bestmodel > 1 && datams[0][0].bestmodel < 5) 
			printf("\nTHE BEST MODEL is M%d. (%d different theta/nt distributed in the loci).\n\n",datams[0][0].bestmodel,datams[0][0].bestmodel);
		if(datams[0][0].bestmodel == 5) printf("\nTHE BEST MODEL is M%d. (Each loci has a different theta/nt).\n\n",datams[0][0].n_loci);
	}

    /*PRINT OUTPUT TO FILE FOR COALESCENT SIMULATIONS*/
	if(file_output) {
		if(data[0].nloci > 1) {
			fprintf(file_output,"\nLikelihood ratio test (LRT) probabilities obtained with %ld coalescent simulations:\n",
				data[0].ncoalsim);
			if(data[0].nloci > 2) {
				if(Pcoal.probc12 <= (double)1.)
					fprintf(file_output,"\nLR(M1 vs M2) = %g. Probability = %g\n",LR12,Pcoal.probc12);
				else
					fprintf(file_output,"\nLR(M1 vs M2) = %g. Probability not calculated\n",LR12);
			}
			/*
			if(data[0].nloci > 3)
				if(Pcoal.probc13 <= (double)1.)
					fprintf(file_output,"LR(M1 vs M3) = %g. Probability = %g\n",LR13,Pcoal.probc13);
				else
					fprintf(file_output,"LR(M1 vs M3) = %g. Probability not calculated\n",LR13);
			if(data[0].nloci > 4)
				if(Pcoal.probc14 <= (double)1.)
					fprintf(file_output,"LR(M1 vs M4) = %g. Probability = %g\n",LR14,Pcoal.probc14);
				else
					fprintf(file_output,"LR(M1 vs M4) = %g. Probability not calculated\n",LR14);
			if(Pcoal.probc1l <= (double)1.)
				fprintf(file_output,"LR(M1 vs M%d) = %g. Probability = %g\n",data[0].nloci,LR1l,Pcoal.probc1l);
			else
				fprintf(file_output,"LR(M1 vs M%d) = %g. Probability not calculated\n",data[0].nloci,LR1l);
			*/
			if(data[0].nloci == 2) {
				if(Pcoal.probc1l <= (double)1.)
					fprintf(file_output,"\nLR(M1 vs M2) = %g. Probability = %g\n",LR1l,Pcoal.probc1l);
				else
					fprintf(file_output,"\nLR(M1 vs M2) = %g. Probability not calculated\n",LR1l);
			}
		}
		if(data[0].nloci > 2) {
			if(data[0].nloci > 3) {
				if(Pcoal.probc23 <= (double)1.)
					fprintf(file_output,"LR(M2 vs M3) = %g. Probability = %g\n",LR23,Pcoal.probc23);
				else
					fprintf(file_output,"LR(M2 vs M3) = %g. Probability not calculated\n",LR23);
			}
			/*
			if(data[0].nloci > 4) {
				if(Pcoal.probc24 <= (double)1.)
					fprintf(file_output,"LR(M2 vs M4) = %g. Probability = %g\n",LR24,Pcoal.probc24);
				else
					fprintf(file_output,"LR(M2 vs M4) = %g. Probability not calculated\n",LR24);
			}
			if(Pcoal.probc2l <= (double)1.)
				fprintf(file_output,"LR(M2 vs M%d) = %g. Probability = %g\n",data[0].nloci,
						LR2l,Pcoal.probc2l);
			else
				fprintf(file_output,"LR(M2 vs M%d) = %g. Probability not calculated\n",data[0].nloci,
						LR2l);
			*/
			if(data[0].nloci == 3) {
				if(Pcoal.probc2l <= (double)1.)
					fprintf(file_output,"LR(M2 vs M3) = %g. Probability = %g\n",LR2l,Pcoal.probc2l);
				else
					fprintf(file_output,"LR(M2 vs M3) = %g. Probability not calculated\n",LR2l);
			}
		}
		if(data[0].nloci > 3) {
			if(data[0].nloci > 4) {
				if(Pcoal.probc34 <= (double)1.)
					fprintf(file_output,"LR(M3 vs M4) = %g. Probability = %g\n",LR34,Pcoal.probc34);
				else
					fprintf(file_output,"LR(M3 vs M4) = %g. Probability not calculated\n",LR34);
			}
			/*
			if(Pcoal.probc3l <= (double)1.)
				fprintf(file_output,"LR(M3 vs M%d) = %g. Probability = %g\n",data[0].nloci,
						LR3l,Pcoal.probc3l);
			else
				fprintf(file_output,"LR(M3 vs M%d) = %g. Probability not calculated\n",data[0].nloci,
						LR3l);
			*/
			if(data[0].nloci == 4) {
				if(Pcoal.probc3l <= (double)1.)
					fprintf(file_output,"LR(M3 vs M4) = %g. Probability = %g\n",LR3l,Pcoal.probc3l);
				else
					fprintf(file_output,"LR(M3 vs M4) = %g. Probability not calculated\n",LR3l);
			}
		}
		if(data[0].nloci > 4) {
			if(Pcoal.probc4l <= (double)1.)
				fprintf(file_output,"LR(M4 vs M%d) = %g. Probability = %g\n",data[0].nloci,
					LR4l,Pcoal.probc4l);
			else
				fprintf(file_output,"LR(M4 vs M%d) = %g. Probability not calculated\n",data[0].nloci,
					LR4l);
		}

		/*THE BEST MODEL*/
		if(datams[0][0].bestmodel == 1) fprintf(file_output,"\nTHE BEST MODEL is M%d. (A single theta/nt per all loci).\n\n",datams[0][0].bestmodel);
		else {
			if(datams[0][0].bestmodel > 1 && datams[0][0].bestmodel < 5) 
				fprintf(file_output,"\nTHE BEST MODEL is M%d. (%d different theta/nt distributed in the loci).\n\n",datams[0][0].bestmodel,datams[0][0].bestmodel);
			if(datams[0][0].bestmodel == 5) fprintf(file_output,"\nTHE BEST MODEL is M%d. (Each loci has a different theta/nt).\n\n",datams[0][0].n_loci);
		}
	}
    
    /*FREE MEMORY AND CLOSE THE FILE*/
    for(j=0;j<=data[0].nloci;j++) free(ml[j]);
    free(ml);
    free(sumlocus);
    free(locusml);
    free(thetalocusml);
    free(sortedlociml);
    free(stillloci);
    free(data[0].nsamA);
    free(data[0].factor_chrn);
    free(data[0].length);
    free(data[0].SA);
    #if EQTAVARE == 0
    for(j=0;j<=maxn;j++) free(Q[j]);
    free(Q);
    for(j=0;j<=maxn;j++) free(P[j]);
    free(P);
    #endif
    free(thetaM2);
    free(mlM2);
    free(thetaM3);
    free(mlM3);
    free(thetaM4);
    free(mlM4);
    free(vector2ml1);
    free(vector2ml2);
    free(vector3ml1);
    free(vector3ml2);
    free(vector3ml3);
    free(vector4ml1);
    free(vector4ml2);
    free(vector4ml3);
    free(vector4ml4);
    if(data[0].nloci > 20) {
        free(relthetaM3);
        free(relmlM3);
        free(relthetaM4);
        free(relmlM4);
    }
	free(data[0].numberloci);
    free(data);

    return 0;
}

/*compare two numbers*/
int compare(const void *i,const void *j) 
{
    return (*((int *)i) - (*(int *)j));
}

/************** calculate Hudson 12 recursive function Oxford surveys in evol. biol. 1-44 **********************/
double dfhudson(int S,double theta, int nsam, double ***Q,double ***P)
{
    int y,i,j;
    double pn;
    double TtalPn = (double)0.;
    double dfthetaSteffi(int,int, double **,double ***);
    
    /*calculate all possible Q: nsam*S values. Init P*/
    for(i=2;i<=nsam;i++) { 
        for(j=0;j<=S;j++) {
            Q[0][i][j] = (double)1.0;
            pn = theta/(theta + i - (double)1.0);
            for(y=0;y<j;y++) Q[0][i][j] *= pn;
            Q[0][i][j] *= (i - (double)1.0)/(theta + i - (double)1.0);
            P[0][i][j] = -1.;
        }
    }
    /*calculate the total probability*/
    TtalPn = dfthetaSteffi(S,nsam,Q[0],&P[0]);
    if(TtalPn == (double)0.) return (double)1E-35;       
    return TtalPn;
}

double dfthetaSteffi(int S, int nsam, double **Q, double ***P)
{
    int x;
    double df = (double)0.0;
    
    if(nsam > 2) {
        for(x=0;x<=S;x++) {
            if(P[0][nsam-1][S-x] != -1.0) df += Q[nsam][x] * P[0][nsam-1][S-x];
            else
                if(Q[nsam][x] != (double)0.)
                    df += Q[nsam][x] * dfthetaSteffi(S-x,nsam-1,Q,&P[0]);
        }
    }
    else df += Q[nsam][S];
    P[0][nsam][S] = df;

    if(df < (double)0) 
		return (double)0;
	if(df > (double)1) 
		return (double)1;

	return df;
}

/************* calculate equation 9.5 from Tavaré 1984 Theor. Pop. Biol. 26: 119-164 *******************************/
double dftavare(int S,double theta,int nsam)
{
    int l;
    double P,sum;
    double bico2(int,int);
    
    double pow1,bin,powe;
/*    
    S=23;
    theta = 57.0429993;
    nsam=12;
*/    
    sum = (double)0;
    for(l=1;l<nsam;l++) {
        pow1 = pow((double)(-1),(double)(l-1));
        bin  = (double)bico2(nsam-2,l-1);
        powe = (double)pow((double)theta/((double)l+(double)theta),(double)(S+1));
        sum += pow1 * bin * powe;
    }
    P = ((double)nsam - (double)1)/(double)theta * sum;
	/**/    
    if(P > (double)1) {
        printf("\nError at P[%d,%g,%d]=%g",S,theta,nsam,P);
	}
	/**/
    if(P <= (double)5e-12) {
        /*printf("\nError at P[%d,%g,%d]=%g",S,theta,nsam,P);*/
		P = 5e-12; /*No more precission*/
	}
	return (double)P;
}

double bico2(int n, int k) /*based in NR in C. page 215*/
{
	int x,y;
	double fact;
	double diffless,diffmore;
	
	if(n-k < k) 
		k = n-k;
	
	fact = (double)1;
	for(x=n,y=k;x>n-k;x--,y--) 
		fact *= (double)x/(double)y;
	
	diffless = fact - floor(fact);
	diffmore = ceil(fact) - fact;
	
	if(diffless < diffmore) 
		fact = floor(fact);
	else
		fact = ceil(fact);
	return (double)fact;
}


