/*
 *  pcoalthetaml.c
 *  thetamlS
 *
 *  Created by sonsins on Wed Apr 28 2004.
 *
 */

#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265359

#define SIG 0.05
#define NLOCISEQ 20
#define OPTION1234l 0
/*option to do only LRT12,LRT23,LRT34 and LRT4l (0 do all)*/
#define PRINTLRTDIST 0

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

int probcoalthhaetasml(struct parameters_thetaml *data,
    double thetaM1,double *thetaM2,double *thetaM3,double *thetaM4,
    int *vector2ml1,int *vector2ml2,
    int *vector3ml1,int *vector3ml2,int *vector3ml3,
    int *vector4ml1,int *vector4ml2,int *vector4ml3,int *vector4ml4,
    int nloci2m1,int nloci2m2,
    int nloci3m1,int nloci3m2,int nloci3m3,
    int nloci4m1,int nloci4m2,int nloci4m3,int nloci4m4,
    double LR12,double LR13,double LR14,double LR1l,
    double LR23,double LR24,double LR2l,
    double LR34,double LR3l,double LR4l,
    struct Prob *Pcoal,
	struct LRTdist **lrdist/*,struct MLthetasim ***mthetasim*/)
{
    int it,j,n;
    double t;
    int   *SC=0;
    double *LRP12=0;
    double *LRP13=0;
    double *LRP14=0;
    double *LRP1l=0;
    double *LRP23=0;
    double *LRP24=0;
    double *LRP2l=0;
    double *LRP34=0;
    double *LRP3l=0;
    double *LRP4l=0;

    double **ml=0;
    double *thetalocusml=0;
    double *locusml=0;
    int   *sortedlociml=0;
    int   *stillloci=0;
    
    double *mlM2=0;
    double *mlM3=0;
    double *relmlM3=0;
    double *mlM4=0;
    double *relmlM4=0;
    
    double ran1();
    int compare_(const void *,const void *);
    int calculateLRM(struct parameters_thetaml *,int *,
        double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int,
        double **,double *,double *,int *,int *,double *,double *,double *,
        double *,double *,int);
    double poisso(double);

    int printdot;
    #if PRINTLRTDIST
    FILE *file_outputLRT;
    #endif
    printdot = (int)ceil((double)data[0].ncoalsim/10.);

    /*define matrix and vectors*/
    if(!(ml = (double **)malloc((data[0].nloci+1)*sizeof(double *)))) return(1);
    for(j=0;j<=data[0].nloci;j++) {
        if(!(ml[j] = (double *)malloc((data[0].steps_user)*sizeof(double)))) return(1);
    }
    if(!(thetalocusml = (double *)malloc((data[0].nloci+1)*sizeof(double)))) return(1);
    if(!(locusml = (double *)malloc((data[0].nloci+1)*sizeof(double)))) return(1);
    if(!(sortedlociml = (int *)malloc((data[0].nloci)*sizeof(int)))) return(1);
    if(!(stillloci = (int *)malloc((data[0].nloci)*sizeof(int)))) return(1);
    /*vectors for M2 model*/
    if(!(mlM2 = (double *)malloc(2*sizeof(double)))) return(1);    
    /*vectors for M3 model*/
    if(!(mlM3 = (double *)malloc(3*sizeof(double)))) return(1);
    /*vectors for M4 model*/
    if(!(mlM4 = (double *)malloc(4*sizeof(double)))) return(1);    
    /*MHMCMC*/
    if(data[0].nloci > NLOCISEQ) {
        if(!(relmlM3 = (double *)malloc(3*sizeof(double)))) return(1);
        if(!(relmlM4 = (double *)malloc(4*sizeof(double)))) return(1);
    }   
    /*SC vector and LRP matrix*/
    if(!(LRP12 = (double *)malloc((data[0].ncoalsim)*sizeof(double)))) return(1);
	memset(LRP12,'\0',(data[0].ncoalsim)*sizeof(double));
    if(!(LRP13 = (double *)malloc((data[0].ncoalsim)*sizeof(double)))) return(1);
	memset(LRP13,'\0',(data[0].ncoalsim)*sizeof(double));
    if(!(LRP14 = (double *)malloc((data[0].ncoalsim)*sizeof(double)))) return(1);
	memset(LRP14,'\0',(data[0].ncoalsim)*sizeof(double));
    if(!(LRP1l = (double *)malloc((data[0].ncoalsim)*sizeof(double)))) return(1);
	memset(LRP1l,'\0',(data[0].ncoalsim)*sizeof(double));
    if(!(LRP23 = (double *)malloc((data[0].ncoalsim)*sizeof(double)))) return(1);
	memset(LRP23,'\0',(data[0].ncoalsim)*sizeof(double));
    if(!(LRP24 = (double *)malloc((data[0].ncoalsim)*sizeof(double)))) return(1);
	memset(LRP24,'\0',(data[0].ncoalsim)*sizeof(double));
    if(!(LRP2l = (double *)malloc((data[0].ncoalsim)*sizeof(double)))) return(1);
	memset(LRP2l,'\0',(data[0].ncoalsim)*sizeof(double));
    if(!(LRP34 = (double *)malloc((data[0].ncoalsim)*sizeof(double)))) return(1);
	memset(LRP34,'\0',(data[0].ncoalsim)*sizeof(double));
    if(!(LRP3l = (double *)malloc((data[0].ncoalsim)*sizeof(double)))) return(1);
	memset(LRP3l,'\0',(data[0].ncoalsim)*sizeof(double));
    if(!(LRP4l = (double *)malloc((data[0].ncoalsim)*sizeof(double)))) return(1);
	memset(LRP4l,'\0',(data[0].ncoalsim)*sizeof(double));
    if(!(SC = (int *)malloc((data[0].nloci)*sizeof(int)))) return(1);
	memset(SC,'\0',(data[0].nloci)*sizeof(int));
    
    /*modeling M1*/
    printf("\n\nCoalescent simulations for model M1 ");
    fflush(stdout);

    for(it=0;it<data[0].ncoalsim;it++) {
        for(n=0;n<data[0].nloci;n++) {
            /*coalescent*/
            t=0.;
            for(j=data[0].nsamA[n];j>1;j--)
                t += -(double)log((double)1.0 - ran1())/((double)j*(j-(double)1.0)/data[0].factor_chrn[n]) * (double)j;
            SC[n] = (int)poisso(thetaM1*(double)data[0].length[n] * t);
        }
        /*calculate ML*/
        if(calculateLRM(data,SC,LRP12+it,LRP13+it,LRP14+it,LRP1l+it,NULL,NULL,NULL,NULL,NULL,NULL,1,
            ml,thetalocusml,locusml,sortedlociml,stillloci,mlM2,mlM3,relmlM3,
            mlM4,relmlM4,OPTION1234l)) return(2);
			
		lrdist[0][it].LRT12 = LRP12[it];
        if((double)it/(double)printdot == (int)it/(int)printdot) {
            printf(".");
            fflush(stdout);
        }
    }
    /*calculate probabilities for model M1*/
    if(data[0].ncoalsim) {
		if(data[0].nloci >= 2) {
			qsort(LRP12,data[0].ncoalsim,sizeof(double),compare_);
			for(j=0;j<data[0].ncoalsim;j++) if(LR12 <= LRP12[j]) break;
			Pcoal[0].probc12 = (double)(data[0].ncoalsim - j)/(double)data[0].ncoalsim;
			#if PRINTLRTDIST
			if (!(file_outputLRT = fopen ("LRT12dist.out","w"))) {
				puts("Error in input/output");
				return(1);
			}
			for(j=0;j<data[0].ncoalsim;j++) fprintf(file_outputLRT,"%g\n",LRP12[j]);
			fclose(file_outputLRT); 
			#endif
		}
		else Pcoal[0].probc12 = (double)10000.;
		
		qsort(LRP1l,data[0].ncoalsim,sizeof(double),compare_);
		for(j=0;j<data[0].ncoalsim;j++) if(LR1l <= LRP1l[j]) break;
		Pcoal[0].probc1l = (double)(data[0].ncoalsim - j)/(double)data[0].ncoalsim;
		#if PRINTLRTDIST
		if (!(file_outputLRT = fopen ("LRT1ldist.out","w"))) {
			puts("Error in input/output");
			return(1);
		}
		for(j=0;j<data[0].ncoalsim;j++) fprintf(file_outputLRT,"%g\n",LRP1l[j]);
		fclose(file_outputLRT); 
		#endif
		if(data[0].nloci >= 3) {
			#if OPTION1234l == 0
			qsort(LRP13,data[0].ncoalsim,sizeof(double),compare_);
			for(j=0;j<data[0].ncoalsim;j++) if(LR13 <= LRP13[j]) break;
			Pcoal[0].probc13 = (double)(data[0].ncoalsim - j)/(double)data[0].ncoalsim;
			#if PRINTLRTDIST
			if (!(file_outputLRT = fopen ("LRT13dist.out","w"))) {
				puts("Error in input/output");
				return(1);
			}
			for(j=0;j<data[0].ncoalsim;j++) fprintf(file_outputLRT,"%g\n",LRP13[j]);
			fclose(file_outputLRT); 
			#endif
		}
		else Pcoal[0].probc13 = (double)10000.;
		
		if(data[0].nloci >= 4) {
			qsort(LRP14,data[0].ncoalsim,sizeof(double),compare_);
			for(j=0;j<data[0].ncoalsim;j++) if(LR14 <= LRP14[j]) break;
			Pcoal[0].probc14 = (double)(data[0].ncoalsim - j)/(double)data[0].ncoalsim;
			#if PRINTLRTDIST
			if (!(file_outputLRT = fopen ("LRT14dist.out","w"))) {
				puts("Error in input/output");
				return(1);
			}
			for(j=0;j<data[0].ncoalsim;j++) fprintf(file_outputLRT,"%g\n",LRP14[j]);
			fclose(file_outputLRT); 
			#endif
			#else 
			Pcoal[0].probc13 = Pcoal[0].probc14 = (double)10000.;
			#endif
		}
		else Pcoal[0].probc14 = (double)10000.;
    }
    else {
       Pcoal[0].probc12 = Pcoal[0].probc13 = Pcoal[0].probc14 = Pcoal[0].probc1l = 10000.;
    }
    if(data[0].nloci > 2) {
        /*modeling M2*/
        if(/*Pcoal[0].probc12 <*/ SIG) {
            printf("\nCoalescent simulations for model M2 ");
            fflush(stdout);
            for(it=0;it<data[0].ncoalsim;it++) {
                for(n=0;n<nloci2m1;n++) {
                    /*coalescent*/
                    t=0.;
                    for(j=data[0].nsamA[vector2ml1[n]];j>1;j--)
                        t += -(double)log((double)1.0 - ran1())/((double)j*(j-(double)1.0)/data[0].factor_chrn[vector2ml1[n]]) * (double)j;
                    SC[vector2ml1[n]] = (int)poisso(thetaM2[0]*(double)data[0].length[vector2ml1[n]] * t);
                }
                for(n=0;n<nloci2m2;n++) {
                    /*coalescent*/
                    t=0.;
                    for(j=data[0].nsamA[vector2ml2[n]];j>1;j--)
                        t += -(double)log((double)1.0 - ran1())/((double)j*(j-(double)1.0)/data[0].factor_chrn[vector2ml2[n]]) * (double)j;
                    SC[vector2ml2[n]] = (int)poisso(thetaM2[1]*(double)data[0].length[vector2ml2[n]] * t);
                }
                /*calculate ML*/
                if(calculateLRM(data,SC,NULL,NULL,NULL,NULL,LRP23+it,LRP24+it,LRP2l+it,NULL,NULL,NULL,2,
                    ml,thetalocusml,locusml,sortedlociml,stillloci,mlM2,mlM3,relmlM3,
                    mlM4,relmlM4,OPTION1234l)) return(2);

				lrdist[0][it].LRT23 = LRP23[it];
                if((double)it/(double)printdot == (int)it/(int)printdot) {
                    printf(".");
                    fflush(stdout);
                }
            }
            /*calculate probabilities for model M2*/
			if(data[0].nloci >= 3) {
				qsort(LRP23,data[0].ncoalsim,sizeof(double),compare_);
				for(j=0;j<data[0].ncoalsim;j++) if(LR23 <= LRP23[j]) break;
				Pcoal[0].probc23 = (double)(data[0].ncoalsim - j)/(double)data[0].ncoalsim;
				#if PRINTLRTDIST
				if (!(file_outputLRT = fopen ("LRT23dist.out","w"))) {
					puts("Error in input/output");
					return(1);
				}
				for(j=0;j<data[0].ncoalsim;j++) fprintf(file_outputLRT,"%g\n",LRP23[j]);
				fclose(file_outputLRT); 
				#endif
			}
			else Pcoal[0].probc23 = (double)10000.;
			
            qsort(LRP2l,data[0].ncoalsim,sizeof(double),compare_);
            for(j=0;j<data[0].ncoalsim;j++) if(LR2l <= LRP2l[j]) break;
            Pcoal[0].probc2l = (double)(data[0].ncoalsim - j)/(double)data[0].ncoalsim;
            #if PRINTLRTDIST
            if (!(file_outputLRT = fopen ("LRT2ldist.out","w"))) {
                puts("Error in input/output");
                return(1);
            }
            for(j=0;j<data[0].ncoalsim;j++) fprintf(file_outputLRT,"%g\n",LRP2l[j]);
            fclose(file_outputLRT); 
            #endif
			if(data[0].nloci >= 4) {
				#if OPTION1234l == 0
				qsort(LRP24,data[0].ncoalsim,sizeof(double),compare_);
				for(j=0;j<data[0].ncoalsim;j++) if(LR24 <= LRP24[j]) break;
				Pcoal[0].probc24 = (double)(data[0].ncoalsim - j)/(double)data[0].ncoalsim;
				#if PRINTLRTDIST
				if (!(file_outputLRT = fopen ("LRT24dist.out","w"))) {
					puts("Error in input/output");
					return(1);
				}
				for(j=0;j<data[0].ncoalsim;j++) fprintf(file_outputLRT,"%g\n",LRP24[j]);
				fclose(file_outputLRT); 
				#endif
				#else
				Pcoal[0].probc24 = (double)10000.;
				#endif
			}
			else Pcoal[0].probc24 = (double)10000.;
        }
        else Pcoal[0].probc23 = Pcoal[0].probc24 = Pcoal[0].probc2l = (double)10000.;
        if(data[0].nloci > 3) {
            /*modeling M3*/
            if(/*Pcoal[0].probc13 < SIG || Pcoal[0].probc23 < */SIG) {
                printf("\nCoalescent simulations for model M3 ");
                fflush(stdout);
                for(it=0;it<data[0].ncoalsim;it++) {
                    for(n=0;n<nloci3m1;n++) {
                        /*coalescent*/
                        t=0.;
                        for(j=data[0].nsamA[vector3ml1[n]];j>1;j--)
                            t += -(double)log((double)1.0 - ran1())/((double)j*(j-(double)1.0)/data[0].factor_chrn[vector3ml1[n]]) * (double)j;
                        SC[vector3ml1[n]] = (int)poisso(thetaM3[0]*(double)data[0].length[vector3ml1[n]] * t);
                    }
                    for(n=0;n<nloci3m2;n++) {
                        /*coalescent*/
                        t=0.;
                        for(j=data[0].nsamA[vector3ml2[n]];j>1;j--)
                            t += -(double)log((double)1.0 - ran1())/((double)j*(j-(double)1.0)/data[0].factor_chrn[vector3ml2[n]]) * (double)j;
                        SC[vector3ml2[n]] = (int)poisso(thetaM3[1]*(double)data[0].length[vector3ml2[n]] * t);
                    }
                    for(n=0;n<nloci3m3;n++) {
                        /*coalescent*/
                        t=0.;
                        for(j=data[0].nsamA[vector3ml3[n]];j>1;j--)
                            t += -(double)log((double)1.0 - ran1())/((double)j*(j-(double)1.0)/data[0].factor_chrn[vector3ml3[n]]) * (double)j;
                        SC[vector3ml3[n]] = (int)poisso(thetaM3[2]*(double)data[0].length[vector3ml3[n]] * t);
                    }
                    /*calculate ML*/
                    if(calculateLRM(data,SC,NULL,NULL,NULL,NULL,NULL,NULL,NULL,LRP34+it,LRP3l+it,NULL,3,
                        ml,thetalocusml,locusml,sortedlociml,stillloci,mlM2,mlM3,relmlM3,
                        mlM4,relmlM4,OPTION1234l)) return(2);

					lrdist[0][it].LRT34 = LRP34[it];
                    if((double)it/(double)printdot == (int)it/(int)printdot) {
                        printf(".");
                        fflush(stdout);
                    }
                }
                /*calculate probabilities for M3*/
				if(data[0].nloci >= 4) {
					qsort(LRP34,data[0].ncoalsim,sizeof(double),compare_);
					for(j=0;j<data[0].ncoalsim;j++) if(LR34 <= LRP34[j]) break;
					Pcoal[0].probc34 = (double)(data[0].ncoalsim - j)/(double)data[0].ncoalsim;
					#if PRINTLRTDIST
					if (!(file_outputLRT = fopen ("LRT34dist.out","w"))) {
						puts("Error in input/output");
						return(1);
					}
					for(j=0;j<data[0].ncoalsim;j++) fprintf(file_outputLRT,"%g\n",LRP34[j]);
					fclose(file_outputLRT); 
					#endif
				}
				else Pcoal[0].probc34 = (double)10000.;
				
                qsort(LRP3l,data[0].ncoalsim,sizeof(double),compare_);
                for(j=0;j<data[0].ncoalsim;j++) if(LR3l <= LRP3l[j]) break;
                Pcoal[0].probc3l = (double)(data[0].ncoalsim - j)/(double)data[0].ncoalsim;
                #if PRINTLRTDIST
                if (!(file_outputLRT = fopen ("LRT3ldist.out","w"))) {
                    puts("Error in input/output");
                    return(1);
                }
                for(j=0;j<data[0].ncoalsim;j++) fprintf(file_outputLRT,"%g\n",LRP3l[j]);
                fclose(file_outputLRT); 
                #endif
            }
            else Pcoal[0].probc34 = Pcoal[0].probc3l = (double)10000.;
            
            if(data[0].nloci > 4) {
                /*modeling M4*/
                if(/*Pcoal[0].probc14 < SIG || Pcoal[0].probc24 < SIG || Pcoal[0].probc34 < */SIG) {
                    printf("\nCoalescent simulations for model M4 ");
                    fflush(stdout);
                    for(it=0;it<data[0].ncoalsim;it++) {
                        for(n=0;n<nloci4m1;n++) {
                            /*coalescent*/
                            t=0.;
                            for(j=data[0].nsamA[vector4ml1[n]];j>1;j--)
                                t += -(double)log((double)1.0 - ran1())/((double)j*(j-(double)1.0)/data[0].factor_chrn[vector4ml1[n]]) * (double)j;
                            SC[vector4ml1[n]] = (int)poisso(thetaM4[0]*(double)data[0].length[vector4ml1[n]] * t);
                        }
                        for(n=0;n<nloci4m2;n++) {
                            /*coalescent*/
                            t=0.;
                            for(j=data[0].nsamA[vector4ml2[n]];j>1;j--)
                                t += -(double)log((double)1.0 - ran1())/((double)j*(j-(double)1.0)/data[0].factor_chrn[vector4ml2[n]]) * (double)j;
                            SC[vector4ml2[n]] = (int)poisso(thetaM4[1]*(double)data[0].length[vector4ml2[n]] * t);
                        }
                        for(n=0;n<nloci4m3;n++) {
                            /*coalescent*/
                            t=0.;
                            for(j=data[0].nsamA[vector4ml3[n]];j>1;j--)
                                t += -(double)log((double)1.0 - ran1())/((double)j*(j-(double)1.0)/data[0].factor_chrn[vector4ml3[n]]) * (double)j;
                            SC[vector4ml3[n]] = (int)poisso(thetaM4[2]*(double)data[0].length[vector4ml3[n]] * t);
                        }
                        for(n=0;n<nloci4m4;n++) {
                            /*coalescent*/
                            t=0.;
                            for(j=data[0].nsamA[vector4ml4[n]];j>1;j--)
                                t += -(double)log((double)1.0 - ran1())/((double)j*(j-(double)1.0)/data[0].factor_chrn[vector4ml4[n]]) * (double)j;
                            SC[vector4ml4[n]] = (int)poisso(thetaM4[3]*(double)data[0].length[vector4ml4[n]] * t);
                        }
                        /*calculate ML*/
                        if(calculateLRM(data,SC,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,LRP4l+it,4,
                            ml,thetalocusml,locusml,sortedlociml,stillloci,mlM2,mlM3,relmlM3,
                            mlM4,relmlM4,OPTION1234l)) return(2);

						lrdist[0][it].LRT4l = LRP4l[it];
                        if((double)it/(double)printdot == (int)it/(int)printdot) {
                            printf(".");
                            fflush(stdout);
                        }
                    }
                    /*calculate probabilities for M4*/
                    qsort(LRP4l,data[0].ncoalsim,sizeof(double),compare_);
                    for(j=0;j<data[0].ncoalsim;j++) if(LR4l <= LRP4l[j]) break;
                    Pcoal[0].probc4l = (double)(data[0].ncoalsim - j)/(double)data[0].ncoalsim;
                    #if PRINTLRTDIST
                    if (!(file_outputLRT = fopen ("LRT4ldist.out","w"))) {
                        puts("Error in input/output");
                        return(1);
                    }
                    for(j=0;j<data[0].ncoalsim;j++) fprintf(file_outputLRT,"%g\n",LRP4l[j]);
                    fclose(file_outputLRT); 
                    #endif
                }
                else Pcoal[0].probc4l = (double)10000.;
            }
        }
    }

    /*FREE MEMORY*/
    free(SC);
    free(LRP12);
    free(LRP13);
    free(LRP14);
    free(LRP1l);
    free(LRP23);
    free(LRP24);
    free(LRP2l);
    free(LRP34);
    free(LRP3l);
    free(LRP4l);
    
    for(j=0;j<=data[0].nloci;j++) free(ml[j]);
    free(ml);
    free(locusml);
    free(thetalocusml);
    free(sortedlociml);
    free(stillloci);
    free(mlM2);
    free(mlM3);
    free(mlM4);
    if(data[0].nloci > 20) {
        free(relmlM3);
        free(relmlM4);
    }

    return 0;                        
}
