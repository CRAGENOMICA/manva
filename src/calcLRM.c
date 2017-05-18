/*
 *  calcLRM.c
 *  thetamlS
 *
 *  Created by sebas on Thu Apr 29 2004.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NLOCISEQ 20

#define RNDSIM 500
#define MHSIM 500
#define EQTAVARE 1

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

int calculateLRM(struct parameters_thetaml *data,int *SC,
    double *LR12,double *LR13,double *LR14,double *LR1l,
    double *LR23,double *LR24,double *LR2l,
    double *LR34,double *LR3l,
    double *LR4l,
    int model,
    double **ml,
    double *thetalocusml,double *locusml,
    int *sortedlociml,int *stillloci,
    double *mlM2,
    double *mlM3,double *relmlM3,
    double *mlM4,double *relmlM4,
    int option1234l)
{
    long int i;
    int j,k,l,m,n;
    double thetaA,mlM1;
    int suml1,suml2,suml3;
    double ran;
    
    double sortthetamlmin;
    int sortlocimin;
    int findml2coalsorted(int,int *,int,int,double *,double **);
    
    int nloci3m1,nloci3m2;
    int relnloci3m1,relnloci3m2;
	double relM3totc;
    int findml3coalsorted(int,int,int *,int,int,double *,double **);
    int findml3coalsortedmh(int,int,int *,int,int,double *,double *,double **,double);
    
    int nloci4m1,nloci4m2,nloci4m3;
    int relnloci4m1,relnloci4m2,relnloci4m3;
	double relM4totc;
    int findml4coalsorted(int,int,int,int *,int,int,double *,double **);
    int findml4coalsortedmh(int,int,int,int *,int,int, double *,double *,double **,double);

    double logranloci;
    double ran1();    
    int compare(const void *,const void *);

    #if EQTAVARE == 0
    double **Q,**P;
    int maxS,maxn;
    double dfhudson(int,double,int,double ***,double ***);
    #else
    double dftavare(int,double,int);
    #endif
	double valuedf;
    
    double bestmlbefore;
    int best1,best2,best3;
    
    #if EQTAVARE == 0
    /*Define max bound for S and n (useful for defining Q and P)*/
    maxS = maxn = 0;
    for(j=0;j<data[0].nloci;j++) {
        if(data[0].nsamA[j] > maxn) maxn = data[0].nsamA[j];
        if(SC[j] > maxS) maxS = SC[j];
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

    /*Calculate the probabilities for the user bounds*/
    for(i=0;i<data[0].steps_user;i++) 
        ml[data[0].nloci][i] = (double)0.;/*init*/
    for(j=0;j<data[0].nloci;j++) {
        for(i=0;i<data[0].steps_user;i++) { 
            thetaA = ((double)i/((double)data[0].steps_user-(double)1.))*
                     (data[0].thetamax_user*(double)data[0].length[j] - data[0].thetamin_user*(double)data[0].length[j]) + 
                      data[0].thetamin_user*(double)data[0].length[j];/*uniform*/
			#if EQTAVARE == 0
			valuedf =  dfhudson(SC[j],thetaA*data[0].factor_chrn[j],data[0].nsamA[j],&Q,&P);
			if(valuedf <= (double)1e-323)
				ml[j][i] = (double)-1000;
			else 
				ml[j][i] = (double)log(valuedf);/*logprobability with SC data*/
			#else
			valuedf =  dftavare(SC[j],thetaA*data[0].factor_chrn[j],data[0].nsamA[j]);
			if(valuedf <= (double)1e-323)
				ml[j][i] = (double)-1000;
			else 
				ml[j][i] = (double)log(valuedf);/*logprobability with SC data*/
			#endif
            /*ml[j][i] = (double)log(dftavare(SC[j],thetaA*data[0].factor_chrn[j],data[0].nsamA[j]));*//*logprobability with SC data*/
            ml[data[0].nloci][i] += ml[j][i];/*sum of likelihood for each iteration of all loci*/
            if(i==0 || locusml[j] < ml[j][i]) {
                locusml[j] = ml[j][i];
                thetalocusml[j] = (double)i;
            }
        }
    }
    /*total maximum likelihood. M1 and Mnloci models*/
    locusml[data[0].nloci] = (double)0.;
    for(j=0;j<data[0].nloci;j++) {
        thetalocusml[j] = (thetalocusml[j]/((double)data[0].steps_user-(double)1.))*
                          (data[0].thetamax_user - data[0].thetamin_user) + data[0].thetamin_user;
        locusml[data[0].nloci] += locusml[j];
    }
    mlM1 = (double)0.;
    for(i=0;i<data[0].steps_user;i++) {
        if(i==0 || mlM1 < ml[data[0].nloci][i]) {
            mlM1 = ml[data[0].nloci][i];
        }
    }
    if(data[0].nloci >= 2) {
        /*loci are sorted by thetaM1 value*/
        for(j=0;j<data[0].nloci;j++) stillloci[j] = j;
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
        if(model < 3) {
            /*M2 model*/
            /*do all lineal combinations (heuristic) and calculate ml*/
            mlM2[0] = mlM2[1] = mlM1;
            for(j=0;j<data[0].nloci-1;j++) {
                if(findml2coalsorted(j+1,sortedlociml,data[0].nloci,data[0].steps_user,mlM2,ml)) return(1);
            }
        }
    }
    if(!(model == 1 && option1234l == 1)) {
        if(data[0].nloci >= 3 && model < 4) {
            /*M3 model*/
            if(data[0].nloci <= NLOCISEQ) {
                /*heuristic but all lineal combinations*/
                mlM3[0] = mlM3[1] = mlM3[2] = mlM1;
                for(j=0;j<data[0].nloci-2;j++) {
                    for(k=j+1,l=0;k<data[0].nloci-1;k++,l++) {
                        if(findml3coalsorted(j+1,l+1,sortedlociml,data[0].nloci,data[0].steps_user,mlM3,ml)) return(1);
                    }
                }
            }
            else {
                /*Find randomly.*/
                /*random using iterations (RNDSIM/MHSIM)*/
                /*random values*/
                do {
					mlM3[0] = mlM3[1] = mlM3[2] = mlM1;
					for(n=0;n<RNDSIM;n++) {
						nloci3m1 = (int)((double)ran1()*(double)(data[0].nloci-2)) + 1;
						nloci3m2 = (int)((double)ran1()*(double)(data[0].nloci-1-nloci3m1)) + 1;
						bestmlbefore = mlM3[0]+mlM3[1]+mlM3[2];
						if(findml3coalsorted(nloci3m1,nloci3m2,sortedlociml,data[0].nloci,data[0].steps_user,mlM3,ml))
							return(1);
						if(bestmlbefore < mlM3[0]+mlM3[1]+mlM3[2]) {
							best1 = nloci3m1;
							best2 = nloci3m2;
						}
					}
					relmlM3[0] = mlM3[0];
					relmlM3[1] = mlM3[1];
					relmlM3[2] = mlM3[2];
					nloci3m1 = best1;
					nloci3m2 = best2;
					relnloci3m1 = suml1 = nloci3m1;
					relnloci3m2 = nloci3m2;
					suml2 = suml1 + nloci3m2;
					/*Now MHMCMC. values close to the best ml*/
					for(n=0;n<MHSIM;n++) {
						if((ran=ran1()) < (double)0.333) {
							if(suml1 == 1) suml1 = 2;
							else suml1 -= 1;
							if((ran=ran1()) < (double)0.333) {
								if(suml2 <= suml1 + 1) suml2 = suml1 + 1;
								else if(suml2 > suml1 + 1) suml2 -= 1;
							}
							else {
								if(ran > (double)0.666)
									if(suml2 < data[0].nloci - 3) suml2 += 1;
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
									}
									if(ran > (double)0.666)
										if(suml2 < data[0].nloci - 3) suml2 += 1;
								}
							}
							else {
								if((ran=ran1()) < (double)0.333) {
									if(suml2 <= suml1 + 1) suml2 = suml1 + 1;
									else suml2 -= 1;
								}
								if(ran > (double)0.666)
									if(suml2 < data[0].nloci - 3) suml2 += 1;
							}
						}
						nloci3m1 = suml1;
						nloci3m2 = suml2 - suml1;

						logranloci = (double)log((double)ran1()/(double)data[0].nloci);
						relM3totc = relmlM3[0] + relmlM3[1] + relmlM3[2];
						if(findml3coalsortedmh(nloci3m1,nloci3m2,sortedlociml,data[0].nloci,
							data[0].steps_user,relmlM3,mlM3,ml,logranloci))
							return(1);
						if(relM3totc == relmlM3[0] + relmlM3[1] + relmlM3[2]) {
							/*rejected*/
							nloci3m1 = relnloci3m1;
							nloci3m2 = relnloci3m2;
						}
						else {
							/*accepted*/
							relnloci3m1 = nloci3m1;
							relnloci3m2 = nloci3m2;
						}
					}
				}while(model<3 && (mlM2[0] + mlM2[1] > mlM3[0]+mlM3[1]+mlM3[2]));
            }
        }
        if(!(model == 2 && option1234l == 1)) {
            if(data[0].nloci >= 4) {
                /*M4 model*/
                if(data[0].nloci <= NLOCISEQ) {
                    /*heuristic lineal all*/
                    mlM4[0] = mlM4[1] = mlM4[2] = mlM4[3] = mlM1;
                    for(j=0;j<data[0].nloci-3;j++) {
                        for(k=j+1,l=0;k<data[0].nloci-2;k++,l++) {
                            for(n=k+1,m=0;n<data[0].nloci-1;n++,m++) {
                                if(findml4coalsorted(j+1,l+1,m+1,sortedlociml,data[0].nloci,data[0].steps_user,mlM4,ml))
                                    return(1);
                            }
                        }
                    }
                }
                else {
                    /*Find randomly*/
                    /*random using iterations (RNDSIM/MHSIM)*/
                    do {
						mlM4[0] = mlM4[1] = mlM4[2] = mlM4[3] = mlM1;
						for(n=0;n<RNDSIM;n++) {/*random values*/
							nloci4m1 = (int)(ran1()*(double)(data[0].nloci-3)) + 1;
							nloci4m2 = (int)(ran1()*(double)(data[0].nloci-2-nloci4m1)) + 1;
							nloci4m3 = (int)(ran1()*(double)(data[0].nloci-1-nloci4m1-nloci4m2)) + 1;
							bestmlbefore = mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3];
							if(findml4coalsorted(nloci4m1,nloci4m2,nloci4m3,sortedlociml,data[0].nloci,
								data[0].steps_user,mlM4,ml))
								return(1);
							if(bestmlbefore < mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]) {
								best1 = nloci4m1;
								best2 = nloci4m2;
								best3 = nloci4m3;
							}
						}
						relmlM4[0] = mlM4[0];
						relmlM4[1] = mlM4[1];
						relmlM4[2] = mlM4[2];
						relmlM4[3] = mlM4[3];
						nloci4m1 = best1;
						nloci4m2 = best2;
						nloci4m3 = best3;
						relnloci4m1 = nloci4m1;
						relnloci4m2 = nloci4m2;
						relnloci4m3 = nloci4m3;            
						suml1 = nloci4m1;
						suml2 = suml1 + nloci4m2;
						suml3 = suml2 + nloci4m3;
						/*Now MHMCMC. values close to the best ml*/
						for(n=0;n<MHSIM;n++) {
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
											if(suml3 < data[0].nloci - 3) suml3 += 1;
									}
								}
								else {
									if(ran > (double)0.666) {
										if(suml2 < data[0].nloci - 4) suml2 += 1;
										if((ran=ran1()) < (double)0.333) {
											if(suml3 <= suml2 + 1) suml3 = suml2 + 1;
											else suml3 -= 1;
										}
										else {
											if(ran > (double)0.666)
												if(suml3 < data[0].nloci - 3) suml3 += 1;
										}
									}
								}
							}
							else {
								if(ran > (double)0.666) {
									if(suml1 >= data[0].nloci - 5) {
										suml1 = data[0].nloci - 4;
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
													if(suml3 < data[0].nloci - 3) suml3 += 1;
											}
										}
										if(ran > (double)0.666) {
											if(suml2 < data[0].nloci - 4) suml2 += 1;
											if((ran=ran1()) < (double)0.333) {
												if(suml3 <= suml2 + 1) suml3 = suml2 + 1;
												else suml3 -= 1;
											}
											else {
												if(ran > (double)0.666)
													if(suml3 < data[0].nloci - 3) suml3 += 1;
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
												if(suml3 < data[0].nloci - 3) suml3 += 1;
										}
									}
									if(ran > (double)0.666) {
										if(suml2 < data[0].nloci - 4) suml2 += 1;
										if((ran=ran1()) < (double)0.333) {
											if(suml3 <= suml2 + 1) suml3 = suml2 + 1;
											else suml3 -= 1;
										}
										else {
											if(ran > (double)0.666)
												if(suml3 < data[0].nloci - 3) suml3 += 1;
										}
									}
								}
							}
							nloci4m1 = suml1;
							nloci4m2 = suml2 - suml1;
							nloci4m3 = suml3 - suml2;
			
							logranloci = (double)log(ran1()/(double)data[0].nloci);
							relM4totc = relmlM4[0] + relmlM4[1] + relmlM4[2] + relmlM4[3];
							if(findml4coalsortedmh(nloci4m1,nloci4m2,nloci4m3,sortedlociml,data[0].nloci,
								data[0].steps_user,relmlM4,mlM4,ml,
								logranloci)) return(1);
							if(relM4totc == relmlM4[0] + relmlM4[1] + relmlM4[2] + relmlM4[3]) {
								/*rejected*/
								nloci4m1 = relnloci4m1;
								nloci4m2 = relnloci4m2;
								nloci4m3 = relnloci4m3;
							}
							else {
								/*accepted*/
								relnloci4m1 = nloci4m1;
								relnloci4m2 = nloci4m2;
								relnloci4m3 = nloci4m3;
							}
						}
					}while(model<4 && (mlM3[0] + mlM3[1] + mlM3[2] > mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]));
                }
            }
        }
    }
    /*likelihood ratio tests*/
    if(model == 1) *(LR1l) = (double)2 * (-mlM1 + locusml[data[0].nloci]);
    if(data[0].nloci > 1) {
        if(model == 1) *(LR12) = (double)2 * (-mlM1 + (mlM2[0]+mlM2[1]));
        if(model == 2) *(LR2l) = (double)2 * (-(mlM2[0]+mlM2[1]) + locusml[data[0].nloci]);
    }
    if(data[0].nloci > 2) {
        if(model == 1) *(LR13) = (double)2 * (-mlM1 + (mlM3[0]+mlM3[1]+mlM3[2]));
        if(model == 2) *(LR23) = (double)2 * (-(mlM2[0]+mlM2[1]) + (mlM3[0]+mlM3[1]+mlM3[2]));
        if(model == 3) *(LR3l) = (double)2 * (-(mlM3[0]+mlM3[1]+mlM3[2]) + locusml[data[0].nloci]);
    }
    if(data[0].nloci > 3) {
        if(model == 1) *(LR14) = (double)2 * (-mlM1 + (mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]));
        if(model == 2) *(LR24) = (double)2 * (-(mlM2[0]+mlM2[1]) + (mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]));
        if(model == 3) *(LR34) = (double)2 * (-(mlM3[0]+mlM3[1]+mlM3[2]) +
                                            (mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]));
        if(model == 4) *(LR4l) = (double)2 * (-(mlM4[0]+mlM4[1]+mlM4[2]+mlM4[3]) + locusml[data[0].nloci]);
    }

    #if EQTAVARE == 0
    for(j=0;j<=maxn;j++) free(Q[j]);
    free(Q);
    for(j=0;j<=maxn;j++) free(P[j]);
    free(P);
    #endif

    return 0;
}

