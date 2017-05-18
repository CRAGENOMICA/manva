/*
 *  displaysimdata2.c
 *  MuLoNeTests
 *
 *  Created by sebas on Tue Mar 18 2003.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ALLRESULTS	1

void open_simdisplayprob(struct var *data,struct var2b *inputms,struct statistics *matrix,struct statmulo *matrixml,struct statistisim **matrixsim,struct statistisimmuloc *matrixmlsim,struct horizontalstatsml *avgstatloci,int *observed_data,int *outgroup/*,int *n_loci,int *montecarlo_sim*/,int onlymulo, FILE *file_output,int dataobsequalsim, int neuttest,int mulo)
{
    /*Display probabilities*/
    char k[1];
    int x;
	int maxnoutgroups;
    long int it,it2;
    long int nunderalfa[30];
    long int nequivalfa[30];
    long int nitertests[30];
    double varobs[30],prob;
    void printfignif(double,FILE *);
    void printfignif1(double,FILE *);
    double *sortvector;
    int comp(const void *,const void *);
    double avg,var,vars;
    double varianceit(double,double,long int);
    double variances(double,double,int);
    /*long int wilcoxonrst_ts(double *,double *,int);*/
    /*double wilcoxontstonormal(int,long int);*/
    /*double table_wilcoxonsr(int,double);*/
    double prob_cumbinomial(int, int,double);
    int *psigned;
    int *nsigned;
	int *qsigned;
	int eq,ps,ns;
	double **ibonf;
	int **ibonfn;
    void print_percentages(struct var *,double *,long int, FILE *);
	void print_ibonf(double,int,int,long int,FILE *);
	void ibonf_sort(double **,int **,int,int);
    
	switch(mulo) {
        case 0: /*all loci*/
            switch(neuttest) {
                case 0: /* statistics, multilocus*/
                    #if COMMAND_LINE
                    if(file_output) fflush(file_output);		
                    printf("\n\nMENU 3. Coalescent Monte Carlo simulations and Statistical inference:");
                    printf("\n     3.1. Display coalescent simulation analysis menu:");
                    printf("\n     3.1.0. Display statistics menu:");
                    printf("\n     3.1.0.1. Display detailed statistics menu:\n");
                    printf("\n     3.1.0.1.1 Display table with detailed statistics for all loci data:\n\n");
                    
                    printf("CHOOSE:\n");
                    printf(" 0 - Display general statistics.\n");
                    printf(" 1 - Display other linkage related statistics.\n");
                    printf(" 2 - Display estimates of total locus variability.\n");
                    printf(" 3 - Display estimates of nucleotide locus variability.\n");
                    printf(" 4 - Back to previous menu.\n");

                    if(file_output) {
                        fprintf(file_output,"\n\n     MENU:\n     3. Coalescent Monte Carlo simulations and Statistical inference:");
                        fprintf(file_output,"\n     3.1. Display coalescent simulation analysis menu:");
                        fprintf(file_output,"\n     3.1.0. Display statistics menu:");
                        fprintf(file_output,"\n     3.1.0.1. Display detailed statistics menu:\n");
                        fprintf(file_output,"\n     3.1.0.1.1 Display table with detailed statistics for all loci data:\n\n");
                        
                        fprintf(file_output," 0 - Display general statistics.\n");
                        fprintf(file_output," 1 - Display other linkage related statistics.\n");
                        fprintf(file_output," 2 - Display estimates of total locus variability.\n");
                        fprintf(file_output," 3 - Display estimates of nucleotide locus variability.\n");
                        fprintf(file_output," 4 - Back to previous menu.\n");
                    }
                    #endif
                    
                    do *k = getchar();
                    while(*k<'0' || *k>'4');

                    if(file_output) fprintf(file_output,"\nOPTION CHOSEN: %c\n\n",*k);

					/*Improved Bonferroni test. Definitions: An Improved Bonferroni Procedure for Multiple Tests of Significance. 
					 R. J. Simes. Biometrika, Vol. 73, No. 3. (Dec., 1986), pp. 751-754.*/
					if((ibonf = (double **)calloc(10,sizeof(double *))) == 0) {
						printf("Error in memory allocation: IBt.");
						if(file_output) 
							fputs("Error in memory allocation: IBt.",file_output);
						break;
					}
					for(x=0;x<10;x++) {
						if((ibonf[x] = (double *)calloc(data[0].n_loci,sizeof(double))) == 0) {
							printf("Error in memory allocation: IBt.");
							if(file_output) 
								fputs("Error in memory allocation: IBt.",file_output);
							break;
						}
					}
					if((ibonfn = (int **)calloc(10,sizeof(int *))) == 0) {
						printf("Error in memory allocation: IBt.");
						if(file_output) 
							fputs("Error in memory allocation: IBt.",file_output);
						break;
					}
					for(x=0;x<10;x++) {
						if((ibonfn[x] = (int *)calloc(data[0].n_loci,sizeof(int))) == 0) {
							printf("Error in memory allocation: IBt.");
							if(file_output) 
								fputs("Error in memory allocation: IBt.",file_output);
							break;
						}
					}
					
                    switch(*k) {
                        case '0': /* 0 - Display general statistics.*/
                            if(*observed_data && dataobsequalsim) {/*compare simulations with observed data.*/
								/*check if more than 1 outgroup*/
								maxnoutgroups = 0;
								if(*outgroup) {
									for(x=0;x<data[0].n_loci;x++) {
										if(matrix[x].noutgroups > maxnoutgroups) 
											maxnoutgroups = matrix[x].noutgroups;
									}
								} 
                                /*Sign tests.*/
                                /*Binomial for each statistic. from median (all data) or average (onlymulo). outgroup and non outgroup.*/
                                if((psigned = (int *)calloc(9,sizeof(int))) == 0) {
                                    printf("Error in memory allocation: Sign tests.");
                                    if(file_output) 
                                        fputs("Error in memory allocation: Sign tests.",file_output);
                                    break;
                                }
                                if((nsigned = (int *)calloc(9,sizeof(int))) == 0) {
                                    printf("Error in memory allocation: Sign tests.");
                                    if(file_output) 
                                        fputs("Error in memory allocation: Sign tests.",file_output);
                                    break;
                                }

                                if(onlymulo) {
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\nNon-parametric SIGN TEST. Assuming that average have equivalent value than median. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl)");
                                    printf("\nComparing observed data with the average of simulated data at each locus.");
                                    printf("\nAdjust to a Binomial.\n");
                                    printf(" The table contains in order the name of the locus,");
                                    printf(" the number of segregating biallelic sites,");
                                    if(data[0].time_spec > 0.) {
										if(maxnoutgroups == 1) printf(" fixed mutations, divergence, ");
										else printf(" divergence");
									}
                                    printf(" the number of haplotypes (not shared, if observed), haplotypes/sample size, haplotype diversity, and the recombination parameter (Hudson 1987) and Rm (Hudson and Kaplan 1985).\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\nNon-parametric SIGN TEST. Assuming that average have equivalent value than median. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl)",file_output);
                                        fputs("\nComparing observed data with the average of simulated data at each locus.",file_output);
                                        fputs("\nAdjust to a Binomial.\n",file_output);
                                        fprintf(file_output," The table contains in order the name of the locus,");
                                        fprintf(file_output," the number of segregating biallelic sites,");
                                        if(data[0].time_spec > 0.) {
											if(maxnoutgroups == 1) fprintf(file_output," fixed mutations, divergence,");
											else fprintf(file_output," divergence");
										}
                                        fprintf(file_output," the number of haplotypes (not shared, if observed), haplotypes/sample size, haplotype diversity, and the recombination parameter (Hudson 1987) and Rm (Hudson and Kaplan 1985).\n");
                                    }
                                }
                                else {
                                    printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									printf("\nNon-parametric SIGN TEST. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl)");
                                    printf("\nComparing observed data with the median of simulated data at each locus.");
                                    printf("\nAdjust to a Binomial.\n");
                                    printf(" The table contains in order the name of the locus,");
                                    printf(" the number of segregating biallelic sites,");
                                    if(data[0].time_spec > 0.) {
										if(maxnoutgroups == 1) printf(" fixed mutations, divergence");
										else printf(" divergence,");
									}
                                    printf(" the number of haplotypes (not shared, if observed), haplotypes/sample size, haplotype diversity, and the recombination parameter (Hudson 1987) and Rm (Hudson and Kaplan 1985).\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\nNon-parametric SIGN TEST . \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl).",file_output);
                                        fputs("\nComparing observed data with the median of simulated data at each locus.",file_output);
                                        fputs("\nAdjust to a Binomial.\n",file_output);
                                        fprintf(file_output," The table contains in order the name of the locus,");
                                        fprintf(file_output," the number of segregating biallelic sites,");
                                        if(data[0].time_spec > 0.) {
											if(maxnoutgroups == 1) fprintf(file_output," fixed mutations, divergence");
											else fprintf(file_output," divergence,");
										}
                                        fprintf(file_output," the number of haplotypes (not shared, if observed), haplotypes/sample size, haplotype diversity, and the recombination parameter (Hudson 1987) and Rm (Hudson and Kaplan 1985).\n");
                                    }
                                }
                                if(data[0].time_spec > 0.) {
                                    if(maxnoutgroups == 1)  puts("\nloci\tObs(S)\tExp(S)\tObs(fix)\tExp(fix)\tObs(Div)\tExp(Div)\tObs(hapl)\tExp(hapl)\tObs(haplsam)\tExp(haplsam)\tObs(hapldiv)\tExp(hapldiv)\tObs(C)\tExp(C)\tObs(Rm)\tExp(Rm)\n");
                                    else puts("\nloci\tObs(S)\tExp(S)\tObs(Div)\tExp(Div)\tObs(hapl)\tExp(hapl)\tObs(haplsam)\tExp(haplsam)\tObs(hapldiv)\tExp(hapldiv)\tObs(Rm)\tExp(Rm)\n");
									if(file_output) {
                                        if(maxnoutgroups == 1)  fputs("\nloci\tObs(S)\tExp(S)\tObs(fix)\tExp(fix)\tObs(Div)\tExp(Div)\tObs(hapl)\tExp(hapl)\tObs(haplsam)\tExp(haplsam)\tObs(hapldiv)\tExp(hapldiv)\tObs(C)\tExp(C)\tObs(Rm)\tExp(Rm)\n",file_output);
										else fputs("\nloci\tObs(S)\tExp(S)\tObs(Div)\tExp(Div)\tObs(hapl)\tExp(hapl)\tObs(haplsam)\tExp(haplsam)\tObs(hapldiv)\tExp(hapldiv)\tObs(Rm)\tExp(Rm)\n",file_output);
									}
									for(x=0;x<data[0].n_loci;x++) {
                                        if(maxnoutgroups == 1) 
											printf("%d:%s\t%d\t%g\t%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%g\n",x,matrix[x].gene,
                                            matrix[x].biallsites,avgstatloci[x].biallsites,
                                            matrix[x].fixed,avgstatloci[x].fixed,
                                            matrix[x].ndivergence,avgstatloci[x].ndivergence,
                                            matrix[x].nhapl,avgstatloci[x].nhapl,
                                            matrix[x].nhaplsam,avgstatloci[x].nhaplsam,
                                            matrix[x].hapldiv,avgstatloci[x].hapldiv,
											matrix[x].Rm,avgstatloci[x].Rm);
										else
											printf("%d:%s\t%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%g\n",x,matrix[x].gene,
                                            matrix[x].biallsites,avgstatloci[x].biallsites,
                                            matrix[x].ndivergence,avgstatloci[x].ndivergence,
                                            matrix[x].nhapl,avgstatloci[x].nhapl,
                                            matrix[x].nhaplsam,avgstatloci[x].nhaplsam,
                                            matrix[x].hapldiv,avgstatloci[x].hapldiv,
											matrix[x].Rm,avgstatloci[x].Rm);
                                        if(file_output) {
                                            if(maxnoutgroups == 1) 
												fprintf(file_output,"%d:%s\t%d\t%g\t%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%g\n",x,matrix[x].gene,
                                                matrix[x].biallsites,avgstatloci[x].biallsites,
                                                matrix[x].fixed,avgstatloci[x].fixed,
                                                matrix[x].ndivergence,avgstatloci[x].ndivergence,
                                                matrix[x].nhapl,avgstatloci[x].nhapl,
                                                matrix[x].nhaplsam,avgstatloci[x].nhaplsam,
                                                matrix[x].hapldiv,avgstatloci[x].hapldiv,
												matrix[x].Rm,avgstatloci[x].Rm);
											else 
												fprintf(file_output,"%d:%s\t%d\t%g\t%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%g\n",x,matrix[x].gene,
                                                matrix[x].biallsites,avgstatloci[x].biallsites,
                                                matrix[x].ndivergence,avgstatloci[x].ndivergence,
                                                matrix[x].nhapl,avgstatloci[x].nhapl,
                                                matrix[x].nhaplsam,avgstatloci[x].nhaplsam,
                                                matrix[x].hapldiv,avgstatloci[x].hapldiv,
												matrix[x].Rm,avgstatloci[x].Rm);
                                        }
										if(matrix[x].biallsites > avgstatloci[x].biallsites) psigned[0] += 1;
                                        if(matrix[x].fixed > avgstatloci[x].fixed) psigned[1] += 1;
                                        if(matrix[x].ndivergence > avgstatloci[x].ndivergence) psigned[2] += 1;
                                        if(matrix[x].nhapl > avgstatloci[x].nhapl) psigned[3] += 1;
                                        if(matrix[x].nhaplsam > avgstatloci[x].nhaplsam) psigned[4] += 1;
                                        if(matrix[x].hapldiv > avgstatloci[x].hapldiv) psigned[5] += 1;
                                        if(matrix[x].Rm > avgstatloci[x].Rm) psigned[6] += 1;

                                        if(matrix[x].biallsites < avgstatloci[x].biallsites) nsigned[0] += 1;
                                        if(matrix[x].fixed < avgstatloci[x].fixed) nsigned[1] += 1;
                                        if(matrix[x].ndivergence < avgstatloci[x].ndivergence) nsigned[2] += 1;
                                        if(matrix[x].nhapl < avgstatloci[x].nhapl) nsigned[3] += 1;
                                        if(matrix[x].nhaplsam < avgstatloci[x].nhaplsam) nsigned[4] += 1;
                                        if(matrix[x].hapldiv < avgstatloci[x].hapldiv) nsigned[5] += 1;
                                        if(matrix[x].Rm < avgstatloci[x].Rm) nsigned[6] += 1;
                                    }
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    for(x=0;x<7;x++) {
										if(x==1 && maxnoutgroups > 1) x = 2; 
                                        switch(x) {
                                            case 0:
                                                printf("\nSign test: Segregating sites.\n");
                                                if(file_output) fputs("\nSign test: Segregating sites.\n",file_output);
                                                break;
                                            case 1:
                                                printf("\nSign test: Fixed mutations.\n");
                                                if(file_output) fputs("\nSign test: Fixed mutations.\n",file_output);
                                                break;
                                            case 2:
                                                printf("\nSign test: Divergence.\n");
                                                if(file_output) fputs("\nSign test: Divergence.\n",file_output);
                                                break;
                                            case 3:
                                                printf("\nSign test: Number of haplotypes.\n");
                                                if(file_output) fputs("\nSign test: Number of haplotypes.\n",file_output);
                                                break;
                                            case 4:
                                                printf("\nSign test: Number of haplotypes / sample size.\n");
                                                if(file_output) fputs("\nSign test: Number of haplotypes / sample size.\n",file_output);
                                                break;
                                            case 5:
                                                printf("\nSign test: Haplotype diversity.\n");
                                                if(file_output) fputs("\nSign test: Haplotype diversity.\n",file_output);
                                                break;
                                            case 6:
                                                printf("\nSign test: Rm parameter (Hudson and Kaplan 1985).\n");
                                                if(file_output) fputs("\nSign test: Rm parameter (Hudson and Kaplan 1985).\n",file_output);
												break;
                                        }
                                        prob = prob_cumbinomial(psigned[x]+nsigned[x],psigned[x],0.5);
										if(prob == -10000) printf("#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
                                        else printf("#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
										if(file_output) {
											if(prob == -10000) fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
											else fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
										}
										printfignif1(prob,file_output);
										
										/*look at the effect of discrete values*/
										if(prob <= 0.05 && prob >= 0.) {
											eq = data[0].n_loci - (psigned[x]+nsigned[x]);
											if(psigned[x] > nsigned[x]) {
												ns = nsigned[x] + eq;
												ps = psigned[x];
											}
											else {
												ns = nsigned[x];
												ps = psigned[x] + eq;
											}
											if((psigned[x] > nsigned[x] && ps < ns) || (psigned[x] < nsigned[x] && ps > ns) || prob_cumbinomial(ps+ns,ps,0.5) > 0.05) {
                                                printf(" WARNING: too many values equal to the median. It might be too liberal.");
                                                if(file_output) fputs(" WARNING: too many values equal to the median. It might be too liberal.",file_output);
											}
										}
                                    }
                                }
                                else {
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\nloci\tObs(S)\tExp(S)\tObs(hapl)\tExp(hapl)\tObs(haplsam)\tExp(haplsam)\tObs(hapldiv)\tExp(hapldiv)\tObs(Rm)\tExp(Rm)\n");
									if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    if(file_output)
                                        fputs("\nloci\tObs(S)\tExp(S)\tObs(hapl)\tExp(hapl)\tObs(haplsam)\tExp(haplsam)\tObs(hapldiv)\tExp(hapldiv)\tObs(Rm)\tExp(Rm)\n",file_output);
                                    for(x=0;x<data[0].n_loci;x++) {
                                        printf("%d:%s\t%d\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%g\n",x,matrix[x].gene,
                                            matrix[x].biallsites,avgstatloci[x].biallsites,
                                            matrix[x].nhapl,avgstatloci[x].nhapl,
                                            matrix[x].nhaplsam,avgstatloci[x].nhaplsam,
                                            matrix[x].hapldiv,avgstatloci[x].hapldiv,
											matrix[x].Rm,avgstatloci[x].Rm);
                                        if(file_output)
                                            fprintf(file_output,"%d:%s\t%d\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%g\n",x,matrix[x].gene,
                                                matrix[x].biallsites,avgstatloci[x].biallsites,
                                                matrix[x].nhapl,avgstatloci[x].nhapl,
                                                matrix[x].nhaplsam,avgstatloci[x].nhaplsam,
                                                matrix[x].hapldiv,avgstatloci[x].hapldiv,
												matrix[x].Rm,avgstatloci[x].Rm);

                                        if(matrix[x].biallsites > avgstatloci[x].biallsites) psigned[0] += 1;
                                        if(matrix[x].nhapl > avgstatloci[x].nhapl) psigned[1] += 1;
                                        if(matrix[x].nhaplsam > avgstatloci[x].nhaplsam) psigned[2] += 1;
                                        if(matrix[x].hapldiv > avgstatloci[x].hapldiv) psigned[3] += 1;
                                        if(matrix[x].Rm > avgstatloci[x].Rm) psigned[4] += 1;

                                        if(matrix[x].biallsites < avgstatloci[x].biallsites) nsigned[0] += 1;
                                        if(matrix[x].nhapl < avgstatloci[x].nhapl) nsigned[1] += 1;
                                        if(matrix[x].nhaplsam < avgstatloci[x].nhaplsam) nsigned[2] += 1;
                                        if(matrix[x].hapldiv < avgstatloci[x].hapldiv) nsigned[3] += 1;
                                        if(matrix[x].Rm < avgstatloci[x].Rm) nsigned[4] += 1;
                                    }
                                    for(x=0;x<5;x++) {
                                        switch(x) {
                                            case 0:
                                                printf("\nSign test: Segregating sites.\n");
                                                if(file_output) fputs("\nSign test: Segregating sites.\n",file_output);
                                                break;
                                            case 1:
                                                printf("\nSign test: Number of haplotypes.\n");
                                                if(file_output) fputs("\nSign test: Number of haplotypes.\n",file_output);
                                                break;
                                            case 2:
                                                printf("\nSign test: Number of haplotypes / sample size.\n");
                                                if(file_output) fputs("\nSign test: Number of haplotypes / sample size.\n",file_output);
                                                break;
                                            case 3:
                                                printf("\nSign test: Haplotype diversity.\n");
                                                if(file_output) fputs("\nSign test: Haplotype diversity.\n",file_output);
                                                break;
                                            case 4:
                                                printf("\nSign test: Rm parameter (Hudson and Kaplan 1985).\n");
                                                if(file_output) fputs("\nSign test: Rm parameter (Hudson and Kaplan 1985).\n",file_output);
                                                break;
                                        }
                                        prob = prob_cumbinomial(psigned[x]+nsigned[x],psigned[x],0.5);
										if(prob == -10000) printf("#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
                                        else printf("#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
										if(file_output) {
											if(prob == -10000) fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
											else fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
										}
                                        printfignif1(prob,file_output);
										
										/*look at the effect of discrete values*/
										if(prob <= 0.05 && prob >= 0.) {
											eq = data[0].n_loci - (psigned[x]+nsigned[x]);
											if(psigned[x] > nsigned[x]) {
												ns = nsigned[x] + eq;
												ps = psigned[x];
											}
											else {
												ns = nsigned[x];
												ps = psigned[x] + eq;
											}
											if((psigned[x] > nsigned[x] && ps < ns) || (psigned[x] < nsigned[x] && ps > ns) || prob_cumbinomial(ps+ns,ps,0.5) > 0.05) {
                                                printf(" WARNING: too many values equal to the median. It might be too liberal.");
                                                if(file_output) fputs(" WARNING: too many values equal to the median. It might be too liberal.",file_output);
											}
										}
                                    }
                                }
                                free(psigned);
                                free(nsigned);

                                if(!onlymulo) {
                                    /*calculate probability each observed value/locus*/
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\n\nDisplay the probability of observed statistics for each locus obtained by coalescent Monte Carlo simulations: ");
                                    printf("\nE.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. Note: & means the ierations would be increased to confirm significance");
                                    printf("\nMultiple hits not included in the analysis.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\n\nDisplay the probability of observed statistics for each locus obtained by coalescent Monte Carlo simulations: ",file_output);
                                        fprintf(file_output,"\nE.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. Note: & means the ierations would be increased to confirm significance");
                                        fputs("\nMultiple hits not included in the analysis.\n",file_output);
                                    }
                                    if(*outgroup) {
                                        if(data[0].time_spec <= 0.) {
                                            printf("\nloci\tObs(S)\tP(S)\tIBonf(S)\tObs(nhapl)\tP(nhapl)\tIBonf(nhapl)\tObs(nhaplsam)\tP(nhaplsam)\tIBonf(nhaplsam)\tObs(hapldiv)\tP(hapldiv)\tIBonf(hapldiv)\tObs(Rm)\tP(Rm)\tIBonf(Rm)\n");
                                            if(file_output)
                                                fputs("\nloci\tObs(S)\tP(S)\tIBonf(S)\tObs(nhapl)\tP(nhapl)\tIBonf(nhapl)\tObs(nhaplsam)\tP(nhaplsam)\tIBonf(nhaplsam)\tObs(hapldiv)\tP(hapldiv)\tIBonf(hapldiv)\tObs(Rm)\tP(Rm)\tIBonf(Rm)\n",file_output);
                                            for(x=0;x<data[0].n_loci;x++) {
                                                nunderalfa[0] = nunderalfa[1] = nunderalfa[2] = nunderalfa[3] = nunderalfa[4]  = nunderalfa[5] = 0;
                                                nequivalfa[0] = nequivalfa[1] = nequivalfa[2] = nequivalfa[3] = nequivalfa[4] = nequivalfa[5] = 0;
                                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                                    if(matrix[x].biallsites < matrixsim[x][it].biallsites)
                                                        nunderalfa[0] += 1;
                                                    if(matrix[x].nhapl < matrixsim[x][it].nhapl)
                                                        nunderalfa[1] += 1;
                                                    if(matrix[x].nhaplsam < matrixsim[x][it].nhaplsam)
                                                        nunderalfa[2] += 1;
                                                    if(matrix[x].hapldiv < matrixsim[x][it].hapldiv)
                                                        nunderalfa[3] += 1;
                                                    if(matrix[x].Rm < matrixsim[x][it].Rm)
                                                        nunderalfa[4] += 1;
                                                    if(matrix[x].biallsites == matrixsim[x][it].biallsites)
                                                        nequivalfa[0] += 1;
                                                    if(matrix[x].nhapl == matrixsim[x][it].nhapl)
                                                        nequivalfa[1] += 1;
                                                    if(matrix[x].nhaplsam == matrixsim[x][it].nhaplsam)
                                                        nequivalfa[2] += 1;
                                                    if(matrix[x].hapldiv == matrixsim[x][it].hapldiv)
                                                        nequivalfa[3] += 1;
                                                    if(matrix[x].Rm == matrixsim[x][it].Rm)
                                                        nequivalfa[4] += 1;
                                                }
                                                ibonf[0][x] = (double)nunderalfa[0]/(double)data[0].n_iter;
												if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                                if(ibonf[0][x] > (double)0.5) ibonf[0][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                                else ibonf[0][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                                
												ibonf[1][x] = (double)nunderalfa[1]/(double)data[0].n_iter;
												if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                                if(ibonf[1][x] > (double)0.5) ibonf[1][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                                else ibonf[1][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                                
												ibonf[2][x] = (double)nunderalfa[2]/(double)data[0].n_iter;
												if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                                if(ibonf[2][x] > (double)0.5) ibonf[2][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                                else ibonf[2][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                                
												ibonf[3][x] = (double)nunderalfa[3]/(double)data[0].n_iter;
												if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                                if(ibonf[3][x] > (double)0.5) ibonf[3][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                                else ibonf[3][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                                
												ibonf[4][x] = (double)nunderalfa[4]/(double)data[0].n_iter;
												if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                                if(ibonf[4][x] > (double)0.5) ibonf[4][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                                else ibonf[4][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];												
                                            }
											ibonf_sort(ibonf,ibonfn,5,data[0].n_loci);
											
											for(x=0;x<data[0].n_loci;x++) {
                                                printf("%d:%s\t",x,matrix[x].gene);
                                                if(file_output) fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
												
                                                printf("%d\t%g\t",matrix[x].biallsites,ibonf[0][x]);
                                                if(file_output) fprintf(file_output,"%d\t%g\t",matrix[x].biallsites,ibonf[0][x]);
												print_ibonf(ibonf[0][x],ibonfn[0][x],data[0].n_loci,data[0].n_iter,file_output);
												
                                                printf("%d\t%g\t",matrix[x].nhapl,ibonf[1][x]);
                                                if(file_output) fprintf(file_output,"%d\t%g\t",matrix[x].nhapl,ibonf[1][x]);
												print_ibonf(ibonf[1][x],ibonfn[1][x],data[0].n_loci,data[0].n_iter,file_output);	
																							
                                                printf("%g\t%g\t",matrix[x].nhaplsam,ibonf[2][x]);
                                                if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].nhaplsam,ibonf[2][x]);
												print_ibonf(ibonf[2][x],ibonfn[2][x],data[0].n_loci,data[0].n_iter,file_output);
																								
                                                printf("%g\t%g\t",matrix[x].hapldiv,ibonf[3][x]);
                                                if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].hapldiv,ibonf[3][x]);
												print_ibonf(ibonf[3][x],ibonfn[3][x],data[0].n_loci,data[0].n_iter,file_output);
																																			
                                                printf("%d\t%g\t",matrix[x].Rm,ibonf[4][x]);
                                                if(file_output) fprintf(file_output,"%d\t%g\t",matrix[x].Rm,ibonf[4][x]);
												print_ibonf(ibonf[5][x],ibonfn[4][x],data[0].n_loci,data[0].n_iter,file_output);	
																							
												printf("\n");
												if(file_output) fprintf(file_output,"\n");
											}
                                        }
                                        else {
											printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                            if(maxnoutgroups == 1)  printf("\nloci\tObs(S)\tP(S)\tIBonf(S)\tObs(Fix)\tP(Fix)\tIBonf(Fix)\tObs(Div)\tP(Div)\tIBonf(Div)\tObs(nhapl)\tP(nhapl)\tIBonf(nhapl)\tObs(nhaplsam)\tP(nhaplsam)\tIBonf(nhaplsam)\tObs(hapldiv)\tP(hapldiv)\tIBonf(hapldiv)\tObs(Rm)\tP(Rm)\tIBonf(Rm)\n");
											else  printf("\nloci\tObs(S)\tP(S)\tIBonf(S)\tObs(Div)\tP(Div)\tIBonf(Div)\tObs(nhapl)\tP(nhapl)\tIBonf(nhapl)\tObs(nhaplsam)\tP(nhaplsam)\tIBonf(nhaplsam)\tObs(hapldiv)\tP(hapldiv)\tIBonf(hapldiv)\tObs(Rm)\tP(Rm)\tIBonf(Rm)\n");
											if(file_output) {
												if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                                if(maxnoutgroups == 1) fputs("\nloci\tObs(S)\tP(S)\tIBonf(S)\tObs(Fix)\tP(Fix)\tIBonf(Fix)\tObs(Div)\tP(Div)\tIBonf(Div)\tObs(nhapl)\tP(nhapl)\tIBonf(nhapl)\tObs(nhaplsam)\tP(nhaplsam)\tIBonf(nhaplsam)\tObs(hapldiv)\tP(hapldiv)\tIBonf(hapldiv)\tObs(Rm)\tP(Rm)\tIBonf(Rm)\n",file_output);
												else fputs("\nloci\tObs(S)\tP(S)\tIBonf(S)\tObs(Div)\tP(Div)\tIBonf(Div)\tObs(nhapl)\tP(nhapl)\tIBonf(nhapl)\tObs(nhaplsam)\tP(nhaplsam)\tIBonf(nhaplsam)\tObs(hapldiv)\tP(hapldiv)\tIBonf(hapldiv)\tObs(Rm)\tP(Rm)\tIBonf(Rm)\n",file_output);
                                            }
											for(x=0;x<data[0].n_loci;x++) {
                                                nunderalfa[0] = nunderalfa[1] = nunderalfa[2] = nunderalfa[3] = nunderalfa[4] = nunderalfa[5] = nunderalfa[6] = nunderalfa[7] = 0;
                                                nequivalfa[0] = nequivalfa[1] = nequivalfa[2] = nequivalfa[3] = nequivalfa[4] = nequivalfa[5] = nequivalfa[6] = nequivalfa[7] = 0;
                                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                                    if(matrix[x].biallsites < matrixsim[x][it].biallsites)
                                                        nunderalfa[0] += 1;
                                                    if(matrix[x].fixed < matrixsim[x][it].fixed)
                                                        nunderalfa[1] += 1;
                                                    if(matrix[x].ndivergence < matrixsim[x][it].ndivergence)
                                                        nunderalfa[2] += 1;
                                                    if(matrix[x].nhapl < matrixsim[x][it].nhapl)
                                                        nunderalfa[3] += 1;
                                                    if(matrix[x].nhaplsam < matrixsim[x][it].nhaplsam)
                                                        nunderalfa[4] += 1;
                                                    if(matrix[x].hapldiv < matrixsim[x][it].hapldiv)
                                                        nunderalfa[5] += 1;
                                                   if(matrix[x].Rm < matrixsim[x][it].Rm)
                                                        nunderalfa[6] += 1;
                                                    if(matrix[x].biallsites == matrixsim[x][it].biallsites)
                                                        nequivalfa[0] += 1;
                                                    if(matrix[x].fixed == matrixsim[x][it].fixed)
                                                        nequivalfa[1] += 1;
                                                    if(matrix[x].ndivergence == matrixsim[x][it].ndivergence)
                                                        nequivalfa[2] += 1;
                                                    if(matrix[x].nhapl == matrixsim[x][it].nhapl)
                                                        nequivalfa[3] += 1;
                                                    if(matrix[x].nhaplsam == matrixsim[x][it].nhaplsam)
                                                        nequivalfa[4] += 1;
                                                    if(matrix[x].hapldiv == matrixsim[x][it].hapldiv)
                                                        nequivalfa[5] += 1;
                                                    if(matrix[x].Rm == matrixsim[x][it].Rm)
                                                        nequivalfa[6] += 1;
                                                }
                                                ibonf[0][x] = (double)nunderalfa[0]/(double)data[0].n_iter;
												if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                                if(ibonf[0][x] > (double)0.5) ibonf[0][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                                else ibonf[0][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
												
												ibonf[1][x] = (double)nunderalfa[1]/(double)data[0].n_iter;
												if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                                if(ibonf[1][x] > (double)0.5) ibonf[1][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                                else ibonf[1][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];

                                                ibonf[2][x] = (double)nunderalfa[2]/(double)data[0].n_iter;
												if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                                if(ibonf[2][x] > (double)0.5) ibonf[2][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                                else ibonf[2][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];

                                                ibonf[3][x] = (double)nunderalfa[3]/(double)data[0].n_iter;
												if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                                if(ibonf[3][x] > (double)0.5) ibonf[3][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                                else ibonf[3][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];

                                                ibonf[4][x] = (double)nunderalfa[4]/(double)data[0].n_iter;
												if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                                if(ibonf[4][x] > (double)0.5) ibonf[4][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                                else ibonf[4][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];

                                                ibonf[5][x] = (double)nunderalfa[5]/(double)data[0].n_iter;
												if((double)1./((double)data[0].n_iter) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                                if(ibonf[5][x] > (double)0.5) ibonf[4][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                                else ibonf[5][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];

                                                ibonf[6][x] = (double)nunderalfa[6]/(double)data[0].n_iter;
												if((double)1./((double)data[0].n_iter) *(double)nequivalfa[6] >= (double)0.5) nequivalfa[6] /= 2;
                                                if(ibonf[6][x] > (double)0.5) ibonf[6][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[6];
                                                else ibonf[6][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[6];      
                                            }
											ibonf_sort(ibonf,ibonfn,7,data[0].n_loci);
											
											for(x=0;x<data[0].n_loci;x++) {
												printf("%d:%s\t",x,matrix[x].gene);
                                                if(file_output) fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
												
                                                printf("%d\t%g\t",matrix[x].biallsites,ibonf[0][x]);
                                                if(file_output) fprintf(file_output,"%d\t%g\t",matrix[x].biallsites,ibonf[0][x]);
                                                print_ibonf(ibonf[0][x],ibonfn[0][x],data[0].n_loci,data[0].n_iter,file_output);
												
												if(maxnoutgroups == 1) printf("%d\t%g\t",matrix[x].fixed,ibonf[1][x]);
                                                if(maxnoutgroups == 1) if(file_output) fprintf(file_output,"%d\t%g\t",matrix[x].fixed,ibonf[1][x]);
                                                if(maxnoutgroups == 1) print_ibonf(ibonf[1][x],ibonfn[1][x],data[0].n_loci,data[0].n_iter,file_output);
												
												printf("%g\t%g\t",matrix[x].ndivergence,ibonf[2][x]);
                                                if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].ndivergence,ibonf[2][x]);
                                                print_ibonf(ibonf[2][x],ibonfn[2][x],data[0].n_loci,data[0].n_iter,file_output);
												
												printf("%d\t%g\t",matrix[x].nhapl,ibonf[3][x]);
                                                if(file_output) fprintf(file_output,"%d\t%g\t",matrix[x].nhapl,ibonf[3][x]);
												print_ibonf(ibonf[3][x],ibonfn[3][x],data[0].n_loci,data[0].n_iter,file_output);
												
												printf("%g\t%g\t",matrix[x].nhaplsam,ibonf[4][x]);
                                                if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].nhaplsam,ibonf[4][x]);
												print_ibonf(ibonf[4][x],ibonfn[4][x],data[0].n_loci,data[0].n_iter,file_output);
												
												printf("%g\t%g\t",matrix[x].hapldiv,ibonf[5][x]);
                                                if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].hapldiv,ibonf[5][x]);
												print_ibonf(ibonf[5][x],ibonfn[5][x],data[0].n_loci,data[0].n_iter,file_output);
																								
												printf("%d\t%g\t",matrix[x].Rm,ibonf[6][x]);
                                                if(file_output) fprintf(file_output,"%d\t%g\t",matrix[x].Rm,ibonf[6][x]);
												print_ibonf(ibonf[6][x],ibonfn[6][x],data[0].n_loci,data[0].n_iter,file_output);
												
												printf("\n");
												if(file_output) fprintf(file_output,"\n");
										}
                                        }
                                    }
                                    else {
										printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        printf("\nloci\tObs(S)\tP(S)\tIBonf(S)\tObs(nhapl)\tP(nhapl)\tIBonf(nhapl)\tObs(nhaplsam)\tP(nhaplsam)\tIBonf(nhaplsam)\tObs(hapldiv)\tP(hapldiv)\tIBonf(hapldiv)\tObs(Rm)\tP(Rm)\tIBonf(Rm)\n");
                                        if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
										if(file_output) fputs("\nloci\tObs(S)\tP(S)\tIBonf(S)\tObs(nhapl)\tP(nhapl)\tIBonf(nhapl)\tObs(nhaplsam)\tP(nhaplsam)\tIBonf(nhaplsam)\tObs(hapldiv)\tP(hapldiv)\tIBonf(hapldiv)\tObs(Rm)\tP(Rm)\tIBonf(Rm)\n",file_output);
                                        for(x=0;x<data[0].n_loci;x++) {
                                            nunderalfa[0] = nunderalfa[1] = nunderalfa[2] = nunderalfa[3] = nunderalfa[4] = nunderalfa[5] = 0;
                                            nequivalfa[0] = nequivalfa[1] = nequivalfa[2] = nequivalfa[3] = nequivalfa[4] = nequivalfa[5] = 0;
                                            for(it=0;it<(long int)data[0].n_iter;it++) {
                                                if(matrix[x].biallsites < matrixsim[x][it].biallsites)
                                                    nunderalfa[0] += 1;
                                                if(matrix[x].nhapl < matrixsim[x][it].nhapl)
                                                    nunderalfa[1] += 1;
                                                if(matrix[x].nhaplsam < matrixsim[x][it].nhaplsam)
                                                    nunderalfa[2] += 1;
                                                if(matrix[x].hapldiv < matrixsim[x][it].hapldiv)
                                                    nunderalfa[3] += 1;
                                                if(matrix[x].Rm < matrixsim[x][it].Rm)
                                                    nunderalfa[4] += 1;
                                                if(matrix[x].biallsites == matrixsim[x][it].biallsites)
                                                    nequivalfa[0] += 1;
                                                if(matrix[x].nhapl == matrixsim[x][it].nhapl)
                                                    nequivalfa[1] += 1;
                                                if(matrix[x].nhaplsam == matrixsim[x][it].nhaplsam)
                                                    nequivalfa[2] += 1;
                                                if(matrix[x].hapldiv == matrixsim[x][it].hapldiv)
                                                    nequivalfa[3] += 1;
                                                 if(matrix[x].Rm == matrixsim[x][it].Rm)
                                                    nequivalfa[4] += 1;
                                            }
											ibonf[0][x] = (double)nunderalfa[0]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
											if(ibonf[0][x] > (double)0.5) ibonf[0][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
											else ibonf[0][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
											
											ibonf[1][x] = (double)nunderalfa[1]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
											if(ibonf[1][x] > (double)0.5) ibonf[1][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
											else ibonf[1][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
											
											ibonf[2][x] = (double)nunderalfa[2]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
											if(ibonf[2][x] > (double)0.5) ibonf[2][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
											else ibonf[2][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
											
											ibonf[3][x] = (double)nunderalfa[3]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
											if(ibonf[3][x] > (double)0.5) ibonf[3][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
											else ibonf[3][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
											
											ibonf[4][x] = (double)nunderalfa[4]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
											if(ibonf[4][x] > (double)0.5) ibonf[4][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
											else ibonf[4][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                        }
										ibonf_sort(ibonf,ibonfn,5,data[0].n_loci);
										
										for(x=0;x<data[0].n_loci;x++) {
											printf("%d:%s\t",x,matrix[x].gene);
											if(file_output) fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
											
											printf("%d\t%g\t",matrix[x].biallsites,ibonf[0][x]);
											if(file_output) fprintf(file_output,"%d\t%g\t",matrix[x].biallsites,ibonf[0][x]);
											print_ibonf(ibonf[0][x],ibonfn[0][x],data[0].n_loci,data[0].n_iter,file_output);
											
											printf("%d\t%g\t",matrix[x].nhapl,ibonf[1][x]);
											if(file_output) fprintf(file_output,"%d\t%g\t",matrix[x].nhapl,ibonf[1][x]);
											print_ibonf(ibonf[1][x],ibonfn[1][x],data[0].n_loci,data[0].n_iter,file_output);		
																					
											printf("%g\t%g\t",matrix[x].nhaplsam,ibonf[2][x]);
											if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].nhaplsam,ibonf[2][x]);
											print_ibonf(ibonf[2][x],ibonfn[2][x],data[0].n_loci,data[0].n_iter,file_output);		
																					
											printf("%g\t%g\t",matrix[x].hapldiv,ibonf[3][x]);
											if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].hapldiv,ibonf[3][x]);
											print_ibonf(ibonf[3][x],ibonfn[3][x],data[0].n_loci,data[0].n_iter,file_output);		
																					
											printf("%d\t%g\t",matrix[x].Rm,ibonf[4][x]);
											if(file_output) fprintf(file_output,"%d\t%g\t",matrix[x].Rm,ibonf[4][x]);
											print_ibonf(ibonf[4][x],ibonfn[4][x],data[0].n_loci,data[0].n_iter,file_output);	
																						
											printf("\n");
											if(file_output) fprintf(file_output,"\n");
                          				}
                                    }
                                }
                                else {
                                    printf("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n");
                                    if(file_output)
                                        fputs("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n",file_output);
                                }
                            }
							#if ALLRESULTS == 0
                            else {
							#endif
                                if(!onlymulo) {
                                    /*calculate median, 10% 5% etc./locus*/
                                    if((sortvector = (double *) malloc(data[0].n_iter * sizeof(double))) == 0) {
                                        printf("Probabilities can not be calculated, sorry.\n");
                                        if(file_output) fputs("Probabilities can not be calculated, sorry.\n",file_output);
                                        break;
                                    }
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n");
                                    printf(" The table contains in order the name of the locus,");
                                    printf(" the number of segregating biallelic sites,");
                                    if(data[0].time_spec > 0.) printf(" fixed mutations, divergence,");
                                    printf(" the number of haplotypes (not shared, if observed), haplotypes/sample size, haplotype diversity, and Rm.\n\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n",file_output);
                                        fprintf(file_output," The table contains in order the name of the locus,");
                                        fprintf(file_output," the number of segregating biallelic sites,");
                                        if(data[0].time_spec > 0.) fprintf(file_output," fixed mutations, divergence,");
                                        fprintf(file_output," the number of haplotypes (not shared, if observed), haplotypes/sample size, haplotype diversity, and Rm.\n\n");
                                    }
									if(*outgroup) {                                                
                                        if(data[0].time_spec <= 0.) {
                                            for(x=0;x<data[0].n_loci;x++) {
                                                printf("\nloci\tstatistic\tmedian\tavg\tvar\t");                                    
                                                if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                                if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                                if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                                if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                                if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                                if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                                if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                                printf("\n");
                                                if(file_output) {
                                                    fputs("\nloci\tstatistic\tmedian\tavg\tvar\t",file_output);                                    
                                                    if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                                    if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                                    if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                                    if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                                    if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                                    if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                                    if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                                    fprintf(file_output,"\n");
                                                }
                                                
                                                /*Sbiallelic sites. */
                                                /*average*/
                                                avg = var = (double)0.0;
                                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                    if(matrixsim[x][it].biallsites != -10000) {
                                                        sortvector[it2] = (double)matrixsim[x][it].biallsites;
                                                        it2 += 1;
                                                        avg += (double)matrixsim[x][it].biallsites;
                                                        var += (double)matrixsim[x][it].biallsites * (double)matrixsim[x][it].biallsites;
                                                    }
                                                }
                                                var = varianceit(var,avg,it2);
                                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                                printf("%d\tbialsites\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(file_output) 
                                                    fprintf(file_output,"%d\tbialsites\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(var != (double)-10000) {
                                                    printf("%g\t",var);
                                                    if(file_output) fprintf(file_output,"%g\t",var);
                                                }
                                                else {
                                                    printf("na\t");
                                                    if(file_output) fprintf(file_output,"na\t");                                                						}
                                                print_percentages(data,sortvector,it2,file_output);
                                                
												/*nhaplotypes*/
                                                /*average*/
                                                avg = var = (double)0.0;
                                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                    if( matrixsim[x][it].nhapl != -10000) {
                                                        sortvector[it2] = (double)matrixsim[x][it].nhapl;
                                                        it2 += 1;
                                                        avg += (double)matrixsim[x][it].nhapl;
                                                        var += (double)matrixsim[x][it].nhapl * (double)matrixsim[x][it].nhapl;
                                                    }
                                                }
                                                var = varianceit(var,avg,it2);
                                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                                printf("%d\tnhapl\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(file_output)
                                                    fprintf(file_output,"%d\tnhapl\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(var != (double)-10000) {
                                                    printf("%g\t",var);
                                                    if(file_output) fprintf(file_output,"%g\t",var);
                                                }
                                                else {
                                                    printf("na\t");
                                                    if(file_output) fprintf(file_output,"na\t");                                                }
                                                print_percentages(data,sortvector,it2,file_output);
                                               
												/*nhaplotypes/samplesize*/
                                                /*average*/
                                                avg = var = (double)0.0;
                                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                    if( matrixsim[x][it].nhaplsam != -10000) {
                                                        sortvector[it2] = (double)matrixsim[x][it].nhaplsam;
                                                        it2 += 1;
                                                        avg += (double)matrixsim[x][it].nhaplsam;
                                                        var += (double)matrixsim[x][it].nhaplsam * (double)matrixsim[x][it].nhaplsam;
                                                    }
                                                }
                                                var = varianceit(var,avg,it2);
                                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                                printf("%d\tnhaplsam\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(file_output)
                                                    fprintf(file_output,"%d\tnhaplsam\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(var != (double)-10000) {
                                                    printf("%g\t",var);
                                                    if(file_output) fprintf(file_output,"%g\t",var);
                                                }
                                                else {
                                                    printf("na\t");
                                                    if(file_output) fprintf(file_output,"na\t");
												}
                                                print_percentages(data,sortvector,it2,file_output);
                                                
												/*nhaplotype diversity*/
                                                /*average*/
                                                avg = var = (double)0.0;
                                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                    if( matrixsim[x][it].hapldiv != -10000) {
                                                        sortvector[it2] = (double)matrixsim[x][it].hapldiv;
                                                        it2 += 1;
                                                        avg += (double)matrixsim[x][it].hapldiv;
                                                        var += (double)matrixsim[x][it].hapldiv * (double)matrixsim[x][it].hapldiv;
                                                    }
                                                }
                                                var = varianceit(var,avg,it2);
                                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                                printf("%d\thapldiv\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(file_output)
                                                    fprintf(file_output,"%d\thapldiv\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(var != (double)-10000) {
                                                    printf("%g\t",var);
                                                    if(file_output) fprintf(file_output,"%g\t",var);
                                                }
                                                else {
                                                    printf("na\t");
                                                    if(file_output) fprintf(file_output,"na\t");
												}
                                                print_percentages(data,sortvector,it2,file_output);
                                                                                                
												/*Rm*/
                                                /*average*/
                                                avg = var = (double)0.0;
                                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                    if( matrixsim[x][it].Rm != -10000) {
                                                        sortvector[it2] = (double)matrixsim[x][it].Rvpi;
                                                        it2 += 1;
                                                        avg += (double)matrixsim[x][it].Rm;
                                                        var += (double)matrixsim[x][it].Rm * (double)matrixsim[x][it].Rm;
                                                    }
                                                }
                                                var = varianceit(var,avg,it2);
                                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                                printf("%d\tRm\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(file_output)
                                                    fprintf(file_output,"%d\tRm\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(var != (double)-10000) {
                                                    printf("%g\t",var);
                                                    if(file_output) fprintf(file_output,"%g\t",var);
                                                }
                                                else {
                                                    printf("na\t");
                                                    if(file_output) fprintf(file_output,"na\t");
												}
                                                print_percentages(data,sortvector,it2,file_output);
                                            }
                                        }
                                        else {
                                            for(x=0;x<data[0].n_loci;x++) {
                                                printf("\nloci\tstatistic\tmedian\tavg\tvar\t");                                    
                                                if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                                if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                                if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                                if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                                if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                                if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                                if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                                printf("\n");
                                                if(file_output) {
                                                    fputs("\nloci\tstatistic\tmedian\tavg\tvar\t",file_output);                                    
                                                    if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                                    if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                                    if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                                    if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                                    if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                                    if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                                    if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                                    fprintf(file_output,"\n");
                                                }
                                                /*Sbiallelic sites. */
                                                /*average*/
                                                avg = var = (double)0.0;
                                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                    if( matrixsim[x][it].biallsites != -10000) {
                                                        sortvector[it2] = (double)matrixsim[x][it].biallsites;
                                                        it2 += 1;
                                                        avg += (double)matrixsim[x][it].biallsites;
                                                        var += (double)matrixsim[x][it].biallsites * (double)matrixsim[x][it].biallsites;
                                                    }
                                                }
                                                var = varianceit(var,avg,it2);
                                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                                printf("%d\tbialsites\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(file_output)
                                                    fprintf(file_output,"%d\tbialsites\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(var != (double)-10000) {
                                                    printf("%g\t",var);
                                                    if(file_output) fprintf(file_output,"%g\t",var);
                                                }
                                                else {
                                                    printf("na\t");
                                                    if(file_output) fprintf(file_output,"na\t");                                                }
                                                print_percentages(data,sortvector,it2,file_output);
                                                
												/*fixed*/
                                                /*average*/
                                                avg = var = (double)0.0;
                                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                    if( matrixsim[x][it].fixed != -10000) {
                                                        sortvector[it2] = (double)matrixsim[x][it].fixed;
                                                        it2 += 1;
                                                        avg += (double)matrixsim[x][it].fixed;
                                                        var += (double)matrixsim[x][it].fixed * (double)matrixsim[x][it].fixed;
                                                    }
                                                }
                                                var = varianceit(var,avg,it2);
                                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                                printf("%d\tfixed\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(file_output)
                                                    fprintf(file_output,"%d\tfixed\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(var != (double)-10000) {
                                                    printf("%g\t",var);
                                                    if(file_output) fprintf(file_output,"%g\t",var);
                                                }
                                                else {
                                                    printf("na\t");
                                                    if(file_output) fprintf(file_output,"na\t");                                                						}
                                                print_percentages(data,sortvector,it2,file_output);
                                                
												/*divergence*/
                                                /*average*/
                                                avg = var = (double)0.0;
                                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                    if( matrixsim[x][it].ndivergence != (double) -10000) {
                                                        sortvector[it2] = (double)matrixsim[x][it].ndivergence;
                                                        it2 += 1;
                                                        avg += (double)matrixsim[x][it].ndivergence;
                                                        var += (double)matrixsim[x][it].ndivergence * (double)matrixsim[x][it].ndivergence;
                                                    }
                                                }
                                                var = varianceit(var,avg,it2);
                                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                                printf("%d\tDiverg\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(file_output)
                                                    fprintf(file_output,"%d\tDiverg\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(var != (double)-10000) {
                                                    printf("%g\t",var);
                                                    if(file_output) fprintf(file_output,"%g\t",var);
                                                }
                                                else {
                                                    printf("na\t");
                                                    if(file_output) fprintf(file_output,"na\t");                                                }
                                                print_percentages(data,sortvector,it2,file_output);
                                                
												/*nhaplotypes*/
                                                /*average*/
                                                avg = var = (double)0.0;
                                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                    if( matrixsim[x][it].nhapl != -10000) {
                                                        sortvector[it2] = (double)matrixsim[x][it].nhapl;
                                                        it2 += 1;
                                                        avg += (double)matrixsim[x][it].nhapl;
                                                        var += (double)matrixsim[x][it].nhapl * (double)matrixsim[x][it].nhapl;
                                                    }
                                                }
                                                var = varianceit(var,avg,it2);
                                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                                printf("%d\tnhapl\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(file_output)
                                                    fprintf(file_output,"%d\tnhapl\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(var != (double)-10000) {
                                                    printf("%g\t",var);
                                                    if(file_output) fprintf(file_output,"%g\t",var);
                                                }
                                                else {
                                                    printf("na\t");
                                                    if(file_output) fprintf(file_output,"na\t");                                                						}
                                                print_percentages(data,sortvector,it2,file_output);
                                                
												/*nhaplotypes/samplesize*/
                                                /*average*/
                                                avg = var = (double)0.0;
                                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                    if( matrixsim[x][it].nhaplsam != -10000) {
                                                        sortvector[it2] = (double)matrixsim[x][it].nhaplsam;
                                                        it2 += 1;
                                                        avg += (double)matrixsim[x][it].nhaplsam;
                                                        var += (double)matrixsim[x][it].nhaplsam * (double)matrixsim[x][it].nhaplsam;
                                                    }
                                                }
                                                var = varianceit(var,avg,it2);
                                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                                printf("%d\tnhaplsam\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(file_output)
                                                    fprintf(file_output,"%d\tnhaplsam\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(var != (double)-10000) {
                                                    printf("%g\t",var);
                                                    if(file_output) fprintf(file_output,"%g\t",var);
                                                }
                                                else {
                                                    printf("na\t");
                                                    if(file_output) fprintf(file_output,"na\t");                                                						}
                                                print_percentages(data,sortvector,it2,file_output);
                                                
												/*nhaplotype diversity*/
                                                /*average*/
                                                avg = var = (double)0.0;
                                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                    if( matrixsim[x][it].hapldiv != -10000) {
                                                        sortvector[it2] = (double)matrixsim[x][it].hapldiv;
                                                        it2 += 1;
                                                        avg += (double)matrixsim[x][it].hapldiv;
                                                        var += (double)matrixsim[x][it].hapldiv * (double)matrixsim[x][it].hapldiv;
                                                    }
                                                }
                                                var = varianceit(var,avg,it2);
                                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                                printf("%d\thapldiv\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(file_output)
                                                    fprintf(file_output,"%d\thapldiv\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(var != (double)-10000) {
                                                    printf("%g\t",var);
                                                    if(file_output) fprintf(file_output,"%g\t",var);
                                                }
                                                else {
                                                    printf("na\t");
                                                    if(file_output) fprintf(file_output,"na\t");                                                						}
                                                print_percentages(data,sortvector,it2,file_output);
                                                                                                
												/*Rm*/
                                                /*average*/
                                                avg = var = (double)0.0;
                                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                    if( matrixsim[x][it].Rm != -10000) {
                                                        sortvector[it2] = (double)matrixsim[x][it].Rm;
                                                        it2 += 1;
                                                        avg += (double)matrixsim[x][it].Rm;
                                                        var += (double)matrixsim[x][it].Rm * (double)matrixsim[x][it].Rm;
                                                    }
                                                }
                                                var = varianceit(var,avg,it2);
                                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                                printf("%d\tRm\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(file_output)
                                                    fprintf(file_output,"%d\tRm\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                                if(var != (double)-10000) {
                                                    printf("%g\t",var);
                                                    if(file_output) fprintf(file_output,"%g\t",var);
                                                }
                                                else {
                                                    printf("na\t");
                                                    if(file_output) fprintf(file_output,"na\t");                                                						}
                                                print_percentages(data,sortvector,it2,file_output);
                                            }
                                        }
                                    }
                                    else {
                                        for(x=0;x<data[0].n_loci;x++) {
                                            printf("\nloci\tstatistic\tmedian\tavg\tvar\t");                                    
                                            if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                            if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                            if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                            if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                            if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                            if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                            if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                            printf("\n");
                                            if(file_output) {
                                                fputs("\nloci\tstatistic\tmedian\tavg\tvar\t",file_output);                                    
                                                if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                                if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                                if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                                if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                                if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                                if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                                if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                                fprintf(file_output,"\n");
                                            }
                                            /*Sbiallelic sites. */
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].biallsites != -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].biallsites;
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].biallsites;
                                                    var += (double)matrixsim[x][it].biallsites * (double)matrixsim[x][it].biallsites;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\tbialsites\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\tbialsites\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*nhapl*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].nhapl != -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].nhapl;
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].nhapl;
                                                    var += (double)matrixsim[x][it].nhapl * (double)matrixsim[x][it].nhapl;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\tnhapl\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\tnhapl\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*nhapl/samsize*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].nhaplsam != -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].nhaplsam;
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].nhaplsam;
                                                    var += (double)matrixsim[x][it].nhaplsam * (double)matrixsim[x][it].nhaplsam;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\tnhaplsam\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\tnhaplsam\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*nhapl diversity*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].hapldiv != -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].hapldiv;
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].hapldiv;
                                                    var += (double)matrixsim[x][it].hapldiv * (double)matrixsim[x][it].hapldiv;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\thapldiv\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\thapldiv\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
																						
											/*Rm*/
											/*average*/
											avg = var = (double)0.0;
											for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
												if( matrixsim[x][it].Rm != -10000) {
													sortvector[it2] = (double)matrixsim[x][it].Rm;
													it2 += 1;
													avg += (double)matrixsim[x][it].Rm;
													var += (double)matrixsim[x][it].Rm * (double)matrixsim[x][it].Rm;
												}
											}
											var = varianceit(var,avg,it2);
											if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
											printf("%d\rRm\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
											if(file_output)
												fprintf(file_output,"%d\tRm\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
											if(var != (double)-10000) {
												printf("%g\t",var);
												if(file_output) fprintf(file_output,"%g\t",var);
											}
											else {
												printf("na\t");
												if(file_output) fprintf(file_output,"na\t");                                                }
											print_percentages(data,sortvector,it2,file_output);
                                        }
                                    }
                                    free(sortvector);
                                }
                                else {
                                    printf("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n");
                                    if(file_output)
                                        fputs("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n",file_output);
                                }
							#if ALLRESULTS == 0
                            }
							#endif
                            break;
                        case '1': /* 1 - Display linkage related statistics.*/
                            if(*observed_data && dataobsequalsim) {/*compare simulations with observed data.*/
                                /*Sign tests.*/
                                /*Binomial for each statistic. from median (all data) or average (onlymulo). outgroup and non outgroup.*/
                                if((psigned = (int *)calloc(3,sizeof(int))) == 0) {
                                    printf("Error in memory allocation: Sign tests are not available.");
                                    if(file_output) fputs("Error in memory allocation: Sign tests are not available.",file_output);
                                    break;
                                }
                                if((nsigned = (int *)calloc(3,sizeof(int))) == 0) {
                                    printf("Error in memory allocation: Sign tests are not available.");
                                    if(file_output) fputs("Error in memory allocation: Sign tests are not available.",file_output);
                                    break;
                                }
    
                                if(onlymulo) {
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\nNon-parametric SIGN TEST. Assuming that average have equivalent value than median. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl)");
                                    printf("\nComparing observed data with the average of simulated data at each locus.");
                                    printf("\nAdjust to a Binomial.\n");
                                    printf(" The table contains in order the name of the locus");
                                    printf(" Wall's not divided statistics B and Q,\n");
                                    printf("and Rozas' et al. not divided ZA (shared excluded in all three statistics, if observed)\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\nNon-parametric SIGN TEST. Assuming that average have equivalent value than median. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl)",file_output);
                                        fputs("\nComparing observed data with the average of simulated data at each locus.",file_output);
                                        fputs("\nAdjust to a Binomial.\n",file_output);
                                        fprintf(file_output," The table contains in order the name of the locus");
                                        fprintf(file_output," Wall's not divided statistics B and Q,\n");
                                        fprintf(file_output,"and Rozas' et al. not divided ZA (shared excluded in all three statistics, if observed)\n");
                                    }
                                }
                                else {
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\nNon-parametric SIGN TEST. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl).");
                                    printf("\nComparing observed data with the median of simulated data at each locus.");
                                    printf("\nAdjust to a Binomial.\n");
                                    printf(" The table contains in order the name of the locus");
                                    printf(" Wall's not divided statistics B and Q,\n");
                                    printf("and Rozas' et al. not divided ZA (shared excluded in all three statistics, if observed)\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\nNon-parametric SIGN TEST. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl).",file_output);
                                        fputs("\nComparing observed data with the median of simulated data at each locus.",file_output);
                                        fputs("\nAdjust to a Binomial.\n",file_output);
                                        fprintf(file_output," The table contains in order the name of the locus");
                                        fprintf(file_output," Wall's not divided statistics B and Q,\n");
                                        fprintf(file_output,"and Rozas' et al. not divided ZA (shared excluded in all three statistics, if observed)\n");
                                    }
                                }
                                printf("\nloci\tObs(b)\tExp(b)\tObs(q)\tExp(q)\tObs(za)\tExp(za)\n");
                                if(file_output)
                                    fputs("\nloci\tObs(b)\tExp(b)\tObs(q)\tExp(q)\tObs(za)\tExp(za)\n",file_output);
                                for(x=0;x<data[0].n_loci;x++) {
                                    printf("%d:%s\t%d\t%g\t%d\t%g\t%g\t%g\n",x,matrix[x].gene,
                                        matrix[x].b,avgstatloci[x].b,
                                        matrix[x].q,avgstatloci[x].q,
                                        matrix[x].za,avgstatloci[x].za);
                                    if(file_output)
                                        fprintf(file_output,"%d:%s\t%d\t%g\t%d\t%g\t%g\t%g\n",x,matrix[x].gene,
                                            matrix[x].b,avgstatloci[x].b,
                                            matrix[x].q,avgstatloci[x].q,
                                            matrix[x].za,avgstatloci[x].za);
    
                                    if(matrix[x].b > avgstatloci[x].b) psigned[0] += 1;
                                    if(matrix[x].q > avgstatloci[x].q) psigned[1] += 1;
                                    if(matrix[x].za > avgstatloci[x].za) psigned[2] += 1;
    
                                    if(matrix[x].b < avgstatloci[x].b) nsigned[0] += 1;
                                    if(matrix[x].q < avgstatloci[x].q) nsigned[1] += 1;
                                    if(matrix[x].za < avgstatloci[x].za) nsigned[2] += 1;
                                }
								printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
								if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                for(x=0;x<3;x++) {
                                    switch(x) {
                                        case 0:
                                            printf("\nSign test: Value of b (Wall's B not divided by S).\n");
                                            if(file_output) fputs("\nSign test: Value of b (Wall's B not divided by S).\n",file_output);
                                            break;
                                        case 1:
                                            printf("\nSign test: Value of q (Wall's Q not divided by S).\n");
                                            if(file_output) fputs("\nSign test: Value of q (Wall's Q not divided by S).\n",file_output);
                                            break;
                                        case 2:
                                            printf("\nSign test: Value of za (Rozas' not divided by S).\n");
                                            if(file_output) fputs("\nSign test: Value of za (Rozas' not divided by S).\n",file_output);
                                            break;
                                    }
                                    prob = prob_cumbinomial(psigned[x]+nsigned[x],psigned[x],0.5);
									if(prob == -10000) printf("#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
                                    else printf("#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
									if(file_output) {
										if(prob == -10000) fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
										else fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
									}
                                    printfignif1(prob,file_output);

									/*look at the effect of discrete values*/
									if(prob <= 0.05 && prob >= 0.) {
										eq = data[0].n_loci - (psigned[x]+nsigned[x]);
										if(psigned[x] > nsigned[x]) {
											ns = nsigned[x] + eq;
											ps = psigned[x];
										}
										else {
											ns = nsigned[x];
											ps = psigned[x] + eq;
										}
										if((psigned[x] > nsigned[x] && ps < ns) || (psigned[x] < nsigned[x] && ps > ns) || prob_cumbinomial(ps+ns,ps,0.5) > 0.05) {
											printf(" WARNING: too many values equal to the median. It might be too liberal.");
											if(file_output) fputs(" WARNING: too many values equal to the median. It might be too liberal.",file_output);
										}
									}
								}
                                free(psigned);
                                free(nsigned);
                                
                                if(!onlymulo) {
                                    /*calculate probability each observed value/locus*/
                                    printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									printf("\n\n Display the probability of observed statistics for each locus obtained by coalescent Monte Carlo simulations: ");
                                    printf("\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. Note: & means the ierations would be increased to confirm significance");
                                    printf("\n Multiple hits not included in the analysis.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\n\n Display the probability of observed statistics for each locus obtained by coalescent Monte Carlo simulations: ",file_output);
                                        fprintf(file_output,"\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. Note: & means the ierations would be increased to confirm significance");
                                        fputs("\n Multiple hits not included in the analysis.\n",file_output);
                                    }
                                    printf("\nloci\tObs(za)\tP(za)\tIBonf(za)\tObs(b)\tP(b)\tIBonf(b)\tObs(q)\tP(q)\tIBonf(q)\n");
                                    if(file_output) fputs("\nloci\tObs(za)\tP(za)\tIBonf(za)\tObs(b)\tP(b)\tIBonf(b)\tObs(q)\tP(q)\tIBonf(q)\n",file_output);
                                    for(x=0;x<data[0].n_loci;x++) {
                                        nunderalfa[0] = nunderalfa[1] = nunderalfa[2] = 0;
                                        nequivalfa[0] = nequivalfa[1] = nequivalfa[2] = 0;
                                        for(it=0;it<(long int)data[0].n_iter;it++) {
                                            if(matrix[x].za < matrixsim[x][it].za)
                                                nunderalfa[0] += 1;
                                            if(matrix[x].b < matrixsim[x][it].b)
                                                nunderalfa[1] += 1;
                                            if(matrix[x].q < matrixsim[x][it].q)
                                                nunderalfa[2] += 1;
                                            if(matrix[x].za == matrixsim[x][it].za)
                                                nunderalfa[0] += 1;
                                            if(matrix[x].b == matrixsim[x][it].b)
                                                nunderalfa[1] += 1;
                                            if(matrix[x].q == matrixsim[x][it].q)
                                                nunderalfa[2] += 1;
                                        }
                                        ibonf[0][x] = (double)nunderalfa[0]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                        if(ibonf[0][x] > (double)0.5) ibonf[0][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                        else ibonf[0][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                        
										ibonf[1][x] = (double)nunderalfa[1]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                        if(ibonf[1][x] > (double)0.5) ibonf[1][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                        else ibonf[1][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                        
										ibonf[2][x] = (double)nunderalfa[2]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                        if(ibonf[2][x] > (double)0.5) ibonf[2][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                        else ibonf[2][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                    }
									ibonf_sort(ibonf,ibonfn,3,data[0].n_loci);

									for(x=0;x<data[0].n_loci;x++) {
                                        printf("%d:%s\t",x,matrix[x].gene);
                                        if(file_output) fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
                                        
										printf("%g\t%g\t",matrix[x].za,ibonf[0][x]);
 										if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].za,ibonf[0][x]);
										print_ibonf(ibonf[0][x],ibonfn[0][x],data[0].n_loci,data[0].n_iter,file_output);
                                        
										printf("%d\t%g\t",matrix[x].b,ibonf[1][x]);
										if(file_output) fprintf(file_output,"%d\t%g\t",matrix[x].b,ibonf[1][x]);
										print_ibonf(ibonf[1][x],ibonfn[1][x],data[0].n_loci,data[0].n_iter,file_output);
                                        
										printf("%d\t%g\t",matrix[x].q,ibonf[2][x]);
  										if(file_output) fprintf(file_output,"%d\t%g\t",matrix[x].q,ibonf[2][x]);
										print_ibonf(ibonf[2][x],ibonfn[2][x],data[0].n_loci,data[0].n_iter,file_output);
										
										printf("\n");
										if(file_output) fprintf(file_output,"\n");
									}
                                }
                                else {
                                    printf("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n");
                                    if(file_output)
                                        fputs("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n",file_output);
                                }
                            }
							#if ALLRESULTS == 0
                            else {
							#endif
                                if(!onlymulo) {
                                    /*calculate median, 10% 5% etc./locus*/
                                    if((sortvector = (double *) malloc(data[0].n_iter * sizeof(double))) == 0) {
                                        printf("Probabilities can not be calculated, sorry.\n");
                                        if(file_output) fputs("Probabilities can not be calculated, sorry.\n",file_output);
                                        break;
                                    }
                                    printf("\n\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n");
                                    if(file_output)
                                        fputs("\n\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n",file_output);
                                    for(x=0;x<data[0].n_loci;x++) {
                                        printf("\nloci\tstatistic\tmedian\tavg\tvar\t");                                    
                                        if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                        if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                        if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                        if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                        if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                        if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                        if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                        printf("\n");
                                        if(file_output) {
                                            fputs("\nloci\tstatistic\tmedian\tavg\tvar\t",file_output);                                    
                                            if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                            if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                            if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                            if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                            if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                            if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                            if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                            fprintf(file_output,"\n");
                                        }
                                        /*za. */
                                        /*average*/
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if( matrixsim[x][it].za != (double) -10000) {
                                                sortvector[it2] = (double)matrixsim[x][it].za;
                                                it2 += 1;
                                                avg += (double)matrixsim[x][it].za;
                                                var += (double)matrixsim[x][it].za * (double)matrixsim[x][it].za;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tza\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tza\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                        
										/*b*/
                                        /*average*/
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if( matrixsim[x][it].b != (double)-10000) {
                                                sortvector[it2] = (double)matrixsim[x][it].b;
                                                it2 += 1;
                                                avg += (double)matrixsim[x][it].b;
                                                var += (double)matrixsim[x][it].b * (double)matrixsim[x][it].b;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tb\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tb\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                        
										/*q*/
                                        /*average*/
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if( matrixsim[x][it].q != (double)-10000) {
                                                sortvector[it2] = (double)matrixsim[x][it].q;
                                                it2 += 1;
                                                avg += (double)matrixsim[x][it].q;
                                                var += (double)matrixsim[x][it].q * (double)matrixsim[x][it].q;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tq\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tq\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    free(sortvector);
                                }
                                else {
                                    printf("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n");
                                    if(file_output)
                                        fputs("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n",file_output);
                                }
							#if ALLRESULTS == 0
                            }
							#endif
                            break;
                        case '2': /* 2 - Display Theta/locus.*/
                            if(*observed_data && dataobsequalsim) {/*compare simulations with observed data.*/
                                /*Sign tests.*/
                                /*Binomial for each statistic. from median (all data) or average (onlymulo). outgroup and non outgroup.*/
                                if((psigned = (int *)calloc(5,sizeof(int))) == 0) {
                                    printf("Error in memory allocation: Sign tests are not available.");
                                    if(file_output) fputs("Error in memory allocation: Sign tests are not available.",file_output);
                                    break;
                                }
                                if((nsigned = (int *)calloc(5,sizeof(int))) == 0) {
                                    printf("Error in memory allocation: Sign tests are not available.");
                                    if(file_output) fputs("Error in memory allocation: Sign tests are not available.",file_output);
                                    break;
                                }

                                if(onlymulo) {
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									printf(" NOTE that all theta estimates are corrected for chromosome population size.\n");
                                    printf("\nNon-parametric SIGN TEST. Assuming that average have equivalent value than median. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl)");
                                    printf("\nComparing observed data with the average of simulated data at each locus.");
                                    printf("\nAdjust to a Binomial.\n");
                                    printf("The table contains in order the name of the locus");
                                    printf(" and locus estimations of variability from Watterson, Tajima, Fu and Li,\n");
                                    printf("Fay and Wu and Zeng et al.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
										fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n",file_output);
                                        fputs("\nNon-parametric SIGN TEST. Assuming that average have equivalent value than median. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl)",file_output);
                                        fputs("\nComparing observed data with the average of simulated data at each locus.",file_output);
                                        fputs("\nAdjust to a Binomial.\n",file_output);
                                        fprintf(file_output,"The table contains in order the name of the locus");
                                        fprintf(file_output," and locus estimations of variability from Watterson, Tajima, Fu and Li,\n");
                                        fprintf(file_output,"Fay and Wu and Zeng et al.\n");
                                    }
                                }
                                else {
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf(" NOTE that all theta estimates are corrected for chromosome population size.\n");
									printf("\nNon-parametric SIGN TEST. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl).");
                                    printf("\nComparing observed data with the median of simulated data at each locus.");
                                    printf("\nAdjust to a Binomial.\n");
                                    printf("The table contains in order the name of the locus");
                                    printf(" and locus estimations of variability from Watterson, Tajima, Fu and Li,\n");
                                    printf("Fay and Wu and Zeng et al.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
										fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n",file_output);
                                        fputs("\nNon-parametric SIGN TEST. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl).",file_output);
                                        fputs("\nComparing observed data with the median of simulated data at each locus.",file_output);
                                        fputs("\nAdjust to a Binomial.\n",file_output);
                                        fprintf(file_output,"The table contains in order the name of the locus");
                                        fprintf(file_output," and locus estimations of variability from Watterson, Tajima, Fu and Li,\n");
                                        fprintf(file_output,"Fay and Wu and Zeng et al.\n");
                                    }
                                }
                                if(*outgroup) {
                                    printf("\nloci\tObs(theta_wat)\tExp(theta_wat)\tObs(theta_taj)\tExp(theta_taj)\tObs(theta_fuli)\tExp(theta_fuli)\tObs(theta_fw)\tExp(theta_fw)\tObs(theta_zeng)\tExp(theta_zeng)\n");
                                    if(file_output)
                                        fputs("\nloci\tObs(theta_wat)\tExp(theta_wat)\tObs(theta_taj)\tExp(theta_taj)\tObs(theta_fuli)\tExp(theta_fuli)\tObs(theta_fw)\tExp(theta_fw)\tObs(theta_zeng)\tExp(theta_zeng)\n",file_output);
                                    for(x=0;x<data[0].n_loci;x++) {
                                        printf("%d:%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",x,matrix[x].gene,
                                            matrix[x].theta_wat*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_wat,
                                            matrix[x].theta_taj*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_taj,
                                            matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_fuli,
                                            matrix[x].theta_fw*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_fw,
											matrix[x].theta_L*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_L);
                                        if(file_output)
                                            fprintf(file_output,"%d:%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",x,matrix[x].gene,
                                                matrix[x].theta_wat*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_wat,
                                                matrix[x].theta_taj*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_taj,
                                                matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_fuli,
                                                matrix[x].theta_fw*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_fw,
												matrix[x].theta_L*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_L);

                                        if(matrix[x].theta_wat*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_wat) psigned[0] += 1;
                                        if(matrix[x].theta_taj*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_taj) psigned[1] += 1;
                                        if(matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_fuli) psigned[2] += 1;
                                        if(matrix[x].theta_fw*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_fw) psigned[3] += 1;
                                        if(matrix[x].theta_L*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_L) psigned[4] += 1;

                                        if(matrix[x].theta_wat*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_wat) nsigned[0] += 1;
                                        if(matrix[x].theta_taj*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_taj) nsigned[1] += 1;
                                        if(matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_fuli) nsigned[2] += 1;
                                        if(matrix[x].theta_fw*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_fw) nsigned[3] += 1;
										if(matrix[x].theta_L*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_L) nsigned[4] += 1;
                                    }
                                    printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									for(x=0;x<5;x++) {
                                        switch(x) {
                                            case 0:
                                                printf("\nSign test: Theta (Watterson) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Watterson) corrected for chromosome population size.\n",file_output);
                                                break;
                                            case 1:
                                                printf("\nSign test: Theta (Tajima) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Tajima) corrected for chromosome population size.\n",file_output);
                                                break;
                                            case 2:
                                                printf("\nSign test: Theta (Fu and Li) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Fu and Li) corrected for chromosome population size.\n",file_output);
                                                break;
                                            case 3:
                                                printf("\nSign test: Theta (Fay and Wu) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Fay and Wu) corrected for chromosome population size.\n",file_output);
                                                break;
                                            case 4:
                                                printf("\nSign test: Theta (Zeng et al.) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Zeng et al.) corrected for chromosome population size.\n",file_output);
                                                break;
                                        }
                                        prob = prob_cumbinomial(psigned[x]+nsigned[x],psigned[x],0.5);
										if(prob == -10000) printf("#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
                                        else printf("#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
										if(file_output) {
											if(prob == -10000) fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
											else fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
										}
                                        printfignif1(prob,file_output);

										/*look at the effect of discrete values*/
										if(prob <= 0.05 && prob >= 0.) {
											eq = data[0].n_loci - (psigned[x]+nsigned[x]);
											if(psigned[x] > nsigned[x]) {
												ns = nsigned[x] + eq;
												ps = psigned[x];
											}
											else {
												ns = nsigned[x];
												ps = psigned[x] + eq;
											}
											if((psigned[x] > nsigned[x] && ps < ns) || (psigned[x] < nsigned[x] && ps > ns) || prob_cumbinomial(ps+ns,ps,0.5) > 0.05) {
                                                printf(" WARNING: too many values equal to the median. It might be too liberal.");
                                                if(file_output) fputs(" WARNING: too many values equal to the median. It might be too liberal.",file_output);
											}
										}
									}
                                }
                                else {
									puts("\nloci\tObs(theta_wat)\tExp(theta_wat)\tObs(theta_taj)\tExp(theta_taj)\tObs(theta_fuli)\tExp(theta_fuli)\n");
                                    if(file_output)
                                        fputs("\nloci\tObs(theta_wat)\tExp(theta_wat)\tObs(theta_taj)\tExp(theta_taj)\tObs(theta_fuli)\tExp(theta_fuli)\n",file_output);
                                    for(x=0;x<data[0].n_loci;x++) {
                                        printf("%d:%s\t%g\t%g\t%g\t%g\t%g\t%g\n",x,matrix[x].gene,
                                            matrix[x].theta_wat*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_wat,
                                            matrix[x].theta_taj*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_taj,
                                            matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_fulin);
                                        if(file_output)
                                            fprintf(file_output,"%d:%s\t%g\t%g\t%g\t%g\t%g\t%g\n",x,matrix[x].gene,
                                                matrix[x].theta_wat*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_wat,
                                                matrix[x].theta_taj*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_taj,
                                                matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_fulin);

                                        if(matrix[x].theta_wat*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_wat) psigned[0] += 1;
                                        if(matrix[x].theta_taj*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_taj) psigned[1] += 1;
                                        if(matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_fulin) psigned[2] += 1;

                                        if(matrix[x].theta_wat*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_wat) nsigned[0] += 1;
                                        if(matrix[x].theta_taj*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_taj) nsigned[1] += 1;
                                        if(matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_fulin) nsigned[2] += 1;
                                    }
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    for(x=0;x<3;x++) {
                                        switch(x) {
                                            case 0:
                                                printf("\nSign test: Theta (Watterson) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Watterson) corrected for chromosome population size.\n",file_output);
                                                break;
                                            case 1:
                                                printf("\nSign test: Theta (Tajima) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Tajima) corrected for chromosome population size.\n",file_output);
                                                break;
                                            case 2:
                                                printf("\nSign test: Theta (Fu and Li) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Fu and Li) corrected for chromosome population size.\n",file_output);
                                                break;
                                        }
                                        prob = prob_cumbinomial(psigned[x]+nsigned[x],psigned[x],0.5);
										if(prob == -10000) printf("#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
                                        else printf("#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
										if(file_output) {
											if(prob == -10000) fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
											else fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
										}
                                        printfignif1(prob,file_output);

										/*look at the effect of discrete values*/
										if(prob <= 0.05 && prob >= 0.) {
											eq = data[0].n_loci - (psigned[x]+nsigned[x]);
											if(psigned[x] > nsigned[x]) {
												ns = nsigned[x] + eq;
												ps = psigned[x];
											}
											else {
												ns = nsigned[x];
												ps = psigned[x] + eq;
											}
											if((psigned[x] > nsigned[x] && ps < ns) || (psigned[x] < nsigned[x] && ps > ns) || prob_cumbinomial(ps+ns,ps,0.5) > 0.05) {
                                                printf(" WARNING: too many values equal to the median. It might be too liberal.");
                                                if(file_output) fputs(" WARNING: too many values equal to the median. It might be too liberal.",file_output);
											}
										}
									}
                                }
                                free(psigned);
                                free(nsigned);
                                
                                if(!onlymulo) {
                                    /*calculate probability each observed value/locus*/
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\n\n Display the probability of observed statistics for each locus obtained by coalescent Monte Carlo simulations: ");
                                    printf("\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. Note: & means the ierations would be increased to confirm significance");
                                    printf("\n Multiple hits not included in the analysis.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\n\n Display the probability of observed statistics for each locus obtained by coalescent Monte Carlo simulations: ",file_output);
                                        fprintf(file_output,"\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. Note: & means the ierations would be increased to confirm significance");
                                        fputs("\n Multiple hits not included in the analysis.\n",file_output);
                                    }
                                    if(*outgroup) {
                                        printf("\nloci\tObs(theta_wat)\tP(theta_wat)\tIBonf(theta_wat)\tObs(theta_taj)\tP(theta_taj)\tIBonf(theta_taj)\tObs(theta_fuli)\tP(theta_fuli)\tIBonf(theta_fuli)\tObs(theta_fw)\tP(theta_fw)\tIBonf(theta_fw)\tP(theta_zeng)\tIBonf(theta_zeng)\n");
                                        if(file_output) 
                                            fputs("\nloci\tObs(theta_wat)\tP(theta_wat)\tIBonf(theta_wat)\tObs(theta_taj)\tP(theta_taj)\tIBonf(theta_taj)\tObs(theta_fuli)\tP(theta_fuli)\tIBonf(theta_fuli)\tObs(theta_fw)\tP(theta_fw)\tIBonf(theta_fw)\tP(theta_zeng)\tIBonf(theta_zeng)\n",file_output);
                                        for(x=0;x<data[0].n_loci;x++) {
                                            nunderalfa[0] = nunderalfa[1] = nunderalfa[2] = nunderalfa[3] = nunderalfa[4] = 0;
                                            nequivalfa[0] = nequivalfa[1] = nequivalfa[2] = nequivalfa[3] = nequivalfa[4] = 0;
                                            for(it=0;it<(long int)data[0].n_iter;it++) {
                                                if(matrix[x].theta_wat*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_wat*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[0] += 1;
                                                 if(matrix[x].theta_taj*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_taj*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[1] += 1;
                                                if(matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_fuli*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[2] += 1;
                                                if(matrix[x].theta_fw*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_fw*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[3] += 1;
                                                if(matrix[x].theta_L*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_L*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[4] += 1;
												if(matrix[x].theta_wat*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_wat*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[0] += 1;
                                                 if(matrix[x].theta_taj*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_taj*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[1] += 1;
                                                if(matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_fuli*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[2] += 1;
                                                if(matrix[x].theta_fw*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_fw*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[3] += 1;
                                                if(matrix[x].theta_L*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_L*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[4] += 1;
                                            }
                                            ibonf[0][x] = (double)nunderalfa[0]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                            if(ibonf[0][x] > (double)0.5) ibonf[0][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                            else ibonf[0][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                            
											ibonf[1][x] = (double)nunderalfa[1]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                            if(ibonf[1][x] > (double)0.5) ibonf[1][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                            else ibonf[1][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                            
											ibonf[2][x] = (double)nunderalfa[2]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                            if(ibonf[2][x] > (double)0.5) ibonf[2][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                            else ibonf[2][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                            
											ibonf[3][x] = (double)nunderalfa[3]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                            if(ibonf[3][x] > (double)0.5) ibonf[3][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                            else ibonf[3][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                            
											ibonf[4][x] = (double)nunderalfa[4]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                            if(ibonf[4][x] > (double)0.5) ibonf[4][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                            else ibonf[4][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                        }
										ibonf_sort(ibonf,ibonfn,5,data[0].n_loci);
										
                                        for(x=0;x<data[0].n_loci;x++) {
                                            printf("%d:%s\t",x,matrix[x].gene);
                                            if(file_output) fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
                                            
											printf("%g\t%g\t",matrix[x].theta_wat*((double)1/matrix[x].factor_chrn),ibonf[0][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_wat*((double)1/matrix[x].factor_chrn),ibonf[0][x]);
											print_ibonf(ibonf[0][x],ibonfn[0][x],data[0].n_loci,data[0].n_iter,file_output);
                                            
											printf("%g\t%g\t",matrix[x].theta_taj*((double)1/matrix[x].factor_chrn),ibonf[1][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_taj*((double)1/matrix[x].factor_chrn),ibonf[1][x]);
											print_ibonf(ibonf[1][x],ibonfn[1][x],data[0].n_loci,data[0].n_iter,file_output);
                                            
											printf("%g\t%g\t",matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn),ibonf[2][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn),ibonf[2][x]);
											print_ibonf(ibonf[2][x],ibonfn[2][x],data[0].n_loci,data[0].n_iter,file_output);
                                            
											printf("%g\t%g\t",matrix[x].theta_fw*((double)1/matrix[x].factor_chrn),ibonf[3][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_fw*((double)1/matrix[x].factor_chrn),ibonf[3][x]);
											print_ibonf(ibonf[3][x],ibonfn[3][x],data[0].n_loci,data[0].n_iter,file_output);
                                            
											printf("%g\t%g\t",matrix[x].theta_L*((double)1/matrix[x].factor_chrn),ibonf[4][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_L*((double)1/matrix[x].factor_chrn),ibonf[4][x]);
											print_ibonf(ibonf[4][x],ibonfn[4][x],data[0].n_loci,data[0].n_iter,file_output);
											
											printf("\n");
											if(file_output) fprintf(file_output,"\n");
										}
                                    }
                                    else {
                                        printf("\nloci\tObs(theta_wat)\tP(theta_wat)\tIBonf(theta_wat)\tObs(theta_taj)\tP(theta_taj)\tIBonf(theta_taj)\tObs(theta_fulin)\tP(theta_fulin)\tIBonf(theta_fulin)\n");
                                        if(file_output) 
                                            fputs("\nloci\tObs(theta_wat)\tP(theta_wat)\tIBonf(theta_wat)\tObs(theta_taj)\tP(theta_taj)\tIBonf(theta_taj)\tObs(theta_fulin)\tP(theta_fulin)\tIBonf(theta_fulin)\n",file_output);
                                        for(x=0;x<data[0].n_loci;x++) {
                                            nunderalfa[0] = nunderalfa[1] = nunderalfa[2] = 0;
                                            nequivalfa[0] = nequivalfa[1] = nequivalfa[2] = 0;
                                            for(it=0;it<(long int)data[0].n_iter;it++) {
                                                if(matrix[x].theta_wat*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_wat*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[0] += 1;
                                                 if(matrix[x].theta_taj*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_taj*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[1] += 1;
                                                if(matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_fulin*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[2] += 1;
                                                if(matrix[x].theta_wat*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_wat*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[0] += 1;
                                                 if(matrix[x].theta_taj*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_taj*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[1] += 1;
                                                if(matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_fulin*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[2] += 1;
                                            }
                                            ibonf[0][x] = (double)nunderalfa[0]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                            if(ibonf[0][x] > (double)0.5) ibonf[0][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                            else ibonf[0][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                            
											ibonf[1][x] = (double)nunderalfa[1]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                            if(ibonf[1][x] > (double)0.5) ibonf[1][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                            else ibonf[1][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                            
											ibonf[2][x] = (double)nunderalfa[2]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                            if(ibonf[2][x] > (double)0.5) ibonf[2][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                            else ibonf[2][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                        }
										ibonf_sort(ibonf,ibonfn,3,data[0].n_loci);
										
										for(x=0;x<data[0].n_loci;x++) {
                                            printf("%d:%s\t",x,matrix[x].gene);
                                            if(file_output) fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
                                            
											printf("%g\t%g\t",matrix[x].theta_wat*((double)1/matrix[x].factor_chrn),ibonf[0][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_wat*((double)1/matrix[x].factor_chrn),ibonf[0][x]);
											print_ibonf(ibonf[0][x],ibonfn[0][x],data[0].n_loci,data[0].n_iter,file_output);
                                            
											printf("%g\t%g\t",matrix[x].theta_taj*((double)1/matrix[x].factor_chrn),ibonf[1][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_taj*((double)1/matrix[x].factor_chrn),ibonf[1][x]);
											print_ibonf(ibonf[1][x],ibonfn[1][x],data[0].n_loci,data[0].n_iter,file_output);
                                            
											printf("%g\t%g\t",matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn),ibonf[2][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn),ibonf[2][x]);
											print_ibonf(ibonf[2][x],ibonfn[2][x],data[0].n_loci,data[0].n_iter,file_output);
											
											printf("\n");
											if(file_output) fprintf(file_output,"\n");
										}
                                    }
                                }
                                else {
                                    printf("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n");
                                    if(file_output)
                                        fputs("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n",file_output);
                                }
                            }
							#if ALLRESULTS == 0
                            else {
							#endif
                                if(!onlymulo) {
                                    /*calculate median, 10% 5% etc./locus*/
                                    if((sortvector = (double *) malloc(data[0].n_iter * sizeof(double))) == 0) {
                                        printf("Probabilities can not be calculated, sorry.\n");
                                        if(file_output) fputs("Probabilities can not be calculated, sorry.\n",file_output);
                                        break;
                                    }
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n");
                                    if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									if(file_output)
                                        fputs("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n",file_output);
                                    if(*outgroup) {    
                                        for(x=0;x<data[0].n_loci;x++) {
                                            printf("\nloci\tstatistic\tmedian\tavg\tvar\t");                                    
                                            if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                            if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                            if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                            if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                            if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                            if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                            if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                            printf("\n");
                                            if(file_output) {
                                                fputs("\nloci\tstatistic\tmedian\tavg\tvar\t",file_output);                                    
                                                if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                                if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                                if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                                if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                                if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                                if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                                if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                                fprintf(file_output,"\n");
                                            }

                                            /*theta_wat */
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_wat != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_wat*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_wat*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_wat * (double)matrixsim[x][it].theta_wat *
														((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(wat)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output) 
                                                fprintf(file_output,"%d\ttheta(wat)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*theta_taj*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_taj != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_taj*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_taj*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_taj * (double)matrixsim[x][it].theta_taj
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(taj)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\ttheta(taj)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*theta_fuli*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_fuli != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_fuli*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_fuli*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_fuli * (double)matrixsim[x][it].theta_fuli
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(fuli)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\ttheta(fuli)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*theta_faywu*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_fw != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_fw*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_fw*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_fw * (double)matrixsim[x][it].theta_fw
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(faywu)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\ttheta(faywu)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*theta_zeng*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_L != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_L*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_L*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_L * (double)matrixsim[x][it].theta_L
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(zeng)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\ttheta(zeng)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                        }
                                    }
                                    else {
                                        for(x=0;x<data[0].n_loci;x++) {
                                            printf("\nloci\tstatistic\tmedian\tavg\tvar\t");                                    
                                            if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                            if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                            if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                            if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                            if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                            if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                            if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                            printf("\n");
                                            if(file_output) {
                                                fputs("\nloci\tstatistic\tmedian\tavg\tvar\t",file_output);                                    
                                                if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                                if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                                if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                                if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                                if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                                if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                                if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                                fprintf(file_output,"\n");
                                            }
											/*theta_wat */
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_wat != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_wat*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_wat*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_wat * (double)matrixsim[x][it].theta_wat
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(wat)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output) 
                                                fprintf(file_output,"%d\ttheta(wat)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*theta_taj*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_taj != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_taj*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_taj*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_taj * (double)matrixsim[x][it].theta_taj
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(taj)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\ttheta(taj)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*theta_fulin*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_fulin != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_fulin*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_fulin*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_fulin * (double)matrixsim[x][it].theta_fulin
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(fulin)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\ttheta(fulin)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                        }
                                    }
                                    free(sortvector);
                                }
                                else {
                                    printf("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n");
                                    if(file_output)
                                        fputs("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n",file_output);
                                }
							#if ALLRESULTS == 0
                            }
							#endif
                            break;
                        case '3': /* 3 - Display Theta/nt.*/
                            if(*observed_data && dataobsequalsim) {/*compare simulations with observed data.*/
                                /*Sign tests.*/
                                /*Binomial for each statistic. from median (all data) or average (onlymulo). outgroup and non outgroup.*/
                                if((psigned = (int *)calloc(5,sizeof(int))) == 0) {
                                    printf("Error in memory allocation: Sign tests are not available.");
                                    if(file_output) fputs("Error in memory allocation: Sign tests are not available.",file_output);
                                    break;
                                }
                                if((nsigned = (int *)calloc(5,sizeof(int))) == 0) {
                                    printf("Error in memory allocation: Sign tests are not available.");
                                    if(file_output) fputs("Error in memory allocation: Sign tests are not available.",file_output);
                                    break;
                                }

                                if(onlymulo) {
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									printf(" NOTE that all theta estimates are corrected for chromosome population size.\n");
                                    printf("\nNon-parametric SIGN TEST. Assuming that average have equivalent value than median. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl)");
                                    printf("\nComparing observed data with the average of simulated data at each locus.");
                                    printf("\nAdjust to a Binomial.\n");
                                    printf("The table contains in order the name of the locus");
                                    printf(" and nucleotide estimations of variability from Watterson, Tajima, Fu and Li,\n");
                                    printf("Fay and Wu and Zeng et al.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
										fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n",file_output);
                                        fputs("\nNon-parametric SIGN TEST. Assuming that average have equivalent value than median. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl)",file_output);
                                        fputs("\nComparing observed data with the average of simulated data at each locus.",file_output);
                                        fputs("\nAdjust to a Binomial.\n",file_output);
                                        fprintf(file_output,"The table contains in order the name of the locus");
                                        fprintf(file_output," and nucleotide estimations of variability from Watterson, Tajima, Fu and Li,\n");
                                        fprintf(file_output,"Fay and Wu and Zeng et al.\n");
                                    }
                                }
                                else {
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									printf(" NOTE that all theta estimates are corrected for chromosome population size.\n");
                                    printf("\nNon-parametric SIGN TEST. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl).");
                                    printf("\nComparing observed data with the median of simulated data at each locus.");
                                    printf("\nAdjust to a Binomial.\n");
                                    printf("The table contains in order the name of the locus");
                                    printf(" and nucleotide estimations of variability from Watterson, Tajima, Fu and Li,\n");
                                    printf("Fay and Wu and Zeng et al.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
										fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n",file_output);
                                        fputs("\nNon-parametric SIGN TEST. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl).",file_output);
                                        fputs("\nComparing observed data with the median of simulated data at each locus.",file_output);
                                        fputs("\nAdjust to a Binomial.\n",file_output);
                                        fprintf(file_output,"The table contains in order the name of the locus");
                                        fprintf(file_output," and nucleotide estimations of variability from Watterson, Tajima, Fu and Li,\n");
                                        fprintf(file_output,"Fay and Wu and Zeng et al.\n");
                                    }
                                }
                                if(*outgroup) {
                                    printf("\nloci\tObs(theta_wat)\tExp(theta_wat)\tObs(theta_taj)\tExp(theta_taj)\tObs(theta_fuli)\tExp(theta_fuli)\tObs(theta_fw)\tExp(theta_fw)\tObs(theta_zeng)\tExp(theta_zeng)\n");
                                    if(file_output)
                                        fputs("\nloci\tObs(theta_wat)\tExp(theta_wat)\tObs(theta_taj)\tExp(theta_taj)\tObs(theta_fuli)\tExp(theta_fuli)\tObs(theta_fw)\tExp(theta_fw)\tObs(theta_zeng)\tExp(theta_zeng)\n",file_output);
                                    for(x=0;x<data[0].n_loci;x++) {
                                        printf("%d:%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",x,matrix[x].gene,
                                            matrix[x].theta_wat/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_wat_nut,
                                            matrix[x].theta_taj/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_taj_nut,
                                            matrix[x].theta_fuli/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_fuli_nut,
                                            matrix[x].theta_fw/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_fw_nut,
											matrix[x].theta_L/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_L_nut);
                                        if(file_output)
                                            fprintf(file_output,"%d:%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",x,matrix[x].gene,
                                                matrix[x].theta_wat/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_wat_nut,
                                                matrix[x].theta_taj/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_taj_nut,
                                                matrix[x].theta_fuli/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_fuli_nut,
                                                matrix[x].theta_fw/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_fw_nut,
												matrix[x].theta_L/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_L_nut);

                                        if(matrix[x].theta_wat/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_wat_nut) psigned[0] += 1;
                                        if(matrix[x].theta_taj/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_taj_nut) psigned[1] += 1;
                                        if(matrix[x].theta_fuli/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_fuli_nut) psigned[2] += 1;
                                        if(matrix[x].theta_fw/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_fw_nut) psigned[3] += 1;
                                        if(matrix[x].theta_L/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_L_nut) psigned[4] += 1;

                                        if(matrix[x].theta_wat/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_wat_nut) nsigned[0] += 1;
                                        if(matrix[x].theta_taj/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_taj_nut) nsigned[1] += 1;
                                        if(matrix[x].theta_fuli/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_fuli_nut) nsigned[2] += 1;
                                        if(matrix[x].theta_fw/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_fw_nut) nsigned[3] += 1;
                                        if(matrix[x].theta_L/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_L_nut) nsigned[4] += 1;
                                    }
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    for(x=0;x<5;x++) {
                                        switch(x) {
                                            case 0:
                                                printf("\nSign test: Theta (Watterson) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Watterson) corrected for chromosome population size.\n",file_output);
                                                break;
                                            case 1:
                                                printf("\nSign test: Theta (Tajima) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Tajima) corrected for chromosome population size.\n",file_output);
                                                break;
                                            case 2:
                                                printf("\nSign test: Theta (Fu and Li) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Fu and Li) corrected for chromosome population size.\n",file_output);
                                                break;
                                            case 3:
                                                printf("\nSign test: Theta (Fay and Wu) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Fay and Wu) corrected for chromosome population size.\n",file_output);
                                                break;
                                            case 4:
                                                printf("\nSign test: Theta (zeng et al.) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (zeng et al.) corrected for chromosome population size.\n",file_output);
                                                break;
                                        }
                                        prob = prob_cumbinomial(psigned[x]+nsigned[x],psigned[x],0.5);
										if(prob == -10000) printf("#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
                                        else printf("#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
										if(file_output) {
											if(prob == -10000) fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
											else fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
										}
                                        printfignif1(prob,file_output);
										
										/*look at the effect of discrete values*/
										if(prob <= 0.05 && prob >= 0.) {
											eq = data[0].n_loci - (psigned[x]+nsigned[x]);
											if(psigned[x] > nsigned[x]) {
												ns = nsigned[x] + eq;
												ps = psigned[x];
											}
											else {
												ns = nsigned[x];
												ps = psigned[x] + eq;
											}
											if((psigned[x] > nsigned[x] && ps < ns) || (psigned[x] < nsigned[x] && ps > ns) || prob_cumbinomial(ps+ns,ps,0.5) > 0.05) {
                                                printf(" WARNING: too many values equal to the median. It might be too liberal.");
                                                if(file_output) fputs(" WARNING: too many values equal to the median. It might be too liberal.",file_output);
											}
										}
                                    }
                                }
                                else {
									puts("\nloci\tObs(theta_wat)\tExp(theta_wat)\tObs(theta_taj)\tExp(theta_taj)\tObs(theta_fuli)\tExp(theta_fuli)\n");
                                    if(file_output)
                                        fputs("\nloci\tObs(theta_wat)\tExp(theta_wat)\tObs(theta_taj)\tExp(theta_taj)\tObs(theta_fuli)\tExp(theta_fuli)\n",file_output);
                                    for(x=0;x<data[0].n_loci;x++) {
                                        printf("%d:%s\t%g\t%g\t%g\t%g\t%g\t%g\n",x,matrix[x].gene,
                                            matrix[x].theta_wat/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_wat_nut,
                                            matrix[x].theta_taj/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_taj_nut,
                                            matrix[x].theta_fulin/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_fulin_nut);
                                        if(file_output)
                                            fprintf(file_output,"%d:%s\t%g\t%g\t%g\t%g\t%g\t%g\n",x,matrix[x].gene,
                                                matrix[x].theta_wat/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_wat_nut,
                                                matrix[x].theta_taj/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_taj_nut,
                                                matrix[x].theta_fulin/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),avgstatloci[x].theta_fulin_nut);

                                        if(matrix[x].theta_wat/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_wat_nut) psigned[0] += 1;
                                        if(matrix[x].theta_taj/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_taj_nut) psigned[1] += 1;
                                        if(matrix[x].theta_fulin/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) > avgstatloci[x].theta_fulin_nut) psigned[2] += 1;

                                        if(matrix[x].theta_wat/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_wat_nut) nsigned[0] += 1;
                                        if(matrix[x].theta_taj/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_taj_nut) nsigned[1] += 1;
                                        if(matrix[x].theta_fulin/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) < avgstatloci[x].theta_fulin_nut) nsigned[2] += 1;
                                    }
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    for(x=0;x<3;x++) {
                                        switch(x) {
                                            case 0:
                                                printf("\nSign test: Theta (Watterson) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Watterson) corrected for chromosome population size.\n",file_output);
                                                break;
                                            case 1:
                                                printf("\nSign test: Theta (Tajima) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Tajima) corrected for chromosome population size.\n",file_output);
                                                break;
                                            case 2:
                                                printf("\nSign test: Theta (Fu and Li) corrected for chromosome population size.\n");
                                                if(file_output) fputs("\nSign test: Theta (Fu and Li) corrected for chromosome population size.\n",file_output);
                                                break;
                                        }
                                        prob = prob_cumbinomial(psigned[x]+nsigned[x],psigned[x],0.5);
										if(prob == -10000) printf("#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
                                        else printf("#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
										if(file_output) {
											if(prob == -10000) fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
											else fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
										}
                                        printfignif1(prob,file_output);

										/*look at the effect of discrete values*/
										if(prob <= 0.05 && prob >= 0.) {
											eq = data[0].n_loci - (psigned[x]+nsigned[x]);
											if(psigned[x] > nsigned[x]) {
												ns = nsigned[x] + eq;
												ps = psigned[x];
											}
											else {
												ns = nsigned[x];
												ps = psigned[x] + eq;
											}
											if((psigned[x] > nsigned[x] && ps < ns) || (psigned[x] < nsigned[x] && ps > ns) || prob_cumbinomial(ps+ns,ps,0.5) > 0.05) {
                                                printf(" WARNING: too many values equal to the median. It might be too liberal.");
                                                if(file_output) fputs(" WARNING: too many values equal to the median. It might be too liberal.",file_output);
											}
										}
									}
                                }
                                free(psigned);
                                free(nsigned);

                                if(!onlymulo) {
                                    /*calculate probability each observed value/locus*/
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\n\n Display the probability of observed statistics for each locus obtained by coalescent Monte Carlo simulations: ");
                                    printf("\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. Note: & means the ierations would be increased to confirm significance");
                                    printf("\n Multiple hits not included in the analysis.\n");
                                    printf("The table contains in order the name of the locus");
                                    printf(" and nucleotide estimations of variability from Watterson, Tajima, Fu and Li,\n");
                                    printf("and normalized Fay and Wu (if outgroup).\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\n\n Display the probability of observed statistics for each locus obtained by coalescent Monte Carlo simulations: ",file_output);
                                        fprintf(file_output,"\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. Note: & means the ierations would be increased to confirm significance");
                                        fputs("\n Multiple hits not included in the analysis.\n",file_output);
                                        fprintf(file_output,"The table contains in order the name of the locus");
                                        fprintf(file_output," and locus estimations of variability from Watterson, Tajima, Fu and Li,\n");
                                        fprintf(file_output,"and normalized Fay and Wu (if outgroup).\n");
                                    }
                                    if(*outgroup) {
                                        printf("\nloci\tObs(theta_wat)\tP(theta_wat)\tIBonf(theta_wat)\tObs(theta_taj)\tP(theta_taj)\tIBonf(theta_taj)\tObs(theta_fuli)\tP(theta_fuli)\tIBonf(theta_fuli)\tObs(theta_fw)\tP(theta_fw)\tIBonf(theta_fw)\tP(theta_zeng)\tIBonf(theta_zeng)\n");
                                        if(file_output) 
                                            fputs("\nloci\tObs(theta_wat)\tP(theta_wat)\tIBonf(theta_wat)\tObs(theta_taj)\tP(theta_taj)\tIBonf(theta_taj)\tObs(theta_fuli)\tP(theta_fuli)\tIBonf(theta_fuli)\tObs(theta_fw)\tP(theta_fw)\tIBonf(theta_fw)\tP(theta_zeng)\tIBonf(theta_zeng)\n",file_output);
                                        for(x=0;x<data[0].n_loci;x++) {
                                            nunderalfa[0] = nunderalfa[1] = nunderalfa[2] = nunderalfa[3] = nunderalfa[4] = 0;
                                            nequivalfa[0] = nequivalfa[1] = nequivalfa[2] = nequivalfa[3] = nequivalfa[4] = 0;
                                            for(it=0;it<(long int)data[0].n_iter;it++) {
                                                if(matrix[x].theta_wat/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_wat_nut*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[0] += 1;
                                                 if(matrix[x].theta_taj/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_taj_nut*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[1] += 1;
                                                if(matrix[x].theta_fuli/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_fuli_nut*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[2] += 1;
                                                if(matrix[x].theta_fw/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_fw_nut*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[3] += 1;
                                                if(matrix[x].theta_L/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_L_nut*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[4] += 1;
                                                if(matrix[x].theta_wat/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_wat_nut*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[0] += 1;
                                                 if(matrix[x].theta_taj/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_taj_nut*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[1] += 1;
                                                if(matrix[x].theta_fuli/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_fuli_nut*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[2] += 1;
                                                if(matrix[x].theta_fw/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_fw_nut*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[3] += 1;
                                                if(matrix[x].theta_L/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_L_nut*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[4] += 1;
                                            }
                                            ibonf[0][x] = (double)nunderalfa[0]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                            if(ibonf[0][x] > (double)0.5) ibonf[0][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                            else ibonf[0][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                            
											ibonf[1][x] = (double)nunderalfa[1]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                            if(ibonf[1][x] > (double)0.5) ibonf[1][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                            else ibonf[1][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                            
											ibonf[2][x] = (double)nunderalfa[2]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                            if(ibonf[2][x] > (double)0.5) ibonf[2][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                            else ibonf[2][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                            
											ibonf[3][x] = (double)nunderalfa[3]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                            if(ibonf[3][x] > (double)0.5) ibonf[3][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                            else ibonf[3][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                            
											ibonf[4][x] = (double)nunderalfa[4]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                            if(ibonf[4][x] > (double)0.5) ibonf[4][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                            else ibonf[4][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                        }
										ibonf_sort(ibonf,ibonfn,5,data[0].n_loci);
										
										for(x=0;x<data[0].n_loci;x++) {
		                                    printf("%d:%s\t",x,matrix[x].gene);
                                            if(file_output) fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
                                            
											printf("%g\t%g\t",matrix[x].theta_wat/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn),ibonf[0][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_wat/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn),ibonf[0][x]);
											print_ibonf(ibonf[0][x],ibonfn[0][x],data[0].n_loci,data[0].n_iter,file_output);
                                            
											printf("%g\t%g\t",matrix[x].theta_taj/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn),ibonf[1][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_taj/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn),ibonf[1][x]);
											print_ibonf(ibonf[1][x],ibonfn[1][x],data[0].n_loci,data[0].n_iter,file_output);
                                            
											printf("%g\t%g\t",matrix[x].theta_fuli/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn),ibonf[2][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_fuli/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn),ibonf[2][x]);
											print_ibonf(ibonf[2][x],ibonfn[2][x],data[0].n_loci,data[0].n_iter,file_output);
                                            
											printf("%g\t%g\t",matrix[x].theta_fw/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn),ibonf[3][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_fw/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn),ibonf[3][x]);
											print_ibonf(ibonf[3][x],ibonfn[3][x],data[0].n_loci,data[0].n_iter,file_output);
                                            
											printf("%g\t%g\t",matrix[x].theta_L/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn),ibonf[4][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_L/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn),ibonf[4][x]);
											print_ibonf(ibonf[4][x],ibonfn[4][x],data[0].n_loci,data[0].n_iter,file_output);
											
											printf("\n");
											if(file_output) fprintf(file_output,"\n");
        								}
                                    }
                                    else {
                                        printf("\nloci\tObs(theta_wat)\tP(theta_wat)\tIBonf(theta_wat)\tObs(theta_taj)\tP(theta_taj)\tIBonf(theta_taj)\tObs(theta_fuli)\tP(theta_fuli)\tIBonf(theta_fuli)\n");
                                        if(file_output) 
                                            fputs("\nloci\tObs(theta_wat)\tP(theta_wat)\tIBonf(theta_wat)\tObs(theta_taj)\tP(theta_taj)\tIBonf(theta_taj)\tObs(theta_fuli)\tP(theta_fuli)\tIBonf(theta_fuli)\n",file_output);
                                        for(x=0;x<data[0].n_loci;x++) {
                                            nunderalfa[0] = nunderalfa[1] = nunderalfa[2] = 0;
                                            nequivalfa[0] = nequivalfa[1] = nequivalfa[2] = 0;
                                            for(it=0;it<(long int)data[0].n_iter;it++) {
                                                if(matrix[x].theta_wat/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_wat_nut*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[0] += 1;
                                                 if(matrix[x].theta_taj/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_taj_nut*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[1] += 1;
                                                if(matrix[x].theta_fulin/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) < matrixsim[x][it].theta_fulin_nut*((double)1/inputms[x].factor_chrn))
                                                    nunderalfa[2] += 1;
                                                if(matrix[x].theta_wat/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_wat_nut*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[0] += 1;
                                                 if(matrix[x].theta_taj/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_taj_nut*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[1] += 1;
                                                if(matrix[x].theta_fulin/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) == matrixsim[x][it].theta_fulin_nut*((double)1/inputms[x].factor_chrn))
                                                    nequivalfa[2] += 1;
                                            }
                                           
                                            ibonf[0][x] = (double)nunderalfa[0]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                            if(ibonf[0][x] > (double)0.5) ibonf[0][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                            else ibonf[0][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                            
                                            ibonf[1][x] = (double)nunderalfa[1]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                            if(ibonf[1][x] > (double)0.5) ibonf[1][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                            else ibonf[1][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                            
                                            ibonf[2][x] = (double)nunderalfa[2]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                            if(ibonf[2][x] > (double)0.5) ibonf[2][x] -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                            else ibonf[2][x] += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                        }
										ibonf_sort(ibonf,ibonfn,3,data[0].n_loci);
										
										for(x=0;x<data[0].n_loci;x++) {
											printf("%d:%s\t",x,matrix[x].gene);
                                            if(file_output) fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
                                            
											printf("%g\t%g\t",matrix[x].theta_wat/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),ibonf[0][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_wat/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),ibonf[0][x]);
                                            print_ibonf(ibonf[0][x],ibonfn[0][x],data[0].n_loci,data[0].n_iter,file_output);
											
											printf("%g\t%g\t",matrix[x].theta_taj/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),ibonf[1][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_taj/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),ibonf[1][x]);
                                            print_ibonf(ibonf[1][x],ibonfn[1][x],data[0].n_loci,data[0].n_iter,file_output);
											
											printf("%g\t%g\t",matrix[x].theta_fulin/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),ibonf[2][x]);
                                            if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].theta_fulin/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn),ibonf[2][x]);
											print_ibonf(ibonf[2][x],ibonfn[2][x],data[0].n_loci,data[0].n_iter,file_output);
											
											printf("\n");
											if(file_output) fprintf(file_output,"\n");
										}
                                    }
                                }
                                else {
                                    printf("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n");
                                    if(file_output)
                                        fputs("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n",file_output);
                                }
                            }
							#if ALLRESULTS == 0
                            else {
							#endif
                                if(!onlymulo) {
                                    /*calculate median, 10% 5% etc./locus*/
                                    if((sortvector = (double *) malloc(data[0].n_iter * sizeof(double))) == 0) {
                                        printf("Probabilities can not be calculated, sorry.\n");
                                        if(file_output) fputs("Probabilities can not be calculated, sorry.\n",file_output);
                                        break;
                                    }
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n");
                                    if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
									if(file_output)
                                        fputs("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n",file_output);
                                    if(*outgroup) {    
                                        for(x=0;x<data[0].n_loci;x++) {
                                            printf("\nloci\tstatistic\tmedian\tavg\tvar\t");                                    
                                            if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                            if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                            if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                            if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                            if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                            if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                            if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                            printf("\n");
                                            if(file_output) {
                                                fputs("\nloci\tstatistic\tmedian\tavg\tvar\t",file_output);                                    
                                                if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                                if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                                if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                                if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                                if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                                if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                                if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                                fprintf(file_output,"\n");
                                            }

                                            /*theta_wat */
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_wat != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_wat_nut*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_wat_nut*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_wat_nut * (double)matrixsim[x][it].theta_wat_nut
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(wat)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output) 
                                                fprintf(file_output,"%d\ttheta(wat)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*theta_taj*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_taj != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_taj_nut*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_taj_nut *((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_taj_nut * (double)matrixsim[x][it].theta_taj_nut
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(taj)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\ttheta(taj)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*theta_fuli*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_fuli != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_fuli_nut*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_fuli_nut*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_fuli_nut * (double)matrixsim[x][it].theta_fuli_nut
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(fuli)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\ttheta(fuli)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*theta_faywu*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_fw != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_fw_nut*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_fw_nut*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_fw_nut * (double)matrixsim[x][it].theta_fw_nut
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(faywu)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\ttheta(faywu)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*theta_zeng*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_L != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_L_nut*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_L_nut*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_L_nut * (double)matrixsim[x][it].theta_L_nut
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(zeng)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\ttheta(zeng)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");
                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                        }
                                    }
                                    else {
                                        for(x=0;x<data[0].n_loci;x++) {
                                            printf("\nloci\tstatistic\tmedian\tavg\tvar\t");                                    
                                            if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                            if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                            if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                            if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                            if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                            if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                            if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                            printf("\n");
                                            if(file_output) {
                                                fputs("\nloci\tstatistic\tmedian\tavg\tvar\t",file_output);                                    
                                                if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                                if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                                if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                                if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                                if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                                if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                                if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                                fprintf(file_output,"\n");
                                            }
                                             /*theta_wat */
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_wat != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_wat_nut*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_wat_nut*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_wat_nut * (double)matrixsim[x][it].theta_wat_nut
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(wat)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output) 
                                                fprintf(file_output,"%d\ttheta(wat)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");
                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*theta_taj*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_taj != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_taj_nut*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_taj_nut*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_taj_nut * (double)matrixsim[x][it].theta_taj_nut
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(taj)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\ttheta(taj)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");
                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                            
											/*theta_fuli*/
                                            /*average*/
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if( matrixsim[x][it].theta_fulin != (double) -10000) {
                                                    sortvector[it2] = (double)matrixsim[x][it].theta_fulin_nut*((double)1/inputms[x].factor_chrn);
                                                    it2 += 1;
                                                    avg += (double)matrixsim[x][it].theta_fulin_nut*((double)1/inputms[x].factor_chrn);
                                                    var += (double)matrixsim[x][it].theta_fulin_nut * (double)matrixsim[x][it].theta_fulin_nut
														*((double)1/inputms[x].factor_chrn)*((double)1/inputms[x].factor_chrn);
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("%d\ttheta(fuli)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"%d\ttheta(fuli)\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");
                                            }
                                            print_percentages(data,sortvector,it2,file_output);
                                        }
                                    }
                                    free(sortvector);
                                }
                                else {
                                    printf("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n");
                                    if(file_output)
                                        fputs("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n",file_output);
                                }
							#if ALLRESULTS == 0
                            }
							#endif
                            break;
                        case '4': /* 4 - Back to the previous menu.*/
                            return;
                            break;
                    }
					for(x=0;x<10;x++) {
						free(ibonf[x]);
						free(ibonfn[x]);
					}
					free(ibonf);
					free(ibonfn);
                    break;
                case 1: /*neutrality tests.*/
					if(file_output) {    
						fprintf(file_output,"\n\n     MENU:\n     3. Statistical inference based on Coalescent Monte Carlo simulations menu:");
						fprintf(file_output,"\n     3.1. Display coalescent simulation analysis menu:");
						fprintf(file_output,"\n     3.1.1. Display neutrality tests menu:");
						fprintf(file_output,"\n     3.1.1.1. Display detailed neutrality tests menu:\n\n");
						fprintf(file_output,"\n     3.1.1.1.1. Display the probability values (from observed data) or confident intervals.\n");
					}

					if((ibonf = (double **)calloc(13,sizeof(double *))) == 0) {
						printf("Error in memory allocation: IBt.");
						if(file_output) 
							fputs("Error in memory allocation: IBt.",file_output);
						break;
					}
					for(x=0;x<13;x++) {
						if((ibonf[x] = (double *)calloc(data[0].n_loci,sizeof(double))) == 0) {
							printf("Error in memory allocation: IBt.");
							if(file_output) 
								fputs("Error in memory allocation: IBt.",file_output);
							break;
						}
					}
					if((ibonfn = (int **)calloc(13,sizeof(int *))) == 0) {
						printf("Error in memory allocation: IBt.");
						if(file_output) 
							fputs("Error in memory allocation: IBt.",file_output);
						break;
					}
					for(x=0;x<13;x++) {
						if((ibonfn[x] = (int *)calloc(data[0].n_loci,sizeof(int))) == 0) {
							printf("Error in memory allocation: IBt.");
							if(file_output) 
								fputs("Error in memory allocation: IBt.",file_output);
							break;
						}
					}
                    if(*observed_data && dataobsequalsim) {/*compare simulations with observed data.*/
                        /*Sign tests.*/
                        /*Binomial for each statistic. from median (all data) or average (onlymulo). outgroup and non outgroup.*/
                        if((psigned = (int *)calloc(13,sizeof(int))) == 0) {
                            printf("Error in memory allocation: Sign tests are not available.");
                            if(file_output) fputs("Error in memory allocation: Sign tests are not available.",file_output);
                            break;
                        }
                        if((nsigned = (int *)calloc(13,sizeof(int))) == 0) {
                            printf("Error in memory allocation: Sign tests are not available.");
                            if(file_output) fputs("Error in memory allocation: Sign tests are not available.",file_output);
                            break;
                        }
                        if((qsigned = (int *)calloc(13,sizeof(int))) == 0) {
                            printf("Error in memory allocation: Sign tests are not available.");
                            if(file_output) fputs("Error in memory allocation: Sign tests are not available.",file_output);
                            break;
                        }
						
                        if(onlymulo) {
							printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                            printf("\nNon-parametric SIGN TEST. Assuming that average have equivalent value than median. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl)");
                            printf("\nComparing observed data with the average of simulated data at each locus.");
                            printf("\nAdjust to a Binomial.");
                            printf("\nThe table contains in order the Observed and expected values of tests of Tajima's D, Fu and Li's D and F");
                            printf("\n(with '*' not outgroup), Fu's Fs, normalized Fay and Wu's H (in case outgroup), Fay and Wu's H (in case outgroup), Rozas's et al. ZA, Wall's B and Q");
                            printf("\nRamos-Onsins and Rozas' R2, Zeng et al. E test and Ewens-Watterson test.\n");
                            if(file_output) {
								if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                fputs("\nNon-parametric SIGN TEST. Assuming that average have equivalent value than median. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl)",file_output);
                                fputs("\nComparing observed data with the average of simulated data at each locus.",file_output);
                                fputs("\nAdjust to a Binomial.",file_output);
                                fprintf(file_output,"\nThe table contains in order the Observed and expected values of tests of Tajima's D, Fu and Li's D and F");
                                fprintf(file_output,"\n(with '*' not outgroup), Fu's Fs, normalized Fay and Wu's H (in case outgroup), Fay and Wu's H (in case outgroup), Rozas's et al. ZA, Wall's B and Q");
                                fprintf(file_output,"\nRamos-Onsins and Rozas' R2, Zeng et al. E test and Ewens-Watterson test.\n");
                            }
                        }
                        else {
							printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                            printf("\nNon-parametric SIGN TEST. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl).");
                            printf("\nComparing observed data with the median of simulated data at each locus.");
                            printf("\nAdjust to a Binomial.");
                            printf("\nThe table contains in order the Observed and expected values of tests of Tajima's D, Fu and Li's D and F");
                            printf("\n(with '*' not outgroup), Fu's Fs, normalized Fay and Wu's H (in case outgroup), Fay and Wu's H (in case outgroup), Rozas's et al. ZA, Wall's B and Q");
                            printf("\nRamos-Onsins and Rozas' R2, Zeng et al. E test and Ewens-Watterson test.\n");
                            if(file_output) {
								if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                fputs("\nNon-parametric SIGN TEST. \nWARNING: Sign test might be falsely significant using discrete statistics (i.e., nhapl).",file_output);
                                fputs("\nComparing observed data with the median of simulated data at each locus.",file_output);
                                fputs("\nAdjust to a Binomial.",file_output);
                                fprintf(file_output,"\nThe table contains in order the Observed and expected values of tests of Tajima's D, Fu and Li's D and F");
                                fprintf(file_output,"\n(with '*' not outgroup), Fu's Fs, normalized Fay and Wu's H (in case outgroup), Fay and Wu's H (in case outgroup), Rozas's et al. ZA, Wall's B and Q");
                                fprintf(file_output,"\nRamos-Onsins and Rozas' R2, Zeng et al. E test and Ewens-Watterson test.\n");
                            }
                        }
                        if(*outgroup) {
                            printf("\n\nloci\tObs(Tajima)\tExp(Tajima)\tObs(Fu&LiD)\tExp(Fu&LiD)\tObs(Fu&LiF)\tExp(Fu&LiF)\tObs(Fs)\tExp(Fs)\tObs(FayWuHn)\tExp(FayWuHn)\tObs(FayWuH)\tExp(FayWuH)\tObs(ZA)\tExp(ZA)\tObs(B)\tExp(B)\tObs(Q)\tExp(Q)\tObs(R2)\tExp(R2)\tObs(ZengE)\tExp(ZengE)\tObs(EW)\tExp(EW)\n");
                            if(file_output)
                                fputs("\n\nloci\tObs(Tajima)\tExp(Tajima)\tObs(Fu&LiD)\tExp(Fu&LiD)\tObs(Fu&LiF)\tExp(Fu&LiF)\tObs(Fs)\tExp(Fs)\tObs(FayWuHn)\tExp(FayWuHn)\tObs(FayWuH)\tExp(FayWuH)\tObs(ZA)\tExp(ZA)\tObs(B)\tExp(B)\tObs(Q)\tExp(Q)\tObs(R2)\tExp(R2)\tObs(ZengE)\tExp(ZengE)\tObs(EW)\tExp(EW)\n",file_output);
                            for(x=0;x<data[0].n_loci;x++) {
                                printf("%d:%s\t",x,matrix[x].gene);
                                if(file_output) fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
                                if(matrix[x].tajimaD != (double) -10000) {
                                    printf("%g\t",matrix[x].tajimaD);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].tajimaD);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].tajimaD != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].tajimaD);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].tajimaD);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].fuliD != (double) -10000) {
                                    printf("%g\t",matrix[x].fuliD);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].fuliD);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].fuliD != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].fuliD);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].fuliD);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].fuliF != (double) -10000) {
                                    printf("%g\t",matrix[x].fuliF);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].fuliF);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].fuliF != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].fuliF);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].fuliF);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].fuFs != (double) -10000) {
                                    printf("%g\t",matrix[x].fuFs);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].fuFs);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].fuFs != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].fuFs);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].fuFs);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].faywuH != (double) -10000) {
                                    printf("%g\t",matrix[x].faywuH);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].faywuH);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].faywuH != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].faywuH);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].faywuH);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].faywuHo != (double) -10000) {
                                    printf("%g\t",matrix[x].faywuHo);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].faywuHo);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].faywuHo != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].faywuHo);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].faywuHo);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].rZA != (double) -10000) {
                                    printf("%g\t",matrix[x].rZA);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].rZA);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].rZA != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].rZA);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].rZA);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].wB != (double) -10000) {
                                    printf("%g\t",matrix[x].wB);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].wB);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].wB != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].wB);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].wB);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].wQ != (double) -10000) {
                                    printf("%g\t",matrix[x].wQ);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].wQ);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].wQ != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].wQ);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].wQ);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].R2 != (double) -10000) {
                                    printf("%g\t",matrix[x].R2);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].R2);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].R2 != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].R2);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].R2);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].zengE != (double) -10000) {
                                    printf("%g\t",matrix[x].zengE);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].zengE);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].zengE != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].zengE);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].zengE);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].ewtest != (double) -10000) {
                                    printf("%g\t",matrix[x].ewtest);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].ewtest);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].ewtest != (double) -10000) {
                                    printf("%g\n",avgstatloci[x].ewtest);
                                    if(file_output) fprintf(file_output,"%g\n",avgstatloci[x].ewtest);
                                }
                                else {
                                    printf("na\n");
                                    if(file_output) fprintf(file_output,"na\n");
                                }
                                if(matrix[x].tajimaD != (double)-10000 && avgstatloci[x].tajimaD != (double)-10000) {
                                    if(matrix[x].tajimaD > avgstatloci[x].tajimaD) psigned[0] += 1;
                                    if(matrix[x].tajimaD < avgstatloci[x].tajimaD) nsigned[0] += 1;
                                    if(matrix[x].tajimaD == avgstatloci[x].tajimaD) qsigned[0] += 1;
                                }
                                if(matrix[x].fuliD != (double)-10000 && avgstatloci[x].fuliD != (double)-10000) {
                                    if(matrix[x].fuliD > avgstatloci[x].fuliD) psigned[1] += 1;
                                    if(matrix[x].fuliD < avgstatloci[x].fuliD) nsigned[1] += 1;
                                    if(matrix[x].fuliD == avgstatloci[x].fuliD) qsigned[1] += 1;
                                }
                                if(matrix[x].fuliF != (double)-10000 && avgstatloci[x].fuliF != (double)-10000) {
                                    if(matrix[x].fuliF > avgstatloci[x].fuliF) psigned[2] += 1;
                                    if(matrix[x].fuliF < avgstatloci[x].fuliF) nsigned[2] += 1;
                                    if(matrix[x].fuliF == avgstatloci[x].fuliF) qsigned[2] += 1;
                                }
                                if(matrix[x].fuFs != (double)-10000 && avgstatloci[x].fuFs != (double)-10000) {
                                    if(matrix[x].fuFs > avgstatloci[x].fuFs) psigned[3] += 1;
                                    if(matrix[x].fuFs < avgstatloci[x].fuFs) nsigned[3] += 1;
                                    if(matrix[x].fuFs == avgstatloci[x].fuFs) qsigned[3] += 1;
                                }
                                if(matrix[x].faywuH != (double)-10000 && avgstatloci[x].faywuH != (double)-10000) {
                                    if(matrix[x].faywuH > avgstatloci[x].faywuH) psigned[4] += 1;
                                    if(matrix[x].faywuH < avgstatloci[x].faywuH) nsigned[4] += 1;
									if(matrix[x].faywuH == avgstatloci[x].faywuH) qsigned[4] += 1;
                                }
                                if(matrix[x].faywuHo != (double)-10000 && avgstatloci[x].faywuHo != (double)-10000) {
                                    if(matrix[x].faywuHo > avgstatloci[x].faywuHo) psigned[5] += 1;
                                    if(matrix[x].faywuHo < avgstatloci[x].faywuHo) nsigned[5] += 1;
                                    if(matrix[x].faywuHo == avgstatloci[x].faywuHo) qsigned[5] += 1;
                                }
                                if(matrix[x].rZA != (double)-10000 && avgstatloci[x].rZA != (double)-10000) {
                                    if(matrix[x].rZA > avgstatloci[x].rZA) psigned[6] += 1;
                                    if(matrix[x].rZA < avgstatloci[x].rZA) nsigned[6] += 1;
                                    if(matrix[x].rZA == avgstatloci[x].rZA) qsigned[6] += 1;
                                }
                                if(matrix[x].wB != (double)-10000 && avgstatloci[x].wB != (double)-10000) {
                                    if(matrix[x].wB > avgstatloci[x].wB) psigned[7] += 1;
                                    if(matrix[x].wB < avgstatloci[x].wB) nsigned[7] += 1;
                                    if(matrix[x].wB == avgstatloci[x].wB) qsigned[7] += 1;
                                }
                                if(matrix[x].wQ != (double)-10000 && avgstatloci[x].wQ != (double)-10000) {
                                    if(matrix[x].wQ > avgstatloci[x].wQ) psigned[8] += 1;
                                    if(matrix[x].wQ < avgstatloci[x].wQ) nsigned[8] += 1;
                                    if(matrix[x].wQ == avgstatloci[x].wQ) qsigned[8] += 1;
                                }
                                if(matrix[x].R2 != (double)-10000 && avgstatloci[x].R2 != (double)-10000) {
                                    if(matrix[x].R2 > avgstatloci[x].R2) psigned[9] += 1;
                                    if(matrix[x].R2 < avgstatloci[x].R2) nsigned[9] += 1;
                                    if(matrix[x].R2 == avgstatloci[x].R2) qsigned[9] += 1;
                                }
                                if(matrix[x].zengE != (double)-10000 && avgstatloci[x].zengE != (double)-10000) {
                                    if(matrix[x].zengE > avgstatloci[x].zengE) psigned[10] += 1;
                                    if(matrix[x].zengE < avgstatloci[x].zengE) nsigned[10] += 1;
                                    if(matrix[x].zengE == avgstatloci[x].zengE) qsigned[10] += 1;
                                }
                                if(matrix[x].ewtest != (double)-10000 && avgstatloci[x].ewtest != (double)-10000) {
                                    if(matrix[x].ewtest > avgstatloci[x].ewtest) psigned[11] += 1;
                                    if(matrix[x].ewtest < avgstatloci[x].ewtest) nsigned[11] += 1;
                                    if(matrix[x].ewtest == avgstatloci[x].ewtest) qsigned[11] += 1;
                                }
                            }
							printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
							if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                            for(x=0;x<12;x++) {
                                switch(x) {
                                    case 0:
                                        printf("\nSign test: Tajima's D.\n");
                                        if(file_output) fputs("\nSign test: Tajima's D.\n",file_output);
                                        break;
                                    case 1:
                                        printf("\nSign test: Fu and Li's D.\n");
                                        if(file_output) fputs("\nSign test: Fu and Li's D.\n",file_output);
                                        break;
                                    case 2:
                                        printf("\nSign test: Fu and Li's F.\n");
                                        if(file_output) fputs("\nSign test: Fu and Li's F.\n",file_output);
                                        break;
                                    case 3:
                                        printf("\nSign test: Fu's Fs.\n");
                                        if(file_output) fputs("\nSign test: Fu's Fs.\n",file_output);
                                        break;
                                    case 4:
                                        printf("\nSign test: normalized Fay and Wu's H.\n");
                                        if(file_output) fputs("\nSign test: normalized Fay and Wu's H.\n",file_output);
                                        break;
                                    case 5:
                                        printf("\nSign test: Fay and Wu's H.\n");
                                        if(file_output) fputs("\nSign test: Fay and Wu's H.\n",file_output);
                                        break;
                                    case 6:
                                        printf("\nSign test: Rozas' et al. ZA.\n");
                                        if(file_output) fputs("\nSign test: Rozas' et al. ZA.\n",file_output);
                                        break;
                                    case 7:
                                        printf("\nSign test: Wall's B.\n");
                                        if(file_output) fputs("\nSign test: Wall's B.\n",file_output);
                                        break;
                                    case 8:
                                        printf("\nSign test: Wall's Q.\n");
                                        if(file_output) fputs("\nSign test: Wall's Q.\n",file_output);
                                        break;
                                    case 9:
                                        printf("\nSign test: R2.\n");
                                        if(file_output) fputs("\nSign test: R2.\n",file_output);
                                        break;
                                    case 10:
                                        printf("\nSign test: Zeng E.\n");
                                        if(file_output) fputs("\nSign test: Zeng E.\n",file_output);
                                        break;
                                    case 11:
                                        printf("\nSign test: EW.\n");
                                        if(file_output) fputs("\nSign test: EW.\n",file_output);
                                        break;
                                }
                                if(psigned[x] == 0 && nsigned[x] == 0) {
                                    printf(" na ");
                                    if(file_output) fprintf(file_output,"na ");
                                }
                                else {
                                    prob = prob_cumbinomial(psigned[x]+nsigned[x],psigned[x],0.5);
									if(prob == -10000) printf("#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
                                    else printf("#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
									if(file_output) {
										if(prob == -10000) fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
										else fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
									}
                                    printfignif1(prob,file_output);

									/*look at the effect of discrete values*/
									if(prob <= 0.05 && prob >= 0.) {
										eq = qsigned[x];
										if(psigned[x] > nsigned[x]) {
											ns = nsigned[x] + eq;
											ps = psigned[x];
										}
										else {
											ns = nsigned[x];
											ps = psigned[x] + eq;
										}
										if((psigned[x] > nsigned[x] && ps < ns) || (psigned[x] < nsigned[x] && ps > ns) || prob_cumbinomial(ps+ns,ps,0.5) > 0.05) {
											printf(" WARNING: too many values equal to the median. It might be too liberal.");
											if(file_output) fputs(" WARNING: too many values equal to the median. It might be too liberal.",file_output);
										}
									}
								}
                            }
                        }
                        else {
							puts("\n\nloci\tObs(Tajima)\tExp(Tajima)\tObs(Fu&LiD*)\tExp(Fu&LiD*)\tObs(Fu&LiF*)\tExp(Fu&LiF*)\tObs(Fs)\tExp(Fs)\tObs(ZA)\tExp(ZA)\tObs(B)\tExp(B)\tObs(Q)\tExp(Q)\tObs(R2)\tExp(R2)\tObs(EW)\tExp(EW)\n");
                            if(file_output)
                                fputs("\n\nloci\tObs(Tajima)\tExp(Tajima)\tObs(Fu&LiD*)\tExp(Fu&LiD*)\tObs(Fu&LiF*)\tExp(Fu&LiF*)\tObs(Fs)\tExp(Fs)\tObs(ZA)\tExp(ZA)\tObs(B)\tExp(B)\tObs(Q)\tExp(Q)\tObs(R2)\tExp(R2)\tObs(EW)\tExp(EW)\n",file_output);
                            for(x=0;x<data[0].n_loci;x++) {
                                printf("%d:%s\t",x,matrix[x].gene);
                                if(file_output) fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
                                if(matrix[x].tajimaD != (double) -10000) {
                                    printf("%g\t",matrix[x].tajimaD);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].tajimaD);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].tajimaD != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].tajimaD);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].tajimaD);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].fuliDn != (double) -10000) {
                                    printf("%g\t",matrix[x].fuliDn);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].fuliDn);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].fuliDn != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].fuliDn);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].fuliDn);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].fuliFn != (double) -10000) {
                                    printf("%g\t",matrix[x].fuliFn);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].fuliFn);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].fuliFn != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].fuliFn);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].fuliFn);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].fuFs != (double) -10000) {
                                    printf("%g\t",matrix[x].fuFs);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].fuFs);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].fuFs != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].fuFs);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].fuFs);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].rZA != (double) -10000) {
                                    printf("%g\t",matrix[x].rZA);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].rZA);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].rZA != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].rZA);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].rZA);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].wB != (double) -10000) {
                                    printf("%g\t",matrix[x].wB);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].wB);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].wB != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].wB);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].wB);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].wQ != (double) -10000) {
                                    printf("%g\t",matrix[x].wQ);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].wQ);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].wQ != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].wQ);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].wQ);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].R2 != (double) -10000) {
                                    printf("%g\t",matrix[x].R2);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].R2);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].R2 != (double) -10000) {
                                    printf("%g\t",avgstatloci[x].R2);
                                    if(file_output) fprintf(file_output,"%g\t",avgstatloci[x].R2);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(matrix[x].ewtest != (double) -10000) {
                                    printf("%g\t",matrix[x].ewtest);
                                    if(file_output) fprintf(file_output,"%g\t",matrix[x].ewtest);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
                                }
                                if(avgstatloci[x].ewtest != (double) -10000) {
                                    printf("%g\n",avgstatloci[x].ewtest);
                                    if(file_output) fprintf(file_output,"%g\n",avgstatloci[x].ewtest);
                                }
                                else {
                                    printf("na\n");
                                    if(file_output) fprintf(file_output,"na\n");
                                }
                                if(matrix[x].tajimaD != (double)-10000 && avgstatloci[x].tajimaD != (double)-10000) {
                                    if(matrix[x].tajimaD > avgstatloci[x].tajimaD) psigned[0] += 1;
                                    if(matrix[x].tajimaD < avgstatloci[x].tajimaD) nsigned[0] += 1;
                                    if(matrix[x].tajimaD == avgstatloci[x].tajimaD) qsigned[0] += 1;
                                }
                                if(matrix[x].fuliDn != (double)-10000 && avgstatloci[x].fuliDn != (double)-10000) {
                                    if(matrix[x].fuliDn > avgstatloci[x].fuliDn) psigned[1] += 1;
                                    if(matrix[x].fuliDn < avgstatloci[x].fuliDn) nsigned[1] += 1;
                                    if(matrix[x].fuliDn == avgstatloci[x].fuliDn) qsigned[1] += 1;
                                }
                                if(matrix[x].fuliFn != (double)-10000 && avgstatloci[x].fuliFn != (double)-10000) {
                                    if(matrix[x].fuliFn > avgstatloci[x].fuliFn) psigned[2] += 1;
                                    if(matrix[x].fuliFn < avgstatloci[x].fuliFn) nsigned[2] += 1;
                                    if(matrix[x].fuliFn == avgstatloci[x].fuliFn) qsigned[2] += 1;
                                }
                                if(matrix[x].fuFs != (double)-10000 && avgstatloci[x].fuFs != (double)-10000) {
                                    if(matrix[x].fuFs > avgstatloci[x].fuFs) psigned[3] += 1;
                                    if(matrix[x].fuFs < avgstatloci[x].fuFs) nsigned[3] += 1;
                                    if(matrix[x].fuFs == avgstatloci[x].fuFs) qsigned[3] += 1;
                                }
                                if(matrix[x].rZA != (double)-10000 && avgstatloci[x].rZA != (double)-10000) {
                                    if(matrix[x].rZA > avgstatloci[x].rZA) psigned[4] += 1;
                                    if(matrix[x].rZA < avgstatloci[x].rZA) nsigned[4] += 1;
                                    if(matrix[x].rZA == avgstatloci[x].rZA) qsigned[4] += 1;
                                }
                                if(matrix[x].wB != (double)-10000 && avgstatloci[x].wB != (double)-10000) {
                                    if(matrix[x].wB > avgstatloci[x].wB) psigned[5] += 1;
                                    if(matrix[x].wB < avgstatloci[x].wB) nsigned[5] += 1;
                                    if(matrix[x].wB == avgstatloci[x].wB) qsigned[5] += 1;
                                }
                                if(matrix[x].wQ != (double)-10000 && avgstatloci[x].wQ != (double)-10000) {
                                    if(matrix[x].wQ > avgstatloci[x].wQ) psigned[6] += 1;
                                    if(matrix[x].wQ < avgstatloci[x].wQ) nsigned[6] += 1;
                                    if(matrix[x].wQ == avgstatloci[x].wQ) qsigned[6] += 1;
                                }
                                if(matrix[x].R2 != (double)-10000 && avgstatloci[x].R2 != (double)-10000) {
                                    if(matrix[x].R2 > avgstatloci[x].R2) psigned[7] += 1;
                                    if(matrix[x].R2 < avgstatloci[x].R2) nsigned[7] += 1;
                                    if(matrix[x].R2 == avgstatloci[x].R2) qsigned[7] += 1;
                                }
                                if(matrix[x].ewtest != (double)-10000 && avgstatloci[x].ewtest != (double)-10000) {
                                    if(matrix[x].ewtest > avgstatloci[x].ewtest) psigned[8] += 1;
                                    if(matrix[x].ewtest < avgstatloci[x].ewtest) nsigned[8] += 1;
                                    if(matrix[x].ewtest == avgstatloci[x].ewtest) qsigned[8] += 1;
                                }
                            }
							printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
							if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                            for(x=0;x<9;x++) {
                                switch(x) {
                                    case 0:
                                        printf("\nSign test: Tajima's D.\n");
                                        if(file_output) fputs("\nSign test: Tajima's D.\n",file_output);
                                        break;
                                    case 1:
                                        printf("\nSign test: Fu and Li's D*.\n");
                                        if(file_output) fputs("\nSign test: Fu and Li's D.\n",file_output);
                                        break;
                                    case 2:
                                        printf("\nSign test: Fu and Li's F*.\n");
                                        if(file_output) fputs("\nSign test: Fu and Li's F.\n",file_output);
                                        break;
                                    case 3:
                                        printf("\nSign test: Fu's Fs.\n");
                                        if(file_output) fputs("\nSign test: Fu's Fs.\n",file_output);
                                        break;
                                    case 4:
                                        printf("\nSign test: Rozas' et al. ZA.\n");
                                        if(file_output) fputs("\nSign test: Rozas' et al. ZA.\n",file_output);
                                        break;
                                    case 5:
                                        printf("\nSign test: Wall's B.\n");
                                        if(file_output) fputs("\nSign test: Wall's B.\n",file_output);
                                        break;
                                    case 6:
                                        printf("\nSign test: Wall's Q.\n");
                                        if(file_output) fputs("\nSign test: Wall's Q.\n",file_output);
                                        break;
                                    case 7:
                                        printf("\nSign test: R2.\n");
                                        if(file_output) fputs("\nSign test: R2.\n",file_output);
                                        break;
                                    case 8:
                                        printf("\nSign test: EW test.\n");
                                        if(file_output) fputs("\nSign test: EW test.\n",file_output);
                                        break;
                                }
                                if(psigned[x] == 0 && nsigned[x] == 0) {
                                    printf(" na ");
                                    if(file_output) fprintf(file_output," na ");
                                }
                                else {
                                    prob = prob_cumbinomial(psigned[x]+nsigned[x],psigned[x],0.5);
									if(prob == -10000) printf("#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
                                    else printf("#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
									if(file_output) {
										if(prob == -10000) fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = na ",psigned[x],nsigned[x]);
										else fprintf(file_output,"#loci above: %d\t #loci below %d\t  Prob = %g ",psigned[x],nsigned[x],prob);
									}
                                    printfignif1(prob,file_output);

									/*look at the effect of discrete values*/
									if(prob <= 0.05 && prob >= 0.) {
										eq = qsigned[x];
										if(psigned[x] > nsigned[x]) {
											ns = nsigned[x] + eq;
											ps = psigned[x];
										}
										else {
											ns = nsigned[x];
											ps = psigned[x] + eq;
										}
										if((psigned[x] > nsigned[x] && ps < ns) || (psigned[x] < nsigned[x] && ps > ns) || prob_cumbinomial(ps+ns,ps,0.5) > 0.05) {
											printf(" WARNING: too many values equal to the median. It might be too liberal.");
											if(file_output) fputs(" WARNING: too many values equal to the median. It might be too liberal.",file_output);
										}
									}
								}
                            }
                        }
                        free(psigned);
                        free(nsigned);
                        free(qsigned);

                        if(!onlymulo) {
                            /*calculate probability each observed value/locus*/
							printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                            printf("\n\n Display the probability of observed statistics for each locus obtained by coalescent Monte Carlo simulations: ");
                            printf("\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. Note: & means the ierations would be increased to confirm significance");
                            printf("\n Multiple hits not included in the analysis.\n");
                            if(file_output) {
								if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                fputs("\n\n Display the probability of observed statistics for each locus obtained by coalescent Monte Carlo simulations: ",file_output);
                                fprintf(file_output,"\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. Note: & means the ierations would be increased to confirm significance");
                                fputs("\n Multiple hits not included in the analysis.\n",file_output);
                            }
                            if(*outgroup) {
                                printf("\n\nloci\tObs(Tajima's_D)\tP(Tajima's_D)\tIBonf(Tajima's_D)\tObs(Fu&Li's_D)\tP(Fu&Li's_D)\tIBonf(Fu&Li's_D)\tObs(Fu&Li's_F)\tP(Fu&Li's_F)\tIBonf(Fu&Li's_F)\tObs(Fu's_Fs)\tP(Fu's_Fs)\tIBonf(Fu's_Fs)\tObs(Fay&Wu's_Hn)\tP(Fay&Wu's_Hn)\tIBonf(Fay&Wu's_Hn)\tObs(Fay&Wu's_H)\tP(Fay&Wu's_H)\tIBonf(Fay&Wu's_H)\tObs(Rozas'_ZA)\tP(Rozas'_ZA)\tIBonf(Rozas'_ZA)\tObs(Wall's_B)\tP(Wall's_B)\tIBonf(Wall's_B)\tObs(Wall's_Q)\tP(Wall's_Q)\tIBonf(Wall's_Q)\tObs(R2)\tP(R2)\tIBonf(R2)\tObs(ZengE)\tP(ZengE)\tIBonf(ZengE)\tObs(EW)\tP(EW)\tIBonf(EW)\n");
                                if(file_output) 
                                    fputs("\n\nloci\tObs(Tajima's_D)\tP(Tajima's_D)\tIBonf(Tajima's_D)\tObs(Fu&Li's_D)\tP(Fu&Li's_D)\tIBonf(Fu&Li's_D)\tObs(Fu&Li's_F)\tP(Fu&Li's_F)\tIBonf(Fu&Li's_F)\tObs(Fu's_Fs)\tP(Fu's_Fs)\tIBonf(Fu's_Fs)\tObs(Fay&Wu's_Hn)\tP(Fay&Wu's_Hn)\tIBonf(Fay&Wu's_Hn)\tObs(Fay&Wu's_H)\tP(Fay&Wu's_H)\tIBonf(Fay&Wu's_H)\tObs(Rozas'_ZA)\tP(Rozas'_ZA)\tIBonf(Rozas'_ZA)\tObs(Wall's_B)\tP(Wall's_B)\tIBonf(Wall's_B)\tObs(Wall's_Q)\tP(Wall's_Q)\tIBonf(Wall's_Q)\tObs(R2)\tP(R2)\tIBonf(R2)\tObs(ZengE)\tP(ZengE)\tIBonf(ZengE)\tObs(EW)\tP(EW)\tIBonf(EW)\n",file_output);
								
                                for(x=0;x<data[0].n_loci;x++) {
                                    nunderalfa[0] = nunderalfa[1] = nunderalfa[2] = nunderalfa[3] = nunderalfa[4] =
                                        nunderalfa[5] = nunderalfa[6] = nunderalfa[7] = nunderalfa[8] = nunderalfa[9] = nunderalfa[10] = nunderalfa[11] = 0;
                                    nequivalfa[0] = nequivalfa[1] = nequivalfa[2] = nequivalfa[3] = nequivalfa[4] = 
                                        nequivalfa[5] = nequivalfa[6] = nequivalfa[7] = nequivalfa[8] = nequivalfa[9] = nequivalfa[10] = nequivalfa[11] = 0;
                                    nitertests[0] = nitertests[1] = nitertests[2] = nitertests[3] = nitertests[4] = 
                                        nitertests[5] = nitertests[6] = nitertests[7] = nitertests[8] = nitertests[9] = nitertests[10] = nitertests[11] = 0;
                                    
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        if(matrix[x].tajimaD != (double)-10000 && matrixsim[x][it].tajimaD != (double)-10000) {
                                            if(matrix[x].tajimaD < matrixsim[x][it].tajimaD)
                                                nunderalfa[0] += 1;
                                            if(matrix[x].tajimaD == matrixsim[x][it].tajimaD)
                                                nequivalfa[0] += 1;
                                            nitertests[0] += 1;
                                        }
                                        if(matrix[x].fuliD != (double)-10000 && matrixsim[x][it].fuliD != (double)-10000) {
                                            if(matrix[x].fuliD < matrixsim[x][it].fuliD)
                                                nunderalfa[1] += 1;
                                            if(matrix[x].fuliD == matrixsim[x][it].fuliD)
                                                nequivalfa[1] += 1;
                                            nitertests[1] += 1;
                                        }
                                        if(matrix[x].fuliF != (double)-10000 && matrixsim[x][it].fuliF != (double)-10000) {
                                            if(matrix[x].fuliF < matrixsim[x][it].fuliF)
                                                nunderalfa[2] += 1;
                                            if(matrix[x].fuliF == matrixsim[x][it].fuliF)
                                                nequivalfa[2] += 1;
                                            nitertests[2] += 1;
                                        }
                                        if(matrix[x].fuFs != (double)-10000 && matrixsim[x][it].fuFs != (double)-10000) {
                                            if(matrix[x].fuFs < matrixsim[x][it].fuFs)
                                                nunderalfa[3] += 1;
                                            if(matrix[x].fuFs == matrixsim[x][it].fuFs)
                                                nequivalfa[3] += 1;
                                            nitertests[3] += 1;
                                        }
                                        if(matrix[x].faywuH != (double)-10000 && matrixsim[x][it].faywuH != (double)-10000) {
                                            if(matrix[x].faywuH < matrixsim[x][it].faywuH)
                                                nunderalfa[4] += 1;
                                            if(matrix[x].faywuH == matrixsim[x][it].faywuH)
                                                nequivalfa[4] += 1;
                                            nitertests[4] += 1;
                                        }
                                        if(matrix[x].faywuHo != (double)-10000 && matrixsim[x][it].faywuHo != (double)-10000) {
                                            if(matrix[x].faywuHo < matrixsim[x][it].faywuHo)
                                                nunderalfa[5] += 1;
                                            if(matrix[x].faywuHo == matrixsim[x][it].faywuHo)
                                                nequivalfa[5] += 1;
                                            nitertests[5] += 1;
                                        }
                                        if(matrix[x].rZA != (double)-10000 && matrixsim[x][it].rZA != (double)-10000) {
                                            if(matrix[x].rZA < matrixsim[x][it].rZA)
                                                nunderalfa[6] += 1;
                                            if(matrix[x].rZA == matrixsim[x][it].rZA)
                                                nequivalfa[6] += 1;
                                            nitertests[6] += 1;
                                        }
                                        if(matrix[x].wB != (double)-10000 && matrixsim[x][it].wB != (double)-10000) {
                                            if(matrix[x].wB < matrixsim[x][it].wB)
                                                nunderalfa[7] += 1;
                                            if(matrix[x].wB == matrixsim[x][it].wB)
                                                nequivalfa[7] += 1;
                                            nitertests[7] += 1;
                                        }
                                        if(matrix[x].wQ != (double)-10000 && matrixsim[x][it].wQ != (double)-10000) {
                                            if(matrix[x].wQ < matrixsim[x][it].wQ)
                                                nunderalfa[8] += 1;
                                            if(matrix[x].wQ == matrixsim[x][it].wQ)
                                                nequivalfa[8] += 1;
                                            nitertests[8] += 1;
                                        }
                                        if(matrix[x].R2 != (double)-10000 && matrixsim[x][it].wQ != (double)-10000) {
                                            if(matrix[x].R2 < matrixsim[x][it].R2)
                                                nunderalfa[9] += 1;
                                            if(matrix[x].R2 == matrixsim[x][it].R2)
                                                nequivalfa[9] += 1;
                                            nitertests[9] += 1;
                                        }
                                        if(matrix[x].zengE != (double)-10000 && matrixsim[x][it].zengE != (double)-10000) {
                                            if(matrix[x].zengE < matrixsim[x][it].zengE)
                                                nunderalfa[10] += 1;
                                            if(matrix[x].zengE == matrixsim[x][it].zengE)
                                                nequivalfa[10] += 1;
                                            nitertests[10] += 1;
                                        }
                                        if(matrix[x].ewtest != (double)-10000 && matrixsim[x][it].ewtest != (double)-10000) {
                                            if(matrix[x].ewtest < matrixsim[x][it].ewtest)
                                                nunderalfa[11] += 1;
                                            if(matrix[x].ewtest == matrixsim[x][it].ewtest)
                                                nequivalfa[11] += 1;
                                            nitertests[11] += 1;
                                        }
                                    }
                                    if(nitertests[0]) {
                                        ibonf[0][x] = (double)nunderalfa[0]/(double)nitertests[0];
										if((double)1./((double)nitertests[0]) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                        if(ibonf[0][x] > (double)0.5) ibonf[0][x] -= (double)1./((double)nitertests[0]) *(double)nequivalfa[0];
                                        else ibonf[0][x] += (double)1./((double)nitertests[0]) *(double)nequivalfa[0];
                                    }
									else ibonf[0][x] = (double)-10000;
                                    if(nitertests[1]) {
                                        ibonf[1][x] = (double)nunderalfa[1]/(double)nitertests[1];
										if((double)1./((double)nitertests[1]) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                        if(ibonf[1][x] > (double)0.5) ibonf[1][x] -= (double)1./((double)nitertests[1]) *(double)nequivalfa[1];
                                        else ibonf[1][x] += (double)1./((double)nitertests[1]) *(double)nequivalfa[1];
                                    }
									else ibonf[1][x] = (double)-10000;
                                    if(nitertests[2]) {
                                        ibonf[2][x] = (double)nunderalfa[2]/(double)nitertests[2];
										if((double)1./((double)nitertests[2]) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                        if(ibonf[2][x] > (double)0.5) ibonf[2][x] -= (double)1./((double)nitertests[2]) *(double)nequivalfa[2];
                                        else ibonf[2][x] += (double)1./((double)nitertests[2]) *(double)nequivalfa[2];
                                    }
									else ibonf[2][x] = (double)-10000;
                                    if(nitertests[3]) {
                                        ibonf[3][x] = (double)nunderalfa[3]/(double)nitertests[3];
										if((double)1./((double)nitertests[3]) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                        if(ibonf[3][x] > (double)0.5) ibonf[3][x] -= (double)1./((double)nitertests[3]) *(double)nequivalfa[3];
                                        else ibonf[3][x] += (double)1./((double)nitertests[3]) *(double)nequivalfa[3];
                                    }
									else ibonf[3][x] = (double)-10000;
                                    if(nitertests[4]) {
                                        ibonf[4][x] = (double)nunderalfa[4]/(double)nitertests[4];
										if((double)1./((double)nitertests[4]) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                        if(ibonf[4][x] > (double)0.5) ibonf[4][x] -= (double)1./((double)nitertests[4]) *(double)nequivalfa[4];
                                        else ibonf[4][x] += (double)1./((double)nitertests[4]) *(double)nequivalfa[4];
                                    }
									else ibonf[4][x] = (double)-10000;
                                    if(nitertests[5]) {
                                        ibonf[5][x] = (double)nunderalfa[5]/(double)nitertests[5];
										if((double)1./((double)nitertests[5]) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                        if(ibonf[5][x] > (double)0.5) ibonf[5][x] -= (double)1./((double)nitertests[5]) *(double)nequivalfa[5];
                                        else ibonf[5][x] += (double)1./((double)nitertests[5]) *(double)nequivalfa[5];
                                    }
									else ibonf[5][x] = (double)-10000;
                                    if(nitertests[6]) {
                                        ibonf[6][x] = (double)nunderalfa[6]/(double)nitertests[6];
										if((double)1./((double)nitertests[6]) *(double)nequivalfa[6] >= (double)0.5) nequivalfa[6] /= 2;
                                        if(ibonf[6][x] > (double)0.5) ibonf[6][x] -= (double)1./((double)nitertests[6]) *(double)nequivalfa[6];
                                        else ibonf[6][x] += (double)1./((double)nitertests[6]) *(double)nequivalfa[6];
                                    }
									else ibonf[6][x] = (double)-10000;
                                    if(nitertests[7]) {
                                        ibonf[7][x] = (double)nunderalfa[7]/(double)nitertests[7];
										if((double)1./((double)nitertests[7]) *(double)nequivalfa[7] >= (double)0.5) nequivalfa[7] /= 2;
                                        if(ibonf[7][x] > (double)0.5) ibonf[7][x] -= (double)1./((double)nitertests[7]) *(double)nequivalfa[7];
                                        else ibonf[7][x] += (double)1./((double)nitertests[7]) *(double)nequivalfa[7];
                                    }
									else ibonf[7][x] = (double)-10000;
                                    if(nitertests[8]) {
                                        ibonf[8][x] = (double)nunderalfa[8]/(double)nitertests[8];
										if((double)1./((double)nitertests[8]) *(double)nequivalfa[8] >= (double)0.5) nequivalfa[8] /= 2;
                                        if(ibonf[8][x] > (double)0.5) ibonf[8][x] -= (double)1./((double)nitertests[8]) *(double)nequivalfa[8];
                                        else ibonf[8][x] += (double)1./((double)nitertests[8]) *(double)nequivalfa[8];
                                    }
									else ibonf[8][x] = (double)-10000;
                                    if(nitertests[9]) {
                                        ibonf[9][x] = (double)nunderalfa[9]/(double)nitertests[9];
										if((double)1./((double)nitertests[9]) *(double)nequivalfa[9] >= (double)0.5) nequivalfa[9] /= 2;
                                        if(ibonf[9][x] > (double)0.5) ibonf[9][x] -= (double)1./((double)nitertests[9]) *(double)nequivalfa[9];
                                        else ibonf[9][x] += (double)1./((double)nitertests[9]) *(double)nequivalfa[9];
                                    }
									else ibonf[9][x] = (double)-10000;
                                    if(nitertests[10]) {
                                        ibonf[10][x] = (double)nunderalfa[10]/(double)nitertests[10];
										if((double)1./((double)nitertests[10]) *(double)nequivalfa[10] >= (double)0.5) nequivalfa[10] /= 2;
                                        if(ibonf[10][x] > (double)0.5) ibonf[10][x] -= (double)1./((double)nitertests[10]) *(double)nequivalfa[10];
                                        else ibonf[10][x] += (double)1./((double)nitertests[10]) *(double)nequivalfa[10];
                                    }
									else ibonf[10][x] = (double)-10000;
                                    if(nitertests[11]) {
                                        ibonf[11][x] = (double)nunderalfa[11]/(double)nitertests[11];
										if((double)1./((double)nitertests[11]) *(double)nequivalfa[11] >= (double)0.5) nequivalfa[11] /= 2;
                                        if(ibonf[11][x] > (double)0.5) ibonf[11][x] -= (double)1./((double)nitertests[11]) *(double)nequivalfa[11];
                                        else ibonf[11][x] += (double)1./((double)nitertests[11]) *(double)nequivalfa[11];
                                    }
									else ibonf[11][x] = (double)-10000;
                                }
								ibonf_sort(ibonf,ibonfn,12,data[0].n_loci);
								
								for(x=0;x<data[0].n_loci;x++) {
                                    printf("%d:%s\t",x,matrix[x].gene);
                                    if(file_output) fprintf(file_output,"%d:%s\t",x,matrix[x].gene);

                                    if(ibonfn[0][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].tajimaD,ibonf[0][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].tajimaD,ibonf[0][x]);
										print_ibonf(ibonf[0][x],ibonfn[0][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[1][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].fuliD,ibonf[1][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].fuliD,ibonf[1][x]);
										print_ibonf(ibonf[1][x],ibonfn[1][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[2][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].fuliF,ibonf[2][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].fuliF,ibonf[2][x]);
										print_ibonf(ibonf[2][x],ibonfn[2][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[3][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].fuFs,ibonf[3][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].fuFs,ibonf[3][x]);
										print_ibonf(ibonf[3][x],ibonfn[3][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[4][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].faywuH,ibonf[4][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].faywuH,ibonf[4][x]);
										print_ibonf(ibonf[4][x],ibonfn[4][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[5][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].faywuHo,ibonf[5][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].faywuHo,ibonf[5][x]);
										print_ibonf(ibonf[5][x],ibonfn[5][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[5][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].rZA,ibonf[6][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].rZA,ibonf[6][x]);
										print_ibonf(ibonf[6][x],ibonfn[6][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[7][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].wB,ibonf[7][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].wB,ibonf[7][x]);
										print_ibonf(ibonf[7][x],ibonfn[7][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[7][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].wQ,ibonf[8][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].wQ,ibonf[8][x]);
										print_ibonf(ibonf[8][x],ibonfn[8][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[9][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].R2,ibonf[9][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].R2,ibonf[9][x]);
										print_ibonf(ibonf[9][x],ibonfn[9][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[10][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].zengE,ibonf[10][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].zengE,ibonf[10][x]);
										print_ibonf(ibonf[10][x],ibonfn[10][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[11][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].ewtest,ibonf[11][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].ewtest,ibonf[11][x]);
										print_ibonf(ibonf[11][x],ibonfn[11][x],data[0].n_loci,data[0].n_iter,file_output);
										printf("\n");
										if(file_output) fprintf(file_output,"\n");
                                    }
                                    else{
                                        printf("na\tna\tna\n");
                                        if(file_output) fprintf(file_output,"na\tna\tna\n");
                                    }
								}
                            }
                            else {
                                printf("\n\nloci\tObs(Tajima's_D)\tP(Tajima's_D)\tIBonf(Tajima's_D)\tObs(Fu&Li's_D*)\tP(Fu&Li's_D*)\tIBonf(Fu&Li's_D*)\tObs(Fu&Li's_F*)\tP(Fu&Li's_F*)\tIBonf(Fu&Li's_F*)\tObs(Fu's_Fs)\tP(Fu's_Fs)\tIBonf(Fu's_Fs)\tObs(Rozas'_ZA)\tP(Rozas'_ZA)\tIBonf(Rozas'_ZA)\tObs(Wall's_B)\tP(Wall's_B)\tIBonf(Wall's_B)\tObs(Wall's_Q)\tP(Wall's_Q)\tIBonf(Wall's_Q)\tObs(R2)\tP(R2)\tIBonf(R2)\tObs(EW)\tP(EW)\tIBonf(EW)\n");
                                if(file_output) 
                                    fputs("\n\nloci\tObs(Tajima's_D)\tP(Tajima's_D)\tIBonf(Tajima's_D)\tObs(Fu&Li's_D*)\tP(Fu&Li's_D*)\tIBonf(Fu&Li's_D*)\tObs(Fu&Li's_F*)\tP(Fu&Li's_F*)\tIBonf(Fu&Li's_F*)\tObs(Fu's_Fs)\tP(Fu's_Fs)\tIBonf(Fu's_Fs)\tObs(Rozas'_ZA)\tP(Rozas'_ZA)\tIBonf(Rozas'_ZA)\tObs(Wall's_B)\tP(Wall's_B)\tIBonf(Wall's_B)\tObs(Wall's_Q)\tP(Wall's_Q)\tIBonf(Wall's_Q)\tObs(R2)\tP(R2)\tIBonf(R2)\tObs(EW)\tP(EW)\tIBonf(EW)\n",file_output);

                                for(x=0;x<data[0].n_loci;x++) {
                                    nunderalfa[0] = nunderalfa[1] = nunderalfa[2] = nunderalfa[3] = nunderalfa[4] =
                                        nunderalfa[5] = nunderalfa[6] = nunderalfa[7]  = nunderalfa[8] = 0;
                                    nequivalfa[0] = nequivalfa[1] = nequivalfa[2] = nequivalfa[3] = nequivalfa[4] =
                                        nequivalfa[5] = nequivalfa[6] = nequivalfa[7] = nequivalfa[8] = 0;
                                    nitertests[0] = nitertests[1] = nitertests[2] = nitertests[3] = nitertests[4] =
                                        nitertests[5] = nitertests[6] = nitertests[7] = nitertests[8] = 0;
                                        
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        if(matrix[x].tajimaD != (double)-10000 && matrixsim[x][it].tajimaD != (double)-10000) {
                                            if(matrix[x].tajimaD < matrixsim[x][it].tajimaD)
                                                nunderalfa[0] += 1;
                                            if(matrix[x].tajimaD == matrixsim[x][it].tajimaD)
                                                nequivalfa[0] += 1;
                                            nitertests[0] += 1;
                                        }
                                        if(matrix[x].fuliDn != (double)-10000 && matrixsim[x][it].fuliDn != (double)-10000) {
                                            if(matrix[x].fuliDn < matrixsim[x][it].fuliDn)
                                                nunderalfa[1] += 1;
                                            if(matrix[x].fuliDn == matrixsim[x][it].fuliDn)
                                                nequivalfa[1] += 1;
                                            nitertests[1] += 1;
                                        }
                                        if(matrix[x].fuliFn != (double)-10000 && matrixsim[x][it].fuliFn != (double)-10000) {
                                            if(matrix[x].fuliFn < matrixsim[x][it].fuliFn)
                                                nunderalfa[2] += 1;
                                            if(matrix[x].fuliFn == matrixsim[x][it].fuliFn)
                                                nequivalfa[2] += 1;
                                            nitertests[2] += 1;
                                        }
                                        if(matrix[x].fuFs != (double)-10000 && matrixsim[x][it].fuFs != (double)-10000) {
                                            if(matrix[x].fuFs < matrixsim[x][it].fuFs)
                                                nunderalfa[3] += 1;
                                            if(matrix[x].fuFs == matrixsim[x][it].fuFs)
                                                nequivalfa[3] += 1;
                                            nitertests[3] += 1;
                                        }
                                        if(matrix[x].rZA != (double)-10000 && matrixsim[x][it].rZA != (double)-10000) {
                                            if(matrix[x].rZA < matrixsim[x][it].rZA)
                                                nunderalfa[4] += 1;
                                            if(matrix[x].rZA == matrixsim[x][it].rZA)
                                                nequivalfa[4] += 1;
                                            nitertests[4] += 1;
                                        }
                                        if(matrix[x].wB != (double)-10000 && matrixsim[x][it].wB != (double)-10000) {
                                            if(matrix[x].wB < matrixsim[x][it].wB)
                                                nunderalfa[5] += 1;
                                            if(matrix[x].wB == matrixsim[x][it].wB)
                                                nequivalfa[5] += 1;
                                            nitertests[5] += 1;
                                        }
                                        if(matrix[x].wQ != (double)-10000 && matrixsim[x][it].wQ != (double)-10000) {
                                            if(matrix[x].wQ < matrixsim[x][it].wQ)
                                                nunderalfa[6] += 1;
                                            if(matrix[x].wQ == matrixsim[x][it].wQ)
                                                nequivalfa[6] += 1;
                                            nitertests[6] += 1;
                                        }
                                        if(matrix[x].R2 != (double)-10000 && matrixsim[x][it].R2 != (double)-10000) {
                                            if(matrix[x].R2 < matrixsim[x][it].R2)
                                                nunderalfa[7] += 1;
                                            if(matrix[x].R2 == matrixsim[x][it].R2)
                                                nequivalfa[7] += 1;
                                            nitertests[7] += 1;
                                        }
                                        if(matrix[x].ewtest != (double)-10000 && matrixsim[x][it].ewtest != (double)-10000) {
                                            if(matrix[x].ewtest < matrixsim[x][it].ewtest)
                                                nunderalfa[8] += 1;
                                            if(matrix[x].ewtest == matrixsim[x][it].ewtest)
                                                nequivalfa[8] += 1;
                                            nitertests[8] += 1;
                                        }
                                    }
                                    if(nitertests[0]) {
                                        ibonf[0][x] = (double)nunderalfa[0]/(double)nitertests[0];
										if((double)1./((double)nitertests[0]) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                        if(ibonf[0][x] > (double)0.5) ibonf[0][x] -= (double)1./((double)nitertests[0]) *(double)nequivalfa[0];
                                        else ibonf[0][x] += (double)1./((double)nitertests[0]) *(double)nequivalfa[0];
                                    }
									else ibonf[0][x] = (double)-10000;
                                    if(nitertests[1]) {
                                        ibonf[1][x] = (double)nunderalfa[1]/(double)nitertests[1];
										if((double)1./((double)nitertests[1]) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                        if(ibonf[1][x] > (double)0.5) ibonf[1][x] -= (double)1./((double)nitertests[1]) *(double)nequivalfa[1];
                                        else ibonf[1][x] += (double)1./((double)nitertests[1]) *(double)nequivalfa[1];
                                    }
									else ibonf[1][x] = (double)-10000;
                                    if(nitertests[2]) {
                                        ibonf[2][x] = (double)nunderalfa[2]/(double)nitertests[2];
										if((double)1./((double)nitertests[2]) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                        if(ibonf[2][x] > (double)0.5) ibonf[2][x] -= (double)1./((double)nitertests[2]) *(double)nequivalfa[2];
                                        else ibonf[2][x] += (double)1./((double)nitertests[2]) *(double)nequivalfa[2];
                                    }
									else ibonf[2][x] = (double)-10000;
                                    if(nitertests[3]) {
                                        ibonf[3][x] = (double)nunderalfa[3]/(double)nitertests[3];
										if((double)1./((double)nitertests[3]) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                        if(ibonf[3][x] > (double)0.5) ibonf[3][x] -= (double)1./((double)nitertests[3]) *(double)nequivalfa[3];
                                        else ibonf[3][x] += (double)1./((double)nitertests[3]) *(double)nequivalfa[3];
                                    }
									else ibonf[3][x] = (double)-10000;
                                    if(nitertests[4]) {
                                        ibonf[4][x] = (double)nunderalfa[4]/(double)nitertests[4];
										if((double)1./((double)nitertests[4]) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                        if(ibonf[4][x] > (double)0.5) ibonf[4][x] -= (double)1./((double)nitertests[4]) *(double)nequivalfa[4];
                                        else ibonf[4][x] += (double)1./((double)nitertests[4]) *(double)nequivalfa[4];
                                    }
									else ibonf[4][x] = (double)-10000;
                                    if(nitertests[5]) {
                                        ibonf[5][x] = (double)nunderalfa[5]/(double)nitertests[5];
										if((double)1./((double)nitertests[5]) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                        if(ibonf[5][x] > (double)0.5) ibonf[5][x] -= (double)1./((double)nitertests[5]) *(double)nequivalfa[5];
                                        else ibonf[5][x] += (double)1./((double)nitertests[5]) *(double)nequivalfa[5];
                                    }
									else ibonf[5][x] = (double)-10000;
                                    if(nitertests[6]) {
                                        ibonf[6][x] = (double)nunderalfa[6]/(double)nitertests[6];
										if((double)1./((double)nitertests[6]) *(double)nequivalfa[6] >= (double)0.5) nequivalfa[6] /= 2;
                                        if(ibonf[6][x] > (double)0.5) ibonf[6][x] -= (double)1./((double)nitertests[6]) *(double)nequivalfa[6];
                                        else ibonf[6][x] += (double)1./((double)nitertests[6]) *(double)nequivalfa[6];
                                    }
									else ibonf[6][x] = (double)-10000;
                                    if(nitertests[7]) {
                                        ibonf[7][x] = (double)nunderalfa[7]/(double)nitertests[7];
										if((double)1./((double)nitertests[7]) *(double)nequivalfa[7] >= (double)0.5) nequivalfa[7] /= 2;
                                        if(ibonf[7][x] > (double)0.5) ibonf[7][x] -= (double)1./((double)nitertests[7]) *(double)nequivalfa[7];
                                        else ibonf[7][x] += (double)1./((double)nitertests[7]) *(double)nequivalfa[7];
                                    }
									else ibonf[7][x] = (double)-10000;
                                    if(nitertests[8]) {
                                        ibonf[8][x] = (double)nunderalfa[8]/(double)nitertests[8];
										if((double)1./((double)nitertests[8]) *(double)nequivalfa[8] >= (double)0.5) nequivalfa[8] /= 2;
                                        if(ibonf[8][x] > (double)0.5) ibonf[8][x] -= (double)1./((double)nitertests[8]) *(double)nequivalfa[8];
                                        else ibonf[8][x] += (double)1./((double)nitertests[8]) *(double)nequivalfa[8];
                                    }
									else ibonf[8][x] = (double)-10000;
                                }
								ibonf_sort(ibonf,ibonfn,9,data[0].n_loci);
								
								for(x=0;x<data[0].n_loci;x++) {
                                    printf("%d:%s\t",x,matrix[x].gene);
                                    if(file_output) fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
                                    
									if(ibonfn[0][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].tajimaD,ibonf[0][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].tajimaD,ibonf[0][x]);
										print_ibonf(ibonf[0][x],ibonfn[0][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[1][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].fuliDn,ibonf[1][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].fuliDn,ibonf[1][x]);
										print_ibonf(ibonf[1][x],ibonfn[1][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[2][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].fuliFn,ibonf[2][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].fuliFn,ibonf[2][x]);
										print_ibonf(ibonf[2][x],ibonfn[2][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[3][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].fuFs,ibonf[3][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].fuFs,ibonf[3][x]);
										print_ibonf(ibonf[3][x],ibonfn[3][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[4][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].rZA,ibonf[4][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].rZA,ibonf[4][x]);
										print_ibonf(ibonf[4][x],ibonfn[4][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[5][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].wB,ibonf[5][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].wB,ibonf[5][x]);
										print_ibonf(ibonf[5][x],ibonfn[5][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[6][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].wQ,ibonf[6][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].wQ,ibonf[6][x]);
										print_ibonf(ibonf[6][x],ibonfn[6][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[7][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].R2,ibonf[7][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].R2,ibonf[7][x]);
										print_ibonf(ibonf[7][x],ibonfn[7][x],data[0].n_loci,data[0].n_iter,file_output);
                                    }
                                    else{
                                        printf("na\tna\tna\t");
                                        if(file_output) fprintf(file_output,"na\tna\tna\t");
                                    }
                                    if(ibonfn[8][x] != (double)-1) {
                                        printf("%g\t%g\t",matrix[x].ewtest,ibonf[8][x]);
                                        if(file_output) fprintf(file_output,"%g\t%g\t",matrix[x].ewtest,ibonf[8][x]);
										print_ibonf(ibonf[8][x],ibonfn[8][x],data[0].n_loci,data[0].n_iter,file_output);
										printf("\n");
										if(file_output) fprintf(file_output,"\n");
                                    }
                                    else{
                                        printf("na\tna\tna\n");
                                        if(file_output) fprintf(file_output,"na\tna\tna\n");
                                    }
								}
							}
                        }
                        else {
                            printf("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n");
                            if(file_output)
                                fputs("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n",file_output);
                        }
                    }
					#if ALLRESULTS == 0
                    else {
					#endif
                        if(!onlymulo) {
                            /*calculate median, 10% 5% etc./locus*/
                            if((sortvector = (double *) malloc(data[0].n_iter * sizeof(double))) == 0) {
                                printf("Probabilities can not be calculated, sorry.\n");
                                if(file_output) fputs("Probabilities can not be calculated, sorry.\n",file_output);
                                break;
                            }
							printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                            printf("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n");
                            if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
							if(file_output)
                                fputs("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n",file_output);
                            if(*outgroup) {    
                                for(x=0;x<data[0].n_loci;x++) {
                                    printf("\nloci\tstatistic\tmedian\tavg\tvar\t");                                    
                                    if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                    if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                    if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                    if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                    if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                    if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                    if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                    printf("\n");
                                    if(file_output) {
                                        fputs("\nloci\tstatistic\tmedian\tavg\tvar\t",file_output);                                    
                                        if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                        if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                        if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                        if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                        if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                        if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                        if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                        fprintf(file_output,"\n");
                                    }

                                    /*Tajima's D */
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].tajimaD != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].tajimaD;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].tajimaD;
                                            var += (double)matrixsim[x][it].tajimaD * (double)matrixsim[x][it].tajimaD;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tTajima's_D\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tTajima's_D\t%g\t%g\t",
                                                x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tTajima's_D\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tTajima's_D\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*Fu&Li'sD*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].fuliD != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].fuliD;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].fuliD;
                                            var += (double)matrixsim[x][it].fuliD * (double)matrixsim[x][it].fuliD;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tFu&Li's_D\t%g\t%g\t",
                                            x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tFu&Li's_D\t%g\t%g\t",
                                                x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tFu&Li's_D\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tFu&Li's_D\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Fu&Li'sF*/
                                    
									/*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].fuliF != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].fuliF;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].fuliF;
                                            var += (double)matrixsim[x][it].fuliF * (double)matrixsim[x][it].fuliF;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tFu&Li's_F\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tFu&Li's_F\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tFu&Li's_F\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tFu&Li's_F\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*Fu'sFs*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].fuFs != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].fuFs;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].fuFs;
                                            var += (double)matrixsim[x][it].fuFs * (double)matrixsim[x][it].fuFs;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tFu's_Fs\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tFu's_Fs\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tFu's_Fs\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tFu's_Fs\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*FayWuH*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].faywuH != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].faywuH;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].faywuH;
                                            var += (double)matrixsim[x][it].faywuH * (double)matrixsim[x][it].faywuH;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tFay&Wu's_Hn\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tFay&Wu's_Hn\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tFay&Wu's_Hn\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tFay&Wu's_Hn\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*FayWuHo*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].faywuHo != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].faywuHo;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].faywuHo;
                                            var += (double)matrixsim[x][it].faywuHo * (double)matrixsim[x][it].faywuHo;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tFay&Wu's_H\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tFay&Wu's_H\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tFay&Wu's_H\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tFay&Wu's_H\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*RozasZA*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].rZA != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].rZA;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].rZA;
                                            var += (double)matrixsim[x][it].rZA * (double)matrixsim[x][it].rZA;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tRozas'_ZA\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tRozas'_ZA\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tRozas'_ZA\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tRozas'_ZA\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*Wall's_B*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].wB != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].wB;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].wB;
                                            var += (double)matrixsim[x][it].wB * (double)matrixsim[x][it].wB;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tWall's_B\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tWall's_B\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tWall's_B\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tWall's_B\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*Wall's_Q*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].wQ != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].wQ;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].wQ;
                                            var += (double)matrixsim[x][it].wQ * (double)matrixsim[x][it].wQ;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tWall's_Q\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tWall's_Q\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tWall's_Q\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tWall's_Q\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*R2*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].R2 != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].R2;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].R2;
                                            var += (double)matrixsim[x][it].R2 * (double)matrixsim[x][it].R2;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tR2\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tR2\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tR2\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tR2\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*ZengE*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].zengE != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].zengE;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].zengE;
                                            var += (double)matrixsim[x][it].zengE * (double)matrixsim[x][it].zengE;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tZeng_E\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tZeng_E\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tZeng_E\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tZeng_E\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                   
									/*EW*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].ewtest != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].ewtest;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].ewtest;
                                            var += (double)matrixsim[x][it].ewtest * (double)matrixsim[x][it].ewtest;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tEW\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tEW\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tEW\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tEW\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                }
                            }
                            else {
                                for(x=0;x<data[0].n_loci;x++) {
                                    printf("\nloci\tstatistic\tmedian\tavg\tvar\t");                                    
                                    if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                    if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                    if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                    if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                    if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                    if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                    if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                    printf("\n");
                                    if(file_output) {
                                        fputs("\nloci\tstatistic\tmedian\tavg\tvar\t",file_output);                                    
                                        if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                        if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                        if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                        if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                        if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                        if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                        if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                        fprintf(file_output,"\n");
                                    }

                                    /*Tajima's D */
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].tajimaD != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].tajimaD;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].tajimaD;
                                            var += (double)matrixsim[x][it].tajimaD * (double)matrixsim[x][it].tajimaD;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tTajima's_D\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tTajima's_D\t%g\t%g\t",
                                                x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tTajima's_D\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tTajima's_D\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*Fu&Li'sD**/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].fuliDn != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].fuliDn;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].fuliDn;
                                            var += (double)matrixsim[x][it].fuliDn * (double)matrixsim[x][it].fuliDn;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tFu&Li's_D*\t%g\t%g\t",
                                            x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tFu&Li's_D*\t%g\t%g\t",
                                                x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tFu&Li's_D*\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tFu&Li's_D*\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*Fu&Li'sF**/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].fuliFn != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].fuliFn;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].fuliFn;
                                            var += (double)matrixsim[x][it].fuliFn * (double)matrixsim[x][it].fuliFn;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tFu&Li's_F*\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tFu&Li's_F*\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tFu&Li's_F*\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tFu&Li's_F*\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                   
									/*Fu'sFs*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].fuFs != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].fuFs;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].fuFs;
                                            var += (double)matrixsim[x][it].fuFs * (double)matrixsim[x][it].fuFs;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tFu's_Fs\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tFu's_Fs\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tFu's_Fs\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tFu's_Fs\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*RozasZA*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].rZA != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].rZA;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].rZA;
                                            var += (double)matrixsim[x][it].rZA * (double)matrixsim[x][it].rZA;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tRozas'_ZA\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tRozas'_ZA\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tRozas'_ZA\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tRozas'_ZA\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*Wall's_B*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].wB != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].wB;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].wB;
                                            var += (double)matrixsim[x][it].wB * (double)matrixsim[x][it].wB;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tWall's_B\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tWall's_B\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tWall's_B\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tWall's_B\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*Wall's_Q*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].wQ != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].wQ;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].wQ;
                                            var += (double)matrixsim[x][it].wQ * (double)matrixsim[x][it].wQ;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tWall's_Q\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tWall's_Q\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tWall's_Q\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tWall's_Q\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*R2*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].R2 != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].R2;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].R2;
                                            var += (double)matrixsim[x][it].R2 * (double)matrixsim[x][it].R2;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tR2\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tR2\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tR2\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tR2\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    
									/*EW*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixsim[x][it].ewtest != (double) -10000) {
                                            sortvector[it2] = (double)matrixsim[x][it].ewtest;
                                            it2 += 1;
                                            avg += (double)matrixsim[x][it].ewtest;
                                            var += (double)matrixsim[x][it].ewtest * (double)matrixsim[x][it].ewtest;
                                        }
                                    }
                                    if(it2 > 0) {
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("%d\tEW\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"%d\tEW\t%g\t%g\t",x,sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    else {
                                        printf("%d\tEW\tna\tna\tna\t",x);
                                        if(file_output) 
                                            fprintf(file_output,"%d\tEW\tna\tna\tna\t",x);
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                }
                            }
                            free(sortvector);
                        }
                        else {
                            printf("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n");
                            if(file_output)
                                fputs("\n Sorry, probability of detailed statistics are not available if data for each loci are not in the memory.\n",file_output);
                        }
					#if ALLRESULTS == 0
                    }
					#endif
					for(x=0;x<13;x++) {
						free(ibonf[x]);
						free(ibonfn[x]);
					}
					free(ibonf);
					free(ibonfn);
                    break;
			}
            break;
        case 1: /*multilocus*/
            switch(neuttest) {
                case 0: /* statistics, multilocus*/
                    #if COMMAND_LINE
                    printf("\n\nMENU 3. Statistical inference based on Coalescent Monte Carlo simulations menu:");
                    printf("\n     3.1. Display coalescent simulation analysis menu:");
                    printf("\n     3.1.0. Display statistics menu:");
                    printf("\n     3.1.0.0. Display multilocus statistics menu:");
                    printf("\n     3.1.0.0.1. Display table with detailed statistics for multilocus data:\n\n");
                    
                    printf("CHOOSE:\n");
                    printf(" 0 - Display general statistics.\n");
                    printf(" 1 - Display other linkage related statistics.\n");
                    printf(" 2 - Display estimates of total locus variability.\n");
                    printf(" 3 - Display estimates of nucleotide locus variability.\n");
                    printf(" 4 - Back to previous menu.\n");

                    if(file_output) {
                        fprintf(file_output,"\n\n     MENU:\n     3. Statistical inference based on Coalescent Monte Carlo simulations menu:");
                        fprintf(file_output,"\n     3.1. Display coalescent simulation analysis menu:");
                        fprintf(file_output,"\n     3.1.0. Display statistics menu:");
                        fprintf(file_output,"\n     3.1.0.0. Display multilocus statistics menu:");
                        fprintf(file_output,"\n     3.1.0.0.1. Display table with detailed statistics for multilocus data:\n\n");
                        
                        fprintf(file_output," 0 - Display general statistics.\n");
                        fprintf(file_output," 1 - Display other linkage related statistics.\n");
                        fprintf(file_output," 2 - Display estimates of total locus variability.\n");
                        fprintf(file_output," 3 - Display estimates of nucleotide locus variability.\n");
                        fprintf(file_output," 4 - Back to previous menu.\n");
                    }
                    #endif
                    
                    do *k = getchar();
                    while(*k<'0' || *k>'4');

                    if(file_output) fprintf(file_output,"\nOPTION CHOSEN: %c\n\n",*k);

                    switch(*k) {
                        case '0': /* 0 - Display general statistics*/
                            if(*observed_data && dataobsequalsim) {/*compare simulations with observed data.*/
                                if(*outgroup) {
									maxnoutgroups = 0;
									for(x=0;x<data[0].n_loci;x++) {
										if(matrix[x].noutgroups > maxnoutgroups) 
											maxnoutgroups = matrix[x].noutgroups;
									} 
                                    /*average and variance comparisons.*/
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ");
                                    printf("\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%.");
                                    printf("\n Multiple hits not included in the analysis.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ",file_output);
                                        fprintf(file_output,"\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%.");
                                        fputs("\n Multiple hits not included in the analysis.\n",file_output);
                                    }
                                    if(data[0].time_spec <= 0.) {
                                        nunderalfa[0] = nunderalfa[1] = nunderalfa[2] = nunderalfa[3] =  nunderalfa[4] =  nunderalfa[5] =  nunderalfa[6] =  nunderalfa[7] =  nunderalfa[8] =  nunderalfa[9] =  nunderalfa[10] =  nunderalfa[11] = 0;
                                        nequivalfa[0] = nequivalfa[1] = nequivalfa[2] = nequivalfa[3] =  nequivalfa[4] =  nequivalfa[5] =  nequivalfa[6] =  nequivalfa[7] =  nequivalfa[8] =  nequivalfa[9] =  nequivalfa[10] =  nequivalfa[11] = 0;
                                        varobs[0] = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,data[0].n_loci);
                                        varobs[1] = variances((double)matrixml[0].S2nhapl,(double)matrixml[0].Snhapl,data[0].n_loci);                                        					
                                        varobs[2] = variances((double)matrixml[0].S2nhaplsam,(double)matrixml[0].Snhaplsam,data[0].n_loci);                                        					
                                        varobs[3] = variances((double)matrixml[0].S2hapldiv,(double)matrixml[0].Shapldiv,data[0].n_loci);                                        					
                                        varobs[4] = variances((double)matrixml[0].S2Rvpi,(double)matrixml[0].SRvpi,data[0].n_loci);                                        					
                                        varobs[5] = variances((double)matrixml[0].S2Rm,(double)matrixml[0].SRm,data[0].n_loci);                                        					
                                        for(it=0;it<(long int)data[0].n_iter;it++) {
                                            if(matrixml[0].Sbiallsites < matrixmlsim[it].Sbiallsites)
                                                nunderalfa[0] += 1;
                                            if(varobs[0] < variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci))
                                                nunderalfa[1] += 1;
                                            if(matrixml[0].Snhapl < matrixmlsim[it].Snhapl)
                                                nunderalfa[2] += 1;
                                            if(varobs[1] < variances((double)matrixmlsim[it].S2nhapl,(double)matrixmlsim[it].Snhapl,data[0].n_loci))
                                                nunderalfa[3] += 1;
                                            if(matrixml[0].Snhaplsam < matrixmlsim[it].Snhaplsam)
                                                nunderalfa[4] += 1;
                                            if(varobs[2] < variances((double)matrixmlsim[it].S2nhaplsam,(double)matrixmlsim[it].Snhaplsam,data[0].n_loci))
                                                nunderalfa[5] += 1;
                                            if(matrixml[0].Shapldiv < matrixmlsim[it].Shapldiv)
                                                nunderalfa[6] += 1;
                                            if(varobs[3] < variances((double)matrixmlsim[it].S2hapldiv,(double)matrixmlsim[it].Shapldiv,data[0].n_loci))
                                                nunderalfa[7] += 1;
                                            if(matrixml[0].SRvpi < matrixmlsim[it].SRvpi)
                                                nunderalfa[8] += 1;
                                            if(varobs[4] < variances((double)matrixmlsim[it].S2Rvpi,(double)matrixmlsim[it].SRvpi,data[0].n_loci))
                                                nunderalfa[9] += 1;
                                            if(matrixml[0].SRm < matrixmlsim[it].SRm)
                                                nunderalfa[10] += 1;
                                            if(varobs[5] < variances((double)matrixmlsim[it].S2Rm,(double)matrixmlsim[it].SRm,data[0].n_loci))
                                                nunderalfa[11] += 1;
                                            if(matrixml[0].Sbiallsites == matrixmlsim[it].Sbiallsites)
                                                nequivalfa[0] += 1;
                                            if(varobs[0] == variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci))
                                                nequivalfa[1] += 1;
                                            if(matrixml[0].Snhapl == matrixmlsim[it].Snhapl)
                                                nequivalfa[2] += 1;
                                            if(varobs[1] == variances((double)matrixmlsim[it].S2nhapl,(double)matrixmlsim[it].Snhapl,data[0].n_loci))
                                                nequivalfa[3] += 1;
                                            if(matrixml[0].Snhaplsam == matrixmlsim[it].Snhaplsam)
                                                nequivalfa[4] += 1;
                                            if(varobs[2] == variances((double)matrixmlsim[it].S2nhaplsam,(double)matrixmlsim[it].Snhaplsam,data[0].n_loci))
                                                nequivalfa[5] += 1;
                                            if(matrixml[0].Shapldiv == matrixmlsim[it].Shapldiv)
                                                nequivalfa[6] += 1;
                                            if(varobs[3] == variances((double)matrixmlsim[it].S2hapldiv,(double)matrixmlsim[it].Shapldiv,data[0].n_loci))
                                                nequivalfa[7] += 1;
                                            if(matrixml[0].SRvpi == matrixmlsim[it].SRvpi)
                                                nequivalfa[8] += 1;
                                            if(varobs[4] == variances((double)matrixmlsim[it].S2Rvpi,(double)matrixmlsim[it].SRvpi,data[0].n_loci))
                                                nequivalfa[9] += 1;
                                            if(matrixml[0].SRm == matrixmlsim[it].SRm)
                                                nequivalfa[10] += 1;
                                            if(varobs[5] == variances((double)matrixmlsim[it].S2Rm,(double)matrixmlsim[it].SRm,data[0].n_loci))
                                                nequivalfa[11] += 1;
                                        }
                                        prob = (double)nunderalfa[0]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                        printf("\n    P(average biallelic sites = %g) = %g ",(double)matrixml[0].Sbiallsites/(double)data[0].n_loci,prob);
                                        if(file_output)
                                            fprintf(file_output,"\n    P(average biallelic sites = %g) = %g ",(double)matrixml[0].Sbiallsites/ (double)data[0].n_loci,prob);
                                        if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[1]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                        if(data[0].n_loci > 2) printf("\n    P(variance of biallelic sites = %g) = %g ",varobs[0],prob);
                                        if(file_output)
                                            if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of biallelic sites = %g) = %g ",varobs[0],prob);
                                        if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[2]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                        printf("\n    P(average number of haplotypes = %g) = %g ",(double)matrixml[0].Snhapl/(double)data[0].n_loci,prob);
                                        if(file_output)
                                           fprintf(file_output,"\n    P(average number of haplotypes = %g) = %g ",(double)matrixml[0].Snhapl/ (double)data[0].n_loci,prob);
                                        if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[3]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                        if(data[0].n_loci > 2) printf("\n    P(variance of number of haplotypes = %g) = %g ",varobs[1],prob);
                                        if(file_output)
                                            if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of number of haplotypes = %g) = %g ",varobs[1],prob);
                                        if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[4]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                        printf("\n    P(average haplotypes / sample size = %g) = %g ",(double)matrixml[0].Snhaplsam/(double)data[0].n_loci,prob);
                                        if(file_output)
                                            fprintf(file_output,"\n    P(average haplotypes / sample size = %g) = %g ",(double)matrixml[0].Snhaplsam/ (double)data[0].n_loci,prob);
                                        if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[5]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                        if(data[0].n_loci > 2) printf("\n    P(variance of haplotypes / sample size = %g) = %g ",varobs[2],prob);
                                        if(file_output)
                                            if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of haplotypes/sample size = %g) = %g ",varobs[2],prob);
                                        if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[6]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[6] >= (double)0.5) nequivalfa[6] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[6];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[6];
                                        printf("\n    P(average haplotype diversity = %g) = %g ",(double)matrixml[0].Shapldiv/(double)data[0].n_loci,prob);
                                        if(file_output)
                                            fprintf(file_output,"\n    P(average haplotype diversity = %g) = %g ",(double)matrixml[0].Shapldiv/ (double)data[0].n_loci,prob);
                                        if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[7]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[7] >= (double)0.5) nequivalfa[7] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[7];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[7];
                                        if(data[0].n_loci > 2) printf("\n    P(variance of haplotype diversity = %g) = %g ",varobs[3],prob);
                                        if(file_output)
                                            if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of haplotype diversity = %g) = %g ",varobs[3],prob);
                                        if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        /*
											prob = (double)nunderalfa[8]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[8] >= (double)0.5) nequivalfa[8] /= 2;
											if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[8];
											else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[8];
											printf("\n    P(average of C = %g) = %g ",(double)matrixml[0].SRvpi/ (double)data[0].n_loci,prob);
											if(file_output)
												fprintf(file_output,"\n    P(average of C = %g) = %g ",(double)matrixml[0].SRvpi/ (double)data[0].n_loci,prob);
											if(data[0].n_iter >=20) printfignif(prob,file_output);
										
											prob = (double)nunderalfa[9]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[9] >= (double)0.5) nequivalfa[9] /= 2;
											if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[9];
											else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[9];
											if(data[0].n_loci > 2) printf("\n    P(variance of C = %g) = %g ",varobs[4],prob);
											if(file_output)
												if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of C = %g) = %g ",varobs[4],prob);
											if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
										*/
										
                                        prob = (double)nunderalfa[10]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[10] >= (double)0.5) nequivalfa[10] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[10];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[10];
                                        printf("\n    P(average of Rm = %g) = %g ",(double)matrixml[0].SRm/ (double)data[0].n_loci,prob);
                                        if(file_output)
                                            fprintf(file_output,"\n    P(average of Rm = %g) = %g ",(double)matrixml[0].SRm/ (double)data[0].n_loci,prob);
                                        if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[11]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[11] >= (double)0.5) nequivalfa[11] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[11];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[11];
                                        if(data[0].n_loci > 2) printf("\n    P(variance of Rm = %g) = %g ",varobs[5],prob);
                                        if(file_output)
                                            if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of Rm = %g) = %g ",varobs[5],prob);
                                        if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                    }
                                    else {
                                        for(x=0;x<20;x++) nunderalfa[x] = 0;
                                        for(x=0;x<20;x++) nequivalfa[x] = 0;
                                        varobs[0] = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,data[0].n_loci);
                                        varobs[1] = variances((double)matrixml[0].S2fixed,(double)matrixml[0].Sfixed,data[0].n_loci);
                                        varobs[2] = variances((double)matrixml[0].S2ndivergence,(double)matrixml[0].Sndivergence,data[0].n_loci);
                                        varobs[3] = variances((double)matrixml[0].S2nhapl,(double)matrixml[0].Snhapl,data[0].n_loci);
                                        varobs[4] = variances((double)matrixml[0].S2nhaplsam,(double)matrixml[0].Snhaplsam,data[0].n_loci);
                                        varobs[5] = variances((double)matrixml[0].S2hapldiv,(double)matrixml[0].Shapldiv,data[0].n_loci);
                                        varobs[6] = variances((double)matrixml[0].S2Rvpi,(double)matrixml[0].SRvpi,data[0].n_loci);
                                        varobs[7] = variances((double)matrixml[0].S2Rm,(double)matrixml[0].SRm,data[0].n_loci);
										for(it=0;it<(long int)data[0].n_iter;it++) {
                                            if(matrixml[0].Sbiallsites < matrixmlsim[it].Sbiallsites)
                                                nunderalfa[0] += 1;
                                            if(varobs[0] < variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci))
                                                nunderalfa[1] += 1;
                                            if(matrixml[0].Sfixed < matrixmlsim[it].Sfixed)
                                                nunderalfa[2] += 1;
                                            if(varobs[1] < variances((double)matrixmlsim[it].S2fixed,(double)matrixmlsim[it].Sfixed,data[0].n_loci))
                                                nunderalfa[3] += 1;
                                            if(matrixml[0].Sndivergence < matrixmlsim[it].Sndivergence)
                                                nunderalfa[4] += 1;
                                            if(varobs[2] < variances((double)matrixmlsim[it].S2ndivergence,(double)matrixmlsim[it].Sndivergence,data[0].n_loci))
                                                nunderalfa[5] += 1;
                                            if(matrixml[0].Snhapl < matrixmlsim[it].Snhapl)
                                                nunderalfa[6] += 1;
                                            if(varobs[3] < variances((double)matrixmlsim[it].S2nhapl,(double)matrixmlsim[it].Snhapl,data[0].n_loci))
                                                nunderalfa[7] += 1;
                                             if(matrixml[0].Snhaplsam < matrixmlsim[it].Snhaplsam)
                                                nunderalfa[8] += 1;
                                            if(varobs[4] < variances((double)matrixmlsim[it].S2nhaplsam,(double)matrixmlsim[it].Snhaplsam,data[0].n_loci))
                                                nunderalfa[9] += 1;
                                             if(matrixml[0].Shapldiv < matrixmlsim[it].Shapldiv)
                                                nunderalfa[10] += 1;
                                            if(varobs[5] < variances((double)matrixmlsim[it].S2hapldiv,(double)matrixmlsim[it].Shapldiv,data[0].n_loci))
                                                nunderalfa[11] += 1;
                                             if(matrixml[0].SRvpi < matrixmlsim[it].SRvpi)
                                                nunderalfa[12] += 1;
                                            if(varobs[6] < variances((double)matrixmlsim[it].S2Rvpi,(double)matrixmlsim[it].SRvpi,data[0].n_loci))
                                                nunderalfa[13] += 1;
                                             if(matrixml[0].SRm < matrixmlsim[it].SRm)
                                                nunderalfa[14] += 1;
                                            if(varobs[7] < variances((double)matrixmlsim[it].S2Rm,(double)matrixmlsim[it].SRm,data[0].n_loci))
                                                nunderalfa[15] += 1;
 
											if(matrixml[0].Sbiallsites == matrixmlsim[it].Sbiallsites)
                                                nequivalfa[0] += 1;
                                            if(varobs[0] == variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci))
                                                nequivalfa[1] += 1;
                                            if(matrixml[0].Sfixed == matrixmlsim[it].Sfixed)
                                                nequivalfa[2] += 1;
                                            if(varobs[1] == variances((double)matrixmlsim[it].S2fixed,(double)matrixmlsim[it].Sfixed,data[0].n_loci))
                                                nequivalfa[3] += 1;
                                            if(matrixml[0].Sndivergence == matrixmlsim[it].Sndivergence)
                                                nequivalfa[4] += 1;
                                            if(varobs[2] == variances((double)matrixmlsim[it].S2ndivergence,(double)matrixmlsim[it].Sndivergence,data[0].n_loci))
                                                nequivalfa[5] += 1;
                                            if(matrixml[0].Snhapl == matrixmlsim[it].Snhapl)
                                                nequivalfa[6] += 1;
                                            if(varobs[3] == variances((double)matrixmlsim[it].S2nhapl,(double)matrixmlsim[it].Snhapl,data[0].n_loci))
                                                nequivalfa[7] += 1;
                                            if(matrixml[0].Snhaplsam == matrixmlsim[it].Snhaplsam)
                                                nequivalfa[8] += 1;
                                            if(varobs[4] == variances((double)matrixmlsim[it].S2nhaplsam,(double)matrixmlsim[it].Snhaplsam,data[0].n_loci))
                                                nequivalfa[9] += 1;
                                            if(matrixml[0].Shapldiv == matrixmlsim[it].Shapldiv)
                                                nequivalfa[10] += 1;
                                            if(varobs[5] == variances((double)matrixmlsim[it].S2hapldiv,(double)matrixmlsim[it].Shapldiv,data[0].n_loci))
                                                nequivalfa[11] += 1;
                                            if(matrixml[0].SRvpi == matrixmlsim[it].SRvpi)
                                                nequivalfa[12] += 1;
                                            if(varobs[6] == variances((double)matrixmlsim[it].S2Rvpi,(double)matrixmlsim[it].SRvpi,data[0].n_loci))
                                                nequivalfa[13] += 1;
                                            if(matrixml[0].SRm == matrixmlsim[it].SRm)
                                                nequivalfa[14] += 1;
                                            if(varobs[7] == variances((double)matrixmlsim[it].S2Rm,(double)matrixmlsim[it].SRm,data[0].n_loci))
                                                nequivalfa[15] += 1;
                                        }
                                        prob = (double)nunderalfa[0]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                        printf("\n    P(average biallelic sites = %g) = %g ",(double)matrixml[0].Sbiallsites/ (double)data[0].n_loci,prob);
                                        if(file_output)
                                            fprintf(file_output,"\n    P(average biallelic sites = %g) = %g ",(double)matrixml[0].Sbiallsites/ (double)data[0].n_loci,prob);
                                        if(data[0].n_iter >=20) printfignif(prob,file_output);
                                        prob = (double)nunderalfa[1]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                        if(data[0].n_loci > 2) printf("\n    P(variance of biallelic sites = %g) = %g ",varobs[0],prob);
                                        if(file_output)
                                            if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of biallelic sites = %g) = %g ",varobs[0],prob);
                                        if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[2]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                        if(maxnoutgroups == 1) printf("\n    P(average of fixed sites = %g) = %g ",(double)matrixml[0].Sfixed/ (double)data[0].n_loci,prob);
                                        if(file_output)
                                            if(maxnoutgroups == 1) fprintf(file_output,"\n    P(average of fixed sites = %g) = %g ",(double)matrixml[0].Sfixed/ (double)data[0].n_loci,prob);
                                        if(maxnoutgroups == 1) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                        prob = (double)nunderalfa[3]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                        if(data[0].n_loci > 2) if(maxnoutgroups == 1) printf("\n    P(variance of fixed sites = %g) = %g ",varobs[1],prob);
                                        if(file_output)
                                            if(data[0].n_loci > 2) if(maxnoutgroups == 1) fprintf(file_output,"\n    P(variance of fixed sites = %g) = %g ",varobs[1],prob);
                                        if(data[0].n_loci > 2) if(maxnoutgroups == 1) if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[4]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                        printf("\n    P(average of total divergence = %g) = %g ",(double)matrixml[0].Sndivergence/ (double)data[0].n_loci,prob);
                                        if(file_output)
                                            fprintf(file_output,"\n    P(average of total divergence = %g) = %g ",(double)matrixml[0].Sndivergence/ (double)data[0].n_loci,prob);
                                        if(data[0].n_iter >=20) printfignif(prob,file_output);
                                        prob = (double)nunderalfa[5]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                        if(data[0].n_loci > 2) printf("\n    P(variance of total divergence = %g) = %g ",varobs[2],prob);
                                        if(file_output)
                                            if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of total divergence = %g) = %g ",varobs[2],prob);
                                        if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[6]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[6] >= (double)0.5) nequivalfa[6] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[6];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[6];
                                        printf("\n    P(average of number of haplotypes = %g) = %g ",(double)matrixml[0].Snhapl/ (double)data[0].n_loci,prob);
                                        if(file_output)
                                            fprintf(file_output,"\n    P(average of number of haplotypes = %g) = %g ",(double)matrixml[0].Snhapl/ (double)data[0].n_loci,prob);
                                        if(data[0].n_iter >=20) printfignif(prob,file_output);
                                        prob = (double)nunderalfa[7]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[7] >= (double)0.5) nequivalfa[7] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[7];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[7];
                                        if(data[0].n_loci > 2) printf("\n    P(variance of number of haplotypes = %g) = %g ",varobs[3],prob);
                                        if(file_output)
                                            if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of number of haplotypes = %g) = %g ",varobs[3],prob);
                                        if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[8]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[8] >= (double)0.5) nequivalfa[8] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[8];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[8];
                                        printf("\n    P(average of number of haplotypes / sample size = %g) = %g ",(double)matrixml[0].Snhaplsam/ (double)data[0].n_loci,prob);
                                        if(file_output)
                                            fprintf(file_output,"\n    P(average of number of haplotypes / sample size = %g) = %g ",(double)matrixml[0].Snhaplsam/ (double)data[0].n_loci,prob);
                                        if(data[0].n_iter >=20) printfignif(prob,file_output);
                                        prob = (double)nunderalfa[9]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[9] >= (double)0.5) nequivalfa[9] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[9];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[9];
                                        if(data[0].n_loci > 2) printf("\n    P(variance of number of haplotypes / sample size = %g) = %g ",varobs[4],prob);
                                        if(file_output)
                                            if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of number of haplotypes / sample size = %g) = %g ",varobs[4],prob);
                                        if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        prob = (double)nunderalfa[10]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[10] >= (double)0.5) nequivalfa[10] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[10];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[10];
                                        printf("\n    P(average of haplotype diversity = %g) = %g ",(double)matrixml[0].Shapldiv / (double)data[0].n_loci,prob);
                                        if(file_output)
                                            fprintf(file_output,"\n    P(average of haplotype diversity = %g) = %g ",(double)matrixml[0].Shapldiv/ (double)data[0].n_loci,prob);
                                        if(data[0].n_iter >=20) printfignif(prob,file_output);
                                        prob = (double)nunderalfa[11]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[11] >= (double)0.5) nequivalfa[11] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[11];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[11];
                                        if(data[0].n_loci > 2) printf("\n    P(variance of haplotype diversity = %g) = %g ",varobs[5],prob);
                                        if(file_output)
                                            if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of haplotype diversity = %g) = %g ",varobs[5],prob);
                                        if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
										
                                        /*
											prob = (double)nunderalfa[12]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[12] >= (double)0.5) nequivalfa[12] /= 2;
											if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[12];
											else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[12];
											printf("\n    P(average of C = %g) = %g ",(double)matrixml[0].SRvpi / (double)data[0].n_loci,prob);
											if(file_output)
												fprintf(file_output,"\n    P(average of C = %g) = %g ",(double)matrixml[0].SRvpi/ (double)data[0].n_loci,prob);
											if(data[0].n_iter >=20) printfignif(prob,file_output);
											prob = (double)nunderalfa[13]/(double)data[0].n_iter;
											if((double)1./((double)data[0].n_iter) *(double)nequivalfa[13] >= (double)0.5) nequivalfa[13] /= 2;
											if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[13];
											else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[13];
											if(data[0].n_loci > 2) printf("\n    P(variance of C = %g) = %g ",varobs[6],prob);
											if(file_output)
												if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of C = %g) = %g ",varobs[6],prob);
											if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
										*/
										
                                        prob = (double)nunderalfa[14]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[14] >= (double)0.5) nequivalfa[14] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[14];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[14];
                                        printf("\n    P(average of Rm = %g) = %g ",(double)matrixml[0].SRm / (double)data[0].n_loci,prob);
                                        if(file_output)
                                            fprintf(file_output,"\n    P(average of Rm = %g) = %g ",(double)matrixml[0].SRm/ (double)data[0].n_loci,prob);
                                        if(data[0].n_iter >=20) printfignif(prob,file_output);
                                        prob = (double)nunderalfa[15]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[15] >= (double)0.5) nequivalfa[15] /= 2;
                                        if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[15];
                                        else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[15];
                                        if(data[0].n_loci > 2) printf("\n    P(variance of Rm = %g) = %g ",varobs[7],prob);
                                        if(file_output)
                                            if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of Rm = %g) = %g ",varobs[7],prob);
                                        if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    }
                                }
                                else {/*not outgroup, observed_data*/
                                    /*average and variance comparisons.*/
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ");
                                    printf("\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%.");
                                    printf("\n Multiple hits not included in the analysis.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ",file_output);
                                        fprintf(file_output,"\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%.");
                                        fputs("\n Multiple hits not included in the analysis.\n",file_output);
                                    }
									for(x=0;x<20;x++) nunderalfa[x] = 0;
									for(x=0;x<20;x++) nequivalfa[x] = 0;
                                    varobs[0] = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,data[0].n_loci);
                                    varobs[1] = variances((double)matrixml[0].S2nhapl,(double)matrixml[0].Snhapl,data[0].n_loci);                                                                          
                                    varobs[2] = variances((double)matrixml[0].S2nhaplsam,(double)matrixml[0].Snhaplsam,data[0].n_loci);                                                                          
                                    varobs[3] = variances((double)matrixml[0].S2hapldiv,(double)matrixml[0].Shapldiv,data[0].n_loci);                                                                          
                                    varobs[4] = variances((double)matrixml[0].S2Rvpi,(double)matrixml[0].SRvpi,data[0].n_loci);                                                                          
                                    varobs[5] = variances((double)matrixml[0].S2Rm,(double)matrixml[0].SRm,data[0].n_loci);                                                                          
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        if(matrixml[0].Sbiallsites < matrixmlsim[it].Sbiallsites)
                                            nunderalfa[0] += 1;
                                        if(varobs[0] < variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci))
                                            nunderalfa[1] += 1;
                                        if(matrixml[0].Snhapl < matrixmlsim[it].Snhapl)
                                            nunderalfa[2] += 1;
                                        if(varobs[1] < variances((double)matrixmlsim[it].S2nhapl,(double)matrixmlsim[it].Snhapl,data[0].n_loci))
                                            nunderalfa[3] += 1;
                                        if(matrixml[0].Snhaplsam < matrixmlsim[it].Snhaplsam)
                                            nunderalfa[4] += 1;
                                        if(varobs[2] < variances((double)matrixmlsim[it].S2nhaplsam,(double)matrixmlsim[it].Snhaplsam,data[0].n_loci))
                                            nunderalfa[5] += 1;
                                        if(matrixml[0].Shapldiv < matrixmlsim[it].Shapldiv)
                                            nunderalfa[6] += 1;
                                        if(varobs[3] < variances((double)matrixmlsim[it].S2hapldiv,(double)matrixmlsim[it].Shapldiv,data[0].n_loci))
                                            nunderalfa[7] += 1;
                                        if(matrixml[0].SRvpi < matrixmlsim[it].SRvpi)
                                            nunderalfa[8] += 1;
                                        if(varobs[4] < variances((double)matrixmlsim[it].S2Rvpi,(double)matrixmlsim[it].SRvpi,data[0].n_loci))
                                            nunderalfa[9] += 1;
                                        if(matrixml[0].SRm < matrixmlsim[it].SRm)
                                            nunderalfa[10] += 1;
                                        if(varobs[5] < variances((double)matrixmlsim[it].S2Rm,(double)matrixmlsim[it].SRm,data[0].n_loci))
                                            nunderalfa[11] += 1;
											
                                        if(matrixml[0].Sbiallsites == matrixmlsim[it].Sbiallsites)
                                            nequivalfa[0] += 1;
                                        if(varobs[0] == variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci))
                                            nequivalfa[1] += 1;
                                        if(matrixml[0].Snhapl == matrixmlsim[it].Snhapl)
                                            nequivalfa[2] += 1;
                                        if(varobs[1] == variances((double)matrixmlsim[it].S2nhapl,(double)matrixmlsim[it].Snhapl,data[0].n_loci))
                                            nequivalfa[3] += 1;
                                        if(matrixml[0].Snhaplsam == matrixmlsim[it].Snhaplsam)
                                            nequivalfa[4] += 1;
                                        if(varobs[2] == variances((double)matrixmlsim[it].S2nhaplsam,(double)matrixmlsim[it].Snhaplsam,data[0].n_loci))
                                            nequivalfa[5] += 1;
                                        if(matrixml[0].Shapldiv == matrixmlsim[it].Shapldiv)
                                            nequivalfa[6] += 1;
                                        if(varobs[3] == variances((double)matrixmlsim[it].S2hapldiv,(double)matrixmlsim[it].Shapldiv,data[0].n_loci))
                                            nequivalfa[7] += 1;
                                        if(matrixml[0].SRvpi == matrixmlsim[it].SRvpi)
                                            nequivalfa[8] += 1;
                                        if(varobs[4] == variances((double)matrixmlsim[it].S2Rvpi,(double)matrixmlsim[it].SRvpi,data[0].n_loci))
                                            nequivalfa[9] += 1;
                                        if(matrixml[0].SRm == matrixmlsim[it].SRm)
                                            nequivalfa[10] += 1;
                                        if(varobs[5] == variances((double)matrixmlsim[it].S2Rm,(double)matrixmlsim[it].SRm,data[0].n_loci))
                                            nequivalfa[11] += 1;
                                    }
                                    prob = (double)nunderalfa[0]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                    printf("\n    P(average biallelic sites = %g) = %g ",(double)matrixml[0].Sbiallsites/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average biallelic sites = %g) = %g ",(double)matrixml[0].Sbiallsites/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    prob = (double)nunderalfa[1]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of biallelic sites = %g) = %g ",varobs[0],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of biallelic sites = %g) = %g ",varobs[0],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
									
                                    prob = (double)nunderalfa[2]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                    printf("\n    P(average number of haplotypes = %g) = %g ",(double)matrixml[0].Snhapl/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average number of haplotypes = %g) = %g ",(double)matrixml[0].Snhapl/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    prob = (double)nunderalfa[3]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of number of haplotypes = %g) = %g ",varobs[1],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of number of haplotypes = %g) = %g ",varobs[1],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
									
                                    prob = (double)nunderalfa[4]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                    printf("\n    P(average number of haplotypes / sample size = %g) = %g ",(double)matrixml[0].Snhaplsam/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average number of haplotypes / sample size = %g) = %g ",(double)matrixml[0].Snhaplsam/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    prob = (double)nunderalfa[5]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of number of haplotypes / sample size = %g) = %g ",varobs[2],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of number of haplotypes / sample size = %g) = %g ",varobs[2],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
									
                                    prob = (double)nunderalfa[6]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[6] >= (double)0.5) nequivalfa[6] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[6];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[6];
                                    printf("\n    P(average of haplotype diversity = %g) = %g ",(double)matrixml[0].Shapldiv/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of haplotype diversity = %g) = %g ",(double)matrixml[0].Shapldiv/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    prob = (double)nunderalfa[7]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[7] >= (double)0.5) nequivalfa[7] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[7];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[7];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of haplotype diversity = %g) = %g ",varobs[3],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of haplotype diversity = %g) = %g ",varobs[3],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
									
									/*
										prob = (double)nunderalfa[8]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[8] >= (double)0.5) nequivalfa[8] /= 2;
										if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[8];
										else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[8];
										printf("\n    P(average of C = %g) = %g ",(double)matrixml[0].SRvpi / (double)data[0].n_loci,prob);
										if(file_output)
											fprintf(file_output,"\n    P(average of C = %g) = %g ",(double)matrixml[0].SRvpi/ (double)data[0].n_loci,prob);
										if(data[0].n_iter >=20) printfignif(prob,file_output);
										prob = (double)nunderalfa[9]/(double)data[0].n_iter;
										if((double)1./((double)data[0].n_iter) *(double)nequivalfa[9] >= (double)0.5) nequivalfa[9] /= 2;
										if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[9];
										else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[9];
										if(data[0].n_loci > 2) printf("\n    P(variance of C = %g) = %g ",varobs[4],prob);
										if(file_output)
											if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of C = %g) = %g ",varobs[4],prob);
										if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
									*/
									
									prob = (double)nunderalfa[10]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[10] >= (double)0.5) nequivalfa[10] /= 2;
									if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[10];
									else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[10];
									printf("\n    P(average of Rm = %g) = %g ",(double)matrixml[0].SRm / (double)data[0].n_loci,prob);
									if(file_output)
										fprintf(file_output,"\n    P(average of Rm = %g) = %g ",(double)matrixml[0].SRm/ (double)data[0].n_loci,prob);
									if(data[0].n_iter >=20) printfignif(prob,file_output);
									prob = (double)nunderalfa[11]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[11] >= (double)0.5) nequivalfa[11] /= 2;
									if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[11];
									else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[11];
									if(data[0].n_loci > 2) printf("\n    P(variance of Rm = %g) = %g ",varobs[5],prob);
									if(file_output)
										if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of Rm = %g) = %g ",varobs[5],prob);
									if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                }
                            }
							#if ALLRESULTS == 0
                            else {/*give confident intervals with simulations. NOT observed_data.*/
							#endif
                                if((sortvector = (double *) malloc(data[0].n_iter * sizeof(double))) == 0) {
                                    printf("Probabilities can not be calculated, sorry.\n");
                                    if(file_output) fputs("Probabilities can not be calculated, sorry.\n",file_output);
                                    break;
                                }
                                if(*outgroup) {
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n");
                                    printf("statistic\tmedian\tavg\tvar\t");                                    
                                    if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                    if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                    if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                    if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                    if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                    if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                    if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                    printf("\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n",file_output);
                                        fputs("statistic\tmedian\tavg\tvar\t",file_output);                                    
                                        if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                        if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                        if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                        if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                        if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                        if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                        if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                        fprintf(file_output,"\n");
                                    }

                                    if(data[0].time_spec <= 0.) {
                                        /*Sbiallelic sites. */
                                        /*average*/
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if( matrixmlsim[it].Sbiallsites != (double) -10000) {
                                                sortvector[it2] = (double)matrixmlsim[it].Sbiallsites/(double)data[0].n_loci;
                                                it2 += 1;
                                                avg += (double)matrixmlsim[it].Sbiallsites/(double)data[0].n_loci;
                                                var += (double)matrixmlsim[it].Sbiallsites/(double)data[0].n_loci * (double)matrixmlsim[it].Sbiallsites/(double)data[0].n_loci;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("avg_bialsites\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"avg_bialsites\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                        /*variance*/
                                        if(data[0].n_loci > 2) {
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if((vars=variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci)) != (double) -10000) {
                                                    sortvector[it2] = vars;
                                                    it2 += 1;
                                                    avg += vars;
                                                    var += vars * vars;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("var_bialsites\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"var_bialsites\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                        						}
                                            print_percentages(data,sortvector,it2,file_output);
                                        }
                                    }
                                    else {
                                        /*Sbiallelic sites. */
                                        /*average*/
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if( matrixmlsim[it].Sbiallsites != (double) -10000) {
                                                sortvector[it2] = (double)matrixmlsim[it].Sbiallsites/(double)data[0].n_loci;
                                                it2 += 1;
                                                avg += (double)matrixmlsim[it].Sbiallsites/(double)data[0].n_loci;
                                                var += (double)matrixmlsim[it].Sbiallsites/(double)data[0].n_loci * (double)matrixmlsim[it].Sbiallsites/(double)data[0].n_loci;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("avg_bialsites\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"avg_bialsites\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                        /*variance*/
                                        if(data[0].n_loci > 2) {
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if((vars=variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci)) != (double) -10000) {
                                                    sortvector[it2] = vars;
                                                    it2 += 1;
                                                    avg += vars;
                                                    var += vars * vars;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("var_bialsites\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"var_bialsites\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                        }
                                            print_percentages(data,sortvector,it2,file_output);
                                        }
                                        /*fixed*/
                                        /*average*/
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if( matrixmlsim[it].Sfixed != (double) -10000) {
                                                sortvector[it2] = (double)matrixmlsim[it].Sfixed/(double)data[0].n_loci;
                                                it2 += 1;
                                                avg += (double)matrixmlsim[it].Sfixed/(double)data[0].n_loci;
                                                var += (double)matrixmlsim[it].Sfixed/(double)data[0].n_loci * (double)matrixmlsim[it].Sfixed/(double)data[0].n_loci;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("avg_fixed\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"avg_fixed\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                        /*variance*/
                                        if(data[0].n_loci > 2) {
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if((vars=variances((double)matrixmlsim[it].S2fixed,(double)matrixmlsim[it].Sfixed,data[0].n_loci)) != (double) -10000) {
                                                    sortvector[it2] = vars;
                                                    it2 += 1;
                                                    avg += vars;
                                                    var += vars * vars;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("var_fixed\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"var_fixed\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                        
											}
                                            print_percentages(data,sortvector,it2,file_output);
                                        }
                                        /*divergence*/
                                        /*average*/
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if( matrixmlsim[it].Sndivergence != (double) -10000) {
                                                sortvector[it2] = (double)matrixmlsim[it].Sndivergence/(double)data[0].n_loci;
                                                it2 += 1;
                                                avg += (double)matrixmlsim[it].Sndivergence/(double)data[0].n_loci;
                                                var += (double)matrixmlsim[it].Sndivergence/(double)data[0].n_loci * (double)matrixmlsim[it].Sndivergence/(double)data[0].n_loci;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("avg_Diverg\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"avg_Diverg\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                        /*variance*/
                                        if(data[0].n_loci > 2) {
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if((vars=variances(matrixmlsim[it].S2ndivergence,matrixmlsim[it].Sndivergence,data[0].n_loci)) != (double) -10000) {
                                                    sortvector[it2] = vars;
                                                    it2 += 1;
                                                    avg += vars;
                                                    var += vars * vars;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("var_Diverg\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"var_Diverg\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                        }
                                            print_percentages(data,sortvector,it2,file_output);
                                        }
                                        /*nhaplotypes*/
                                        /*average*/
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if( matrixmlsim[it].Snhapl != (double) -10000) {
                                                sortvector[it2] = (double)matrixmlsim[it].Snhapl/(double)data[0].n_loci;
                                                it2 += 1;
                                                avg += (double)matrixmlsim[it].Snhapl/(double)data[0].n_loci;
                                                var += (double)matrixmlsim[it].Snhapl/(double)data[0].n_loci * (double)matrixmlsim[it].Snhapl/(double)data[0].n_loci;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("avg_nhapl\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"avg_nhapl\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                        /*variance*/
                                        if(data[0].n_loci > 2) {
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if((vars=variances(matrixmlsim[it].S2nhapl,matrixmlsim[it].Snhapl,data[0].n_loci)) != (double) -10000) {
                                                    sortvector[it2] = vars;
                                                    it2 += 1;
                                                    avg += vars;
                                                    var += vars * vars;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("var_nhapl\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"var_nhapl\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                        }
                                            print_percentages(data,sortvector,it2,file_output);
                                        }
                                         /*nhaplotypes/samsize*/
                                        /*average*/
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if( matrixmlsim[it].Snhaplsam != (double) -10000) {
                                                sortvector[it2] = (double)matrixmlsim[it].Snhaplsam/(double)data[0].n_loci;
                                                it2 += 1;
                                                avg += (double)matrixmlsim[it].Snhaplsam/(double)data[0].n_loci;
                                                var += (double)matrixmlsim[it].Snhaplsam/(double)data[0].n_loci * (double)matrixmlsim[it].Snhaplsam/(double)data[0].n_loci;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("avg_nhaplsam\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"avg_nhaplsam\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                        /*variance*/
                                        if(data[0].n_loci > 2) {
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if((vars=variances(matrixmlsim[it].S2nhaplsam,matrixmlsim[it].Snhaplsam,data[0].n_loci)) != (double) -10000) {
                                                    sortvector[it2] = vars;
                                                    it2 += 1;
                                                    avg += vars;
                                                    var += vars * vars;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("var_nhaplsam\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"var_nhaplsam\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                        }
                                            print_percentages(data,sortvector,it2,file_output);
                                        }
                                        /*nhaplotype diversity*/
                                        /*average*/
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if( matrixmlsim[it].Shapldiv != (double) -10000) {
                                                sortvector[it2] = (double)matrixmlsim[it].Shapldiv/(double)data[0].n_loci;
                                                it2 += 1;
                                                avg += (double)matrixmlsim[it].Shapldiv/(double)data[0].n_loci;
                                                var += (double)matrixmlsim[it].Shapldiv/(double)data[0].n_loci * (double)matrixmlsim[it].Shapldiv/(double)data[0].n_loci;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("avg_hapldiv\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"avg_hapldiv\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                        /*variance*/
                                        if(data[0].n_loci > 2) {
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if((vars=variances(matrixmlsim[it].S2hapldiv,matrixmlsim[it].Shapldiv,data[0].n_loci)) != (double) -10000) {
                                                    sortvector[it2] = vars;
                                                    it2 += 1;
                                                    avg += vars;
                                                    var += vars * vars;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("var_hapldiv\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"var_hapldiv\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                        }
                                            print_percentages(data,sortvector,it2,file_output);
                                        }
                                        /*C*/
                                        /*average*//*
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if( matrixmlsim[it].SRvpi != (double) -10000) {
                                                sortvector[it2] = (double)matrixmlsim[it].SRvpi/(double)data[0].n_loci;
                                                it2 += 1;
                                                avg += (double)matrixmlsim[it].SRvpi/(double)data[0].n_loci;
                                                var += (double)matrixmlsim[it].SRvpi/(double)data[0].n_loci * (double)matrixmlsim[it].SRvpi/(double)data[0].n_loci;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("avg_C\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"avg_C\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                        *//*variance*//*
                                        if(data[0].n_loci > 2) {
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if((vars=variances(matrixmlsim[it].S2Rvpi,matrixmlsim[it].SRvpi,data[0].n_loci)) != (double) -10000) {
                                                    sortvector[it2] = vars;
                                                    it2 += 1;
                                                    avg += vars;
                                                    var += vars * vars;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("var_C\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"var_C\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                        }
                                            print_percentages(data,sortvector,it2,file_output);
                                        }*/
                                        /*Rm*/
                                        /*average*/
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if( matrixmlsim[it].SRm != (double) -10000) {
                                                sortvector[it2] = (double)matrixmlsim[it].SRm/(double)data[0].n_loci;
                                                it2 += 1;
                                                avg += (double)matrixmlsim[it].SRm/(double)data[0].n_loci;
                                                var += (double)matrixmlsim[it].SRm/(double)data[0].n_loci * (double)matrixmlsim[it].SRm/(double)data[0].n_loci;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("avg_Rm\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"avg_Rm\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                        }
                                        print_percentages(data,sortvector,it2,file_output);
                                        /*variance*/
                                        if(data[0].n_loci > 2) {
                                            avg = var = (double)0.0;
                                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                                if((vars=variances((double)matrixmlsim[it].S2Rm,(double)matrixmlsim[it].SRm,data[0].n_loci)) != (double) -10000) {
                                                    sortvector[it2] = vars;
                                                    it2 += 1;
                                                    avg += vars;
                                                    var += vars * vars;
                                                }
                                            }
                                            var = varianceit(var,avg,it2);
                                            if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                            printf("var_Rm\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(file_output)
                                                fprintf(file_output,"var_Rm\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                            if(var != (double)-10000) {
                                                printf("%g\t",var);
                                                if(file_output) fprintf(file_output,"%g\t",var);
                                            }
                                            else {
                                                printf("na\t");
                                                if(file_output) fprintf(file_output,"na\t");                                        }
                                            print_percentages(data,sortvector,it2,file_output);
                                        }
                                   }
                                }
                                else {
                                    /*Do a NEW vector with all the values, SORT.*/
                                    /* Take the median (50%), average and variance, and 10%, 5%, 2.5%, 1%, 0.5%, 0.1% and 0.05% for 1 tails.*/
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n");
                                    printf("statistic\tmedian\tavg\tvar\t");                                    
                                    if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                    if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                    if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                    if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                    if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                    if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                    if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                    printf("\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n",file_output);
                                        fputs("statistic\tmedian\tavg\tvar\t",file_output);                                    
                                        if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                        if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                        if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                        if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                        if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                        if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                        if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\n",file_output);
                                        fprintf(file_output,"\n");
                                    }
                                    /*Sbiallelic sites. */
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Sbiallsites != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Sbiallsites/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Sbiallsites/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Sbiallsites/(double)data[0].n_loci * (double)matrixmlsim[it].Sbiallsites/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_bialsites\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_bialsites\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");
                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_bialsites\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_bialsites\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*nhapl*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Snhapl != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Snhapl/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Snhapl/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Snhapl/(double)data[0].n_loci * (double)matrixmlsim[it].Snhapl/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_nhapl\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_nhapl\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2nhapl,matrixmlsim[it].Snhapl,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_nhapl\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_nhapl\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*nhapl/sam*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Snhaplsam != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Snhaplsam/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Snhaplsam/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Snhaplsam/(double)data[0].n_loci * (double)matrixmlsim[it].Snhaplsam/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_nhaplsam\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_nhaplsam\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2nhaplsam,matrixmlsim[it].Snhaplsam,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_nhaplsam\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_nhaplsam\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*nhapl div*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Shapldiv != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Shapldiv/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Shapldiv/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Shapldiv/(double)data[0].n_loci * (double)matrixmlsim[it].Shapldiv/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_hapldiv\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_hapldiv\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2hapldiv,matrixmlsim[it].Shapldiv,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_hapldiv\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_hapldiv\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
									/*C*/
									/*average*//*
									avg = var = (double)0.0;
									for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
										if( matrixmlsim[it].SRvpi != (double) -10000) {
											sortvector[it2] = (double)matrixmlsim[it].SRvpi/(double)data[0].n_loci;
											it2 += 1;
											avg += (double)matrixmlsim[it].SRvpi/(double)data[0].n_loci;
											var += (double)matrixmlsim[it].SRvpi/(double)data[0].n_loci * (double)matrixmlsim[it].SRvpi/(double)data[0].n_loci;
										}
									}
									var = varianceit(var,avg,it2);
									if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
									printf("avg_C\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
									if(file_output)
										fprintf(file_output,"avg_C\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
									if(var != (double)-10000) {
										printf("%g\t",var);
										if(file_output) fprintf(file_output,"%g\t",var);
									}
									else {
										printf("na\t");
										if(file_output) fprintf(file_output,"na\t");                                        }
									print_percentages(data,sortvector,it2,file_output);
									*//*variance*//*
									if(data[0].n_loci > 2) {
										avg = var = (double)0.0;
										for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
											if((vars=variances(matrixmlsim[it].S2Rvpi,matrixmlsim[it].SRvpi,data[0].n_loci)) != (double) -10000) {
												sortvector[it2] = vars;
												it2 += 1;
												avg += vars;
												var += vars * vars;
											}
										}
										var = varianceit(var,avg,it2);
										if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
										printf("var_C\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
										if(file_output)
											fprintf(file_output,"var_C\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
										if(var != (double)-10000) {
											printf("%g\t",var);
											if(file_output) fprintf(file_output,"%g\t",var);
										}
										else {
											printf("na\t");
											if(file_output) fprintf(file_output,"na\t");                                        }
										print_percentages(data,sortvector,it2,file_output);
									}*/
									/*Rm*/
									/*average*/
									avg = var = (double)0.0;
									for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
										if( matrixmlsim[it].SRm != (double) -10000) {
											sortvector[it2] = (double)matrixmlsim[it].SRm/(double)data[0].n_loci;
											it2 += 1;
											avg += (double)matrixmlsim[it].SRm/(double)data[0].n_loci;
											var += (double)matrixmlsim[it].SRm/(double)data[0].n_loci * (double)matrixmlsim[it].SRm/(double)data[0].n_loci;
										}
									}
									var = varianceit(var,avg,it2);
									if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
									printf("avg_Rm\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
									if(file_output)
										fprintf(file_output,"avg_Rm\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
									if(var != (double)-10000) {
										printf("%g\t",var);
										if(file_output) fprintf(file_output,"%g\t",var);
									}
									else {
										printf("na\t");
										if(file_output) fprintf(file_output,"na\t");                                        }
									print_percentages(data,sortvector,it2,file_output);
									/*variance*/
									if(data[0].n_loci > 2) {
										avg = var = (double)0.0;
										for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
											if((vars=variances((double)matrixmlsim[it].S2Rm,(double)matrixmlsim[it].SRm,data[0].n_loci)) != (double) -10000) {
												sortvector[it2] = vars;
												it2 += 1;
												avg += vars;
												var += vars * vars;
											}
										}
										var = varianceit(var,avg,it2);
										if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
										printf("var_Rm\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
										if(file_output)
											fprintf(file_output,"var_Rm\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
										if(var != (double)-10000) {
											printf("%g\t",var);
											if(file_output) fprintf(file_output,"%g\t",var);
										}
										else {
											printf("na\t");
											if(file_output) fprintf(file_output,"na\t");                                        }
										print_percentages(data,sortvector,it2,file_output);
									}
                                }
                                free(sortvector);
							#if ALLRESULTS == 0
                            }
							#endif
                            break;
                        case '1': /* 1 - Display linkage related statistics*/
                            if(*observed_data && dataobsequalsim) {/*compare simulations with observed data.*/
                                /*average and variance comparisons.*/
								printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                printf("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ");
                                printf("\n E.g., Probabilities lower than 0.025 and higher than 0.975 would be significant at 5%%. ");
                                printf("\n Multiple hits not included in the analysis.\n");
                                if(file_output) {
									if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    fputs("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ",file_output);
                                    fprintf(file_output,"\n E.g., Probabilities lower than 0.025 and higher than 0.975 would be significant at 5%%. ");
                                    fputs("\n Multiple hits not included in the analysis.\n",file_output);
                                }
                                for(x=0;x<20;x++) nunderalfa[x] = 0;
                                for(x=0;x<20;x++) nequivalfa[x] = 0;
                                varobs[0] = variances((double)matrixml[0].S2b,(double)matrixml[0].Sb,data[0].n_loci);
                                varobs[1] = variances((double)matrixml[0].S2q,(double)matrixml[0].Sq,data[0].n_loci);
                                varobs[2] = variances((double)matrixml[0].S2za,(double)matrixml[0].Sza,data[0].n_loci);
                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                    if(matrixml[0].Sb < matrixmlsim[it].Sb)
                                        nunderalfa[0] += 1;
                                    if(varobs[0] < variances((double)matrixmlsim[it].S2b,(double)matrixmlsim[it].Sb,data[0].n_loci))
                                        nunderalfa[1] += 1;
                                    if(matrixml[0].Sq < matrixmlsim[it].Sq)
                                        nunderalfa[2] += 1;
                                    if(varobs[1] < variances((double)matrixmlsim[it].S2q,(double)matrixmlsim[it].Sq,data[0].n_loci))
                                        nunderalfa[3] += 1;
                                    if(matrixml[0].Sza < matrixmlsim[it].Sza)
                                        nunderalfa[4] += 1;
                                    if(varobs[2] < variances((double)matrixmlsim[it].S2za,(double)matrixmlsim[it].Sza,data[0].n_loci))
                                        nunderalfa[5] += 1;
                                    if(matrixml[0].Sb == matrixmlsim[it].Sb)
                                        nequivalfa[0] += 1;
                                    if(varobs[0] == variances((double)matrixmlsim[it].S2b,(double)matrixmlsim[it].Sb,data[0].n_loci))
                                        nequivalfa[1] += 1;
                                    if(matrixml[0].Sq == matrixmlsim[it].Sq)
                                        nequivalfa[2] += 1;
                                    if(varobs[1] == variances((double)matrixmlsim[it].S2q,(double)matrixmlsim[it].Sq,data[0].n_loci))
                                        nequivalfa[3] += 1;
                                    if(matrixml[0].Sza == matrixmlsim[it].Sza)
                                        nequivalfa[4] += 1;
                                    if(varobs[2] == variances((double)matrixmlsim[it].S2za,(double)matrixmlsim[it].Sza,data[0].n_loci))
                                        nequivalfa[5] += 1;
                                }
                                prob = (double)nunderalfa[0]/(double)data[0].n_iter;
								if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                printf("\n    P(average of b = %g) = %g ",(double)matrixml[0].Sb/ (double)data[0].n_loci,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of b = %g) = %g ",(double)matrixml[0].Sb/ (double)data[0].n_loci,prob);
                                if(data[0].n_iter >=20) printfignif(prob,file_output);
                               
								prob = (double)nunderalfa[1]/(double)data[0].n_iter;
								if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                if(data[0].n_loci > 2) printf("\n    P(variance of b = %g) = %g ",varobs[0],prob);
                                if(file_output)
                                    if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of b = %g) = %g ",varobs[0],prob);
                                if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                
								prob = (double)nunderalfa[2]/(double)data[0].n_iter;
								if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                printf("\n    P(average of q = %g) = %g ",(double)matrixml[0].Sq/ (double)data[0].n_loci,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of q = %g) = %g ",(double)matrixml[0].Sq/ (double)data[0].n_loci,prob);
                                if(data[0].n_iter >=20) printfignif(prob,file_output);
                                
								prob = (double)nunderalfa[3]/(double)data[0].n_iter;
								if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                if(data[0].n_loci > 2) printf("\n    P(variance of q = %g) = %g ",varobs[1],prob);
                                if(file_output)
                                    if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of q = %g) = %g ",varobs[1],prob);
                                if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                
								prob = (double)nunderalfa[4]/(double)data[0].n_iter;
								if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                printf("\n    P(average of za = %g) = %g ",(double)matrixml[0].Sza/ (double)data[0].n_loci,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of za = %g) = %g ",(double)matrixml[0].Sza/ (double)data[0].n_loci,prob);
                                if(data[0].n_iter >=20) printfignif(prob,file_output);
                                
								prob = (double)nunderalfa[5]/(double)data[0].n_iter;
								if((double)1./((double)data[0].n_iter) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                if(data[0].n_loci > 2) printf("\n    P(variance of za = %g) = %g ",varobs[2],prob);
                                if(file_output)
                                    if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of za = %g) = %g ",varobs[2],prob);
                                if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
								
								printf("\n");
								if(file_output) fputs("\n",file_output);
                            }
							#if ALLRESULTS == 0
                            else {/*give confident intervals with simulations. NOT observed_data.*/
							#endif
                                if((sortvector = (double *) malloc(data[0].n_iter * sizeof(double))) == 0) {
                                    printf("Probabilities can not be calculated, sorry.\n");
                                    if(file_output) fputs("Probabilities can not be calculated, sorry.\n",file_output);
                                    break;
                                }
                                printf("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n");
                                printf("statistic\tmedian\tavg\tvar\t");                                    
                                if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                printf("\n");
                                if(file_output) {
                                    fputs("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n",file_output);
                                    fputs("statistic\tmedian\tavg\tvar\t",file_output);                                    
                                    if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                    if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                    if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                    if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                    if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                    if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                    if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                    fprintf(file_output,"\n");
                                }

                                /*b*/
                                /*average*/
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if( matrixmlsim[it].Sb != (double) -10000) {
                                        sortvector[it2] = (double)matrixmlsim[it].Sb/(double)data[0].n_loci;
                                        it2 += 1;
                                        avg += (double)matrixmlsim[it].Sb/(double)data[0].n_loci;
                                        var += (double)matrixmlsim[it].Sb/(double)data[0].n_loci * (double)matrixmlsim[it].Sb/(double)data[0].n_loci;
                                    }
                                }
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_b\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_b\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                                /*variance*/
                                if(data[0].n_loci > 2) {
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if((vars=variances((double)matrixmlsim[it].S2b,(double)matrixmlsim[it].Sb,data[0].n_loci)) != (double) -10000) {
                                            sortvector[it2] = vars;
                                            it2 += 1;
                                            avg += vars;
                                            var += vars * vars;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_b\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_b\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");
                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                /*q*/
                                /*average*/
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if( matrixmlsim[it].Sq != (double) -10000) {
                                        sortvector[it2] = (double)matrixmlsim[it].Sq/(double)data[0].n_loci;
                                        it2 += 1;
                                        avg += (double)matrixmlsim[it].Sq/(double)data[0].n_loci;
                                        var += (double)matrixmlsim[it].Sq/(double)data[0].n_loci * (double)matrixmlsim[it].Sq/(double)data[0].n_loci;
                                    }
                                }
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_q\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_q\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                                /*variance*/
                                if(data[0].n_loci > 2) {
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if((vars=variances((double)matrixmlsim[it].S2q,(double)matrixmlsim[it].Sq,data[0].n_loci)) != (double) -10000) {
                                            sortvector[it2] = vars;
                                            it2 += 1;
                                            avg += vars;
                                            var += vars * vars;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_q\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_q\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");
                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                /*za*/
                                /*average*/
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if( matrixmlsim[it].Sza != (double) -10000) {
                                        sortvector[it2] = (double)matrixmlsim[it].Snhapl/(double)data[0].n_loci;
                                        it2 += 1;
                                        avg += (double)matrixmlsim[it].Sza/(double)data[0].n_loci;
                                        var += (double)matrixmlsim[it].Sza/(double)data[0].n_loci * (double)matrixmlsim[it].Sza/(double)data[0].n_loci;
                                    }
                                }
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_za\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_za\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                                /*variance*/
                                if(data[0].n_loci > 2) {
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if((vars=variances(matrixmlsim[it].S2za,matrixmlsim[it].Sza,data[0].n_loci)) != (double) -10000)
                                        {
                                            sortvector[it2] = vars;
                                            it2 += 1;
                                            avg += vars;
                                            var += vars * vars;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_za\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_za\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                free(sortvector);
							#if ALLRESULTS == 0
                            }
							#endif
                            break;
                        case '2': /* 2 - Display estimates of total locus variability*/
                            if(*observed_data && dataobsequalsim) {/*compare simulations with observed data.*/
                                if(*outgroup) {
                                    /*average and variance comparisons.*/
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ");
                                    printf("\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. ");
                                    printf("\n Multiple hits not included in the analysis.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ",file_output);
                                        fprintf(file_output,"\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. ");
                                        fputs("\n Multiple hits not included in the analysis.\n",file_output);
                                    }

                                    for(x=0;x<20;x++) nunderalfa[x] = 0;
                                    for(x=0;x<20;x++) nequivalfa[x] = 0;
                                    varobs[0] = variances((double)matrixml[0].S2theta_wat,(double)matrixml[0].Stheta_wat,data[0].n_loci);
                                    varobs[1] = variances((double)matrixml[0].S2theta_taj,(double)matrixml[0].Stheta_taj,data[0].n_loci);
                                    varobs[2] = variances((double)matrixml[0].S2theta_fuli,(double)matrixml[0].Stheta_fuli,data[0].n_loci);
                                    varobs[3] = variances((double)matrixml[0].S2theta_fw,(double)matrixml[0].Stheta_fw,data[0].n_loci);
                                    varobs[4] = variances((double)matrixml[0].S2theta_L,(double)matrixml[0].Stheta_L,data[0].n_loci);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        if(matrixml[0].Stheta_wat < matrixmlsim[it].Stheta_wat)
                                            nunderalfa[0] += 1;
                                        if(varobs[0] < variances((double)matrixmlsim[it].S2theta_wat,(double)matrixmlsim[it].Stheta_wat,data[0].n_loci))
                                            nunderalfa[1] += 1;
                                        if(matrixml[0].Stheta_taj < matrixmlsim[it].Stheta_taj)
                                            nunderalfa[2] += 1;
                                        if(varobs[1] < variances((double)matrixmlsim[it].S2theta_taj,(double)matrixmlsim[it].Stheta_taj,data[0].n_loci))
                                            nunderalfa[3] += 1;
                                        if(matrixml[0].Stheta_fuli < matrixmlsim[it].Stheta_fuli)
                                            nunderalfa[4] += 1;
                                        if(varobs[2] < variances((double)matrixmlsim[it].S2theta_fuli,(double)matrixmlsim[it].Stheta_fuli,data[0].n_loci))
                                            nunderalfa[5] += 1;
                                        if(matrixml[0].Stheta_fw < matrixmlsim[it].Stheta_fw)
                                            nunderalfa[6] += 1;
                                        if(varobs[3] < variances((double)matrixmlsim[it].S2theta_fw,(double)matrixmlsim[it].Stheta_fw,data[0].n_loci))
                                            nunderalfa[7] += 1;
                                        if(matrixml[0].Stheta_L < matrixmlsim[it].Stheta_L)
                                            nunderalfa[8] += 1;
                                        if(varobs[4] < variances((double)matrixmlsim[it].S2theta_L,(double)matrixmlsim[it].Stheta_L,data[0].n_loci))
                                            nunderalfa[9] += 1;
                                        if(matrixml[0].Stheta_wat == matrixmlsim[it].Stheta_wat)
                                            nequivalfa[0] += 1;
                                        if(varobs[0] == variances((double)matrixmlsim[it].S2theta_wat,(double)matrixmlsim[it].Stheta_wat,data[0].n_loci))
                                            nequivalfa[1] += 1;
                                        if(matrixml[0].Stheta_taj == matrixmlsim[it].Stheta_taj)
                                            nequivalfa[2] += 1;
                                        if(varobs[1] == variances((double)matrixmlsim[it].S2theta_taj,(double)matrixmlsim[it].Stheta_taj,data[0].n_loci))
                                            nequivalfa[3] += 1;
                                        if(matrixml[0].Stheta_fuli == matrixmlsim[it].Stheta_fuli)
                                            nequivalfa[4] += 1;
                                        if(varobs[2] == variances((double)matrixmlsim[it].S2theta_fuli,(double)matrixmlsim[it].Stheta_fuli,data[0].n_loci))
                                            nequivalfa[5] += 1;
                                        if(matrixml[0].Stheta_fw == matrixmlsim[it].Stheta_fw)
                                            nequivalfa[6] += 1;
                                        if(varobs[3] == variances((double)matrixmlsim[it].S2theta_fw,(double)matrixmlsim[it].Stheta_fw,data[0].n_loci))
                                            nequivalfa[7] += 1;
                                        if(matrixml[0].Stheta_L == matrixmlsim[it].Stheta_L)
                                            nequivalfa[8] += 1;
                                        if(varobs[4] == variances((double)matrixmlsim[it].S2theta_L,(double)matrixmlsim[it].Stheta_L,data[0].n_loci))
                                            nequivalfa[9] += 1;
                                    }
                                    prob = (double)nunderalfa[0]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                    printf("\n    P(average of theta_wat = %g) = %g ",(double)matrixml[0].Stheta_wat/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_wat = %g) = %g ",(double)matrixml[0].Stheta_wat/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[1]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_wat = %g) = %g ",varobs[0],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_wat = %g) = %g ",varobs[0],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[2]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                    printf("\n    P(average of theta_taj = %g) = %g ",(double)matrixml[0].Stheta_taj/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_taj = %g) = %g ",(double)matrixml[0].Stheta_taj/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[3]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_taj = %g) = %g ",varobs[1],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_taj = %g) = %g ",varobs[1],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[4]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                    printf("\n    P(average of theta_fuli = %g) = %g ",(double)matrixml[0].Stheta_fuli/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_fuli = %g) = %g ",(double)matrixml[0].Stheta_fuli/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[5]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_fuli = %g) = %g ",varobs[2],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_fuli = %g) = %g ",varobs[2],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[6]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[6] >= (double)0.5) nequivalfa[6] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[6];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[6];
                                    printf("\n    P(average of theta_fw = %g) = %g ",(double)matrixml[0].Stheta_fw/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_fw = %g) = %g ",(double)matrixml[0].Stheta_fw/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[7]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[7] >= (double)0.5) nequivalfa[7] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[7];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[7];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_fw = %g) = %g ",varobs[3],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_fw = %g) = %g ",varobs[3],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);

                                    prob = (double)nunderalfa[8]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[8] >= (double)0.5) nequivalfa[8] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[8];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[8];
                                    printf("\n    P(average of theta_Zeng = %g) = %g ",(double)matrixml[0].Stheta_L/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_Zeng = %g) = %g ",(double)matrixml[0].Stheta_L/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[9]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[9] >= (double)0.5) nequivalfa[9] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[9];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[9];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_Zeng = %g) = %g ",varobs[4],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_Zeng = %g) = %g ",varobs[4],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                }
                                else {/*not outgroup, observed_data*/
                                    /*average and variance comparisons.*/
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ");
                                    printf("\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. ");
                                    printf("\n Multiple hits not included in the analysis.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ",file_output);
                                        fprintf(file_output,"\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. ");
                                        fputs("\n Multiple hits not included in the analysis.\n",file_output);
                                    }
                                    for(x=0;x<20;x++) nunderalfa[x] = 0;
                                    for(x=0;x<20;x++) nequivalfa[x] = 0;
                                    varobs[0] = variances((double)matrixml[0].S2theta_wat,(double)matrixml[0].Stheta_wat,data[0].n_loci);
                                    varobs[1] = variances((double)matrixml[0].S2theta_taj,(double)matrixml[0].Stheta_taj,data[0].n_loci);
                                    varobs[2] = variances((double)matrixml[0].S2theta_fulin,(double)matrixml[0].Stheta_fulin,data[0].n_loci);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        if(matrixml[0].Stheta_wat < matrixmlsim[it].Stheta_wat)
                                            nunderalfa[0] += 1;
                                        if(varobs[0] < variances((double)matrixmlsim[it].S2theta_wat,(double)matrixmlsim[it].Stheta_wat,data[0].n_loci))
                                            nunderalfa[1] += 1;
                                        if(matrixml[0].Stheta_taj < matrixmlsim[it].Stheta_taj)
                                            nunderalfa[2] += 1;
                                        if(varobs[1] < variances((double)matrixmlsim[it].S2theta_taj,(double)matrixmlsim[it].Stheta_taj,data[0].n_loci))
                                            nunderalfa[3] += 1;
                                        if(matrixml[0].Stheta_fulin < matrixmlsim[it].Stheta_fulin)
                                            nunderalfa[4] += 1;
                                        if(varobs[2] < variances((double)matrixmlsim[it].S2theta_fulin,(double)matrixmlsim[it].Stheta_fulin,data[0].n_loci))
                                            nunderalfa[5] += 1;
                                        if(matrixml[0].Stheta_wat == matrixmlsim[it].Stheta_wat)
                                            nequivalfa[0] += 1;
                                        if(varobs[0] == variances((double)matrixmlsim[it].S2theta_wat,(double)matrixmlsim[it].Stheta_wat,data[0].n_loci))
                                            nequivalfa[1] += 1;
                                        if(matrixml[0].Stheta_taj == matrixmlsim[it].Stheta_taj)
                                            nequivalfa[2] += 1;
                                        if(varobs[1] == variances((double)matrixmlsim[it].S2theta_taj,(double)matrixmlsim[it].Stheta_taj,data[0].n_loci))
                                            nequivalfa[3] += 1;
                                        if(matrixml[0].Stheta_fulin == matrixmlsim[it].Stheta_fulin)
                                            nequivalfa[4] += 1;
                                        if(varobs[2] == variances((double)matrixmlsim[it].S2theta_fulin,(double)matrixmlsim[it].Stheta_fulin,data[0].n_loci))
                                            nequivalfa[5] += 1;
                                    }
                                    prob = (double)nunderalfa[0]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                    printf("\n    P(average of theta_wat = %g) = %g ",(double)matrixml[0].Stheta_wat/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_wat = %g) = %g ",(double)matrixml[0].Stheta_wat/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[1]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_wat = %g) = %g ",varobs[0],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_wat = %g) = %g ",varobs[0],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[2]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                    printf("\n    P(average of theta_taj = %g) = %g ",(double)matrixml[0].Stheta_taj/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_taj = %g) = %g ",(double)matrixml[0].Stheta_taj/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[3]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_taj = %g) = %g ",varobs[1],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_taj = %g) = %g ",varobs[1],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[4]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                    printf("\n    P(average of theta_fuli = %g) = %g ",(double)matrixml[0].Stheta_fulin/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_fuli = %g) = %g ",(double)matrixml[0].Stheta_fulin/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[5]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_fuli = %g) = %g ",varobs[2],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_fuli = %g) = %g ",varobs[2],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    if(file_output) {
                                        
                                    }
                                }
                            }
							#if ALLRESULTS == 0
                            else {/*give confident intervals with simulations. NOT observed_data.*/
							#endif
                                if((sortvector = (double *) malloc(data[0].n_iter * sizeof(double))) == 0) {
                                    printf("Probabilities can not be calculated, sorry.\n");
                                    if(file_output) fputs("Probabilities can not be calculated, sorry.\n",file_output);
                                    break;
                                }
                                if(*outgroup) {
                                    printf("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n");
                                    printf("statistic\tmedian\tavg\tvar\t");                                    
                                    if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                    if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                    if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                    if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                    if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                    if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                    if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                    printf("\n");
                                    if(file_output) {
                                        fputs("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n",file_output);
                                        fputs("statistic\tmedian\tavg\tvar\t",file_output);                                    
                                        if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                        if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                        if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                        if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                        if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                        if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                        if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                        fprintf(file_output,"\n");
                                    }

                                    /*Stheta_wat. */
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_wat != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_wat/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_wat/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_wat/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_wat/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_wat\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_wat\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_wat,matrixmlsim[it].Stheta_wat,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_wat\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_wat\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Stheta_taj*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_taj != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_taj/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_taj/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_taj/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_taj/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_taj\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_taj\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_taj,matrixmlsim[it].Stheta_taj,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_taj\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_taj\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Stheta_fuli*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_fuli != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_fuli/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_fuli/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_fuli/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_fuli/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_fuli\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_fuli\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_fuli,matrixmlsim[it].Stheta_fuli,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_fuli\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_fuli\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Stheta_fw*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_fw != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_fw/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_fw/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_fw/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_fw/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_fw\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_fw\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_fw,matrixmlsim[it].Stheta_fw,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_fw\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_fw\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Stheta_zeng*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_L != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_L/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_L/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_L/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_L/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_zeng\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_zeng\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_L,matrixmlsim[it].Stheta_L,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_zeng\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_zeng\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                }
                                else {
                                    /*Do a NEW vector with all the values, SORT.*/
                                    /* Take the median (50%), average and variance, and 10%, 5%, 2.5%, 1%, 0.5%, 0.1% and 0.05% for 1 tails.*/
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n");
                                    printf("statistic\tmedian\tavg\tvar\t");                                    
                                    if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                    if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                    if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                    if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                    if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                    if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                    if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                    printf("\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n",file_output);
                                        fputs("statistic\tmedian\tavg\tvar\t",file_output);                                    
                                        if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                        if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                        if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                        if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                        if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                        if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                        if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                        fprintf(file_output,"\n");
                                    }
                                    
                                    /*Stheta_wat. */
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_wat != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_wat/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_wat/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_wat/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_wat/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_wat\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_wat\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_wat,matrixmlsim[it].Stheta_wat,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_wat\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_wat\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Stheta_taj*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_taj != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_taj/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_taj/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_taj/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_taj/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_taj\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_taj\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_taj,matrixmlsim[it].Stheta_taj,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_taj\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_taj\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Stheta_fulin*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_fulin != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_fulin/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_fulin/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_fulin/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_fulin/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_fulin\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_fulin\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_fulin,matrixmlsim[it].Stheta_fulin,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_fulin\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_fulin\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                }
                                free(sortvector);
							#if ALLRESULTS == 0
                            }
							#endif
                            break;
                        case '3': /* 3 - Display estimates of nucleotide locus variability*/
                            if(*observed_data && dataobsequalsim) {/*compare simulations with observed data.*/
                                if(*outgroup) {
                                    /*average and variance comparisons.*/
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ");
                                    printf("\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. ");
                                    printf("\n Multiple hits not included in the analysis.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ",file_output);
                                        fprintf(file_output,"\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. ");
                                        fputs("\n Multiple hits not included in the analysis.\n",file_output);
                                    }

                                    for(x=0;x<20;x++) nunderalfa[x] = 0;
                                    for(x=0;x<20;x++) nequivalfa[x] = 0;
                                    varobs[0] = variances((double)matrixml[0].S2theta_wat_nut,(double)matrixml[0].Stheta_wat_nut,data[0].n_loci);
                                    varobs[1] = variances((double)matrixml[0].S2theta_taj_nut,(double)matrixml[0].Stheta_taj_nut,data[0].n_loci);
                                    varobs[2] = variances((double)matrixml[0].S2theta_fuli_nut,(double)matrixml[0].Stheta_fuli_nut,data[0].n_loci);
                                    varobs[3] = variances((double)matrixml[0].S2theta_fw_nut,(double)matrixml[0].Stheta_fw_nut,data[0].n_loci);
									varobs[4] = variances((double)matrixml[0].S2theta_L_nut,(double)matrixml[0].Stheta_L_nut,data[0].n_loci);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        if(matrixml[0].Stheta_wat_nut < matrixmlsim[it].Stheta_wat_nut)
                                            nunderalfa[0] += 1;
                                        if(varobs[0] < variances((double)matrixmlsim[it].S2theta_wat_nut,(double)matrixmlsim[it].Stheta_wat_nut,data[0].n_loci))
                                            nunderalfa[1] += 1;
                                        if(matrixml[0].Stheta_taj_nut < matrixmlsim[it].Stheta_taj_nut)
                                            nunderalfa[2] += 1;
                                        if(varobs[1] < variances((double)matrixmlsim[it].S2theta_taj_nut,(double)matrixmlsim[it].Stheta_taj_nut,data[0].n_loci))
                                            nunderalfa[3] += 1;
                                        if(matrixml[0].Stheta_fuli_nut < matrixmlsim[it].Stheta_fuli_nut)
                                            nunderalfa[4] += 1;
                                        if(varobs[2] < variances((double)matrixmlsim[it].S2theta_fuli_nut,(double)matrixmlsim[it].Stheta_fuli_nut,data[0].n_loci))
                                            nunderalfa[5] += 1;
                                        if(matrixml[0].Stheta_fw_nut < matrixmlsim[it].Stheta_fw_nut)
                                            nunderalfa[6] += 1;
                                        if(varobs[3] < variances((double)matrixmlsim[it].S2theta_fw_nut,(double)matrixmlsim[it].Stheta_fw_nut,data[0].n_loci))
                                            nunderalfa[7] += 1;
                                        if(matrixml[0].Stheta_L_nut < matrixmlsim[it].Stheta_L_nut)
                                            nunderalfa[8] += 1;
                                        if(varobs[4] < variances((double)matrixmlsim[it].S2theta_L_nut,(double)matrixmlsim[it].Stheta_L_nut,data[0].n_loci))
                                            nunderalfa[9] += 1;
                                        if(matrixml[0].Stheta_wat_nut == matrixmlsim[it].Stheta_wat_nut)
                                            nequivalfa[0] += 1;
                                        if(varobs[0] == variances((double)matrixmlsim[it].S2theta_wat_nut,(double)matrixmlsim[it].Stheta_wat_nut,data[0].n_loci))
                                            nequivalfa[1] += 1;
                                        if(matrixml[0].Stheta_taj_nut == matrixmlsim[it].Stheta_taj_nut)
                                            nequivalfa[2] += 1;
                                        if(varobs[1] == variances((double)matrixmlsim[it].S2theta_taj_nut,(double)matrixmlsim[it].Stheta_taj_nut,data[0].n_loci))
                                            nequivalfa[3] += 1;
                                        if(matrixml[0].Stheta_fuli_nut == matrixmlsim[it].Stheta_fuli_nut)
                                            nequivalfa[4] += 1;
                                        if(varobs[2] == variances((double)matrixmlsim[it].S2theta_fuli_nut,(double)matrixmlsim[it].Stheta_fuli_nut,data[0].n_loci))
                                            nequivalfa[5] += 1;
                                        if(matrixml[0].Stheta_fw_nut == matrixmlsim[it].Stheta_fw_nut)
                                            nequivalfa[6] += 1;
                                        if(varobs[3] == variances((double)matrixmlsim[it].S2theta_fw_nut,(double)matrixmlsim[it].Stheta_fw_nut,data[0].n_loci))
                                            nequivalfa[7] += 1;
                                        if(matrixml[0].Stheta_L_nut == matrixmlsim[it].Stheta_L_nut)
                                            nequivalfa[8] += 1;
                                        if(varobs[4] == variances((double)matrixmlsim[it].S2theta_L_nut,(double)matrixmlsim[it].Stheta_L_nut,data[0].n_loci))
                                            nequivalfa[9] += 1;
                                    }
                                    prob = (double)nunderalfa[0]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                    printf("\n    P(average of theta_wat/nt = %g) = %g ",(double)matrixml[0].Stheta_wat_nut/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_wat/nt = %g) = %g ",(double)matrixml[0].Stheta_wat_nut/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[1]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_wat/nt = %g) = %g ",varobs[0],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_wat/nt = %g) = %g ",varobs[0],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[2]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                    printf("\n    P(average of theta_taj/nt = %g) = %g ",(double)matrixml[0].Stheta_taj_nut/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_taj/nt = %g) = %g ",(double)matrixml[0].Stheta_taj_nut/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[3]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_taj/nt = %g) = %g ",varobs[1],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_taj/nt = %g) = %g ",varobs[1],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[4]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                    printf("\n    P(average of theta_fuli/nt = %g) = %g ",(double)matrixml[0].Stheta_fuli_nut/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_fuli/nt = %g) = %g ",(double)matrixml[0].Stheta_fuli_nut/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[5]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_fuli/nt = %g) = %g ",varobs[2],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_fuli/nt = %g) = %g ",varobs[2],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[6]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[6] >= (double)0.5) nequivalfa[6] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[6];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[6];
                                    printf("\n    P(average of theta_fw/nt = %g) = %g ",(double)matrixml[0].Stheta_fw_nut/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_fw/nt = %g) = %g ",(double)matrixml[0].Stheta_fw_nut/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[7]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[7] >= (double)0.5) nequivalfa[7] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[7];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[7];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_fw/nt = %g) = %g ",varobs[3],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_fw/nt = %g) = %g ",varobs[3],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
                                    prob = (double)nunderalfa[8]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[8] >= (double)0.5) nequivalfa[8] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[8];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[8];
                                    printf("\n    P(average of theta_zeng/nt = %g) = %g ",(double)matrixml[8].Stheta_L_nut/ (double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_zeng/nt = %g) = %g ",(double)matrixml[0].Stheta_L_nut/ (double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[9]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[9] >= (double)0.5) nequivalfa[9] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[9];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[9];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_zeng/nt = %g) = %g ",varobs[4],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_zeng/nt = %g) = %g ",varobs[4],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                }
                                else {/*not outgroup, observed_data*/
                                    /*average and variance comparisons.*/
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ");
                                    printf("\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. ");
                                    printf("\n Multiple hits not included in the analysis.\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\n\n Display the probability of observed multilocus statistics obtained by coalescent Monte Carlo simulations: ",file_output);
                                        fprintf(file_output,"\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. ");
                                        fputs("\n Multiple hits not included in the analysis.\n",file_output);
                                    }
                                    nunderalfa[0] = nunderalfa[1] = nunderalfa[2] = nunderalfa[3] = nunderalfa[4] = nunderalfa[5] = 0;
                                    nequivalfa[0] = nequivalfa[1] = nequivalfa[2] = nequivalfa[3] = nequivalfa[4] = nequivalfa[5] = 0;
                                    varobs[0] = variances((double)matrixml[0].S2theta_wat_nut,(double)matrixml[0].Stheta_wat_nut,data[0].n_loci);
                                    varobs[1] = variances((double)matrixml[0].S2theta_taj_nut,(double)matrixml[0].Stheta_taj_nut,data[0].n_loci);
                                    varobs[2] = variances((double)matrixml[0].S2theta_fulin_nut,(double)matrixml[0].Stheta_fulin_nut,data[0].n_loci);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        if(matrixml[0].Stheta_wat_nut < matrixmlsim[it].Stheta_wat_nut)
                                            nunderalfa[0] += 1;
                                        if(varobs[0] < variances((double)matrixmlsim[it].S2theta_wat_nut,(double)matrixmlsim[it].Stheta_wat_nut,data[0].n_loci))
                                            nunderalfa[1] += 1;
                                        if(matrixml[0].Stheta_taj_nut < matrixmlsim[it].Stheta_taj_nut)
                                            nunderalfa[2] += 1;
                                        if(varobs[1] < variances((double)matrixmlsim[it].S2theta_taj_nut,(double)matrixmlsim[it].Stheta_taj_nut,data[0].n_loci))
                                            nunderalfa[3] += 1;
                                        if(matrixml[0].Stheta_fulin_nut < matrixmlsim[it].Stheta_fulin_nut)
                                            nunderalfa[4] += 1;
                                        if(varobs[2] < variances((double)matrixmlsim[it].S2theta_fulin_nut,(double)matrixmlsim[it].Stheta_fulin_nut,data[0].n_loci))
                                            nunderalfa[5] += 1;
                                        if(matrixml[0].Stheta_wat_nut == matrixmlsim[it].Stheta_wat_nut)
                                            nequivalfa[0] += 1;
                                        if(varobs[0] == variances((double)matrixmlsim[it].S2theta_wat_nut,(double)matrixmlsim[it].Stheta_wat_nut,data[0].n_loci))
                                            nequivalfa[1] += 1;
                                        if(matrixml[0].Stheta_taj_nut == matrixmlsim[it].Stheta_taj_nut)
                                            nequivalfa[2] += 1;
                                        if(varobs[1] == variances((double)matrixmlsim[it].S2theta_taj_nut,(double)matrixmlsim[it].Stheta_taj_nut,data[0].n_loci))
                                            nequivalfa[3] += 1;
                                        if(matrixml[0].Stheta_fulin_nut == matrixmlsim[it].Stheta_fulin_nut)
                                            nequivalfa[4] += 1;
                                        if(varobs[2] == variances((double)matrixmlsim[it].S2theta_fulin_nut,(double)matrixmlsim[it].Stheta_fulin_nut,data[0].n_loci))
                                            nequivalfa[5] += 1;
                                    }
                                    prob = (double)nunderalfa[0]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[0];
                                    printf("\n    P(average of theta_wat/nt = %g) = %g ",(double)matrixml[0].Stheta_wat_nut/(double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_wat/nt = %g) = %g ",(double)matrixml[0].Stheta_wat_nut/(double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[1]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[1];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_wat/nt = %g) = %g ",varobs[0],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_wat/nt = %g) = %g ",varobs[0],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[2]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[2];
                                    printf("\n    P(average of theta_taj/nt = %g) = %g ",(double)matrixml[0].Stheta_taj_nut/(double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_taj/nt = %g) = %g ",(double)matrixml[0].Stheta_taj_nut/(double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[3]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[3];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_taj/nt = %g) = %g ",varobs[1],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_taj/nt = %g) = %g ",varobs[1],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[4]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[4];
                                    printf("\n    P(average of theta_fuli/nt = %g) = %g ",(double)matrixml[0].Stheta_fulin_nut/(double)data[0].n_loci,prob);
                                    if(file_output)
                                        fprintf(file_output,"\n    P(average of theta_fuli/nt = %g) = %g ",(double)matrixml[0].Stheta_fulin_nut/(double)data[0].n_loci,prob);
                                    if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
									prob = (double)nunderalfa[5]/(double)data[0].n_iter;
									if((double)1./((double)data[0].n_iter) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                    if(prob > (double)0.5) prob -= (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                    else prob += (double)1./((double)data[0].n_iter) *(double)nequivalfa[5];
                                    if(data[0].n_loci > 2) printf("\n    P(variance of theta_fuli/nt = %g) = %g ",varobs[2],prob);
                                    if(file_output)
                                        if(data[0].n_loci > 2) fprintf(file_output,"\n    P(variance of theta_fuli/nt = %g) = %g ",varobs[2],prob);
                                    if(data[0].n_loci > 2) if(data[0].n_iter >=20) printfignif(prob,file_output);
                                    
                                }
                            }
							#if ALLRESULTS == 0
                            else {/*give confident intervals with simulations. NOT observed_data.*/
							#endif
                                if((sortvector = (double *) malloc(data[0].n_iter * sizeof(double))) == 0) {
                                    printf("Probabilities can not be calculated, sorry.\n");
                                    if(file_output) fputs("Probabilities can not be calculated, sorry.\n",file_output);
                                    break;
                                }
                                if(*outgroup) {
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n");
                                    printf("statistic\tmedian\tavg\tvar\t");                                    
                                    if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                    if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                    if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                    if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                    if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                    if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                    if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                    printf("\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n",file_output);
                                        fputs("statistic\tmedian\tavg\tvar\t",file_output);                                    
                                        if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                        if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                        if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                        if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                        if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                        if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                        if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                        fprintf(file_output,"\n");
                                    }

                                    /*Stheta_wat. */
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_wat_nut != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_wat_nut/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_wat_nut/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_wat_nut/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_wat_nut/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_wat/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_wat/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_wat_nut,matrixmlsim[it].Stheta_wat_nut,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_wat/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_wat/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Stheta_taj*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_taj_nut != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_taj_nut/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_taj_nut/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_taj_nut/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_taj_nut/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_taj/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_taj/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_taj_nut,matrixmlsim[it].Stheta_taj_nut,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_taj/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_taj/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Stheta_fuli*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_fuli_nut != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_fuli_nut/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_fuli_nut/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_fuli_nut/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_fuli_nut/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_fuli/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_fuli/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_fuli_nut,matrixmlsim[it].Stheta_fuli_nut,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_fuli/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_fuli/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Stheta_fw*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_fw_nut != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_fw_nut/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_fw_nut/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_fw_nut/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_fw_nut/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_fw/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_fw/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_fw_nut,matrixmlsim[it].Stheta_fw_nut,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_fw/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_fw/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");
										}
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Stheta_zeng*/
                                    /*average*/
                                    avg = var = (double)0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_L_nut != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_L_nut/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_L_nut/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_L_nut/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_L_nut/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_zeng/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_zeng/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_L_nut,matrixmlsim[it].Stheta_L_nut,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_zeng/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_zeng/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");
										}
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                }
                                else {
                                    /*Do a NEW vector with all the values, SORT.*/
                                    /* Take the median (50%), average and variance, and 10%, 5%, 2.5%, 1%, 0.5%, 0.1% and 0.05% for 1 tails.*/
									printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                    printf("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n");
                                    printf("statistic\tmedian\tavg\tvar\t");                                    
                                    if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                                    if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                                    if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                                    if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                                    if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                                    if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                                    if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                                    printf("\n");
                                    if(file_output) {
										if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                        fputs("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n",file_output);
                                        fputs("statistic\tmedian\tavg\tvar\t",file_output);                                    
                                        if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                        if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                        if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                        if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                        if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                        if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                        if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                        fprintf(file_output,"\n");
                                    }
                                    /*Stheta_wat. */
                                    /*average*/
                                    avg = var = 0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_wat_nut != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_wat_nut/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_wat_nut/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_wat_nut/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_wat_nut/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_wat/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_wat/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = 0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_wat_nut,matrixmlsim[it].Stheta_wat_nut,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_wat/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_wat/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Stheta_taj*/
                                    /*average*/
                                    avg = var = 0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_taj_nut != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_taj_nut/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_taj_nut/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_taj_nut/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_taj_nut/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_taj/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_taj/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = 0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_taj_nut,matrixmlsim[it].Stheta_taj_nut,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_taj/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_taj/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                    /*Stheta_fulin*/
                                    /*average*/
                                    avg = var = 0.0;
                                    for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                        if( matrixmlsim[it].Stheta_fulin_nut != (double) -10000) {
                                            sortvector[it2] = (double)matrixmlsim[it].Stheta_fulin_nut/(double)data[0].n_loci;
                                            it2 += 1;
                                            avg += (double)matrixmlsim[it].Stheta_fulin_nut/(double)data[0].n_loci;
                                            var += (double)matrixmlsim[it].Stheta_fulin_nut/(double)data[0].n_loci * (double)matrixmlsim[it].Stheta_fulin_nut/(double)data[0].n_loci;
                                        }
                                    }
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("avg_theta_fuli/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"avg_theta_fuli/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                    /*variance*/
                                    if(data[0].n_loci > 2) {
                                        avg = var = (double)0.0;
                                        for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                            if((vars=variances(matrixmlsim[it].S2theta_fulin_nut,matrixmlsim[it].Stheta_fulin_nut,data[0].n_loci)) != (double) -10000) {
                                                sortvector[it2] = vars;
                                                it2 += 1;
                                                avg += vars;
                                                var += vars * vars;
                                            }
                                        }
                                        var = varianceit(var,avg,it2);
                                        if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                        printf("var_theta_fulin/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(file_output)
                                            fprintf(file_output,"var_theta_fulin/nt\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                        if(var != (double)-10000) {
                                            printf("%g\t",var);
                                            if(file_output) fprintf(file_output,"%g\t",var);
                                        }
                                        else {
                                            printf("na\t");
                                            if(file_output) fprintf(file_output,"na\t");                                    }
                                        print_percentages(data,sortvector,it2,file_output);
                                    }
                                }
                                free(sortvector);
							#if ALLRESULTS == 0
							}
							#endif                            
							break;
                        case '4': /* 4 - Back to previous menu*/
                            return;
                            break;
                    }
                    break;
               case 1:	/* neutrality tests, multilocus*/
					if(file_output) {    
						fprintf(file_output,"\n\n     MENU:\n     3. Statistical inference based on Coalescent Monte Carlo simulations menu:");
						fprintf(file_output,"\n     3.1. Display coalescent simulation analysis menu:");
						fprintf(file_output,"\n     3.1.1. Display neutrality tests menu:");
						fprintf(file_output,"\n     3.1.1.0. Display multilocus neutrality tests menu:\n\n");
						fprintf(file_output,"\n     3.1.1.0.1. Display the probability values (from observed data) or confident intervals.\n");
					}

                    if(*observed_data && dataobsequalsim) {/*compare simulations with observed data.*/
                        if(*outgroup) {
                            /*average and variance comparisons.*/
							printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                            printf("\n\n Display the probability of observed multilocus neutrality tests obtained by coalescent Monte Carlo simulations: ");
                            printf("\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%, except for HKA test (one tail). ");
                            printf("\n Multiple hits not included in the analysis.\n");
                            printf("\n Shared polymorphisms (if they are) are not included in the analysis, except in HKA test (if performed).\n");
                            if(file_output) {
								if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                fputs("\n\n Display the probability of observed multilocus neutrality tests obtained by coalescent Monte Carlo simulations: ",file_output);
                                fprintf(file_output,"\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%, except for HKA test (one tail). ");
                                fputs("\n Multiple hits not included in the analysis.\n",file_output);
                                fputs("\n Shared polymorphisms (if they are) are not included in the analysis, except in HKA test (if performed).\n",file_output);
                            }

                            for(x=0;x<30;x++) nunderalfa[x] = 0;
                            for(x=0;x<30;x++) nequivalfa[x] = 0;
                            for(x=0;x<30;x++) nitertests[x] = 0;
                            
                            varobs[0] = variances((double)matrixml[0].S2tajimaD,(double)matrixml[0].StajimaD,matrixml[0].nltajD);
                            varobs[1] = variances((double)matrixml[0].S2fuliD,(double)matrixml[0].SfuliD,matrixml[0].nlflD);
                            varobs[2] = variances((double)matrixml[0].S2fuliF,(double)matrixml[0].SfuliF,matrixml[0].nlflF);
                            varobs[3] = variances((double)matrixml[0].S2fuFs,(double)matrixml[0].SfuFs,matrixml[0].nlFs);
                            varobs[4] = variances((double)matrixml[0].S2faywuH,(double)matrixml[0].SfaywuH,matrixml[0].nlH);
                            varobs[5] = variances((double)matrixml[0].S2faywuHo,(double)matrixml[0].SfaywuHo,matrixml[0].nlHo);
                            varobs[6] = variances((double)matrixml[0].S2rZA,(double)matrixml[0].SrZA,matrixml[0].nlZ);
                            varobs[7] = variances((double)matrixml[0].S2wB,(double)matrixml[0].SwB,matrixml[0].nlB);
                            varobs[8] = variances((double)matrixml[0].S2wQ,(double)matrixml[0].SwQ,matrixml[0].nlQ);
                            varobs[9] = variances((double)matrixml[0].S2R2,(double)matrixml[0].SR2,matrixml[0].nlR2);
                            varobs[10] = variances((double)matrixml[0].S2zengE,(double)matrixml[0].SzengE,matrixml[0].nlE);
                            varobs[11] = variances((double)matrixml[0].S2ewtest,(double)matrixml[0].Sewtest,data[0].n_loci);
                            
                            for(it=0;it<(long int)data[0].n_iter;it++) {
                                if(matrixml[0].StajimaD != (double)-10000 && matrixmlsim[it].StajimaD != (double)-10000) {
                                    if(matrixml[0].StajimaD/(double)matrixml[0].nltajD <
                                        matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD)
                                        nunderalfa[0] += 1;
                                    if(matrixml[0].StajimaD/(double)matrixml[0].nltajD ==
                                        matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD)
                                        nequivalfa[0] += 1;
                                    nitertests[0] += 1;
                                    if(matrixmlsim[it].nltajD > 2) {
                                        if(varobs[0] < variances((double)matrixmlsim[it].S2tajimaD,(double)matrixmlsim[it].StajimaD,matrixmlsim[it].nltajD))
                                            nunderalfa[1] += 1;
                                        if(varobs[0] == variances((double)matrixmlsim[it].S2tajimaD,(double)matrixmlsim[it].StajimaD,matrixmlsim[it].nltajD))
                                            nequivalfa[1] += 1;
                                        nitertests[1] += 1;
                                    }
                                }
                                if(matrixml[0].SfuliD != (double)-10000 && matrixmlsim[it].SfuliD != (double)-10000) {
                                    if(matrixml[0].SfuliD/(double)matrixml[0].nlflD <
                                        matrixmlsim[it].SfuliD/(double)matrixmlsim[it].nlflD)
                                        nunderalfa[2] += 1;
                                    if(matrixml[0].SfuliD/(double)matrixml[0].nlflD == 
                                        matrixmlsim[it].SfuliD/(double)matrixmlsim[it].nlflD)
                                        nequivalfa[2] += 1;
                                    nitertests[2] += 1;
                                    if(matrixmlsim[it].nlflD > 2) {
                                        if(varobs[1] < variances((double)matrixmlsim[it].S2fuliD,(double)matrixmlsim[it].SfuliD,matrixmlsim[it].nlflD))
                                            nunderalfa[3] += 1;
                                        if(varobs[1] == variances((double)matrixmlsim[it].S2fuliD,(double)matrixmlsim[it].SfuliD,matrixmlsim[it].nlflD))
                                            nequivalfa[3] += 1;
                                        nitertests[3] += 1;
                                    }
                                }
                                if(matrixml[0].SfuliF != (double)-10000 && matrixmlsim[it].SfuliF != (double)-10000) {
                                    if(matrixml[0].SfuliF/(double)matrixml[0].nlflF < 
                                        matrixmlsim[it].SfuliF/(double)matrixmlsim[it].nlflF)
                                        nunderalfa[4] += 1;
                                    if(matrixml[0].SfuliF/(double)matrixml[0].nlflF == 
                                        matrixmlsim[it].SfuliF/(double)matrixmlsim[it].nlflF)
                                        nequivalfa[4] += 1;
                                    nitertests[4] += 1;
                                    if(matrixmlsim[it].nlflF > 2) {
                                        if(varobs[2] < variances((double)matrixmlsim[it].S2fuliF,(double)matrixmlsim[it].SfuliF,matrixmlsim[it].nlflF))
                                            nunderalfa[5] += 1;
                                        if(varobs[2] == variances((double)matrixmlsim[it].S2fuliF,(double)matrixmlsim[it].SfuliF,matrixmlsim[it].nlflF))
                                            nequivalfa[5] += 1;
                                        nitertests[5] += 1;
                                    }
                                }
                                if(matrixml[0].SfuFs != (double)-10000 && matrixmlsim[it].SfuFs != (double)-10000) {
                                    if(matrixml[0].SfuFs/(double)matrixml[0].nlFs < 
                                        matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs)
                                        nunderalfa[6] += 1;
                                    if(matrixml[0].SfuFs/(double)matrixml[0].nlFs == 
                                        matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs)
                                        nequivalfa[6] += 1;
                                    nitertests[6] += 1;
                                    if(matrixmlsim[it].nlFs > 2) {
                                        if(varobs[3] < variances((double)matrixmlsim[it].S2fuFs,(double)matrixmlsim[it].SfuFs,matrixmlsim[it].nlFs))
                                            nunderalfa[7] += 1;
                                        if(varobs[3] == variances((double)matrixmlsim[it].S2fuFs,(double)matrixmlsim[it].SfuFs,matrixmlsim[it].nlFs))
                                            nequivalfa[7] += 1;
                                        nitertests[7] += 1;
                                    }
                                }
                                if(matrixml[0].SfaywuH != (double)-10000 && matrixmlsim[it].SfaywuH != (double)-10000) {
                                    if(matrixml[0].SfaywuH/(double)matrixml[0].nlH <
                                        matrixmlsim[it].SfaywuH/(double)matrixmlsim[it].nlH)
                                        nunderalfa[8] += 1;
                                    if(matrixml[0].SfaywuH/(double)matrixml[0].nlH ==
                                        matrixmlsim[it].SfaywuH/(double)matrixmlsim[it].nlH)
                                        nequivalfa[8] += 1;
                                    nitertests[8] += 1;
                                    if(matrixmlsim[it].nlH > 2) {
                                        if(varobs[4] < variances((double)matrixmlsim[it].S2faywuH,(double)matrixmlsim[it].SfaywuH,matrixmlsim[it].nlH))
                                            nunderalfa[9] += 1;
                                        if(varobs[4] == variances((double)matrixmlsim[it].S2faywuH,(double)matrixmlsim[it].SfaywuH,matrixmlsim[it].nlH))
                                            nequivalfa[9] += 1;
                                        nitertests[9] += 1;
                                    }
                                }
                                if(matrixml[0].SfaywuHo != (double)-10000 && matrixmlsim[it].SfaywuHo != (double)-10000) {
                                    if(matrixml[0].SfaywuHo/(double)matrixml[0].nlHo <
                                        matrixmlsim[it].SfaywuHo/(double)matrixmlsim[it].nlHo)
                                        nunderalfa[10] += 1;
                                    if(matrixml[0].SfaywuHo/(double)matrixml[0].nlHo ==
                                        matrixmlsim[it].SfaywuHo/(double)matrixmlsim[it].nlHo)
                                        nequivalfa[10] += 1;
                                    nitertests[10] += 1;
                                    if(matrixmlsim[it].nlHo > 2) {
                                        if(varobs[5] < variances((double)matrixmlsim[it].S2faywuHo,(double)matrixmlsim[it].SfaywuHo,matrixmlsim[it].nlHo))
                                            nunderalfa[11] += 1;
                                        if(varobs[5] == variances((double)matrixmlsim[it].S2faywuHo,(double)matrixmlsim[it].SfaywuHo,matrixmlsim[it].nlHo))
                                            nequivalfa[11] += 1;
                                        nitertests[11] += 1;
                                    }
                                }
                                if(matrixml[0].SrZA != (double)-10000 && matrixmlsim[it].SrZA != (double)-10000) {
                                    if(matrixml[0].SrZA/(double)matrixml[0].nlZ < matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ)
                                        nunderalfa[12] += 1;
                                    if(matrixml[0].SrZA/(double)matrixml[0].nlZ == matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ)
                                        nequivalfa[12] += 1;
                                    nitertests[12] += 1;
                                    if(matrixmlsim[it].nlZ > 2) {
                                        if(varobs[6] < variances((double)matrixmlsim[it].S2rZA,(double)matrixmlsim[it].SrZA,matrixmlsim[it].nlZ))
                                            nunderalfa[13] += 1;
                                        if(varobs[6] == variances((double)matrixmlsim[it].S2rZA,(double)matrixmlsim[it].SrZA,matrixmlsim[it].nlZ))
                                            nequivalfa[13] += 1;
                                        nitertests[13] += 1;
                                    }
                                }
                                if(matrixml[0].SwB != (double)-10000 && matrixmlsim[it].SwB != (double)-10000) {
                                    if(matrixml[0].SwB/(double)matrixml[0].nlB < matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB)
                                        nunderalfa[14] += 1;
                                    if(matrixml[0].SwB/(double)matrixml[0].nlB == matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB)
                                        nequivalfa[14] += 1;
                                    nitertests[14] += 1;
                                    if(matrixmlsim[it].nlB > 2) {
                                        if(varobs[7] < variances((double)matrixmlsim[it].S2wB,(double)matrixmlsim[it].SwB,matrixmlsim[it].nlB))
                                            nunderalfa[15] += 1;
                                        if(varobs[7] == variances((double)matrixmlsim[it].S2wB,(double)matrixmlsim[it].SwB,matrixmlsim[it].nlB))
                                            nequivalfa[15] += 1;
                                        nitertests[15] += 1;
                                    }
                                }
                                if(matrixml[0].SwQ != (double)-10000 && matrixmlsim[it].SwQ != (double)-10000) {
                                    if(matrixml[0].SwQ/(double)matrixml[0].nlQ <
                                        matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ)
                                        nunderalfa[16] += 1;
                                    if(matrixml[0].SwQ/(double)matrixml[0].nlQ ==
                                        matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ)
                                        nequivalfa[16] += 1;
                                    nitertests[16] += 1;
                                    if(matrixmlsim[it].nlQ > 2) {
                                        if(varobs[8] < variances((double)matrixmlsim[it].S2wQ,(double)matrixmlsim[it].SwQ,matrixmlsim[it].nlQ))
                                            nunderalfa[17] += 1;
                                        if(varobs[8] == variances((double)matrixmlsim[it].S2wQ,(double)matrixmlsim[it].SwQ,matrixmlsim[it].nlQ))
                                            nequivalfa[17] += 1;
                                        nitertests[17] += 1;
                                    }
                                }
                                if(matrixml[0].SR2 != (double)-10000 && matrixmlsim[it].SR2 != (double)-10000) {
                                    if(matrixml[0].SR2/(double)matrixml[0].nlR2 <
                                        matrixmlsim[it].SR2/matrixmlsim[it].nlR2)
                                        nunderalfa[18] += 1;
                                    if(matrixml[0].SR2/(double)matrixml[0].nlR2 == matrixmlsim[it].SR2/matrixmlsim[it].nlR2)
                                        nequivalfa[18] += 1;
                                    nitertests[18] += 1;
                                    if( matrixmlsim[it].nlR2 > 2) {
                                        if(varobs[9] < variances((double)matrixmlsim[it].S2R2,(double)matrixmlsim[it].SR2,matrixmlsim[it].nlR2))
                                            nunderalfa[19] += 1;
                                        if(varobs[9] == variances((double)matrixmlsim[it].S2R2,(double)matrixmlsim[it].SR2,matrixmlsim[it].nlR2))
                                            nequivalfa[19] += 1;
                                        nitertests[19] += 1;
                                    }
                                }

                                if(matrixml[0].SzengE != (double)-10000 && matrixmlsim[it].SzengE != (double)-10000) {
                                    if(matrixml[0].SzengE/(double)matrixml[0].nlE <
                                        matrixmlsim[it].SzengE/matrixmlsim[it].nlE)
                                        nunderalfa[20] += 1;
                                    if(matrixml[0].SzengE/(double)matrixml[0].nlE == matrixmlsim[it].SzengE/matrixmlsim[it].nlE)
                                        nequivalfa[20] += 1;
                                    nitertests[20] += 1;
                                    if( matrixmlsim[it].nlE > 2) {
                                        if(varobs[10] < variances((double)matrixmlsim[it].S2zengE,(double)matrixmlsim[it].SzengE,matrixmlsim[it].nlE))
                                            nunderalfa[21] += 1;
                                        if(varobs[10] == variances((double)matrixmlsim[it].S2zengE,(double)matrixmlsim[it].SzengE,matrixmlsim[it].nlE))
                                            nequivalfa[21] += 1;
                                        nitertests[21] += 1;
                                    }
                                }
                                if(matrixml[0].Sewtest != (double)-10000 && matrixmlsim[it].Sewtest != (double)-10000) {
                                    if(matrixml[0].Sewtest/(double)data[0].n_loci <
                                        matrixmlsim[it].Sewtest/(double)data[0].n_loci)
                                        nunderalfa[22] += 1;
                                    if(matrixml[0].Sewtest/(double)data[0].n_loci == matrixmlsim[it].Sewtest/(double)data[0].n_loci)
                                        nequivalfa[22] += 1;
                                    nitertests[22] += 1;
                                    if( data[0].n_loci > 2) {
                                        if(varobs[11] < variances((double)matrixmlsim[it].S2ewtest,(double)matrixmlsim[it].Sewtest,data[0].n_loci))
                                            nunderalfa[23] += 1;
                                        if(varobs[11] == variances((double)matrixmlsim[it].S2ewtest,(double)matrixmlsim[it].Sewtest,data[0].n_loci))
                                            nequivalfa[23] += 1;
                                        nitertests[23] += 1;
                                    }
                                }

                                if(matrixml[0].Shka != (double)-10000 && matrixmlsim[it].Shka != (double)-10000 && data[0].n_loci > 1) {
                                    if(matrixml[0].Shka < matrixmlsim[it].Shka)
                                        nunderalfa[24] += 1;
                                    if(matrixml[0].Shka == matrixmlsim[it].Shka)
                                        nequivalfa[24] += 1;
                                    nitertests[24] += 1;
                                }
                            }
                            /*print*/
                            if(nitertests[0]) {
                                prob = (double)nunderalfa[0]/(double)nitertests[0];
								if((double)1./((double)nitertests[0]) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[0]) *(double)nequivalfa[0];
                                else prob += (double)1./((double)nitertests[0]) *(double)nequivalfa[0];
                                printf("\n    P(average of Tajima's D = %g) = %g ",(double)matrixml[0].StajimaD/(double)matrixml[0].nltajD,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Tajima's D = %g) = %g ",(double)matrixml[0].StajimaD/ (double)matrixml[0].nltajD,prob);
                                if(nitertests[0] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Tajima's D) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Tajima's D) = na\t");
                            }
                            if(nitertests[1] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[1]/(double)nitertests[1];
								if((double)1./((double)nitertests[1]) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[1]) *(double)nequivalfa[1];
                                else prob += (double)1./((double)nitertests[1]) *(double)nequivalfa[1];
                                printf("\n    P(variance of Tajima's D = %g) = %g ",varobs[0],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Tajima's D = %g) = %g ",varobs[0],prob);
                                if(nitertests[1] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Tajima's D) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Tajima's D) = na\t");
                            }
                            if(nitertests[2]) {
                                prob = (double)nunderalfa[2]/(double)nitertests[2];
								if((double)1./((double)nitertests[2]) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[2]) *(double)nequivalfa[2];
                                else prob += (double)1./((double)nitertests[2]) *(double)nequivalfa[2];
                                printf("\n    P(average of Fu and Li's D = %g) = %g ",(double)matrixml[0].SfuliD/ (double)matrixml[0].nlflD,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Fu and Li's D = %g) = %g ",(double)matrixml[0].SfuliD/ (double)matrixml[0].nlflD,prob);
                                if(nitertests[2] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Fu and Li's D) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Fu and Li's D) = na\t");
                            }
                            if(nitertests[3] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[3]/(double)nitertests[3];
								if((double)1./((double)nitertests[3]) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[3]) *(double)nequivalfa[3];
                                else prob += (double)1./((double)nitertests[3]) *(double)nequivalfa[3];
                                printf("\n    P(variance of Fu and Li's D = %g) = %g ",varobs[1],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Fu and Li's D = %g) = %g ",varobs[1],prob);
                                if(data[0].n_loci > 2) if(nitertests[3] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Fu and Li's D) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Fu and Li's D) = na\t");
                            }
                            if(nitertests[4]) {
                                prob = (double)nunderalfa[4]/(double)nitertests[4];
								if((double)1./((double)nitertests[4]) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[4]) *(double)nequivalfa[4];
                                else prob += (double)1./((double)nitertests[4]) *(double)nequivalfa[4];
                                printf("\n    P(average of Fu and Li's F = %g) = %g ",(double)matrixml[0].SfuliF/ (double)matrixml[0].nlflF,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Fu and Li's F = %g) = %g ",(double)matrixml[0].SfuliF/ (double)matrixml[0].nlflF,prob);
                                if(nitertests[4] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Fu and Li's F) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Fu and Li's F) = na\t");
                            }
                            if(nitertests[5] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[5]/(double)nitertests[5];
								if((double)1./((double)nitertests[5]) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[5]) *(double)nequivalfa[5];
                                else prob += (double)1./((double)nitertests[5]) *(double)nequivalfa[5];
                                printf("\n    P(variance of Fu and Li's F = %g) = %g ",varobs[2],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Fu and Li's F = %g) = %g ",varobs[2],prob);
                                if(data[0].n_loci > 2) if(nitertests[5] >=20) printfignif(prob,file_output);
                             }
                            else{
                                printf("\n    P(variance of Fu and Li's F) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Fu and Li's F) = na\t");
                            }
                            if(nitertests[6]) {
								prob = (double)nunderalfa[6]/(double)nitertests[6];
								if((double)1./((double)nitertests[6]) *(double)nequivalfa[6] >= (double)0.5) nequivalfa[6] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[6]) *(double)nequivalfa[6];
                                else prob += (double)1./((double)nitertests[6]) *(double)nequivalfa[6];
                                printf("\n    P(average of Fu's Fs = %g) = %g ",(double)matrixml[0].SfuFs/ (double)matrixml[0].nlFs,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Fu's Fs = %g) = %g ",(double)matrixml[0].SfuFs/ (double)matrixml[0].nlFs,prob);
                                if(nitertests[6] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Fu's Fs) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Fu's Fs) = na\t");
                            }
                            if(nitertests[7] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[7]/(double)nitertests[7];
								if((double)1./((double)nitertests[7]) *(double)nequivalfa[7] >= (double)0.5) nequivalfa[7] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[7]) *(double)nequivalfa[7];
                                else prob += (double)1./((double)nitertests[7]) *(double)nequivalfa[7];
                                printf("\n    P(variance of Fu's Fs = %g) = %g ",varobs[3],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Fu's Fs = %g) = %g ",varobs[3],prob);
                                if(data[0].n_loci > 2) if(nitertests[7] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Fu's Fs) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Fu's Fs) = na\t");
                            }
                            if(nitertests[8]) {
                                prob = (double)nunderalfa[8]/(double)nitertests[8];
								if((double)1./((double)nitertests[8]) *(double)nequivalfa[8] >= (double)0.5) nequivalfa[8] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[8]) *(double)nequivalfa[8];
                                else prob += (double)1./((double)nitertests[8]) *(double)nequivalfa[8];
                                printf("\n    P(average of normalized Fay and Wu's H = %g) = %g ",(double)matrixml[0].SfaywuH/ (double)matrixml[0].nlH,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of normalized Fay and Wu's H = %g) = %g ",(double)matrixml[0].SfaywuH/ (double)matrixml[0].nlH,prob);
                                if(nitertests[8] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of normalized Fay and Wu's H) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of normalized Fay and Wu's H) = na\t");
                            }
                            if(nitertests[9] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[9]/(double)nitertests[9];
								if((double)1./((double)nitertests[9]) *(double)nequivalfa[9] >= (double)0.5) nequivalfa[9] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[9]) *(double)nequivalfa[9];
                                else prob += (double)1./((double)nitertests[9]) *(double)nequivalfa[9];
                                printf("\n    P(variance of normalized Fay and Wu's H = %g) = %g ",varobs[4],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of normalized Fay and Wu's H = %g) = %g ",varobs[4],prob);
                                if(data[0].n_loci > 2) if(nitertests[9] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of normalized Fay and Wu's H) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of normalized Fay and Wu's H) = na\t");
                            }
                            if(nitertests[10]) {
                                prob = (double)nunderalfa[10]/(double)nitertests[10];
								if((double)1./((double)nitertests[10]) *(double)nequivalfa[10] >= (double)0.5) nequivalfa[10] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[10]) *(double)nequivalfa[10];
                                else prob += (double)1./((double)nitertests[10]) *(double)nequivalfa[10];
                                printf("\n    P(average of Fay and Wu's H = %g) = %g ",(double)matrixml[0].SfaywuHo/ (double)matrixml[0].nlHo,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Fay and Wu's H = %g) = %g ",(double)matrixml[0].SfaywuHo/ (double)matrixml[0].nlHo,prob);
                                if(nitertests[10] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Fay and Wu's H) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Fay and Wu's H) = na\t");
                            }
                            if(nitertests[11] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[11]/(double)nitertests[11];
								if((double)1./((double)nitertests[11]) *(double)nequivalfa[11] >= (double)0.5) nequivalfa[11] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[11]) *(double)nequivalfa[11];
                                else prob += (double)1./((double)nitertests[11]) *(double)nequivalfa[11];
                                printf("\n    P(variance of Fay and Wu's H = %g) = %g ",varobs[5],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Fay and Wu's H = %g) = %g ",varobs[5],prob);
                                if(data[0].n_loci > 2) if(nitertests[11] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Fay and Wu's H) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Fay and Wu's H) = na\t");
                            }
                            if(nitertests[12]) {
                                prob = (double)nunderalfa[12]/(double)nitertests[12];
								if((double)1./((double)nitertests[12]) *(double)nequivalfa[12] >= (double)0.5) nequivalfa[12] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[12]) *(double)nequivalfa[12];
                                else prob += (double)1./((double)nitertests[12]) *(double)nequivalfa[12];
                                printf("\n    P(average of Rozas' et al. ZA = %g) = %g ",(double)matrixml[0].SrZA/ (double)matrixml[0].nlZ,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Rozas' et al. ZA = %g) = %g ",(double)matrixml[0].SrZA/ (double)matrixml[0].nlZ,prob);
                                if(nitertests[12] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Rozas' et al. ZA) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of rozas' et al. ZA) = na\t");
                            }
                            if(nitertests[13] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[13]/(double)nitertests[13];
								if((double)1./((double)nitertests[13]) *(double)nequivalfa[13] >= (double)0.5) nequivalfa[13] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[13]) *(double)nequivalfa[13];
                                else prob += (double)1./((double)nitertests[13]) *(double)nequivalfa[13];
                                printf("\n    P(variance of Rozas' et al. ZA = %g) = %g ",varobs[6],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Rozas' et al. ZA = %g) = %g ",varobs[6],prob);
                                if(data[0].n_loci > 2) if(nitertests[13] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Rozas' et al. ZA) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of rozas' et al. ZA) = na\t");
                            }
                            if(nitertests[14]) {
                                prob = (double)nunderalfa[14]/(double)nitertests[14];
								if((double)1./((double)nitertests[14]) *(double)nequivalfa[14] >= (double)0.5) nequivalfa[14] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[14]) *(double)nequivalfa[14];
                                else prob += (double)1./((double)nitertests[14]) *(double)nequivalfa[14];
                                printf("\n    P(average of Wall's B = %g) = %g ",(double)matrixml[0].SwB/ (double)matrixml[0].nlB,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Wall's B = %g) = %g ",(double)matrixml[0].SwB/ (double)matrixml[0].nlB,prob);
                                if(nitertests[14] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Wall's B) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Wall's B) = na\t");
                            }
                            if(nitertests[15] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[15]/(double)nitertests[15];
								if((double)1./((double)nitertests[15]) *(double)nequivalfa[15] >= (double)0.5) nequivalfa[15] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[15]) *(double)nequivalfa[15];
                                else prob += (double)1./((double)nitertests[15]) *(double)nequivalfa[15];
                                printf("\n    P(variance of Wall's B = %g) = %g ",varobs[7],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Wall's B = %g) = %g ",varobs[7],prob);
                                if(data[0].n_loci > 2) if(nitertests[15] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Wall's B) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Wall's B) = na\t");
                            }
                            if(nitertests[16]) {
                                prob = (double)nunderalfa[16]/(double)nitertests[16];
								if((double)1./((double)nitertests[16]) *(double)nequivalfa[16] >= (double)0.5) nequivalfa[16] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[16]) *(double)nequivalfa[16];
                                else prob += (double)1./((double)nitertests[16]) *(double)nequivalfa[16];
                                printf("\n    P(average of Wall's Q = %g) = %g ",(double)matrixml[0].SwQ/ (double)matrixml[0].nlQ,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Wall's Q = %g) = %g ",(double)matrixml[0].SwQ/ (double)matrixml[0].nlQ,prob);
                                if(nitertests[16] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Wall's Q) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Wall's Q) = na\t");
                            }
                            if(nitertests[17] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[17]/(double)nitertests[17];
								if((double)1./((double)nitertests[17]) *(double)nequivalfa[17] >= (double)0.5) nequivalfa[17] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[17]) *(double)nequivalfa[17];
                                else prob += (double)1./((double)nitertests[17]) *(double)nequivalfa[17];
                                printf("\n    P(variance of Wall's Q = %g) = %g ",varobs[8],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Wall's Q = %g) = %g ",varobs[8],prob);
                                if(data[0].n_loci > 2) if(nitertests[17] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Wall's Q) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Wall's Q) = na\t");
                            }
                            if(nitertests[18]) {
                                prob = (double)nunderalfa[18]/(double)nitertests[18];
								if((double)1./((double)nitertests[18]) *(double)nequivalfa[18] >= (double)0.5) nequivalfa[18] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[18]) *(double)nequivalfa[18];
                                else prob += (double)1./((double)nitertests[18]) *(double)nequivalfa[18];
                                printf("\n    P(average of R2 = %g) = %g ",(double)matrixml[0].SR2/ (double)matrixml[0].nlR2,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of R2 = %g) = %g ",(double)matrixml[0].SR2/ (double)matrixml[0].nlR2,prob);
                                if(nitertests[18] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Ramos and Rozas' R2) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Ramos and Rozas' R2) = na\t");
                            }
                            if(nitertests[19] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[19]/(double)nitertests[19];
								if((double)1./((double)nitertests[19]) *(double)nequivalfa[19] >= (double)0.5) nequivalfa[19] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[19]) *(double)nequivalfa[19];
                                else prob += (double)1./((double)nitertests[19]) *(double)nequivalfa[19];
                                printf("\n    P(variance of R2 = %g) = %g ",varobs[9],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of R2 = %g) = %g ",varobs[9],prob);
                                if(data[0].n_loci > 2) if(nitertests[19] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Ramos and Rozas' R2) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Ramos and Rozas' R2) = na\t");
                            }

                            if(nitertests[20]) {
                                prob = (double)nunderalfa[20]/(double)nitertests[20];
								if((double)1./((double)nitertests[20]) *(double)nequivalfa[20] >= (double)0.5) nequivalfa[20] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[20]) *(double)nequivalfa[20];
                                else prob += (double)1./((double)nitertests[20]) *(double)nequivalfa[20];
                                printf("\n    P(average of Zeng et al. E = %g) = %g ",(double)matrixml[0].SzengE/ (double)matrixml[0].nlE,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Zeng et al. E = %g) = %g ",(double)matrixml[0].SzengE/ (double)matrixml[0].nlE,prob);
                                if(nitertests[20] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Zeng et al. E) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Zeng et al. E) = na\t");
                            }
                            if(nitertests[21] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[21]/(double)nitertests[21];
								if((double)1./((double)nitertests[21]) *(double)nequivalfa[21] >= (double)0.5) nequivalfa[21] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[21]) *(double)nequivalfa[21];
                                else prob += (double)1./((double)nitertests[21]) *(double)nequivalfa[21];
                                printf("\n    P(variance of Zeng et al. E = %g) = %g ",varobs[10],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Zeng et al. E = %g) = %g ",varobs[10],prob);
                                if(data[0].n_loci > 2) if(nitertests[21] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Zeng et al. E) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Zeng et al. E) = na\t");
                            }
                            if(nitertests[22]) {
                                prob = (double)nunderalfa[22]/(double)nitertests[22];
								if((double)1./((double)nitertests[22]) *(double)nequivalfa[22] >= (double)0.5) nequivalfa[22] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[22]) *(double)nequivalfa[22];
                                else prob += (double)1./((double)nitertests[22]) *(double)nequivalfa[22];
                                printf("\n    P(average of Ewens-Watterson test = %g) = %g ",(double)matrixml[0].Sewtest/ (double)data[0].n_loci,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Ewens-Watterson test = %g) = %g ",(double)matrixml[0].Sewtest/ (double)data[0].n_loci,prob);
                                if(nitertests[20] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Ewens-Watterson test) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Ewens-Watterson test) = na\t");
                            }
                            if(nitertests[23] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[23]/(double)nitertests[23];
								if((double)1./((double)nitertests[23]) *(double)nequivalfa[23] >= (double)0.5) nequivalfa[23] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[23]) *(double)nequivalfa[23];
                                else prob += (double)1./((double)nitertests[23]) *(double)nequivalfa[23];
                                printf("\n    P(variance of Ewens-Watterson test = %g) = %g ",varobs[11],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Ewens-Watterson test = %g) = %g ",varobs[11],prob);
                                if(data[0].n_loci > 2) if(nitertests[23] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Ewens-Watterson test) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Ewens-Watterson test) = na\t");
                            }
                            if(nitertests[24] && data[0].n_loci > 1) {
                                prob = (double)nunderalfa[24]/(double)nitertests[24];
								/*if((double)1./((double)nitertests[24]) *(double)nequivalfa[24] >= (double)0.5) nequivalfa[24] /= 2;*/
                                /*
								if(prob > (double)0.5) prob -= (double)1./((double)nitertests[24]) *(double)nequivalfa[24];
                                else */prob += (double)1./((double)nitertests[24]) *(double)nequivalfa[24];
                                printf("\n    P(HKA test X = %g) = %g ",(double)matrixml[0].Shka,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(HKA test X = %g) = %g ",(double)matrixml[0].Shka,prob);
                                if(nitertests[24] >= 2) printfignif1(prob,file_output);
                            }
                            else{
                                printf("\n    P(HKA test X) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(HKA test X) = na\t");
                            }
                        }
                        else {/*not outgroup, observed_data*/
                            /*average and variance comparisons.*/
							printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                            printf("\n\n Display the probability of observed multilocus neutrality tests obtained by coalescent Monte Carlo simulations: \n");
                            printf("\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. \n");
                            printf("\n Multiple hits not included in the analysis.\n");
                            if(file_output) {
								if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                fputs("\n\n Display the probability of observed multilocus neutrality tests obtained by coalescent Monte Carlo simulations: \n",file_output);
                                fprintf(file_output,"\n E.g., Probabilities of 0.025 and 0.975 would be significant at 5%%. \n");
                                fputs("\n Multiple hits not included in the analysis.\n",file_output);
                            }
                            for(x=0;x<20;x++) nunderalfa[x] = 0;
                            for(x=0;x<20;x++) nequivalfa[x] = 0;
                            for(x=0;x<20;x++) nitertests[x] = 0;
                            
                            varobs[0] = variances((double)matrixml[0].S2tajimaD,(double)matrixml[0].StajimaD,matrixml[0].nltajD);
                            varobs[1] = variances((double)matrixml[0].S2fuliDn,(double)matrixml[0].SfuliDn,matrixml[0].nlflDn);
                            varobs[2] = variances((double)matrixml[0].S2fuliFn,(double)matrixml[0].SfuliFn,matrixml[0].nlflFn);
                            varobs[3] = variances((double)matrixml[0].S2fuFs,(double)matrixml[0].SfuFs,matrixml[0].nlFs);
                            varobs[4] = variances((double)matrixml[0].S2rZA,(double)matrixml[0].SrZA,matrixml[0].nlZ);
                            varobs[5] = variances((double)matrixml[0].S2wB,(double)matrixml[0].SwB,matrixml[0].nlB);
                            varobs[6] = variances((double)matrixml[0].S2wQ,(double)matrixml[0].SwQ,matrixml[0].nlQ);
                            varobs[7] = variances((double)matrixml[0].S2R2,(double)matrixml[0].SR2,matrixml[0].nlR2);
                            varobs[8] = variances((double)matrixml[0].S2ewtest,(double)matrixml[0].Sewtest,data[0].n_loci);
                            
                            for(it=0;it<(long int)data[0].n_iter;it++) {
                                if(matrixml[0].StajimaD != (double)-10000 && matrixmlsim[it].StajimaD != (double)-10000) {
                                    if(matrixml[0].StajimaD/(double)matrixml[0].nltajD <
                                        matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD)
                                        nunderalfa[0] += 1;
                                    if(matrixml[0].StajimaD/(double)matrixml[0].nltajD ==
                                        matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD)
                                        nequivalfa[0] += 1;
                                    nitertests[0] += 1;
                                    if(matrixmlsim[it].nltajD > 2) {
                                        if(varobs[0] < variances((double)matrixmlsim[it].S2tajimaD,(double)matrixmlsim[it].StajimaD,matrixmlsim[it].nltajD))
                                            nunderalfa[1] += 1;
                                        if(varobs[0] == variances((double)matrixmlsim[it].S2tajimaD,(double)matrixmlsim[it].StajimaD,matrixmlsim[it].nltajD))
                                            nequivalfa[1] += 1;
                                        nitertests[1] += 1;
                                    }
                                }
                                if(matrixml[0].SfuliDn != (double)-10000 && matrixmlsim[it].SfuliDn != (double)-10000) {
                                    if(matrixml[0].SfuliDn/(double)matrixml[0].nlflDn <
                                        matrixmlsim[it].SfuliDn/(double)matrixmlsim[it].nlflDn)
                                        nunderalfa[2] += 1;
                                    if(matrixml[0].SfuliDn/(double)matrixml[0].nlflDn ==
                                        matrixmlsim[it].SfuliDn/(double)matrixmlsim[it].nlflDn)
                                        nequivalfa[2] += 1;
                                    nitertests[2] += 1;
                                    if(matrixmlsim[it].nlflDn > 2) {
                                        if(varobs[1] < variances((double)matrixmlsim[it].S2fuliDn,(double)matrixmlsim[it].SfuliDn,matrixmlsim[it].nlflDn))
                                            nunderalfa[3] += 1;
                                        if(varobs[1] == variances((double)matrixmlsim[it].S2fuliDn,(double)matrixmlsim[it].SfuliDn,matrixmlsim[it].nlflDn))
                                            nequivalfa[3] += 1;
                                        nitertests[3] += 1;
                                    }
                                }
                                if(matrixml[0].SfuliFn != (double)-10000 && matrixmlsim[it].SfuliFn != (double)-10000) {
                                    if(matrixml[0].SfuliFn/(double)matrixml[0].nlflFn < 
                                        matrixmlsim[it].SfuliFn/(double)matrixmlsim[it].nlflFn)
                                        nunderalfa[4] += 1;
                                    if(matrixml[0].SfuliFn/(double)matrixml[0].nlflFn == 
                                        matrixmlsim[it].SfuliFn/(double)matrixmlsim[it].nlflFn)
                                        nequivalfa[4] += 1;
                                    nitertests[4] += 1;
                                    if(matrixmlsim[it].nlflFn > 2) {
                                        if(varobs[2] < variances((double)matrixmlsim[it].S2fuliFn,(double)matrixmlsim[it].SfuliFn,matrixmlsim[it].nlflFn))
                                            nunderalfa[5] += 1;
                                        if(varobs[2] == variances((double)matrixmlsim[it].S2fuliFn,(double)matrixmlsim[it].SfuliFn,matrixmlsim[it].nlflFn))
                                            nequivalfa[5] += 1;
                                        nitertests[5] += 1;
                                    }
                                }
                                if(matrixml[0].SfuFs != (double)-10000 && matrixmlsim[it].SfuFs != (double)-10000) {
                                    if(matrixml[0].SfuFs/(double)matrixml[0].nlFs < 
                                        matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs)
                                        nunderalfa[6] += 1;
                                    if(matrixml[0].SfuFs/(double)matrixml[0].nlFs == 
                                        matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs)
                                        nequivalfa[6] += 1;
                                    nitertests[6] += 1;
                                    if(matrixmlsim[it].nlFs > 2) {
                                        if(varobs[3] < variances((double)matrixmlsim[it].S2fuFs,(double)matrixmlsim[it].SfuFs,matrixmlsim[it].nlFs))
                                            nunderalfa[7] += 1;
                                        if(varobs[3] == variances((double)matrixmlsim[it].S2fuFs,(double)matrixmlsim[it].SfuFs,matrixmlsim[it].nlFs))
                                            nequivalfa[7] += 1;
                                        nitertests[7] += 1;
                                    }
                                }
                                if(matrixml[0].SrZA != (double)-10000 && matrixmlsim[it].SrZA != (double)-10000) {
                                    if(matrixml[0].SrZA/(double)matrixml[0].nlZ < matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ)
                                        nunderalfa[8] += 1;
                                    if(matrixml[0].SrZA/(double)matrixml[0].nlZ == matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ)
                                        nequivalfa[8] += 1;
                                    nitertests[8] += 1;
                                    if(matrixmlsim[it].nlZ > 2) {
                                        if(varobs[4] < variances((double)matrixmlsim[it].S2rZA,(double)matrixmlsim[it].SrZA,matrixmlsim[it].nlZ))
                                            nunderalfa[9] += 1;
                                        if(varobs[4] == variances((double)matrixmlsim[it].S2rZA,(double)matrixmlsim[it].SrZA,matrixmlsim[it].nlZ))
                                            nequivalfa[9] += 1;
                                        nitertests[9] += 1;
                                    }
                                }
                                if(matrixml[0].SwB != (double)-10000 && matrixmlsim[it].SwB != (double)-10000) {
                                    if(matrixml[0].SwB/(double)matrixml[0].nlB < matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB)
                                        nunderalfa[10] += 1;
                                    if(matrixml[0].SwB/(double)matrixml[0].nlB == matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB)
                                        nequivalfa[10] += 1;
                                    nitertests[10] += 1;
                                    if(matrixmlsim[it].nlB > 2) {
                                        if(varobs[5] < variances((double)matrixmlsim[it].S2wB,(double)matrixmlsim[it].SwB,matrixmlsim[it].nlB))
                                            nunderalfa[11] += 1;
                                        if(varobs[5] == variances((double)matrixmlsim[it].S2wB,(double)matrixmlsim[it].SwB,matrixmlsim[it].nlB))
                                            nequivalfa[11] += 1;
                                        nitertests[11] += 1;
                                    }
                                }
                                if(matrixml[0].SwQ != (double)-10000 && matrixmlsim[it].SwQ != (double)-10000) {
                                    if(matrixml[0].SwQ/(double)matrixml[0].nlQ < matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ)
                                       nunderalfa[12] += 1;
                                    if(matrixml[0].SwQ/(double)matrixml[0].nlQ == matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ)
                                        nequivalfa[12] += 1;
                                    nitertests[12] += 1;
                                    if(matrixmlsim[it].nlQ > 2) {
                                        if(varobs[6] < variances((double)matrixmlsim[it].S2wQ,(double)matrixmlsim[it].SwQ,matrixmlsim[it].nlQ))
                                            nunderalfa[13] += 1;
                                        if(varobs[6] == variances((double)matrixmlsim[it].S2wQ,(double)matrixmlsim[it].SwQ,matrixmlsim[it].nlQ))
                                            nequivalfa[13] += 1;
                                        nitertests[13] += 1;
                                    }
                                }
                                if(matrixml[0].SR2 != (double)-10000 && matrixmlsim[it].SR2 != (double)-10000) {
                                    if(matrixml[0].SR2/(double)matrixml[0].nlR2 <
                                        matrixmlsim[it].SR2/matrixmlsim[it].nlR2)
                                        nunderalfa[14] += 1;
                                    if(matrixml[0].SR2/(double)matrixml[0].nlR2 ==
                                        matrixmlsim[it].SR2/matrixmlsim[it].nlR2)
                                        nequivalfa[14] += 1;
                                    nitertests[14] += 1;
                                    if( matrixmlsim[it].nlR2 > 2) {
                                        if(varobs[7] < variances((double)matrixmlsim[it].S2R2,(double)matrixmlsim[it].SR2,matrixmlsim[it].nlR2))
                                            nunderalfa[15] += 1;
                                        if(varobs[7] == variances((double)matrixmlsim[it].S2R2,(double)matrixmlsim[it].SR2,matrixmlsim[it].nlR2))
                                            nequivalfa[15] += 1;
                                        nitertests[15] += 1;
                                    }
                                }
                                if(matrixml[0].Sewtest != (double)-10000 && matrixmlsim[it].Sewtest != (double)-10000) {
                                    if(matrixml[0].Sewtest/(double)data[0].n_loci <
                                        matrixmlsim[it].Sewtest/(double)data[0].n_loci)
                                        nunderalfa[16] += 1;
                                    if(matrixml[0].Sewtest/(double)data[0].n_loci ==
                                        matrixmlsim[it].Sewtest/(double)data[0].n_loci)
                                        nequivalfa[16] += 1;
                                    nitertests[16] += 1;
                                    if( data[0].n_loci > 2) {
                                        if(varobs[8] < variances((double)matrixmlsim[it].S2ewtest,(double)matrixmlsim[it].Sewtest,data[0].n_loci))
                                            nunderalfa[17] += 1;
                                        if(varobs[8] == variances((double)matrixmlsim[it].S2ewtest,(double)matrixmlsim[it].Sewtest,data[0].n_loci))
                                            nequivalfa[17] += 1;
                                        nitertests[17] += 1;
                                    }
                                }
                            }
                            /*print*/
                            if(nitertests[0]) {
                                prob = (double)nunderalfa[0]/(double)nitertests[0];
								if((double)1./((double)nitertests[0]) *(double)nequivalfa[0] >= (double)0.5) nequivalfa[0] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[0]) *(double)nequivalfa[0];
                                else prob += (double)1./((double)nitertests[0]) *(double)nequivalfa[0];
                                printf("\n    P(average of Tajima's D = %g) = %g ",(double)matrixml[0].StajimaD/(double)matrixml[0].nltajD,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Tajima's D = %g) = %g ",(double)matrixml[0].StajimaD/ (double)matrixml[0].nltajD,prob);
                                if(nitertests[0] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Tajima's D) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Tajima's D) = na\t");
                            }
                            if(nitertests[1] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[1]/(double)nitertests[1];
								if((double)1./((double)nitertests[1]) *(double)nequivalfa[1] >= (double)0.5) nequivalfa[1] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[1]) *(double)nequivalfa[1];
                                else prob += (double)1./((double)nitertests[1]) *(double)nequivalfa[1];
                                printf("\n    P(variance of Tajima's D = %g) = %g ",varobs[0],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Tajima's D = %g) = %g ",varobs[0],prob);
                                if(data[0].n_loci > 2) if(nitertests[1] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Tajima's D) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Tajima's D) = na\t");
                            }
                            if(nitertests[2]) {
                                prob = (double)nunderalfa[2]/(double)nitertests[2];
								if((double)1./((double)nitertests[2]) *(double)nequivalfa[2] >= (double)0.5) nequivalfa[2] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[2]) *(double)nequivalfa[2];
                                else prob += (double)1./((double)nitertests[2]) *(double)nequivalfa[2];
                                printf("\n    P(average of Fu and Li's D* = %g) = %g ",(double)matrixml[0].SfuliDn/ (double)matrixml[0].nlflDn,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Fu and Li's D* = %g) = %g ",(double)matrixml[0].SfuliDn/ (double)matrixml[0].nlflDn,prob);
                                if(nitertests[2] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Fu and Li's D*) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Fu and Li's D*) = na\t");
                            }
                            if(nitertests[3] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[3]/(double)nitertests[3];
								if((double)1./((double)nitertests[3]) *(double)nequivalfa[3] >= (double)0.5) nequivalfa[3] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[3]) *(double)nequivalfa[3];
                                else prob += (double)1./((double)nitertests[3]) *(double)nequivalfa[3];
                                printf("\n    P(variance of Fu and Li's D* = %g) = %g ",varobs[1],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Fu and Li's D* = %g) = %g ",varobs[1],prob);
                                if(data[0].n_loci > 2) if(nitertests[3] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Fu and Li's D*) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Fu and Li's D*) = na\t");
                            }
                            if(nitertests[4]) {
                                prob = (double)nunderalfa[4]/(double)nitertests[4];
								if((double)1./((double)nitertests[4]) *(double)nequivalfa[4] >= (double)0.5) nequivalfa[4] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[4]) *(double)nequivalfa[4];
                                else prob += (double)1./((double)nitertests[4]) *(double)nequivalfa[4];
                                printf("\n    P(average of Fu and Li's F* = %g) = %g ",(double)matrixml[0].SfuliFn/ (double)matrixml[0].nlflFn,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Fu and Li's F* = %g) = %g ",(double)matrixml[0].SfuliFn/ (double)matrixml[0].nlflFn,prob);
                                if(nitertests[4] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Fu and Li's F*) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Fu and Li's F*) = na\t");
                            }
                            if(nitertests[5] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[5]/(double)nitertests[5];
								if((double)1./((double)nitertests[5]) *(double)nequivalfa[5] >= (double)0.5) nequivalfa[5] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[5]) *(double)nequivalfa[5];
                                else prob += (double)1./((double)nitertests[5]) *(double)nequivalfa[5];
                                printf("\n    P(variance of Fu and Li's F* = %g) = %g ",varobs[2],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Fu and Li's F* = %g) = %g ",varobs[2],prob);
                                if(data[0].n_loci > 2) if(nitertests[5] >=20) printfignif(prob,file_output);
                             }
                            else{
                                printf("\n    P(variance of Fu and Li's F*) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Fu and Li's F*) = na\t");
                            }
                            if(nitertests[6]) {
								prob = (double)nunderalfa[6]/(double)nitertests[6];
								if((double)1./((double)nitertests[6]) *(double)nequivalfa[6] >= (double)0.5) nequivalfa[6] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[6]) *(double)nequivalfa[6];
                                else prob += (double)1./((double)nitertests[6]) *(double)nequivalfa[6];
                                printf("\n    P(average of Fu's Fs = %g) = %g ",(double)matrixml[0].SfuFs/ (double)matrixml[0].nlFs,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Fu's Fs = %g) = %g ",(double)matrixml[0].SfuFs/ (double)matrixml[0].nlFs,prob);
                                if(nitertests[6] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Fu's Fs) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Fu's Fs) = na\t");
                            }
                            if(nitertests[7] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[7]/(double)nitertests[7];
								if((double)1./((double)nitertests[7]) *(double)nequivalfa[7] >= (double)0.5) nequivalfa[7] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[7]) *(double)nequivalfa[7];
                                else prob += (double)1./((double)nitertests[7]) *(double)nequivalfa[7];
                                printf("\n    P(variance of Fu's Fs = %g) = %g ",varobs[3],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Fu's Fs = %g) = %g ",varobs[3],prob);
                                if(data[0].n_loci > 2) if(nitertests[7] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Fu's Fs) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Fu's Fs) = na\t");
                            }
                            if(nitertests[8]) {
                                prob = (double)nunderalfa[8]/(double)nitertests[8];
								if((double)1./((double)nitertests[8]) *(double)nequivalfa[8] >= (double)0.5) nequivalfa[8] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[8]) *(double)nequivalfa[8];
                                else prob += (double)1./((double)nitertests[8]) *(double)nequivalfa[8];
                                printf("\n    P(average of Rozas' et al. ZA = %g) = %g ",(double)matrixml[0].SrZA/ (double)matrixml[0].nlZ,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Rozas' et al. ZA = %g) = %g ",(double)matrixml[0].SrZA/ (double)matrixml[0].nlZ,prob);
                                if(nitertests[8] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Rozas' et al. ZA) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of rozas' et al. ZA) = na\t");
                            }
                            if(nitertests[9] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[9]/(double)nitertests[9];
								if((double)1./((double)nitertests[9]) *(double)nequivalfa[9] >= (double)0.5) nequivalfa[9] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[9]) *(double)nequivalfa[9];
                                else prob += (double)1./((double)nitertests[9]) *(double)nequivalfa[9];
                                printf("\n    P(variance of Rozas' et al. ZA = %g) = %g ",varobs[4],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Rozas' et al. ZA = %g) = %g ",varobs[4],prob);
                                if(data[0].n_loci > 2) if(nitertests[9] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Rozas' et al. ZA) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of rozas' et al. ZA) = na\t");
                            }
                            if(nitertests[10]) {
                                prob = (double)nunderalfa[10]/(double)nitertests[10];
								if((double)1./((double)nitertests[10]) *(double)nequivalfa[10] >= (double)0.5) nequivalfa[10] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[10]) *(double)nequivalfa[10];
                                else prob += (double)1./((double)nitertests[10]) *(double)nequivalfa[10];
                                printf("\n    P(average of Wall's B = %g) = %g ",(double)matrixml[0].SwB/ (double)matrixml[0].nlB,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Wall's B = %g) = %g ",(double)matrixml[0].SwB/ (double)matrixml[0].nlB,prob);
                                if(nitertests[10] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Wall's B) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Wall's B) = na\t");
                            }
                            if(nitertests[11] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[11]/(double)nitertests[11];
								if((double)1./((double)nitertests[11]) *(double)nequivalfa[11] >= (double)0.5) nequivalfa[11] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[11]) *(double)nequivalfa[11];
                                else prob += (double)1./((double)nitertests[11]) *(double)nequivalfa[11];
                                printf("\n    P(variance of Wall's B = %g) = %g ",varobs[5],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Wall's B = %g) = %g ",varobs[5],prob);
                                if(data[0].n_loci > 2) if(nitertests[11] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Wall's B) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Wall's B) = na\t");
                            }
                            if(nitertests[12]) {
                                prob = (double)nunderalfa[12]/(double)nitertests[12];
								if((double)1./((double)nitertests[12]) *(double)nequivalfa[12] >= (double)0.5) nequivalfa[12] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[12]) *(double)nequivalfa[12];
                                else prob += (double)1./((double)nitertests[12]) *(double)nequivalfa[12];
                                printf("\n    P(average of Wall's Q = %g) = %g ",(double)matrixml[0].SwQ/ (double)matrixml[0].nlQ,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Wall's Q = %g) = %g ",(double)matrixml[0].SwQ/ (double)matrixml[0].nlQ,prob);
                                if(nitertests[12] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Wall's Q) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Wall's Q) = na\t");
                            }
                            if(nitertests[13] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[13]/(double)nitertests[13];
								if((double)1./((double)nitertests[13]) *(double)nequivalfa[13] >= (double)0.5) nequivalfa[13] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[13]) *(double)nequivalfa[13];
                                else prob += (double)1./((double)nitertests[13]) *(double)nequivalfa[13];
                                printf("\n    P(variance of Wall's Q = %g) = %g ",varobs[6],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Wall's Q = %g) = %g ",varobs[6],prob);
                                if(data[0].n_loci > 2) if(nitertests[13] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Wall's Q) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Wall's Q) = na\t");
                            }
                            if(nitertests[14]) {
                                prob = (double)nunderalfa[14]/(double)nitertests[14];
								if((double)1./((double)nitertests[14]) *(double)nequivalfa[14] >= (double)0.5) nequivalfa[14] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[14]) *(double)nequivalfa[14];
                                else prob += (double)1./((double)nitertests[14]) *(double)nequivalfa[14];
                                printf("\n    P(average of R2 = %g) = %g ",(double)matrixml[0].SR2/ (double)matrixml[0].nlR2,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of R2 = %g) = %g ",(double)matrixml[0].SR2/ (double)matrixml[0].nlR2,prob);
                                if(nitertests[14] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average of Ramos and Rozas' R2) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Ramos and Rozas' R2) = na\t");
                            }
                            if(nitertests[15] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[15]/(double)nitertests[15];
								if((double)1./((double)nitertests[15]) *(double)nequivalfa[15] >= (double)0.5) nequivalfa[15] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[15]) *(double)nequivalfa[15];
                                else prob += (double)1./((double)nitertests[15]) *(double)nequivalfa[15];
                                printf("\n    P(variance of R2 = %g) = %g ",varobs[7],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of R2 = %g) = %g ",varobs[7],prob);
                                if(data[0].n_loci > 2) if(nitertests[15] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Ramos and Rozas' R2) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Ramos and Rozas' R2) = na\t");
                            }
                            if(nitertests[16]) {
                                prob = (double)nunderalfa[16]/(double)nitertests[16];
								if((double)1./((double)nitertests[16]) *(double)nequivalfa[16] >= (double)0.5) nequivalfa[16] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[16]) *(double)nequivalfa[16];
                                else prob += (double)1./((double)nitertests[16]) *(double)nequivalfa[16];
                                printf("\n    P(average of Ewens-Watterson test = %g) = %g ",(double)matrixml[0].Sewtest/ (double)data[0].n_loci,prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(average of Ewens-Watterson test = %g) = %g ",(double)matrixml[0].Sewtest/ (double)data[0].n_loci,prob);
                                if(nitertests[16] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(average ofEwens-Watterson test) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(average of Ewens-Watterson test) = na\t");
                            }
                            if(nitertests[17] && data[0].n_loci > 2) {
                                prob = (double)nunderalfa[17]/(double)nitertests[17];
								if((double)1./((double)nitertests[17]) *(double)nequivalfa[17] >= (double)0.5) nequivalfa[17] /= 2;
                                if(prob > (double)0.5) prob -= (double)1./((double)nitertests[17]) *(double)nequivalfa[17];
                                else prob += (double)1./((double)nitertests[17]) *(double)nequivalfa[17];
                                printf("\n    P(variance of Ewens-Watterson test = %g) = %g ",varobs[8],prob);
                                if(file_output)
                                    fprintf(file_output,"\n    P(variance of Ewens-Watterson test = %g) = %g ",varobs[8],prob);
                                if(data[0].n_loci > 2) if(nitertests[17] >=20) printfignif(prob,file_output);
                            }
                            else{
                                printf("\n    P(variance of Ewens-Watterson test) = na\t");
                                if(file_output) fprintf(file_output,"\n    P(variance of Ewens-Watterson test) = na\t");
                            }
                        }
                    }
					#if ALLRESULTS == 0
                    else {/*give confident intervals with simulations. NOT observed_data.*/
					#endif
                        if((sortvector = (double *) malloc(data[0].n_iter * sizeof(double))) == 0) {
                            printf("Probabilities can not be calculated, sorry.\n");
                            if(file_output) fputs("Probabilities can not be calculated, sorry.\n",file_output);
                            break;
                        }
                        if(*outgroup) {
							printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                            printf("\nTable containing the name of neutrality tests, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n");
                            printf("statistic\tmedian\tavg\tvar\t");                                    
                            if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                            if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                            if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                            if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                            if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                            if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                            if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                            printf("\n");
                            if(file_output) {
								if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                fputs("\nTable containing the name of neutrality tests, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n",file_output);
                                fputs("statistic\tmedian\tavg\tvar\t",file_output);                                    
                                if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                fprintf(file_output,"\n");
                            }
                            
                            /*StajimaD. */
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if(matrixmlsim[it].StajimaD != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD;
                                    var += (double)matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD * (double)matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_TajimaD\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_TajimaD\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_Tajima's_D\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_Tajima's_D\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2tajimaD,matrixmlsim[it].StajimaD,matrixmlsim[it].nltajD)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_TajimaD\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_TajimaD\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_Tajima's_D\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_Tajima's_D\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SfuliD*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SfuliD != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SfuliD/(double)matrixmlsim[it].nlflDn;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SfuliD/(double)matrixmlsim[it].nlflDn;
                                    var += (double)matrixmlsim[it].SfuliD/(double)matrixmlsim[it].nlflDn * (double)matrixmlsim[it].SfuliD/(double)matrixmlsim[it].nlflDn;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_FuLiD\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_FuLiD\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_FuLiD\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_FuLiD\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2fuliD,matrixmlsim[it].SfuliD,matrixmlsim[it].nlflD)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_FuLiD\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_FuLiD\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_FuLiD\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_FuLiD\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SfuliF*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SfuliF != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SfuliF/(double)matrixmlsim[it].nlflFn;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SfuliF/(double)matrixmlsim[it].nlflFn;
                                    var += (double)matrixmlsim[it].SfuliF/(double)matrixmlsim[it].nlflFn * (double)matrixmlsim[it].SfuliF/(double)matrixmlsim[it].nlflFn;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_FuLiF\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_FuLiF\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_FuLiF\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_FuLiF\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2fuliF,matrixmlsim[it].SfuliF,matrixmlsim[it].nlflF)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_FuLiF\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_FuLiF\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_FuLiF\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_FuLiF\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SfuFs*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SfuFs != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs;
                                    var += (double)matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs * (double)matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_FuFs\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_FuFs\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_FuFs\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_FuFs\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2fuFs,matrixmlsim[it].SfuFs,matrixmlsim[it].nlFs)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_FuFs\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_FuFs\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_FuFs\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_FuFs\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SfaywuH*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SfaywuH != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SfaywuH/(double)matrixmlsim[it].nlH;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SfaywuH/(double)matrixmlsim[it].nlH;
                                    var += (double)matrixmlsim[it].SfaywuH/(double)matrixmlsim[it].nlH * (double)matrixmlsim[it].SfaywuH/(double)matrixmlsim[it].nlH;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_FayWuHn\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_FaynWuH\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_FayWuHn\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_FayWuHn\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2faywuH,matrixmlsim[it].SfaywuH,matrixmlsim[it].nlH)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_FayWuHn\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_FayWuHn\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_FayWuHn\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_FayWuHn\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SfaywuHo*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SfaywuHo != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SfaywuHo/(double)data[0].n_loci;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SfaywuHo/(double)data[0].n_loci;
                                    var += (double)matrixmlsim[it].SfaywuHo/(double)data[0].n_loci * (double)matrixmlsim[it].SfaywuH/(double)data[0].n_loci;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_FayWuH\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_FayWuH\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_FayWuH\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_FayWuH\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2faywuHo,matrixmlsim[it].SfaywuHo,matrixmlsim[it].nlHo)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_FayWuH\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_FayWuH\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_FayWuH\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_FayWuH\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SrZA*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SrZA != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ;
                                    var += (double)matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ * (double)matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_RozasZA\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_RozasZA\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_RozasZA\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_RozasZA\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2rZA,matrixmlsim[it].SrZA,matrixmlsim[it].nlZ)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_RozasZA\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_RozasZA\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_RozasZA\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_RozasZA\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SwB*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SwB != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB;
                                    var += (double)matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB * (double)matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB;
                                }
                            }
                            if(it2 >0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_WallB\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_WallB\t%g\t%g\t",
                                        sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_WallB\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_WallB\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2wB,
                                        matrixmlsim[it].SwB,matrixmlsim[it].nlB)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_WallB\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_WallB\t%g\t%g\t",
                                            sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_WallB\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_WallB\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SwQ*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SwQ != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ;
                                    var += (double)matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ * (double)matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_WallQ\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_WallQ\t%g\t%g\t",
                                        sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_WallQ\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_WallQ\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2wQ,matrixmlsim[it].SwQ,matrixmlsim[it].nlQ)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_WallQ\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_WallQ\t%g\t%g\t",
                                            sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_WallQ\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_WallQ\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SR2*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SR2 != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SR2/(double)matrixmlsim[it].nlR2;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SR2/(double)matrixmlsim[it].nlR2;
                                    var += (double)matrixmlsim[it].SR2/(double)matrixmlsim[it].nlR2 * (double)matrixmlsim[it].SR2/(double)matrixmlsim[it].nlR2;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_R2\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_R2\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_R2\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_R2\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2R2,matrixmlsim[it].SR2,matrixmlsim[it].nlR2)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_R2\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_R2\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_R2\tna\tna\tna\n");
                                    if(file_output) 
                                        fprintf(file_output,"var_R2\tna\tna\tna\n");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
							}
                            /*SzengE*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SzengE != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SzengE/(double)matrixmlsim[it].nlE;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SzengE/(double)matrixmlsim[it].nlE;
                                    var += (double)matrixmlsim[it].SzengE/(double)matrixmlsim[it].nlE * (double)matrixmlsim[it].SzengE/(double)matrixmlsim[it].nlE;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_zengE\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_zengE\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_zengE\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_zengE\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2zengE,matrixmlsim[it].SzengE,matrixmlsim[it].nlE)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_ZengE\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_ZengE\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_zengE\tna\tna\tna\n");
                                    if(file_output) 
                                        fprintf(file_output,"var_zengE\tna\tna\tna\n");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SEW*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].Sewtest != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].Sewtest/(double)data[0].n_loci;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].Sewtest/(double)data[0].n_loci;
                                    var += (double)matrixmlsim[it].Sewtest/(double)data[0].n_loci * (double)matrixmlsim[it].Sewtest/(double)data[0].n_loci;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_EW\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_EW\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_EW\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_EW\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2ewtest,matrixmlsim[it].Sewtest,data[0].n_loci)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_EW\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_EW\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_EW\tna\tna\tna\n");
                                    if(file_output) 
                                        fprintf(file_output,"var_EW\tna\tna\tna\n");
                                    print_percentages(data,sortvector,it2,file_output);
                                }

								/*HKA*/
								/*Shka*/
								avg = var = (double)0.0;
								for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
									if(matrixmlsim[it].Shka != (double) -10000) {
										sortvector[it2] = (double)matrixmlsim[it].Shka;
										it2 += 1;
										avg += (double)matrixmlsim[it].Shka;
										var += (double)matrixmlsim[it].Shka * (double)matrixmlsim[it].Shka;
									}
								}
								if(it2 > 0) {
									var = varianceit(var,avg,it2);
									if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
									printf("HKA\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
									if(file_output)
										fprintf(file_output,"HKA\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
									if(var != (double)-10000) {
										printf("%g\t",var);
										if(file_output) fprintf(file_output,"%g\t",var);
									}
									else {
										printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
									print_percentages(data,sortvector,it2,file_output);
								}
								else {
									printf("HKA\tna\tna\tna\t");
									if(file_output) 
										fprintf(file_output,"HKA\tna\tna\tna\t");
									print_percentages(data,sortvector,it2,file_output);
								}
								/*hka_T*/
								avg = var = (double)0.0;
								for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
									if(matrixmlsim[it].hka_T != (double) -10000) {
										sortvector[it2] = (double)matrixmlsim[it].hka_T;
										it2 += 1;
										avg += (double)matrixmlsim[it].hka_T;
										var += (double)matrixmlsim[it].hka_T * (double)matrixmlsim[it].hka_T;
									}
								}
								if(it2 > 0) {
									var = varianceit(var,avg,it2);
									if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
									printf("t(HKA)\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
									if(file_output)
										fprintf(file_output,"t(HKA)\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
									if(var != (double)-10000) {
										printf("%g\t",var);
										if(file_output) fprintf(file_output,"%g\t",var);
									}
									else {
										printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
									print_percentages(data,sortvector,it2,file_output);
								}
								else {
									printf("t(HKA)\tna\tna\tna\t");
									if(file_output) 
										fprintf(file_output,"t(HKA)\tna\tna\tna\t");
									print_percentages(data,sortvector,it2,file_output);
								}
							}
                        }
                        else {
                            /*Do a NEW vector with all the values, SORT.*/
                            /* Take the median (50%), average and variance, and 10%, 5%, 2.5%, 1%, 0.5%, 0.1% and 0.05% for 1 tails.*/
							printf("\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                            printf("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n");
                            printf("statistic\tmedian\tavg\tvar\t");                                    
                            if(data[0].n_iter >= 10) printf("-10%%\t+10%%\t");
                            if(data[0].n_iter >= 20) printf("-5%%\t+5%%\t");
                            if(data[0].n_iter >= 40) printf("-2.5%%\t+2.5%%\t");
                            if(data[0].n_iter >= 100) printf("-1%%\t+1%%\t");
                            if(data[0].n_iter >= 200) printf("-0.5%%\t+0.5%%\t");
                            if(data[0].n_iter >= 1000) printf("-0.1%%\t+0.1%%\t");
                            if(data[0].n_iter >= 2000) printf("-0.05%%\t+0.005%%\t");
                            printf("\n");
                            if(file_output) {
								if(file_output) fprintf(file_output,"\nAnalysis based on %ld iterations for each locus.\n\n",data[0].n_iter);
                                fputs("\nTable containing the name of statistic, median, average and variance, and values at different probabilities for one tail \n (for two tails the value of the probability is doubled).\n\n",file_output);
                                fputs("statistic\tmedian\tavg\tvar\t",file_output);                                    
                                if(data[0].n_iter >= 10) fputs("-10%\t+10%\t",file_output);
                                if(data[0].n_iter >= 20) fputs("-5%\t+5%\t",file_output);
                                if(data[0].n_iter >= 40) fputs("-2.5%\t+2.5%\t",file_output);
                                if(data[0].n_iter >= 100) fputs("-1%\t+1%\t",file_output);
                                if(data[0].n_iter >= 200) fputs("-0.5%\t+0.5%\t",file_output);
                                if(data[0].n_iter >= 1000) fputs("-0.1%\t+0.1%\t",file_output);
                                if(data[0].n_iter >= 2000) fputs("-0.05%\t+0.005%\t",file_output);
                                fprintf(file_output,"\n");
                            }
                                
                            /*StajimaD. */
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].StajimaD != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD;
                                    var += (double)matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD * (double)matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_TajimaD\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_TajimaD\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");
								}
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_Tajima's_D\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_Tajima's_D\tna\tna\tna\t");
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2tajimaD,matrixmlsim[it].StajimaD,matrixmlsim[it].nltajD)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_TajimaD\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_TajimaD\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_Tajima's_D\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_Tajima's_D\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SfuliD**/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SfuliDn != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SfuliDn/(double)matrixmlsim[it].nlflDn;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SfuliDn/(double)matrixmlsim[it].nlflDn;
                                    var += (double)matrixmlsim[it].SfuliDn/(double)matrixmlsim[it].nlflDn * (double)matrixmlsim[it].SfuliDn/(double)matrixmlsim[it].nlflDn;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_FuLiD*\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_FuLiD*\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_FuLiD*\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_FuLiD*\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2fuliDn,matrixmlsim[it].SfuliDn,matrixmlsim[it].nlflDn)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_FuLiD*\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_FuLiD*\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_FuLiD*\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_FuLiD*\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SfuliF*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SfuliFn != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SfuliFn/(double)matrixmlsim[it].nlflFn;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SfuliFn/(double)matrixmlsim[it].nlflFn;
                                    var += (double)matrixmlsim[it].SfuliFn/(double)matrixmlsim[it].nlflFn * (double)matrixmlsim[it].SfuliFn/(double)matrixmlsim[it].nlflFn;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_FuLiF*\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_FuLiF*\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_FuLiF*\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_FuLiF*\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2fuliFn,matrixmlsim[it].SfuliFn,matrixmlsim[it].nlflFn)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_FuLiF*\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_FuLiF*\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_FuLiF*\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_FuLiF*\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SfuFs*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SfuFs != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs;
                                    var += (double)matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs * (double)matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_FuFs\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_FuFs\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_FuFs\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_FuFs\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2fuFs,matrixmlsim[it].SfuFs,matrixmlsim[it].nlFs)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_FuFs\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_FuFs\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_FuFs\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_FuFs\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SrZA*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SrZA != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ;
                                    var += (double)matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ * (double)matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_RozasZA\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_RozasZA\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_RozasZA\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_RozasZA\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2rZA,matrixmlsim[it].SrZA,matrixmlsim[it].nlZ)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_RozasZA\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_RozasZA\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_RozasZA\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_RozasZA\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SwB*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SwB != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB;
                                    var += (double)matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB * (double)matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB;
                                }
                            }
                            if(it2 >0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_WallB\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_WallB\t%g\t%g\t",
                                        sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_WallB\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_WallB\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2wB,
                                        matrixmlsim[it].SwB,matrixmlsim[it].nlB)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_WallB\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_WallB\t%g\t%g\t",
                                            sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_WallB\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_WallB\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SwQ*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SwQ != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ;
                                    var += (double)matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ * (double)matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_WallQ\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_WallQ\t%g\t%g\t",
                                        sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_WallQ\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_WallQ\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2wQ,matrixmlsim[it].SwQ,
                                        matrixmlsim[it].nlQ)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_WallQ\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_WallQ\t%g\t%g\t",
                                            sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");                                }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_WallQ\tna\tna\tna\t");
                                    if(file_output) 
                                        fprintf(file_output,"var_WallQ\tna\tna\tna\t");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SR2*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].SR2 != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].SR2/(double)matrixmlsim[it].nlR2;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].SR2/(double)matrixmlsim[it].nlR2;
                                    var += (double)matrixmlsim[it].SR2/(double)matrixmlsim[it].nlR2 * (double)matrixmlsim[it].SR2/(double)matrixmlsim[it].nlR2;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_R2\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_R2\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_R2\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_R2\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2R2,matrixmlsim[it].SR2,matrixmlsim[it].nlR2)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_R2\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_R2\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");
                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_R2\tna\tna\tna\n");
                                    if(file_output) 
                                        fprintf(file_output,"var_R2\tna\tna\tna\n");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                            /*SEW*/
                            /*average*/
                            avg = var = (double)0.0;
                            for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                if( matrixmlsim[it].Sewtest != (double) -10000) {
                                    sortvector[it2] = (double)matrixmlsim[it].Sewtest/(double)data[0].n_loci;
                                    it2 += 1;
                                    avg += (double)matrixmlsim[it].Sewtest/(double)data[0].n_loci;
                                    var += (double)matrixmlsim[it].Sewtest/(double)data[0].n_loci * (double)matrixmlsim[it].Sewtest/(double)data[0].n_loci;
                                }
                            }
                            if(it2 > 0) {
                                var = varianceit(var,avg,it2);
                                if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                printf("avg_EW\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(file_output)
                                    fprintf(file_output,"avg_EW\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                if(var != (double)-10000) {
                                    printf("%g\t",var);
                                    if(file_output) fprintf(file_output,"%g\t",var);
                                }
                                else {
                                    printf("na\t");
                                    if(file_output) fprintf(file_output,"na\t");                                }
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            else {
                                printf("avg_EW\tna\tna\tna\t");
                                if(file_output) 
                                    fprintf(file_output,"avg_EW\tna\tna\tna\t");
                                print_percentages(data,sortvector,it2,file_output);
                            }
                            /*variance*/
                            if(data[0].n_loci > 2) {
                                avg = var = (double)0.0;
                                for(it=0,it2=0;it<(long int)data[0].n_iter;it++) {
                                    if((vars=variances(matrixmlsim[it].S2ewtest,matrixmlsim[it].Sewtest,data[0].n_loci)) != (double) -10000) {
                                        sortvector[it2] = vars;
                                        it2 += 1;
                                        avg += vars;
                                        var += vars * vars;
                                    }
                                }
                                if(it2 > 0) {
                                    var = varianceit(var,avg,it2);
                                    if(it2 > 1) qsort(sortvector,it2,sizeof(double),comp);
                                    printf("var_EW\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(file_output)
                                        fprintf(file_output,"var_EW\t%g\t%g\t",sortvector[(long int)floor(((double)it2*0.5))],avg/(double)it2);
                                    if(var != (double)-10000) {
                                        printf("%g\t",var);
                                        if(file_output) fprintf(file_output,"%g\t",var);
                                    }
                                    else {
                                        printf("na\t");
                                        if(file_output) fprintf(file_output,"na\t");
                                    }
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                                else {
                                    printf("var_EW\tna\tna\tna\n");
                                    if(file_output) 
                                        fprintf(file_output,"var_EW\tna\tna\tna\n");
                                    print_percentages(data,sortvector,it2,file_output);
                                }
                            }
                        }
                        free(sortvector);
					#if ALLRESULTS == 0
                    }
					#endif
                    break;
                default: /*error*/
                    return;
                    break;
            }
            break;
        default: /*error*/
            return;
            break;
    }
    return;
}

void printfignif(double prob,FILE *file_output)
{
    if(prob == -10000) {
		printf("na");
		if(file_output) fputs("na",file_output);
		return;
	}
	if(prob <= (double)0.001/2. || prob >= (double)1. - (double)(0.001/2.)) {
        printf("***");
        if(file_output) fputs("***",file_output);
    }
    else if(prob <= (double)0.01/2. || prob >= (double)1. - (double)(0.01/2.)) {
            printf("**");
            if(file_output) fputs("**",file_output);
        }
        else if(prob <= (double)0.05/2. || prob >= (double)1. - (double)(0.05/2.)) {
            printf("*");
            if(file_output) fputs("*",file_output);
        }
    
    return;
}
void printfignif1(double prob,FILE *file_output)
{
    if(prob == -10000) {
		printf("na");
		if(file_output) fputs("na",file_output);
		return;
	}
    if(prob <= (double)0.001) {
        printf("***");
        if(file_output) fputs("***",file_output);
    }
    else if(prob <= (double)0.01) {
            printf("**");
            if(file_output) fputs("**",file_output);
        }
        else if(prob <= (double)0.05) {
            printf("*");
            if(file_output) fputs("*",file_output);
        }
    
    return;
}
/*Sort data in qsort function*/
int comp(const void *i, const void *j)
{
    if(*(double *)i < *(double *)j) return -1;
    if(*(double *)i > *(double *)j) return  1;
    return 0;
}

void print_percentages(struct var *data,double *sortvector,long int it2, FILE *file_output)
{
    long int pcount,ncount;
    double valuep,valuen;
    
    if(data[0].n_iter >= 10) {
        if(it2 >= 10) {
            ncount = (long int)(it2*0.1);
            valuen = sortvector[ncount];
            while(sortvector[ncount] == sortvector[(long int)(it2*0.1)+1] && ncount > -1)  ncount--;
            pcount = (long int)(it2-it2*0.1);
            valuep = sortvector[pcount];
            while(sortvector[pcount] == sortvector[(long int)(it2*0.1)-1] && pcount < it2) pcount++;
            if(ncount == -1 && pcount < it2) {
                printf("na\t%g\t",sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"na\t%g\t",sortvector[(long int)pcount]);
            }
            if(ncount > -1 && pcount == it2) {
                printf("%g\tna\t",sortvector[(long int)ncount]);
                if(file_output) fprintf(file_output,"%g\tna\t",sortvector[(long int)ncount]);
            }
            if(ncount == -1 && pcount == it2) {
                printf("na\tna\t");
                if(file_output) fprintf(file_output,"na\tna\t");
            }
            if(ncount > -1 && pcount < it2) {
                printf("%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
            }
        }
        else {
            printf("na\tna\t");
            if(file_output) fprintf(file_output,"na\tna\t");
        }
    }

    if(data[0].n_iter >= 20) {
        if(it2 >= 20) {
            ncount = (long int)(it2*0.05);
            valuen = sortvector[ncount];
            while(sortvector[ncount] == sortvector[(long int)(it2*0.05)+1] && ncount > -1)  ncount--;
            pcount = (long int)(it2-it2*0.05);
            valuep = sortvector[pcount];
            while(sortvector[pcount] == sortvector[(long int)(it2*0.05)-1] && pcount < it2) pcount++;
            if(ncount == -1 && pcount < it2) {
                printf("na\t%g\t",sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"na\t%g\t",sortvector[(long int)pcount]);
            }
            if(ncount > -1 && pcount == it2) {
                printf("%g\tna\t",sortvector[(long int)ncount]);
                if(file_output) fprintf(file_output,"%g\tna\t",sortvector[(long int)ncount]);
            }
            if(ncount == -1 && pcount == it2) {
                printf("na\tna\t");
                if(file_output) fprintf(file_output,"na\tna\t");
            }
            if(ncount > -1 && pcount < it2) {
                printf("%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
            }
        }
        else {
            printf("na\tna\t");
            if(file_output) fprintf(file_output,"na\tna\t");
        }
	}
    if(data[0].n_iter >= 40) {
        if(it2 >= 40) {
            ncount = (long int)(it2*0.025);
            valuen = sortvector[ncount];
            while(sortvector[ncount] == sortvector[(long int)(it2*0.025)+1] && ncount > -1)  ncount--;
            pcount = (long int)(it2-it2*0.025);
            valuep = sortvector[pcount];
            while(sortvector[pcount] == sortvector[(long int)(it2*0.025)-1] && pcount < it2) pcount++;
            if(ncount == -1 && pcount < it2) {
                printf("na\t%g\t",sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"na\t%g\t",sortvector[(long int)pcount]);
            }
            if(ncount > -1 && pcount == it2) {
                printf("%g\tna\t",sortvector[(long int)ncount]);
                if(file_output) fprintf(file_output,"%g\tna\t",sortvector[(long int)ncount]);
            }
            if(ncount == -1 && pcount == it2) {
                printf("na\tna\t");
                if(file_output) fprintf(file_output,"na\tna\t");
            }
            if(ncount > -1 && pcount < it2) {
                printf("%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
            }
        }
        else {
            printf("na\tna\t");
            if(file_output) fprintf(file_output,"na\tna\t");
        }
	}
    if(data[0].n_iter >= 100) {
        if(it2 >= 100) {
            ncount = (long int)(it2*0.01);
            valuen = sortvector[ncount];
            while(sortvector[ncount] == sortvector[(long int)(it2*0.01)+1] && ncount > -1)  ncount--;
            pcount = (long int)(it2-it2*0.01);
            valuep = sortvector[pcount];
            while(sortvector[pcount] == sortvector[(long int)(it2*0.01)-1] && pcount < it2) pcount++;
            if(ncount == -1 && pcount < it2) {
                printf("na\t%g\t",sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"na\t%g\t",sortvector[(long int)pcount]);
            }
            if(ncount > -1 && pcount == it2) {
                printf("%g\tna\t",sortvector[(long int)ncount]);
                if(file_output) fprintf(file_output,"%g\tna\t",sortvector[(long int)ncount]);
            }
            if(ncount == -1 && pcount == it2) {
                printf("na\tna\t");
                if(file_output) fprintf(file_output,"na\tna\t");
            }
            if(ncount > -1 && pcount < it2) {
                printf("%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
            }
        }
        else {
            printf("na\tna\t");
            if(file_output) fprintf(file_output,"na\tna\t");
        }
	}
    if(data[0].n_iter >= 200) {
        if(it2 >= 200) {
            ncount = (long int)(it2*0.005);
            valuen = sortvector[ncount];
            while(sortvector[ncount] == sortvector[(long int)(it2*0.005)+1] && ncount > -1)  ncount--;
            pcount = (long int)(it2-it2*0.005);
            valuep = sortvector[pcount];
            while(sortvector[pcount] == sortvector[(long int)(it2*0.005)-1] && pcount < it2) pcount++;
            if(ncount == -1 && pcount < it2) {
                printf("na\t%g\t",sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"na\t%g\t",sortvector[(long int)pcount]);
            }
            if(ncount > -1 && pcount == it2) {
                printf("%g\tna\t",sortvector[(long int)ncount]);
                if(file_output) fprintf(file_output,"%g\tna\t",sortvector[(long int)ncount]);
            }
            if(ncount == -1 && pcount == it2) {
                printf("na\tna\t");
                if(file_output) fprintf(file_output,"na\tna\t");
            }
            if(ncount > -1 && pcount < it2) {
                printf("%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
            }
        }
        else {
            printf("na\tna\t");
            if(file_output) fprintf(file_output,"na\tna\t");
        }
	}
    if(data[0].n_iter >= 1000) {
        if(it2 >= 1000) {
            ncount = (long int)(it2*0.001);
            valuen = sortvector[ncount];
            while(sortvector[ncount] == sortvector[(long int)(it2*0.001)+1] && ncount > -1)  ncount--;
            pcount = (long int)(it2-it2*0.001);
            valuep = sortvector[pcount];
            while(sortvector[pcount] == sortvector[(long int)(it2*0.001)-1] && pcount < it2) pcount++;
            if(ncount == -1 && pcount < it2) {
                printf("na\t%g\t",sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"na\t%g\t",sortvector[(long int)pcount]);
            }
            if(ncount > -1 && pcount == it2) {
                printf("%g\tna\t",sortvector[(long int)ncount]);
                if(file_output) fprintf(file_output,"%g\tna\t",sortvector[(long int)ncount]);
            }
            if(ncount == -1 && pcount == it2) {
                printf("na\tna\t");
                if(file_output) fprintf(file_output,"na\tna\t");
            }
            if(ncount > -1 && pcount < it2) {
                printf("%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
            }
        }
        else {
            printf("na\tna\t");
            if(file_output) fprintf(file_output,"na\tna\t");
        }
	}
    if(data[0].n_iter >= 2000) {
        if(it2 >= 2000) {
            ncount = (long int)(it2*0.0005);
            valuen = sortvector[ncount];
            while(sortvector[ncount] == sortvector[(long int)(it2*0.0005)+1] && ncount > -1)  ncount--;
            pcount = (long int)(it2-it2*0.0005);
            valuep = sortvector[pcount];
            while(sortvector[pcount] == sortvector[(long int)(it2*0.0005)-1] && pcount < it2) pcount++;
            if(ncount == -1 && pcount < it2) {
                printf("na\t%g\t",sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"na\t%g\t",sortvector[(long int)pcount]);
            }
            if(ncount > -1 && pcount == it2) {
                printf("%g\tna\t",sortvector[(long int)ncount]);
                if(file_output) fprintf(file_output,"%g\tna\t",sortvector[(long int)ncount]);
            }
            if(ncount == -1 && pcount == it2) {
                printf("na\tna\t");
                if(file_output) fprintf(file_output,"na\tna\t");
            }
            if(ncount > -1 && pcount < it2) {
                printf("%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
                if(file_output) fprintf(file_output,"%g\t%g\t",sortvector[(long int)ncount],sortvector[(long int)pcount]);
            }
        }
        else {
            printf("na\tna\t");
            if(file_output) fprintf(file_output,"na\tna\t");
        }
	}
    printf("\n");
    if(file_output) fprintf(file_output,"\n");
    
    return;
}

void print_ibonf(double ibonf,int ibonfn,int n_loci,long int niter,FILE *file_output)
{
	double maxbon;
	if(ibonf > (double)0.5) {
		if((maxbon =(((double)1-ibonf)*(double)n_loci/(double)ibonfn)) <= (double)0.025) {
			if(niter < 40) {
				printf("*&\t");
				if(file_output) fprintf(file_output,"*&\t");
				return;
			}
			printf("*");
			if(file_output) fprintf(file_output,"*");
		}
		if(maxbon <= (double)0.005) {
			if(niter < 200) {
				printf("*&\t");
				if(file_output) fprintf(file_output,"*&\t");
				return;
			}
			printf("*");
			if(file_output) fprintf(file_output,"*");
		}
		if(maxbon <= (double)0.0005) {
			if(niter < 2000) {
				printf("*&\t");
				if(file_output) fprintf(file_output,"*&\t");
				return;
			}
			printf("*");
			if(file_output) fprintf(file_output,"*");
		}
		if(maxbon > (double)0.025) {
			printf("ns");
			if(file_output) fprintf(file_output,"ns");
		}
		printf("\t");
		if(file_output) fprintf(file_output,"\t");
	}
	else {
		if((maxbon =(ibonf*(double)n_loci/(double)ibonfn)) <= (double)0.025) {
			if(niter < 40) {
				printf("*&\t");
				if(file_output) fprintf(file_output,"*&\t");
				return;
			}
			printf("*");
			if(file_output) fprintf(file_output,"*");
		}
		if(maxbon <= (double)0.005) {
			if(niter < 200) {
				printf("*&\t");
				if(file_output) fprintf(file_output,"*&\t");
				return;
			}
			printf("*");
			if(file_output) fprintf(file_output,"*");
		}
		if(maxbon <= (double)0.0005) {
			if(niter < 2000) {
				printf("*&\t");
				if(file_output) fprintf(file_output,"*&\t");
				return;
			}
			printf("*");
			if(file_output) fprintf(file_output,"*");
		}
		if(maxbon > (double)0.025) {
			printf("ns");
			if(file_output) fprintf(file_output,"ns");
		}
		printf("\t");
		if(file_output) fprintf(file_output,"\t");
	}
	return;
}

void ibonf_sort(double **ibonf,int **ibonfn,int number,int loci)
{
	int count,x,y,z;
	double maxbon;
	int maxbon_n;
	double value;
	
	for(x=0;x<number;x++) 
		for(y=0;y<loci;y++) 
			ibonfn[x][y] = -1;
	
	for(x=0;x<number;x++) {
		count = 1;
		for(z=0;z<loci;z++) {
			maxbon = (double)-1;
			for(y=0;y<loci;y++) {
				if(ibonf[x][y] != (double)-10000) {
					value = ((double)0.5 - ibonf[x][y]);
					if(value < (double)0) value = -value;
					if(value >= maxbon) {
						if(ibonfn[x][y] == -1) {
							maxbon = value; 
							maxbon_n = y;
						}
					}
				}
			}
			if(maxbon != (double)-1) {
				ibonfn[x][maxbon_n] = count;
				count ++;
			}
		}
	}
	
	return;
}


