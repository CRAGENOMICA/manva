/*
 *  displaysimdata.c
 *  MuLoNeTests
 *
 *  Created by sonsins on Fri Mar 14 2003.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void open_simdispldetgrid(struct var *data,struct var2b *inputms,struct statistisim **matrixsim,struct statistisimmuloc *matrixmlsim,struct horizontalstatsml *avgstatloci,int *outgroup/*,int *n_loci*/,int onlymulo, FILE *file_output, int neuttest,int mulo,int dataobsequalsim,struct statistics *matrix, int observed_data)
{
    /*display a Table.*/
    char k[1];
    long int it;
    int x;
    double var,jc,is;
    double variances(double,double,int);
	FILE *file_output2;
	char fileout2[512];
    
    switch(mulo) {
        case 0: /*all loci*/
            switch(neuttest) {
                case 0: /* statistics, all loci, puf..*/
                    #if COMMAND_LINE
                    if(file_output) fflush(file_output);		
                    printf("\n\nMENU 3. Coalescent Monte Carlo simulations and Statistical inference:");
                    printf("\n     3.1. Display coalescent simulation analysis menu:");
                    printf("\n     3.1.0. Display statistics menu:");
                    printf("\n     3.1.0.1. Display detailed statistics menu:");
                    printf("\n     3.1.0.1.0. Display table with detailed statistics for each locus per iteration:\n\n");
                    
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
                        fprintf(file_output,"\n     3.1.0.1. Display detailed statistics menu:");
                        fprintf(file_output,"\n     3.1.0.1.0. Display table with detailed statistics for each locus per iteration:\n\n");
                        
                        fprintf(file_output," 0 - Display general statistics.\n");
                        fprintf(file_output," 1 - Display other linkage related statistics.\n");
                        fprintf(file_output," 2 - Display estimates of total locus variability.\n");
                        fprintf(file_output," 3 - Display estimates of nucleotide locus variability.\n");
                        fprintf(file_output," 4 - Back to previous menu.\n");
                    }
                    #endif
                    
                    do *k = getchar();
                    while(*k<'0' || *k>'4');

                    if(file_output) fprintf(file_output,"OPTION CHOSEN: %c\n\n",*k);
					
					if(*k<'4') {
						printf("\n\n This information will be sent into a new file. Please introduce the name of the file: ");
						scanf("%s",fileout2);

						if((file_output2 = fopen (fileout2,"w")) == 0) {
							printf("\n  It is not possible to open the file %s",fileout2);
							if(file_output) fprintf(file_output,"\n  It is not possible to open the file %s",fileout2);
							return;
						}
						if(file_output) fprintf(file_output," Data results are sent to file %s.\n",fileout2);
					}
                    
                    switch(*k) {
                        case '0': /* 0 - Display general statistics*/
                            if(*outgroup) {
                                /*Display statistics*//*
                                printf("\n\n Display statistics: \n");
                                printf(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n");
                                printf("samples, number of sites, segregating biallelic sites, \n");
                                if(data[0].time_spec > 0.) printf(" fixed mutations, averaged divergence,");
                                printf(" the number of haplotypes (not shared, if observed), haplotypes/sample size, haplotype diversity,\n");
                                printf(" the recombination parameter per loci (Hudson 1987) and Rm.\n");
                                for(x=0;x<data[0].n_loci;x++) {
                                    printf("\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
									printf("\n#samples(%d): %d\t#sites(%d): %ld\n",x,inputms[x].nsam,x,inputms[x].nsites);
                                    printf("\n#iter\t");
                                    printf("#S(%d)\t",x);
                                    if(data[0].time_spec > 0.) printf("#fixed(%d)\tDiverg(%d)\t",x,x);
                                    printf("#haplot(%d)\thaplsam(%d)\thapldiv(%d)\tC(%d)\tRm(%d)\n",x,x,x,x,x);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        printf("%ld\t",it);
                                        printf("%d\t",matrixsim[x][it].biallsites);
                                        if(data[0].time_spec > 0.) printf("%d\t",matrixsim[x][it].fixed);
                                        if(data[0].time_spec > 0.) printf("%g\t",matrixsim[x][it].ndivergence);
                                        printf("%d\t",matrixsim[x][it].nhapl);
										printf("%g\t",matrixsim[x][it].nhaplsam);
										printf("%g\t",matrixsim[x][it].hapldiv);
                                        printf("%g\n",matrixsim[x][it].Rvpi);
                                        printf("%d\n",matrixsim[x][it].Rm);
                                    }
                                    printf("\n\n");                            
                                }*/
                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display statistics: \n",file_output2);
                                    fputs(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n",file_output2);
                                    fputs("samples, number of sites, segregating biallelic sites, \n",file_output2);
                                    if(data[0].time_spec > 0.) fputs(" fixed mutations, averaged divergence,",file_output2);
                                    fputs(" the number of haplotypes (not shared, if observed), haplotypes/sample size, haplotype diversity,\n",file_output2);
                                    fputs(" and Rm.\n",file_output2);
                                    for(x=0;x<data[0].n_loci;x++) {
										fprintf(file_output2,"\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                        fprintf(file_output2,"\n#samples(%d): %d\t#sites(%d): %ld\n",x,inputms[x].nsam,x,inputms[x].nsites);
                                        fprintf(file_output2,"\n#iter\t");
                                        fprintf(file_output2,"#sites_nomh(%d)\t",x);
                                        fprintf(file_output2,"#S(%d)\t",x);
                                        if(data[0].time_spec > 0.) fprintf(file_output2,"#fixed(%d)\tDiverg(%d)\t",x,x);
                                        fprintf(file_output2,"#haplot(%d)\thaplsam(%d)\thapldiv(%d)\tRm(%d)\n",x,x,x,x);
                                        for(it=0;it<(long int)data[0].n_iter;it++) {
                                            fprintf(file_output2,"%ld\t",it);
                                            fprintf(file_output2,"%ld\t",matrixsim[x][it].nsites);
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].biallsites);
											if(data[0].time_spec > 0.) {
                                                fprintf(file_output2,"%d\t",matrixsim[x][it].fixed);
                                                fprintf(file_output2,"%g\t",matrixsim[x][it].ndivergence);
											}
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].nhapl);
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].nhaplsam);
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].hapldiv);
                                            fprintf(file_output2,"%d\n",matrixsim[x][it].Rm);
                                        }
                                        fputs("\n\n",file_output2);                            
                                    }
                                }
                            }
                            else {
                                /*not outgroup*/
                                /*Display statistics*//*
                                printf("\n\n Display statistics: \n");
                                printf(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n");
                                printf("samples, number of sites, segregating biallelic sites,\n");
								printf(" the number of haplotypes (not shared, if observed), haplotypes/sample size, haplotype diversity,\n");
                                printf(" the recombination parameter per loci (Hudson 1987) and Rm.\n");
                                for(x=0;x<data[0].n_loci;x++) {
                                    printf("\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
									printf("\n#samples(%d): %d\t#sites(%d): %ld\n",x,inputms[x].nsam,x,inputms[x].nsites);
                                    printf("\n#iter\t");
                                    printf("#S(%d)\t",x);
                                    printf("#haplot(%d)\thaplsam(%d)\thapldiv(%d)\tC(%d)\tRm(%d)\n",x,x,x,x,x);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        printf("%ld\t",it);
                                        printf("%d\t",matrixsim[x][it].biallsites);
                                        printf("%d\t",matrixsim[x][it].nhapl);
										printf("%g\t",matrixsim[x][it].nhaplsam);
										printf("%g\t",matrixsim[x][it].hapldiv);
                                        printf("%g\n",matrixsim[x][it].Rvpi);
                                        printf("%d\n",matrixsim[x][it].Rm);
                                    }
                                    printf("\n\n");                            
                                }*/
                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display statistics: \n",file_output2);
                                    fputs(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n",file_output2);
                                    fputs("samples, number of sites, segregating biallelic sites,\n",file_output2);
									fputs(" the number of haplotypes (not shared, if observed), haplotypes/sample size, haplotype diversity,\n",file_output2);
                                    fputs(" and Rm.\n",file_output2);
                                    for(x=0;x<data[0].n_loci;x++) {
										fprintf(file_output2,"\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                        fprintf(file_output2,
                                            "\n#samples(%d): %d\t#sites(%d): %ld\n",
                                            x,inputms[x].nsam,x,inputms[x].nsites);
                                        fprintf(file_output2,"\n#iter\t");
                                        fprintf(file_output2,"#sites_nomh(%d)\t",x);
                                        fprintf(file_output2,"#S(%d)\t",x);
                                        fprintf(file_output2,"#haplot(%d)\thaplsam(%d)\thapldiv(%d)\tRm(%d)\n",x,x,x,x);
                                        for(it=0;it<(long int)data[0].n_iter;it++) {
                                            fprintf(file_output2,"%ld\t",it);
                                            fprintf(file_output2,"%ld\t",matrixsim[x][it].nsites);
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].biallsites);
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].nhapl);
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].nhaplsam);
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].hapldiv);
                                            fprintf(file_output2,"%d\n",matrixsim[x][it].Rm);
                                        }
                                        fputs("\n\n",file_output2);                            
                                    }
                                }
                            }
							fclose(file_output2);
                            break;
                        case '1': /* 1 - Display linkage related statistics*/
                            if(*outgroup) {
                                /*Display statistics*//*
                                printf("\n\n Display statistics: \n");
                                printf(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n");
                                printf("samples, number of sites, segregating biallelic sites,\n");
                                printf("Wall's not divided statistics B and Q,\n");
                                printf(" and Rozas' et al. not divided ZA (shared excluded in all three statistics, if observed)\n");
                                for(x=0;x<data[0].n_loci;x++) {
									printf("\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                    printf("\n#samples(%d): %d\t#sites(%d): %ld\n",x,inputms[x].nsam,x,inputms[x].nsites);
                                    printf("\n#iter\t");
                                    printf("#S(%d)\t",x);
                                    printf("b(%d)\tq(%d)\tza(%d)\n",x,x,x);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        printf("%ld\t",it);
                                        printf("%d\t",matrixsim[x][it].biallsites);
                                        printf("%d\t",matrixsim[x][it].b);
                                        printf("%d\t",matrixsim[x][it].q);
                                        printf("%g\n",matrixsim[x][it].za);
                                    }
                                    printf("\n\n");                            
                                }*/
                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display statistics: \n",file_output2);
                                    fputs(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n",file_output2);
                                    fputs("samples, number of sites, segregating biallelic sites,\n",file_output2);
                                    fputs("Wall's not divided statistics B and Q,\n",file_output2);
                                    fputs(" and Rozas' et al. not divided ZA (shared excluded in all three statistics, if observed)\n",file_output2);
                                    for(x=0;x<data[0].n_loci;x++) {
										fprintf(file_output2,"\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                        fprintf(file_output2,
                                            "\n#samples(%d): %d\t#sites(%d): %ld\n",
                                            x,inputms[x].nsam,x,inputms[x].nsites);
                                        fprintf(file_output2,"\n#iter\t");
                                        fprintf(file_output2,"#sites_nomh(%d)\t",x);
                                        fprintf(file_output2,"#S(%d)\t",x);
                                        fprintf(file_output2,"b(%d)\tq(%d)\tza(%d)\n",x,x,x);
                                        for(it=0;it<(long int)data[0].n_iter;it++) {
                                            fprintf(file_output2,"%ld\t",it);
                                            fprintf(file_output2,"%ld\t",matrixsim[x][it].nsites);
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].biallsites);
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].b);
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].q);
                                            fprintf(file_output2,"%g\n",matrixsim[x][it].za);
                                        }
                                        fputs("\n\n",file_output2);                            
                                    }
                                }
                            }
                            else {
                                /*not outgroup*/
                                /*Display statistics*//*
                                printf("\n\n Display statistics: \n");
                                printf(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n");
                                printf("samples, number of sites, segregating biallelic sites,\n");
                                printf("Wall's not divided statistics B and Q,\n");
                                printf(" and Rozas' et al. not divided ZA\n");
                                for(x=0;x<data[0].n_loci;x++) {
									printf("\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                    printf("\n#samples(%d): %d\t#sites(%d): %ld\n",x,inputms[x].nsam,x,inputms[x].nsites);
                                    printf("\n#iter\t");
                                    printf("#S(%d)\t",x);
                                    printf("b(%d)\tq(%d)\tza(%d)\n",x,x,x);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        printf("%ld\t",it);
                                        printf("%d\t",matrixsim[x][it].biallsites);
                                        printf("%d\t",matrixsim[x][it].b);
                                        printf("%d\t",matrixsim[x][it].q);
                                        printf("%g\n",matrixsim[x][it].za);
                                    }
                                    printf("\n\n");                            
                                }*/
                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display statistics: \n",file_output2);
                                    fputs(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n",file_output2);
                                    fputs("samples, number of sites, segregating biallelic sites,\n",file_output2);
                                    fputs("Wall's not divided statistics B and Q,\n",file_output2);
                                    fputs(" and Rozas' et al. not divided ZA\n",file_output2);
                                    for(x=0;x<data[0].n_loci;x++) {
										fprintf(file_output2,"\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                        fprintf(file_output2,
                                            "\n#samples(%d): %d\t#sites(%d): %ld\n",
                                            x,inputms[x].nsam,x,inputms[x].nsites);
                                        fprintf(file_output2,"\n#iter\t");
                                        fprintf(file_output2,"#sites_nomh(%d)\t",x);
                                        fprintf(file_output2,"#S(%d)\t",x);
                                        fprintf(file_output2,"b(%d)\tq(%d)\tza(%d)\n",x,x,x);
                                        for(it=0;it<(long int)data[0].n_iter;it++) {
                                            fprintf(file_output2,"%ld\t",it);
                                            fprintf(file_output2,"%ld\t",matrixsim[x][it].nsites);
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].biallsites);
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].b);
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].q);
                                            fprintf(file_output2,"%g\n",matrixsim[x][it].za);
                                        }
                                        fputs("\n\n",file_output2);                            
                                    }
                                }
                            }
							fclose(file_output2);
                            break;
                        case '2': /* 2 - Display estimates of total locus variability*/
                            if(*outgroup) {
                                /*Display statistics*//*
                                printf("\n\n Display statistics: \n");
                                printf(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n");
                                printf("samples, number of sites, segregating biallelic sites, \n");
                                if(data[0].time_spec > 0.) 
                                    printf(" averaged divergence (uncorrected and corrected by Jukes and Cantor) \n");
                                printf(" and locus estimations of variability from Watterson, Tajima, Fu and Li,\n");
                                printf(" Fay and Wu and Zeng et al.\n")
                                printf(" NOTE that all theta estimates are corrected for chromosome population size.\n");
								for(x=0;x<data[0].n_loci;x++) {
									printf("\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                    printf("\n#samples(%d): %d\t#sites(%d): %ld\n",x,inputms[x].nsam,x,inputms[x].nsites);
                                    printf("\n#iter\t");
                                    printf("#S(%d)\t",x);
                                    if(data[0].time_spec > 0.) 
                                        printf("Diverg(%d)\tDivJC(%d)\t",x,x);
                                    printf("theta_wat(%d)\ttheta_taj(%d)\t",x,x);
                                    printf("theta_fuli(%d)\ttheta_faywu(%d)\ttheta_zeng(%d)\n",x,x,x);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        printf("%ld\t",it);
                                        printf("%d\t",matrixsim[x][it].biallsites);
                                        if(data[0].time_spec > 0.) {
                                            printf("%g\t",matrixsim[x][it].ndivergence);
                                            if((jc = matrixsim[x][it].ndivergence/(double)inputms[x].nsites) < (double)0.75) {
                                                jc = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * jc) * (double)inputms[x].nsites;
                                                printf("%g\t",jc);
                                            }
                                            else printf("na\t");
                                        }
                                        printf("%g\t",matrixsim[x][it].theta_wat*((double)1/inputms[x].factor_chrn));
                                        printf("%g\t",matrixsim[x][it].theta_taj*((double)1/inputms[x].factor_chrn));
                                        printf("%g\t",matrixsim[x][it].theta_fuli*((double)1/inputms[x].factor_chrn));
                                        printf("%g\n",matrixsim[x][it].theta_fw*((double)1/inputms[x].factor_chrn));
                                        printf("%g\n",matrixsim[x][it].theta_L*((double)1/inputms[x].factor_chrn));
                                    }
                                    printf("\n\n");                            
                                }*/
                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display statistics: \n",file_output2);
                                    fputs(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n",file_output2);
                                    fputs("samples, number of sites, segregating biallelic sites, \n",file_output2);
                                    if(data[0].time_spec > 0.) 
                                        fputs(" averaged divergence (uncorrected and corrected by Jukes and Cantor) \n",file_output2);
                                    fputs(" and locus estimations of variability from Watterson, Tajima, Fu and Li,\n",file_output2);
                                    fputs(" Fay and Wu and Zeng et al. Also from HKA test.\n",file_output2);
									fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n",file_output2);
                                    for(x=0;x<data[0].n_loci;x++) {
    									fprintf(file_output2,"\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
										fprintf(file_output2,
                                            "\n#samples(%d): %d\t#sites(%d): %ld\n",
                                            x,inputms[x].nsam,x,inputms[x].nsites);
                                        fprintf(file_output2,"\n#iter\t");
                                        fprintf(file_output2,"#sites_nomh(%d)\t",x);
                                        fprintf(file_output2,"#S(%d)\t",x);
                                        if(data[0].time_spec > 0.) fprintf(file_output2,"Diverg(%d)\tDivJC(%d)\t",x,x);
                                        fprintf(file_output2,"theta_wat(%d)\ttheta_taj(%d)\t",x,x);
                                        fprintf(file_output2,"theta_fuli(%d)\ttheta_faywu(%d)\ttheta_zeng(%d)\ttheta_HKA(%d)\n",x,x,x,x);
                                        for(it=0;it<(long int)data[0].n_iter;it++) {
                                            fprintf(file_output2,"%ld\t",it);
                                            fprintf(file_output2,"%ld\t",matrixsim[x][it].nsites);
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].biallsites);
                                            if(data[0].time_spec > 0.) {
                                                fprintf(file_output2,"%g\t",matrixsim[x][it].ndivergence);
                                                if((jc = matrixsim[x][it].ndivergence/(double)inputms[x].nsites) < (double)0.75) {
                                                    jc = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * jc) * (double)inputms[x].nsites;
                                                    fprintf(file_output2,"%g\t",jc);
                                                }
                                                else fputs("na\t",file_output2);
                                            }
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_wat*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_taj*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_fuli*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_fw*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_L*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\n",matrixsim[x][it].hka_theta*((double)1/inputms[x].factor_chrn));
                                        }
                                        fputs("\n\n",file_output2);                            
                                    }
                                }
                            }
                            else {
                                /*not outgroup*/
                                /*Display statistics*//*
                                printf("\n\n Display statistics: \n");
                                printf(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n");
                                printf("samples, number of sites, segregating biallelic sites,\n");
                                printf(" locus estimations of variability from Watterson, Tajima and Fu and Li.\n");
                                printf(" NOTE that all theta estimates are corrected for chromosome population size.\n");
								for(x=0;x<data[0].n_loci;x++) {
									printf("\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                    printf("\n#samples(%d): %d\t#sites(%d): %ld\n",x,inputms[x].nsam,x,inputms[x].nsites);
                                    printf("\n#iter\t");
                                    printf("#S(%d)\t",x);
                                    printf("theta_wat(%d)\ttheta_taj(%d)\t",x,x);
                                    printf("theta_fuli(%d)\n",x);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        printf("%ld\t",it);
                                        printf("%d\t",matrixsim[x][it].biallsites);
                                        printf("%g\t",matrixsim[x][it].theta_wat*((double)1/inputms[x].factor_chrn));
                                        printf("%g\t",matrixsim[x][it].theta_taj*((double)1/inputms[x].factor_chrn));
                                        printf("%g\n",matrixsim[x][it].theta_fulin*((double)1/inputms[x].factor_chrn));
                                    }
                                    printf("\n\n");                            
                                }*/
                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display statistics: \n",file_output2);
                                    fputs(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n",file_output2);
                                    fputs("samples, number of sites, segregating biallelic sites,\n",file_output2);
                                    fputs("\n#iter\t",file_output2);
                                    fputs(" locus estimations of variability from Watterson, Tajima and Fu and Li.\n",file_output2);
                                    fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n",file_output2);
									for(x=0;x<data[0].n_loci;x++) {
										fprintf(file_output2,"\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                        fprintf(file_output2,
                                            "\n#samples(%d): %d\t#sites(%d): %ld\n",
                                            x,inputms[x].nsam,x,inputms[x].nsites);
                                        fprintf(file_output2,"\n#iter\t");
                                        /*fprintf(file_output2,"#sites_nomh(%d)\t",x);*/
                                        fprintf(file_output2,"#S(%d)\t",x);
                                        fprintf(file_output2,"theta_wat(%d)\ttheta_taj(%d)\t",x,x);
                                        fprintf(file_output2,"theta_fulin(%d)\n",x);
                                        for(it=0;it<(long int)data[0].n_iter;it++) {
                                            fprintf(file_output2,"%ld\t",it);
                                            /*fprintf(file_output2,"%ld\t",matrixsim[x][it].nsites);*/
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].biallsites);
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_wat*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_taj*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\n",matrixsim[x][it].theta_fulin*((double)1/inputms[x].factor_chrn));
                                        }
                                        fputs("\n\n",file_output2);                            
                                    }
                                }
                            }
							fclose(file_output2);
                            break;
                        case '3': /* 3 - Display estimates of nucleotide locus variability*/
                            if(*outgroup) {
                                /*Display statistics*//*
                                printf("\n\n Display statistics: \n");
                                printf(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n");
                                printf("samples, number of sites, segregating biallelic sites,\n");
                                if(data[0].time_spec > 0.) 
                                    printf(" nucleotide divergence (uncorrected and corrected by Jukes and Cantor) \n");
                                printf(" and nucleotide estimations of variability from Watterson, Tajima, Fu and Li,\n");
                                printf(" Fay and Wu and Zeng et al.\n");
                                printf(" NOTE that all theta estimates are corrected for chromosome population size.\n");
								for(x=0;x<data[0].n_loci;x++) {
									printf("\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                    printf("\n#samples(%d): %d\t#sites(%d): %ld\n",x,inputms[x].nsam,x,inputms[x].nsites);
                                    printf("\n#iter\t");
                                    printf("#S(%d)\t",x);
                                    if(data[0].time_spec > 0.) printf("Diverg(%d)\tDivJC(%d)\t",x,x);
                                    printf("theta_wat(%d)\ttheta_taj(%d)\t",x,x);
                                    printf("theta_fuli(%d)\ttheta_faywu(%d)\ttheta_zeng(%d)\n",x,x,x);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        printf("%ld\t",it);
                                        printf("%d\t",matrixsim[x][it].biallsites);
                                        if(data[0].time_spec > 0.) {
                                            printf("%g\t",matrixsim[x][it].ndivergence/(double)inputms[x].nsites);
                                            if((jc = matrixsim[x][it].ndivergence/(double)inputms[x].nsites) < (double)0.75) {
                                                jc = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * jc);
                                                printf("%g\t",jc);
                                            }
                                            else printf("na\t");
                                        }
                                        printf("%g\t",matrixsim[x][it].theta_wat/
                                            (inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                        printf("%g\t",matrixsim[x][it].theta_taj/
                                            ((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                        printf("%g\t",matrixsim[x][it].theta_fuli/
                                            ((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                        printf("%g\n",matrixsim[x][it].theta_fw/
                                            ((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                        printf("%g\n",matrixsim[x][it].theta_L/
                                            ((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                    }
                                    printf("\n\n");                            
                                }*/
                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display statistics: \n",file_output2);
                                    fputs(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n",file_output2);
                                    fputs("samples, number of sites, segregating biallelic sites,\n",file_output2);
                                    if(data[0].time_spec > 0.) 
                                        fputs(" nucleotide divergence (uncorrected and corrected by Jukes and Cantor) \n",file_output2);
                                    fputs(" and nucleotide estimations of variability from Watterson, Tajima, Fu and Li.\n",file_output2);
                                    fputs(" Fay and Wu and Zeng et al. Also theta from HKA test.\n",file_output2);
                                    fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n",file_output2);
									for(x=0;x<data[0].n_loci;x++) {
										fprintf(file_output2,"\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                        fprintf(file_output2,
                                            "\n#samples(%d): %d\t#sites(%d): %ld\n",
                                            x,inputms[x].nsam,x,inputms[x].nsites);
                                        fprintf(file_output2,"\n#iter\t");
                                        fprintf(file_output2,"#S(%d)\t",x);
                                        if(data[0].time_spec > 0.) fprintf(file_output2,"Diverg(%d)\tDivJC(%d)\t",x,x);
                                        fprintf(file_output2,"theta_wat(%d)\ttheta_taj(%d)\t",x,x);
                                        fprintf(file_output2,"theta_fuli(%d)\ttheta_faywu(%d)\ttheta_zeng(%d)\ttheta_HKA(%d)\n",x,x,x,x);
                                        for(it=0;it<(long int)data[0].n_iter;it++) {
                                            fprintf(file_output2,"%ld\t",it);
											fprintf(file_output2,"#sites_nomh(%d)\t",x);
                                            fprintf(file_output2,"%ld\t",matrixsim[x][it].nsites);
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].biallsites);
                                            if(data[0].time_spec > 0.) {
                                                fprintf(file_output2,"%g\t",matrixsim[x][it].ndivergence/(double)inputms[x].nsites);
                                                if((jc = matrixsim[x][it].ndivergence/(double)inputms[x].nsites) < (double)0.75) {
                                                    jc = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * jc);
                                                    fprintf(file_output2,"%g\t",jc);
                                                }
                                                else fputs("na\t",file_output2);
                                            }
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_wat/(inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_taj/((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_fuli/((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_fw/((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_L/((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
											fprintf(file_output2,"%g\n",matrixsim[x][it].hka_theta/((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));                                    }
                                        fputs("\n\n",file_output2);                            
                                    }
                                }
                            }
                            else {
                                /*not outgroup*/
                                /*Display statistics*//*
                                printf("\n\n Display statistics: \n");
                                printf(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n");
                                printf("samples, number of sites, segregating biallelic sites,\n");
                                printf("nucleotide estimations of variability from Watterson, Tajima and Fu and Li.\n");
                                printf(" NOTE that all theta estimates are corrected for chromosome population size.\n");
								for(x=0;x<data[0].n_loci;x++) {
									printf("\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                    printf("\n#samples(%d): %d\t#sites(%d): %ld\n",x,inputms[x].nsam,x,inputms[x].nsites);
                                    printf("\n#iter\t");
                                    printf("#S(%d)\t",x);
                                    printf("theta_wat(%d)\ttheta_taj(%d)\t",x,x);
                                    printf("theta_fuli(%d)\n",x);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        printf("%ld\t",it);
                                        printf("%d\t",inputms[x].nsam);
                                        printf("%g\t",matrixsim[x][it].theta_wat/((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                        printf("%g\t",matrixsim[x][it].theta_taj/((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                        printf("%g\n",matrixsim[x][it].theta_fulin/((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                    }
                                    printf("\n\n");                            
                                }*/
                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display statistics: \n",file_output2);
                                    fputs(" The table contains in order the number of iteration, and for each locus the chromosome population size, the number of\n",file_output2);
                                    fputs("samples, number of sites, segregating biallelic sites,\n",file_output2);
                                    fputs("nucleotide estimations of variability from Watterson, Tajima and Fu and Li.\n",file_output2);
                                    fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n",file_output2);
									for(x=0;x<data[0].n_loci;x++) {
										fprintf(file_output2,"\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                        fprintf(file_output2,
                                            "\n#samples(%d): %d\t#sites(%d): %ld\n",
                                            x,inputms[x].nsam,x,inputms[x].nsites);
                                        fprintf(file_output2,"\n#iter\t");
                                        /*fprintf(file_output2,"#sites_nomh(%d)\t",x);*/
                                        fprintf(file_output2,"#S(%d)\t",x);
                                        fprintf(file_output2,"theta_wat(%d)\ttheta_taj(%d)\t",x,x);
                                        fprintf(file_output2,"theta_fulin(%d)\n",x);
                                        for(it=0;it<(long int)data[0].n_iter;it++) {
                                            fprintf(file_output2,"%ld\t",it);
                                            /*fprintf(file_output2,"%ld\t",matrixsim[x][it].nsites);*/
                                            fprintf(file_output2,"%d\t",matrixsim[x][it].biallsites);
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_wat/((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\t",matrixsim[x][it].theta_taj/((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                            fprintf(file_output2,"%g\n",matrixsim[x][it].theta_fulin/((double)inputms[x].nsites)*((double)1/inputms[x].factor_chrn));
                                        }
                                        fputs("\n\n",file_output2);                            
                                    }
                                }
                            }
							fclose(file_output2);
                            break;
                        case '4': /* 4 - Back to previous menu*/
                            return;
                            break;
                    }
                    break;
                case 1: /* neutrality tests, all loci, puf..*/
					if(file_output) {    
						fprintf(file_output,"\n\n     MENU:\n     3. Statistical inference based on Coalescent Monte Carlo simulations menu:");
						fprintf(file_output,"\n     3.1. Display coalescent simulation analysis menu:");
						fprintf(file_output,"\n     3.1.1. Display neutrality tests menu:");
						fprintf(file_output,"\n     3.1.1.1. Display detailed neutrality tests menu:\n\n");
						fprintf(file_output,"\n     3.1.1.1.0. Display the results in a table (VERY LONG OUTPUT!).\n");
					}
					
					printf("\n\n This information will be sent into a new file. Please introduce the name of the file: ");
					scanf("%s",fileout2);

					if((file_output2 = fopen (fileout2,"w")) == 0) {
						printf("\n  It is not possible to open the file %s",fileout2);
						if(file_output) fprintf(file_output,"\n  It is not possible to open the file %s",fileout2);
						return;
					}
					fprintf(file_output," Data results are sent to file %s.\n",fileout2);

                    if(*outgroup) {/*
                        printf("\n\n Display detailed neutral tests: \n\n Multiple hits not included in the analysis.\n");
                        printf(" The table contains the factor correction for chromosome population size, the number of samples for each locus,\n");
						printf(" Tajima's D test, Fu and Li's D and F tests, Fay and Wu's H test, normalized Fay and Wu H test, Fu's Fs test,\n");
						printf(" Rozas' et a. ZA test, Wall's B and Q tests, Ramos-Onsins & Rozas (R2) test, Zeng E test and Ewens-Watterson test. \n");
                        printf(" Shared polymorphisms (if observed) are excluded from aanlysis in the above tests. \n");
                        if(data[0].n_loci > 1 && data[0].time_spec > 0.0) 
                            printf(" Also the result of the partial Chi-square value for HKA test in each locus (values corrected by Jukes and Cantor and correcting for factor_chrn).\n");
                        printf(" 'na' means the test is not calculated because the number of samples or the number of segregating sites\n");
                        printf("were insufficient.\n");
                        
                        for(x=0;x<data[0].n_loci;x++) {
							printf("\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                            printf("\n#samples(%d): %d\t#sites(%d): %ld\n",x,inputms[x].nsam,x,inputms[x].nsites);
                            printf("\n#iter\t");
                            printf("TajD(%d)\tFuLiD(%d)\tFuLiF(%d)\tFayWuH(%d)\tFayWuHn(%d)\tFuFs(%d)\t",x,x,x,x,x,x);
                            printf("RozasZA(%d)\tWallB(%d)\tWallQ(%d)\tR2(%d)\tE(%d)\tEW(%d)\t",x,x,x,x,x,x); 
                            if(data[0].n_loci > 1 && data[0].time_spec > 0.0)  printf("partialhka X\t");
                            printf("\n");  
                            for(it=0;it<(long int)data[0].n_iter;it++) {
                                printf("%ld\t",it);
                                if(matrixsim[x][it].tajimaD < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].tajimaD);
                                if(matrixsim[x][it].fuliD < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].fuliD);
                                if(matrixsim[x][it].fuliF < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].fuliF);
                                if(matrixsim[x][it].faywuHo == (double) -10000) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].faywuHo);
                                if(matrixsim[x][it].faywuH == (double) -10000) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].faywuH);
                                if(matrixsim[x][it].fuFs == (double) -10000) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].fuFs);
                                if(matrixsim[x][it].rZA < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].rZA);
                                if(matrixsim[x][it].wB < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].wB);
                                if(matrixsim[x][it].wQ < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].wQ);
                                if(matrixsim[x][it].R2 < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].R2);
                                if(matrixsim[x][it].zengE < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].zengE);
                                if(matrixsim[x][it].ewtest < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].ewtest);
                                if(data[0].n_loci > 1 && data[0].time_spec > 0.0) {
									printf("%g\t",matrixsim[x][it].hka);
									fprintf(file_output2,"%g\t",matrixsim[x][it].hka_theta);
								}
                                printf("\n");
                            }
                            printf("\n\n");
                        }*/
                        if(file_output2) {
                            /*print to file*/
                            fputs("\n\n Display detailed neutral tests: \n\n Multiple hits not included in the analysis.\n",file_output2);
                            fputs(" The table contains the factor correction for chromosome population size, the number of samples for each locus,\n",file_output2);
                            fputs(" Tajima's D test, Fu and Li's D and F tests, Fay and Wu's H test, normalized Fay and Wu H test, Fu's Fs test,\n",file_output2);
                            fputs(" Rozas' et a. ZA test, Wall's B and Q tests, Ramos-Onsins & Rozas (R2) test, Zeng E test and Ewens-Watterson test. \n",file_output2);
                            fputs(" Shared polymorphisms (if observed) are excluded from aanlysis in the above tests. \n",file_output2);
                            if(data[0].n_loci > 1 && data[0].time_spec > 0.0) fputs(" Also the result of the partial Chi-square value for HKA test in each locus and the theta value \n(values corrected by Jukes and Cantor and correcting for factor_chrn).\n",file_output2);
                            fputs(" 'na' means the test is not calculated because the number of samples or the number of segregating sites\n",file_output2);
                            fputs("were insufficient.\n",file_output2);
                            
                            for(x=0;x<data[0].n_loci;x++) {
								fprintf(file_output2,"\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                fprintf(file_output2,
                                    "\n#samples(%d): %d\t#sites(%d): %ld\n",
                                    x,inputms[x].nsam,x,inputms[x].nsites);
                                fprintf(file_output2,"\n#iter\t");
								fprintf(file_output2,"#sites_nomh(%d)\t",x);
                                fprintf(file_output2,"TajD(%d)\tFuLiD(%d)\tFuLiF(%d)\tFayWuH(%d)\tFayWuHn(%d)\tFuFs(%d)\t",x,x,x,x,x,x);
                                fprintf(file_output2,"RozasZA(%d)\tWallB(%d)\tWallQ(%d)\tR2(%d)\tE(%d)\tEW(%d)\t",x,x,x,x,x,x); 
                                if(data[0].n_loci > 1 && data[0].time_spec > 0.0)  fprintf(file_output2,"partialhka_X(%d)\t",x);
                                fputs("\n",file_output2);  
                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                    fprintf(file_output2,"%ld\t",it);
									fprintf(file_output2,"%ld\t",matrixsim[x][it].nsites);
                                    if(matrixsim[x][it].tajimaD < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].tajimaD);
                                    if(matrixsim[x][it].fuliD < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].fuliD);
                                    if(matrixsim[x][it].fuliF < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].fuliF);
                                    if(matrixsim[x][it].faywuH == (double) -10000) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].faywuHo);
                                    if(matrixsim[x][it].faywuHo == (double) -10000) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].faywuH);
                                    if(matrixsim[x][it].fuFs == (double) -10000) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].fuFs);
                                    if(matrixsim[x][it].rZA < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].rZA);
                                    if(matrixsim[x][it].wB < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].wB);
                                    if(matrixsim[x][it].wQ < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].wQ);
                                    if(matrixsim[x][it].R2 < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].R2);
                                    if(matrixsim[x][it].zengE < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].zengE);
                                    if(matrixsim[x][it].ewtest < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].ewtest);
                                    if(data[0].n_loci > 1 && data[0].time_spec > 0.0) {
                                        fprintf(file_output2,"%g\t",matrixsim[x][it].hka);
										/*fprintf(file_output2,"%g\t",matrixsim[x][it].hka_theta);*/
									}
                                    fputs("\n",file_output2);
                                }
                                fputs("\n\n",file_output2);
                            }
                        }
                    }
                    else {/*
                        printf("\n\n Display detailed neutral tests: \n\n Multiple hits not included in the analysis.\n");
                        printf(" The table contains the factor correction for chromosome population size, the number of samples for each locus,\n");
                        printf(" Tajima's D test, Fu and Li's D and F tests, normalized Fay and Wu H test, Fu's Fs test,\n");
                        printf(" Rozas' et a. ZA test, Wall's B, Q tests, Ramos-Onsins & Rozas (R2) test and Ewens-Watterson test.\n");
                        printf(" 'na' means the test is not calculated because the number of samples or the number of segregating sites\n");
                        printf("were insufficient.\n");
                        
                        for(x=0;x<data[0].n_loci;x++) {
							printf("\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                            printf("\n#samples(%d): %d\t#sites(%d): %ld\n",x,inputms[x].nsam,x,inputms[x].nsites);
                            printf("\n#iter\t");
                            printf("TajD(%d)\tFuLiD*(%d)\tFuLiF*(%d)\tFuFs(%d)\t",x,x,x,x);
                            printf("RozasZA(%d)\tWallB(%d)\tWallQ(%d)\tR2(%d)\tEW(%d)\n",x,x,x,x,x);    
                            for(it=0;it<(long int)data[0].n_iter;it++) {
                                printf("%ld\t",it);
                                if(matrixsim[x][it].tajimaD < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].tajimaD);
                                if(matrixsim[x][it].fuliDn < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].fuliDn);
                                if(matrixsim[x][it].fuliFn < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].fuliFn);
                                if(matrixsim[x][it].fuFs == (double) -10000) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].fuFs);
                                if(matrixsim[x][it].rZA < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].rZA);
                                if(matrixsim[x][it].wB < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].wB);
                                if(matrixsim[x][it].wQ < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].wQ);
                                if(matrixsim[x][it].R2 < -9999) printf("na\t");
                                else printf("%g\t",matrixsim[x][it].R2);
                                if(matrixsim[x][it].ewtest < -9999) printf("na\n");
                                else printf("%g\n",matrixsim[x][it].ewtest);
                            }
                            printf("\n\n");
                        }*/
                        if(file_output2) {
                            /*print to file*/
                            fputs("\n\n Display detailed neutral tests: \n\n Multiple hits not included in the analysis.\n",file_output2);
                            fputs(" The table contains the factor correction for chromosome population size, the number of samples for each locus,\n",file_output2);
                            fputs(" Tajima's D test, Fu and Li's D and F tests, normalized Fay and Wu H test, Fu's Fs test,\n",file_output2);
                            fputs(" Rozas' et a. ZA test, Wall's B and Q tests, Ramos-Onsins & Rozas (R2) test and Ewens-Watterson test.\n",file_output2);
                            fputs(" 'na' means the test is not calculated because the number of samples or the number of segregating sites\n",file_output2);
                            fputs("were insufficient.\n",file_output2);
                            
                            for(x=0;x<data[0].n_loci;x++) {
								fprintf(file_output2,"\nfactor_chrn(%d): %g\n",x,inputms[x].factor_chrn);
                                fprintf(file_output2,
                                    "\n#samples(%d): %d\t#sites(%d): %ld\n",
                                    x,inputms[x].nsam,x,inputms[x].nsites);
                                fprintf(file_output2,"\n#iter\t");
								/*fprintf(file_output2,"#sites_nomh(%d)\t",x);*/
                                fprintf(file_output2,"TajD(%d)\tFuLiD*(%d)\tFuLiF*(%d)\tFuFs(%d)\t",x,x,x,x);
                                fprintf(file_output2,"RozasZA(%d)\tWallB(%d)\tWallQ(%d)\tR2(%d)\tEW(%d)\n",x,x,x,x,x);    
                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                    fprintf(file_output2,"%ld\t",it);
									/*fprintf(file_output2,"%ld\t",matrixsim[x][it].nsites);*/
                                    if(matrixsim[x][it].tajimaD < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].tajimaD);
                                    if(matrixsim[x][it].fuliDn < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].fuliDn);
                                    if(matrixsim[x][it].fuliFn < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].fuliFn);
                                    if(matrixsim[x][it].fuFs == (double) -10000) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].fuFs);
                                    if(matrixsim[x][it].rZA < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].rZA);
                                    if(matrixsim[x][it].wB < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].wB);
                                    if(matrixsim[x][it].wQ < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].wQ);
                                    if(matrixsim[x][it].R2 < -9999) fputs("na\t",file_output2);
                                    else fprintf(file_output2,"%g\t",matrixsim[x][it].R2);
                                    if(matrixsim[x][it].ewtest < -9999) fputs("na\n",file_output2);
                                    else fprintf(file_output2,"%g\n",matrixsim[x][it].ewtest);
                                }
                                fputs("\n\n",file_output2);
                            }
                        }
                    }
					fclose(file_output2);
                    break;
                default: /*error*/
                    return;
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
                    printf("\n     3.1.0.0.0. Display table with detailed statistics for multilocus data per iteration:\n\n");
                    
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
                        fprintf(file_output,"\n     3.1.0.0.0. Display table with detailed statistics for multilocus data per iteration:\n\n");
                        
                        fprintf(file_output," 0 - Display general statistics.\n");
                        fprintf(file_output," 1 - Display other linkage related statistics.\n");
                        fprintf(file_output," 2 - Display estimates of total locus variability.\n");
                        fprintf(file_output," 3 - Display estimates of nucleotide locus variability.\n");
                        fprintf(file_output," 4 - Back to previous menu.\n");
                    }
                    #endif
                    
                    do *k = getchar();
                    while(*k<'0' || *k>'4');
                    
                    if(file_output) fprintf(file_output,"OPTION CHOSEN: %c\n\n",*k);

					if(*k<'4') {
						printf("\n\n This information will be sent into a new file. Please introduce the name of the file: ");
						scanf("%s",fileout2);

						if((file_output2 = fopen (fileout2,"w")) == 0) {
							printf("\n  It is not possible to open the file %s",fileout2);
							if(file_output) fprintf(file_output,"\n  It is not possible to open the file %s",fileout2);
							return;
						}
						if(file_output) fprintf(file_output," Data results are sent to file %s.\n",fileout2);
					}
                    switch(*k) {
                        case '0': /* 0 - Display general statistics*/                            
							if(*outgroup) {
                                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
                                /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                /*
								printf(" Total number of loci: %d\n\n",data[0].n_loci);
                                it = 0;
                                for(x=0;x<data[0].n_loci;x++) it += (long int)inputms[x].nsites;
                                printf(" Total number of sites: %ld\n",it);
                                printf(" Average number of sites per loci: %g\n\n",(double)it/(double)data[0].n_loci);
                                if(data[0].n_loci > 2) printf(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) for \n");
                                else printf(" Table containing in order, the observed Total (Tot) and average (Avg) for \n");
                                printf(" biallelic sites,");
                                if(data[0].time_spec > 0.) 
                                    printf(" fixed positions,\ndivergence (uncorrected and Jukes and Cantor correction)\n");
                                printf( "haplotypes (excluded shared polymorphisms, if observed) haplotypes/sample size, haplotype diversity \n ");                                                 
                                printf( "and the recombination parameter C (Hudson 1987) and Rm.\n ");                                                 
                                if(data[0].n_loci > 2) {
                                    printf("\nniter\tTotS\tAvgS\tVarS\t");
                                    if(data[0].time_spec > 0.) {
                                        printf("TotFix\tAvgFix\tVarFix\t");
                                        printf("TotDiv\tAvgDiv\tVarDiv\t");
                                        printf("TotDivjc\tAvgDivjc\tVarDivjc\t");
                                    }
                                    printf("Tothap\tAvghap\tVarhap\t");
                                    printf("Tothapsam\tAvghapsam\tVarhapsam\t");
                                    printf("Tothapdiv\tAvghapdiv\tVarhapdiv\t");
                                    printf("TotC\tAvgC\tVarC\t");
                                    printf("TotRm\tAvgRm\tVarRm\t");
                                }
                                else {
                                    printf("\nniter\tTotS\tAvgS\t");
                                    if(data[0].time_spec > 0.) {
                                        printf("TotFix\tAvgFix\t");
                                        printf("TotDiv\tAvgDiv\t");
                                        printf("TotDivjc\tAvgDivjc\t");
                                    }
                                    printf("Tothap\tAvghap\t");
                                    printf("Tothapsam\tAvghapsam\t");
                                    printf("Tothapdiv\tAvghapdiv\t");
                                    printf("TotC\tAvgC\tVarC\t");
                                    printf("TotRm\tAvgRm\tVarRm\t");
                                }
                                printf("\n");
                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                    printf("%ld\t",it);
                                    printf("%ld\t",matrixmlsim[it].Sbiallsites);
                                    printf("%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    if(data[0].time_spec > 0.) {
                                        printf("%ld\t",matrixmlsim[it].Sfixed);
                                        printf("%g\t",(double)matrixmlsim[it].Sfixed/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2fixed,(double)matrixmlsim[it].Sfixed,data[0].n_loci);
                                        if(var > (double)-10000) printf("%g\t",var);
                                        printf("%g\t",matrixmlsim[it].Sndivergence);
                                        printf("%g\t",(double)matrixmlsim[it].Sndivergence/(double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2ndivergence,(double)matrixmlsim[it].Sndivergence,data[0].n_loci);
                                        if(var > (double)-10000) printf("%g\t",var);
                                        printf("%g\t",matrixmlsim[it].Sndivergencejc);
                                        printf("%g\t",(double)matrixmlsim[it].Sndivergencejc/(double)matrixmlsim[it].nldivjc);
                                        var = variances((double)matrixmlsim[it].S2ndivergencejc,(double)matrixmlsim[it].Sndivergencejc,matrixmlsim[it].nldivjc);
                                        if(var > (double)-10000) printf("%g\t",var);
                                    }
                                    printf("%ld\t",matrixmlsim[it].Snhapl);
                                    printf("%g\t",(double)matrixmlsim[it].Snhapl/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2nhapl,(double)matrixmlsim[it].Snhapl,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    printf("%ld\t",matrixmlsim[it].Snhaplsam);
                                    printf("%g\t",(double)matrixmlsim[it].Snhaplsam/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2nhaplsam,(double)matrixmlsim[it].Snhaplsam,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    printf("%ld\t",matrixmlsim[it].Shapldiv);
                                    printf("%g\t",(double)matrixmlsim[it].Shapldiv/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2hapldiv,(double)matrixmlsim[it].Shapldiv,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    printf("%g\t",matrixmlsim[it].SRvpi);
                                    printf("%g\t",(double)matrixmlsim[it].SRvpi/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2Rvpi,(double)matrixmlsim[it].SRvpi,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    printf("%ld\t",matrixmlsim[it].SRm);
                                    printf("%g\t",(double)matrixmlsim[it].SRm/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2Rm,(double)matrixmlsim[it].SRm,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    else printf("\n");
                                }
                                */
                                
                                /* print the median/average for all statistic in each loci.*/
                                if(onlymulo)
                                    printf("\n Value of the average obtained by simulations for statistics in each locus:\n\n");
                                else 
                                    printf("\n Value of the median obtained by simulations for statistics in each locus:\n\n");

                                printf("\nlocus\t");
                                printf("S\t");
                                if(data[0].time_spec > 0.) {
                                    printf("Fixed\t");
                                    printf("Div\t");
                                }
                                printf("nhapl\t");
                                printf("nhaplsam\t");
                                printf("hapldiv\t");
                                printf("C\tRm\n");
                                
                                for(x=0;x<data[0].n_loci;x++) {
                                    if(observed_data && dataobsequalsim) printf("%d:%s\t",x,matrix[x].gene);
                                    else printf("%d\t",x);
                                    printf("%g\t",avgstatloci[x].biallsites);
                                    if(data[0].time_spec > 0.) {
                                        printf("%g\t",avgstatloci[x].fixed);
                                        printf("%g\t",avgstatloci[x].ndivergence);
                                    }
                                    printf("%g\t",avgstatloci[x].nhapl);
                                    printf("%g\t",avgstatloci[x].nhaplsam);
                                    printf("%g\t",avgstatloci[x].hapldiv);
                                    if(avgstatloci[x].Rvpi >= 0.) printf("%g\t",avgstatloci[x].Rvpi);
									else printf("na\t");
                                    printf("%g\n",avgstatloci[x].Rm);
                                }
                                
                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output2);
                                    /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                    fprintf(file_output2," Total number of loci: %d\n\n",data[0].n_loci);
                                    is = (double)0;
                                    for(x=0;x<data[0].n_loci;x++) for(it=0;it<(long int)data[0].n_iter;it++) is += (double)matrixsim[x][it].nsites;
									is /= (double)data[0].n_iter;
                                    fprintf(file_output2," Total number of sites: %ld\n",(long int)is);
                                    fprintf(file_output2," Average number of sites per loci: %g\n\n",(double)is/(double)data[0].n_loci);
                                    if(data[0].n_loci > 2) fputs(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) for \n",file_output2);
                                    else fputs(" Table containing in order, the observed Total (Tot) and average (Avg) for \n",file_output2);
                                    fputs(" biallelic sites,",file_output2);
                                    if(data[0].time_spec > 0.) fputs(" fixed positions,\ndivergence (uncorrected and Jukes and Cantor correction)\n",file_output2);
                                    fputs( "haplotypes (excluded shared polymorphisms, if observed) haplotypes/sample size, haplotype diversity \n ",file_output2);                                
                                    fputs( "and Rm.\n ",file_output2);                               
                                    if(data[0].n_loci > 2) {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\tVarS\t");
                                        if(data[0].time_spec > 0.) {
                                            fprintf(file_output2,"TotFix\tAvgFix\tVarFix\t");
                                            fprintf(file_output2,"TotDiv\tAvgDiv\tVarDiv\t");
                                            fprintf(file_output2,"TotDivjc\tAvgDivjc\tVarDivjc\t");
                                        }
                                        fprintf(file_output2,"Tothap\tAvghap\tVarhap\t");
                                        fprintf(file_output2,"Tothapsam\tAvghapsam\tVarhapsam\t");
                                        fprintf(file_output2,"Tothapdiv\tAvghapdiv\tVarhapdiv\t");
                                        fprintf(file_output2,"TotRm\tAvgRm\tVarRm\t");
                                    }
                                    else {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\t");
                                        if(data[0].time_spec > 0.) {
                                            fprintf(file_output2,"TotFix\tAvgFix\t");
                                            fprintf(file_output2,"TotDiv\tAvgDiv\t");
                                            fprintf(file_output2,"TotDivjc\tAvgDivjc\t");
                                        }
                                        fprintf(file_output2,"Tothap\tAvghap\t");
                                        fprintf(file_output2,"Tothapsam\tAvghapsam\t");
                                        fprintf(file_output2,"Tothapdiv\tAvghapdiv\t");
                                        fprintf(file_output2,"TotRm\tAvgRm\tVarRm\t");
                                    }

                                    fputs("\n",file_output2);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        fprintf(file_output2,"%ld\t",it);
                                        
										fprintf(file_output2,"%ld\t",matrixmlsim[it].Sbiallsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
										else if(data[0].n_loci > 2) fputs("na\t",file_output2);
                                        
										if(data[0].time_spec > 0.) {
                                            fprintf(file_output2,"%ld\t",matrixmlsim[it].Sfixed);
                                            fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sfixed/ (double)data[0].n_loci);
                                            var = variances((double)matrixmlsim[it].S2fixed,(double)matrixmlsim[it].Sfixed,data[0].n_loci);
                                            if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
											else if(data[0].n_loci > 2) fputs("na\t",file_output2);
                                            
											fprintf(file_output2,"%g\t",matrixmlsim[it].Sndivergence);
                                            fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sndivergence/ (double)data[0].n_loci);
                                            var = variances((double)matrixmlsim[it].S2ndivergence,(double)matrixmlsim[it].Sndivergence,data[0].n_loci);
                                            if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
											else if(data[0].n_loci > 2) fputs("na\t",file_output2);
                                            
											fprintf(file_output2,"%g\t",matrixmlsim[it].Sndivergencejc);
                                            fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sndivergencejc/(double)matrixmlsim[it].nldivjc);
                                            var = variances((double)matrixmlsim[it].S2ndivergencejc,(double)matrixmlsim[it].Sndivergencejc,matrixmlsim[it].nldivjc);
                                            if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
											else if(data[0].n_loci > 2) fputs("na\t",file_output2);
                                        }
                                        fprintf(file_output2,"%g\t",matrixmlsim[it].Snhapl);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Snhapl/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2nhapl,(double)matrixmlsim[it].Snhapl,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
										else if(data[0].n_loci > 2) fputs("na\t",file_output2);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Snhaplsam);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Snhaplsam/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2nhaplsam,(double)matrixmlsim[it].Snhaplsam,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
										else if(data[0].n_loci > 2) fputs("na\t",file_output2);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Shapldiv);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Shapldiv/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2hapldiv,(double)matrixmlsim[it].Shapldiv,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
										else if(data[0].n_loci > 2) fputs("na\t",file_output2);
										
										fprintf(file_output2,"%ld\t",matrixmlsim[it].SRm);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].SRm/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2Rm,(double)matrixmlsim[it].SRm,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\n",var);
										else if(data[0].n_loci > 2) fputs("na\n",file_output2);
                                    }
                                    
                                    
                                    /* print the median/average for all statistic in each loci.*/
                                    if(onlymulo)
                                        fputs("\n Value of the average obtained by simulations for statistics in each locus:\n\n",file_output2);
                                    else 
                                        fputs("\n Value of the median obtained by simulations for statistics in each locus:\n\n",file_output2);
    
                                    fprintf(file_output2,"\nlocus\t");
                                    fputs("S\t",file_output2);
                                    if(data[0].time_spec > 0.) {
                                        fputs("Fixed\t",file_output2);
                                        fputs("Div\t",file_output2);
                                    }
                                    fputs("nhapl\t",file_output2);
                                    fputs("nhaplsam\t",file_output2);
                                    fputs("hapldiv\t",file_output2);
                                    fputs("Rm\n",file_output2);
                                    
                                    for(x=0;x<data[0].n_loci;x++) {
                                        if(observed_data && dataobsequalsim) fprintf(file_output2,"%d:%s\t",x,matrix[x].gene);
                                        else fprintf(file_output2,"%d\t",x);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].biallsites);
                                        if(data[0].time_spec > 0.) {
                                            fprintf(file_output2,"%g\t",avgstatloci[x].fixed);
                                            fprintf(file_output2,"%g\t",avgstatloci[x].ndivergence);
                                        }
                                        fprintf(file_output2,"%g\t",avgstatloci[x].nhapl);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].nhaplsam);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].hapldiv);
                                         fprintf(file_output2,"%g\n",avgstatloci[x].Rm);
                                    }
                                }
                            }
                            else {
                                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
                                /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                /*
								printf(" Total number of loci: %d\n\n",data[0].n_loci);
                                it = 0;
                                for(x=0;x<data[0].n_loci;x++) it += (long int)inputms[x].nsites;
                                printf(" Total number of sites: %ld\n",it);
                                printf(" Average number of sites per loci: %g\n\n",(double)it/(double)data[0].n_loci);
                                if(data[0].n_loci > 2) printf(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) for \n");
                                else printf(" Table containing in order, the observed Total (Tot) and average (Avg) for \n");
                                printf(" biallelic sites,");
                                printf(" haplotypes, haplotypes/sample size, haplotype diversity and recombination parameter C (Husdon 1987) and Rm.\n ");                                

                                if(data[0].n_loci > 2) {
                                    printf("\nniter\tTotS\tAvgS\tVarS\t");
                                    printf("Tothap\tAvghap\tVarhap\t");
                                    printf("Tothapsam\tAvghapsam\tVarhapsam\t");
                                    printf("Tothapdiv\tAvghapdiv\tVarhapdiv\t");
                                    printf("TotC\tAvgC\tVarC\t");
                                    printf("TotRm\tAvgRm\tVarRm\t");
                                }
                                else {
                                    printf("\nniter\tTotS\tAvgS\t");
                                    printf("Tothap\tAvghap\t");
                                    printf("Tothapsam\tAvghapsam\t");
                                    printf("Tothapdiv\tAvghapdiv\t");
                                    printf("TotC\tAvgC\tVarC\t");
                                    printf("TotRm\tAvgRm\tVarRm\t");
                                }

                                printf("\n");
                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                    printf("%ld\t",it);
                                    printf("%ld\t",matrixmlsim[it].Sbiallsites);
                                    printf("%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%ld\t",matrixmlsim[it].Snhapl);
                                    printf("%g\t",(double)matrixmlsim[it].Snhapl/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2nhapl,(double)matrixmlsim[it].Snhapl,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    printf("%ld\t",matrixmlsim[it].Snhaplsam);
                                    printf("%g\t",(double)matrixmlsim[it].Snhaplsam/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2nhaplsam,(double)matrixmlsim[it].Snhaplsam,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    printf("%ld\t",matrixmlsim[it].Shapldiv);
                                    printf("%g\t",(double)matrixmlsim[it].Shapldiv/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2hapldiv,(double)matrixmlsim[it].Shapldiv,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    printf("%g\t",matrixmlsim[it].SRvpi);
                                    printf("%g\t",(double)matrixmlsim[it].SRvpi/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2Rvpi,(double)matrixmlsim[it].SRvpi,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
									printf("%ld\t",matrixmlsim[it].SRm);
									printf("%g\t",(double)matrixmlsim[it].SRm/ (double)data[0].n_loci);
									var = variances((double)matrixmlsim[it].S2Rm,(double)matrixmlsim[it].SRm,data[0].n_loci);
									if(var > (double)-10000) printf("%g\n",var);
                                    else printf("\n");
                                }
                                */              
                                /* print the median/average for all statistic in each loci.*/
                                if(onlymulo)
                                    printf("\n Value of the average obtained by simulations for statistics in each locus:\n\n");
                                else 
                                    printf("\n Value of the median obtained by simulations for statistics in each locus:\n\n");

                                printf("\nlocus\t");
                                printf("S\t");
                                printf("nhapl\tnhaplsam\thapldiv\tRm\n");
                                
                                for(x=0;x<data[0].n_loci;x++) {
                                    if(observed_data && dataobsequalsim) printf("%d:%s\t",x,matrix[x].gene);
                                    else printf("%d\t",x);
                                    printf("%g\t",avgstatloci[x].biallsites);
                                    printf("%g\t",avgstatloci[x].nhapl);
                                    printf("%g\t",avgstatloci[x].nhaplsam);
                                    printf("%g\t",avgstatloci[x].hapldiv);
                                    printf("%g\n",avgstatloci[x].Rm);
                                }

                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output2);
                                    /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                    fprintf(file_output2," Total number of loci: %d\n\n",data[0].n_loci);
                                    is = (double)0;
                                    for(x=0;x<data[0].n_loci;x++) for(it=0;it<(long int)data[0].n_iter;it++) is += (double)matrixsim[x][it].nsites;
									is /= (double)data[0].n_iter;
                                    fprintf(file_output2," Total number of sites: %ld\n",(long int)is);
                                    fprintf(file_output2," Average number of sites per loci: %g\n\n",(double)is/(double)data[0].n_loci);
                                    if(data[0].n_loci > 2) fputs(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) for \n", file_output2);
                                    else fputs(" Table containing in order, the observed Total (Tot) and average (Avg) for \n", file_output2);
                                    fputs(" biallelic sites,",file_output2);
                                    fputs(" haplotypes, haplotypes/sample size and haplotype diversity. \n",file_output2);                                
                                    fputs(" and Rm.\n ",file_output2);                               

                                    if(data[0].n_loci > 2) {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\tVarS\t");
                                        fprintf(file_output2,"Tothap\tAvghap\tVarhap\t");
                                        fprintf(file_output2,"Tothapsam\tAvghapsam\tVarhapsam\t");
                                        fprintf(file_output2,"Tothapdiv\tAvghapdiv\tVarhapdiv\t");
                                        fprintf(file_output2,"TotRm\tAvgRm\tVarRm\t");
                                    }
                                    else {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\t");
                                        fprintf(file_output2,"Tothap\tAvghap\t");
                                        fprintf(file_output2,"Tothapsam\tAvghapsam\t");
                                        fprintf(file_output2,"Tothapdiv\tAvghapdiv\t");
                                         fprintf(file_output2,"TotRm\tAvgRm\tVarRm\t");
                                    }

                                    fputs("\n",file_output2);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        fprintf(file_output2,"%ld\t",it);
                                        
										fprintf(file_output2,"%ld\t",matrixmlsim[it].Sbiallsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
										else if(data[0].n_loci > 2) fputs("na\t",file_output2);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Snhapl);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Snhapl/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2nhapl,(double)matrixmlsim[it].Snhapl,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
										else if(data[0].n_loci > 2) fputs("na\t",file_output2);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Snhaplsam);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Snhaplsam/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2nhaplsam,(double)matrixmlsim[it].Snhaplsam,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
										else if(data[0].n_loci > 2) fputs("na\t",file_output2);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Shapldiv);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Shapldiv/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2hapldiv,(double)matrixmlsim[it].Shapldiv,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
										else if(data[0].n_loci > 2) fputs("na\t",file_output2);
										
										fprintf(file_output2,"%ld\t",matrixmlsim[it].SRm);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].SRm/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2Rm,(double)matrixmlsim[it].SRm,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\n",var);
										else if(data[0].n_loci > 2) fputs("na\n",file_output2);
                                    }
                                                
                                    /* print the median/average for all statistic in each loci.*/
                                    if(onlymulo)
                                        fputs("\n Value of the average obtained by simulations for statistics in each locus:\n\n",file_output2);
                                    else 
                                        fputs("\n Value of the median obtained by simulations for statistics in each locus:\n\n",file_output2);
    
                                    fprintf(file_output2,"\nlocus\t");
                                    fputs("S\t",file_output2);
                                    fputs("nhapl\tnhaplsam\thapldiv\tRm\n",file_output2);
                                    
                                    for(x=0;x<data[0].n_loci;x++) {
                                        if(observed_data && dataobsequalsim) fprintf(file_output2,"%d:%s\t",x,matrix[x].gene);
                                        else fprintf(file_output2,"%d\t",x);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].biallsites);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].nhapl);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].nhaplsam);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].hapldiv);
										fprintf(file_output2,"%g\n",avgstatloci[x].Rm);
                                    }    
                                }
                            }
							fclose(file_output2);
                            break;
                        case '1': /* 1 - Display linkage related statistics*/
                            if(*outgroup) {
                                printf("\n\n Display multilocus statistics. Linkage related statistics: \n\n Multiple hits not included in the analysis.\n");
                                /* nombre de loci, then SUM, AVERAGE and VARIANCE */
								/*
                                printf(" Total number of loci: %d\n\n",data[0].n_loci);
                                it = 0;
                                for(x=0;x<data[0].n_loci;x++) it += (long int)inputms[x].nsites;
                                printf(" Total number of sites: %ld\n",it);
                                printf(" Average number of sites per loci: %g\n\n",(double)it/(double)data[0].n_loci);
                                if(data[0].n_loci > 2) printf(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) for \n");
                                else printf(" Table containing in order, the observed Total (Tot) and average (Avg) for \n");
                                printf(" biallelic sites,");
                                printf( " Wall's b and q statistics, and Rozas' et al. za statistic (excluded shared polymorphisms, if observed).\n ");                                
                                if(data[0].n_loci > 2) {
                                    printf("\nniter\tTotS\tAvgS\tVarS\t");
                                    printf("Totb\tAvgb\tVarb\tTotq\tAvgq\tVarq\tTotza\tAvgza\tVarza\t");
                                }
                                else {
                                    printf("\nniter\tTotS\tAvgS\t");
                                    printf("Totb\tAvgb\tTotq\tAvgq\tTotza\tAvgza\t");
                                }

                                printf("\n");
                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                    printf("%ld\t",it);
                                    printf("%ld\t",matrixmlsim[it].Sbiallsites);
                                    printf("%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%ld\t",matrixmlsim[it].Sb);
                                    printf("%g\t",(double)matrixmlsim[it].Sb/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2b,(double)matrixmlsim[it].Sb,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%ld\t",matrixmlsim[it].Sq);
                                    printf("%g\t",(double)matrixmlsim[it].Sq/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2q,(double)matrixmlsim[it].Sq,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Sza);
                                    printf("%g\t",(double)matrixmlsim[it].Sza/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2za,(double)matrixmlsim[it].Sza,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    else printf("\n");
                                }
								*/
                                /* print the median/average for all statistic in each loci.*/
                                if(onlymulo)
                                    printf("\n Value of the average obtained by simulations in each locus:\n\n");
                                else 
                                    printf("\n Value of the median obtained by simulations in each locus:\n\n");

                                printf("\nlocus\t");
                                printf("S\t");
                                printf("b\t");
                                printf("q\t");
                                printf("za\n");
                                
                                for(x=0;x<data[0].n_loci;x++) {
                                    if(observed_data && dataobsequalsim) printf("%d:%s\t",x,matrix[x].gene);
                                    else printf("%d\t",x);
                                    printf("%g\t",avgstatloci[x].biallsites);
                                    printf("%g\t",avgstatloci[x].b);
                                    printf("%g\t",avgstatloci[x].q);
                                    printf("%g\n",avgstatloci[x].za);
                                }

                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display multilocus statistics. Linkage related statistics: \n\n Multiple hits not included in the analysis.\n",file_output2);
                                    /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                    fprintf(file_output2," Total number of loci: %d\n\n",data[0].n_loci);
                                    is = (double)0;
                                    for(x=0;x<data[0].n_loci;x++) for(it=0;it<(long int)data[0].n_iter;it++) is += (double)matrixsim[x][it].nsites;
									is /= (double)data[0].n_iter;
                                    fprintf(file_output2," Total number of sites: %ld\n",(long int)is);
                                    fprintf(file_output2," Average number of sites per loci: %g\n\n",(double)is/(double)data[0].n_loci);
                                    if(data[0].n_loci > 2) fputs(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) for \n", file_output2);
                                    else fputs(" Table containing in order, the observed Total (Tot) and average (Avg) for \n", file_output2);
                                    fputs(" biallelic sites,",file_output2);
                                    fputs( " Wall's b and q statistics, and Rozas' et al. za statistic (excluded shared polymorphisms, if observed).\n ",file_output2);                                
                                    if(data[0].n_loci > 2) {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\tVarS\t");
                                        fprintf(file_output2,"Totb\tAvgb\tVarb\tTotq\tAvgq\tVarq\tTotza\tAvgza\tVarza\t");
                                    }
                                    else {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\t");
                                        fprintf(file_output2,"Totb\tAvgb\tTotq\tAvgq\tTotza\tAvgza\t");
                                    }

                                    fputs("\n",file_output2);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        fprintf(file_output2,"%ld\t",it);
                                        
										fprintf(file_output2,"%ld\t",matrixmlsim[it].Sbiallsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%ld\t",matrixmlsim[it].Sb);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sb/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2b,(double)matrixmlsim[it].Sb,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%ld\t",matrixmlsim[it].Sq);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sq/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2q,(double)matrixmlsim[it].Sq,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Sza);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sza/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2za,(double)matrixmlsim[it].Sza,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\n",var);
                                        else fputs("\n",file_output2);
                                    }
                                                    
                                    /* print the median/average for all statistic in each loci.*/
                                    if(onlymulo)
                                        fputs("\n Value of the average obtained by simulations in each locus:\n\n",file_output2);
                                    else 
                                        fputs("\n Value of the median obtained by simulations in each locus:\n\n",file_output2);
    
                                    fprintf(file_output2,"\nlocus\t");
                                    fputs("S\t",file_output2);
                                    fputs("b\t",file_output2);
                                    fputs("q\t",file_output2);
                                    fputs("za\n",file_output2);
                                    
                                    for(x=0;x<data[0].n_loci;x++) {
                                        if(observed_data && dataobsequalsim) fprintf(file_output2,"%d:%s\t",x,matrix[x].gene);
                                        else fprintf(file_output2,"%d\t",x);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].biallsites);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].b);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].q);
                                        fprintf(file_output2,"%g\n",avgstatloci[x].za);
                                    }
    
                                }
                            }
                            else {
                                printf("\n\n Display multilocus statistics. Linkage related statistics: \n\n Multiple hits not included in the analysis.\n");
                                /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                /*
								printf(" Total number of loci: %d\n\n",data[0].n_loci);
                                it = 0;
                                for(x=0;x<data[0].n_loci;x++) it += (long int)inputms[x].nsites;
                                printf(" Total number of sites: %ld\n",it);
                                printf(" Average number of sites per loci: %g\n\n",(double)it/(double)data[0].n_loci);
                                if(data[0].n_loci > 2) printf(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) for \n");
                                else printf(" Table containing in order, the observed Total (Tot) and average (Avg) for \n");
                                printf(" biallelic sites,");
                                printf( " Wall's b and q statistics, and Rozas' et al. za statistic.\n ");                                

                                if(data[0].n_loci > 2)  {
                                    printf("\nniter\tTotS\tAvgS\tVarS\t");
                                    printf("Totb\tAvgb\tVarb\tTotq\tAvgq\tVarq\tTotza\tAvgza\tVarza\t");
                                }
                                else {
                                    printf("\nniter\tTotS\tAvgS\t");
                                    printf("Totb\tAvgb\tTotq\tAvgq\tTotza\tAvgza\t");
                                }

                                printf("\n");
                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                    printf("%ld\t",it);
                                    printf("%ld\t",matrixmlsim[it].Sbiallsites);
                                    printf("%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%ld\t",matrixmlsim[it].Sb);
                                    printf("%g\t",(double)matrixmlsim[it].Sb/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2b,(double)matrixmlsim[it].Sb,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%ld\t",matrixmlsim[it].Sq);
                                    printf("%g\t",(double)matrixmlsim[it].Sq/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2q,(double)matrixmlsim[it].Sq,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Sza);
                                    printf("%g\t",(double)matrixmlsim[it].Sza/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2za,(double)matrixmlsim[it].Sza,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    else printf("\n");
                                }                                
								*/
                                /* print the median/average for all statistic in each loci.*/
                                if(onlymulo)
                                    printf("\n Value of the average obtained by simulations in each locus:\n\n");
                                else 
                                    printf("\n Value of the median obtained by simulations in each locus:\n\n");

                                printf("\nlocus\t");
                                printf("S\t");
                                printf("b\t");
                                printf("q\t");
                                printf("za\n");
                                
                                for(x=0;x<data[0].n_loci;x++) {
                                    if(observed_data && dataobsequalsim) printf("%d:%s\t",x,matrix[x].gene);
                                    else printf("%d\t",x);
                                    printf("%g\t",avgstatloci[x].biallsites);
                                    printf("%g\t",avgstatloci[x].b);
                                    printf("%g\t",avgstatloci[x].q);
                                    printf("%g\n",avgstatloci[x].za);
                                }

                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display multilocus statistics. Linkage related statistics: \n\n Multiple hits not included in the analysis.\n",file_output2);
                                    /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                    fprintf(file_output2," Total number of loci: %d\n\n",data[0].n_loci);
                                    is = (double)0;
                                    for(x=0;x<data[0].n_loci;x++) for(it=0;it<(long int)data[0].n_iter;it++) is += (double)matrixsim[x][it].nsites;
									is /= (double)data[0].n_iter;
                                    fprintf(file_output2," Total number of sites: %ld\n",(long int)is);
                                    fprintf(file_output2," Average number of sites per loci: %g\n\n",(double)is/(double)data[0].n_loci);
                                    if(data[0].n_loci > 2) fputs(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) for \n", file_output2);
                                    else fputs(" Table containing in order, the observed Total (Tot) and average (Avg) for \n", file_output2);
                                    fputs(" biallelic sites,",file_output2);
                                    fputs( " Wall's b and q statistics, and Rozas' et al. za statistic.\n ",file_output2);                                

                                    if(data[0].n_loci > 2) {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\tVarS\t");
                                        fprintf(file_output2,"Totb\tAvgb\tVarb\tTotq\tAvgq\tVarq\tTotza\tAvgza\tVarza\t");
                                    }
                                    else {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\t");
                                        fprintf(file_output2,"Totb\tAvgb\tTotq\tAvgq\tTotza\tAvgza\t");
                                    }

                                    fputs("\n",file_output2);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        fprintf(file_output2,"%ld\t",it);
                                        
										fprintf(file_output2,"%ld\t",matrixmlsim[it].Sbiallsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%ld\t",matrixmlsim[it].Sb);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sb/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2b,(double)matrixmlsim[it].Sb,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%ld\t",matrixmlsim[it].Sq);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sq/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2q,(double)matrixmlsim[it].Sq,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Sza);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sza/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2za,(double)matrixmlsim[it].Sza,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\n",var);
                                        else fputs("\n",file_output2);
                                    }                                
                                                    
                                    /* print the median/average for all statistic in each loci.*/
                                    if(onlymulo)
                                        fputs("\n Value of the average obtained by simulations in each locus:\n\n",file_output2);
                                    else 
                                        fputs("\n Value of the median obtained by simulations in each locus:\n\n",file_output2);
    
                                    fprintf(file_output2,"\nlocus\t");
                                    fputs("S\t",file_output2);
                                    fputs("b\t",file_output2);
                                    fputs("q\t",file_output2);
                                    fputs("za\n",file_output2);
                                    
                                    for(x=0;x<data[0].n_loci;x++) {
                                        if(observed_data && dataobsequalsim) fprintf(file_output2,"%d:%s\t",x,matrix[x].gene);
                                        else fprintf(file_output2,"%d\t",x);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].biallsites);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].b);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].q);
                                        fprintf(file_output2,"%g\n",avgstatloci[x].za);
                                    }
                                }
                            }
							fclose(file_output2);
                            break;
                        case '2': /* 2 - Display estimates of total locus variability*/
                            if(*outgroup) {
                                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
                                /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                /*
								printf(" Total number of loci: %d\n\n",data[0].n_loci);
                                it = 0;
                                for(x=0;x<data[0].n_loci;x++) it += (long int)inputms[x].nsites;
                                printf(" Total number of sites: %ld\n",it);
                                printf(" Average number of sites per loci: %g\n\n",(double)it/(double)data[0].n_loci);
                                if(data[0].n_loci > 2) printf(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) per locus for \n");
                                else printf(" Table containing in order, the observed Total (Tot) and average (Avg) for \n");
                                printf(" biallelic sites,");
                                if(data[0].time_spec > 0.) printf(" averaged divergence (uncorrected and corrected by Jukes and Cantor), \n ");
                                printf( "and estimations of variability from Watterson, Tajima, Fu and Li,\n ");                                
                                printf(" Fay and Wu and Zeng E.\n");

                                if(data[0].n_loci > 2) {
                                    printf("\nniter\tTotS\tAvgS\tVarS\t");
                                    if(data[0].time_spec > 0.) {
                                        printf("TotDiv\tAvgDiv\tVarDiv\tTotDivJC\tAvgDivJC\tVarDivJC\t");
                                    }
                                    printf("TotTheta_wat\tAvgTheta_wat\tVarTheta_wat\t");
                                    printf("TotTheta_taj\tAvgTheta_taj\tVarTheta_taj\t");
                                    printf("TotTheta_fuli\tAvgTheta_fuli\tVarTheta_fuli\t");
                                    printf("TotTheta_fw\tAvgTheta_fw\tVarTheta_fw\t");
                                    printf("TotTheta_Z\tAvgTheta_Z\tVarTheta_Z\t");
                                }
                                else {
                                    printf("\nniter\tTotS\tAvgS\t");
                                    if(data[0].time_spec > 0.) {
                                        printf("TotDiv\tAvgDiv\tTotDivJC\tAvgDivJC\t");
                                    }
                                    printf("TotTheta_wat\tAvgTheta_wat\t");
                                    printf("TotTheta_taj\tAvgTheta_taj\t");
                                    printf("TotTheta_fuli\tAvgTheta_fuli\t");
                                    printf("TotTheta_fw\tAvgTheta_fw\t");
                                    printf("TotTheta_Z\tAvgTheta_Z\t");
                                }

                                printf("\n");
                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                    printf("%ld\t",it);
                                    printf("%ld\t",matrixmlsim[it].Sbiallsites);
                                    printf("%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    if(data[0].time_spec > 0.) {
                                        printf("%g\t",matrixmlsim[it].Sndivergence);
                                        printf("%g\t",(double)matrixmlsim[it].Sndivergence/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2ndivergence,(double)matrixmlsim[it].Sndivergence,data[0].n_loci);
                                        if(var > (double)-10000) printf("%g\t",var);
                                        printf("%g\t",matrixmlsim[it].Sndivergencejc);
                                        printf("%g\t",(double)matrixmlsim[it].Sndivergencejc/(double)matrixmlsim[it].nldivjc);
                                        var = variances((double)matrixmlsim[it].S2ndivergencejc,(double)matrixmlsim[it].Sndivergencejc,matrixmlsim[it].nldivjc);
                                        if(var > (double)-10000) printf("%g\t",var);
                                    }
                                    printf("%g\t",matrixmlsim[it].Stheta_wat);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_wat/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_wat,(double)matrixmlsim[it].Stheta_wat,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_taj);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_taj/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_taj,(double)matrixmlsim[it].Stheta_taj,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_fuli);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_fuli/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_fuli,(double)matrixmlsim[it].Stheta_fuli,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_fw);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_fw/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_fw,(double)matrixmlsim[it].Stheta_fw,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_L);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_L/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_L,(double)matrixmlsim[it].Stheta_L,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    else printf("\n");
                                }                                
								*/
                                /* print the median/average for all statistic in each loci.*/
                                if(onlymulo)
                                    printf("\n Value of the average obtained by simulations in each locus:\n\n");
                                else 
                                    printf("\n Value of the median obtained by simulations in each locus:\n\n");

                                printf("\nlocus\t");
                                printf("S\t");
                                if(data[0].time_spec > 0.) printf("Div\tDivJC\t");
                                printf("Theta_wat\t");
                                printf("Theta_taj\t");
                                printf("Theta_fuli\t");
                                printf("Theta_fw\t");
                                printf("Theta_Z\n");
                                
                                for(x=0;x<data[0].n_loci;x++) {
                                    if(observed_data && dataobsequalsim) printf("%d:%s\t",x,matrix[x].gene);
                                    else printf("%d\t",x);
                                    printf("%g\t",avgstatloci[x].biallsites);
                                    if(data[0].time_spec > 0.) {
                                        printf("%g\t",avgstatloci[x].fixed);
                                        printf("%g\t",avgstatloci[x].ndivergence);
                                    }
                                    printf("%g\t",avgstatloci[x].theta_wat);
                                    printf("%g\t",avgstatloci[x].theta_taj);
                                    printf("%g\t",avgstatloci[x].theta_fuli);
                                    printf("%g\t",avgstatloci[x].theta_fw);
                                    printf("%g\n",avgstatloci[x].theta_L);
                                }
                                
                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output2);
                                    /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                    fprintf(file_output2," Total number of loci: %d\n\n",data[0].n_loci);
                                    is = (double)0;
                                    for(x=0;x<data[0].n_loci;x++) for(it=0;it<(long int)data[0].n_iter;it++) is += (double)matrixsim[x][it].nsites;
									is /= (double)data[0].n_iter;
                                    fprintf(file_output2," Total number of sites: %ld\n",(long int)is);
                                    fprintf(file_output2," Average number of sites per loci: %g\n\n",(double)is/(double)data[0].n_loci);
                                    if(data[0].n_loci > 2) fputs(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) per locus for \n",file_output2);
                                    else fputs(" Table containing in order, the observed Total (Tot) and average (Avg) for \n",file_output2);
                                    fputs(" biallelic sites,",file_output2);
                                    if(data[0].time_spec > 0.) 
                                        fputs(" averaged divergence (uncorrected and corrected by Jukes and Cantor), \n ",file_output2);
                                    fputs( "and estimations of variability from Watterson, Tajima, Fu and Li,\n ",file_output2);                                
                                    fputs(" Fay and Wu and Zeng et al..\n",file_output2);

                                    if(data[0].n_loci > 2) {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\tVarS\t");
                                        if(data[0].time_spec > 0.) {
                                            fprintf(file_output2,"TotDiv\tAvgDiv\tVarDiv\tTotDivJC\tAvgDivJC\tVarDivJC\t");
                                        }
                                        fprintf(file_output2,"TotTheta_wat\tAvgTheta_wat\tVarTheta_wat\t");
                                        fprintf(file_output2,"TotTheta_taj\tAvgTheta_taj\tVarTheta_taj\t");
                                        fprintf(file_output2,"TotTheta_fuli\tAvgTheta_fuli\tVarTheta_fuli\t");
                                        fprintf(file_output2,"TotTheta_fw\tAvgTheta_fw\tVarTheta_fw\t");
                                        fprintf(file_output2,"TotTheta_Z\tAvgTheta_Z\tVarTheta_Z\t");
                                    }
                                    else {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\t");
                                        if(data[0].time_spec > 0.) {
                                            fprintf(file_output2,"TotDiv\tAvgDiv\tTotDivJC\tAvgDivJC\t");
                                        }
                                        fprintf(file_output2,"TotTheta_wat\tAvgTheta_wat\t");
                                        fprintf(file_output2,"TotTheta_taj\tAvgTheta_taj\t");
                                        fprintf(file_output2,"TotTheta_fuli\tAvgTheta_fuli\t");
                                        fprintf(file_output2,"TotTheta_fw\tAvgTheta_fw\t");
                                        fprintf(file_output2,"TotTheta_Z\tAvgTheta_Z");
                                    }

                                    fputs("\n",file_output2);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        fprintf(file_output2,"%ld\t",it);
                                        
										fprintf(file_output2,"%ld\t",matrixmlsim[it].Sbiallsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										if(data[0].time_spec > 0.) {
                                            fprintf(file_output2,"%g\t",matrixmlsim[it].Sndivergence);
                                            fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sndivergence/ (double)data[0].n_loci);
                                            var = variances((double)matrixmlsim[it].S2ndivergence,(double)matrixmlsim[it].Sndivergence,data[0].n_loci);
                                            if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                            
											fprintf(file_output2,"%g\t",matrixmlsim[it].Sndivergencejc);
                                            fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sndivergencejc/ (double)matrixmlsim[it].nldivjc);
                                            var = variances((double)matrixmlsim[it].S2ndivergencejc,(double)matrixmlsim[it].Sndivergencejc,matrixmlsim[it].nldivjc);
                                            if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        }
                                        fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_wat);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_wat/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_wat,(double)matrixmlsim[it].Stheta_wat,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_taj);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_taj/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_taj,(double)matrixmlsim[it].Stheta_taj,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_fuli);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_fuli/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_fuli,(double)matrixmlsim[it].Stheta_fuli,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_fw);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_fw/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_fw,(double)matrixmlsim[it].Stheta_fw,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_L);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_L/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_L,(double)matrixmlsim[it].Stheta_L,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\n",var);
                                        else fputs("\n",file_output2);
                                    }                                
                                                    
                                    /* print the median/average for all statistic in each loci.*/
                                    if(onlymulo)
                                        fputs("\n Value of the average obtained by simulations in each locus:\n\n",file_output2);
                                    else 
                                        fputs("\n Value of the median obtained by simulations in each locus:\n\n",file_output2);
    
                                    fprintf(file_output2,"\nlocus\t");
                                    fputs("S\t",file_output2);
                                    if(data[0].time_spec > 0.) fputs("Div\tDivJC\t",file_output2);
                                    fputs("Theta_wat\t",file_output2);
                                    fputs("Theta_taj\t",file_output2);
                                    fputs("Theta_fuli\t",file_output2);
                                    fputs("Theta_fw\t",file_output2);
                                    fputs("Theta_Z\n",file_output2);
                                    
                                    for(x=0;x<data[0].n_loci;x++) {
                                        if(observed_data && dataobsequalsim) fprintf(file_output2,"%d:%s\t",x,matrix[x].gene);
                                        else fprintf(file_output2,"%d\t",x);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].biallsites);
                                        if(data[0].time_spec > 0.) {
                                            fprintf(file_output2,"%g\t",avgstatloci[x].fixed);
                                            fprintf(file_output2,"%g\t",avgstatloci[x].ndivergence);
                                        }
                                        fprintf(file_output2,"%g\t",avgstatloci[x].theta_wat);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].theta_taj);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].theta_fuli);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].theta_fw);
                                        fprintf(file_output2,"%g\n",avgstatloci[x].theta_L);
                                    }
                                }
                            }
                            else {
                                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
                                /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                /*
								printf(" Total number of loci: %d\n\n",data[0].n_loci);
                                it = 0;
                                for(x=0;x<data[0].n_loci;x++) it += (long int)inputms[x].nsites;
                                printf(" Total number of sites: %ld\n",it);
                                printf(" Average number of sites per loci: %g\n\n",(double)it/(double)data[0].n_loci);
                                if(data[0].n_loci > 2) printf(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) per locus for \n");
                                else printf(" Table containing in order, the observed Total (Tot) and average (Avg) per locus for \n");
                                printf(" biallelic sites,");
                                printf( " and estimations of variability from Watterson, Tajima and Fu and Li.\n ");                                
                                if(data[0].n_loci > 2) {
                                    printf("\nniter\tTotS\tAvgS\tVarS\t");
                                    printf("TotTheta_wat\tAvgTheta_wat\tVarTheta_wat\t");
                                    printf("TotTheta_taj\tAvgTheta_taj\tVarTheta_taj\t");
                                    printf("TotTheta_fuli\tAvgTheta_fuli\tVarTheta_fuli\t");
                                }
                                else {
                                    printf("\nniter\tTotS\tAvgS\t");
                                    printf("TotTheta_wat\tAvgTheta_wat\t");
                                    printf("TotTheta_taj\tAvgTheta_taj\t");
                                    printf("TotTheta_fuli\tAvgTheta_fuli\t");
                                }

                                printf("\n");
                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                    printf("%ld\t",it);
                                    printf("%ld\t",matrixmlsim[it].Sbiallsites);
                                    printf("%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_wat);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_wat/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_wat,(double)matrixmlsim[it].Stheta_wat,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_taj);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_taj/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_taj,(double)matrixmlsim[it].Stheta_taj,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_fulin);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_fulin/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_fulin,(double)matrixmlsim[it].Stheta_fulin,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    else printf("\n");
                                }                                
								*/
                                /* print the median/average for all statistic in each loci.*/
                                if(onlymulo)
                                    printf("\n Value of the average obtained by simulations in each locus:\n\n");
                                else 
                                    printf("\n Value of the median obtained by simulations in each locus:\n\n");

                                printf("\nlocus\t");
                                printf("S\t");
                                printf("Theta_wat\t");
                                printf("Theta_taj\t");
                                printf("Theta_fulin\n");
                                
                                for(x=0;x<data[0].n_loci;x++) {
                                    if(observed_data && dataobsequalsim) printf("%d:%s\t",x,matrix[x].gene);
                                    else printf("%d\t",x);
                                    printf("%g\t",avgstatloci[x].biallsites);
                                    printf("%g\t",avgstatloci[x].theta_wat);
                                    printf("%g\t",avgstatloci[x].theta_taj);
                                    printf("%g\n",avgstatloci[x].theta_fulin);
                                }

                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output2);
                                    /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                    fprintf(file_output2," Total number of loci: %d\n\n",data[0].n_loci);
                                    is = (double)0;
                                    for(x=0;x<data[0].n_loci;x++) for(it=0;it<(long int)data[0].n_iter;it++) is += (double)matrixsim[x][it].nsites;
									is /= (double)data[0].n_iter;
                                    fprintf(file_output2," Total number of sites: %ld\n",(long int)is);
                                    fprintf(file_output2," Average number of sites per loci: %g\n\n",(double)is/(double)data[0].n_loci);
                                    if(data[0].n_loci > 2) fputs(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) per locus for \n",file_output2);
                                    else fputs(" Table containing in order, the observed Total (Tot) and average (Avg) per locus for \n",file_output2);
                                    fputs(" biallelic sites,",file_output2);
                                    fputs( " and estimations of variability from Watterson, Tajima and Fu and Li.\n ",file_output2);                                
                                    if(data[0].n_loci > 2) {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\tVarS\t");
                                        fprintf(file_output2,"TotTheta_wat\tAvgTheta_wat\tVarTheta_wat\t");
                                        fprintf(file_output2,"TotTheta_taj\tAvgTheta_taj\tVarTheta_taj\t");
                                        fprintf(file_output2,"TotTheta_fulin\tAvgTheta_fulin\tVarTheta_fulin\t");
                                    }
                                    else {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\t");
                                        fprintf(file_output2,"TotTheta_wat\tAvgTheta_wat\t");
                                        fprintf(file_output2,"TotTheta_taj\tAvgTheta_taj\t");
                                        fprintf(file_output2,"TotTheta_fulin\tAvgTheta_fulin\t");
                                    }

                                    fputs("\n",file_output2);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        fprintf(file_output2,"%ld\t",it);
                                        
										fprintf(file_output2,"%ld\t",matrixmlsim[it].Sbiallsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_wat);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_wat/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_wat,(double)matrixmlsim[it].Stheta_wat,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_taj);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_taj/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_taj,(double)matrixmlsim[it].Stheta_taj,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_fulin);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_fulin/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_fulin,(double)matrixmlsim[it].Stheta_fulin,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\n",var);
                                        else fputs("\n",file_output2);
                                    }                                
                                                    
                                    /* print the median/average for all statistic in each loci.*/
                                    if(onlymulo)
                                        fputs("\n Value of the average obtained by simulations in each locus:\n\n",file_output2);
                                    else 
                                        fputs("\n Value of the median obtained by simulations in each locus:\n\n",file_output2);
    
                                    fprintf(file_output2,"\nlocus\t");
                                    fputs("S\t",file_output2);
                                    fputs("Theta_wat\t",file_output2);
                                    fputs("Theta_taj\t",file_output2);
                                    fputs("Theta_fuli\n",file_output2);
                                    
                                    for(x=0;x<data[0].n_loci;x++) {
                                        if(observed_data && dataobsequalsim) fprintf(file_output2,"%d:%s\t",x,matrix[x].gene);
                                        else fprintf(file_output2,"%d\t",x);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].biallsites);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].theta_wat);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].theta_taj);
                                        fprintf(file_output2,"%g\n",avgstatloci[x].theta_fuli);
                                    }
                                }
                            }
							fclose(file_output2);
                            break;
                        case '3': /* 3 - Display estimates of nucleotide locus variability*/
                            if(*outgroup) {
                                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
                                /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                /*
								printf(" Total number of loci: %d\n\n",data[0].n_loci);
                                it = 0;
                                for(x=0;x<data[0].n_loci;x++) it += (long int)inputms[x].nsites;
                                printf(" Total number of sites: %ld\n",it);
                                printf(" Average number of sites per loci: %g\n\n",(double)it/(double)data[0].n_loci);
                                if(data[0].n_loci > 2) printf(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) per locus for \n");
                                else printf(" Table containing in order, the observed Total (Tot) and average (Avg) per locus for \n");
                                printf(" biallelic sites,");
                                if(data[0].time_spec > 0.) printf(" nucleotide divergence (uncorrected and corrected by Jukes and Cantor), \n ");
                                printf( "and nucleotide estimations of variability from Watterson, Tajima, Fu and Li,\n ");                                
                                printf(" Fay and Wu and Zeng et al.\n");

                                if(data[0].n_loci > 2) {
                                    printf("\nniter\tTotS\tAvgS\tVarS\t");
                                    if(data[0].time_spec > 0.) {
                                        printf("TotDiv\tAvgDiv\tVarDiv\tTotDivJC\tAvgDivJC\tVarDivJC\t");
                                    }
                                    printf("TotTheta_wat\tAvgTheta_wat\tVarTheta_wat\t");
                                    printf("TotTheta_taj\tAvgTheta_taj\tVarTheta_taj\t");
                                    printf("TotTheta_fuli\tAvgTheta_fuli\tVarTheta_fuli\t");
                                    printf("TotTheta_fw\tAvgTheta_fw\tVarTheta_fw\t");
                                    printf("TotTheta_z\tAvgTheta_z\tVarTheta_z\t");
                                }
                                else {
                                    printf("\nniter\tTotS\tAvgS\t");
                                    if(data[0].time_spec > 0.) {
                                        printf("TotDiv\tAvgDiv\tTotDivJC\tAvgDivJC\t");
                                    }
                                    printf("TotTheta_wat\tAvgTheta_wat\t");
                                    printf("TotTheta_taj\tAvgTheta_taj\t");
                                    printf("TotTheta_fuli\tAvgTheta_fuli\t");
                                    printf("TotTheta_fw\tAvgTheta_fw\t");
                                    printf("TotTheta_z\tAvgTheta_z\tVarTheta_z\t");
                                }

                                printf("\n");
                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                    printf("%ld\t",it);
                                    printf("%ld\t",matrixmlsim[it].Sbiallsites);
                                    printf("%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    if(data[0].time_spec > 0.) {
                                        printf("%g\t",matrixmlsim[it].Sndivergence/(double)matrixmlsim[it].Snsites);
                                        printf("%g\t",(double)matrixmlsim[it].Sndivergence_nut/(double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2ndivergence_nut,(double)matrixmlsim[it].Sndivergence_nut,data[0].n_loci);
                                        if(var > (double)-10000) printf("%g\t",var);
                                        printf("%g\t",matrixmlsim[it].Sndivergencejc/(double)matrixmlsim[it].Snsites);
                                        printf("%g\t",(double)matrixmlsim[it].Sndivergencejc_nut/ (double)matrixmlsim[it].nldivjc);
                                        var = variances((double)matrixmlsim[it].S2ndivergencejc_nut,(double)matrixmlsim[it].Sndivergencejc_nut,matrixmlsim[it].nldivjc);
                                        if(var > (double)-10000) printf("%g\t",var);
                                    }
                                    printf("%g\t",matrixmlsim[it].Stheta_wat/(double)matrixmlsim[it].Snsites);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_wat_nut/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_wat_nut,(double)matrixmlsim[it].Stheta_wat_nut,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_taj/(double)matrixmlsim[it].Snsites);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_taj_nut/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_taj_nut,(double)matrixmlsim[it].Stheta_taj_nut,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_fuli/(double)matrixmlsim[it].Snsites);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_fuli_nut/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_fuli_nut,(double)matrixmlsim[it].Stheta_fuli_nut,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_fw/(double)matrixmlsim[it].Snsites);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_fw_nut/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_fw_nut,(double)matrixmlsim[it].Stheta_fw_nut,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_L/(double)matrixmlsim[it].Snsites);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_L_nut/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_L_nut,(double)matrixmlsim[it].Stheta_L_nut,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    else printf("\n");
                                }
								*/
                                /* print the median/average for all statistic in each loci.*/
                                if(onlymulo)
                                    printf("\n Value of the average obtained by simulations in each locus:\n\n");
                                else 
                                    printf("\n Value of the median obtained by simulations in each locus:\n\n");

                                printf("\nlocus\t");
                                printf("S\t");
                                if(data[0].time_spec > 0.) printf("Fix\tDiv_nt\t");
                                printf("Theta_wat_nt\t");
                                printf("Theta_taj_nt\t");
                                printf("Theta_fuli_nt\t");
                                printf("Theta_fw_nt\t");
                                printf("Theta_Z_nt\n");
                                
                                for(x=0;x<data[0].n_loci;x++) {
                                    if(observed_data && dataobsequalsim) printf("%d:%s\t",x,matrix[x].gene);
                                    else printf("%d\t",x);
                                    printf("%g\t",avgstatloci[x].biallsites);
                                    if(data[0].time_spec > 0.) {
                                        printf("%g\t",avgstatloci[x].fixed);
                                        printf("%g\t",avgstatloci[x].ndivergence_nut);
                                    }
                                    printf("%g\t",avgstatloci[x].theta_wat_nut);
                                    printf("%g\t",avgstatloci[x].theta_taj_nut);
                                    printf("%g\t",avgstatloci[x].theta_fuli_nut);
                                    printf("%g\t",avgstatloci[x].theta_fw_nut);
                                    printf("%g\n",avgstatloci[x].theta_L_nut);
                                }
                                    
                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output2);
                                    /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                    fprintf(file_output2," Total number of loci: %d\n\n",data[0].n_loci);
                                    is = (double)0;
                                    for(x=0;x<data[0].n_loci;x++) for(it=0;it<(long int)data[0].n_iter;it++) is += (double)matrixsim[x][it].nsites;
									is /= (double)data[0].n_iter;
                                    fprintf(file_output2," Total number of sites: %ld\n",(long int)is);
                                    fprintf(file_output2," Average number of sites per loci: %g\n\n",(double)is/(double)data[0].n_loci);
                                    if(data[0].n_loci > 2) fputs(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) per locus for \n",file_output2);
                                    else fputs(" Table containing in order, the observed Total (Tot) and average (Avg) per locus for \n",file_output2);
                                    fputs(" biallelic sites,",file_output2);
                                    if(data[0].time_spec > 0.) fputs(" nucleotide divergence (uncorrected and corrected by Jukes and Cantor), \n ",file_output2);
                                    fputs( "and nucleotide estimations of variability from Watterson, Tajima, Fu and Li,\n ",file_output2);                                
                                    fputs(" Fay and Wu and Zeng et al.\n",file_output2);

                                    if(data[0].n_loci > 2) {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\tVarS\t");
                                        if(data[0].time_spec > 0.) {
                                            fprintf(file_output2,"TotDiv\tAvgDiv\tVarDiv\tTotDivJC\tAvgDivJC\tVarDivJC\t");
                                        }
                                        fprintf(file_output2,"TotTheta_wat\tAvgTheta_wat\tVarTheta_wat\t");
                                        fprintf(file_output2,"TotTheta_taj\tAvgTheta_taj\tVarTheta_taj\t");
                                        fprintf(file_output2,"TotTheta_fuli\tAvgTheta_fuli\tVarTheta_fuli\t");
                                        fprintf(file_output2,"TotTheta_fw\tAvgTheta_fw\tVarTheta_fw\t");
                                        fprintf(file_output2,"TotTheta_z\tAvgTheta_z\tVarTheta_z\t");
                                    }
                                    else {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\t");
                                        if(data[0].time_spec > 0.) {
                                            fprintf(file_output2,"TotDiv\tAvgDiv\tTotDivJC\tAvgDivJC\t");
                                        }
                                        fprintf(file_output2,"TotTheta_wat\tAvgTheta_wat\t");
                                        fprintf(file_output2,"TotTheta_taj\tAvgTheta_taj\t");
                                        fprintf(file_output2,"TotTheta_fuli\tAvgTheta_fuli\t");
                                        fprintf(file_output2,"TotTheta_fw\tAvgTheta_fw\t");
                                        fprintf(file_output2,"TotTheta_z\tAvgTheta_z\t");
                                    }

                                    fputs("\n",file_output2);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        fprintf(file_output2,"%ld\t",it);
                                        
										fprintf(file_output2,"%ld\t",matrixmlsim[it].Sbiallsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										if(data[0].time_spec > 0.) {
                                            fprintf(file_output2,"%g\t",matrixmlsim[it].Sndivergence/(double)matrixmlsim[it].Snsites);
                                            fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sndivergence_nut/ (double)data[0].n_loci);
                                            var = variances((double)matrixmlsim[it].S2ndivergence_nut,(double)matrixmlsim[it].Sndivergence_nut,data[0].n_loci);
                                            if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                            
											fprintf(file_output2,"%g\t",matrixmlsim[it].Sndivergencejc/(double)matrixmlsim[it].Snsites);
                                            fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sndivergencejc_nut/ (double)matrixmlsim[it].nldivjc);
                                            var = variances((double)matrixmlsim[it].S2ndivergencejc_nut,(double)matrixmlsim[it].Sndivergencejc_nut,matrixmlsim[it].nldivjc);
                                            if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        }
                                        fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_wat/(double)matrixmlsim[it].Snsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_wat_nut/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_wat_nut,(double)matrixmlsim[it].Stheta_wat_nut,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_taj/(double)matrixmlsim[it].Snsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_taj_nut/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_taj_nut,(double)matrixmlsim[it].Stheta_taj_nut,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_fuli/(double)matrixmlsim[it].Snsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_fuli_nut/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_fuli_nut,(double)matrixmlsim[it].Stheta_fuli_nut,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_fw/(double)matrixmlsim[it].Snsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_fw_nut/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_fw_nut,(double)matrixmlsim[it].Stheta_fw_nut,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_L/(double)matrixmlsim[it].Snsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_L_nut/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_L_nut,(double)matrixmlsim[it].Stheta_L_nut,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\n",var);
                                        else fputs("\n",file_output2);
                                    }
                                    /*print median/averages*/
                                    if(onlymulo)
                                        fputs("\n Value of the average obtained by simulations in each locus:\n\n",file_output2);
                                    else 
                                        fputs("\n Value of the median obtained by simulations in each locus:\n\n",file_output2);
    
                                    fprintf(file_output2,"\nlocus\t");
                                    fputs("S\t",file_output2);
                                    if(data[0].time_spec > 0.) fputs("Fix\tDiv\t",file_output2);
                                    fputs("Theta_wat_nt\t",file_output2);
                                    fputs("Theta_taj_nt\t",file_output2);
                                    fputs("Theta_fuli_nt\t",file_output2);
                                    fputs("Theta_fw_nt\t",file_output2);
                                    fputs("Theta_Z_nt\n",file_output2);
                                    
                                    for(x=0;x<data[0].n_loci;x++) {
                                        if(observed_data && dataobsequalsim) fprintf(file_output2,"%d:%s\t",x,matrix[x].gene);
                                        else fprintf(file_output2,"%d\t",x);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].biallsites);
                                        if(data[0].time_spec > 0.) {
                                            fprintf(file_output2,"%g\t",avgstatloci[x].fixed);
                                            fprintf(file_output2,"%g\t",avgstatloci[x].ndivergence_nut);
                                        }
                                        fprintf(file_output2,"%g\t",avgstatloci[x].theta_wat_nut);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].theta_taj_nut);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].theta_fuli_nut);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].theta_fw_nut);
                                        fprintf(file_output2,"%g\n",avgstatloci[x].theta_L_nut);
                                    }
                                }
                            }
                            else {
                                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
                                /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                /*
								printf(" Total number of loci: %d\n\n",data[0].n_loci);
                                it = 0;
                                for(x=0;x<data[0].n_loci;x++) it += (long int)inputms[x].nsites;
                                printf(" Total number of sites: %ld\n",it);
                                printf(" Average number of sites per loci: %g\n\n",(double)it/(double)data[0].n_loci);
                                if(data[0].n_loci > 2) printf(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) per locus for \n");
                                else printf(" Table containing in order, the observed Total (Tot) and average (Avg) per locus for \n");
                                printf(" biallelic sites,");
                                printf( " and nucleotide estimations of variability from Watterson, Tajima and Fu and Li.\n ");                                

                                if(data[0].n_loci > 2) {
                                    printf("\nniter\tTotS\tAvgS\tVarS\t");
                                    printf("TotTheta_wat\tAvgTheta_wat\tVarTheta_wat\t");
                                    printf("TotTheta_taj\tAvgTheta_taj\tVarTheta_taj\t");
                                    printf("TotTheta_fuli\tAvgTheta_fuli\tVarTheta_fuli\t");
                                }
                                else {
                                    printf("\nniter\tTotS\tAvgS\t");
                                    printf("TotTheta_wat\tAvgTheta_wat\t");
                                    printf("TotTheta_taj\tAvgTheta_taj\t");
                                    printf("TotTheta_fuli\tAvgTheta_fuli\t");
                                }

                                printf("\n");
                                for(it=0;it<(long int)data[0].n_iter;it++) {
                                    printf("%ld\t",it);
                                    printf("%ld\t",matrixmlsim[it].Sbiallsites);
                                    printf("%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_wat/(double)matrixmlsim[it].Snsites);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_wat_nut/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_wat_nut,(double)matrixmlsim[it].Stheta_wat_nut,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_taj/(double)matrixmlsim[it].Snsites);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_taj_nut/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_taj_nut,(double)matrixmlsim[it].Stheta_taj_nut,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\t",var);
                                    printf("%g\t",matrixmlsim[it].Stheta_fulin/(double)matrixmlsim[it].Snsites);
                                    printf("%g\t",(double)matrixmlsim[it].Stheta_fulin_nut/ (double)data[0].n_loci);
                                    var = variances((double)matrixmlsim[it].S2theta_fulin_nut,(double)matrixmlsim[it].Stheta_fulin_nut,data[0].n_loci);
                                    if(var > (double)-10000) printf("%g\n",var);
                                    else printf("\n");
                                }
								*/
                                /* print the median for all statistic in each loci.*/
                                if(onlymulo)
                                    printf("\n Value of the average obtained by simulations in each locus:\n\n");
                                else 
                                    printf("\n Value of the median obtained by simulations in each locus:\n\n");

                                printf("\nlocus\t");
                                printf("S\t");
                                printf("Theta_wat_nt\t");
                                printf("Theta_taj_nt\t");
                                printf("Theta_fulin_nt\n");
                                
                                for(x=0;x<data[0].n_loci;x++) {
                                    if(observed_data && dataobsequalsim) printf("%d:%s\t",x,matrix[x].gene);
                                    else printf("%d\t",x);
                                    printf("%g\t",avgstatloci[x].biallsites);
                                    printf("%g\t",avgstatloci[x].theta_wat_nut);
                                    printf("%g\t",avgstatloci[x].theta_taj_nut);
                                    printf("%g\n",avgstatloci[x].theta_fulin_nut);
                                }
                                    
                                if(file_output2) {
                                    /*print to file*/
                                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output2);
                                    /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                                    fprintf(file_output2," Total number of loci: %d\n\n",data[0].n_loci);
                                    is = (double)0;
                                    for(x=0;x<data[0].n_loci;x++) for(it=0;it<(long int)data[0].n_iter;it++) is += (double)matrixsim[x][it].nsites;
									is /= (double)data[0].n_iter;
                                    fprintf(file_output2," Total number of sites: %ld\n",(long int)is);
                                    fprintf(file_output2," Average number of sites per loci: %g\n\n",(double)is/(double)data[0].n_loci);
                                    if(data[0].n_loci > 2) fputs(" Table containing in order, the observed Total (Tot), average (Avg) and variance (Var) per locus for \n",file_output2);
                                    else fputs(" Table containing in order, the observed Total (Tot) and average (Avg) per locus for \n",file_output2);
                                    fputs(" biallelic sites,",file_output2);
                                    fputs( " and nucleotide estimations of variability from Watterson, Tajima and Fu and Li.\n ",file_output2);                                
                                    if(data[0].n_loci > 2) {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\tVarS\t");
                                        fprintf(file_output2,"TotTheta_wat\tAvgTheta_wat\tVarTheta_wat\t");
                                        fprintf(file_output2,"TotTheta_taj\tAvgTheta_taj\tVarTheta_taj\t");
                                        fprintf(file_output2,"TotTheta_fulin\tAvgTheta_fulin\tVarTheta_fulin\t");
                                    }
                                    else {
                                        fprintf(file_output2,"\nniter\tTotS\tAvgS\t");
                                        fprintf(file_output2,"TotTheta_wat\tAvgTheta_wat\t");
                                        fprintf(file_output2,"TotTheta_taj\tAvgTheta_taj\t");
                                        fprintf(file_output2,"TotTheta_fulin\tAvgTheta_fulin\t");
                                    }

                                    fputs("\n",file_output2);
                                    for(it=0;it<(long int)data[0].n_iter;it++) {
                                        fprintf(file_output2,"%ld\t",it);
                                        
										fprintf(file_output2,"%ld\t",matrixmlsim[it].Sbiallsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Sbiallsites/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2biallsites,(double)matrixmlsim[it].Sbiallsites,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_wat/(double)matrixmlsim[it].Snsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_wat_nut/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_wat_nut,(double)matrixmlsim[it].Stheta_wat_nut,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_taj/(double)matrixmlsim[it].Snsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_taj_nut/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_taj_nut,(double)matrixmlsim[it].Stheta_taj_nut,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\t",var);
                                        
										fprintf(file_output2,"%g\t",matrixmlsim[it].Stheta_fulin/(double)matrixmlsim[it].Snsites);
                                        fprintf(file_output2,"%g\t",(double)matrixmlsim[it].Stheta_fulin_nut/ (double)data[0].n_loci);
                                        var = variances((double)matrixmlsim[it].S2theta_fulin_nut,(double)matrixmlsim[it].Stheta_fulin_nut,data[0].n_loci);
                                        if(var > (double)-10000) fprintf(file_output2,"%g\n",var);
                                        else fputs("\n",file_output2);
                                    }
                                    /*print the median/average*/
                                    if(onlymulo)
                                        fputs("\n Value of the average obtained by simulations in each locus:\n\n",file_output2);
                                    else 
                                        fputs("\n Value of the median obtained by simulations in each locus:\n\n",file_output2);
    
                                    fprintf(file_output2,"\nlocus\t");
                                    fputs("S\t",file_output2);
                                    fputs("Theta_wat_nt\t",file_output2);
                                    fputs("Theta_taj_nt\t",file_output2);
                                    fputs("Theta_fulin_nt\n",file_output2);
                                    
                                    for(x=0;x<data[0].n_loci;x++) {
                                        if(observed_data && dataobsequalsim) fprintf(file_output2,"%d:%s\t",x,matrix[x].gene);
                                        else fprintf(file_output2,"%d\t",x);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].biallsites);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].theta_wat_nut);
                                        fprintf(file_output2,"%g\t",avgstatloci[x].theta_taj_nut);
                                        fprintf(file_output2,"%g\n",avgstatloci[x].theta_fulin_nut);
                                    }
                                }
                            }
							fclose(file_output2);
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
						fprintf(file_output,"\n     3.1.1.0.0. Display all the results in a table (VERY LONG OUTPUT!).\n");
					}

					printf("\n\n This information will be sent into a new file. Please introduce the name of the file: ");
					scanf("%s",fileout2);

					if((file_output2 = fopen (fileout2,"w")) == 0) {
						printf("\n  It is not possible to open the file %s",fileout2);
						if(file_output) fprintf(file_output,"\n  It is not possible to open the file %s",fileout2);
						return;
					}
					fprintf(file_output," Data results are sent to file %s.\n",fileout2);

                    if(*outgroup) {
                        /*
						printf("\n\n Display multilocus neutrality tests: \n\n Multiple hits not included in the analysis.\n");
                        printf(" The table contains the Total, averages and variances of\n");
                        printf(" Tajima's D test, Fu and Li's D and F tests, normalized Fay and Wu H test, Fay and Wu's H test, Fu's Fs test,\n");
                        printf(" Rozas' et a. ZA test, Wall's B and Q tests, Ramos-Onsins and Rozas R2 test, Zeng et al. E test and Ewens-Watterson test. \n");
                        if(data[0].n_loci > 1 && data[0].time_spec > 0.0) printf(" Also the result of the total Chi-square value and Time and Theta estimates for HKA test (values corrected by Jukes and Cantor).\n");
                        printf(" 'na' means the test is not calculated because the number of samples or the number of segregating sites\n");
                        printf("were insufficient.\n");
                        
                        printf("\nniter\tTotTajD\tAvgTajD\tVarTajD\t");
                        printf("TotFuLiD\tAvgFuLiD\tVarFuLiD\t");
                        printf("TotFuLiF\tAvgFuLiF\tVarFuLiF\t");
                        printf("TotFayWuHn\tAvgFayWuHn\tVarFayWuHn\t");
                        printf("TotFayWuH\tAvgFayWuH\tVarFayWuH\t");
                        printf("TotFuFs\tAvgFuFs\tVarFuFs\t");
                        printf("TotRozZA\tAvgRozZA\tVarRozZA\t");
                        printf("TotWallB\tAvgWallB\tVarWallB\t");
                        printf("TotWallQ\tAvgWallQ\tVarWallQ\t");
                        printf("TotR2\tAvgR2\tVarR2\t");
                        printf("TotZE\tAvgZE\tVarZE\t");
                        printf("TotEW\tAvgEW\tVarEW\t");
                        if(data[0].n_loci > 1 && data[0].time_spec > 0.0) {
                            printf("hka\tE(T)\t");
                            for(x=0;x<data[0].n_loci;x++) printf("E(theta%d)\t",x);
                        }
                        printf("\n");

                        for(it=0;it<(long int)data[0].n_iter;it++) {
                            printf("%ld\t",it);
                            if(matrixmlsim[it].nltajD > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].StajimaD,matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD);
                                if(matrixmlsim[it].nltajD > 2) {
                                    var = variances((double)matrixmlsim[it].S2tajimaD,(double)matrixmlsim[it].StajimaD,matrixmlsim[it].nltajD);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlflD > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SfuliD,matrixmlsim[it].SfuliD/(double)matrixmlsim[it].nlflD);
                                if(matrixmlsim[it].nlflD > 2) {
                                    var = variances((double)matrixmlsim[it].S2fuliD,(double)matrixmlsim[it].SfuliD,matrixmlsim[it].nlflD);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlflF > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SfuliF,matrixmlsim[it].SfuliF/(double)matrixmlsim[it].nlflF);
                                if(matrixmlsim[it].nlflF > 2) {
                                    var = variances((double)matrixmlsim[it].S2fuliF,(double)matrixmlsim[it].SfuliF,matrixmlsim[it].nlflF);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlH > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SfaywuH,matrixmlsim[it].SfaywuH/(double)matrixmlsim[it].nlH);
                                if(matrixmlsim[it].nlH > 2) {
                                    var = variances((double)matrixmlsim[it].S2faywuH,(double)matrixmlsim[it].SfaywuH,matrixmlsim[it].nlH);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlHo > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SfaywuHo,matrixmlsim[it].SfaywuHo/(double)matrixmlsim[it].nlHo);
                                if(matrixmlsim[it].nlHo > 2) {
                                    var = variances((double)matrixmlsim[it].S2faywuHo,(double)matrixmlsim[it].SfaywuHo,matrixmlsim[it].nlHo);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlFs > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SfuFs,matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs);
                                if(matrixmlsim[it].nlFs > 2) {
                                    var = variances((double)matrixmlsim[it].S2fuFs,(double)matrixmlsim[it].SfuFs,matrixmlsim[it].nlFs);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlZ > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SrZA,matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ);
                                if(matrixmlsim[it].nlZ > 2) {
                                    var = variances((double)matrixmlsim[it].S2rZA,(double)matrixmlsim[it].SrZA,matrixmlsim[it].nlZ);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlB > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SwB,matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB);
                                if(matrixmlsim[it].nlB > 2) {
                                    var = variances((double)matrixmlsim[it].S2wB,(double)matrixmlsim[it].SwB,matrixmlsim[it].nlB);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlQ > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SwQ,matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ);
                                if(matrixmlsim[it].nlQ > 2) {
                                    var = variances((double)matrixmlsim[it].S2wQ,(double)matrixmlsim[it].SwQ,matrixmlsim[it].nlQ);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        
                            
                            if(matrixmlsim[it].nlR2 > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SR2,matrixmlsim[it].SR2/(double)matrixmlsim[it].nlR2);
                                if(matrixmlsim[it].nlR2 > 2) {
                                    var = variances((double)matrixmlsim[it].S2R2,(double)matrixmlsim[it].SR2,matrixmlsim[it].nlR2);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlE > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SzengE,matrixmlsim[it].SzengE/(double)matrixmlsim[it].nlE);
                                if(matrixmlsim[it].nlE > 2) {
                                    var = variances((double)matrixmlsim[it].S2zengE,(double)matrixmlsim[it].SzengE,matrixmlsim[it].nlE);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

							if(matrixmlsim[it].n_loci > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].Sewtest,matrixmlsim[it].Sewtest/(double)matrixmlsim[it].n_loci);
                                if(matrixmlsim[it].n_loci > 2) {
                                    var = variances((double)matrixmlsim[it].S2ewtest,(double)matrixmlsim[it].Sewtest,matrixmlsim[it].n_loci);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

							if(data[0].n_loci > 1 && data[0].time_spec > 0.0) {
                                printf("%g\t",matrixmlsim[it].Shka);
                                printf("%g\t",matrixmlsim[it].hka_T);
                                for(x=0;x<data[0].n_loci;x++) {
                                    printf("%g\t",matrixsim[x][it].hka_theta);
                                }
                            }
                            printf("\n");
                        }
                        */
                        /* print the median/average for all statistic in each loci.*/
                        if(onlymulo)
                            printf("\n Value of the average obtained by simulations in each locus:\n\n");
                        else 
                            printf("\n Value of the median obtained by simulations in each locus:\n\n");

                        printf("\nlocus\t");
                        printf("TajD\t");
                        printf("FuLiD\t");
                        printf("FuLiF\t");
                        printf("FayWuHn\t");
                        printf("FayWuH\t");
                        printf("FuFs\t");
                        printf("RozZA\t");
                        printf("WallB\t");
                        printf("WallQ\t");
                        printf("R2\t");
                        printf("ZengE\t");
                        printf("EW\n");
                        
                        for(x=0;x<data[0].n_loci;x++) {
                            if(observed_data && dataobsequalsim) printf("%d:%s\t",x,matrix[x].gene);
                            else printf("%d\t",x);
                            printf("%g\t",avgstatloci[x].tajimaD);
                            printf("%g\t",avgstatloci[x].fuliD);
                            printf("%g\t",avgstatloci[x].fuliF);
                            printf("%g\t",avgstatloci[x].faywuH);
                            printf("%g\t",avgstatloci[x].faywuHo);
                            printf("%g\t",avgstatloci[x].fuFs);
                            printf("%g\t",avgstatloci[x].rZA);
                            printf("%g\t",avgstatloci[x].wB);
                            printf("%g\t",avgstatloci[x].wQ);
                            printf("%g\t",avgstatloci[x].R2);
                            printf("%g\t",avgstatloci[x].zengE);
                            printf("%g\n",avgstatloci[x].ewtest);
                        }

                        if(file_output2) {
                            /*print to file*/
                            fputs("\n\n Display multilocus neutrality tests: \n\n Multiple hits not included in the analysis.\n",file_output2);
                            fputs(" The table contains the Total, averages and variances of\n",file_output2);
                            fputs(" Tajima's D test, Fu and Li's D and F tests, normalized Fay and Wu H test, Fay and Wu's H test, Fu's Fs test,\n",file_output2);
                            fputs("  Rozas' et a. ZA test, Wall's B and Q tests, Ramos-Onsins and Rozas R2 test, Zeng et al. E test and Ewens-Watterson test.\n",file_output2);
                            if(data[0].n_loci > 1 && data[0].time_spec > 0.0) fputs(" Also the result of the total Chi-square value and Time and Theta estimates for HKA test (values corrected by Jukes and Cantor).\n",file_output2);
                            fputs(" 'na' means the test is not calculated because the number of samples or the number of segregating sites\n",file_output2);
                            fputs("were insufficient.\n",file_output2);
                            
                            fprintf(file_output2,"\nniter\tTotTajD\tAvgTajD\tVarTajD\t");
                            fprintf(file_output2,"TotFuLiD\tAvgFuLiD\tVarFuLiD\t");
                            fprintf(file_output2,"TotFuLiF\tAvgFuLiF\tVarFuLiF\t");
                            fprintf(file_output2,"TotFayWuHn\tAvgFayWuHn\tVarFayWuHn\t");
                            fprintf(file_output2,"TotFayWuH\tAvgFayWuH\tVarFayWuH\t");
                            fprintf(file_output2,"TotFuFs\tAvgFuFs\tVarFuFs\t");
                            fprintf(file_output2,"TotRozZA\tAvgRozZA\tVarRozZA\t");
                            fprintf(file_output2,"TotWallB\tAvgWallB\tVarWallB\t");
                            fprintf(file_output2,"TotWallQ\tAvgWallQ\tVarWallQ\t");
                            fprintf(file_output2,"TotR2\tAvgR2\tVarR2\t");
                            fprintf(file_output2,"TotZengE\tAvgZengE\tVarZengE\t");
                            fprintf(file_output2,"TotEW\tAvgEW\tVarEW\t");
                            if(data[0].n_loci > 1 && data[0].time_spec > 0.0) {
                                fputs("hka\tE(T)\t",file_output2);
                                /*for(x=0;x<data[0].n_loci;x++) fprintf(file_output2,"E(theta%d)\t",x);*/
                            }
                            fputs("\n",file_output2);
    
                            for(it=0;it<(long int)data[0].n_iter;it++) {
                                fprintf(file_output2,"%ld\t",it);
                                if(matrixmlsim[it].nltajD > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].StajimaD,matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD);
                                    if(matrixmlsim[it].nltajD > 2) {
                                        var = variances((double)matrixmlsim[it].S2tajimaD,(double)matrixmlsim[it].StajimaD,matrixmlsim[it].nltajD);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlflD > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SfuliD,matrixmlsim[it].SfuliD/(double)matrixmlsim[it].nlflD);
                                    if(matrixmlsim[it].nlflD > 2) {
                                        var = variances((double)matrixmlsim[it].S2fuliD,(double)matrixmlsim[it].SfuliD,matrixmlsim[it].nlflD);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlflF > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SfuliF,matrixmlsim[it].SfuliF/(double)matrixmlsim[it].nlflF);
                                    if(matrixmlsim[it].nlflF > 2) {
                                        var = variances((double)matrixmlsim[it].S2fuliF,(double)matrixmlsim[it].SfuliF,matrixmlsim[it].nlflF);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlH > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SfaywuH,matrixmlsim[it].SfaywuH/(double)matrixmlsim[it].nlH);
                                    if(matrixmlsim[it].nlH > 2) {
                                        var = variances((double)matrixmlsim[it].S2faywuH,(double)matrixmlsim[it].SfaywuH,matrixmlsim[it].nlH);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
								
                                if(matrixmlsim[it].nlHo > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SfaywuHo,matrixmlsim[it].SfaywuHo/(double)matrixmlsim[it].nlHo);
                                    if(matrixmlsim[it].nlHo > 2) {
                                        var = variances((double)matrixmlsim[it].S2faywuHo,(double)matrixmlsim[it].SfaywuHo,matrixmlsim[it].nlHo);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlFs > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SfuFs,matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs);
                                    if(matrixmlsim[it].nlFs > 2) {
                                        var = variances((double)matrixmlsim[it].S2fuFs,(double)matrixmlsim[it].SfuFs,matrixmlsim[it].nlFs);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlZ > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SrZA,matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ);
                                    if(matrixmlsim[it].nlZ > 2) {
                                        var = variances((double)matrixmlsim[it].S2rZA,(double)matrixmlsim[it].SrZA,matrixmlsim[it].nlZ);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlB > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SwB,matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB);
                                    if(matrixmlsim[it].nlB > 2) {
                                        var = variances((double)matrixmlsim[it].S2wB,(double)matrixmlsim[it].SwB,matrixmlsim[it].nlB);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlQ > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SwQ,matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ);
                                    if(matrixmlsim[it].nlQ > 2) {
                                        var = variances((double)matrixmlsim[it].S2wQ,(double)matrixmlsim[it].SwQ,matrixmlsim[it].nlQ);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
                                
                                if(matrixmlsim[it].nlR2 > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SR2,matrixmlsim[it].SR2/(double)matrixmlsim[it].nlR2);
                                    if(matrixmlsim[it].nlR2 > 2) {
                                        var = variances((double)matrixmlsim[it].S2R2,(double)matrixmlsim[it].SR2,matrixmlsim[it].nlR2);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        

                                if(matrixmlsim[it].nlE > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SzengE,matrixmlsim[it].SzengE/(double)matrixmlsim[it].nlE);
                                    if(matrixmlsim[it].nlE > 2) {
                                        var = variances((double)matrixmlsim[it].S2zengE,(double)matrixmlsim[it].SzengE,matrixmlsim[it].nlE);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        

                                if(data[0].n_loci > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].Sewtest,matrixmlsim[it].Sewtest/(double)data[0].n_loci);
                                    if(data[0].n_loci > 2) {
                                        var = variances((double)matrixmlsim[it].S2ewtest,(double)matrixmlsim[it].Sewtest,data[0].n_loci);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        

                                if(data[0].n_loci > 1 && data[0].time_spec > 0.0) {
                                    fprintf(file_output2,"%g\t",matrixmlsim[it].Shka);
                                    fprintf(file_output2,"%g\t",matrixmlsim[it].hka_T);
                                    /*
									 for(x=0;x<data[0].n_loci;x++) {
                                        fprintf(file_output2,"%g\t",matrixsim[x][it].hka_theta);
                                    }
									*/
                                }
                                fputs("\n",file_output2);
                            }
                            
                            /* print the median/average for all statistic in each loci.*/
                            if(onlymulo)
                                fputs("\n Value of the average obtained by simulations in each locus:\n\n",file_output2);
                            else 
                                fputs("\n Value of the median obtained by simulations in each locus:\n\n",file_output2);
    
                            fprintf(file_output2,"\nlocus\t");
                            fputs("TajD\t",file_output2);
                            fputs("FuLiD\t",file_output2);
                            fputs("FuLiF\t",file_output2);
                            fputs("FayWuHn\t",file_output2);
                            fputs("FayWuH\t",file_output2);
                            fputs("FuFs\t",file_output2);
                            fputs("RozZA\t",file_output2);
                            fputs("WallB\t",file_output2);
                            fputs("WallQ\t",file_output2);
                            fputs("R2\t",file_output2);
                            fputs("ZengE\t",file_output2);
                            fputs("EW\n",file_output2);
                            
                            for(x=0;x<data[0].n_loci;x++) {
                                if(observed_data && dataobsequalsim) fprintf(file_output2,"%d:%s\t",x,matrix[x].gene);
                                else fprintf(file_output2,"%d\t",x);
                                fprintf(file_output2,"%g\t",avgstatloci[x].tajimaD);
                                fprintf(file_output2,"%g\t",avgstatloci[x].fuliD);
                                fprintf(file_output2,"%g\t",avgstatloci[x].fuliF);
                                fprintf(file_output2,"%g\t",avgstatloci[x].faywuH);
                                fprintf(file_output2,"%g\t",avgstatloci[x].faywuHo);
                                fprintf(file_output2,"%g\t",avgstatloci[x].fuFs);
                                fprintf(file_output2,"%g\t",avgstatloci[x].rZA);
                                fprintf(file_output2,"%g\t",avgstatloci[x].wB);
                                fprintf(file_output2,"%g\t",avgstatloci[x].wQ);
                                fprintf(file_output2,"%g\t",avgstatloci[x].R2);
                               fprintf(file_output2,"%g\t",avgstatloci[x].zengE);
                               fprintf(file_output2,"%g\n",avgstatloci[x].ewtest);
                            }
                        }
                    }
                    else {
                        /*no outgroup*/
                        printf("\n\n Display multilocus neutrality tests: \n\n Multiple hits not included in the analysis.\n");
                        /*
						printf(" The table contains the Total, averages and variances of\n");
                        printf(" Tajima's D test, Fu and Li's D and F tests, Fu's Fs test,\n");
                        printf(" Rozas' et a. ZA test, Wall's B and Q tests, Ramos-Onsins and Rozas R2 test and Ewens-Watterson test.\n");
                        printf(" 'na' means the test is not calculated because the number of samples or the number of segregating sites\n");
                        printf(" were insufficient.\n");
                        
                        printf("\nniter\tTotTajD\tAvgTajD\tVarTajD\t");
                        printf("TotFuLiD*\tAvgFuLiD*\tVarFuLiD*\t");
                        printf("TotFuLiF*\tAvgFuLiF*\tVarFuLiF*\t");
                        printf("TotFuFs\tAvgFuFs\tVarFuFs\t");
                        printf("TotRozZA\tAvgRozZA\tVarRozZA\t");
                        printf("TotWallB\tAvgWallB\tVarWallB\t");
                        printf("TotWallQ\tAvgWallQ\tVarWallQ\t");
                        printf("TotR2\tAvgR2\tVarR2\t");
                        printf("TotEW\tAvgEW\tVarEW\t");
                        printf("\n");

                        for(it=0;it<(long int)data[0].n_iter;it++) {
                            printf("%ld\t",it);
                            if(matrixmlsim[it].nltajD > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].StajimaD,matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD);
                                if(matrixmlsim[it].nltajD > 2) {
                                    var = variances((double)matrixmlsim[it].S2tajimaD,(double)matrixmlsim[it].StajimaD,matrixmlsim[it].nltajD);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlflDn > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SfuliDn,matrixmlsim[it].SfuliDn/(double)matrixmlsim[it].nlflDn);
                                if(matrixmlsim[it].nlflDn > 2) {
                                    var = variances((double)matrixmlsim[it].S2fuliDn,(double)matrixmlsim[it].SfuliDn,matrixmlsim[it].nlflDn);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlflFn > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SfuliFn,matrixmlsim[it].SfuliFn/(double)matrixmlsim[it].nlflFn);
                                if(matrixmlsim[it].nlflFn > 2) {
                                    var = variances((double)matrixmlsim[it].S2fuliFn,(double)matrixmlsim[it].SfuliFn,matrixmlsim[it].nlflFn);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlFs > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SfuFs,matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs);
                                if(matrixmlsim[it].nlFs > 2) {
                                    var = variances((double)matrixmlsim[it].S2fuFs,(double)matrixmlsim[it].SfuFs,matrixmlsim[it].nlFs);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlZ > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SrZA,matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ);
                                if(matrixmlsim[it].nlZ > 2) {
                                    var = variances((double)matrixmlsim[it].S2rZA,(double)matrixmlsim[it].SrZA,matrixmlsim[it].nlZ);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlB > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SwB,matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB);
                                if(matrixmlsim[it].nlB > 2) {
                                    var = variances((double)matrixmlsim[it].S2wB,(double)matrixmlsim[it].SwB,matrixmlsim[it].nlB);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlQ > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SwQ,matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ);
                                if(matrixmlsim[it].nlQ > 2) {
                                    var = variances((double)matrixmlsim[it].S2wQ,(double)matrixmlsim[it].SwQ,matrixmlsim[it].nlQ);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            if(matrixmlsim[it].nlR2 > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].SR2,matrixmlsim[it].SR2/(double)matrixmlsim[it].nlR2);
                                if(matrixmlsim[it].nlR2 > 2) {
                                    var = variances((double)matrixmlsim[it].S2R2,(double)matrixmlsim[it].SR2,matrixmlsim[it].nlR2);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

							if(matrixmlsim[it].n_loci > 0) {
                                printf("%g\t%g\t",matrixmlsim[it].Sewtest,matrixmlsim[it].Sewtest/(double)matrixmlsim[it].n_loci);
                                if(matrixmlsim[it].n_loci > 2) {
                                    var = variances((double)matrixmlsim[it].S2ewtest,(double)matrixmlsim[it].Sewtest,matrixmlsim[it].n_loci);
                                    printf("%g\t",var);
                                }
                                else printf("na\t");
                            }
                            else printf("na\tna\tna\t");        

                            printf("\n");
                        }
                        */
                        /* print the median/average for all statistic in each loci.*/
                        if(onlymulo)
                            printf("\n Value of the average obtained by simulations in each locus:\n\n");
                        else 
                            printf("\n Value of the median obtained by simulations in each locus:\n\n");

                        printf("locus\t");
                        printf("TajD\t");
                        printf("FuLiD*\t");
                        printf("FuLiF*\t");
                        printf("FuFs\t");
                        printf("RozZA\t");
                        printf("WallB\t");
                        printf("WallQ\t");
                        printf("R2\t");
                        printf("EW\n");
                        
                        for(x=0;x<data[0].n_loci;x++) {
                            if(observed_data && dataobsequalsim) printf("%d:%s\t",x,matrix[x].gene);
                            else printf("%d\t",x);
                            printf("%g\t",avgstatloci[x].tajimaD);
                            printf("%g\t",avgstatloci[x].fuliDn);
                            printf("%g\t",avgstatloci[x].fuliFn);
                            printf("%g\t",avgstatloci[x].fuFs);
                            printf("%g\t",avgstatloci[x].rZA);
                            printf("%g\t",avgstatloci[x].wB);
                            printf("%g\t",avgstatloci[x].wQ);
                            printf("%g\t",avgstatloci[x].R2);
                            printf("%g\n",avgstatloci[x].ewtest);
                        }
                        
                        if(file_output2) {
                            /*print to file*/
                            fputs("\n\n Display multilocus neutrality tests: \n\n Multiple hits not included in the analysis.\n",file_output2);
                            fputs(" The table contains the Total, averages and variances of\n",file_output2);
                            fputs(" Tajima's D test, Fu and Li's D and F tests, Fu's Fs test,\n",file_output2);
                            fputs(" Rozas' et a. ZA test, Wall's B and Q tests, Ramos-Onsins and Rozas R2 test and Ewens-Watterson test.\n",file_output2);
                            fputs(" 'na' means the test is not calculated because the number of samples or the number of segregating sites\n",file_output2);
                            fputs(" were insufficient.\n",file_output2);
                            
                            fprintf(file_output2,"\nniter\tTotTajD\tAvgTajD\tVarTajD\t");
                            fprintf(file_output2,"TotFuLiD*\tAvgFuLiD*\tVarFuLiD*\t");
                            fprintf(file_output2,"TotFuLiF*\tAvgFuLiF*\tVarFuLiF*\t");
                            fprintf(file_output2,"TotFuFs\tAvgFuFs\tVarFuFs\t");
                            fprintf(file_output2,"TotRozZA\tAvgRozZA\tVarRozZA\t");
                            fprintf(file_output2,"TotWallB\tAvgWallB\tVarWallB\t");
                            fprintf(file_output2,"TotWallQ\tAvgWallQ\tVarWallQ\t");
                            fprintf(file_output2,"TotR2\tAvgR2\tVarR2\t");
                            fprintf(file_output2,"TotEW\tAvgEW\tVarEW\t");
                            fputs("\n",file_output2);
    
                            for(it=0;it<(long int)data[0].n_iter;it++) {
                                fprintf(file_output2,"%ld\t",it);
                                if(matrixmlsim[it].nltajD > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].StajimaD,matrixmlsim[it].StajimaD/(double)matrixmlsim[it].nltajD);
                                    if(matrixmlsim[it].nltajD > 2) {
                                        var = variances((double)matrixmlsim[it].S2tajimaD,(double)matrixmlsim[it].StajimaD,matrixmlsim[it].nltajD);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlflDn > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SfuliDn,matrixmlsim[it].SfuliDn/(double)matrixmlsim[it].nlflDn);
                                    if(matrixmlsim[it].nlflDn > 2) {
                                        var = variances((double)matrixmlsim[it].S2fuliDn,(double)matrixmlsim[it].SfuliDn,matrixmlsim[it].nlflDn);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlflFn > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SfuliFn,matrixmlsim[it].SfuliFn/(double)matrixmlsim[it].nlflFn);
                                    if(matrixmlsim[it].nlflFn > 2) {
                                        var = variances((double)matrixmlsim[it].S2fuliFn,(double)matrixmlsim[it].SfuliFn,matrixmlsim[it].nlflFn);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlFs > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SfuFs,matrixmlsim[it].SfuFs/(double)matrixmlsim[it].nlFs);
                                    if(matrixmlsim[it].nlFs > 2) {
                                        var = variances((double)matrixmlsim[it].S2fuFs,(double)matrixmlsim[it].SfuFs,matrixmlsim[it].nlFs);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlZ > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SrZA,matrixmlsim[it].SrZA/(double)matrixmlsim[it].nlZ);
                                    if(matrixmlsim[it].nlZ > 2) {
                                        var = variances((double)matrixmlsim[it].S2rZA,(double)matrixmlsim[it].SrZA,matrixmlsim[it].nlZ);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlB > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SwB,matrixmlsim[it].SwB/(double)matrixmlsim[it].nlB);
                                    if(matrixmlsim[it].nlB > 2) {
                                        var = variances((double)matrixmlsim[it].S2wB,(double)matrixmlsim[it].SwB,matrixmlsim[it].nlB);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlQ > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SwQ,matrixmlsim[it].SwQ/(double)matrixmlsim[it].nlQ);
                                    if(matrixmlsim[it].nlQ > 2) {
                                        var = variances((double)matrixmlsim[it].S2wQ,(double)matrixmlsim[it].SwQ,matrixmlsim[it].nlQ);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        
    
                                if(matrixmlsim[it].nlR2 > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].SR2,matrixmlsim[it].SR2/(double)matrixmlsim[it].nlR2);
                                    if(matrixmlsim[it].nlR2 > 2) {
                                        var = variances((double)matrixmlsim[it].S2R2,(double)matrixmlsim[it].SR2,matrixmlsim[it].nlR2);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        

                                if(data[0].n_loci > 0) {
                                    fprintf(file_output2,"%g\t%g\t",matrixmlsim[it].Sewtest,matrixmlsim[it].Sewtest/(double)data[0].n_loci);
                                    if(data[0].n_loci > 2) {
                                        var = variances((double)matrixmlsim[it].S2ewtest,(double)matrixmlsim[it].Sewtest,data[0].n_loci);
                                        fprintf(file_output2,"%g\t",var);
                                    }
                                    else fputs("na\t",file_output2);
                                }
                                else fputs("na\tna\tna\t",file_output2);        

                                fputs("\n",file_output2);
                            }
                            
                            /* print the median/average for all statistic in each loci.*/
                            if(onlymulo)
                                fputs("\n Value of the average obtained by simulations in each locus:\n\n",file_output2);
                            else 
                                fputs("\n Value of the median obtained by simulations in each locus:\n\n",file_output2);
    
                            fprintf(file_output2,"\nlocus\t");
                            fputs("TajD\t",file_output2);
                            fputs("FuLiD*\t",file_output2);
                            fputs("FuLiF*\t",file_output2);
                            fputs("FuFs\t",file_output2);
                            fputs("RozZA\t",file_output2);
                            fputs("WallB\t",file_output2);
                            fputs("WallQ\t",file_output2);
                            fputs("R2\t",file_output2);
                            fputs("EW\n",file_output2);
                            
                            for(x=0;x<data[0].n_loci;x++) {
                                if(observed_data && dataobsequalsim) fprintf(file_output2,"%d:%s\t",x,matrix[x].gene);
                                else fprintf(file_output2,"%d\t",x);
                                fprintf(file_output2,"%g\t",avgstatloci[x].tajimaD);
                                fprintf(file_output2,"%g\t",avgstatloci[x].fuliDn);
                                fprintf(file_output2,"%g\t",avgstatloci[x].fuliFn);
                                fprintf(file_output2,"%g\t",avgstatloci[x].fuFs);
                                fprintf(file_output2,"%g\t",avgstatloci[x].rZA);
                                fprintf(file_output2,"%g\t",avgstatloci[x].wB);
                                fprintf(file_output2,"%g\t",avgstatloci[x].wQ);
                                fprintf(file_output2,"%g\t",avgstatloci[x].R2);
                                fprintf(file_output2,"%g\n",avgstatloci[x].ewtest);
                            }
                        }
                    }
					fclose(file_output2);
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
	if(file_output) fflush(file_output);
    return;
}

