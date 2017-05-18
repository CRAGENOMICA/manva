/*
 *  displayobsdata.c
 *  MuLoNeTests
 *
 *  Created by sebas on Wed Feb 26 2003.
 *
 */

#include "MuLoNeTests.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void open_displstdetgrid(struct statistics *matrix/*,int *observed_data*/,int *outgroup,int *n_loci,char *subset_positions,FILE *file_output)
{
    char k[1];
    int x,y;
    long int ss;
    double jc;

    #if COMMAND_LINE
    if(file_output) fflush(file_output);
    
    printf("\n\nMENU 1. Observed analysis menu:");
    printf("\n     1.0. Display statistics menu:");
    printf("\n     1.0.1. Display detailed statistics menu:");
    printf("\n     1.0.1.0. Display table with detailed statistics for each locus:\n\n");
    
    printf("CHOOSE:\n");
    printf(" 0 - Display general statistics.\n");
    printf(" 1 - Display other linkage related statistics.\n");
    printf(" 2 - Display estimates of total locus variability.\n");
    printf(" 3 - Display estimates of nucleotide locus variability.\n");
    printf(" 4 - Back to previous menu.\n");

    if(file_output) {
        
        fprintf(file_output,"\n\n     MENU:\n     1. Observed analysis menu:");
        fprintf(file_output,"\n     1.0. Display statistics menu:");
        fprintf(file_output,"\n     1.0.1. Display detailed statistics menu:");
        fprintf(file_output,"\n     1.0.1.0. Display table with detailed statistics for each locus:\n\n");
        
        fprintf(file_output," 0 - Display general statistics.\n");
        fprintf(file_output," 1 - Display other linkage related statistics.\n");
        fprintf(file_output," 2 - Display estimates of total locus variability.\n");
        fprintf(file_output," 3 - Display estimates of nucleotide locus variability.\n");
        fprintf(file_output," 4 - Back to previous menu.\n");
    }
    #endif
    
    do *k = getchar();
    while(*k<'0' || *k>'4');
    
    switch(*k) {
        case '0':
            if(*outgroup) {
                ss = 0;
                for(y=0;y<*n_loci;y++) ss += matrix[y].shared;
                printf("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n");
                printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n");
				if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
                printf(" The table contains the name of the locus, correction factor for chromosome population size, the transition/transversion ratio(s/v) and the number of samples, outgroups, sites, multiple hits, ");
                printf("segregating biallelic sites\n");
                if(ss > 0) printf("not including shared polymorphisms (nsh) or including them (sh), number of shared polymorphisms, \n");
                printf(" number of fixed mutations, averaged divergence, the number of haplotypes, haplotypes/sample size, haplotype diversity\n");
                printf(" and the recombination parameter per loci (Hudson 1987).\n");
                if(ss > 0) printf(" Shared polymorphisms are not included in counting the number of haplotypes.\n");
                printf("\n\nname_loci\tfactor_chrn\ts/v\t#samples\t#outg\t#sites\t#mhits\t");
                if(ss > 0) printf("S(nsh)\tS(sh)\t#shared\t");
                else printf("Seg_sites\t");
                printf("#fixed\tDiverg\t#haplot\thapl/sam\thapldiv\tC\tRm\n");
                for(x=0;x<*n_loci;x++) {
                    printf("\n%d:%-20s\t",x,matrix[x].gene);
                    printf("%g\t",matrix[x].factor_chrn);
                    if(matrix[x].transversions) printf("%.2f\t",(double)matrix[x].transitions/matrix[x].transversions);
					else printf("na\t");
					printf("%d\t",matrix[x].nsamples);
                    printf("%d\t",matrix[x].noutgroups);
                    printf("%.2f\t",matrix[x].nsites);
                    printf("%d\t",matrix[x].nmhits);
                    printf("%d\t",matrix[x].biallsites);
                    if(ss > 0) {
                        printf("%d\t",matrix[x].biallsitesn);
                        printf("%d\t",matrix[x].shared);
                    }
                    printf("%d\t",matrix[x].fixed);
                    printf("%g\t",matrix[x].ndivergence);
                    printf("%d\t",matrix[x].nhapl);
                    printf("%g\t",matrix[x].nhaplsam);
                    printf("%g\t",matrix[x].hapldiv);
                    printf("%g\t",matrix[x].Rvpi);
                    printf("%d\t",matrix[x].Rm);
                }
                printf("\n\n");
                
                if(file_output) {
                    ss = 0;
                    for(y=0;y<*n_loci;y++) ss += matrix[y].shared;
                    fputs("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					fputs(" The table contains the name of the locus, correction factor for chromosome population size, the transition/transversion ratio(s/v) and the number of samples, outgroups, sites, multiple hits, ",file_output);
                    fputs("segregating biallelic sites\n",file_output);
                    if(ss > 0) fputs("not including shared polymorphisms (nsh) or including them (sh), number of shared polymorphisms,\n",file_output);
                    fputs(" number of fixed mutations, averaged divergence, the number of haplotypes, haplotypes/sample size, haplotype diversity\n",file_output);
					fputs(" and the recombination parameter per loci (Hudson 1987).\n",file_output);
                    if(ss > 0) fputs(" Shared polymorphisms are not included in counting the number of haplotypes.\n",file_output);
                    fputs("\n\nname_loci\tfactor_chrn\ts/v\t#samples\t#outg\t#sites\t#mhits\t",file_output);
                    if(ss > 0) fputs("S(nsh)\tS(sh)\t#shared\t",file_output);
                    else fputs("Seg_sites\t",file_output);
                    fputs("#fixed\tDiverg\t#haplot\thapl/sam\thapldiv\tC\tRm\n",file_output);
                    for(x=0;x<*n_loci;x++) {
                        fprintf(file_output,"\n%d:%-20s\t",x,matrix[x].gene);
                        fprintf(file_output,"%g\t",matrix[x].factor_chrn);
						if(matrix[x].transversions) fprintf(file_output,"%.2f\t",(double)matrix[x].transitions/matrix[x].transversions);
						else fprintf(file_output,"na\t");
                        fprintf(file_output,"%d\t",matrix[x].nsamples);
                        fprintf(file_output,"%d\t",matrix[x].noutgroups);
                        fprintf(file_output,"%.2f\t",matrix[x].nsites);
                        fprintf(file_output,"%d\t",matrix[x].nmhits);
                        fprintf(file_output,"%d\t",matrix[x].biallsites);
                        if(ss > 0) {
                            fprintf(file_output,"%d\t",matrix[x].biallsitesn);
                            fprintf(file_output,"%d\t",matrix[x].shared);
                        }
                        fprintf(file_output,"%d\t",matrix[x].fixed);
                        fprintf(file_output,"%g\t",matrix[x].ndivergence);
                        fprintf(file_output,"%d\t",matrix[x].nhapl);
						fprintf(file_output,"%g\t",matrix[x].nhaplsam);
						fprintf(file_output,"%g\t",matrix[x].hapldiv);
                        fprintf(file_output,"%g\t",matrix[x].Rvpi);
                        fprintf(file_output,"%d\t",matrix[x].Rm);
                    }
                    fputs("\n\n",file_output);
                }
            }
            else {
                printf("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				printf(" The table contains the name of the locus, correction factor for chromosome population size, the transition/transversion ratio(s/v) and the number of samples, outgroups, sites, multiple hits, \n ");
                printf(" segregating biallelic sites, the number of haplotypes, haplotypes/sample size, haplotype diversity\n");
                printf(" and the recombination parameter per loci (Hudson 1987).\n");
                printf("\n\nname_loci\tfactor_chrn\ts/v\t#samples\t#outg\t#sites\t#mhits\t");
                printf("Seg_sites\t");
                printf("#haplot\thapl/sam\thapldiv\t");
                printf("C\tRm\n");
                
                for(x=0;x<*n_loci;x++) {
                    printf("\n%d:%-20s\t",x,matrix[x].gene);
                    printf("%g\t",matrix[x].factor_chrn);
                    if(matrix[x].transversions) printf("%.2f\t",(double)matrix[x].transitions/matrix[x].transversions);
					else printf("na\t");
                    printf("%d\t",matrix[x].nsamples);
                    printf("%d\t",matrix[x].noutgroups);
                    printf("%.2f\t",matrix[x].nsites);
                    printf("%d\t",matrix[x].nmhits);
                    printf("%d\t",matrix[x].biallsites);
                    printf("%d\t",matrix[x].nhapl);
                    printf("%g\t",matrix[x].nhaplsam);
                    printf("%g\t",matrix[x].hapldiv);
                    printf("%g\t",matrix[x].Rvpi);
                    printf("%d\t",matrix[x].Rm);
                }
                printf("\n\n");
                
                if(file_output) {
                    fputs("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					fputs(" The table contains the name of the locus, correction factor for chromosome population size, the transition/transversion ratio(s/v) and the number of samples, outgroups, sites, multiple hits, \n ",file_output);
                    fputs(" segregating biallelic sites, the number of haplotypes, haplotypes/sample size, haplotype diversity\n",file_output);
					fputs(" and the recombination parameter per loci (Hudson 1987) corrected for chromosome population size (1/f).\n",file_output);
                    fputs("\n\nname_loci\tfactor_chrn\ts/v\t#samples\t#outg\t#sites\t#mhits\t",file_output);
                    fputs("Seg_sites\t",file_output);
                    fputs("#haplot\thapl/sam\thapldiv\t",file_output);
                    fputs("C\tRm\n",file_output);
                    
                    for(x=0;x<*n_loci;x++) {
                        fprintf(file_output,"\n%d:%-20s\t",x,matrix[x].gene);
                        fprintf(file_output,"%g\t",matrix[x].factor_chrn);
						if(matrix[x].transversions) fprintf(file_output,"%.2f\t",(double)matrix[x].transitions/matrix[x].transversions);
						else fprintf(file_output,"na\t");
                        fprintf(file_output,"%d\t",matrix[x].nsamples);
                        fprintf(file_output,"%d\t",matrix[x].noutgroups);
                        fprintf(file_output,"%.2f\t",matrix[x].nsites);
                        fprintf(file_output,"%d\t",matrix[x].nmhits);
                        fprintf(file_output,"%d\t",matrix[x].biallsites);
                        fprintf(file_output,"%d\t",matrix[x].nhapl);
						fprintf(file_output,"%g\t",matrix[x].nhaplsam);
						fprintf(file_output,"%g\t",matrix[x].hapldiv);
                        fprintf(file_output,"%g\t",matrix[x].Rvpi);
                        fprintf(file_output,"%d\t",matrix[x].Rm);
                    }
                    fputs("\n\n",file_output);
                }
            }    
            break;
        case '1':
            if(*outgroup) {
                ss = 0;
                for(y=0;y<*n_loci;y++) ss += matrix[y].shared;
                printf("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n");
                printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				printf(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples, outgroups, sites, ");
                printf("segregating biallelic sites\n");
                if(ss > 0) printf("not including shared polymorphisms (nsh) or including them (sh), \n");
                printf("Wall's not divided statistics B and Q,\n");
                printf(" and Rozas' et al. not divided ZA\n");
                if(ss > 0) printf(" Shared polymorphisms are not included in counting B, Q and ZA.\n");
                printf("\n\nname_loci\tfactor_chrn\t#samples\t#outg\t#sites\t");
                if(ss > 0) printf("S(nsh)\tS(sh)\t");
                else printf("Seg_sites\t");
                printf("b\tq\tza\n");
                for(x=0;x<*n_loci;x++) {
                    printf("\n%d:%-20s\t",x,matrix[x].gene);
                    printf("%g\t",matrix[x].factor_chrn);
                    printf("%d\t",matrix[x].nsamples);
                    printf("%d\t",matrix[x].noutgroups);
                    printf("%.2f\t",matrix[x].nsites);
                    printf("%d\t",matrix[x].biallsites);
                    if(ss > 0) {
                        printf("%d\t",matrix[x].biallsitesn);
                    }
                    printf("%d\t",matrix[x].b);
                    printf("%d\t",matrix[x].q);
                    printf("%g\t",matrix[x].za);
                }
                printf("\n\n");
                
                if(file_output) {
                    ss = 0;
                    for(y=0;y<*n_loci;y++) ss += matrix[y].shared;
                    fputs("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					fputs(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples, outgroups, sites, ",file_output);
                    fputs("segregating biallelic sites\n",file_output);
                    if(ss > 0) fputs("not including shared polymorphisms (nsh) or including them (sh), \n",file_output);
                    fputs(" Wall's not divided statistics B and Q,",file_output);
                    fputs(" and Rozas' et al. not divided ZA\n",file_output);
                    fputs("\n\nname_loci\tfactor_chrn\t#samples\t#outg\t#sites\t",file_output);
                    if(ss > 0) fputs(" Shared polymorphisms are not included in counting B, Q and ZA.\n",file_output);
                    if(ss > 0) fputs("S(nsh)\tS(sh)\t",file_output);
                    else fputs("Seg_sites\t",file_output);
                    fputs("b\tq\tza\n",file_output);
                    for(x=0;x<*n_loci;x++) {
                        fprintf(file_output,"\n%d:%-20s\t",x,matrix[x].gene);
                        fprintf(file_output,"%g\t",matrix[x].factor_chrn);
                        fprintf(file_output,"%d\t",matrix[x].nsamples);
                        fprintf(file_output,"%d\t",matrix[x].noutgroups);
                        fprintf(file_output,"%.2f\t",matrix[x].nsites);
                        fprintf(file_output,"%d\t",matrix[x].biallsites);
                        if(ss > 0) {
                            fprintf(file_output,"%d\t",matrix[x].biallsitesn);
                        }
                        fprintf(file_output,"%d\t",matrix[x].b);
                        fprintf(file_output,"%d\t",matrix[x].q);
                        fprintf(file_output,"%g\t",matrix[x].za);
                    }
                    fputs("\n\n",file_output);
                }
            }
            else {
                printf("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				printf(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples, outgroups, sites,\n ");
                printf(" segregating biallelic sites, Wall's not divided statistics B and Q,\n");
                printf(" and Rozas' et al. not divided ZA\n");
                printf("\n\nname_loci\tfactor_chrn\t#samples\t#outg\t#sites\t");
                printf("Seg_sites\t");
                printf("b\tq\tza\n");
                
                for(x=0;x<*n_loci;x++) {
                    printf("\n%d:%-20s\t",x,matrix[x].gene);
                    printf("%g\t",matrix[x].factor_chrn);
                    printf("%d\t",matrix[x].nsamples);
                    printf("%d\t",matrix[x].noutgroups);
                    printf("%.2f\t",matrix[x].nsites);
                    printf("%d\t",matrix[x].biallsites);
                    printf("%d\t",matrix[x].b);
                    printf("%d\t",matrix[x].q);
                    printf("%g\t",matrix[x].za);
                }
                printf("\n\n");
                
                if(file_output) {
                    fputs("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					fputs(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples, outgroups, sites,\n ",file_output);
                    fputs(" segregating biallelic sites, Wall's not divided statistics B and Q,\n",file_output);
                    fputs(" and Rozas' et al. not divided ZA\n",file_output);
                    fputs("\n\nname_loci\tfactor_chrn\t#samples\t#outg\t#sites\t",file_output);
                    fputs("Seg_sites\t",file_output);
                    fputs("b\tq\tza\n",file_output);
                    
                    for(x=0;x<*n_loci;x++) {
                        fprintf(file_output,"\n%d:%-20s\t",x,matrix[x].gene);
                        fprintf(file_output,"%g\t",matrix[x].factor_chrn);
                        fprintf(file_output,"%d\t",matrix[x].nsamples);
                        fprintf(file_output,"%d\t",matrix[x].noutgroups);
                        fprintf(file_output,"%.2f\t",matrix[x].nsites);
                        fprintf(file_output,"%d\t",matrix[x].biallsites);
                        fprintf(file_output,"%d\t",matrix[x].b);
                        fprintf(file_output,"%d\t",matrix[x].q);
                        fprintf(file_output,"%g\t",matrix[x].za);
                    }
                    fputs("\n\n",file_output);
                }
            }    
            break;
        case '2':
            if(*outgroup) {
                ss = 0;
                for(y=0;y<*n_loci;y++) ss += matrix[y].shared;
                printf("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n");
                printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				printf(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples, outgroups, sites, ");
                printf("segregating biallelic sites\n");
                if(ss > 0) printf("not including shared polymorphisms (nsh) or including them (sh), \n");
                printf(" averaged divergence (not corrected and with Jukes and Cantor correction)\n");
                printf(" and locus estimations of variability from Watterson, Tajima, Fu and Li,\n");
                printf(" Fay and Wu and Zeng et al.\n");
				printf(" NOTE that all theta estimates are corrected for chromosome population size.\n");
                printf("\n\nname_loci\tfactor_chrn\t#samples\t#outg\t#sites\t");
                if(ss > 0) printf("S(nsh)\tS(sh)\t");
                else printf("Seg_sites\tdiv\tdivJC\t");
                if(ss > 0) printf("theta_wat(nsh)\ttheta_wat(sh)\ttheta_taj(nsh)\ttheta_taj(sh)\t");
                else printf("theta_wat\ttheta_taj\t");
                printf("theta_fuli\ttheta_faywu\ttheta_zeng\n");

                for(x=0;x<*n_loci;x++) {
                    printf("\n%d:%-20s\t",x,matrix[x].gene);
                    printf("%g\t",matrix[x].factor_chrn);
                    printf("%d\t",matrix[x].nsamples);
                    printf("%d\t",matrix[x].noutgroups);
                    printf("%.2f\t",matrix[x].nsites);
                    printf("%d\t",matrix[x].biallsites);
                    if(ss > 0) {
                        printf("%d\t",matrix[x].biallsitesn);
                    }
                    printf("%g\t",matrix[x].ndivergence);
                    if((jc = matrix[x].ndivergence/(double)matrix[x].nsites) < (double)0.75) {
                        jc = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * jc) * (double)matrix[x].nsites;
                        printf("%g\t",jc);
                    }
                    else printf("na\t");
                    printf("%g\t",matrix[x].theta_wat*((double)1/matrix[x].factor_chrn));
                    if(ss > 0) printf("%g\t",matrix[x].theta_watn*((double)1/matrix[x].factor_chrn));
                    printf("%g\t",matrix[x].theta_taj*((double)1/matrix[x].factor_chrn));
                    if(ss > 0) printf("%g\t",matrix[x].theta_tajn*((double)1/matrix[x].factor_chrn));
                    printf("%g\t",matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn));
                    printf("%g\t",matrix[x].theta_fw*((double)1/matrix[x].factor_chrn));
                    printf("%g\t",matrix[x].theta_L*((double)1/matrix[x].factor_chrn));
                }
                printf("\n\n");
                
                if(file_output) {
                    ss = 0;
                    for(y=0;y<*n_loci;y++) ss += matrix[y].shared;
                    fputs("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n", file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					fputs(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples, outgroups, sites, ",file_output);
                    fputs("segregating biallelic sites\n",file_output);
                    if(ss > 0) fputs("not including shared polymorphisms (nsh) or including them (sh), \n",file_output);
                    fputs(" averaged divergence (not corrected and with Jukes and Cantor correction)\n",file_output);
                    fputs(" and locus estimations of variability from Watterson, Tajima, Fu and Li,\n",file_output);
                    fputs(" Fay and Wu and Zeng et al.\n",file_output);
					fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n",file_output);
                    fputs("\n\nname_loci\tfactor_chrn\t#samples\t#outg\t#sites\t",file_output);
                    if(ss > 0) fputs("S(nsh)\tS(sh)\t",file_output);
                    else fputs("Seg_sites\tdiv\tdivJC\t",file_output);
                    if(ss > 0) fputs("theta_wat(nsh)\ttheta_wat(sh)\ttheta_taj(nsh)\ttheta_taj(sh)\t",file_output);
                    else fputs("theta_wat\ttheta_taj\t",file_output);
                    fputs("theta_fuli\ttheta_faywu\ttheta_zeng\n",file_output);
                    for(x=0;x<*n_loci;x++) {
                        fprintf(file_output,"\n%d:%-20s\t",x,matrix[x].gene);
                        fprintf(file_output,"%g\t",matrix[x].factor_chrn);
                        fprintf(file_output,"%d\t",matrix[x].nsamples);
                        fprintf(file_output,"%d\t",matrix[x].noutgroups);
                        fprintf(file_output,"%.2f\t",matrix[x].nsites);
                        fprintf(file_output,"%d\t",matrix[x].biallsites);
                        if(ss > 0) {
                            fprintf(file_output,"%d\t",matrix[x].biallsitesn);
                        }
                        fprintf(file_output,"%g\t",matrix[x].ndivergence);
                        if((jc = matrix[x].ndivergence/(double)matrix[x].nsites) < (double)0.75) {
                            jc = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * jc) * (double)matrix[x].nsites;
                            fprintf(file_output,"%g\t",jc);
                        }
                        else fputs("na\t",file_output);
                        fprintf(file_output,"%g\t",matrix[x].theta_wat*((double)1/matrix[x].factor_chrn));
                        if(ss > 0) fprintf(file_output,"%g\t",matrix[x].theta_watn*((double)1/matrix[x].factor_chrn));
                        fprintf(file_output,"%g\t",matrix[x].theta_taj*((double)1/matrix[x].factor_chrn));
                        if(ss > 0) fprintf(file_output,"%g\t",matrix[x].theta_tajn*((double)1/matrix[x].factor_chrn));
                        fprintf(file_output,"%g\t",matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn));
                        fprintf(file_output,"%g\t",matrix[x].theta_fw*((double)1/matrix[x].factor_chrn));
                        fprintf(file_output,"%g\t",matrix[x].theta_L*((double)1/matrix[x].factor_chrn));
                    }
                    fputs("\n\n",file_output);
                }
            }
            else {
                printf("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				printf(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples, outgroups, sites, ");
                printf("segregating biallelic sites,\n");
                printf(" and locus estimations of nucleotide from Watterson, Tajima and Fu and Li.\n");
				printf(" NOTE that all theta estimates are corrected for chromosome population size.\n");
                printf("\n\nname_loci\tfactor_chrn\t#samples\t#outg\t#sites\t");
                printf("Seg_sites\t");
                printf("theta_wat\ttheta_taj\t");
                printf("theta_fuli\n");
                for(x=0;x<*n_loci;x++) {
                    printf("\n%d:%-20s\t",x,matrix[x].gene);
                    printf("%g\t",matrix[x].factor_chrn);
                    printf("%d\t",matrix[x].nsamples);
                    printf("%d\t",matrix[x].noutgroups);
                    printf("%.2f\t",matrix[x].nsites);
                    printf("%d\t",matrix[x].biallsites);
                    printf("%g\t",matrix[x].theta_wat*((double)1/matrix[x].factor_chrn));
                    printf("%g\t",matrix[x].theta_taj*((double)1/matrix[x].factor_chrn));
                    printf("%g\t",matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn));
                }
                printf("\n\n");
                
                if(file_output) {
                    fputs("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					fputs(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples, outgroups, sites, ",file_output);
                    fputs("segregating biallelic sites,\n",file_output);
                    fputs(" and locus estimations of variability from Watterson, Tajima and Fu and Li.\n",file_output);
					fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n",file_output);
                    fputs("\n\nname_loci\tfactor_chrn\t#samples\t#outg\t#sites\t",file_output);
                    fputs("Seg_sites\t",file_output);
                    fputs("theta_wat\ttheta_taj\t",file_output);
                    fputs("theta_fuli\n",file_output);
                    for(x=0;x<*n_loci;x++) {
                        fprintf(file_output,"\n%d:%-20s\t",x,matrix[x].gene);
                        fprintf(file_output,"%g\t",matrix[x].factor_chrn);
                        fprintf(file_output,"%d\t",matrix[x].nsamples);
                        fprintf(file_output,"%d\t",matrix[x].noutgroups);
                        fprintf(file_output,"%.2f\t",matrix[x].nsites);
                        fprintf(file_output,"%d\t",matrix[x].biallsites);
                        fprintf(file_output,"%g\t",matrix[x].theta_wat*((double)1/matrix[x].factor_chrn));
                        fprintf(file_output,"%g\t",matrix[x].theta_taj*((double)1/matrix[x].factor_chrn));
                        fprintf(file_output,"%g\t",matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn));
                    }
                    fputs("\n\n",file_output);
                }
            }    
            break;
        case '3':
            if(*outgroup) {
                ss = 0;
                for(y=0;y<*n_loci;y++) ss += matrix[y].shared;
                printf("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n");
                printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				printf(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples, outgroups, sites, ");
                printf("segregating biallelic sites\n");
                if(ss > 0) printf(" not including shared polymorphisms (nsh) or including them (sh), \n");
                printf(" nucleotide divergence (not corrected and with Jukes and Cantor correction)\n");
                printf(" and nucleotide estimations of variability from Watterson, Tajima, Fu and Li,\n");
                printf(" Fay and Wu and Zeng et al.\n");
                printf(" NOTE that all theta estimates are corrected for chromosome population size.\n");
                printf("\n\nname_loci\tfactor_chrn\t#samples\t#outg\t#sites\t");
                if(ss > 0) printf("S(nsh)\tS(sh)\t");
                else printf("Seg_sites\tdiv\tdivJC\t");
                if(ss > 0) printf("theta_wat(nsh)\ttheta_wat(sh)\ttheta_taj(nsh)\ttheta_taj(sh)\t");
                else printf("theta_wat\ttheta_taj\t");
                printf("theta_fuli\ttheta_faywu\ttheta_zeng\n");

                for(x=0;x<*n_loci;x++) {
                    printf("\n%d:%-20s\t",x,matrix[x].gene);
                    printf("%g\t",matrix[x].factor_chrn);
                    printf("%d\t",matrix[x].nsamples);
                    printf("%d\t",matrix[x].noutgroups);
                    printf("%.2f\t",matrix[x].nsites);
                    printf("%d\t",matrix[x].biallsites);
                    if(ss > 0) {
                        printf("%d\t",matrix[x].biallsitesn);
                    }
                    printf("%g\t",matrix[x].ndivergence/(double)matrix[x].nsites);
                    if((jc = matrix[x].ndivergence/(double)matrix[x].nsites) < (double)0.75) {
                        jc = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * jc);
                        printf("%g\t",jc);
                    }
                    else printf("na\t");
                    printf("%g\t",matrix[x].theta_wat/(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn));
                    if(ss > 0) printf("%g\t",matrix[x].theta_watn/matrix[x].nsites*((double)1/matrix[x].factor_chrn));
                    printf("%g\t",matrix[x].theta_taj/(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn));
                    if(ss > 0) printf("%g\t",matrix[x].theta_tajn/matrix[x].nsites*((double)1/matrix[x].factor_chrn));
                    printf("%g\t",matrix[x].theta_fuli/(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn));
                    printf("%g\t",matrix[x].theta_fw/(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn));
                    printf("%g\t",matrix[x].theta_L/(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn));
                }
                printf("\n\n");
                
                if(file_output) {
                    ss = 0;
                    for(y=0;y<*n_loci;y++) ss += matrix[y].shared;
                    fputs("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					fputs(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples, outgroups, sites, ",file_output);
                    fputs("segregating biallelic sites\n",file_output);
                    if(ss > 0) fputs(" not including shared polymorphisms (nsh) or including them (sh), \n",file_output);
                    fputs(" nucleotide divergence (not corrected and with Jukes and Cantor correction)\n",file_output);
                    fputs(" and nucleotide estimations of variability from Watterson, Tajima, Fu and Li,\n",file_output);
                    fputs(" Fay and Wu anf Zeng et al.\n",file_output);
					fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n",file_output);
                    fputs("\n\nname_loci\tfactor_chrn\t#samples\t#outg\t#sites\t",file_output);
                    if(ss > 0) fputs("S(nsh)\tS(sh)\t",file_output);
                    else fputs("Seg_sites\tdiv\tdivJC\t",file_output);
                    if(ss > 0) fputs("theta_wat(nsh)\ttheta_wat(sh)\ttheta_taj(nsh)\ttheta_taj(sh)\t",file_output);
                    else fputs("theta_wat\ttheta_taj\t",file_output);
                    fputs("theta_fuli\ttheta_faywu\ttheta_zeng\n",file_output);
                    for(x=0;x<*n_loci;x++) {
                        fprintf(file_output,"\n%d:%-20s\t",x,matrix[x].gene);
                        fprintf(file_output,"%g\t",matrix[x].factor_chrn);
                        fprintf(file_output,"%d\t",matrix[x].nsamples);
                        fprintf(file_output,"%d\t",matrix[x].noutgroups);
                        fprintf(file_output,"%.2f\t",matrix[x].nsites);
                        fprintf(file_output,"%d\t",matrix[x].biallsites);
                        if(ss > 0) {
                            fprintf(file_output,"%d\t",matrix[x].biallsitesn);
                            fprintf(file_output,"%d\t",matrix[x].shared);
                        }
                        fprintf(file_output,"%g\t",matrix[x].ndivergence/(double)matrix[x].nsites);
                        if((jc = matrix[x].ndivergence/(double)matrix[x].nsites) < (double)0.75) {
                            jc = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * jc);
                            fprintf(file_output,"%g\t",jc);
                        }
                        else fputs("na\t",file_output);
                        fprintf(file_output,"%g\t",matrix[x].theta_wat/(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn));
                        if(ss > 0) fprintf(file_output,"%g\t",matrix[x].theta_watn/matrix[x].nsites*((double)1/matrix[x].factor_chrn));
                        fprintf(file_output,"%g\t",matrix[x].theta_taj/(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn));
                        if(ss > 0) fprintf(file_output,"%g\t",matrix[x].theta_tajn/matrix[x].nsites*((double)1/matrix[x].factor_chrn));
                        fprintf(file_output,"%g\t",matrix[x].theta_fuli/(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn));
                        fprintf(file_output,"%g\t",matrix[x].theta_fw/(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn));
                        fprintf(file_output,"%g\t",matrix[x].theta_L/(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn));
                    }
                    fputs("\n\n",file_output);
                }
            }
            else {
                printf("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				printf(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples, outgroups, sites, ");
                printf("segregating biallelic sites,\n");
                printf(" and nucleotide estimations of nucleotide from Watterson, Tajima and Fu and Li.\n");
				printf(" NOTE that all theta estimates are corrected for chromosome population size.\n");
                printf("\n\nname_loci\tfactor_chrn\t#samples\t#outg\t#sites\t");
                printf("Seg_sites\t");
                printf("theta_wat\ttheta_taj\t");
                printf("theta_fuli\n");
                for(x=0;x<*n_loci;x++) {
                    printf("\n%d:%-20s\t",x,matrix[x].gene);
                    printf("%g\t",matrix[x].factor_chrn);
                    printf("%d\t",matrix[x].nsamples);
                    printf("%d\t",matrix[x].noutgroups);
                    printf("%.2f\t",matrix[x].nsites);
                    printf("%d\t",matrix[x].biallsites);
                    printf("%g\t",matrix[x].theta_wat/matrix[x].nsites*((double)1/matrix[x].factor_chrn));
                    printf("%g\t",matrix[x].theta_taj/matrix[x].nsites*((double)1/matrix[x].factor_chrn));
                    printf("%g\t",matrix[x].theta_fulin/matrix[x].nsites*((double)1/matrix[x].factor_chrn));
                }
                printf("\n\n");
                
                if(file_output) {
                    fputs("\n\n Display detailed statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					fputs(" The table contains the name of the locus and the number of samples, outgroups, sites, ",file_output);
                    fputs("segregating biallelic sites,\n",file_output);
                    fputs(" and nucleotide estimations of variability from Watterson, Tajima and Fu and Li.\n",file_output);
					fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n",file_output);
                    fputs("\n\nname_loci\tfactor_chrn\t#samples\t#outg\t#sites\t",file_output);
                    fputs("Seg_sites\t",file_output);
                    fputs("theta_wat\ttheta_taj\t",file_output);
                    fputs("theta_fuli\n",file_output);
                    for(x=0;x<*n_loci;x++) {
                        fprintf(file_output,"\n%d:%-20s\t",x,matrix[x].gene);
                        fprintf(file_output,"%g\t",matrix[x].factor_chrn);
                        fprintf(file_output,"%d\t",matrix[x].nsamples);
                        fprintf(file_output,"%d\t",matrix[x].noutgroups);
                        fprintf(file_output,"%.2f\t",matrix[x].nsites);
                        fprintf(file_output,"%d\t",matrix[x].biallsites);
                        fprintf(file_output,"%g\t",matrix[x].theta_wat/matrix[x].nsites*((double)1/matrix[x].factor_chrn));
                        fprintf(file_output,"%g\t",matrix[x].theta_taj/matrix[x].nsites*((double)1/matrix[x].factor_chrn));
                        fprintf(file_output,"%g\t",matrix[x].theta_fulin/matrix[x].nsites*((double)1/matrix[x].factor_chrn));
                    }
                    fputs("\n\n",file_output);
                }
            }    
            break;
        case '4':
            return;
            break;
    }
    return;
}

void open_displntdetgrid(struct statistics *matrix/*,int *observed_data*/,int *outgroup,int *n_loci,char *subset_positions,FILE *file_output)
{
    int x;
    
	if(file_output) {		
		fprintf(file_output,"\n\n     MENU:\n     1. Observed analysis menu:");
		fprintf(file_output,"\n     1.1. Display neutrality tests menu:");
		fprintf(file_output,"\n     1.1.1. Display detailed neutrality tests menu:\n\n");
		fprintf(file_output,"\n     1.1.1.0. Display a table with the tests of neutrality for each locus.\n");
	}

    if(*outgroup) {
        printf("\n\n Display detailed neutral tests: \n\n Multiple hits not included in the analysis.\n");
        printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n");
        if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
		printf(" The table contains the name of the locus, correction factor for chromosome population size, the number of samples for each locus,\n");
        printf(" The number of biallelic sites (not including shared polymorphisms, if observed),\n");
        printf(" Tajima's D test, Fu and Li's D and F tests, Fay and Wu's H test, normalized Fay and Wu H test, Fu's Fs test,\n");
        printf(" Rozas' et a. ZA test, Wall's B and Q tests, Ramos-Onsins & Rozas (R2) test, Zeng E test and Ewens-Watterson test. \n");
        printf(" Shared polymorphisms (if observed) are not included in the above tests of neutrality.\n");
        if(*n_loci > 1 && matrix[0].hka >= (double)0.0) printf(" Also the result of the partial Chi-square value for HKA test in each locus (values corrected by Jukes and Cantor and correcting for factor_chrn).\n");
        printf(" 'na' means the test is not calculated because the number of samples or the number of segregating sites\n");
        printf("were insufficient.\n");
        
        printf("\n\nname_loci\tfactor_chrn\t#samples\tS\t");
        printf("TajD\tFuLiD\tFuLiF\tFayWuH\tFayWuHn\tFuFs\t");
        printf("RozasZA\tWallB\tWallQ\tR2\tE\tEW\t");
        if(*n_loci > 1 && matrix[0].hka >= (double)0.0) printf("partialhka\t");
        printf("\n");
        
        for(x=0;x<*n_loci;x++) {
            printf("\n%d:%-20s\t",x,matrix[x].gene);
			printf("%g\t",matrix[x].factor_chrn);
            printf("%d\t",matrix[x].nsamples);
            printf("%d\t",matrix[x].biallsites);
            if(matrix[x].tajimaD < -9999) printf("na\t");
            else printf("%g\t",matrix[x].tajimaD);
            if(matrix[x].fuliD < -9999) printf("na\t");
            else printf("%g\t",matrix[x].fuliD);
            if(matrix[x].fuliF < -9999) printf("na\t");
            else printf("%g\t",matrix[x].fuliF);
            if(matrix[x].faywuHo == (double) -10000) printf("na\t");
            else printf("%g\t",matrix[x].faywuHo);
            if(matrix[x].faywuH < -9999) printf("na\t");
            else printf("%g\t",matrix[x].faywuH);
            if(matrix[x].fuFs == (double) -10000) printf("na\t");
            else printf("%g\t",matrix[x].fuFs);
            if(matrix[x].rZA < -9999) printf("na\t");
            else printf("%g\t",matrix[x].rZA);
            if(matrix[x].wB < -9999) printf("na\t");
            else printf("%g\t",matrix[x].wB);
            if(matrix[x].wQ < -9999) printf("na\t");
            else printf("%g\t",matrix[x].wQ);
            if(matrix[x].R2 < -9999) printf("na\t");
            else printf("%g\t",matrix[x].R2);
            if(matrix[x].zengE < -9999) printf("na\t");
            else printf("%g\t",matrix[x].zengE);
            if(matrix[x].ewtest < -9999) printf("na\t");
            else printf("%g\t",matrix[x].ewtest);
            if(*n_loci > 1 && matrix[0].hka >= (double)0.0) printf("%g\t",matrix[x].hka);
        }
        printf("\n\n");
        
        if(file_output) {
            fputs("\n\n Display detailed neutral tests: \n\n Multiple hits not included in the analysis.\n",file_output);
            fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n",file_output);
            if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
			fputs(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples for each locus,\n",file_output);
            fputs(" The number of biallelic sites (not including shared polymorphisms, if observed),\n",file_output);
            fputs(" Tajima's D test, Fu and Li's D and F tests, Fay and Wu H test, normalized Fay and Wu H test, Fu's Fs test,\n",file_output);
            fputs(" Rozas' et a. ZA test, Wall's B and Q tests, Ramos-Onsins & Rozas (R2) test, Zeng E test and Ewens-Watterson test.\n",file_output);
            fputs(" Shared polymorphisms (if observed) are not included in the above tests of neutrality.\n",file_output);
            if(*n_loci > 1 && matrix[0].hka >= (double)0.0) fputs(" Also the result of the partial Chi-square value for HKA test in each locus (values corrected by Jukes and Cantor and correcting for factor_chrn).\n",file_output);
            fputs(" 'na' means the test is not calculated because the number of samples or the number of segregating sites\n",file_output);
            fputs(" were insufficient.\n",file_output);
            
            fputs("\n\nname_loci\tfactor_chrn\t#samples\tS\t",file_output);
            fputs("TajD\tFuLiD\tFuLiF\tFayWuH\tFayWuHn\tFuFs\t",file_output);
            fputs("RozasZA\tWallB\tWallQ\tR2\tE\tEW\t",file_output);
            if(*n_loci > 1 && matrix[0].hka >= (double)0.0) fputs("partialhka\t",file_output);
            fputs("\n",file_output);
            
            for(x=0;x<*n_loci;x++) {
                fprintf(file_output,"\n%d:%-20s\t",x,matrix[x].gene);
				fprintf(file_output,"%g\t",matrix[x].factor_chrn);
                fprintf(file_output,"%d\t",matrix[x].nsamples);
                fprintf(file_output,"%d\t",matrix[x].biallsites);
                if(matrix[x].tajimaD < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].tajimaD);
                if(matrix[x].fuliD < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].fuliD);
                if(matrix[x].fuliF < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].fuliF);
                if(matrix[x].faywuHo == (double) -10000) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].faywuHo);
                if(matrix[x].faywuH < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].faywuH);
                if(matrix[x].fuFs == (double) -10000) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].fuFs);
                if(matrix[x].rZA < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].rZA);
                if(matrix[x].wB < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].wB);
                if(matrix[x].wQ < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].wQ);
                if(matrix[x].R2 < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].R2);
                if(matrix[x].zengE < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].zengE);
                if(matrix[x].ewtest < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].ewtest);
                if(*n_loci > 1 && matrix[0].hka >= (double)0.0) fprintf(file_output,"%g\t",matrix[x].hka);
            }
            fputs("\n\n",file_output);
        }
    }
    else {
        printf("\n\n Display detailed neutral tests: \n\n Multiple hits not included in the analysis.\n");
        if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
		printf(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples for each locus,\n");
        printf(" The number of biallelic sites,\n");
        printf(" Tajima's D test, Fu and Li's D* and F* tests, Fu's Fs test,\n");
        printf(" Rozas' et a. ZA test, Wall's B and Q tests, Ramos-Onsins & Rozas (R2) test and Ewens-Watterson test.\n");
        printf(" 'na' means the test is not calculated because the number of samples or the number of segregating sites\n");
        printf(" were insufficient.\n");
        
        printf("\n\nname_loci\tfactor_chrn\t#samples\tS\t");
        printf("TajD\tFuLiD*\tFuLiF*\tFuFs\t");
        printf("RozasZA\tWallB\tWallQ\tR2\tEW\n");
        
        for(x=0;x<*n_loci;x++) {
            printf("\n%d:%-20s\t",x,matrix[x].gene);
			printf("%g\t",matrix[x].factor_chrn);
            printf("%d\t",matrix[x].nsamples);
            printf("%d\t",matrix[x].biallsites);
            if(matrix[x].tajimaD < -9999) printf("na\t");
            else printf("%g\t",matrix[x].tajimaD);
            if(matrix[x].fuliDn < -9999) printf("na\t");
            else printf("%g\t",matrix[x].fuliDn);
            if(matrix[x].fuliFn < -9999) printf("na\t");
            else printf("%g\t",matrix[x].fuliFn);
            if(matrix[x].fuFs == (double) -10000) printf("na\t");
            else printf("%g\t",matrix[x].fuFs);
            if(matrix[x].rZA < -9999) printf("na\t");
            else printf("%g\t",matrix[x].rZA);
            if(matrix[x].wB < -9999) printf("na\t");
            else printf("%g\t",matrix[x].wB);
            if(matrix[x].wQ < -9999) printf("na\t");
            else printf("%g\t",matrix[x].wQ);
            if(matrix[x].R2 < -9999) printf("na\t");
            else printf("%g\t",matrix[x].R2);
            if(matrix[x].ewtest < -9999) printf("na\t");
            else printf("%g\t",matrix[x].ewtest);
        }
        printf("\n\n");
        
        if(file_output) {
            fputs("\n\n Display detailed neutral tests: \n\n Multiple hits not included in the analysis.\n",file_output);
            if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
			fputs(" The table contains the name of the locus, correction factor for chromosome population size and the number of samples for each locus,\n",file_output);
            fputs(" The number of biallelic sites,\n",file_output);
            fputs(" Tajima's D test, Fu and Li's D* and F* tests, Fu's Fs test,\n",file_output);
            fputs(" Rozas' et a. ZA test, Wall's B and Q tests, Ramos-Onsins & Rozas (R2) test and Ewens-Watterson test\n",file_output);
            fputs(" 'na' means the tests is not calculated because the number of samples or the number of segregating sites\n",file_output);
            fputs(" were insufficient.\n",file_output);
            
            fputs("\n\nname_loci\tfactor_chrn\t#samples\tS\t",file_output);
            fputs("TajD\tFuLiD*\tFuLiF*\tFuFs\t",file_output);
            fputs("RozasZA\tWallB\tWallQ\tR2\tEW\n",file_output);
            
            for(x=0;x<*n_loci;x++) {
                fprintf(file_output,"\n%d:%-20s\t",x,matrix[x].gene);
				fprintf(file_output,"%g\t",matrix[x].factor_chrn);
                fprintf(file_output,"%d\t",matrix[x].nsamples);
                fprintf(file_output,"%d\t",matrix[x].biallsites);
                if(matrix[x].tajimaD < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].tajimaD);
                if(matrix[x].fuliDn < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].fuliDn);
                if(matrix[x].fuliFn < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].fuliFn);
                if(matrix[x].fuFs == (double) -10000) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].fuFs);
                if(matrix[x].rZA < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].rZA);
                if(matrix[x].wB < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].wB);
                if(matrix[x].wQ < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].wQ);
                if(matrix[x].R2 < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].R2);
                if(matrix[x].ewtest < -9999) fputs("na\t",file_output);
                else fprintf(file_output,"%g\t",matrix[x].ewtest);
            }
            fputs("\n\n",file_output);
        }
    }    
    return;
}

void open_displstmulo(struct statmulo *matrixml/*,int *observed_data*/,int *outgroup/*,int *n_loci */,char *subset_positions,FILE *file_output)
{
    char k[1];
    double var;
    double variances(double,double,int);

    #if COMMAND_LINE        
    if(file_output) fflush(file_output);
    
    printf("\n\nMENU 1. Observed analysis menu:");
    printf("\n     1.0. Display statistics menu:");
    printf("\n     1.0.0. Display multilocus statistics:\n\n");
    
    printf("CHOOSE:\n");
    printf(" 0 - Display general statistics.\n");
    printf(" 1 - Display other linkage related statistics.\n");
    printf(" 2 - Display estimates of total variability.\n");
    printf(" 3 - Display estimates of nucleotide variability.\n");
    printf(" 4 - Back to Main menu.\n");

    if(file_output) {
        
        fprintf(file_output,"\n\n     MENU:\n     1. Observed analysis menu:");
        fprintf(file_output,"\n     1.0. Display statistics menu:");
        fprintf(file_output,"\n     1.0.0. Display multilocus statistics:\n\n");
        
        fprintf(file_output," 0 - Display general statistics.\n");
        fprintf(file_output," 1 - Display other linkage related statistics.\n");
        fprintf(file_output," 2 - Display estimates of total variability.\n");
        fprintf(file_output," 3 - Display estimates of nucleotide variability.\n");
        fprintf(file_output," 4 - Back to Main menu.\n");
    }
    #endif
    
    do *k = getchar();
    while(*k<'0' || *k>'4');
    
    switch(*k) {
        case '0':
            if(*outgroup) {
                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
                printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				/* nombre de loci, then SUM, AVERAGE and VARIANCE */
                printf("Total number of loci: %d\n\n",matrixml[0].nloci);
                printf("Total number of sites: %.2f\n",matrixml[0].Snsites);
                printf("Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                printf("Total number of multiple hits: %ld\n",matrixml[0].Snmhits);
                printf("Average number of multiple hits: %g\n\n",(double)matrixml[0].Snmhits/(double)matrixml[0].nloci);
                if(matrixml[0].Stransversions) {
					printf("Total s/v ratio: %.2f\n\n",(double)matrixml[0].Stransitions/(double)matrixml[0].Stransversions);
				}
                printf("Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                printf("Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                if(var > (double) -10000) printf("Variance number of biallelic sites per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) {
                    printf("Total number of biallelic sites (incl. shared): %ld\n",matrixml[0].Sbiallsitesn);
                    printf("Average number of biallelic sites per loci (incl. shared): %g\n",(double)matrixml[0].Sbiallsitesn/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2biallsitesn,(double)matrixml[0].Sbiallsitesn,matrixml[0].nloci);
                    if(var > (double)-10000) printf("Variance number of biallelic sites per loci (incl. shared): %g\n\n",var);
                    printf("Total number of shared polymorphisms: %ld\n",matrixml[0].Sshared);
                    printf("Average number of shared polymorphisms per loci: %g\n",(double)matrixml[0].Sshared/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2shared,(double)matrixml[0].Sshared,matrixml[0].nloci);
                    if(var > (double)-10000) printf("Variance number of shared polymorphisms per loci: %g\n\n",var);
                }
                printf("Total number of fixed sites: %ld\n",matrixml[0].Sfixed);
                printf("Average number of fixed sites per loci: %g\n",(double)matrixml[0].Sfixed/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2fixed,(double)matrixml[0].Sfixed,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance number of fixed sites per loci: %g\n\n",var);
                printf("Total divergence: %g\n",matrixml[0].Sndivergence);
                printf("Averaged divergence per loci: %g\n",(double)matrixml[0].Sndivergence/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2ndivergence,(double)matrixml[0].Sndivergence,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of divergence per loci: %g\n\n",var);
                printf("Average number of haplotypes per loci: %g\n",(double)matrixml[0].Snhapl/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2nhapl,(double)matrixml[0].Snhapl,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of the number of haplotypes per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) printf("Not shared polymorphisms included in counting the number of haplotypes");
                printf("Average number of haplotypes/sample size per loci: %g\n",(double)matrixml[0].Snhaplsam/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2nhaplsam,(double)matrixml[0].Snhaplsam,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of the number of haplotypes/sample size per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) printf("Not shared polymorphisms included in counting the number of haplotypes");
                printf("Average of haplotype diversity per loci: %g\n",(double)matrixml[0].Shapldiv/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2hapldiv,(double)matrixml[0].Shapldiv,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of haplotype diversity per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) printf("Not shared polymorphisms included in counting the number of haplotypes");
                /*
				 printf("Average of C per loci: %g\n",(double)matrixml[0].SRvpi/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2Rvpi,(double)matrixml[0].SRvpi,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of C per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) printf("Not shared polymorphisms included in counting C");
                */
				printf("Average of Rm per loci: %g\n",(double)matrixml[0].SRm/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2Rm,(double)matrixml[0].SRm,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of Rm per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) printf("Not shared polymorphisms included in counting Rm");

                if(file_output) {
                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					/* nombre de loci, then SUM, AVERAGE and VARIANCE */
                    fprintf(file_output,"Total number of loci: %d\n\n",matrixml[0].nloci);
                    fprintf(file_output,"Total number of sites: %.2f\n",matrixml[0].Snsites);
                    fprintf(file_output,"Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
					fprintf(file_output,"Total number of multiple hits: %ld\n",matrixml[0].Snmhits);
					fprintf(file_output,"Average number of multiple hits: %g\n\n",(double)matrixml[0].Snmhits/(double)matrixml[0].nloci);
					if(matrixml[0].Stransversions) {
						fprintf(file_output,"Total s/v ratio: %.2f\n\n",(double)matrixml[0].Stransitions/(double)matrixml[0].Stransversions);
					}
                    fprintf(file_output,"Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                    fprintf(file_output,"Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance number of biallelic sites per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0) {
                        fprintf(file_output,"Total number of biallelic sites (incl. shared): %ld\n",matrixml[0].Sbiallsitesn);
                        fprintf(file_output,"Average number of biallelic sites per loci (incl. shared): %g\n",(double)matrixml[0].Sbiallsitesn/(double)matrixml[0].nloci);
                        var = variances((double)matrixml[0].S2biallsitesn,(double)matrixml[0].Sbiallsitesn,matrixml[0].nloci);
                        if(var > (double)-10000) fprintf(file_output,"Variance number of biallelic sites per loci (incl. shared): %g\n\n",var);
                        fprintf(file_output,"Total number of shared polymorphisms: %ld\n",matrixml[0].Sshared);
                        fprintf(file_output,"Average number of shared polymorphisms per loci: %g\n",(double)matrixml[0].Sshared/(double)matrixml[0].nloci);
                        var = variances((double)matrixml[0].S2shared,(double)matrixml[0].Sshared,matrixml[0].nloci);
                        if(var > (double)-10000) fprintf(file_output,"Variance number of shared polymorphisms per loci: %g\n\n",var);
                    }
                    fprintf(file_output,"Total number of fixed sites: %ld\n",matrixml[0].Sfixed);
                    fprintf(file_output,"Average number of fixed sites per loci: %g\n",(double)matrixml[0].Sfixed/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2fixed,(double)matrixml[0].Sfixed,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance number of fixed sites per loci: %g\n\n",var);
                    fprintf(file_output,"Total divergence: %g\n",matrixml[0].Sndivergence);
                    fprintf(file_output,"Averaged divergence per loci: %g\n",(double)matrixml[0].Sndivergence/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2ndivergence,(double)matrixml[0].Sndivergence,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of divergence per loci: %g\n\n",var);
                    fprintf(file_output,"Average number of haplotypes per loci: %g\n",(double)matrixml[0].Snhapl/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2nhapl,(double)matrixml[0].Snhapl,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of the number of haplotypes per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0)
                        fputs("Not shared polymorphisms included in counting the number of haplotypes",file_output);
                    fprintf(file_output,"Average number of haplotypes/sample size per loci: %g\n",(double)matrixml[0].Snhaplsam/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2nhaplsam,(double)matrixml[0].Snhaplsam,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of the number of haplotypes/sample size per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0)
                        fputs("Not shared polymorphisms included in counting the number of haplotypes",file_output);
                    fprintf(file_output,"Average of haplotype diversity per loci: %g\n",(double)matrixml[0].Shapldiv/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2hapldiv,(double)matrixml[0].Shapldiv,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of haplotype diversity per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0)
                        fputs("Not shared polymorphisms included in counting the number of haplotypes",file_output);
					fprintf(file_output,"Average of C per loci: %g\n",(double)matrixml[0].SRvpi/(double)matrixml[0].nloci);
					var = variances((double)matrixml[0].S2Rvpi,(double)matrixml[0].SRvpi,matrixml[0].nloci);
					if(var > (double)-10000) fprintf(file_output,"Variance of C per loci: %g\n\n",var);
					if(matrixml[0].Sshared > 0) fprintf(file_output,"Not shared polymorphisms included in counting C");
					fprintf(file_output,"Average of Rm per loci: %g\n",(double)matrixml[0].SRm/(double)matrixml[0].nloci);
					var = variances((double)matrixml[0].S2Rm,(double)matrixml[0].SRm,matrixml[0].nloci);
					if(var > (double)-10000) fprintf(file_output,"Variance of Rm per loci: %g\n\n",var);
					if(matrixml[0].Sshared > 0) fprintf(file_output,"Not shared polymorphisms included in counting Rm");
                }
            }
            else {
                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
                printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				/* nombre de loci, then SUM, AVERAGE and VARIANCE */
                printf("Total number of loci: %d\n\n",matrixml[0].nloci);
                printf("Total number of sites: %.2f\n",matrixml[0].Snsites);
                printf("Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                printf("Total number of multiple hits: %ld\n",matrixml[0].Snmhits);
                printf("Average number of multiple hits: %g\n\n",(double)matrixml[0].Snmhits/(double)matrixml[0].nloci);
                if(matrixml[0].Stransversions) {
					printf("Total s/v ratio: %.2f\n\n",(double)matrixml[0].Stransitions/(double)matrixml[0].Stransversions);
				}
                printf("Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                printf("Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance number of biallelic sites per loci: %g\n\n",var);
                printf("Average number of haplotypes per loci: %g\n",(double)matrixml[0].Snhapl/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2nhapl,(double)matrixml[0].Snhapl,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of the number of haplotypes per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) printf("Not shared polymorphisms included in counting the number of haplotypes");
                else printf("\n");
                printf("Average number of haplotypes/sample size per loci: %g\n",(double)matrixml[0].Snhaplsam/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2nhaplsam,(double)matrixml[0].Snhaplsam,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of the number of haplotypes/sample size per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) printf("Not shared polymorphisms included in counting the number of haplotypes");
                printf("Average of haplotype diversity per loci: %g\n",(double)matrixml[0].Shapldiv/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2hapldiv,(double)matrixml[0].Shapldiv,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of haplotype diversity per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) printf("Not shared polymorphisms included in counting the number of haplotypes");
                printf("Average of C per loci: %g\n",(double)matrixml[0].SRvpi/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2Rvpi,(double)matrixml[0].SRvpi,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of C per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) printf("Not shared polymorphisms included in counting C");
                printf("Average of Rm per loci: %g\n",(double)matrixml[0].SRm/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2Rm,(double)matrixml[0].SRm,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of Rm per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) printf("Not shared polymorphisms included in counting Rm");
                else printf("\n");
            
                if(file_output) {
                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					/* nombre de loci, then SUM, AVERAGE and VARIANCE */
                    fprintf(file_output,"Total number of loci: %d\n\n",matrixml[0].nloci);
                    fprintf(file_output,"Total number of sites: %.2f\n",matrixml[0].Snsites);
                    fprintf(file_output,"Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
					fprintf(file_output,"Total number of multiple hits: %ld\n",matrixml[0].Snmhits);
					fprintf(file_output,"Average number of multiple hits: %g\n\n",(double)matrixml[0].Snmhits/(double)matrixml[0].nloci);
					if(matrixml[0].Stransversions) {
						fprintf(file_output,"Total s/v ratio: %.2f\n\n",(double)matrixml[0].Stransitions/(double)matrixml[0].Stransversions);
					}
                    fprintf(file_output,"Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                    fprintf(file_output,"Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance number of biallelic sites per loci: %g\n\n",var);
                    fprintf(file_output,"Average number of haplotypes per loci: %g\n",(double)matrixml[0].Snhapl/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2nhapl,(double)matrixml[0].Snhapl,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of the number of haplotypes per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0)
                        fputs("Not shared polymorphisms included in counting the number of haplotypes",file_output);
                    else fputs("\n",file_output);
                    fprintf(file_output,"Average number of haplotypes/sample size per loci: %g\n",(double)matrixml[0].Snhaplsam/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2nhaplsam,(double)matrixml[0].Snhaplsam,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of the number of haplotypes/sample size per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0)
                        fputs("Not shared polymorphisms included in counting the number of haplotypes",file_output);
                    fprintf(file_output,"Average of haplotype diversity per loci: %g\n",(double)matrixml[0].Shapldiv/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2hapldiv,(double)matrixml[0].Shapldiv,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of haplotype diversity per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0)
                        fputs("Not shared polymorphisms included in counting the number of haplotypes",file_output);
                    else fputs("\n",file_output);
					fprintf(file_output,"Average of C per loci: %g\n",(double)matrixml[0].SRvpi/(double)matrixml[0].nloci);
					var = variances((double)matrixml[0].S2Rvpi,(double)matrixml[0].SRvpi,matrixml[0].nloci);
					if(var > (double)-10000) fprintf(file_output,"Variance of C per loci: %g\n\n",var);
					if(matrixml[0].Sshared > 0) fprintf(file_output,"Not shared polymorphisms included in counting C");
					fprintf(file_output,"Average of Rm per loci: %g\n",(double)matrixml[0].SRm/(double)matrixml[0].nloci);
					var = variances((double)matrixml[0].S2Rm,(double)matrixml[0].SRm,matrixml[0].nloci);
					if(var > (double)-10000) fprintf(file_output,"Variance of Rm per loci: %g\n\n",var);
					if(matrixml[0].Sshared > 0) fprintf(file_output,"Not shared polymorphisms included in counting Rm");
                }
            }
            break;
        case '1':
            if(*outgroup) {
                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
                printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				if(matrixml[0].Sshared > 0)
                    printf("Not shared polymorphisms included in the statistics B, Q and ZA.\n\n");
                /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                printf("Total number of loci: %d\n\n",matrixml[0].nloci);
                printf("Total number of sites: %.2f\n",matrixml[0].Snsites);
                printf("Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                printf("Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                printf("Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance number of biallelic sites per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) {
                    printf("Total number of biallelic sites (incl. shared): %ld\n",matrixml[0].Sbiallsitesn);
                    printf("Average number of biallelic sites per loci (incl. shared): %g\n",(double)matrixml[0].Sbiallsitesn/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2biallsitesn,(double)matrixml[0].Sbiallsitesn,matrixml[0].nloci);
                    if(var > (double)-10000) printf("Variance number of biallelic sites per loci (incl. shared): %g\n\n",var);
                }
                printf("Total value of b (Wall's B not divided): %ld\n",matrixml[0].Sb);
                printf("Average number of b (Wall's B not divided) per loci: %g\n",(double)matrixml[0].Sb/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2b,(double)matrixml[0].Sb,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of b (Wall's B not divided) per loci: %g\n\n",var);
                printf("Total value of q (Wall's Q not divided): %ld\n",matrixml[0].Sq);
                printf("Average number of q (Wall's Q not divided) per loci: %g\n",(double)matrixml[0].Sq/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2q,(double)matrixml[0].Sq,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of q (Wall's Q not divided) per loci: %g\n\n",var);
                printf("Total value of za (Rozas' et al. ZA but not divided): %g\n",matrixml[0].Sza);
                printf("Average number of za (Rozas' et al. ZA but not divided) per loci: %g\n",(double)matrixml[0].Sza/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2za,(double)matrixml[0].Sza,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of za (Rozas' et al. ZA but not divided) per loci: %g\n\n",var);
                else printf("\n");

                if(file_output) {
                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					if(matrixml[0].Sshared > 0)
                        fputs("Not shared polymorphisms included in the statistics B, Q and ZA.\n\n",file_output);
                    /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                    fprintf(file_output,"Total number of loci: %d\n\n",matrixml[0].nloci);
                    fprintf(file_output,"Total number of sites: %.2f\n",matrixml[0].Snsites);
                    fprintf(file_output,"Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                    fprintf(file_output,"Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                    fprintf(file_output,"Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance number of biallelic sites per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0) {
                        fprintf(file_output,"Total number of biallelic sites (incl. shared): %ld\n",matrixml[0].Sbiallsitesn);
                        fprintf(file_output,"Average number of biallelic sites per loci (incl. shared): %g\n",(double)matrixml[0].Sbiallsitesn/(double)matrixml[0].nloci);
                        var = variances((double)matrixml[0].S2biallsitesn,(double)matrixml[0].Sbiallsitesn,matrixml[0].nloci);
                        if(var > (double)-10000) fprintf(file_output,"Variance number of biallelic sites per loci (incl. shared): %g\n\n",var);
                    }
                    fprintf(file_output,"Total value of b (Wall's B not divided): %ld\n",matrixml[0].Sb);
                    fprintf(file_output,"Average number of b (Wall's B not divided) per loci: %g\n",(double)matrixml[0].Sb/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2b,(double)matrixml[0].Sb,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of b (Wall's B not divided) per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of q (Wall's Q not divided): %ld\n",matrixml[0].Sq);
                    fprintf(file_output,"Average number of q (Wall's Q not divided) per loci: %g\n",(double)matrixml[0].Sq/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2q,(double)matrixml[0].Sq,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of q (Wall's Q not divided) per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of za (Rozas' et al. ZA but not divided): %g\n",matrixml[0].Sza);
                    fprintf(file_output,"Average number of za (Rozas' et al. ZA but not divided) per loci: %g\n",(double)matrixml[0].Sza/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2za,(double)matrixml[0].Sza,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of za (Rozas' et al. ZA but not divided) per loci: %g\n\n",var);
                    else fputs("\n",file_output); 
                }
            }
            else {
                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
                printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				/* nombre de loci, then SUM, AVERAGE and VARIANCE */
                printf("Total number of loci: %d\n\n",matrixml[0].nloci);
                printf("Total number of sites: %.2f\n",matrixml[0].Snsites);
                printf("Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                printf("Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                printf("Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance number of biallelic sites per loci: %g\n\n",var);
                printf("Total value of b (Wall's B not divided): %ld\n",matrixml[0].Sb);
                printf("Average number of b (Wall's B not divided) per loci: %g\n",(double)matrixml[0].Sb/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2b,(double)matrixml[0].Sb,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of b (Wall's B not divided) per loci: %g\n\n",var);
                printf("Total value of q (Wall's Q not divided): %ld\n",matrixml[0].Sq);
                printf("Average number of q (Wall's Q not divided) per loci: %g\n",(double)matrixml[0].Sq/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2q,(double)matrixml[0].Sq,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of q (Wall's Q not divided) per loci: %g\n\n",var);
                printf("Total value of za (Rozas' et al. ZA but not divided): %g\n",matrixml[0].Sza);
                printf("Average number of za (Rozas' et al. ZA but not divided) per loci: %g\n",(double)matrixml[0].Sza/(double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2za,(double)matrixml[0].Sza,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of za (Rozas' et al. ZA but not divided) per loci: %g\n\n",var);
                else printf("\n");
            
                if(file_output) {
                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					/* nombre de loci, then SUM, AVERAGE and VARIANCE */
                    fprintf(file_output,"Total number of loci: %d\n\n",matrixml[0].nloci);
                    fprintf(file_output,"Total number of sites: %.2f\n",matrixml[0].Snsites);
                    fprintf(file_output,"Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                    fprintf(file_output,"Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                    fprintf(file_output,"Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance number of biallelic sites per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of b (Wall's B not divided): %ld\n",matrixml[0].Sb);
                    fprintf(file_output,"Average number of b (Wall's B not divided) per loci: %g\n",(double)matrixml[0].Sb/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2b,(double)matrixml[0].Sb,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of b (Wall's B not divided) per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of q (Wall's Q not divided): %ld\n",matrixml[0].Sq);
                    fprintf(file_output,"Average number of q (Wall's Q not divided) per loci: %g\n",(double)matrixml[0].Sq/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2q,(double)matrixml[0].Sq,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of q (Wall's Q not divided) per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of za (Rozas' et al. ZA but not divided): %g\n",matrixml[0].Sza);
                    fprintf(file_output,"Average number of za (Rozas' et al. ZA but not divided) per loci: %g\n",(double)matrixml[0].Sza/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2za,(double)matrixml[0].Sza,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of za (Rozas' et al. ZA but not divided) per loci: %g\n\n",var);
                    else fputs("\n",file_output);
                }
            }
            break;
        case '2':
            if(*outgroup) {
                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
                printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n");
				if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				printf(" NOTE that all theta estimates are corrected for chromosome population size.\n\n");
                /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                printf("Total number of loci: %d\n\n",matrixml[0].nloci);
                printf("Total number of sites: %.2f\n",matrixml[0].Snsites);
                printf("Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                printf("Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                printf("Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance number of biallelic sites per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) {
                    printf("Total number of biallelic sites (incl. shared): %ld\n",matrixml[0].Sbiallsitesn);
                    printf("Average number of biallelic sites per loci (incl. shared): %g\n",(double)matrixml[0].Sbiallsitesn/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2biallsitesn,(double)matrixml[0].Sbiallsitesn,matrixml[0].nloci);
                    if(var > (double)-10000) printf("Variance number of biallelic sites per loci (incl. shared): %g\n\n",var);
                }
                printf("Total divergence: %g\n",matrixml[0].Sndivergence);
                printf("Averaged divergence per loci: %g\n",(double)matrixml[0].Sndivergence/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2ndivergence,(double)matrixml[0].Sndivergence,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of divergence per loci: %g\n\n",var);
                printf("Total divergence corrected by Jukes and Cantor: %g\n",matrixml[0].Sndivergencejc);
                printf("Averaged divergence per loci corrected by Jukes and Cantor: %g\n",(double)matrixml[0].Sndivergencejc / (double)matrixml[0].nldivjc);
                var = variances((double)matrixml[0].S2ndivergencejc,(double)matrixml[0].Sndivergencejc,matrixml[0].nldivjc);
                if(var > (double)-10000) printf("Variance of divergence per loci corrected by Jukes and Cantor: %g\n\n",var);
                printf("Total value of theta (Watterson): %g\n",matrixml[0].Stheta_wat);
                printf("Average value of theta (Watterson) per loci: %g\n",(double)matrixml[0].Stheta_wat/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_wat,(double)matrixml[0].Stheta_wat,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Watterson) per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) {
                    printf("Total value of theta (Watterson) including shared pol.: %g\n",matrixml[0].Stheta_watn);
                    printf("Average value of theta (Watterson) per loci including shared pol.: %g\n",(double)matrixml[0].Stheta_watn/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_watn,(double)matrixml[0].Stheta_watn,matrixml[0].nloci);
                    if(var > (double)-10000) printf("Variance of theta (Watterson) per loci including shared pol.: %g\n\n",var);
                }
                printf("Total value of theta (Tajima): %g\n",matrixml[0].Stheta_taj);
                printf("Average value of theta (Tajima) per loci: %g\n",(double)matrixml[0].Stheta_taj/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_taj,(double)matrixml[0].Stheta_taj,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Tajima) per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) {
                    printf("Total value of theta (Tajima) including shared pol.: %g\n",matrixml[0].Stheta_tajn);
                    printf("Average value of theta (Tajima) per loci including shared pol.: %g\n",(double)matrixml[0].Stheta_tajn/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_tajn,(double)matrixml[0].Stheta_tajn,matrixml[0].nloci);
                    if(var > (double)-10000) printf("Variance of theta (Tajima) per loci including shared pol.: %g\n\n",var);
                }
                printf("Total value of theta (Fu & Li): %g\n",matrixml[0].Stheta_fuli);
                printf("Average value of theta (Fu & Li) per loci: %g\n",(double)matrixml[0].Stheta_fuli/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_fuli,(double)matrixml[0].Stheta_fuli,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Fu & Li) per loci: %g\n\n",var);
                printf("Total value of theta (Fay & Wu): %g\n",matrixml[0].Stheta_fw);
                printf("Average value of theta (Fay & Wu) per loci: %g\n",(double)matrixml[0].Stheta_fw/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_fw,(double)matrixml[0].Stheta_fw,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Fay & Wu) per loci: %g\n\n",var);
                printf("Total value of theta (Zeng et al.): %g\n",matrixml[0].Stheta_L);
                printf("Average value of theta (Zeng et al.) per loci: %g\n",(double)matrixml[0].Stheta_L/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_L,(double)matrixml[0].Stheta_L,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Zeng et al.) per loci: %g\n\n",var);
                else printf("\n");

                if(file_output) {
                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n\n",file_output);
					/* nombre de loci, then SUM, AVERAGE and VARIANCE */
                    fprintf(file_output,"Total number of loci: %d\n\n",matrixml[0].nloci);
                    fprintf(file_output,"Total number of sites: %.2f\n",matrixml[0].Snsites);
                    fprintf(file_output,"Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                    fprintf(file_output,"Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                    fprintf(file_output,"Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance number of biallelic sites per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0) {
                        fprintf(file_output,"Total number of biallelic sites (incl. shared): %ld\n",matrixml[0].Sbiallsitesn);
                        fprintf(file_output,"Average number of biallelic sites per loci (incl. shared): %g\n",(double)matrixml[0].Sbiallsitesn/(double)matrixml[0].nloci);
                        var = variances((double)matrixml[0].S2biallsitesn,(double)matrixml[0].Sbiallsitesn,matrixml[0].nloci);
                        if(var > (double)-10000) fprintf(file_output,"Variance number of biallelic sites per loci (incl. shared): %g\n\n",var);
                    }
                    fprintf(file_output,"Total divergence: %g\n",matrixml[0].Sndivergence);
                    fprintf(file_output,"Averaged divergence per loci: %g\n",(double)matrixml[0].Sndivergence/ (double)matrixml[0].nldivjc);
                    var = variances((double)matrixml[0].S2ndivergence,(double)matrixml[0].Sndivergence,matrixml[0].nldivjc);
                    if(var > (double)-10000) fprintf(file_output,"Variance of divergence per loci: %g\n\n",var);
                    fprintf(file_output,"Total divergence corrected by Jukes and Cantor: %g\n",matrixml[0].Sndivergencejc);
                    fprintf(file_output,"Averaged divergence per loci corrected by Jukes and Cantor: %g\n",(double)matrixml[0].Sndivergencejc / (double)matrixml[0].nldivjc);
                    var = variances((double)matrixml[0].S2ndivergencejc,(double)matrixml[0].Sndivergencejc,matrixml[0].nldivjc);
                    if(var > (double)-10000) fprintf(file_output,"Variance of divergence per loci corrected by Jukes and Cantor: %g\n\n",var);
                    fprintf(file_output,"Total value of theta (Watterson): %g\n",matrixml[0].Stheta_wat);
                    fprintf(file_output,"Average value of theta (Watterson) per loci: %g\n",(double)matrixml[0].Stheta_wat/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_wat,(double)matrixml[0].Stheta_wat,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Watterson) per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0) {
                        fprintf(file_output,"Total value of theta (Watterson) including shared pol.: %g\n",matrixml[0].Stheta_watn);
                        fprintf(file_output,"Average value of theta (Watterson) per loci including shared pol.: %g\n",(double)matrixml[0].Stheta_watn/ (double)matrixml[0].nloci);
                        var = variances((double)matrixml[0].S2theta_watn,(double)matrixml[0].Stheta_watn,matrixml[0].nloci);
                        if(var > (double)-10000) fprintf(file_output,"Variance of theta (Watterson) per loci including shared pol.: %g\n\n",var);
                    }
                    fprintf(file_output,"Total value of theta (Tajima): %g\n",matrixml[0].Stheta_taj);
                    fprintf(file_output,"Average value of theta (Tajima) per loci: %g\n",(double)matrixml[0].Stheta_taj/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_taj,(double)matrixml[0].Stheta_taj,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Tajima) per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0) {
                        fprintf(file_output,"Total value of theta (Tajima) including shared pol.: %g\n",matrixml[0].Stheta_tajn);
                        fprintf(file_output,"Average value of theta (Tajima) per loci including shared pol.: %g\n",(double)matrixml[0].Stheta_tajn/ (double)matrixml[0].nloci);
                        var = variances((double)matrixml[0].S2theta_tajn,(double)matrixml[0].Stheta_tajn,matrixml[0].nloci);
                        if(var > (double)-10000) fprintf(file_output,"Variance of theta (Tajima) per loci including shared pol.: %g\n\n",var);
                    }
                    fprintf(file_output,"Total value of theta (Fu & Li): %g\n",matrixml[0].Stheta_fuli);
                    fprintf(file_output,"Average value of theta (Fu & Li) per loci: %g\n",(double)matrixml[0].Stheta_fuli/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_fuli,(double)matrixml[0].Stheta_fuli,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Fu & Li) per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of theta (Fay & Wu): %g\n",matrixml[0].Stheta_fw);
                    fprintf(file_output,"Average value of theta (Fay & Wu) per loci: %g\n",(double)matrixml[0].Stheta_fw/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_fw,(double)matrixml[0].Stheta_fw,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Fay & Wu) per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of theta (Zeng et al.): %g\n",matrixml[0].Stheta_L);
                    fprintf(file_output,"Average value of theta (Zeng et al.) per loci: %g\n",(double)matrixml[0].Stheta_L/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_L,(double)matrixml[0].Stheta_L,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Zeng et al.) per loci: %g\n\n",var);
                    else fputs("\n",file_output);
                }
            }
            else {
                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
				printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n");
				if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				printf(" NOTE that all theta estimates are corrected for chromosome population size.\n\n");
                /* nombre de loci, then SUM, AVERAGE and VARIANCE */
                printf("Total number of loci: %d\n\n",matrixml[0].nloci);
                printf("Total number of sites: %.2f\n",matrixml[0].Snsites);
                printf("Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                printf("Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                printf("Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                printf("Variance number of biallelic sites per loci: %g\n\n",var);
                printf("Total value of theta (Watterson): %g\n",matrixml[0].Stheta_wat);
                printf("Average value of theta (Watterson) per loci: %g\n",(double)matrixml[0].Stheta_wat/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_wat,(double)matrixml[0].Stheta_wat,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Watterson) per loci: %g\n\n",var);
                printf("Total value of theta (Tajima): %g\n",matrixml[0].Stheta_taj);
                printf("Average value of theta (Tajima) per loci: %g\n",(double)matrixml[0].Stheta_taj/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_taj,(double)matrixml[0].Stheta_taj,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Tajima) per loci: %g\n\n",var);
                printf("Total value of theta (Fu & Li): %g\n",matrixml[0].Stheta_fulin);
                printf("Average value of theta (Fu & Li) per loci: %g\n",(double)matrixml[0].Stheta_fulin/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_fulin,(double)matrixml[0].Stheta_fulin,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Fu & Li) per loci: %g\n\n",var);
                else printf("\n");
            
                if(file_output) {
                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n\n",file_output);
					/* nombre de loci, then SUM, AVERAGE and VARIANCE */
                    fprintf(file_output,"Total number of loci: %d\n\n",matrixml[0].nloci);
                    fprintf(file_output,"Total number of sites: %.2f\n",matrixml[0].Snsites);
                    fprintf(file_output,"Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                    fprintf(file_output,"Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                    fprintf(file_output,"Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance number of biallelic sites per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of theta (Watterson): %g\n",matrixml[0].Stheta_wat);
                    fprintf(file_output,"Average value of theta (Watterson) per loci: %g\n",(double)matrixml[0].Stheta_wat/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_wat,(double)matrixml[0].Stheta_wat,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Watterson) per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of theta (Tajima): %g\n",matrixml[0].Stheta_taj);
                    fprintf(file_output,"Average value of theta (Tajima) per loci: %g\n",(double)matrixml[0].Stheta_taj/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_taj,(double)matrixml[0].Stheta_taj,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Tajima) per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of theta (Fu & Li): %g\n",matrixml[0].Stheta_fulin);
                    fprintf(file_output,"Average value of theta (Fu & Li) per loci: %g\n",(double)matrixml[0].Stheta_fulin/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_fulin,(double)matrixml[0].Stheta_fulin,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Fu & Li) per loci: %g\n\n",var);
                    else fputs("\n",file_output);
                }
            }
            break;
         case '3':
            if(*outgroup) {
                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
                printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				printf(" NOTE that all theta estimates are corrected for chromosome population size.\n\n");
				/* nombre de loci, then SUM, AVERAGE and VARIANCE */
                printf("Total number of loci: %d\n\n",matrixml[0].nloci);
                printf("Total number of sites: %.2f\n",matrixml[0].Snsites);
                printf("Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                printf("Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                printf("Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance number of biallelic sites per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) {
                    printf("Total number of biallelic sites (incl. shared): %ld\n",matrixml[0].Sbiallsitesn);
                    printf("Average number of biallelic sites per loci (incl. shared): %g\n",(double)matrixml[0].Sbiallsitesn/(double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2biallsitesn,(double)matrixml[0].Sbiallsitesn,matrixml[0].nloci);
                    if(var > (double)-10000) printf("Variance number of biallelic sites per loci (incl. shared): %g\n\n",var);
                }
                printf("Total divergence per nucleotide: %g\n",matrixml[0].Sndivergence/(double)matrixml[0].Snsites);
                printf("Averaged divergence per loci and per nucleotide: %g\n",(double)matrixml[0].Sndivergence_nut/ (double)matrixml[0].nldivjc);
                var = variances((double)matrixml[0].S2ndivergence_nut,(double)matrixml[0].Sndivergence_nut,matrixml[0].nldivjc);
                if(var > (double)-10000) printf("Variance of divergence per loci and per nucleotide: %g\n\n",var);
                printf("Total divergence per nucleotide corrected by Jukes and Cantor: %g\n",matrixml[0].Sndivergencejc/ (double)matrixml[0].Snsites);
                printf("Averaged divergence per nucleotide per loci corrected by Jukes and Cantor: %g\n",(double)matrixml[0].Sndivergencejc_nut / (double)matrixml[0].nldivjc);
                var = variances((double)matrixml[0].S2ndivergencejc_nut,(double)matrixml[0].Sndivergencejc_nut,matrixml[0].nldivjc);
                if(var > (double)-10000) printf("Variance of divergence per nucleotide per loci corrected by Jukes and Cantor: %g\n\n",var);
                printf("Total value of theta (Watterson) per nucleotide: %g\n",matrixml[0].Stheta_wat/(double)matrixml[0].Snsites);
                printf("Average value of theta (Watterson) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_wat_nut/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_wat_nut,(double)matrixml[0].Stheta_wat_nut,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Watterson) per nucleotide and per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) {
                    printf("Total value of theta (Watterson) per nucleotide including shared pol.: %g\n",matrixml[0].Stheta_watn/(double)matrixml[0].Snsites);
                    printf("Average value of theta (Watterson) per nucleotide and per loci including shared pol.: %g\n",(double)matrixml[0].Stheta_watn_nut/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_watn_nut,(double)matrixml[0].Stheta_watn_nut,matrixml[0].nloci);
                    if(var > (double)-10000) printf("Variance of theta (Watterson) per nucleotide per loci including shared pol.: %g\n\n",var);
                }
                printf("Total value of theta (Tajima) per nucleotide: %g\n",matrixml[0].Stheta_taj/(double)matrixml[0].Snsites);
                printf("Average value of theta (Tajima) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_taj_nut/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_taj_nut,(double)matrixml[0].Stheta_taj_nut,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Tajima) per nucleotide and per loci: %g\n\n",var);
                if(matrixml[0].Sshared > 0) {
                    printf("Total value of theta (Tajima) per nucleotide including shared pol.: %g\n",matrixml[0].Stheta_tajn/(double)matrixml[0].Snsites);
                    printf("Average value of theta (Tajima) per nucleotide and per loci including shared pol.: %g\n",(double)matrixml[0].Stheta_tajn_nut/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_tajn_nut,(double)matrixml[0].Stheta_tajn_nut,matrixml[0].nloci);
                    if(var > (double)-10000) printf("Variance of theta (Tajima) per nucleotide and per loci including shared pol.: %g\n\n",var);
                }
                printf("Total value of theta (Fu & Li) per nucleotide: %g\n",matrixml[0].Stheta_fuli/(double)matrixml[0].Snsites);
                printf("Average value of theta (Fu & Li) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_fuli_nut/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_fuli_nut,(double)matrixml[0].Stheta_fuli_nut,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Fu & Li) per nucleotide and per loci: %g\n\n",var);
                printf("Total value of theta (Fay & Wu) per nucleotide: %g\n",matrixml[0].Stheta_fw/(double)matrixml[0].Snsites);
                printf("Average value of theta (Fay & Wu) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_fw_nut/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_fw_nut,(double)matrixml[0].Stheta_fw_nut,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Fay & Wu) per nucleotide and per loci: %g\n\n",var);
                printf("Total value of theta (Zeng et al.) per nucleotide: %g\n",matrixml[0].Stheta_L/(double)matrixml[0].Snsites);
                printf("Average value of theta (Zeng et al.) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_L_nut/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_L_nut,(double)matrixml[0].Stheta_L_nut,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Zeng et al.) per nucleotide and per loci: %g\n\n",var);
                else printf("\n");

                if(file_output) {
                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n\n",file_output);
					/* nombre de loci, then SUM, AVERAGE and VARIANCE */
                    fprintf(file_output,"Total number of loci: %d\n\n",matrixml[0].nloci);
                    fprintf(file_output,"Total number of sites: %.2f\n",matrixml[0].Snsites);
                    fprintf(file_output,"Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                    fprintf(file_output,"Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                    fprintf(file_output,"Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance number of biallelic sites per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0) {
                        fprintf(file_output,"Total number of biallelic sites (incl. shared): %ld\n",matrixml[0].Sbiallsitesn);
                        fprintf(file_output,"Average number of biallelic sites per loci (incl. shared): %g\n",(double)matrixml[0].Sbiallsitesn/(double)matrixml[0].nloci);
                        var = variances((double)matrixml[0].S2biallsitesn,(double)matrixml[0].Sbiallsitesn,matrixml[0].nloci);
                        if(var > (double)-10000) fprintf(file_output,"Variance number of biallelic sites per loci (incl. shared): %g\n\n",var);
                    }
                    fprintf(file_output,"Total divergence per nucleotide: %g\n",matrixml[0].Sndivergence/(double)matrixml[0].Snsites);
                    fprintf(file_output,"Averaged divergence per loci and per nucleotide: %g\n",(double)matrixml[0].Sndivergence_nut/ (double)matrixml[0].nldivjc);
                    var = variances((double)matrixml[0].S2ndivergence_nut,(double)matrixml[0].Sndivergence_nut,matrixml[0].nldivjc);
                    if(var > (double)-10000) fprintf(file_output,"Variance of divergence per loci and per nucleotide: %g\n\n",var);
                    fprintf(file_output,"Total divergence per nucleotide corrected by Jukes and Cantor: %g\n",matrixml[0].Sndivergencejc/(double)matrixml[0].Snsites);
                    fprintf(file_output,"Averaged divergence per nucleotide per loci corrected by Jukes and Cantor: %g\n",(double)matrixml[0].Sndivergencejc_nut/(double)matrixml[0].nldivjc);
                    var = variances((double)matrixml[0].S2ndivergencejc_nut,(double)matrixml[0].Sndivergencejc_nut,matrixml[0].nldivjc);
                    if(var > (double)-10000) fprintf(file_output,"Variance of divergence per nucleotide per loci corrected by Jukes and Cantor: %g\n\n",var);
                    fprintf(file_output,"Total value of theta (Watterson) per nucleotide: %g\n",matrixml[0].Stheta_wat/(double)matrixml[0].Snsites);
                    fprintf(file_output,"Average value of theta (Watterson) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_wat_nut/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_wat_nut,(double)matrixml[0].Stheta_wat_nut,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Watterson) per nucleotide and per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0) {
                        fprintf(file_output,"Total value of theta (Watterson) per nucleotide including shared pol.: %g\n",matrixml[0].Stheta_watn/(double)matrixml[0].Snsites);
                        fprintf(file_output,"Average value of theta (Watterson) per nucleotide and per loci including shared pol.: %g\n",(double)matrixml[0].Stheta_watn_nut/ (double)matrixml[0].nloci);
                        var = variances((double)matrixml[0].S2theta_watn_nut,(double)matrixml[0].Stheta_watn_nut,matrixml[0].nloci);
                        if(var > (double)-10000) fprintf(file_output,"Variance of theta (Watterson) per nucleotide and per loci including shared pol.: %g\n\n",var);
                    }
                    fprintf(file_output,"Total value of theta (Tajima) per nucleotide: %g\n",matrixml[0].Stheta_taj/(double)matrixml[0].Snsites);
                    fprintf(file_output,"Average value of theta (Tajima) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_taj_nut/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_taj_nut,(double)matrixml[0].Stheta_taj_nut,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Tajima) per nucleotide and per loci: %g\n\n",var);
                    if(matrixml[0].Sshared > 0) {
                        fprintf(file_output,"Total value of theta (Tajima) per nucleotide and including shared pol.: %g\n",matrixml[0].Stheta_tajn/(double)matrixml[0].Snsites);
                        fprintf(file_output,"Average value of theta (Tajima) per nucleotide and per loci including shared pol.: %g\n",(double)matrixml[0].Stheta_tajn_nut/ (double)matrixml[0].nloci);
                        var = variances((double)matrixml[0].S2theta_tajn_nut,(double)matrixml[0].Stheta_tajn_nut,matrixml[0].nloci);
                        if(var > (double)-10000) fprintf(file_output,"Variance of theta (Tajima) per nucleotide and per loci including shared pol.: %g\n\n",var);
                    }
                    fprintf(file_output,"Total value of theta (Fu & Li) per nucleotide: %g\n",matrixml[0].Stheta_fuli/(double)matrixml[0].Snsites);
                    fprintf(file_output,"Average value of theta (Fu & Li) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_fuli_nut/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_fuli_nut,(double)matrixml[0].Stheta_fuli_nut,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Fu & Li) per nucleotide and per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of theta (Fay & Wu) per nucleotide: %g\n",matrixml[0].Stheta_fw/(double)matrixml[0].Snsites);
                    fprintf(file_output,"Average value of theta (Fay & Wu) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_fw_nut/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_fw_nut,(double)matrixml[0].Stheta_fw_nut,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Fay & Wu) per nucleotide and per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of theta (Zeng et al.) per nucleotide: %g\n",matrixml[0].Stheta_L/(double)matrixml[0].Snsites);
                    fprintf(file_output,"Average value of theta (Zeng et al.) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_L_nut/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_L_nut,(double)matrixml[0].Stheta_L_nut,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Zeng et al.) per nucleotide and per loci: %g\n\n",var);
                    else fputs("\n",file_output);
                }
            }
            else {
                printf("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n");
				printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n");
                if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
				printf(" NOTE that all theta estimates are corrected for chromosome population size.\n\n");
				/* nombre de loci, then SUM, AVERAGE and VARIANCE */
                printf("Total number of loci: %d\n\n",matrixml[0].nloci);
                printf("Total number of sites: %.2f\n",matrixml[0].Snsites);
                printf("Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                printf("Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                printf("Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance number of biallelic sites per loci: %g\n\n",var);
                printf("Total value of theta (Watterson) per nucleotide: %g\n",matrixml[0].Stheta_wat/(double)matrixml[0].Snsites);
                printf("Average value of theta (Watterson) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_wat_nut/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_wat_nut,(double)matrixml[0].Stheta_wat_nut,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Watterson) per nucleotide and per loci: %g\n\n",var);
                printf("Total value of theta (Tajima) per nucleotide: %g\n",matrixml[0].Stheta_taj/(double)matrixml[0].Snsites);
                printf("Average value of theta (Tajima) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_taj_nut/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_taj_nut,(double)matrixml[0].Stheta_taj_nut,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Tajima) per nucleotide and per loci: %g\n\n",var);
                printf("Total value of theta (Fu & Li) per nucleotide: %g\n",matrixml[0].Stheta_fulin/(double)matrixml[0].Snsites);
                printf("Average value of theta (Fu & Li) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_fulin_nut/ (double)matrixml[0].nloci);
                var = variances((double)matrixml[0].S2theta_fulin_nut,(double)matrixml[0].Stheta_fulin_nut,matrixml[0].nloci);
                if(var > (double)-10000) printf("Variance of theta (Fu & Li) per nucleotide and per loci: %g\n\n",var);
                else printf("\n");
            
                if(file_output) {
                    fputs("\n\n Display multilocus statistics: \n\n Multiple hits not included in the analysis.\n",file_output);
                    fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n",file_output);
                    if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
					fputs(" NOTE that all theta estimates are corrected for chromosome population size.\n\n",file_output);
					/* nombre de loci, then SUM, AVERAGE and VARIANCE */
                    fprintf(file_output,"Total number of loci: %d\n\n",matrixml[0].nloci);
                    fprintf(file_output,"Total number of sites: %.2f\n",matrixml[0].Snsites);
                    fprintf(file_output,"Average number of sites per loci: %g\n\n",(double)matrixml[0].Snsites/(double)matrixml[0].nloci);
                    fprintf(file_output,"Total number of biallelic sites: %ld\n",matrixml[0].Sbiallsites);
                    fprintf(file_output,"Average number of biallelic sites per loci: %g\n",(double)matrixml[0].Sbiallsites/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance number of biallelic sites per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of theta (Watterson) per nucleotide: %g\n",matrixml[0].Stheta_wat/(double)matrixml[0].Snsites);
                    fprintf(file_output,"Average value of theta (Watterson) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_wat_nut/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_wat_nut,(double)matrixml[0].Stheta_wat_nut,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Watterson) per nucleotide and per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of theta (Tajima) per nucleotide: %g\n",matrixml[0].Stheta_taj/(double)matrixml[0].Snsites);
                    fprintf(file_output,"Average value of theta (Tajima) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_taj_nut/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_taj_nut,(double)matrixml[0].Stheta_taj_nut,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Tajima) per nucleotide and per loci: %g\n\n",var);
                    fprintf(file_output,"Total value of theta (Fu & Li) per nucleotide: %g\n",matrixml[0].Stheta_fulin/(double)matrixml[0].Snsites);
                    fprintf(file_output,"Average value of theta (Fu & Li) per nucleotide and per loci: %g\n",(double)matrixml[0].Stheta_fulin_nut/ (double)matrixml[0].nloci);
                    var = variances((double)matrixml[0].S2theta_fulin_nut,(double)matrixml[0].Stheta_fulin_nut,matrixml[0].nloci);
                    if(var > (double)-10000) fprintf(file_output,"Variance of theta (Fu & Li) per nucleotide and per loci: %g\n\n",var);
                    else fputs("\n",file_output);
                }
            }
            break;
        case '4':
            return;
            break;
    }
    return;
}

void open_displntmulo(struct statistics *matrix,struct statmulo *matrixml/*,int *observed_data*/,int *outgroup,int *n_loci,char *subset_positions,FILE *file_output)
{
    int x;
    double var;
    double variances(double,double,int);
    void print_prob_obshka(struct statmulo *, int *, FILE *);
       
	if(file_output) {		
		fprintf(file_output,"\n\n     MENU:\n     1. Observed analysis menu:");
		fprintf(file_output,"\n     1.1. Display neutrality tests menu:\n\n");		
		fprintf(file_output,"\n     1.1.0. Display MULTILOCUS Neutrality tests.\n");
	}

    if(*outgroup) {
        printf("\n\n Display multilocus neutral tests: \n\n Multiple hits not included in the analysis.\n");
        printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n\n");
        if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
		
        if(matrixml[0].nltajD > 0) {
            printf("Average multilocus Tajima's D test for %d loci: %g\n",matrixml[0].nltajD,matrixml[0].StajimaD/(double)matrixml[0].nltajD);
            if(matrixml[0].nltajD > 2) {
                var = variances((double)matrixml[0].S2tajimaD,(double)matrixml[0].StajimaD,matrixml[0].nltajD);
                printf("Variance multilocus Tajima's D test: %g\n\n",var);
            }
        }
        else printf("multilocus Tajima's D test not available.\n");        
        if(matrixml[0].nlflD > 0) {
            printf("Average multilocus Fu and Li's D test for %d loci: %g\n",matrixml[0].nlflD,matrixml[0].SfuliD/(double)matrixml[0].nlflD);
            if(matrixml[0].nlflD > 2) {
                var = variances((double)matrixml[0].S2fuliD,(double)matrixml[0].SfuliD,matrixml[0].nlflD);
                printf("Variance multilocus Fu and Li's D test: %g\n\n",var);
            }
        }
        else  printf("multilocus Fu and Li's D test not available.\n");
        if(matrixml[0].nlflF > 0) {
            printf("Average multilocus Fu and Li's F test for %d loci: %g\n",matrixml[0].nlflF,matrixml[0].SfuliF/(double)matrixml[0].nlflF);
            if(matrixml[0].nlflF > 2) {
                var = variances((double)matrixml[0].S2fuliF,(double)matrixml[0].SfuliF,matrixml[0].nlflF);
                printf("Variance multilocus Fu and Li's F test: %g\n\n",var);
            }
        }
        else  printf("multilocus Fu and Li's F test not available.\n");
        if(matrixml[0].nlHo > 0) {
            printf("Average multilocus Fay and Wu's H test for %d loci: %g\n",matrixml[0].nlHo,matrixml[0].SfaywuHo/(double)matrixml[0].nlHo);
            if(matrixml[0].nlHo > 2) {
                var = variances((double)matrixml[0].S2faywuHo,(double)matrixml[0].SfaywuHo,matrixml[0].nlHo);
                printf("Variance multilocus Fay and Wu's H test: %g\n\n",var);
            }
        }
        else  printf("multilocus Fay and Wu's H test not available.\n");
        if(matrixml[0].nlH > 0) {
            printf("Average multilocus normalized Fay and Wu's H test for %d loci: %g\n",matrixml[0].nlH,matrixml[0].SfaywuH/(double)matrixml[0].nlH);
            if(matrixml[0].nlH > 2) {
                var = variances((double)matrixml[0].S2faywuH,(double)matrixml[0].SfaywuH,matrixml[0].nlH);
                printf("Variance multilocus normalized Fay and Wu's H test: %g\n\n",var);
            }
        }
        else  printf("multilocus normalized Fay annd Wu's H test not available.\n");
        if(matrixml[0].nlFs > 0) {
            printf("Average multilocus Fu's Fs test for %d loci: %g\n",matrixml[0].nlFs,matrixml[0].SfuFs/(double)matrixml[0].nlFs);
            if(matrixml[0].nlFs > 2) {
                var = variances((double)matrixml[0].S2fuFs,(double)matrixml[0].SfuFs,matrixml[0].nlFs);
                printf("Variance multilocus Fu's Fs test: %g\n\n",var);
            }
        }
        else  printf("multilocus Fu's Fs test not available.\n");
        if(matrixml[0].nlZ > 0) {
            printf("Average multilocus Rozas's et al. ZA test for %d loci: %g\n",matrixml[0].nlZ,matrixml[0].SrZA/(double)matrixml[0].nlZ);
            if(matrixml[0].nlZ > 2) {
                var = variances((double)matrixml[0].S2rZA,(double)matrixml[0].SrZA,matrixml[0].nlZ);
                printf("Variance multilocus Rozas's et al. ZA test: %g\n\n",var);
            }
        }
        else  printf("multilocus Rozas' et al. ZA test not available.\n");
        if(matrixml[0].nlB > 0) {
            printf("Average multilocus Wall's B test for %d loci: %g\n",matrixml[0].nlB,matrixml[0].SwB/(double)matrixml[0].nlB);
            if(matrixml[0].nlB > 2) {
                var = variances((double)matrixml[0].S2wB,(double)matrixml[0].SwB,matrixml[0].nlB);
                printf("Variance multilocus Wall's B test: %g\n\n",var);
            }
        }
        else  printf("multilocus Wall's B test not available.\n");
        if(matrixml[0].nlQ > 0) {
            printf("Average multilocus Wall's Q test for %d loci: %g\n",matrixml[0].nlQ,matrixml[0].SwQ/(double)matrixml[0].nlQ);
            if(matrixml[0].nlQ > 2) {
                var = variances((double)matrixml[0].S2wQ,(double)matrixml[0].SwQ,matrixml[0].nlQ);
                printf("Variance multilocus Wall's Q test: %g\n\n",var);
            }
        }
        else  printf("multilocus Wall's Q test not available.\n");
        
        if(matrixml[0].nlR2 > 0) {
            printf("Average multilocus R2 test for %d loci: %g\n",
                matrixml[0].nlR2,matrixml[0].SR2/(double)matrixml[0].nlR2);
            if(matrixml[0].nlR2 > 2) {
                var = variances((double)matrixml[0].S2R2,(double)matrixml[0].SR2,matrixml[0].nlR2);
                printf("Variance multilocus R2 test: %g\n\n",var);
            }
        }
        else  printf("multilocus R2 test not available.\n");
        if(matrixml[0].nlE > 0) {
            printf("Average multilocus Zeng et al. E test for %d loci: %g\n",matrixml[0].nlE,matrixml[0].SzengE/(double)matrixml[0].nlE);
            if(matrixml[0].nlE > 2) {
                var = variances((double)matrixml[0].S2zengE,(double)matrixml[0].SzengE,matrixml[0].nlE);
                printf("Variance multilocus Zeng et al. E test: %g\n\n",var);
            }
        }
        else  printf("multilocus Zeng et al.E test not available.\n");

		printf("Average multilocus Ewens-Watterson test for %d loci: %g\n",matrixml[0].nloci,(double)matrixml[0].Sewtest/(double)matrixml[0].nloci);
		var = variances((double)matrixml[0].S2ewtest,(double)matrixml[0].Sewtest,matrixml[0].nloci);
		if(var > (double)-10000) printf("Variance multilocus Ewens-Watterson test: %g\n\n",var);

        if(matrixml[0].Sshared > 0)
            printf("\n\nNot shared polymorphisms included in the above neutrality tests.\n\n");

        /*HKA test.*/
        if(*n_loci > 1 && matrix[0].hka >= (double)0.0) {
            printf("HKA test X (corrected by Jukes and Cantor): %g\n",matrixml[0].Shka);
            print_prob_obshka(matrixml,n_loci,0);
            /*parameters.*/
            printf("Estimated parameters using the HKA model:\n");
            printf("\n\tT(in 2N)\t%g\n",matrixml[0].hka_T);
            if(matrixml[0].hka_T <= 0) printf("Warning!! A negative T indicates no difference between ingroup and outgroup!.\n");
			printf("\nValues given per locus and corrected for chromosome population size: value of expected theta corrected by Jukes and Cantor (considering mhits are excluded), ");
			printf("\nnumber of observed segregating sites and observed divergence, same but corrected by Jukes and Cantor, ");
			printf("\nand expected number of segregating sites and divergence corrected by Jukes and Cantor.\n Also the partial value of Chi-square for each locus.\n");
			printf("\nname_locus\ttheta_exp(JC)\tSobs\tDiv_obs\tSobs(JC)\tDiv_obs(JC)\tSexp(JC)\tDiv_exp(JC)\tpartial_hka\n\n");
			for(x = 0; x < *n_loci; x++) {
				printf("%d:%-20s",x,matrix[x].gene);
                printf("\t%g",matrix[x].hka_theta*((double)1/matrix[x].factor_chrn));
                printf("\t%d",matrix[x].biallsitesn);
                printf("\t%g",matrix[x].ndivergence);
                printf("\t%g",matrix[x].Sobshka);
                printf("\t%g",matrix[x].Dobshka);
                printf("\t%g",matrix[x].Sexphka);
                printf("\t%g",matrix[x].Dexphka);
                printf("\t%g",matrix[x].hka);
                printf("\n");
            }
			printf("\n");
        }
        else  fputs("HKA test not available.\n",file_output);
        
        if(file_output) {
            fputs("\n\n Display multilocus neutral tests: \n\n Multiple hits not included in the analysis.\n",file_output);
            fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n\n",file_output);
            if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
			
            if(matrixml[0].nltajD > 0) {
                fprintf(file_output,"Average multilocus Tajima's D test for %d loci: %g\n",matrixml[0].nltajD,matrixml[0].StajimaD/(double)matrixml[0].nltajD);
                if(matrixml[0].nltajD > 2) {
                    var = variances((double)matrixml[0].S2tajimaD,(double)matrixml[0].StajimaD,matrixml[0].nltajD);
                    fprintf(file_output,"Variance multilocus Tajima's D test: %g\n\n",var);
                }
            }
            else fputs("multilocus Tajima's D test not available.\n",file_output);        
            if(matrixml[0].nlflD > 0) {
                fprintf(file_output,"Average multilocus Fu and Li's D test for %d loci: %g\n",matrixml[0].nlflD,matrixml[0].SfuliD/(double)matrixml[0].nlflD);
                if(matrixml[0].nlflD > 2) {
                    var = variances((double)matrixml[0].S2fuliD,(double)matrixml[0].SfuliD,matrixml[0].nlflD);
                    fprintf(file_output,"Variance multilocus Fu and Li's D test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Fu and Li's D test not available.\n",file_output);
            if(matrixml[0].nlflF > 0) {
                fprintf(file_output,"Average multilocus Fu and Li's F test for %d loci: %g\n",matrixml[0].nlflF,matrixml[0].SfuliF/(double)matrixml[0].nlflF);
                if(matrixml[0].nlflF > 2) {
                    var = variances((double)matrixml[0].S2fuliF,(double)matrixml[0].SfuliF,matrixml[0].nlflF);
                    fprintf(file_output,"Variance multilocus Fu and Li's F test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Fu and Li's F test not available.\n",file_output);
            if(matrixml[0].nlHo > 0) {
                fprintf(file_output,"Average multilocus Fay and Wu's H test for %d loci: %g\n",matrixml[0].nlHo,matrixml[0].SfaywuHo/(double)matrixml[0].nlHo);
                if(matrixml[0].nlHo > 2) {
                    var = variances((double)matrixml[0].S2faywuHo,(double)matrixml[0].SfaywuHo,matrixml[0].nlHo);
                    fprintf(file_output,"Variance multilocus Fay and Wu's H test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Fay annd Wu's H test not available.\n",file_output);
            if(matrixml[0].nlH > 0) {
                fprintf(file_output,"Average multilocus normalized Fay and Wu's H test for %d loci: %g\n",matrixml[0].nlH,matrixml[0].SfaywuH/(double)matrixml[0].nlH);
                if(matrixml[0].nlH > 2) {
                    var = variances((double)matrixml[0].S2faywuH,(double)matrixml[0].SfaywuH,matrixml[0].nlH);
                    fprintf(file_output,"Variance multilocus normalized Fay and Wu's H test: %g\n\n",var);
                }
            }
            else  fputs("multilocus normalized Fay annd Wu's H test not available.\n",file_output);
            if(matrixml[0].nlFs > 0) {
                fprintf(file_output,"Average multilocus Fu's Fs test for %d loci: %g\n",matrixml[0].nlFs,matrixml[0].SfuFs/(double)matrixml[0].nlFs);
                if(matrixml[0].nlFs > 2) {
                    var = variances((double)matrixml[0].S2fuFs,(double)matrixml[0].SfuFs,matrixml[0].nlFs);
                    fprintf(file_output,"Variance multilocus Fu's Fs test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Fu's Fs test not available.\n",file_output);
            if(matrixml[0].nlZ > 0) {
                fprintf(file_output,"Average multilocus Rozas's et al. ZA test for %d loci: %g\n",matrixml[0].nlZ,matrixml[0].SrZA/(double)matrixml[0].nlZ);
                if(matrixml[0].nlZ > 2) {
                    var = variances((double)matrixml[0].S2rZA,(double)matrixml[0].SrZA,matrixml[0].nlZ);
                    fprintf(file_output,"Variance multilocus Rozas's et al. ZA test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Rozas' et al. ZA test not available.\n",file_output);
            if(matrixml[0].nlB > 0) {
                fprintf(file_output,"Average multilocus Wall's B test for %d loci: %g\n",matrixml[0].nlB,matrixml[0].SwB/(double)matrixml[0].nlB);
                if(matrixml[0].nlB > 2) {
                    var = variances((double)matrixml[0].S2wB,(double)matrixml[0].SwB,matrixml[0].nlB);
                    fprintf(file_output,"Variance multilocus Wall's B test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Wall's B test not available.\n",file_output);
            if(matrixml[0].nlQ > 0) {
                fprintf(file_output,"Average multilocus Wall's Q test for %d loci: %g\n",matrixml[0].nlQ,matrixml[0].SwQ/(double)matrixml[0].nlQ);
                if(matrixml[0].nlQ > 2) {
                    var = variances((double)matrixml[0].S2wQ,(double)matrixml[0].SwQ,matrixml[0].nlQ);
                    fprintf(file_output,"Variance multilocus Wall's Q test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Wall's Q test not available.\n",file_output);

            if(matrixml[0].nlR2 > 0) {
                fprintf(file_output,"Average multilocus R2 test for %d loci: %g\n",
                    matrixml[0].nlR2,matrixml[0].SR2/(double)matrixml[0].nlR2);
                if(matrixml[0].nlR2 > 2) {
                    var = variances((double)matrixml[0].S2R2,(double)matrixml[0].SR2,matrixml[0].nlR2);
                    fprintf(file_output,"Variance multilocus R2 test: %g\n\n",var);
                }
            }
            else  fputs("multilocus R2 test not available.\n",file_output);
            if(matrixml[0].nlE > 0) {
                fprintf(file_output,"Average multilocus Zeng et al. E test for %d loci: %g\n",matrixml[0].nlE,matrixml[0].SzengE/(double)matrixml[0].nlE);
                if(matrixml[0].nlH > 2) {
                    var = variances((double)matrixml[0].S2zengE,(double)matrixml[0].SzengE,matrixml[0].nlE);
                    fprintf(file_output,"Variance multilocus Zeng et al. E test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Zeng et al. E test not available.\n",file_output);

			fprintf(file_output,"Average multilocus Ewens-Watterson test for %d loci: %g\n",matrixml[0].nloci,(double)matrixml[0].Sewtest/(double)matrixml[0].nloci);
			var = variances((double)matrixml[0].S2ewtest,(double)matrixml[0].Sewtest,matrixml[0].nloci);
			if(var > (double)-10000) fprintf(file_output,"Variance multilocus Ewens-Watterson test: %g\n\n",var);

            if(matrixml[0].Sshared > 0)
                printf("\n\nNot shared polymorphisms included in the above neutrality tests.\n\n");

            if(*n_loci > 1 && matrix[0].hka >= (double)0.0) {
                fprintf(file_output,"HKA test (values corrected by Jukes and Cantor): %g\n",matrixml[0].Shka);
                print_prob_obshka(matrixml,n_loci,file_output);
                fputs("Estimated parameters using the HKA model:\n",file_output);
                fprintf(file_output,"\n\tT(in 2N)\t%g\n\n",matrixml[0].hka_T);
				if(matrixml[0].hka_T <= 0) fprintf(file_output,"Warning!! A negative T indicates no difference between ingroup and outgroup!.\n");
				fprintf(file_output,"\nValues given per locus and corrected for chromosome population size: value of expected theta corrected by Jukes and Cantor (considering mhits are excluded), ");
				fprintf(file_output,"\nnumber of observed segregating sites and observed divergence, same but corrected by Jukes and Cantor, ");
				fprintf(file_output,"\nand expected number of segregating sites and divergence corrected by Jukes and Cantor.\n Also the partial value of Chi-square for each locus.\n");
				fprintf(file_output,"\nname_locus\ttheta_exp(JC)\tSobs\tDiv_obs\tSobs(JC)\tDiv_obs(JC)\tSexp(JC)\tDiv_exp(JC)\tpartial_hka\n\n");
				for(x = 0; x < *n_loci; x++) {
					fprintf(file_output,"%d:%-20s",x,matrix[x].gene);
					fprintf(file_output,"\t%g",matrix[x].hka_theta*((double)1/matrix[x].factor_chrn));
					fprintf(file_output,"\t%d",matrix[x].biallsitesn);
					fprintf(file_output,"\t%g",matrix[x].ndivergence);
					fprintf(file_output,"\t%g",matrix[x].Sobshka);
					fprintf(file_output,"\t%g",matrix[x].Dobshka);
					fprintf(file_output,"\t%g",matrix[x].Sexphka);
					fprintf(file_output,"\t%g",matrix[x].Dexphka);
					fprintf(file_output,"\t%g",matrix[x].hka);
					fprintf(file_output,"\n");
				}
				fprintf(file_output,"\n");
            }
            else  fputs("HKA test not available.\n",file_output);
        }
    }
    else {
        printf("\n\n Display multilocus neutral tests: \n\n Multiple hits not included in the analysis.\n");
        printf(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n\n");
        if(subset_positions[0]) printf(" Subset of positions analyzed: %s.\n",subset_positions);
		
        if(matrixml[0].nltajD > 0) {
            printf("Average multilocus Tajima's D test for %d loci: %g\n",matrixml[0].nltajD,matrixml[0].StajimaD/(double)matrixml[0].nltajD);
            if(matrixml[0].nltajD > 2) {
                var = variances((double)matrixml[0].S2tajimaD,(double)matrixml[0].StajimaD,matrixml[0].nltajD);
                printf("Variance multilocus Tajima's D test: %g\n\n",var);
            }
        }
        else printf("multilocus Tajima's D test not available.\n");        
        
		if(matrixml[0].nlflDn > 0) {
            printf("Average multilocus Fu and Li's D* test for %d loci: %g\n",matrixml[0].nlflDn,matrixml[0].SfuliDn/(double)matrixml[0].nlflDn);
            if(matrixml[0].nlflDn > 2) {
                var = variances((double)matrixml[0].S2fuliDn,(double)matrixml[0].SfuliDn,matrixml[0].nlflDn);
                printf("Variance multilocus Fu and Li's D* test: %g\n\n",var);
            }
        }
        else  printf("multilocus Fu and Li's D* test not available.\n");
        
		if(matrixml[0].nlflFn > 0) {
            printf("Average multilocus Fu and Li's F* test for %d loci: %g\n",matrixml[0].nlflFn,matrixml[0].SfuliFn/(double)matrixml[0].nlflFn);
            if(matrixml[0].nlflFn > 2) {
                var = variances((double)matrixml[0].S2fuliFn,(double)matrixml[0].SfuliFn,matrixml[0].nlflFn);
                printf("Variance multilocus Fu and Li's F* test: %g\n\n",var);
            }
        }
        else  printf("multilocus Fu and Li's F* test not available.\n");
        
		if(matrixml[0].nlFs > 0) {
            printf("Average multilocus Fu's Fs test for %d loci: %g\n",matrixml[0].nlFs,matrixml[0].SfuFs/(double)matrixml[0].nlFs);
            if(matrixml[0].nlFs > 2) {
                var = variances((double)matrixml[0].S2fuFs,(double)matrixml[0].SfuFs,matrixml[0].nlFs);
                printf("Variance multilocus Fu's Fs test: %g\n\n",var);
            }
        }
        else  printf("multilocus Fu's Fs test not available.\n");
        
		if(matrixml[0].nlZ > 0) {
            printf("Average multilocus Rozas's et al. ZA test for %d loci: %g\n",matrixml[0].nlZ,matrixml[0].SrZA/(double)matrixml[0].nlZ);
            if(matrixml[0].nlZ > 2) {
                var = variances((double)matrixml[0].S2rZA,(double)matrixml[0].SrZA,matrixml[0].nlZ);
                printf("Variance multilocus Rozas's et al. ZA test: %g\n\n",var);
            }
        }
        else  printf("multilocus Rozas' et al. ZA test not available.\n");
        
		if(matrixml[0].nlB > 0) {
            printf("Average multilocus Wall's B test for %d loci: %g\n",matrixml[0].nlB,matrixml[0].SwB/(double)matrixml[0].nlB);
            if(matrixml[0].nlB > 2) {
                var = variances((double)matrixml[0].S2wB,(double)matrixml[0].SwB,matrixml[0].nlB);
                printf("Variance multilocus Wall's B test: %g\n\n",var);
            }
        }
        else  printf("multilocus Wall's B test not available.\n");
        
		if(matrixml[0].nlQ > 0) {
            printf("Average multilocus Wall's Q test for %d loci: %g\n",matrixml[0].nlQ,matrixml[0].SwQ/(double)matrixml[0].nlQ);
            if(matrixml[0].nlQ > 2) {
                var = variances((double)matrixml[0].S2wQ,(double)matrixml[0].SwQ,matrixml[0].nlQ);
                printf("Variance multilocus Wall's Q test: %g\n\n",var);
            }
        }
        else  printf("multilocus Wall's Q test not available.\n");
        
        if(matrixml[0].nlR2 > 0) {
            printf("Average multilocus R2 test for %d loci: %g\n",
                matrixml[0].nlR2,matrixml[0].SR2/(double)matrixml[0].nlR2);
            if(matrixml[0].nlR2 > 2) {
                var = variances((double)matrixml[0].S2R2,(double)matrixml[0].SR2,matrixml[0].nlR2);
                printf("Variance multilocus R2 test: %g\n\n",var);
            }
        }
        else  printf("multilocus R2 test not available.\n");

		printf("Average multilocus Ewens-Watterson test for %d loci: %g\n",matrixml[0].nloci,(double)matrixml[0].Sewtest/(double)matrixml[0].nloci);
		var = variances((double)matrixml[0].S2ewtest,(double)matrixml[0].Sewtest,matrixml[0].nloci);
		if(var > (double)-10000) printf("Variance multilocus Ewens-Watterson test: %g\n\n",var);

        if(file_output) {
            fputs("\n\n Display multilocus neutral tests: \n\n Multiple hits not included in the analysis.\n",file_output);
            fputs(" Insertion-deletions are not considered in the ENTIRE alignment (including outgroup lines).\n\n",file_output);
            if(subset_positions[0]) fprintf(file_output," Subset of positions analyzed: %s.\n",subset_positions);
			
            if(matrixml[0].nltajD > 0) {
                fprintf(file_output,"Average multilocus Tajima's D test for %d loci: %g\n",matrixml[0].nltajD,matrixml[0].StajimaD/(double)matrixml[0].nltajD);
                if(matrixml[0].nltajD > 2) {
                    var = variances((double)matrixml[0].S2tajimaD,(double)matrixml[0].StajimaD,matrixml[0].nltajD);
                    fprintf(file_output,"Variance multilocus Tajima's D test: %g\n\n",var);
                }
            }
            else fputs("multilocus Tajima's D test not available.\n",file_output);        
            
			if(matrixml[0].nlflDn > 0) {
                fprintf(file_output,"Average multilocus Fu and Li's D* test for %d loci: %g\n",matrixml[0].nlflDn,matrixml[0].SfuliDn/(double)matrixml[0].nlflDn);
                if(matrixml[0].nlflDn > 2) {
                    var = variances((double)matrixml[0].S2fuliDn,(double)matrixml[0].SfuliDn,matrixml[0].nlflDn);
                    fprintf(file_output,"Variance multilocus Fu and Li's D* test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Fu and Li's D* test not available.\n",file_output);
            
			if(matrixml[0].nlflFn > 0) {
                fprintf(file_output,"Average multilocus Fu and Li's F* test for %d loci: %g\n",matrixml[0].nlflFn,matrixml[0].SfuliFn/(double)matrixml[0].nlflFn);
                if(matrixml[0].nlflFn > 2) {
                    var = variances((double)matrixml[0].S2fuliFn,(double)matrixml[0].SfuliFn,matrixml[0].nlflFn);
                    fprintf(file_output,"Variance multilocus Fu and Li's F* test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Fu and Li's F* test not available.\n",file_output);
            
			if(matrixml[0].nlFs > 0) {
                fprintf(file_output,"Average multilocus Fu's Fs test for %d loci: %g\n",matrixml[0].nlFs,matrixml[0].SfuFs/(double)matrixml[0].nlFs);
                if(matrixml[0].nlFs > 2) {
                    var = variances((double)matrixml[0].S2fuFs,(double)matrixml[0].SfuFs,matrixml[0].nlFs);
                    fprintf(file_output,"Variance multilocus Fu's Fs test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Fu's Fs test not available.\n",file_output);
            
			if(matrixml[0].nlZ > 0) {
                fprintf(file_output,"Average multilocus Rozas's et al. ZA test for %d loci: %g\n",matrixml[0].nlZ,matrixml[0].SrZA/(double)matrixml[0].nlZ);
                if(matrixml[0].nlZ > 2) {
                    var = variances((double)matrixml[0].S2rZA,(double)matrixml[0].SrZA,matrixml[0].nlZ);
                    fprintf(file_output,"Variance multilocus Rozas's et al. ZA test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Rozas' et al. ZA test not available.\n",file_output);
            
			if(matrixml[0].nlB > 0) {
                fprintf(file_output,"Average multilocus Wall's B test for %d loci: %g\n",matrixml[0].nlB,matrixml[0].SwB/(double)matrixml[0].nlB);
                if(matrixml[0].nlB > 2) {
                    var = variances((double)matrixml[0].S2wB,(double)matrixml[0].SwB,matrixml[0].nlB);
                    fprintf(file_output,"Variance multilocus Wall's B test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Wall's B test not available.\n",file_output);
            
			if(matrixml[0].nlQ > 0) {
                fprintf(file_output,"Average multilocus Wall's Q test for %d loci: %g\n",matrixml[0].nlQ,matrixml[0].SwQ/(double)matrixml[0].nlQ);
                if(matrixml[0].nlQ > 2) {
                    var = variances((double)matrixml[0].S2wQ,(double)matrixml[0].SwQ,matrixml[0].nlQ);
                    fprintf(file_output,"Variance multilocus Wall's Q test: %g\n\n",var);
                }
            }
            else  fputs("multilocus Wall's Q test not available.\n",file_output);

            if(matrixml[0].nlR2 > 0) {
                fprintf(file_output,"Average multilocus R2 test for %d loci: %g\n",
                    matrixml[0].nlR2,matrixml[0].SR2/(double)matrixml[0].nlR2);
                if(matrixml[0].nlR2 > 2) {
                    var = variances((double)matrixml[0].S2R2,(double)matrixml[0].SR2,matrixml[0].nlR2);
                    fprintf(file_output,"Variance multilocus R2 test: %g\n\n",var);
                }
            }
            else  fputs("multilocus R2 test not available.\n",file_output);

			fprintf(file_output,"Average multilocus Ewens-Watterson test for %d loci: %g\n",matrixml[0].nloci,(double)matrixml[0].Sewtest/(double)matrixml[0].nloci);
			var = variances((double)matrixml[0].S2ewtest,(double)matrixml[0].Sewtest,matrixml[0].nloci);
			if(var > (double)-10000) fprintf(file_output,"Variance multilocus Ewens-Watterson test: %g\n\n",var);
        }
    }    
    return;
}

void open_displayhistobs(struct statistics *matrix/*,int *observed_data*/,int *outgroup,int *n_loci/*,char *subset_positions*/,FILE *file_output,int test)
{
    char k[1];
    int h,x,y,z;
    long int ss;
    double *vector_obs;
    char *name_obs;
    int start_histogram(double *,double *,int,long int,char *,char *,FILE *);
    
	if(file_output) {
		if(test == 1) {
			fprintf(file_output,"\n\n     MENU:\n     1. Observed analysis menu:");
			fprintf(file_output,"\n     1.1. Display neutrality tests menu:");
			fprintf(file_output,"\n     1.1.1. Display detailed neutrality tests menu:\n\n");
			fprintf(file_output,"\n     1.1.1.1  Display histograms with tests of neutrality.\n");
		}
		else {
            fprintf(file_output,"\n\n     MENU:\n     1. Observed analysis menu:");
            fprintf(file_output,"\n     1.0. Display statistics menu:");
            fprintf(file_output,"\n     1.0.1. Display detailed statistics menu:\n\n");
            fprintf(file_output,"\n     1.0.1.1. Display histograms with the statistics.\n");

		}
	}
	
    /*We need calculate a vector with the data and also give the name of the data we are using*/
    if((vector_obs = (double *) malloc((*n_loci)*sizeof(double))) == 0) {
        printf("\nError: memory not reallocated. open_displayhistobs.1 \n");
        if(file_output) fputs("\nError: memory not reallocated. Histograms not displayed. \n",file_output);
        return;
    }

    if((name_obs = (char *) malloc(128*sizeof(char))) == 0) {
        printf("\nError: memory not reallocated. open_displayhistobs.1 \n");
        if(file_output) fputs("\nError: memory not reallocated. Histograms not displayed. \n",file_output);
        free(vector_obs);
		return;
    }
        
    if(test) {
        for(z=0;z<13;z++) {
            if(!(*outgroup)) if(z == 3) z = 4;
            if(!(*outgroup)) if(z == 10) break;
			if(*n_loci == 1) if(z == 10) break;
			if(z == 12) if(*n_loci <= 1 || matrix[0].hka < (double)0.00) break;
            switch(z) {
                case 0:
                    /*Tajima's D*/
                    y = 0;
                    for(x=0;x<*n_loci;x++) {
                        if(matrix[x].tajimaD != (double) -10000) {
                            vector_obs[y] = matrix[x].tajimaD;
                            y++;
                        }
                    }
                    strcpy(name_obs,"Obs. Tajima's D");
                    break;
                case 1:
                    /*FuLi's D*/
                    y = 0;
                    if(*outgroup) {
                        for(x=0;x<*n_loci;x++) {
                            if(matrix[x].fuliD != (double) -10000) {
                                vector_obs[y] = matrix[x].fuliD;
                                y++;
                            }
                        }
                        strcpy(name_obs,"Obs. Fu and Li's D");
                    }
                    else {
                        for(x=0;x<*n_loci;x++) {
                            if(matrix[x].fuliDn != (double) -10000) {
                                vector_obs[y] = matrix[x].fuliDn;
                                y++;
                            }
                        }
                        strcpy(name_obs,"Obs. Fu and Li's D*");
                    }
                    break;
                case 2:
                    /*FuLi's F*/
                    y = 0;
                    if(*outgroup) {
                        for(x=0;x<*n_loci;x++) {
                            if(matrix[x].fuliF != (double) -10000) {
                                vector_obs[y] = matrix[x].fuliF;
                                y++;
                            }
                        }
                        strcpy(name_obs,"Obs. Fu and Li's F");
                    }
                    else {
                        for(x=0;x<*n_loci;x++) {
                            if(matrix[x].fuliFn != (double) -10000) {
                                vector_obs[y] = matrix[x].fuliFn;
                                y++;
                            }
                        }
                        strcpy(name_obs,"Obs. Fu and Li's F*");
                    }
                    break;
                case 3:
                    /*FayWu H*/
                    y = 0;
                    for(x=0;x<*n_loci;x++) {
                        if(matrix[x].faywuH != (double) -10000) {
                            vector_obs[y] = matrix[x].faywuH;
                            y++;
                        }
                    }
                    strcpy(name_obs,"Obs. normalized Fay and Wu's H");
                    break;
                case 4:
                    /*Fu Fs*/
                    y = 0;
                    for(x=0;x<*n_loci;x++) {
                        if(matrix[x].fuFs != (double) -10000) {
                            vector_obs[y] = matrix[x].fuFs;
                            y++;
                        }
                    }
                    strcpy(name_obs,"Obs. Fu's Fs");
                    break;
                case 5:
                    /*Rozas ZA*/
                    y = 0;
                    for(x=0;x<*n_loci;x++) {
                        if(matrix[x].rZA != (double) -10000) {
                            vector_obs[y] = matrix[x].rZA;
                            y++;
                        }
                    }
                    strcpy(name_obs,"Obs. Rozas' et al. ZA");
                    break;
                case 6:
                    /*Wall's B*/
                    y = 0;
                    for(x=0;x<*n_loci;x++) {
                        if(matrix[x].wB != (double) -10000) {
                            vector_obs[y] = matrix[x].wB;
                            y++;
                        }
                    }
                    strcpy(name_obs,"Obs. Wall's B");
                    break;
                case 7:
                    /*Wall's Q*/
                    y = 0;
                    for(x=0;x<*n_loci;x++) {
                        if(matrix[x].wQ != (double) -10000) {
                            vector_obs[y] = matrix[x].wQ;
                            y++;
                        }
                    }
                    strcpy(name_obs,"Obs. Wall's Q");
                    break;
                case 8:
                    /*R2*/
                    y = 0;
                    for(x=0;x<*n_loci;x++) {
                        if(matrix[x].R2 != (double) -10000) {
                            vector_obs[y] = matrix[x].R2;
                            y++;
                        }
                    }
                    strcpy(name_obs,"Obs. Ramos-Onsins & Rozas' R2");
                    break;
                case 9:
                    /*EWtest*/
                    y = 0;
                    for(x=0;x<*n_loci;x++) {
                        if(matrix[x].ewtest != (double) -10000) {
                            vector_obs[y] = matrix[x].ewtest;
                            y++;
                        }
                    }
                    strcpy(name_obs,"Obs. Ewens-Watterson test");
                    break;
                case 10:
                    /*Zeng E*/
                    y = 0;
                    for(x=0;x<*n_loci;x++) {
                        if(matrix[x].zengE != (double) -10000) {
                            vector_obs[y] = matrix[x].zengE;
                            y++;
                        }
                    }
                    strcpy(name_obs,"Obs. Zeng et al. E");
                    break;
                case 11:
                    /*FW*/
                    y = 0;
                    for(x=0;x<*n_loci;x++) {
                        if(matrix[x].faywuHo != (double) -10000) {
                            vector_obs[y] = matrix[x].faywuHo;
                            y++;
                        }
                    }
                    strcpy(name_obs,"Obs. Fay and Wu's H");
                    break;
                case 12:
                    /*Partial HKA tests*/
                    y = 0;
                    for(x=0;x<*n_loci;x++) {
                        if(matrix[x].hka >= 0) {
                            vector_obs[y] = matrix[x].hka;
                            y++;
                        }
                    }
                    strcpy(name_obs,"Obs. Chi-sq HKA(JC)/locus");
                    break;
            }
            if((h = start_histogram(vector_obs,0,y,0,name_obs,0,file_output)) == 0) {
                printf("\nHistogram for %s is not available.\n",name_obs);
                if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
            }
            else if(h == 2) z = 13;
            strcpy(name_obs,"");
        }
    }
    else {
        #if COMMAND_LINE
        
        printf("\n\nMENU 1. Observed analysis menu:");
        printf("\n     1.0. Display statistics menu:");
        printf("\n     1.0.1. Display detailed statistics menu:");
        printf("\n     1.0.1.1 Display histograms with detailed statistics for each locus:\n\n");
        
        printf("CHOOSE:\n");
        printf(" 0 - Display general statistics.\n");
        printf(" 1 - Display other linkage related statistics.\n");
        printf(" 2 - Display estimates of total locus variability.\n");
        printf(" 3 - Display estimates of nucleotide locus variability.\n");
        printf(" 4 - Back to Display detailed neutrality tests menu.\n");

        if(file_output) {
            
            fprintf(file_output,"\n\n     MENU:\n     1. Observed analysis menu:");
            fprintf(file_output,"\n     1.0. Display statistics menu:");
            fprintf(file_output,"\n     1.0.1. Display detailed statistics menu:");
            fprintf(file_output,"\n     1.0.1.1 Display histograms with detailed statistics for each locus:\n\n");
            
            fprintf(file_output," 0 - Display general statistics.\n");
            fprintf(file_output," 1 - Display other linkage related statistics.\n");
            fprintf(file_output," 2 - Display estimates of total locus variability.\n");
            fprintf(file_output," 3 - Display estimates of nucleotide locus variability.\n");
            fprintf(file_output," 4 - Back to Display detailed neutrality tests menu.\n");
        }
        #endif
        
        do *k = getchar();
        while(*k<'0' || *k>'4');
        
        ss = 0;
        for(y=0;y<*n_loci;y++) ss += matrix[y].shared;

        switch(*k) {
            case '0':
                for(z=0;z<10;z++) {
                    if(*outgroup) if((z == 1 || z == 2) && ss == 0) z = 3;
                    if(!(*outgroup)) if(z == 1) z = 5;
                    switch(z) {
                        case 0:
                            /*biallsites*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].biallsites;
                            strcpy(name_obs,"Biallelic sites");
                            break;
                        case 1:
                            /*biallsitesn*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].biallsitesn;
                            strcpy(name_obs,"Biallelic + shared");
                            break;
                        case 2:
                            /*shared*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].shared;
                            strcpy(name_obs,"Shared pol");
                            break;
                        case 3:
                            /*fixed*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].fixed;
                            strcpy(name_obs,"Fixed sites");
                            break;
                        case 4:
                            /*ndivergence*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].ndivergence;
                            strcpy(name_obs,"Divergence");
                            break;
                        case 5:
                            /*nhapl*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].nhapl;
                            strcpy(name_obs,"Haplotypes");
                            break;
                        case 6:
                            /*nhapl/sam*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].nhaplsam;
                            strcpy(name_obs,"Haplotypes/sample size");
                            break;
                        case 7:
                            /*hapldiv*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].hapldiv;
                            strcpy(name_obs,"Haplotype diversity");
                            break;
                        case 8:
                            /*Rvpi*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].Rvpi;
                            strcpy(name_obs,"Recombination size");
                            break;
                        case 9:
                            /*Rm*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].Rm;
                            strcpy(name_obs,"Minimum Recombination events");
                            break;
                    }
                    if((h = start_histogram(vector_obs,0,x,0,name_obs,0,file_output)) == 0) {
                        printf("\nHistogram for %s is not available.\n",name_obs);
                        if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                    }
                    else if(h == 2) z = 9;
                    strcpy(name_obs,"");
                }
                break;
            case '1':
                for(z=0;z<3;z++) {
                    switch(z) {
                        case 0:
                            /*za*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].za;
                            strcpy(name_obs,"Rozas' za (not div. by S)");
                            break;
                        case 1:
                            /*b*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].b;
                            strcpy(name_obs,"Wall's b (not div. by S)");
                            break;
                        case 2:
                            /*q*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].q;
                            strcpy(name_obs,"Wall's q (not div. by S)");
                            break;
                    }
                    if((h = start_histogram(vector_obs,0,x,0,name_obs,0,file_output)) == 0) {
                        printf("\nHistogram for %s is not available.\n",name_obs);
                        if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                    }
                    else if(h == 2) z = 3;
                    strcpy(name_obs,"");
                }
                break;
            case '2':
                for(z=0;z<9;z++) {
                    if(*outgroup) if((z == 1) && ss == 0) z = 2;
                    if(*outgroup) if((z == 3) && ss == 0) z = 4;
                    if(*outgroup) if(z == 6) z = 7;
                    if(!(*outgroup)) if(z == 1) z = 2;
                    if(!(*outgroup)) if(z == 3) z = 4;
                    if(!(*outgroup)) if(z == 4) z = 6;
                    if(!(*outgroup)) if(z == 7) break;
                    switch(z) {
                        case 0:
                            /*theta_wat*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_wat*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta (Watterson) corrected for chromosome population size");
                            break;
                        case 1:
                            /*theta_watn*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_watn*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta (Watterson) incl. shared corrected for chromosome population size");
                            break;
                        case 2:
                            /*theta_taj*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_taj*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta (Tajima) corrected for chromosome population size");
                            break;
                        case 3:
                            /*theta_tajn*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_tajn*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta (Tajima) incl. shared corrected for chromosome population size");
                            break;
                        case 4:
                            /*theta_fw*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_fw*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta (Fay and Wu) corrected for chromosome population size");
                            break;
                        case 5:
                            /*theta_fuli*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta (Fu and Li) corrected for chromosome population size");
                            break;
                        case 6:
                            /*theta_fulin*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta (Fu and Li) corrected for chromosome population size");
                            break;
                        case 7:
                            /*theta_L*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_L*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta (Zeng et al.) corrected for chromosome population size");
                            break;
                        case 8:
                            /*ndivergence*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].ndivergence;
                            strcpy(name_obs,"Divergence");
                            break;
                    }
                    if((h = start_histogram(vector_obs,0,x,0,name_obs,0,file_output)) == 0) {
                        printf("\nHistogram for %s is not available.\n",name_obs);
                        if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                    }
                    else if(h == 2) z = 9;
                    strcpy(name_obs,"");
                }
                break;
            case '3':
                for(z=0;z<9;z++) {
                    if(*outgroup) if((z == 1) && ss == 0) z = 2;
                    if(*outgroup) if((z == 3) && ss == 0) z = 4;
                    if(*outgroup) if(z == 6) z = 7;
                    if(!(*outgroup)) if(z == 1) z = 2;
                    if(!(*outgroup)) if(z == 3) z = 4;
                    if(!(*outgroup)) if(z == 4) z = 6;
                    if(!(*outgroup)) if(z == 7) break;
                    switch(z) {
                        case 0:
                            /*theta_wat*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_wat/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta/nt (Watterson) corrected for chromosome population size");
                            break;
                        case 1:
                            /*theta_watn*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_watn/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta/nt (Watterson) incl. shared corrected for chromosome population size");
                            break;
                        case 2:
                            /*theta_taj*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_taj/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta/nt (Tajima) corrected for chromosome population size");
                            break;
                        case 3:
                            /*theta_tajn*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_tajn/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta/nt (Tajima) incl. shared corrected for chromosome population size");
                            break;
                        case 4:
                            /*theta_fw*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_fw/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta/nt (Fay and Wu) corrected for chromosome population size");
                            break;
                        case 5:
                            /*theta_fuli*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_fuli/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta/nt (Fu and Li) corrected for chromosome population size");
                            break;
                        case 6:
                            /*theta_fulin*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_fulin/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta/nt (Fu and Li) corrected for chromosome population size");
                            break;
                        case 7:
                            /*theta_L*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].theta_L/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
                            strcpy(name_obs,"Theta/nt (Zeng et al.) corrected for chromosome population size");
                            break;
                        case 8:
                            /*ndivergence*/
                            for(x=0;x<*n_loci;x++) vector_obs[x] = matrix[x].ndivergence/(double)matrix[x].nsites;
                            strcpy(name_obs,"Divergence/nt");
                            break;
                    }
                    if((h = start_histogram(vector_obs,0,x,0,name_obs,0,file_output)) == 0) {
                        printf("\nHistogram for %s is not available.\n",name_obs);
                        if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                    }
                    else if(h == 2) z = 9;
                    strcpy(name_obs,"");
                }
                break;
            case '4':
				free(name_obs);
				free(vector_obs);
                return;
                break;
        }
    }    
    free(name_obs);
    free(vector_obs);
    return;
}

