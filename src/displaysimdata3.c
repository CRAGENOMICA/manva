/*
 *  displaysimdata3.c
 *  MuLoNeTests
 *
 *  Created by sebas on Tue Mar 18 2003.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void open_simdisplayhist(struct var *data,struct var2b *inputms,struct statistics *matrix,struct statmulo *matrixml,struct statistisim **matrixsim,struct statistisimmuloc *matrixmlsim/*,struct horizontalstatsml *avgstatloci,int *observed_data*/,int *outgroup,int *n_loci/*,int *montecarlo_sim,int onlymulo*/, FILE *file_output,int dataobsequalsim,int neuttest,int mulo)
{
    /*Display histogram*/
    char k[1];
    int h,x,y,z,v;
    long int xs,ys,zs;
    double *vector_obs=0;
	double *vector_sim=0;
    char name_obs[250],name_sim[250];
    int start_histogram(double *,double *,int,long int,char *,char *,FILE *);
    double varianceit(double,double,long int);
    double variances(double,double,int);

    
    name_obs[0]='\0';
    name_sim[0]='\0';
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
                    printf("\n     3.1.0.1.2. Display histograms with detailed statistics for each locus:\n\n");
                    
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
                        fprintf(file_output,"\n     3.1.0.1.2. Display histograms with detailed statistics for each locus:\n\n");
                        
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

                    if(*k < '4') {
                        vector_obs = 0;
                        name_obs[0] = '\0';
                        if(dataobsequalsim)
                            if((vector_obs = (double *) malloc((*n_loci)*sizeof(double))) == 0) {
                                printf("\nError: memory not reallocated. Histograms not displayed. \n");
                                if(file_output) fputs("\nError: memory not reallocated. Histograms not displayed. \n",file_output);
                                return;
                            }
                        if((vector_sim = (double *) malloc(((data[0].n_iter)*(data[0].n_loci))*sizeof(double))) == 0) {
                            printf("\nError: memory not reallocated. Histograms not displayed. \n");
                            if(file_output) fputs("\nError: memory not reallocated. Histograms not displayed. \n",file_output);
                            return;
                        }
                    }

                    switch(*k) {
                        case '0':
                            for(z=0;z<7;z++) {
                                if(!(data[0].time_spec > 0.)) if(z == 1) z = 3;
                                x = 0;
                                switch(z) {
                                    case 0:
                                        /*biallsites*/
                                        if(dataobsequalsim) 
                                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].biallsites;
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].biallsites;

                                        if(dataobsequalsim) strcpy(name_obs,"Obs. biallelic sites (not shared)");
                                        strcpy(name_sim,"Exp. biallelic sites (not shared)");
                                        break;
                                    case 1:
                                        /*fixed*/
                                        if(dataobsequalsim) 
                                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].fixed;
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].fixed;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. fixed sites");
                                        strcpy(name_sim,"Exp. fixed sites");
                                        break;
                                    case 2:
                                        /*ndivergence*/
                                        if(dataobsequalsim) 
                                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].ndivergence;
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].ndivergence;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. divergence");
                                        strcpy(name_sim,"Exp. divergence");
                                        break;
                                    case 3:
                                        /*nhapl*/
                                        if(dataobsequalsim) 
                                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].nhapl;
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].nhapl;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. number of haplotypes");
                                        strcpy(name_sim,"Exp. number of haplotypes");
                                        break;
                                    case 4:
                                        /*nhaplsam*/
                                        if(dataobsequalsim) 
                                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].nhaplsam;
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].nhaplsam;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. haplotypes / sample size");
                                        strcpy(name_sim,"Exp. haplotypes / sample size");
                                        break;
                                    case 5:
                                        /*hapldiv*/
                                        if(dataobsequalsim) 
                                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].hapldiv;
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].hapldiv;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. haplotype diversity");
                                        strcpy(name_sim,"Exp. haplotype diversity");
                                        break;
 									case 6:
                                        /*Rm*/
                                        if(dataobsequalsim) 
                                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].Rm;
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].Rm;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Rm (Hudson and Kaplan 1985)");
                                        strcpy(name_sim,"Exp. Rm (Hudson and Kaplan 1985)");
                                        break;
                                }
                                if((h = start_histogram(vector_obs,vector_sim,x,xs,name_obs,name_sim,file_output)) == 0) {
                                    printf("\nHistogram for %s is not available.\n",name_obs);
                                    if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                                }
                                else if(h == 2) z = 7;
                                strcpy(name_obs,"");
                                strcpy(name_sim,"");
                            }
                            free(vector_obs);
                            free(vector_sim);
                            break;
                        case '1':
                            for(z=0;z<3;z++) {
                                x = 0;
                                switch(z) {
                                    case 0:
                                        /*za*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].za;
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].za;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Rozas' za (not divided by S)");
                                        strcpy(name_sim,"Exp. Rozas' za (not divided by S)");
                                        break;
                                    case 1:
                                        /*b*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].b;
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].b;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Wall's b (not divided by S)");
                                        strcpy(name_sim,"Exp. Wall's b (not divided by S)");
                                        break;
                                    case 2:
                                        /*q*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].q;
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].q;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Wall's q (not divided by S)");
                                        strcpy(name_sim,"Exp. Wall's q (not divided by S)");
                                        break;
                                }
                                if((h = start_histogram(vector_obs,vector_sim,x,xs,name_obs,name_sim,file_output)) == 0) {
                                    printf("\nHistogram for %s is not available.\n",name_obs);
                                    if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                                }
                                else if(h == 2) z = 3;
                                strcpy(name_obs,"");
                                strcpy(name_sim,"");
                            }
                            free(vector_obs);
                            free(vector_sim);
                            break;
                        case '2':
                            for(z=0;z<7;z++) {
                                if(data[0].time_spec > 0.) if(z == 5) z = 6;
                                if(!(data[0].time_spec > 0.)) if(z == 2) z = 5;
                                if(!(data[0].time_spec > 0.)) if(z == 6) break;
                                x = 0;
                                switch(z) {
                                    case 0:
                                        /*theta_wat*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].theta_wat*((double)1/matrix[x].factor_chrn);
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].theta_wat*((double)1/inputms[v].factor_chrn);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Theta (Watterson) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. Theta (Watterson) corrected for chromosome population size");
                                        break;
                                    case 1:
                                        /*theta_taj*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].theta_taj*((double)1/matrix[x].factor_chrn);
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].theta_taj*((double)1/inputms[v].factor_chrn);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Theta (Tajima) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. Theta (Tajima) corrected for chromosome population size");
                                        break;
                                    case 2:
                                        /*theta_fw*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].theta_fw*((double)1/matrix[x].factor_chrn);
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].theta_fw*((double)1/inputms[v].factor_chrn);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Theta (normalized Fay and Wu) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. Theta (normalized Fay and Wu) corrected for chromosome population size");
                                        break;
                                    case 3:
                                        /*theta_fuli*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn);
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].theta_fuli*((double)1/inputms[v].factor_chrn);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Theta (Fu and Li) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. Theta (Fu and Li) corrected for chromosome population size");
                                        break;
                                    case 4:
                                        /*theta_L*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].theta_L*((double)1/matrix[x].factor_chrn);
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].theta_L*((double)1/inputms[v].factor_chrn);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Theta (Zeng et al.) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. Theta (Zeng et al.) corrected for chromosome population size");
                                        break;
                                    case 5:
                                        /*theta_fulin*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn);
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].theta_fulin*((double)1/inputms[v].factor_chrn);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Theta (Fu and Li) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. Theta (Fu and Li) corrected for chromosome population size");
                                        break;
                                    case 6:
                                        /*ndivergence*/
                                        if(dataobsequalsim) 
                                            for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].ndivergence;
                                        xs = 0;
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].ndivergence;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. divergence");
                                        strcpy(name_sim,"Exp. divergence");
                                        break;
                                }
                                if((h = start_histogram(vector_obs,vector_sim,x,xs,name_obs,name_sim,file_output)) == 0) {
                                    printf("\nHistogram for %s is not available.\n",name_obs);
                                    if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                                }
                                else if(h == 2) z = 7;
                                strcpy(name_obs,"");
                                strcpy(name_sim,"");
                            }
                            free(vector_obs);
                            free(vector_sim);
                            break;
                        case '3':
                            for(z=0;z<7;z++) {
                                if(data[0].time_spec > 0.) if(z == 5) z = 6;
                                if(!(data[0].time_spec > 0.)) if(z == 2) z = 5;
                                if(!(data[0].time_spec > 0.)) if(z == 6) break;
                                x = 0;
                                xs = 0;
                                switch(z) {
                                    case 0:
                                        /*theta_wat*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].theta_wat/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn);
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].theta_wat/(double)(matrixsim[v][zs].nsites)*((double)1/inputms[v].factor_chrn);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Theta/nt (Watterson) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. Theta/nt (Watterson) corrected for chromosome population size");
                                        break;
                                    case 1:
                                        /*theta_taj*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].theta_taj/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn);
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].theta_taj/(double)(matrixsim[v][zs].nsites)*((double)1/inputms[v].factor_chrn);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Theta/nt (Tajima) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. Theta/nt (Tajima) corrected for chromosome population size");
                                        break;
                                    case 2:
                                        /*theta_fw*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].theta_fw/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn);
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].theta_fw/(double)(matrixsim[v][zs].nsites)*((double)1/inputms[v].factor_chrn);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Theta/nt (normalized Fay and Wu) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. Theta/nt (normalized Fay and Wu) corrected for chromosome population size");
                                        break;
                                    case 3:
                                        /*theta_fuli*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].theta_fuli/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn);
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].theta_fuli/(double)(matrixsim[v][zs].nsites)*((double)1/inputms[v].factor_chrn);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Theta/nt (Fu and Li) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. Theta/nt (Fu and Li) corrected for chromosome population size");
                                        break;
                                    case 4:
                                        /*theta_L*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].theta_L/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn);
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].theta_L/(double)(matrixsim[v][zs].nsites)*((double)1/inputms[v].factor_chrn);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Theta/nt (Zeng et al.) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. Theta/nt (Zeng et al.) corrected for chromosome population size");
                                        break;
                                    case 5:
                                        /*theta_fulin*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].theta_fulin/(double)(matrix[x].nsites-matrix[x].shared)*((double)1/matrix[x].factor_chrn);
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].theta_fulin/(double)(matrixsim[v][zs].nsites)*((double)1/inputms[v].factor_chrn);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. Theta/nt (Fu and Li) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. Theta/nt (Fu and Li) corrected for chromosome population size");
                                        break;
                                    case 6:
                                        /*ndivergence*/
                                        if(dataobsequalsim) for(x=0;x<*n_loci;x++) vector_obs[x] = (double)matrix[x].ndivergence/(double)(matrix[x].nsites);
                                        for(zs=0;zs<(long int)data[0].n_iter;zs++)
                                            for(v=0;v<data[0].n_loci;v++)
                                                vector_sim[xs++] = (double)matrixsim[v][zs].ndivergence/(double)(matrixsim[v][zs].nsites);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. divergence/nt");
                                        strcpy(name_sim,"Exp. divergence/nt");
                                        break;
                                }
                                if((h = start_histogram(vector_obs,vector_sim,x,xs,name_obs,name_sim,file_output)) == 0) {
                                    printf("\nHistogram for %s is not available.\n",name_obs);
                                    if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                                }
                                else if(h == 2) z = 7;
                                strcpy(name_obs,"");
                                strcpy(name_sim,"");
                            }
                            free(vector_obs);
                            free(vector_sim);
                            break;
                        case '4':
                            return;
                            break;
                    }
                    break;
                case 1: /* neutrality tests, all loci, puf..*/
					if(file_output) {    
						fprintf(file_output,"\n\n     MENU:\n     3. Coalescent Monte Carlo simulations and Statistical inference:");
						fprintf(file_output,"\n     3.1. Display coalescent simulation analysis menu:");
						fprintf(file_output,"\n     3.1.1. Display neutrality tests menu:");
						fprintf(file_output,"\n     3.1.1.1. Display detailed neutrality tests menu:\n\n");
						fprintf(file_output,"\n     3.1.1.1.2. Display histograms of simulated results.\n");
					}

                    vector_obs = 0;
                    name_sim[0]='\0';
                    if(dataobsequalsim) 
                        if((vector_obs = (double *) malloc((*n_loci)*sizeof(double))) == 0) {
                            printf("\nError: memory not reallocated. Histograms not displayed. \n");
                            if(file_output) fputs("\nError: memory not reallocated. Histograms not displayed. \n",file_output);
                            return;
                        }
                    if((vector_sim = (double *) malloc(((data[0].n_iter)*(data[0].n_loci))*sizeof(double))) == 0) {
                        printf("\nError: memory not reallocated. Histograms not displayed. \n");
                        if(file_output) fputs("\nError: memory not reallocated. Histograms not displayed. \n",file_output);
                        return;
                    }
                    for(z=0;z<13;z++) {
                        if(!(data[0].time_spec > 0.)) if(z == 3) z = 5;
                        if(!(data[0].time_spec > 0.)) if(z == 11) break;
						if(data[0].n_loci == 1 || matrix[0].hka < (double)0.000) if(z == 12) break;
                        switch(z) {
                            case 0:
                                /*Tajima's D*/
                                y = 0;
                                if(dataobsequalsim) {
                                    for(x=0;x<*n_loci;x++) {
                                        if(matrix[x].tajimaD != (double) -10000) {
                                            vector_obs[y] = matrix[x].tajimaD;
                                            y++;
                                        }
                                    }
                                    strcpy(name_obs,"Obs. Tajima's D");
                                }
                                ys = 0;
                                for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                    for(v=0;v<data[0].n_loci;v++) {
                                        if(matrixsim[v][zs].tajimaD != (double) -10000) {
                                            vector_sim[ys] = matrixsim[v][zs].tajimaD;
                                            ys++;
                                        }
                                    }
                                }
                                strcpy(name_sim,"Exp. Tajima's D");
                                break;
                            case 1:
                                /*FuLi's D*/
                                y = 0;
                                if(data[0].time_spec > 0.) {
                                    if(dataobsequalsim) 
                                        for(x=0;x<*n_loci;x++) {
                                            if(matrix[x].fuliD != (double) -10000) {
                                                vector_obs[y] = matrix[x].fuliD;
                                                y++;
                                            }
                                        }
                                    ys = 0;
                                    for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                        for(v=0;v<data[0].n_loci;v++) {
                                            if(matrixsim[v][zs].fuliD != (double) -10000) {
                                                vector_sim[ys] = matrixsim[v][zs].fuliD;
                                                ys++;
                                            }
                                        }
                                    }
                                    if(dataobsequalsim) strcpy(name_obs,"Obs. Fu and Li's D");
                                    strcpy(name_sim,"Exp. Fu and Li's D");
                                }
                                else {
                                    if(dataobsequalsim) 
                                        for(x=0;x<*n_loci;x++) {
                                            if(matrix[x].fuliDn != (double) -10000) {
                                                vector_obs[y] = matrix[x].fuliDn;
                                                y++;
                                            }
                                        }
                                    ys = 0;
                                    for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                        for(v=0;v<data[0].n_loci;v++) {
                                            if(matrixsim[v][zs].fuliDn != (double) -10000) {
                                                vector_sim[ys] = matrixsim[v][zs].fuliDn;
                                                ys++;
                                            }
                                        }
                                    }
                                    if(dataobsequalsim) strcpy(name_obs,"Obs. Fu and Li's D*");
                                    strcpy(name_sim,"Exp. Fu and Li's D*");
                                }
                                break;
                            case 2:
                                /*FuLi's F*/
                                y = 0;
                                if(data[0].time_spec > 0.) {
                                    if(dataobsequalsim) 
                                        for(x=0;x<*n_loci;x++) {
                                            if(matrix[x].fuliF != (double) -10000) {
                                                vector_obs[y] = matrix[x].fuliF;
                                                y++;
                                            }
                                        }
                                    ys = 0;
                                    for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                        for(v=0;v<data[0].n_loci;v++) {
                                            if(matrixsim[v][zs].fuliF != (double) -10000) {
                                                vector_sim[ys] = matrixsim[v][zs].fuliF;
                                                ys++;
                                            }
                                        }
                                    }
                                    if(dataobsequalsim) strcpy(name_obs,"Obs. Fu and Li's F");
                                        strcpy(name_sim,"Exp. Fu and Li's F");
                                }
                                else {
                                    if(dataobsequalsim) 
                                        for(x=0;x<*n_loci;x++) {
                                            if(matrix[x].fuliFn != (double) -10000) {
                                                vector_obs[y] = matrix[x].fuliFn;
                                                y++;
                                            }
                                        }
                                    ys = 0;
                                    for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                        for(v=0;v<data[0].n_loci;v++) {
                                            if(matrixsim[v][zs].fuliFn != (double) -10000) {
                                                vector_sim[ys] = matrixsim[v][zs].fuliFn;
                                                ys++;
                                            }
                                        }
                                    }
                                    if(dataobsequalsim) strcpy(name_obs,"Obs. Fu and Li's F*");
                                        strcpy(name_sim,"Exp. Fu and Li's F*");
                                    }
                                break;
                            case 3:
                                /*FayWu Hn*/
                                y = 0;
                                if(dataobsequalsim) {
                                    for(x=0;x<*n_loci;x++) {
                                        if(matrix[x].faywuH != (double) -10000) {
                                            vector_obs[y] = matrix[x].faywuH;
                                            y++;
                                        }
                                    }
                                    strcpy(name_obs,"Obs. normalized Fay and Wu's H");
                                }
                                ys = 0;
                                for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                    for(v=0;v<data[0].n_loci;v++) {
                                        if(matrixsim[v][zs].faywuH != (double) -10000) {
                                            vector_sim[ys] = matrixsim[v][zs].faywuH;
                                            ys++;
                                        }
                                    }
                                }
                                strcpy(name_sim,"Exp. normalized Fay and Wu's H");
                                break;
                            case 4:
                                /*FayWu H*/
                                y = 0;
                                if(dataobsequalsim) {
                                    for(x=0;x<*n_loci;x++) {
                                        if(matrix[x].faywuHo != (double) -10000) {
                                            vector_obs[y] = matrix[x].faywuHo;
                                            y++;
                                        }
                                    }
                                    strcpy(name_obs,"Obs. Fay and Wu's H");
                                }
                                ys = 0;
                                for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                    for(v=0;v<data[0].n_loci;v++) {
                                        if(matrixsim[v][zs].faywuHo != (double) -10000) {
                                            vector_sim[ys] = matrixsim[v][zs].faywuHo;
                                            ys++;
                                        }
                                    }
                                }
                                strcpy(name_sim,"Exp. Fay and Wu's H");
                                break;
                            case 5:
                                /*Fu Fs*/
                                y = 0;
                                if(dataobsequalsim) {
                                    for(x=0;x<*n_loci;x++) {
                                        if(matrix[x].fuFs != (double) -10000) {
                                            vector_obs[y] = matrix[x].fuFs;
                                            y++;
                                        }
                                    }
                                    strcpy(name_obs,"Obs. Fu's Fs");
                                }
                                ys = 0;
                                for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                    for(v=0;v<data[0].n_loci;v++) {
                                        if(matrixsim[v][zs].fuFs != (double) -10000) {
                                            vector_sim[ys] = matrixsim[v][zs].fuFs;
                                            ys++;
                                        }
                                    }
                                }
                                strcpy(name_sim,"Exp. Fu's Fs");
                                break;
                            case 6:
                                /*Rozas ZA*/
                                y = 0;
                                if(dataobsequalsim) {
                                    for(x=0;x<*n_loci;x++) {
                                        if(matrix[x].rZA != (double) -10000) {
                                            vector_obs[y] = matrix[x].rZA;
                                            y++;
                                        }
                                    }
                                    strcpy(name_obs,"Obs. Rozas' et al. ZA");
                                }
                                ys = 0;
                                for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                    for(v=0;v<data[0].n_loci;v++) {
                                        if(matrixsim[v][zs].rZA != (double) -10000) {
                                            vector_sim[ys] = matrixsim[v][zs].rZA;
                                            ys++;
                                        }
                                    }
                                }
                                strcpy(name_sim,"Exp. Roza's et al. ZA");
                                break;
                            case 7:
                                /*Wall's B*/
                                y = 0;
                                if(dataobsequalsim) {
                                    for(x=0;x<*n_loci;x++) {
                                        if(matrix[x].wB != (double) -10000) {
                                            vector_obs[y] = matrix[x].wB;
                                            y++;
                                        }
                                    }
                                    strcpy(name_obs,"Obs. Wall's B");
                                }
                                ys = 0;
                                for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                    for(v=0;v<data[0].n_loci;v++) {
                                        if(matrixsim[v][zs].wB != (double) -10000) {
                                            vector_sim[ys] = matrixsim[v][zs].wB;
                                            ys++;
                                        }
                                    }
                                }
                                strcpy(name_sim,"Exp. Wall's B");
                                break;
                            case 8:
                                /*Wall's Q*/
                                y = 0;
                                if(dataobsequalsim) {
                                    for(x=0;x<*n_loci;x++) {
                                        if(matrix[x].wQ != (double) -10000) {
                                            vector_obs[y] = matrix[x].wQ;
                                            y++;
                                        }
                                    }
                                    strcpy(name_obs,"Obs. Wall's Q");
                                }
                                ys = 0;
                                for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                    for(v=0;v<data[0].n_loci;v++) {
                                        if(matrixsim[v][zs].wQ != (double) -10000) {
                                            vector_sim[ys] = matrixsim[v][zs].wQ;
                                            ys++;
                                        }
                                    }
                                }
                               strcpy(name_sim,"Exp. Wall's Q");
                                break;
                            case 9:
                                /*R2*/
                                y = 0;
                                if(dataobsequalsim) {
                                    for(x=0;x<*n_loci;x++) {
                                        if(matrix[x].R2 != (double) -10000) {
                                            vector_obs[y] = matrix[x].R2;
                                            y++;
                                        }
                                    }
                                    strcpy(name_obs,"Obs. Ramos-Onsins & Rozas' R2");
                                }
                                ys = 0;
                                for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                    for(v=0;v<data[0].n_loci;v++) {
                                        if(matrixsim[v][zs].R2 != (double) -10000) {
                                            vector_sim[ys] = matrixsim[v][zs].R2;
                                            ys++;
                                        }
                                    }
                                }
                                strcpy(name_sim,"Exp. Ramos-Onsins & Rozas' R2");
                                break;
                            case 10:
                                /*EW*/
                                y = 0;
                                if(dataobsequalsim) {
                                    for(x=0;x<*n_loci;x++) {
                                        if(matrix[x].ewtest != (double) -10000) {
                                            vector_obs[y] = matrix[x].ewtest;
                                            y++;
                                        }
                                    }
                                    strcpy(name_obs,"Obs. Ewens-Watterson test");
                                }
                                ys = 0;
                                for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                    for(v=0;v<data[0].n_loci;v++) {
                                        if(matrixsim[v][zs].ewtest != (double) -10000) {
                                            vector_sim[ys] = matrixsim[v][zs].ewtest;
                                            ys++;
                                        }
                                    }
                                }
                                strcpy(name_sim,"Exp. Ewens-Watterson test");
                                break;
                            case 11:
                                /*Zeng E*/
                                y = 0;
                                if(dataobsequalsim) {
                                    for(x=0;x<*n_loci;x++) {
                                        if(matrix[x].zengE != (double) -10000) {
                                            vector_obs[y] = matrix[x].zengE;
                                            y++;
                                        }
                                    }
                                    strcpy(name_obs,"Obs. Zeng et al. E test");
                                }
                                ys = 0;
                                for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                    for(v=0;v<data[0].n_loci;v++) {
                                        if(matrixsim[v][zs].zengE != (double) -10000) {
                                            vector_sim[ys] = matrixsim[v][zs].zengE;
                                            ys++;
                                        }
                                    }
                                }
                                strcpy(name_sim,"Exp. Zeng et al. E test");
                                break;
                            case 12:
                                /*Partial HKA tests*/
                                y = 0;
                                if(dataobsequalsim) {
                                    for(x=0;x<*n_loci;x++) {
                                        if(matrix[x].hka >= 0) {
                                            vector_obs[y] = matrix[x].hka;
                                            y++;
                                        }
                                    }
                                    strcpy(name_obs,"Obs. Chi-sq HKA(JC)/locus");
                                }
                                ys = 0;
                                for(zs=0;zs<(long int)data[0].n_iter;zs++) {
                                    for(v=0;v<data[0].n_loci;v++) {
                                        if(matrixsim[v][zs].hka != (double) -10000) {
                                            vector_sim[ys] = matrixsim[v][zs].hka;
                                            ys++;
                                        }
                                    }
                                }
                                strcpy(name_sim,"Exp. Chi-sq HKA(JC)/locus");
                                break;
                        }
                        if((h = start_histogram(vector_obs,vector_sim,y,ys,name_obs,name_sim,file_output)) == 0) {
                            printf("\nHistogram for %s is not available.\n",name_obs);
                            if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                        }
                        else if(h == 2) z = 13;
                        strcpy(name_obs,"");
                        strcpy(name_sim,"");
                    }
                    free(vector_obs);
                    free(vector_sim);
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
                    printf("\n     3.1.0.0.2. Display histograms with detailed statistics for multilocus data:\n\n");
                    
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
                        fprintf(file_output,"\n     3.1.0.0.2. Display histograms with detailed statistics for multilocus data:\n\n");
                        
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

                    if(*k < '4') {
                        vector_obs = 0;
                        name_obs[0] = '\0';
                        if(dataobsequalsim) 
                            if((vector_obs = (double *) malloc(sizeof(double))) == 0) {
                                printf("\nError: memory not reallocated. Histograms not displayed. \n");
                                if(file_output) fputs("\nError: memory not reallocated. Histograms not displayed. \n",file_output);
                                return;
                            }
                        if((vector_sim = (double *) malloc((data[0].n_iter)*sizeof(double))) == 0) {
                            printf("\nError: memory not reallocated. Histograms not displayed. \n");
                            if(file_output) fputs("\nError: memory not reallocated. Histograms not displayed. \n",file_output);
                            return;
                        }
                    }

                    switch(*k) {
                        case '0':
                            for(z=0;z<14;z++) {
                                if(z == 7 && data[0].n_loci < 3) break;
                                if(!(data[0].time_spec > 0.)) if(z == 1) z = 3;
                                if(z > 6 && data[0].n_loci > 2) {
                                    if(!(data[0].time_spec > 0.)) if(z == 8) z = 10;
                                }
                                x = 0;
                                switch(z) {
                                    case 0:
                                        /*average biallsites*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Sbiallsites/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Sbiallsites/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average biallelic sites (not shared)");
                                        strcpy(name_sim,"Exp. average biallelic sites (not shared)");
                                        break;
                                    case 7:
                                        /*variance biallsites*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances((double)matrixml[0].S2biallsites,(double)matrixml[0].Sbiallsites,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit((double)matrixmlsim[xs].S2biallsites,(double)matrixmlsim[xs].Sbiallsites,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance biallelic sites (not shared)");
                                        strcpy(name_sim,"Exp. variance biallelic sites (not shared)");
                                        break;
                                    case 1:
                                        /*average fixed*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Sfixed/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Sfixed/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average fixed sites");
                                        strcpy(name_sim,"Exp. average fixed sites");
                                        break;
                                    case 8:
                                        /*variance fixed*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances((double)matrixml[0].S2fixed,(double)matrixml[0].Sfixed,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit((double)matrixmlsim[xs].S2fixed,(double)matrixmlsim[xs].Sfixed,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance fixed sites");
                                        strcpy(name_sim,"Exp. variance fixed sites");
                                        break;
                                    case 2:
                                        /*average ndivergence*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Sndivergence/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Sndivergence/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average divergence");
                                        strcpy(name_sim,"Exp. average divergence");
                                        break;
                                    case 9:
                                        /*variance ndivergence*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2ndivergence,matrixml[0].Sndivergence,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2ndivergence,matrixmlsim[xs].Sndivergence,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance divergence");
                                        strcpy(name_sim,"Exp. variance divergence");
                                        break;
                                    case 3:
                                        /*average nhapl*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Snhapl/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Snhapl/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average number of haplotypes");
                                        strcpy(name_sim,"Exp. average number of haplotypes");
                                        break;
                                    case 10:
                                        /*variance nhapl*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances((double)matrixml[0].S2nhapl,(double)matrixml[0].Snhapl,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit((double)matrixmlsim[xs].S2nhapl,(double)matrixmlsim[xs].Snhapl,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance number of haplotypes");
                                        strcpy(name_sim,"Exp. variance number of haplotypes");
                                        break;
                                    case 4:
                                        /*average nhaplsam*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Snhaplsam/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Snhaplsam/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average number of haplotypes/nsam");
                                        strcpy(name_sim,"Exp. average number of haplotypes/nsam");
                                        break;
                                    case 11:
                                        /*variance nhaplsam*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances((double)matrixml[0].S2nhaplsam,(double)matrixml[0].Snhaplsam,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit((double)matrixmlsim[xs].S2nhaplsam,(double)matrixmlsim[xs].Snhaplsam,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance number of haplotypes/nsam");
                                        strcpy(name_sim,"Exp. variance number of haplotypes/nsam");
                                        break;
                                    case 5:
                                        /*average hapldiv*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Shapldiv/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Shapldiv/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average of haplotype diversity");
                                        strcpy(name_sim,"Exp. average of haplotype diversity");
                                        break;
                                    case 12:
                                        /*variance hapldiv*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances((double)matrixml[0].S2hapldiv,(double)matrixml[0].Shapldiv,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit((double)matrixmlsim[xs].S2hapldiv,(double)matrixmlsim[xs].Shapldiv,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance of haplotype diversity");
                                        strcpy(name_sim,"Exp. variance of haplotype diversity");
                                        break;
                                    case 6:
                                        /*average Rm*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].SRm/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].SRm/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average of Rm");
                                        strcpy(name_sim,"Exp. average of Rm");
                                        break;
                                    case 13:
                                        /*variance Rm*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances((double)matrixml[0].S2Rm,(double)matrixml[0].SRm,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit((double)matrixmlsim[xs].S2Rm,(double)matrixmlsim[xs].SRm,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance of Rm");
                                        strcpy(name_sim,"Exp. variance of Rm");
                                        break;
                                }
                                if((h = start_histogram(vector_obs,vector_sim,x+1,xs,name_obs,name_sim,file_output)) == 0) {
                                    printf("\nHistogram for %s is not available.\n",name_obs);
                                    if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                                }
                                else if(h == 2) z = 14;
                                strcpy(name_obs,"");
                                strcpy(name_sim,"");
                            }
                            free(vector_obs);
                            free(vector_sim);
                            break;
                        case '1':
                            for(z=0;z<6;z++) {
                                x = 0;
                                if(z == 3 && data[0].n_loci < 3) break;
                                switch(z) {
                                    case 0:
                                        /*average za*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Sza/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Sza/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Rozas' za (not divided by S)");
                                        strcpy(name_sim,"Exp. average Rozas' za (not divided by S)");
                                        break;
                                    case 3:
                                        /*variance za*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2za,matrixml[0].Sza,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2za,matrixmlsim[xs].Sza,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Rozas' za (not divided by S)");
                                        strcpy(name_sim,"Exp. variance Rozas' za (not divided by S)");
                                        break;
                                    case 1:
                                        /*average b*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Sb/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Sb/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Wall's b (not divided by S)");
                                        strcpy(name_sim,"Exp. average Wall's b (not divided by S)");
                                        break;
                                    case 4:
                                        /*variance b*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances((double)matrixml[0].S2b,(double)matrixml[0].Sb,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit((double)matrixmlsim[xs].S2b,(double)matrixmlsim[xs].Sb,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Wall's b (not divided by S)");
                                        strcpy(name_sim,"Exp. variance Wall's b (not divided by S)");
                                        break;
                                    case 2:
                                        /*average q*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Sq/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Sq/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Wall's q (not divided by S)");
                                        strcpy(name_sim,"Exp. average Wall's q (not divided by S)");
                                        break;
                                    case 5:
                                        /*variance q*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances((double)matrixml[0].S2q,(double)matrixml[0].Sq,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit((double)matrixmlsim[xs].S2q,(double)matrixmlsim[xs].Sq,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Wall's q (not divided by S)");
                                        strcpy(name_sim,"Exp. variance Wall's q (not divided by S)");
                                        break;
                                }
                                if((h = start_histogram(vector_obs,vector_sim,x+1,xs,name_obs,name_sim,file_output)) == 0) {
                                    printf("\nHistogram for %s is not available.\n",name_obs);
                                    if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                                }
                                else if(h == 2) z = 6;
                                strcpy(name_obs,"");
                                strcpy(name_sim,"");
                            }
                            free(vector_obs);
                            free(vector_sim);
                            break;
                        case '2':
                            for(z=0;z<14;z++) {
                                if(z == 6 && data[0].n_loci < 3) break;
                                if(data[0].time_spec > 0.) if(z == 4) z = 5;
                                if(data[0].time_spec > 0.) if(z == 11) z = 12;
                                if(!(data[0].time_spec > 0.)) if(z == 2) z = 4;
                                if(!(data[0].time_spec > 0.)) 
                                    if(z == 5) {
                                        if(data[0].n_loci > 2) z = 7;
                                        else break;
                                    }
                                if(!(data[0].time_spec > 0.)) if(z == 9) z = 11;
                                if(!(data[0].time_spec > 0.)) if(z == 12) break;
                                x = 0;
                                switch(z) {
                                    case 0:
                                        /*average theta_wat*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Stheta_wat/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Stheta_wat/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Theta (Watterson) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. average Theta (Watterson) corrected for chromosome population size");
                                        break;
                                    case 7:
                                        /*variance theta_wat*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2theta_wat,matrixml[0].Stheta_wat,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2theta_wat,matrixmlsim[xs].Stheta_wat,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Theta (Watterson) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. variance Theta (Watterson) corrected for chromosome population size");
                                        break;
                                    case 1:
                                        /*average theta_taj*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Stheta_taj/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Stheta_taj/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Theta (Tajima) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. average Theta (Tajima) corrected for chromosome population size");
                                        break;
                                    case 8:
                                        /*variance theta_taj*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2theta_taj,matrixml[0].Stheta_taj,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2theta_taj,matrixmlsim[xs].Stheta_taj,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Theta (Tajima) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. variance Theta (Tajima) corrected for chromosome population size");
                                        break;
                                    case 2:
                                        /*average theta_fw*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Stheta_fw/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Stheta_fw/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Theta (normalized Fay and Wu) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. average Theta (normalized Fay and Wu) corrected for chromosome population size");
                                        break;
                                    case 9:
                                        /*variance theta_fw*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2theta_fw,matrixml[0].Stheta_fw,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2theta_fw,matrixmlsim[xs].Stheta_fw,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Theta (normalized Fay and Wu) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. variance Theta (normalized Fay and Wu) corrected for chromosome population size");
                                        break;
                                    case 3:
                                        /*average theta_fuli*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Stheta_fuli/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Stheta_fuli/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Theta (Fu and Li) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. average Theta (Fu and Li) corrected for chromosome population size");
                                        break;
                                    case 10:
                                        /*variance theta_fuli*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2theta_fuli,matrixml[0].Stheta_fuli,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2theta_fuli,matrixmlsim[xs].Stheta_fuli,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Theta (Fu and Li) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. variance Theta (Fu and Li) corrected for chromosome population size");
                                        break;
                                    case 4:
                                        /*average theta_fulin*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Stheta_fulin/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Stheta_fulin/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Theta (Fu and Li) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. average Theta (Fu and Li) corrected for chromosome population size");
                                        break;
                                    case 11:
                                        /*variance theta_fulin*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2theta_fulin,matrixml[0].Stheta_fulin,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2theta_fulin,matrixmlsim[xs].Stheta_fulin,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Theta (Fu and Li) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. variance Theta (Fu and Li) corrected for chromosome population size");
                                        break;
                                    case 5:
                                        /*average theta_L*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Stheta_L/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Stheta_L/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Theta (Zeng et al.) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. average Theta (Zeng et al.) corrected for chromosome population size");
                                        break;
                                    case 12:
                                        /*variance theta_L*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2theta_L,matrixml[0].Stheta_L,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2theta_L,matrixmlsim[xs].Stheta_L,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Theta (Zeng et al.) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. variance Theta (Zeng et al.) corrected for chromosome population size");
                                        break;
                                    case 6:
                                        /*average ndivergence*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Sndivergence/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Sndivergence/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average divergence");
                                        strcpy(name_sim,"Exp. average divergence");
                                        break;
                                    case 13:
                                        /*variance ndivergence*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2ndivergence,matrixml[0].Sndivergence,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2ndivergence,matrixmlsim[xs].Sndivergence,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance divergence");
                                        strcpy(name_sim,"Exp. variance divergence");
                                        break;
                                }
                                if((h = start_histogram(vector_obs,vector_sim,x+1,xs,name_obs,name_sim,file_output)) == 0) {
                                    printf("\nHistogram for %s is not available.\n",name_obs);
                                    if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                                }
                                else if(h == 2) z = 14;
                                strcpy(name_obs,"");
                                strcpy(name_sim,"");
                            }
                            free(vector_obs);
                            free(vector_sim);
                            break;
                        case '3':
                            for(z=0;z<14;z++) {
                                if(z == 6 && data[0].n_loci < 3) break;
                                if(data[0].time_spec > 0.) if(z == 4) z = 5;
                                if(data[0].time_spec > 0.) if(z == 11) z = 12;
                                if(!(data[0].time_spec > 0.)) if(z == 2) z = 4;
                                if(!(data[0].time_spec > 0.)) 
                                    if(z == 5) {
                                        if(data[0].n_loci > 2) z = 7;
                                        else break;
                                    }
                                if(!(data[0].time_spec > 0.)) if(z == 9) z = 11;
                                if(!(data[0].time_spec > 0.)) if(z == 12) break;
                                x = 0;
                                switch(z) {
                                    case 0:
                                        /*average theta_wat_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Stheta_wat_nut/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Stheta_wat_nut/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Theta/nt (Watterson) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. average Theta/nt (Watterson) corrected for chromosome population size");
                                        break;
                                    case 7:
                                        /*variance theta_wat_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2theta_wat_nut,matrixml[0].Stheta_wat_nut,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2theta_wat_nut,matrixmlsim[xs].Stheta_wat_nut,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Theta/nt (Watterson) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. variance Theta/nt (Watterson) corrected for chromosome population size");
                                        break;
                                    case 1:
                                        /*average theta_taj_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Stheta_taj_nut/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Stheta_taj_nut/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Theta/nt (Tajima) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. average Theta/nt (Tajima) corrected for chromosome population size");
                                        break;
                                    case 8:
                                        /*variance theta_taj_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2theta_taj_nut,matrixml[0].Stheta_taj_nut,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2theta_taj_nut,matrixmlsim[xs].Stheta_taj_nut,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Theta/nt (Tajima) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. variance Theta/nt (Tajima) corrected for chromosome population size");
                                        break;
                                    case 2:
                                        /*average theta_fw_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Stheta_fw_nut/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Stheta_fw_nut/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Theta/nt (normalized Fay and Wu) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. average Theta/nt (normalized Fay and Wu) corrected for chromosome population size");
                                        break;
                                    case 9:
                                        /*variance theta_fw_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2theta_fw_nut,matrixml[0].Stheta_fw_nut,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2theta_fw_nut,matrixmlsim[xs].Stheta_fw_nut,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Theta/nt (normalized Fay and Wu) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. variance Theta/nt (normalized Fay and Wu) corrected for chromosome population size");
                                        break;
                                    case 3:
                                        /*average theta_fuli_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Stheta_fuli_nut/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Stheta_fuli_nut/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Theta/nt (Fu and Li) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. average Theta/nt (Fu and Li) corrected for chromosome population size");
                                        break;
                                    case 10:
                                        /*variance theta_fuli_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2theta_fuli_nut,matrixml[0].Stheta_fuli_nut,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2theta_fuli_nut,matrixmlsim[xs].Stheta_fuli_nut,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Theta/nt (Fu and Li) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. variance Theta/nt (Fu and Li) corrected for chromosome population size");
                                        break;
                                    case 4:
                                        /*average theta_fulin_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Stheta_fulin_nut/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Stheta_fulin_nut/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Theta/nt (Fu and Li) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. average Theta/nt (Fu and Li) corrected for chromosome population size");
                                        break;
                                    case 11:
                                        /*variance theta_fulin_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2theta_fulin_nut,matrixml[0].Stheta_fulin_nut,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2theta_fulin_nut,matrixmlsim[xs].Stheta_fulin_nut,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Theta/nt (Fu and Li) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. variance Theta/nt (Fu and Li) corrected for chromosome population size");
                                        break;
                                    case 5:
                                        /*average theta_L_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Stheta_L_nut/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Stheta_L_nut/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average Theta/nt (Zeng et al.) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. average Theta/nt (Zeng et al.) corrected for chromosome population size");
                                        break;
                                    case 12:
                                        /*variance theta_L_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2theta_L_nut,matrixml[0].Stheta_L_nut,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2theta_L_nut,matrixmlsim[xs].Stheta_L_nut,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance Theta/nt (Zeng et al.) corrected for chromosome population size");
                                        strcpy(name_sim,"Exp. variance Theta/nt (Zeng et al.) corrected for chromosome population size");
                                        break;
                                    case 6:
                                        /*average ndivergence_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = (double)matrixml[0].Sndivergence_nut/(double)matrixml[0].nloci;
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = (double)matrixmlsim[xs].Sndivergence_nut/(double)data[0].n_loci;
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. average divergence/nt");
                                        strcpy(name_sim,"Exp. average divergence/nt");
                                        break;
                                    case 13:
                                        /*variance ndivergence_nut*/
                                        if(dataobsequalsim) vector_obs[x=0] = variances(matrixml[0].S2ndivergence_nut,matrixml[0].Sndivergence_nut,matrixml[0].nloci);
                                        for(xs=0;xs<(long int)data[0].n_iter;xs++) vector_sim[xs] = varianceit(matrixmlsim[xs].S2ndivergence_nut,matrixmlsim[xs].Sndivergence_nut,data[0].n_loci);
                                        if(dataobsequalsim) strcpy(name_obs,"Obs. variance divergence/nt");
                                        strcpy(name_sim,"Exp. variance divergence/nt");
                                        break;
                                }
                                if((h = start_histogram(vector_obs,vector_sim,x+1,xs,name_obs,name_sim,file_output)) == 0) {
                                    printf("\nHistogram for %s is not available.\n",name_obs);
                                    if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                                }
                                else if(h == 2) z = 12;
                                strcpy(name_obs,"");
                                strcpy(name_sim,"");
                            }
                            free(vector_obs);
                            free(vector_sim);
                            break;
                        case '4':
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
						fprintf(file_output,"\n     3.1.1.0.2. Display histograms of simulated results.\n");
					}

                    vector_obs = 0;
                    name_obs[0] = '\0';
                    if(dataobsequalsim) 
                        if((vector_obs = (double *) malloc(sizeof(double))) == 0) {
                            printf("\nError: memory not reallocated. Histograms not displayed. \n");
                            if(file_output) fputs("\nError: memory not reallocated. Histograms not displayed. \n",file_output);
                            return;
                        }
                    if((vector_sim = (double *) malloc((data[0].n_iter)*sizeof(double))) == 0) {
                        printf("\nError: memory not reallocated. Histograms not displayed. \n");
                        if(file_output) fputs("\nError: memory not reallocated. Histograms not displayed. \n",file_output);
                        return;
                    }
                    for(z=0;z<25;z++) {
                        if(z == 13 && data[0].n_loci < 3) break;
						if(data[0].n_loci == 1) if(z==12) break;
                        if(!(data[0].time_spec > 0.)) if(z == 3) z = 5;
                        if(!(data[0].time_spec > 0.)) if(z == 10) z = 11;
						if(!(data[0].time_spec > 0.)) if(z == 12) z = 13;
                        if(!(data[0].time_spec > 0.)) 
                            if(z == 11) {
                                if(data[0].n_loci > 2) z = 12;
                                else break;
                            }
                        if(!(data[0].time_spec > 0.)) if(z == 16) z = 18;
                        if(!(data[0].time_spec > 0.)) if(z == 23) z = 24;
                        switch(z) {
                            case 0:
                                /*average Tajima's D*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].StajimaD != (double) -10000) {
                                        vector_obs[y] = matrixml[0].StajimaD/(double)matrixml[0].nltajD;
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. average Tajima's D");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].StajimaD != (double) -10000) {
                                        vector_sim[ys] = matrixmlsim[xs].StajimaD/(double)matrixmlsim[xs].nltajD;
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. average Tajima's D");
                                break;
                            case 13:
                                /*variance Tajima's D*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].S2tajimaD != (double) -10000) {
                                        vector_obs[y] = variances(matrixml[0].S2tajimaD,matrixml[0].StajimaD,matrixml[0].nltajD);
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. variance Tajima's D");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].S2tajimaD != (double) -10000) {
                                        vector_sim[ys] = varianceit(matrixmlsim[xs].S2tajimaD,matrixmlsim[xs].StajimaD,matrixmlsim[xs].nltajD);
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. variance Tajima's D");
                                break;
                            case 1:
                                /*average FuLi's D*/
                                y = 0;
                                if(data[0].time_spec > 0.) {
                                    if(dataobsequalsim)
                                        if(matrixml[0].SfuliD != (double) -10000) {
                                            vector_obs[y] = matrixml[0].SfuliD/(double)matrixml[0].nlflD;
                                            y++;
                                        }
                                    ys = 0;
                                    for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                        if(matrixmlsim[xs].SfuliD != (double) -10000) {
                                            vector_sim[ys] = matrixmlsim[xs].SfuliD/(double)matrixmlsim[xs].nlflD;
                                            ys++;
                                        }
                                    }
                                    if(dataobsequalsim) strcpy(name_obs,"Obs. average Fu and Li's D");
                                    strcpy(name_sim,"Exp. average Fu and Li's D");
                                }
                                else {
                                    if(dataobsequalsim) 
                                        if(matrixml[0].SfuliDn != (double) -10000) {
                                            vector_obs[y] = matrixml[0].SfuliDn/(double)matrixml[0].nlflDn;
                                            y++;
                                        }
                                    ys = 0;
                                    for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                        if(matrixmlsim[xs].SfuliDn != (double) -10000) {
                                            vector_sim[ys] = matrixmlsim[xs].SfuliDn/(double)matrixmlsim[xs].nlflDn;
                                            ys++;
                                        }
                                    }
                                    if(dataobsequalsim) strcpy(name_obs,"Obs. average Fu and Li's D*");
                                    strcpy(name_sim,"Exp. average Fu and Li's D*");
                                }
                                break;
                            case 14:
                                /*variance FuLi's D*/
                                y = 0;
                                if(data[0].time_spec > 0.) {
                                    if(dataobsequalsim) 
                                        if(matrixml[0].S2fuliD != (double) -10000) {
                                            vector_obs[y] = variances(matrixml[0].S2fuliD,matrixml[0].SfuliD,matrixml[0].nlflD);
                                            y++;
                                        }
                                    ys = 0;
                                    for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                        if(matrixmlsim[xs].S2fuliD != (double) -10000) {
                                            vector_sim[ys] = varianceit(matrixmlsim[xs].S2fuliD,matrixmlsim[xs].SfuliD,matrixmlsim[xs].nlflD);
                                            ys++;
                                        }
                                    }
                                    if(dataobsequalsim) strcpy(name_obs,"Obs. variance Fu and Li's D");
                                    strcpy(name_sim,"Exp. variance Fu and Li's D");
                                }
                                else {
                                    if(dataobsequalsim) 
                                        if(matrixml[0].S2fuliDn != (double) -10000) {
                                            vector_obs[y] = variances(matrixml[0].S2fuliDn,matrixml[0].SfuliDn,matrixml[0].nlflDn);
                                            y++;
                                        }
                                    ys = 0;
                                    for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                        if(matrixmlsim[xs].S2fuliDn != (double) -10000) {
                                            vector_sim[ys] = varianceit(matrixmlsim[xs].S2fuliDn,matrixmlsim[xs].SfuliDn,matrixmlsim[xs].nlflDn);
                                            ys++;
                                        }
                                    }
                                    if(dataobsequalsim) strcpy(name_obs,"Obs. variance Fu and Li's D*");
                                    strcpy(name_sim,"Exp. variance Fu and Li's D*");
                                }
                                break;
                            case 2:
                                /*average FuLi's F*/
                                y = 0;
                                if(data[0].time_spec > 0.) {
                                    if(dataobsequalsim) 
                                        if(matrixml[0].SfuliF != (double) -10000) {
                                            vector_obs[y] = matrixml[0].SfuliF/(double)matrixml[0].nlflF;
                                            y++;
                                        }
                                    ys = 0;
                                    for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                        if(matrixmlsim[xs].SfuliF != (double) -10000) {
                                            vector_sim[ys] = matrixmlsim[xs].SfuliF/(double)matrixmlsim[xs].nlflF;
                                            ys++;
                                        }
                                    }
                                    if(dataobsequalsim) strcpy(name_obs,"Obs. average Fu and Li's F");
                                    strcpy(name_sim,"Exp. average Fu and Li's F");
                                }
                                else {
                                    if(dataobsequalsim) 
                                        if(matrixml[0].SfuliFn != (double) -10000) {
                                            vector_obs[y] = matrixml[0].SfuliFn/(double)matrixml[0].nlflFn;
                                            y++;
                                        }
                                    ys = 0;
                                    for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                        if(matrixmlsim[xs].SfuliFn != (double) -10000) {
                                            vector_sim[ys] = matrixmlsim[xs].SfuliFn/(double)matrixmlsim[xs].nlflFn;
                                            ys++;
                                        }
                                    }
                                    if(dataobsequalsim) strcpy(name_obs,"Obs. average Fu and Li's F*");
                                    strcpy(name_sim,"Exp. average Fu and Li's F*");
                                }
                                break;
                            case 15:
                                /*variance FuLi's F*/
                                y = 0;
                                if(data[0].time_spec > 0.) {
                                    if(dataobsequalsim) 
                                        if(matrixml[0].S2fuliF != (double) -10000) {
                                            vector_obs[y] = variances(matrixml[0].S2fuliF,matrixml[0].SfuliF,matrixml[0].nlflF);
                                            y++;
                                        }
                                    ys = 0;
                                    for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                        if(matrixmlsim[xs].S2fuliF != (double) -10000) {
                                            vector_sim[ys] = varianceit(matrixmlsim[xs].S2fuliF,matrixmlsim[xs].SfuliF,matrixmlsim[xs].nlflF);
                                            ys++;
                                        }
                                    }
                                    if(dataobsequalsim) strcpy(name_obs,"Obs. variance Fu and Li's F");
                                    strcpy(name_sim,"Exp. variance Fu and Li's F");
                                }
                                else {
                                    if(dataobsequalsim) 
                                        if(matrixml[0].S2fuliFn != (double) -10000) {
                                            vector_obs[y] = variances(matrixml[0].S2fuliFn,matrixml[0].SfuliFn,matrixml[0].nlflFn);
                                            y++;
                                        }
                                    ys = 0;
                                    for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                        if(matrixmlsim[xs].S2fuliFn != (double) -10000) {
                                            vector_sim[ys] = varianceit(matrixmlsim[xs].S2fuliFn,matrixmlsim[xs].SfuliFn,matrixmlsim[xs].nlflFn);
                                            ys++;
                                        }
                                    }
                                    if(dataobsequalsim) strcpy(name_obs,"Obs. variance Fu and Li's F*");
                                    strcpy(name_sim,"Exp. variance Fu and Li's F*");
                                }
                                break;
                            case 3:
                                /*average FayWu H*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].SfaywuH != (double) -10000) {
                                        vector_obs[y] = matrixml[0].SfaywuH/(double)matrixml[0].nlH;
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. average normalized Fay and Wu's H");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].SfaywuH != (double) -10000) {
                                        vector_sim[ys] = matrixmlsim[xs].SfaywuH/(double)matrixmlsim[xs].nlH;
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. average normalized Fay and Wu's H");
                                break;
                            case 16:
                                /*variance FayWu H*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].S2faywuH != (double) -10000) {
                                        vector_obs[y] = variances(matrixml[0].S2faywuH,matrixml[0].SfaywuH,matrixml[0].nlH);
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. variance normalized Fay and Wu's H");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].S2faywuH != (double) -10000) {
                                        vector_sim[ys] = varianceit(matrixmlsim[xs].S2faywuH, matrixmlsim[xs].SfaywuH,matrixmlsim[xs].nlH);
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. variance normalized Fay and Wu's H");
                                break;
                            case 4:
                                /*average FayWu Ho*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].SfaywuHo != (double) -10000) {
                                        vector_obs[y] = matrixml[0].SfaywuHo/(double)matrixml[0].nlHo;
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. average Fay and Wu's H");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].SfaywuHo != (double) -10000) {
                                        vector_sim[ys] = matrixmlsim[xs].SfaywuHo/(double)matrixmlsim[xs].nlHo;
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. average Fay and Wu's H");
                                break;
                            case 17:
                                /*variance FayWu Ho*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].S2faywuHo != (double) -10000) {
                                        vector_obs[y] = variances(matrixml[0].S2faywuHo,matrixml[0].SfaywuHo,matrixml[0].nlHo);
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. variance Fay and Wu's H");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].S2faywuHo != (double) -10000) {
                                        vector_sim[ys] = varianceit(matrixmlsim[xs].S2faywuHo, matrixmlsim[xs].SfaywuHo,matrixmlsim[xs].nlHo);
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. variance Fay and Wu's H");
                                break;
                            case 5:
                                /*average Fu Fs*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].SfuFs != (double) -10000) {
                                        vector_obs[y] = matrixml[0].SfuFs/(double)matrixml[0].nlFs;
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. average Fu's Fs");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].SfuFs != (double) -10000) {
                                        vector_sim[ys] = matrixmlsim[xs].SfuFs/(double)matrixmlsim[xs].nlFs;
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. average Fu's Fs");
                                break;
                            case 18:
                                /*variance Fu Fs*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].S2fuFs != (double) -10000) {
                                        vector_obs[y] = variances(matrixml[0].S2fuFs,matrixml[0].SfuFs,matrixml[0].nlFs);
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. variance Fu's Fs");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].S2fuFs != (double) -10000) {
                                        vector_sim[ys] = varianceit(matrixmlsim[xs].S2fuFs, matrixmlsim[xs].SfuFs,matrixmlsim[xs].nlFs);
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. variance Fu's Fs");
                                break;
                            case 6:
                                /*average Rozas ZA*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].SrZA != (double) -10000) {
                                        vector_obs[y] = matrixml[0].SrZA/(double)matrixml[0].nlZ;
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. average Rozas' ZA");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].SrZA != (double) -10000) {
                                        vector_sim[ys] = matrixmlsim[xs].SrZA/(double)matrixmlsim[xs].nlZ;
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. average Rozas' ZA");
                                break;
                            case 19:
                                /*variance Rozas ZA*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].S2rZA != (double) -10000) {
                                        vector_obs[y] = variances(matrixml[0].S2rZA,matrixml[0].SrZA,matrixml[0].nlZ);
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. variance Rozas' ZA");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].S2rZA != (double) -10000) {
                                        vector_sim[ys] = varianceit(matrixmlsim[xs].S2rZA, matrixmlsim[xs].SrZA,matrixmlsim[xs].nlZ);
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. variance Rozas' ZA");
                                break;
                            case 7:
                                /*average Wall's B*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].SwB != (double) -10000) {
                                        vector_obs[y] = matrixml[0].SwB/(double)matrixml[0].nlB;
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. average Wall's B");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].SwB != (double) -10000) {
                                        vector_sim[ys] = matrixmlsim[xs].SwB/(double)matrixmlsim[xs].nlB;
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. average Wall's B");
                                break;
                            case 20:
                                /*variance Wall's B*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].S2wB != (double) -10000) {
                                        vector_obs[y] = variances(matrixml[0].S2wB,matrixml[0].SwB,matrixml[0].nlB);
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. variance Wall's B");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].S2wB != (double) -10000) {
                                        vector_sim[ys] = varianceit(matrixmlsim[xs].S2wB, matrixmlsim[xs].SwB,matrixmlsim[xs].nlB);
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. variance Wall's B");
                                break;
                            case 8:
                                /*average Wall's Q*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].SwQ != (double) -10000) {
                                        vector_obs[y] = matrixml[0].SwQ/(double)matrixml[0].nlQ;
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. average Wall's Q");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].SwQ != (double) -10000) {
                                        vector_sim[ys] = matrixmlsim[xs].SwQ/(double)matrixmlsim[xs].nlQ;
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. average Wall's Q");
                                break;
                            case 21:
                                /*variance Wall's Q*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].S2wQ != (double) -10000) {
                                        vector_obs[y] = variances(matrixml[0].S2wQ,matrixml[0].SwQ,matrixml[0].nlQ);
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. variance Wall's Q");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].S2wQ != (double) -10000) {
                                        vector_sim[ys] = varianceit(matrixmlsim[xs].S2wQ, matrixmlsim[xs].SwQ,matrixmlsim[xs].nlQ);
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. variance Wall's Q");
                                break;
                            case 9:
                                /*average R2*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].SR2 != (double) -10000) {
                                        vector_obs[y] = matrixml[0].SR2/(double)matrixml[0].nlR2;
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. average R2");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].SR2 != (double) -10000) {
                                        vector_sim[ys] = matrixmlsim[xs].SR2/(double)matrixmlsim[xs].nlR2;
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. average R2");
                                break;
                            case 22:
                                /*variance R2*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].S2R2 != (double) -10000) {
                                        vector_obs[y] = variances(matrixml[0].S2R2,matrixml[0].SR2,matrixml[0].nlR2);
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. variance R2");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].S2R2 != (double) -10000) {
                                        vector_sim[ys] = varianceit(matrixmlsim[xs].S2R2, matrixmlsim[xs].SR2,matrixmlsim[xs].nlR2);
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. variance R2");
                                break;
                            case 10:
                                /*average ZengE*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].SzengE != (double) -10000) {
                                        vector_obs[y] = matrixml[0].SzengE/(double)matrixml[0].nlE;
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. average Zeng et al. E");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].SzengE != (double) -10000) {
                                        vector_sim[ys] = matrixmlsim[xs].SzengE/(double)matrixmlsim[xs].nlE;
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. average Zeng et al. E");
                                break;
                            case 23:
                                /*variance ZengE*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].S2zengE != (double) -10000) {
                                        vector_obs[y] = variances(matrixml[0].S2zengE,matrixml[0].SzengE,matrixml[0].nlE);
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. variance Zeng et al. E");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].S2zengE != (double) -10000) {
                                        vector_sim[ys] = varianceit(matrixmlsim[xs].S2zengE, matrixmlsim[xs].SzengE,matrixmlsim[xs].nlE);
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. variance Zeng et al. E");
                                break;
                            case 11:
                                /*average EW*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].Sewtest != (double) -10000) {
                                        vector_obs[y] = matrixml[0].Sewtest/(double)*n_loci;
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. average EW test");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].Sewtest != (double) -10000) {
                                        vector_sim[ys] = matrixmlsim[xs].Sewtest/(double)*n_loci;
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. average EW test");
                                break;
                            case 24:
                                /*variance EW*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].S2ewtest != (double) -10000) {
                                        vector_obs[y] = variances(matrixml[0].S2ewtest,matrixml[0].Sewtest,(int)*n_loci);
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. variance EW test");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].S2ewtest != (double) -10000) {
                                        vector_sim[ys] = varianceit(matrixmlsim[xs].S2ewtest, matrixmlsim[xs].Sewtest,(int)*n_loci);
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. variance EW test");
                                break;
                            case 12:
                                /*HKA tests*/
                                y = 0;
                                if(dataobsequalsim) {
                                    if(matrixml[0].Shka != (double) -10000) {
                                        vector_obs[y] = matrixml[0].Shka/(double)*n_loci;
                                        y++;
                                    }
                                    strcpy(name_obs,"Obs. Chi-sq HKA(JC)");
                                }
                                ys = 0;
                                for(xs=0;xs<(long int)data[0].n_iter;xs++) {
                                    if(matrixmlsim[xs].Shka != (double) -10000) {
                                        vector_sim[ys] = matrixmlsim[xs].Shka/(double)data[0].n_loci;
                                        ys++;
                                    }
                                }
                                strcpy(name_sim,"Exp. Chi-sq HKA(JC)");
                                break;
                        }
                        if((h = start_histogram(vector_obs,vector_sim,y,ys,name_obs,name_sim,file_output)) == 0) {
                            printf("\nHistogram for %s is not available.\n",name_obs);
                            if(file_output) fprintf(file_output,"\nHistogram for %s is not available.\n",name_obs);
                        }
                        else if(h == 2) z = 24;
                        strcpy(name_obs,"");
                        strcpy(name_sim,"");
                    }
                    free(vector_obs);
                    free(vector_sim);
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

