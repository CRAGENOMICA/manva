/*
 *  open_main.c
 *  MuLoNeTests
 *
 *  Created by sebas on Fri Feb 21 2003.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>

int main (/*int argc,char *argv[]*/)
{	
	int x;
    static struct globvar *global = 0;/*global variables*/
    static struct statistics *matrix = 0;/*matrix with all information for observed data*/
    static struct statmulo *matrixml = 0;/*matrix with information on multilocus observed data*/
    static struct var *datams = 0;/*vector with global data for simulations*/
    static struct var2b *inputms = 0;/*matrix with data for simulations in each locus*/
    static struct statistisim **matrixsim = 0;/*simulations (lot of memory. recommended: nloci x niterations < 10to5)*/
    static struct statistisimmuloc *matrixmlsim = 0;/*multilocus simulation data*/
    static struct horizontalstatsml *avgstatloci = 0; /*matrix for avg or median for sign test analysis in each locus*/    
    static struct LRTdist *lrdist = 0; /*vector with LRT distributions for theta models*/
	/*static struct MLthetasim **mthetasim = 0; *//*matrix with simulated theta under all models*/
	FILE *file_output = 0;
    

    int open_mainmenu(struct globvar **,struct statistics **,struct statmulo **,struct var **, struct  var2b  **,
        struct statistisim ***, struct statistisimmuloc **, struct horizontalstatsml **,struct LRTdist ** /*,struct MLthetasim ****/,FILE **);
    
    if((global = (struct globvar *) calloc(1,sizeof(struct globvar))) == 0) {
        puts("\nError: memory not reallocated. main.0 \n");
        exit(0);
    }
    if((matrix = (struct statistics *) calloc(1,sizeof(struct statistics))) == 0) {
        puts("\nError: memory not reallocated. main.1 \n");
        exit(0);
    }
    if((matrixml = (struct statmulo *) calloc(1,sizeof(struct statmulo))) == 0) {
        puts("\nError: memory not reallocated. main.2 \n");
        exit(0);
    }
    if((inputms = (struct var2b *) calloc(1,sizeof(struct var2b))) == 0) {
        puts("\nError: memory not reallocated. main.4 \n");
        exit(0);
    }
    if((matrixsim = (struct statistisim **) calloc(1,sizeof(struct statistisim *))) == 0) {
        puts("\nError: memory not reallocated. main.5 \n");
        exit(0);
    }
    if((matrixmlsim = (struct statistisimmuloc *) calloc(1,sizeof(struct statistisimmuloc))) == 0) {
        puts("\nError: memory not reallocated. main.6 \n");
        exit(0);
    }
    if((avgstatloci = (struct horizontalstatsml *) calloc(1,sizeof(struct horizontalstatsml))) == 0) {
        puts("\nError: memory not reallocated. main.7 \n");
        exit(0);
    }
    if((datams = (struct var *) calloc(1,sizeof(struct var))) == 0) {
        puts("\nError: memory not reallocated. main.3 \n");
        exit(0);
    }
    if((inputms[0].factor_pop = (double *) calloc(1,sizeof(double))) == 0) {
        puts("\nError: memory not reallocated. main.8 \n");
        exit(0);
    }
    if((inputms[0].nrec = (double *) calloc(2,sizeof(double))) == 0) {
        puts("\nError: memory not reallocated. main.9 \n");
        exit(0);
    }
    if((inputms[0].npast = (double *) calloc(2,sizeof(double))) == 0) {
        puts("\nError: memory not reallocated. main.10 \n");
        exit(0);
    }
    if((inputms[0].tpast = (double *) calloc(2,sizeof(double))) == 0) {
        puts("\nError: memory not reallocated. main.11 \n");
        exit(0);
    }
    if((inputms[0].freq = (double *) calloc(1,sizeof(double))) == 0) {
        puts("\nError: memory not reallocated. main.12 \n");
        exit(0);
    }
    if((lrdist = (struct LRTdist *) calloc(1,sizeof(struct LRTdist))) == 0) {
        puts("\nError: memory not reallocated. main.12 \n");
        exit(0);
    }
	/*
    if((mthetasim = (struct MLthetasim **) calloc(1,sizeof(struct MLthetasim *))) == 0) {
        puts("\nError: memory not reallocated. main.12 \n");
        exit(0);
    }
	*/

    global[0].observed_data = 0;
    global[0].montecarlo_sim = 0;
    global[0].dataforsim = 0;
    global[0].dataforsimt = 0;
    global[0].outgroup = -1;
    global[0].n_loci = 0;
	/*global[0].maxloc = 0;*/
	global[0].dataindata = 0;
    global[0].nlocimat = 0;
    global[0].nlocimatt = 0;
    global[0].dataequalsim = 0;
    global[0].ml0done = 0;
    global[0].nitermat = 0;
    global[0].nitermatt = 0;
    
	global[0].gfffiles = 0;
	global[0].ifgencode = 0;
	
	inputms[0].nintn = 0;
	
    puts("\n\n\n\n");
    puts(MANVa);
    puts("\n\n\n\n");
    if((open_mainmenu(&global,&matrix,&matrixml,&datams,&inputms,&matrixsim,&matrixmlsim,&avgstatloci,&lrdist/*,&mthetasim*/,&file_output)) == 0) {
		free(lrdist);
		/*free(mthetasim);*/
		free(inputms[0].freq);
		free(inputms[0].tpast);
		free(inputms[0].npast);
		free(inputms[0].nrec);
		free(inputms[0].factor_pop);
		if(global[0].dataforsimt) {
			for(x=1;x<global[0].nlocimat;x++) {
				free(inputms[x].freq);
				free(inputms[x].tpast);
				free(inputms[x].npast);
				free(inputms[x].nrec);
				free(inputms[x].factor_pop);
			}
		}
		free(datams);
		free(avgstatloci);
		free(matrixmlsim);
		if(global[0].montecarlo_sim) 
			for(x=0;x<global[0].nlocimat;x++) 
				free(matrixsim[x]);
		free(matrixsim);
		free(inputms);
		free(matrixml);
		free(matrix);
		free(global);
		
		if(file_output) {
			fprintf(file_output,"\nExit of MANVa program.\n\n");
			fclose(file_output);
		}
	}
	return 0;
}

