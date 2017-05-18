/*
 *  open_coalsimmenu.c
 *  MuLoNeTests
 *
 *  Created by sonsins on Fri Mar 14 2003.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void menu_coalescent_sim(struct globvar **global,struct statistics *matrix,struct statmulo *matrixml,
    struct statistisim ***matrixsim,struct statistisimmuloc **matrixmlsim,struct horizontalstatsml **avgstatloci,
    struct var **data,struct var2b **inputms,struct LRTdist **lrdist/*,struct MLthetasim ***mthetasim*/,FILE *file_output)
{
    char k[1],l[1];
    long int valueld;
    double valuef;
    int x,y;
    long int recom;

    int open_menu_seqinfdata_sim(struct globvar *,struct statistics *,struct statmulo *,struct var **,struct var2b **,
        int * /*,int * */,int *,struct LRTdist **,FILE *);
    int open_menu_evolcond_sim(struct var **,struct var2b **,int * /*,struct globvar * */,struct statistics *,FILE *);
    
    void print_seqinfodata_sim(int,struct var *,struct var2b *, FILE *,int ,struct statistics *);
    void print_evocond_sim(struct var *,struct var2b *, FILE * /*,int */,struct statistics *);
    void print_datasimulations(int,struct var *,struct var2b *, FILE *,int ,struct statistics *);

    int freeallocvectors(int * /*,int * */,struct statistisim *** /*,long int */,int,int,int * /*,FILE * */);
    int reallocvectors(int * /*,int * */,struct statistisim ***,struct statistisimmuloc **,struct horizontalstatsml **,
        long int, int, FILE *);
    int mainms(struct var **,struct var2b **,struct statistisim ***,struct statistisimmuloc **,
        struct horizontalstatsml **,int,
        FILE *);
    void save_simmatrix(struct globvar *,struct statistics *,struct statmulo *,struct var *, struct var2b *,
        struct statistisim **, struct statistisimmuloc *, struct horizontalstatsml *,struct LRTdist * /*,struct MLthetasim ***/,FILE *);
	void erase_results(struct globvar **,struct statistisim ***, struct statistisimmuloc **,struct horizontalstatsml ** /*,struct var2b **,struct LRTdist ** */,struct var **data,FILE *);

    if(global[0][0].montecarlo_sim == 1) {
        printf(" Do you want to overwrite the results from memory (y/n)? \n");
        do *k = getchar();
        while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');
        
        if(*k == 'y' || *k == 'Y') {
            erase_results(global,matrixsim,matrixmlsim,avgstatloci/*,inputms,lrdist*/,data,file_output);
        }
        else return;
    }    
    while(1) {
        #if COMMAND_LINE
        printf("\n\nMENU 3. Statistical inference based on Coalescent Monte Carlo simulations menu:");
        printf("\n     3.0. Perform Coalescent Monte Carlo simulations menu:\n\n");
        
        printf("CHOOSE:\n");
        printf(" 0 - Introduce data and parameters to perform Coalescent simulations.\n");
        printf(" 1 - Run Monte Carlo Coalescent simulations.\n");
        printf(" 2 - Back to Main menu.\n\n");

        if(file_output) {
            fprintf(file_output,"\n\n     MENU:\n     3. Statistical inference based on Coalescent Monte Carlo simulations menu:");
            fprintf(file_output,"\n     3.0. Perform Coalescent Monte Carlo simulations menu:\n\n");
            
            fprintf(file_output," 0 - Introduce data and parameters to perform Coalescent simulations.\n");
            fprintf(file_output," 1 - Run Monte Carlo Coalescent simulations.\n");
            fprintf(file_output," 2 - Back to Main menu.\n\n");
        }
        #endif
    
        do *k = getchar();
        while(*k<'0' || *k>'2');
        
        if(file_output) fprintf(file_output,"OPTION CHOSEN: %c\n\n",*k);

        switch(*k) {
            case '0':
                /*menu for sequence information data*/
                if(global[0][0].dataindata == 1 || global[0][0].dataindata == 3) {
                    print_seqinfodata_sim(global[0][0].dataindata,*data,*inputms,file_output,global[0][0].dataequalsim,matrix);
                    printf("Do you want to use these current data AND the parameter values (data for loci, i.e., nsamples, theta, recombination, etc.) introduced before (y/n)? ");
                    do *l = getchar();
                    while(*l!='y' && *l!='n' && *l!='Y' && *l!='N');
					if(*l == 'n' || *l == 'n') {
						if(freeallocvectors(&(global[0][0].onlymulo)/*,&(global[0][0].maxloc)*/,matrixsim,data[0][0].n_iter/*,data[0][0].n_loci*/,global[0][0].dataindata,&(global[0][0].dataforsim)/*,file_output */) == 1) return;
						if(open_menu_seqinfdata_sim(*global,matrix,matrixml,data,inputms,&(global[0][0].dataequalsim),
							/*&(global[0][0].dataforsim),*/&(global[0][0].dataindata),lrdist,
							file_output) == 0) {
							return;
						}
						print_seqinfodata_sim(global[0][0].dataindata,*data,*inputms,file_output,global[0][0].dataequalsim,matrix);
					}
                }
				else {
					if(freeallocvectors(&(global[0][0].onlymulo)/*,&(global[0][0].maxloc)*/,matrixsim,data[0][0].n_iter/*,data[0][0].n_loci*/,global[0][0].dataindata,&(global[0][0].dataforsim)/*,file_output*/) == 1) return;
					if(open_menu_seqinfdata_sim(*global,matrix,matrixml,data,inputms,&(global[0][0].dataequalsim),
						/*&(global[0][0].dataforsim),*/&(global[0][0].dataindata),lrdist,
						file_output) == 0) {
						return;
					}
					print_seqinfodata_sim(global[0][0].dataindata,*data,*inputms,file_output,global[0][0].dataequalsim,matrix);
				}
                /*menu for evolutionary conditions of simulations*/
                if(open_menu_evolcond_sim(data,inputms,&(global[0][0].dataforsim)/*,*global*/,matrix,file_output) == 0) {
                    return;
                }
                print_evocond_sim(*data,*inputms,file_output/*,global[0][0].dataequalsim*/,matrix);
                break;
            case '1':
                if(global[0][0].dataforsim == 0) {
                    printf("Coalescent simulations only run when all parameters (data + conditions) are included.\n");
                    if(file_output)
                        fputs("Coalescent simulations only run when all parameters (data + conditions) are included.\n",
                            file_output);
                    return;
                }
                /*ask for parameters: n_iter,tlimit,seed1,time_spec*/
                do {
                    if(data[0][0].n_loci <= 10) recom = 10000;
                    else if(data[0][0].n_loci <= 100) recom = 1000;
                    else if(data[0][0].n_loci <= 1000) recom = 1000;
                    else recom = 100;
                    printf("\nHow many iterations per loci (recommended <= %ld)? ",recom);
                    scanf(" %ld",&valueld);
                    data[0][0].n_iter = valueld;
                    if(data[0][0].n_iter <= 0) puts("Error: the value must be positive.");
                }while(data[0][0].n_iter <= 0);
                do {
                    printf("\nSeed for the pseudorandom numbers (positive value)? ");
                    scanf(" %ld",&valueld);
                    data[0][0].seed1 = valueld;
                    if(data[0][0].seed1 <= 0) puts("Error: the value must be positive.");
                }while(data[0][0].seed1 <= 0);
                
                y = 0;
                for(x=0;x<data[0][0].n_loci;x++) {
                    if(inputms[0][x].r > (double)0.) {
                        y = 1;
                        break;
                    }
                }
                if(y == 1) {
                    printf("\nIn order to increase the speed of coalescent simulations,");
                    printf("\n the recombination process can be inactivated after a period of time T in the tree.\n");
                    printf("Do you want to inactivate the recombination process after a time T (y/n)? ");
                    do *l = getchar();
                    while(*l!='y' && *l!='n' && *l!='Y' && *l!='N');
                    if(*l == 'y' || *l == 'Y') {
                        do {
                            puts("\nIndicate Time (in 4N generations) at what recombination is inactivated: ");
                            scanf(" %lg",&valuef);
                            data[0][0].tlimit = valuef;
                            if(data[0][0].tlimit <= (double)0.) puts("Error: the value must be positive.");
                        }while(data[0][0].tlimit <= (double)0.);
                    }
                    else data[0][0].tlimit = 10000.;
                }
                else data[0][0].tlimit = 10000.;
                
                /*reallocate matrix for sim*/
                if(reallocvectors(&(global[0][0].onlymulo)/*,&(global[0][0].maxloc)*/,matrixsim,matrixmlsim,avgstatloci,
                    data[0][0].n_iter,data[0][0].n_loci,file_output) == 1) {
					return;
				}
                if(global[0][0].onlymulo == 1) {
                    puts("\n\nWARNING: ONLY MULTILOCUS STATISTICS ARE AVAILABLE!");
                    if(file_output)
                        fputs("\n\nWARNING: ONLY MULTILOCUS STATISTICS ARE AVAILABLE!",file_output);
                    printf("\n Press 0 to continue.\n\n");
                    
                    do {
                        *l = getchar();
                    }while(*l!='0');
                }
                /*print conditions*/
                print_datasimulations(global[0][0].dataindata,*data,*inputms,file_output,global[0][0].dataequalsim,matrix);

                /*link to Hudson's modified program: COALESCENT PROGRAM*/                
                if(mainms(data,inputms,matrixsim,matrixmlsim,avgstatloci,global[0][0].onlymulo,file_output) == 0) {
                    puts("Error: Coalescent simulations could not be calculated.\n");
                    if(file_output) fputs("Error: Coalescent simulations  could not be calculated.\n",file_output);
                    return;
                }
                /*successful simulations*/
                global[0][0].nitermat = data[0][0].n_iter;
                global[0][0].nlocimat = global[0][0].nlocimatt = data[0][0].n_loci;
                global[0][0].montecarlo_sim = 1;
                
                printf("\n\n Do you want to save all the loaded information in a MANVa file (y/n)? ");
                do {
                    *l = getc(stdin);
                }while(*l !='y' && *l !='Y' && *l !='n' && *l!='N');
                if(*l == 'y' || *l == 'Y') {
                    save_simmatrix(*global,matrix,matrixml,*data,*inputms,*matrixsim,*matrixmlsim,*avgstatloci,*lrdist/*,*mthetasim*/,file_output);
                }
                return;
                break;
            case '2':
                return;
                break;
        }
    }
    return;
}

int open_menu_seqinfdata_sim(struct globvar *global,struct statistics *matrix,struct statmulo *matrixml,
    struct var **data,struct var2b **inputms,int *dataobsequalsim/*,int *dataforsim*/,int *dataindata,struct LRTdist **lrdist,
    FILE *file_output)
{
    int x,y,xx,yy,bestm;
    char k[1],l[1],kk[1],ll[1],mm[1];
    double value,valuef;
    int valuei;
    long int valueli;
	int correctedjc;
    void print_seqinfodata_sim(int,struct var *,struct var2b *, FILE *,int ,struct statistics *);
	void print_seqinfodataobs_sim(/*int,struct var *,struct var2b *,*/ FILE * /*,int */,struct statistics *, struct globvar *);
	
	if(global[0].dataindata) {
		print_seqinfodata_sim(0,*data,*inputms,file_output,global[0].dataequalsim,matrix);
		printf("Do you want to use these current data values to make coalescent simulations (y/n)? ");
		do *mm = getchar();
		while(*mm!='y' && *mm!='n' && *mm!='Y' && *mm!='N');
	}
	else {
		*mm = 'n';
	}
	if(*mm == 'n' || *mm == 'N') {
		/*erase also data from thetaml estimation*/
		if(global[0].dataforsimt == 1) {
			if(global[0].ml0done == 1) {
				if((*lrdist = (struct LRTdist *)realloc(*lrdist,1*sizeof(struct LRTdist))) == 0) {
					puts(" Matrix not reallocated! NOT SUCCESSFUL.");
					if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
					return 0;
				}
				memset(*lrdist,'\0',sizeof(struct LRTdist));
			}
		}
		global[0].ml0done = 0;
		global[0].nlocimatt = 0;
		global[0].nitermatt = 0;
		global[0].dataforsimt = 0;
		global[0].mltheta = 0;
	}
	/*if considering observed values*/
	if(global[0].observed_data == 1) {
		if(*mm == 'n' || *mm == 'N') {
			print_seqinfodataobs_sim(file_output,matrix,global);
			printf("Do you want to use the values of the OBSERVED data to make coalescent simulations (y/n)? ");
			do *ll = getchar();
			while(*ll!='y' && *ll!='n' && *ll!='Y' && *ll!='N');
		}
		else {
			if(*dataindata == 1 || *dataindata == 3) {
				printf("Do you want to modify the parameters of trans/transv, recombination, theta, etc. (y/n)? ");
				do *ll = getchar();
				while(*ll!='y' && *ll!='n' && *ll!='Y' && *ll!='N');
				if(*ll == 'n' || *ll == 'N') *ll = '\0';
				if((*ll == 'y' || *ll == 'Y') &&  *dataobsequalsim == 0) *ll = 'n';
			}
			else {
				if(*dataobsequalsim == 0) *ll = 'n';
				else *ll = 'y';
			}
		}
        if(*ll == 'y' || *ll == 'Y') {
			for(x=1;x<data[0][0].n_loci;x++) {
				free(inputms[0][x].factor_pop);
				free(inputms[0][x].nrec);
				free(inputms[0][x].npast);
				free(inputms[0][x].tpast);
				free(inputms[0][x].freq);
			}
			
            if((inputms[0] = (struct var2b *)realloc(inputms[0],global[0].n_loci*sizeof(struct var2b))) == 0) {
                puts("\nError: memory not reallocated. main.4 \n");
                exit(0);
            }
            if(global[0].outgroup == 1) {
                puts("\nCoalescent simulations assuming outgroup take into account the multiple hit process.");
                do {
                    printf("\nTime (in 2N generations) to the ancestral species between the ingroup and the ancestral population? \n");
                    if(global[0].n_loci > 1) {
                        printf(" Note: HKA estimated time to the ancestral species in 2N generations is %g.\n",matrixml[0].hka_T);
						if(matrixml[0].hka_T < (double)0.) printf(" (Warning: only a positive or zero value is allowed)\n");
					}
                    scanf(" %lg",&valuef);
                    data[0][0].time_spec = valuef; /*time in 2N to the ancestor*/
                    if(valuef < (double)0.) puts("Error: T must be positive or zero.");
                    if(data[0][0].time_spec < (double)0. && valuef > (double)0.) data[0][0].time_spec = (double)0;
                }while(valuef < (double)0.);
                printf("\n The time (in 2N generations) to the ancestral species between ingroup and the ancestral population is %g\n",valuef);
            }
            else {
                data[0][0].time_spec = 0.;/*-1.0*/
            }
            data[0][0].n_loci = global[0].n_loci;

			x=0;
			if((inputms[0][x].nrec = (double *) realloc(inputms[0][x].nrec,2*sizeof(double))) == 0) {
				puts("\nError: memory not reallocated. main.9 \n");
				return(0);
			}
			memset(inputms[0][x].nrec,'\0',2*sizeof(double));
			if((inputms[0][x].npast = (double *) realloc(inputms[0][x].npast,2*sizeof(double))) == 0) {
				puts("\nError: memory not reallocated. main.10 \n");
				return(0);
			}
			memset(inputms[0][x].npast,'\0',2*sizeof(double));
			if((inputms[0][x].tpast = (double *) realloc(inputms[0][x].tpast,2*sizeof(double))) == 0) {
				puts("\nError: memory not reallocated. main.11 \n");
				return(0);
			}
			memset(inputms[0][x].tpast,'\0',2*sizeof(double));
			if((inputms[0][x].factor_pop = (double *) realloc(inputms[0][x].factor_pop,1*sizeof(double))) == 0) {
				puts("\nError: memory not reallocated. main.8 \n");
				return(0);
			}
			memset(inputms[0][x].factor_pop,'\0',1*sizeof(double));
			if((inputms[0][x].freq = (double *) realloc(inputms[0][x].freq,1*sizeof(double))) == 0) {
				puts("\nError: memory not reallocated. main.12 \n");
				return(0);
			}
			memset(inputms[0][x].freq,'\0',1*sizeof(double));
			
			for(x=1;x<data[0][0].n_loci;x++) {
				if((inputms[0][x].nrec = (double *) calloc(2,sizeof(double))) == 0) {
					puts("\nError: memory not reallocated. main.9 \n");
					if(x>1) data[0][0].n_loci = x-1;
					else data[0][0].n_loci = 1;
					if((inputms[0] = (struct var2b *)realloc(inputms[0],data[0][0].n_loci*sizeof(struct var2b))) == 0) {
						puts("\nError: memory not reallocated. seqinfodata.4 \n");
						return(0);
					}			
					return(0);
				}
				if((inputms[0][x].npast = (double *) calloc(2,sizeof(double))) == 0) {
					puts("\nError: memory not reallocated. main.10 \n");
					free(inputms[0][x].nrec);
					if(x>1) data[0][0].n_loci = x-1;
					else data[0][0].n_loci = 1;
					if((inputms[0] = (struct var2b *)realloc(inputms[0],data[0][0].n_loci*sizeof(struct var2b))) == 0) {
						puts("\nError: memory not reallocated. seqinfodata.4 \n");
						return(0);
					}			
					return(0);
				}
				if((inputms[0][x].tpast = (double *) calloc(2,sizeof(double))) == 0) {
					puts("\nError: memory not reallocated. main.11 \n");
					free(inputms[0][x].npast);
					free(inputms[0][x].nrec);
					if(x>1) data[0][0].n_loci = x-1;
					else data[0][0].n_loci = 1;
					if((inputms[0] = (struct var2b *)realloc(inputms[0],data[0][0].n_loci*sizeof(struct var2b))) == 0) {
						puts("\nError: memory not reallocated. seqinfodata.4 \n");
						return(0);
					}			
					return(0);
				}
				if((inputms[0][x].freq = (double *) calloc(1,sizeof(double))) == 0) {
					puts("\nError: memory not reallocated. main.12 \n");
					free(inputms[0][x].tpast);
					free(inputms[0][x].npast);
					free(inputms[0][x].nrec);
					if(x>1) data[0][0].n_loci = x-1;
					else data[0][0].n_loci = 1;
					if((inputms[0] = (struct var2b *)realloc(inputms[0],data[0][0].n_loci*sizeof(struct var2b))) == 0) {
						puts("\nError: memory not reallocated. seqinfodata.4 \n");
						return(0);
					}			
					return(0);
				}
				if((inputms[0][x].factor_pop = (double *) calloc(1,sizeof(double))) == 0) {
					puts("\nError: memory not reallocated. main.8 \n");
					free(inputms[0][x].freq);
					free(inputms[0][x].tpast);
					free(inputms[0][x].npast);
					free(inputms[0][x].nrec);
					if(x>1) data[0][0].n_loci = x-1;
					else data[0][0].n_loci = 1;
					if((inputms[0] = (struct var2b *)realloc(inputms[0],data[0][0].n_loci*sizeof(struct var2b))) == 0) {
						puts("\nError: memory not reallocated. seqinfodata.4 \n");
						return(0);
					}			
					return(0);
				}
			}
			
            for(x=0;x<data[0][0].n_loci;x++) {
				inputms[0][x].npop = 1;
				inputms[0][x].nintn = 0;
                inputms[0][x].factor_chrn = matrix[x].factor_chrn;
                inputms[0][x].nsam = matrix[x].nsamples;
                inputms[0][x].nsites = (long int)matrix[x].nsites;
				inputms[0][x].nmhits = (unsigned long)matrix[x].nmhits;
            }
            puts("\nSimulations are run given the level of variation for each loci.\n");
            y = 0;
            
            do{
                #if COMMAND_LINE
                printf("\n\n     Indicate the level of variation theta (4Nu) chosen:\n");
				printf("     NOTE: The level of variation can be later corrected by a factor in case using Non-standard neutral model\n");
				if(file_output) {
					fprintf(file_output,"\n\n     Indicate the level of variation theta (4Nu) chosen:\n");
					fprintf(file_output,"     NOTE: The level of variation can be later corrected by a factor in case using Non-standard neutral model\n");
				}
                if(data[0][0].time_spec != -1.) {
                    printf("     (In case outgroup, the estimations include the Jukes and Cantor correction).\n");
                    if(file_output) 
                        fprintf(file_output,"\n In case outgroup, the estimations include the Jukes and Cantor correction.\n");
                }
                printf("\n");
                
                printf(" 0 - Watterson's estimation.\n");
                printf(" 1 - Tajima's estimation.\n");
                printf(" 2 - Fu and Li's estimation.\n");
                printf(" 3 - normalized Fay and Wu's estimation (only in case data with outgroup).\n");
                printf(" 4 - HKA test estimation (only in case data with outgroup).\n");
                printf(" 5 - ML estimation (if it was calculated).\n");
                printf(" 6 - Other levels of variation.\n\n");
                #endif
            
                do *k = getchar();
                while(*k<'0' || *k>'6');
				
				if(file_output) {
					fprintf(file_output,"\n");					
					fprintf(file_output," 0 - Watterson's estimation.\n");
					fprintf(file_output," 1 - Tajima's estimation.\n");
					fprintf(file_output," 2 - Fu and Li's estimation.\n");
					fprintf(file_output," 3 - normalized Fay and Wu's estimation (only in case data with outgroup).\n");
					fprintf(file_output," 4 - HKA test estimation (only in case data with outgroup).\n");
					fprintf(file_output," 5 - ML estimation (if it was calculated).\n");
					fprintf(file_output," 6 - Other levels of variation.\n\n");
					fprintf(file_output,"OPTION CHOSEN: %d\n\n",*k);
				}

                switch(*k) {
                    case '0':
                        for(x=0;x<data[0][0].n_loci;x++) {
                            if(matrix[x].noutgroups==0) inputms[0][x].theta = (double)-.75*(double)log((double)1.-(double)4./(double)3.*matrix[x].theta_wat/(double)matrix[x].nsites) * ((double)matrix[x].nsites+(double)inputms[0][x].nmhits);
                            else  {
                                if(matrix[x].theta_wat*(data[0][0].time_spec/2.0+1.0)/(double)matrix[x].nsites >= (double).75) {
                                    printf("\nVariation for locus %d is too high or divergence too far. Outgroup inactivated.",x);
                                    if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far. Outgroup inactivated",x);
                                    inputms[0][x].theta = matrix[x].theta_wat;
									data[0][0].time_spec = 0.;
									return 0;
                                }
                                else inputms[0][x].theta = (double)-.75*(double)log((double)1.-(double)4./(double)3.*matrix[x].theta_wat*(data[0][0].time_spec/2.0+(double)1.0)/(double)matrix[x].nsites)
									* ((double)matrix[x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
                            }
                        }
                        y = 1;
                        break;
                    case '1':
                        for(x=0;x<data[0][0].n_loci;x++) {
                            if(matrix[x].noutgroups==0) inputms[0][x].theta = (double)-.75*(double)log((double)1.-(double)4./(double)3.*matrix[x].theta_taj/(double)matrix[x].nsites) * ((double)matrix[x].nsites+(double)inputms[0][x].nmhits);
                            else  {
                                if(matrix[x].theta_taj*(data[0][0].time_spec/2.0+1.0)/(double)matrix[x].nsites >= (double).75) {
                                    printf("\nVariation for locus %d is too high or divergence too far. Outgroup inactivated.",x);
                                    if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far. Outgroup inactivated",x);
                                    inputms[0][x].theta = matrix[x].theta_taj;
									data[0][0].time_spec = 0.;
									return 0;
                                }
                                else inputms[0][x].theta = (double)-.75*(double)log((double)1.-(double)4./(double)3.*matrix[x].theta_wat*(data[0][0].time_spec/2.0+(double)1.0)/(double)matrix[x].nsites)
									* ((double)matrix[x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
                            }
                        }
                        y = 1;
                        break;
                    case '2':
                        for(x=0;x<data[0][0].n_loci;x++) {
                            if(matrix[x].noutgroups==0) inputms[0][x].theta = (double)-.75*(double)log((double)1.-(double)4./(double)3.*matrix[x].theta_fuli/(double)matrix[x].nsites) * ((double)matrix[x].nsites+(double)inputms[0][x].nmhits);
                            else  {
                                if(matrix[x].theta_fuli*(data[0][0].time_spec/2.0+1.0)/(double)matrix[x].nsites >= (double).75) {
                                    printf("\nVariation for locus %d is too high or divergence too far. Outgroup inactivated.",x);
                                    if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far. Outgroup inactivated",x);
                                    inputms[0][x].theta = matrix[x].theta_fuli;
									data[0][0].time_spec = 0.;
									return 0;
                                }
                                else inputms[0][x].theta = (double)-.75*(double)log((double)1.-(double)4./(double)3.*matrix[x].theta_fuli*(data[0][0].time_spec/2.0+(double)1.0)/(double)matrix[x].nsites)
									* ((double)matrix[x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
                            }
                        }
                        y = 1;
                        break;
                    case '3':
                        for(x=0;x<data[0][0].n_loci;x++) {
                            if(matrix[x].noutgroups) {
                                if(matrix[x].theta_fw*(data[0][0].time_spec/2.0+1.0)/(double)matrix[x].nsites >= (double).75) {
                                    printf("\nVariation for locus %d is too high or divergence too far",x);
                                    if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
                                    inputms[0][x].theta = matrix[x].theta_fw;
									data[0][0].time_spec = (double)-1.;
									return 0;
                                }
                                else inputms[0][x].theta = (double)-.75*(double)log((double)1.-(double)4./(double)3.*matrix[x].theta_fw*(data[0][0].time_spec/2.0+(double)1.0)/(double)matrix[x].nsites)
									* ((double)matrix[x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
                            }
                        }
                        if(matrix[x].noutgroups) y = 1;
                        else {
                            puts("\nnormalized Fay and Wu's theta estimation is not possible without outgroup.\n");
                            y = 0;
                        }
                        break;
                    case '4':
						/*HKA theta estimations are already corrected by recurrent mutations*/
						if(global[0].outgroup) {
							for(x=0;x<data[0][0].n_loci;x++) {
                                    inputms[0][x].theta = matrix[x].hka_theta;
							}
							y = 1;
                        }
                        else {
                            puts("\nHKA theta estimation is not possible without outgroup.\n");
                            y = 0;
                        }
                        break;
                    case '5':
						/*ML theta estimations:*/
						if(global[0].ml0done == 1 && global[0].dataequalsim == 1) {
							/*menu: take 1,2,3,4, all or the best model of ML thetas*/
							#if COMMAND_LINE
							printf("\n\n     Indicate the ML estimate theta (4Nu) chosen:\n");
							if(file_output) {
								fprintf(file_output,"\n\n     Indicate the ML estimate theta (4Nu) chosen:\n");
							}
							printf("\n");
							
							if(global[0].mltheta == 5) bestm = global[0].nlocimatt;
							else bestm = global[0].mltheta;
							
							printf(" 0 - Best model (recommended): %d different theta/s.\n",bestm);
							printf(" 1 - A single theta for all loci.\n");
							printf(" 2 - Two different thetas\n");
							printf(" 3 - Three different thetas.\n");
							printf(" 4 - Four different thetas\n");
							printf(" 5 - Each locus has a different theta.\n\n");
							#endif
						
							do {
								*kk = getchar();
								if((global[0].nlocimat == 1 && (*kk =!'1' && *kk !='0' && *kk !='5' )) && (*kk>='0' && *kk<'6')) {
									printf("Sorry, only option '0', '1' or '5' are allowed. Try again.\n");
									*kk = 0;
								}
								if((global[0].nlocimat == 2 && (*kk >'2' && *kk !='5' )) && (*kk>='0' && *kk<'6')) {
									printf("Sorry, only option '0', '1', '2' or '5' are allowed. Try again.\n");
									*kk = 0;
								}
								if((global[0].nlocimat == 3 && (*kk >'3' && *kk !='5' )) && (*kk>='0' && *kk<'6')) {
									printf("Sorry, only option '0', '1', '2', '3' or '5' are allowed. Try again.\n");
									*kk = 0;
								}
							} while(*kk<'0' || *kk>'5');
							
							if(file_output) {
								fprintf(file_output,"\n");					
								fprintf(file_output," 0 - Best model (recommended): %d different theta/s.\n",bestm);
								fprintf(file_output," 1 - A single theta for all loci.\n");
								fprintf(file_output," 2 - Two different thetas for all loci\n");
								fprintf(file_output," 3 - Three different thetas for all loci.\n");
								fprintf(file_output," 4 - Four different thetas for all loci.\n");
								fprintf(file_output," 5 - Each locus has a different theta.\n\n");
								fprintf(file_output,"OPTION CHOSEN: %d\n\n",*kk);
							}
							
							switch(*kk) {
								case '0':
									if(global[0].nitermatt==0) global[0].mltheta = 1;/*simplest model*/
									printf("Best model: %d different theta/s.\n\n",bestm);
									if(file_output) fprintf(file_output,"Best model: %d different theta/s.\n\n",global[0].mltheta);
									for(x=0;x<data[0][0].n_loci;x++) {
										inputms[0][x].theta = inputms[0][x].factor_chrn;
										if(matrix[x].noutgroups) {
											if(inputms[0][x].thetaml[global[0].mltheta]*(data[0][0].time_spec/2.0+1.0) >= (double).75) {
												printf("\nVariation for locus %d is too high or divergence too far",x);
												if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
												inputms[0][x].theta *= inputms[0][x].thetaml[global[0].mltheta] * (double)inputms[0][x].nsites;
												data[0][0].time_spec = 0.;
												return 0;
											}
											else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[global[0].mltheta]*(data[0][0].time_spec/2.0+(double)1.0))
												* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
										}
										else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[global[0].mltheta]) * ((double)inputms[0][x].nsites+(double)inputms[0][x].nmhits);
									}
									break;
								case '1':
									for(x=0;x<data[0][0].n_loci;x++) {
										inputms[0][x].theta = inputms[0][x].factor_chrn;
										if(matrix[x].noutgroups) {
											if(inputms[0][x].thetaml[1]*(data[0][0].time_spec/2.0+1.0) >= (double).75) {
												printf("\nVariation for locus %d is too high or divergence too far",x);
												if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
												inputms[0][x].theta *= inputms[0][x].thetaml[1] * (double)inputms[0][x].nsites;
												data[0][0].time_spec = 0.;
												return 0;
											}
											else inputms[0][x].theta *= -(double).75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[1]*(data[0][0].time_spec/2.0+(double)1.0))
												* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
										}
										else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[1]) * ((double)inputms[0][x].nsites+(double)inputms[0][x].nmhits);
									}
									break;
								case '2':
									for(x=0;x<data[0][0].n_loci;x++) {
										inputms[0][x].theta = inputms[0][x].factor_chrn;
										if(matrix[x].noutgroups) {
											if(inputms[0][x].thetaml[2]*(data[0][0].time_spec/2.0+1.0) >= (double).75) {
												printf("\nVariation for locus %d is too high or divergence too far",x);
												if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
												inputms[0][x].theta *= inputms[0][x].thetaml[2] * (double)inputms[0][x].nsites;
												data[0][0].time_spec = 0.;
												return 0;
											}
											else inputms[0][x].theta *= -(double).75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[2]*(data[0][0].time_spec/2.0+(double)1.0))
												* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
										}
										else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[2]) * ((double)inputms[0][x].nsites+(double)inputms[0][x].nmhits);
									}
									break;
								case '3':
									for(x=0;x<data[0][0].n_loci;x++) {
										inputms[0][x].theta = inputms[0][x].factor_chrn;
										if(matrix[x].noutgroups) {
											if(inputms[0][x].thetaml[3]*(data[0][0].time_spec/2.0+1.0) >= (double).75) {
												printf("\nVariation for locus %d is too high or divergence too far",x);
												if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
												inputms[0][x].theta *= inputms[0][x].thetaml[3] * (double)inputms[0][x].nsites;
												data[0][0].time_spec = 0.;
												return 0;
											}
											else inputms[0][x].theta *= -(double).75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[2]*(data[0][0].time_spec/2.0+(double)1.0))
												* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
										}
										else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[3]) * ((double)inputms[0][x].nsites+(double)inputms[0][x].nmhits);
									}
									break;
								case '4':
									for(x=0;x<data[0][0].n_loci;x++) {
										inputms[0][x].theta = inputms[0][x].factor_chrn;
										if(matrix[x].noutgroups) {
											if(inputms[0][x].thetaml[4]*(data[0][0].time_spec/2.0+1.0) >= (double).75) {
												printf("\nVariation for locus %d is too high or divergence too far",x);
												if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
												inputms[0][x].theta *= inputms[0][x].thetaml[4] * (double)inputms[0][x].nsites;
												data[0][0].time_spec = 0.;
												return 0;
											}
											else inputms[0][x].theta *= -(double).75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[4]*(data[0][0].time_spec/2.0+(double)1.0))
												* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
										}
										else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[4]) * ((double)inputms[0][x].nsites+(double)inputms[0][x].nmhits);;
									}
									break;
								case '5':
									for(x=0;x<data[0][0].n_loci;x++) {
										inputms[0][x].theta = inputms[0][x].factor_chrn;
										if(matrix[x].noutgroups) {
											if(inputms[0][x].thetaml[5]*(data[0][0].time_spec/2.0+1.0) >= (double).75) {
												printf("\nVariation for locus %d is too high or divergence too far",x);
												if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
												inputms[0][x].theta *= inputms[0][x].thetaml[5] * (double)inputms[0][x].nsites;
												data[0][0].time_spec = 0.;
												return 0;
											}
											else inputms[0][x].theta *= -(double).75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[5]*(data[0][0].time_spec/2.0+(double)1.0))
												* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
										}
										else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[5]) * ((double)inputms[0][x].nsites+(double)inputms[0][x].nmhits);;
									}
									break;
							}							
							y = 1;
						}
                        else {
                            puts("\nSorry, ML theta estimation is not calculated. Do it from Main menu.\n");
                            y = 0;
                        }
                        break;
                    case '6':
						yy = 0;
						if(data[0][0].time_spec >= (double)0) {
							#if COMMAND_LINE
							printf("\n\n     Indicate the level of variation theta (4Nu) chosen:\n");
							printf("     NOTE: The level of variation can be later corrected by a factor in case using Non-standard neutral model\n");
  							printf("\n     Are you including theta values corrected for recurrent variation (y/n)?\n");
							#endif
						
							do {
								*kk = getchar();
							}while(*kk !='y' && *kk !='Y' && *kk !='n' && *kk !='N');
							if(*kk == 'y' || *kk == 'Y') correctedjc = 1;
							else correctedjc = 0;
						}
						else correctedjc = 1;

						do{
							#if COMMAND_LINE
							printf("\n\n     Indicate the parameter of theta:\n\n");
							
							printf(" 0 - Theta per nucleotide.\n");
							printf(" 1 - Theta per locus.\n");
							#endif
						
							do *kk = getchar();
							while(*kk<'0' || *kk>'1');
							
							switch(*kk) {
								case '0':
									xx = 0;
									do{
										#if COMMAND_LINE
										printf("\n\n     Values of theta per nucleotide:\n\n");
										
										printf(" 0 - Same theta value for all loci.\n");
										printf(" 1 - Different theta values per loci.\n");
										#endif
									
										do *l = getchar();
										while(*l<'0' || *l>'1');
										
										switch(*l) {
											case '0':
												do {
													printf("\nIndicate the value of theta/nt: ");
													scanf(" %lg",&value);
													if(value < (double)0. || value > (double)1.0) puts("Error: the value must be between one and zero.");
												}while(value < (double)0. || value > (double)1.0);
												
												for(x=0;x<data[0][0].n_loci;x++) {
													if(correctedjc == 1) inputms[0][x].theta = value * ((double)matrix[x].nsites+(double)inputms[0][x].nmhits);
													else {
														if(value*(data[0][0].time_spec/2.0+1.0) > (double)0.75) {
															printf("\nVariation for locus %d is too high or divergence too far. Outgroup inactivated.",x);
															if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far. Outgroup inactivated",x);
															return 0;
														} 
														inputms[0][x].theta = -(double).75*(double)log((double)1.-(double)4./(double)3.*value*(data[0][0].time_spec/2.0+(double)1.0))
														* ((double)matrix[x].nsites+(double)matrix[x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
													}
												}
												xx = 1;
												break;
											case '1':
												for(x=0;x<data[0][0].n_loci;x++) {
													do {
														printf("\nIndicate the value of theta/nt for the locus %d: ",x);
														scanf(" %lg",&value);
														if(value < (double)0.) puts("Error: the value must be postive or zero.");
														else {
															if(correctedjc == 1) inputms[0][x].theta = value * ((double)matrix[x].nsites+(double)inputms[0][x].nmhits);
															else {
																if(value*(data[0][0].time_spec/2.0+1.0) > (double)0.75) {
																	printf("\nVariation for locus %d is too high or divergence too far. Outgroup inactivated.",x);
																	if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far. Outgroup inactivated",x);
																	value = (double)0.;
																} 
																inputms[0][x].theta = -(double).75*(double)log((double)1.-(double)4./(double)3.*value*(data[0][0].time_spec/2.0+(double)1.0))
																* ((double)matrix[x].nsites+(double)matrix[x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
															}
														}
													}while(value < (double)0.);
												}
												xx = 1;
												break;
										}
									}while(xx==0);
									yy = 1;
									y = 1;
									break;
								case '1':
									for(x=0;x<data[0][0].n_loci;x++) {
										do {
											printf("\nIndicate the value of the level of variation for the locus %d: ",x);
											scanf(" %lg",&value);
											if(inputms[0][x].theta < (double)0.) puts("Error: the value must be postive or zero.");
											else {
												if(correctedjc == 1) inputms[0][x].theta = value;
												else {
													if(value*(data[0][0].time_spec/2.0+1.0)/(double)matrix[x].nsites > (double)0.75) {
														printf("\nVariation for locus %d is too high or divergence too far. Outgroup inactivated.",x);
														if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far. Outgroup inactivated",x);
														return 0;
													} 
													inputms[0][x].theta = -(double).75*(double)log((double)1.-(double)4./(double)3.*value*(data[0][0].time_spec/2.0+(double)1.0)/(double)matrix[x].nsites)
													* ((double)matrix[x].nsites+(double)matrix[x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
												}
											}
										}while(inputms[0][x].theta < (double)0.);
									}
									yy = 1;
									y = 1;
									break;
							}
						}while(yy==0);
                        break;
                }
            }while(y==0);
            y = 0;
            do{
                #if COMMAND_LINE
                printf("\n\n     Indicate the level of recombination (4Nr) per locus:\n\n");
				printf("     NOTE: The level of recombination can be later corrected by a factor in case using Non-standard neutral model\n");
				
                printf(" 0 - Null recombination for all loci.\n");
				printf(" 1 - Hudson (1987) estimations (Warning: these estimations have high variance).\n");
                printf(" 2 - Other values.\n");
                #endif
            
                do *k = getchar();
                while(*k<'0' || *k>'2');
                
                switch(*k) {
                    case '0':
                        for(x=0;x<data[0][0].n_loci;x++) {
                            inputms[0][x].r = (double)0.;
                        }
                        y = 1;
                        break;
					case '1':
						printf("\nThe recombination value might be too high in some loci.\n");
                        for(x=0;x<data[0][0].n_loci;x++) {
							if(matrix[x].Rvpi > (double)matrix[x].nsites) {
								printf("Recombination parameter in locus #%d is too high. Replaced by %.2f (=nsites).",x,matrix[x].nsites);
								inputms[0][x].r = (double)matrix[x].nsites;
							}
							else inputms[0][x].r = matrix[x].Rvpi * inputms[0][x].factor_chrn;
                        }
                        y = 1;
                        break;
                    case '2':						
						do{
							#if COMMAND_LINE
							printf("\n\n     Indicate the parameter of 4Nr:\n\n");
							
							printf(" 0 - 4Nr per nucleotide.\n");
							printf(" 1 - 4Nr per locus.\n");
							#endif
						
							do *kk = getchar();
							while(*kk<'0' || *kk>'1');
							
							switch(*kk) {
								case '0':
									xx = 0;
									do{
										#if COMMAND_LINE
										printf("\n\n     Values of 4Nr per nucleotide:\n\n");
										
										printf(" 0 - Same 4Nr value for all loci.\n");
										printf(" 1 - Different 4Nr values per loci.\n");
										#endif
									
										do *l = getchar();
										while(*l<'0' || *l>'1');
										
										switch(*l) {
											case '0':
												do {
													printf("\nIndicate the value of 4Nr/nt: ");
													scanf(" %lg",&value);
													if(value < (double)0.) puts("Error: the value must be postive or zero.");
													inputms[0][x].r = value;
												}while(inputms[0][x].r < (double)0.);
												
												for(x=0;x<data[0][0].n_loci;x++) {
													inputms[0][x].r = value *(double)matrix[x].nsites;
												}
												xx = 1;
												break;
											case '1':
												for(x=0;x<data[0][0].n_loci;x++) {
													do {
														printf("\nIndicate the value of 4Nr/nt for the locus %d: ",x);
														scanf(" %lg",&value);
														inputms[0][x].r = value;
														if(inputms[0][x].r < (double)0.) puts("Error: the value must be postive or zero.");
													}while(inputms[0][x].r < (double)0.);
												}
												xx = 1;
												break;
										}
									}while(xx==0);
									yy = 1;
									y = 1;
									break;
								case '1':
									for(x=0;x<data[0][0].n_loci;x++) {
										do {
											printf("\nIndicate the value of the 4Nr for the locus %d: ",x);
											scanf(" %lg",&value);
											if(value < (double)0.) puts("Error, the value must be positive or zero.");
											inputms[0][x].r = value;
										}while(inputms[0][x].r < (double)0.);
									}
									yy = 1;
									y = 1;
									break;
							}
						}while(yy==0);
                        break;
                }
            }while(y==0);

			for(x=0;x<data[0][0].n_loci;x++) if(inputms[0][x].r > 0.) break;
			if(x<data[0][0].n_loci) {
				puts("\n    Is Recombination in males supressed (e.g., Drosophila)? (y/n)");
				do *l = getchar();
				while(*l!='y' && *l!='n' && *l!='Y' && *l!='N');
				if(*l=='y' || *l=='Y') data[0][0].no_rec_males = 1;
				else data[0][0].no_rec_males = 0;
			}

            y = 0;
            do{
                #if COMMAND_LINE
                printf("\n\n     Indicate the ratio of transition/transversion per each locus:\n");
                printf("     (Only affects the variation when they are many multiple hits events (when assuming outgroup))\n\n");
                
                printf(" 0 - All ratios are equal.\n");
                printf(" 1 - Use the observed ratios s/v for each locus (if available).\n");
                printf(" 2 - Other values.\n");
                #endif
            
                do *k = getchar();
                while(*k<'0' || *k>'2');
                
                switch(*k) {
                    case '0':
                        do {
                            printf("\nIndicate the ratio s/v for all loci (s/v = 0.5 is not biased. s/v = %.2f is the average multilocus s/v ratio): ",(double)matrixml[0].Stransitions/(double)matrixml[0].Stransversions);
                            scanf(" %lg",&value);
                            inputms[0][0].ratio_sv = value;
                            if(inputms[0][0].ratio_sv <= (double)0.) puts("Error: the value must be positive.");
                        }while(inputms[0][0].ratio_sv <= (double)0.);
                        for(x=1;x<data[0][0].n_loci;x++) {
                            inputms[0][x].ratio_sv = inputms[0][0].ratio_sv;
                        }
                        y = 1;
                        break;
                    case '1':
                        for(x=0;x<data[0][0].n_loci;x++) {
							if((double)matrix[x].transversions)
								inputms[0][x].ratio_sv = (double)matrix[x].transitions/(double)matrix[x].transversions;
							else {
								do {
									printf("\nloci[%d] has not a valid ratio s/v.\nIndicate the ratio s/v for loci[%d] (s/v = 0.5 is not biased. s/v = %.2f is the average multilocus s/v ratio): ",x,x,(double)matrixml[0].Stransitions/(double)matrixml[0].Stransversions);
									scanf(" %lg",&value);
									inputms[0][x].ratio_sv = value;
									if(inputms[0][x].ratio_sv <= (double)0.) puts("Error: the value must be positive.");
								}while(inputms[0][x].ratio_sv <= (double)0.);
							}
                        }
                        y = 1;
						break;
                    case '2':
                        for(x=0;x<data[0][0].n_loci;x++) {
                            do {
                                printf("\nIndicate the value of the ratio s/v for the locus %d: ",x);
                                scanf(" %lg",&value);
                                if(value <= (double)0.) puts("Error: the value must be positive.");
                            }while(value <= (double)0.);
                            inputms[0][x].ratio_sv = value;
                        }
                        y = 1;
                        break;
                }
            }while(y==0);
            *dataobsequalsim = 1;
        }
    }
	else {
		if((*mm == 'y' || *mm == 'Y') && (*dataindata == 1 || *dataindata == 3)) {
			printf("Do you want to modify the parameters of trans/transv, recombination, theta, etc. (y/n)? ");
			do *ll = getchar();
			while(*ll!='y' && *ll!='n' && *ll!='Y' && *ll!='N');
			if(*ll == 'n' || *ll == 'N') *ll = '\0';
			if(*ll == 'y' || *ll == 'Y') *ll = 'n';
		}
		else *ll = 'n';
	}
	
	/*Not considering observed values (defined or not)*/
	/*include different options: change the dataindata and include the ML theta values*/ 
	/*options: no current data (all new) and no observed data OR use current data */
    if(*ll == 'n' || *ll == 'N') {
		if(*dataindata > 0 && (*mm == 'n' || *mm == 'N')) {
			for(x=1;x<data[0][0].n_loci;x++) {
				free(inputms[0][x].factor_pop);
				free(inputms[0][x].nrec);
				free(inputms[0][x].npast);
				free(inputms[0][x].tpast);
				free(inputms[0][x].freq);
			}
		}
		if(*mm == 'n' || *mm == 'N') {
			do {
				puts("\nSimulations run conditioned to the level of variation for each locus.\n");
				puts("\nIndicate the number of loci:\n");
				scanf(" %d",&valuei);
				data[0][0].n_loci = valuei;
				if(data[0][0].n_loci <= 0) puts("Error. The value must be positive.");
			}while(data[0][0].n_loci <= 0);

			if((inputms[0] = (struct var2b *)realloc(inputms[0],data[0][0].n_loci*sizeof(struct var2b))) == 0) {
				puts("\nError: memory not reallocated. main.4 \n");
				exit(0);
			}

			x=0;
			if((inputms[0][x].nrec = (double *) realloc(inputms[0][x].nrec,2*sizeof(double))) == 0) {
				puts("\nError: memory not reallocated. main.9 \n");
				return(0);
			}
			memset(inputms[0][x].nrec,'\0',2*sizeof(double));
			if((inputms[0][x].npast = (double *) realloc(inputms[0][x].npast,2*sizeof(double))) == 0) {
				puts("\nError: memory not reallocated. main.10 \n");
				return(0);
			}
			memset(inputms[0][x].npast,'\0',2*sizeof(double));
			if((inputms[0][x].tpast = (double *) realloc(inputms[0][x].tpast,2*sizeof(double))) == 0) {
				puts("\nError: memory not reallocated. main.11 \n");
				return(0);
			}
			memset(inputms[0][x].tpast,'\0',2*sizeof(double));
			if((inputms[0][x].factor_pop = (double *) realloc(inputms[0][x].factor_pop,1*sizeof(double))) == 0) {
				puts("\nError: memory not reallocated. main.8 \n");
				return(0);
			}
			memset(inputms[0][x].factor_pop,'\0',1*sizeof(double));
			if((inputms[0][x].freq = (double *) realloc(inputms[0][x].freq,1*sizeof(double))) == 0) {
				puts("\nError: memory not reallocated. main.12 \n");
				return(0);
			}
			memset(inputms[0][x].freq,'\0',1*sizeof(double));
			
			for(x=1;x<data[0][0].n_loci;x++) {
				if((inputms[0][x].nrec = (double *) calloc(2,sizeof(double))) == 0) {
					puts("\nError: memory not reallocated. main.9 \n");
					if(x>1) data[0][0].n_loci = x-1;
					else data[0][0].n_loci = 1;
					if((inputms[0] = (struct var2b *)realloc(inputms[0],data[0][0].n_loci*sizeof(struct var2b))) == 0) {
						puts("\nError: memory not reallocated. seqinfodata.4 \n");
						return(0);
					}			
					return(0);
				}
				if((inputms[0][x].npast = (double *) calloc(2,sizeof(double))) == 0) {
					puts("\nError: memory not reallocated. main.10 \n");
					free(inputms[0][x].nrec);
					if(x>1) data[0][0].n_loci = x-1;
					else data[0][0].n_loci = 1;
					if((inputms[0] = (struct var2b *)realloc(inputms[0],data[0][0].n_loci*sizeof(struct var2b))) == 0) {
						puts("\nError: memory not reallocated. seqinfodata.4 \n");
						return(0);
					}			
					return(0);
				}
				if((inputms[0][x].tpast = (double *) calloc(2,sizeof(double))) == 0) {
					puts("\nError: memory not reallocated. main.11 \n");
					free(inputms[0][x].npast);
					free(inputms[0][x].nrec);
					if(x>1) data[0][0].n_loci = x-1;
					else data[0][0].n_loci = 1;
					if((inputms[0] = (struct var2b *)realloc(inputms[0],data[0][0].n_loci*sizeof(struct var2b))) == 0) {
						puts("\nError: memory not reallocated. seqinfodata.4 \n");
						return(0);
					}			
					return(0);
				}
				if((inputms[0][x].freq = (double *) calloc(1,sizeof(double))) == 0) {
					puts("\nError: memory not reallocated. main.12 \n");
					free(inputms[0][x].tpast);
					free(inputms[0][x].npast);
					free(inputms[0][x].nrec);
					if(x>1) data[0][0].n_loci = x-1;
					else data[0][0].n_loci = 1;
					if((inputms[0] = (struct var2b *)realloc(inputms[0],data[0][0].n_loci*sizeof(struct var2b))) == 0) {
						puts("\nError: memory not reallocated. seqinfodata.4 \n");
						return(0);
					}			
					return(0);
				}
				if((inputms[0][x].factor_pop = (double *) calloc(1,sizeof(double))) == 0) {
					puts("\nError: memory not reallocated. main.8 \n");
					free(inputms[0][x].freq);
					free(inputms[0][x].tpast);
					free(inputms[0][x].npast);
					free(inputms[0][x].nrec);
					if(x>1) data[0][0].n_loci = x-1;
					else data[0][0].n_loci = 1;
					if((inputms[0] = (struct var2b *)realloc(inputms[0],data[0][0].n_loci*sizeof(struct var2b))) == 0) {
						puts("\nError: memory not reallocated. seqinfodata.4 \n");
						return(0);
					}			
					return(0);
				}
			}
			
            for(x=0;x<data[0][0].n_loci;x++) {
				inputms[0][x].npop = 1;
				inputms[0][x].nintn = 0;
                inputms[0][x].factor_chrn = matrix[x].factor_chrn;
                inputms[0][x].nsam = matrix[x].nsamples;
                inputms[0][x].nsites = (long int)matrix[x].nsites;
				inputms[0][x].nmhits = 0;
            }

			puts("\nIndicate the correction factor for chromosome population size:\n");
			for(x=0;x<data[0][0].n_loci;x++) {
				do{
					printf("\n   Indicate the correction factor for chromosome population size for the locus %d: ",x);
					scanf(" %lf",&valuef);
					inputms[0][x].factor_chrn = valuef;
					if(inputms[0][x].factor_chrn <= (double)0) puts("Error. The value must be higher than (double)0.");
				}while(inputms[0][x].factor_chrn <= (double)0);
			}
			puts("\nIndicate the number of samples:\n");
			for(x=0;x<data[0][0].n_loci;x++) {
				do{
					printf("\n   Indicate the number of samples for the locus %d: ",x);
					scanf(" %d",&valuei);
					inputms[0][x].nsam = valuei;
					if(inputms[0][x].nsam <= 1) puts("Error. The value must be higher than 1.");
				}while(inputms[0][x].nsam <= 1);
			}
			puts("\nIndicate the number of sites (including mhits):\n");
			for(x=0;x<data[0][0].n_loci;x++) {
				do {
					printf("\n   Indicate the number of sites for the locus %d: ",x);
					scanf(" %ld",&valueli);
					inputms[0][x].nsites = valueli;
					if(inputms[0][x].nsites <= 0) puts("Error. The value must be positive.");
				}while(inputms[0][x].nsites <= 0);
				inputms[0][x].nmhits = 0;
			}
		}
        puts("Do you want to include an outgroup (y/n)? ");
        do *l = getchar();
        while(*l!='y' && *l!='n' && *l!='Y' && *l!='N');
        if(*l == 'y' || *l == 'Y') {
            do {
                puts("\nTime (in 2N generations) of divergence between the ingroup and the ancestral population? ");
                scanf(" %lg",&valuef);
                data[0][0].time_spec = valuef;
                if(data[0][0].time_spec < (double)0.) puts("Error: the value must be positive.");
            }while(data[0][0].time_spec < (double)0.);
        }
        else data[0][0].time_spec = 0.;/*-1.0*/
		
		y = 0;
		while(y==0) {
			#if COMMAND_LINE
			printf("\nIndicate the level of variation per locus:\n");
			if(file_output) {
				fprintf(file_output,"\nIndicate the level of variation per locus:\n");
			}
			printf("\n");
			
			printf(" 0 - Introduce manually the level of variation for each locus.\n");
			printf(" 1 - Maximum likelihood estimate (if previosly calculated).\n");
			#endif

			do *k = getchar();
			while(*k<'0' || *k>'1');
			
			switch(*k) {
				case '0':
					/*Indicate manually the level of variation*/
					yy = 0;
					if(data[0][0].time_spec >= (double)0) {
						#if COMMAND_LINE
						printf("\n\n     Indicate the parameter of theta:\n");					
						printf("\n     Are you including theta values corrected for recurrent variation (y/n)?\n");
						#endif
					
						do {
							*kk = getchar();
						}while(*kk !='y' && *kk !='Y' && *kk !='n' && *kk !='N');
						if(*kk == 'y' || *kk == 'Y') correctedjc = 1;
						else correctedjc = 0;
					}
					else correctedjc = 1;

					do{
						#if COMMAND_LINE
						printf("\n\n     Indicate the level of variation theta (4Nu) chosen:\n");
						printf("     NOTE: The level of variation can be later corrected by a factor in case using Non-standard neutral model\n");
  						
						printf(" 0 - Theta per nucleotide.\n");
						printf(" 1 - Theta per locus.\n");
						#endif
					
						do *kk = getchar();
						while(*kk<'0' || *kk>'1');
						
						switch(*kk) {
							case '0':
								xx = 0;
								do{
									#if COMMAND_LINE
									printf("\n\n     Values of theta per nucleotide:\n\n");
									
									printf(" 0 - Same theta value for all loci.\n");
									printf(" 1 - Different theta values per loci.\n");
									#endif
								
									do *l = getchar();
									while(*l<'0' || *l>'1');
									
									switch(*l) {
										case '0':
											do {
												printf("\nIndicate the value of theta/nt: ");
												scanf(" %lg",&value);
												if(value < (double)0. || value > (double)1.0) puts("Error: the value must be between one and zero.");
											}while(value < (double)0. || value > (double)1.0);
											
											for(x=0;x<data[0][0].n_loci;x++) {
												if(correctedjc == 1) inputms[0][x].theta = value * ((double)matrix[x].nsites+(double)inputms[0][x].nmhits);
												else {
													if(value*(data[0][0].time_spec/2.0+1.0) > (double)0.75) {
														printf("\nVariation for locus %d is too high or divergence too far. Outgroup inactivated.",x);
														if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far. Outgroup inactivated",x);
														return 0;
													} 
													inputms[0][x].theta = -(double).75*(double)log((double)1.-(double)4./(double)3.*value*(data[0][0].time_spec/2.0+(double)1.0))
													* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
												}
											}
											xx = 1;
											break;
										case '1':
											for(x=0;x<data[0][0].n_loci;x++) {
												do {
													printf("\nIndicate the value of theta/nt for the locus %d: ",x);
													scanf(" %lg",&value);
													if(value < (double)0.) puts("Error: the value must be postive or zero.");
													else {
														if(correctedjc == 1) inputms[0][x].theta = value * ((double)matrix[x].nsites+(double)inputms[0][x].nmhits);
														else {
															if(value*(data[0][0].time_spec/2.0+1.0) > (double)0.75) {
																printf("\nVariation for locus %d is too high or divergence too far. Outgroup inactivated.",x);
																if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far. Outgroup inactivated",x);
																value = (double)0.;
															} 
															inputms[0][x].theta = -(double).75*(double)log((double)1.-(double)4./(double)3.*value*(data[0][0].time_spec/2.0+(double)1.0))
															* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
														}
													}
												}while(value < (double)0.);
											}
											xx = 1;
											break;
									}
								}while(xx==0);
								yy = 1;
								y = 1;
								break;
							case '1':
								for(x=0;x<data[0][0].n_loci;x++) {
									do {
										printf("\nIndicate the value of the level of variation for the locus %d: ",x);
										scanf(" %lg",&value);
										if(inputms[0][x].theta < (double)0.) puts("Error: the value must be postive or zero.");
										else {
											if(correctedjc == 1) inputms[0][x].theta = value;
											else {
												if(value*(data[0][0].time_spec/2.0+1.0)/(double)matrix[x].nsites > (double)0.75) {
													printf("\nVariation for locus %d is too high or divergence too far. Outgroup inactivated.",x);
													if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far. Outgroup inactivated",x);
													return 0;
												} 
												inputms[0][x].theta = -(double).75*(double)log((double)1.-(double)4./(double)3.*value*(data[0][0].time_spec/2.0+(double)1.0)/(double)inputms[0][x].nsites)
												* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
											}
										}
									}while(inputms[0][x].theta < (double)0.);
								}
								yy = 1;
								y = 1;
								break;
						}
					}while(yy==0);
					break;
				case '1':
					/*MLtheta*/
					if(*dataindata == 2 && global[0].ml0done == 1) {
						/*menu: take 1,2,3,4, all or the best model of ML thetas*/
						#if COMMAND_LINE
						printf("\n\n     Indicate the ML estimate theta (4Nu) chosen:\n");
						if(file_output) {
							fprintf(file_output,"\n\n     Indicate the ML estimate theta (4Nu) chosen:\n");
						}
						printf("\n");
						
						if(global[0].mltheta == 5) bestm = global[0].nlocimatt;
						else bestm = global[0].mltheta;
							
						printf(" 0 - Best model (recommended): %d different theta/s.\n",bestm);
						printf(" 1 - A single theta for all loci.\n");
						printf(" 2 - Two different thetas\n");
						printf(" 3 - Three different thetas.\n");
						printf(" 4 - Four different thetas\n");
						printf(" 5 - Each locus has a different theta.\n\n");
						#endif
					
						if(file_output) {
							fprintf(file_output,"\n");					
							fprintf(file_output," 0 - Best model (recommended): %d different theta/s.\n",bestm);
							fprintf(file_output," 1 - A single theta for all loci.\n");
							fprintf(file_output," 2 - Two different thetas for all loci\n");
							fprintf(file_output," 3 - Three different thetas for all loci.\n");
							fprintf(file_output," 4 - Four different thetas for all loci.\n");
							fprintf(file_output," 5 - Each locus has a different theta.\n\n");
							fprintf(file_output,"OPTION CHOSEN: %d\n\n",*kk);
						}

						do {
							*kk = getchar();
							if((global[0].nlocimat == 1 && (*kk =!'1' && *kk !='0' && *kk !='5' )) && (*kk>='0' && *kk<'6')) {
								printf("Sorry, only option '0', '1' or '5' are allowed. Try again.\n");
								*kk = 0;
							}
							if((global[0].nlocimat == 2 && (*kk >'2' && *kk !='5' )) && (*kk>='0' && *kk<'6')) {
								printf("Sorry, only option '0', '1', '2' or '5' are allowed. Try again.\n");
								*kk = 0;
							}
							if((global[0].nlocimat == 3 && (*kk >'3' && *kk !='5' )) && (*kk>='0' && *kk<'6')) {
								printf("Sorry, only option '0', '1', '2', '3' or '5' are allowed. Try again.\n");
								*kk = 0;
							}
						} while(*kk<'0' || *kk>'5');
						
						switch(*kk) {
							case '0':
								if(global[0].nitermatt==0) global[0].mltheta = 1;/*simplest model*/
								printf("Best model: %d different theta/s.\n\n",bestm);
								if(file_output) fprintf(file_output,"Best model: %d different theta/s.\n\n",global[0].mltheta);
								for(x=0;x<data[0][0].n_loci;x++) {
									inputms[0][x].theta = inputms[0][x].factor_chrn;
									if(data[0][0].time_spec > (double)-1) {
										if(inputms[0][x].thetaml[global[0].mltheta]*(data[0][0].time_spec/2.0+(double)1.0) >= (double).75) {
											printf("\nVariation for locus %d is too high or divergence too far",x);
											if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
											inputms[0][x].theta *= inputms[0][x].thetaml[global[0].mltheta] * (double)inputms[0][x].nsites;
											data[0][0].time_spec = 0.;
											return 0;
										}
										else inputms[0][x].theta *= -(double).75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[global[0].mltheta]*(data[0][0].time_spec/2.0+(double)1.0))
											* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
									}
									else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[global[0].mltheta]) * ((double)inputms[0][x].nsites+(double)inputms[0][x].nmhits);
								}
								break;
							case '1':
								for(x=0;x<data[0][0].n_loci;x++) {
									inputms[0][x].theta = inputms[0][x].factor_chrn;
									if(data[0][0].time_spec > (double)-1) {
										if(inputms[0][x].thetaml[1]*(data[0][0].time_spec/2.0+(double)1.0) >= (double).75) {
											printf("\nVariation for locus %d is too high or divergence too far",x);
											if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
											inputms[0][x].theta *= inputms[0][x].thetaml[1] * (double)inputms[0][x].nsites;
											data[0][0].time_spec = 0.;
											return 0;
										}
										else inputms[0][x].theta *= -(double).75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[1]*(data[0][0].time_spec/2.0+(double)1.0))
											* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
									}
									else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[1]) * ((double)inputms[0][x].nsites+(double)inputms[0][x].nmhits);
								}
								break;
							case '2':
								for(x=0;x<data[0][0].n_loci;x++) {
									inputms[0][x].theta = inputms[0][x].factor_chrn;
									if(data[0][0].time_spec > (double)-1) {
										if(inputms[0][x].thetaml[2]*(data[0][0].time_spec/2.0+1.0) >= (double).75) {
											printf("\nVariation for locus %d is too high or divergence too far",x);
											if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
											inputms[0][x].theta *= inputms[0][x].thetaml[2] * (double)inputms[0][x].nsites;
											data[0][0].time_spec = 0.;
											return 0;
										}
										else inputms[0][x].theta *= -(double).75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[2]*(data[0][0].time_spec/2.0+(double)1.0))
											* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
									}
									else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[2]) * ((double)inputms[0][x].nsites+(double)inputms[0][x].nmhits);
								}
								break;
							case '3':
								for(x=0;x<data[0][0].n_loci;x++) {
									inputms[0][x].theta = inputms[0][x].factor_chrn;
									if(data[0][0].time_spec > (double)-1) {
										if(inputms[0][x].thetaml[3]*(data[0][0].time_spec/2.0+1.0) >= (double).75) {
											printf("\nVariation for locus %d is too high or divergence too far",x);
											if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
											inputms[0][x].theta *= inputms[0][x].thetaml[3] * (double)inputms[0][x].nsites;
											data[0][0].time_spec = 0.;
											return 0;
										}
										else inputms[0][x].theta *= -(double).75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[2]*(data[0][0].time_spec/2.0+(double)1.0))
											* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
									}
									else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[3]) * ((double)inputms[0][x].nsites+(double)inputms[0][x].nmhits);
								}
								break;
							case '4':
								for(x=0;x<data[0][0].n_loci;x++) {
									inputms[0][x].theta = inputms[0][x].factor_chrn;
									if(data[0][0].time_spec > (double)-1) {
										if(inputms[0][x].thetaml[4]*(data[0][0].time_spec/2.0+1.0) >= (double).75) {
											printf("\nVariation for locus %d is too high or divergence too far",x);
											if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
											inputms[0][x].theta *= inputms[0][x].thetaml[4] * (double)inputms[0][x].nsites;
											data[0][0].time_spec = 0.;
											return 0;
										}
										else inputms[0][x].theta *= -(double).75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[4]*(data[0][0].time_spec/2.0+(double)1.0))
											* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
									}
									else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[4]) * ((double)inputms[0][x].nsites+(double)inputms[0][x].nmhits);
								}
								break;
							case '5':
								for(x=0;x<data[0][0].n_loci;x++) {
									inputms[0][x].theta = inputms[0][x].factor_chrn;
									if(data[0][0].time_spec > (double)-1) {
										if(inputms[0][x].thetaml[5]*(data[0][0].time_spec/2.0+1.0) >= (double).75) {
											printf("\nVariation for locus %d is too high or divergence too far",x);
											if(file_output) fprintf(file_output,"\nVariation for locus %d is too high or divergence too far",x);
											inputms[0][x].theta *= inputms[0][x].thetaml[5] * (double)inputms[0][x].nsites;
											data[0][0].time_spec = 0.;
											return 0;
										}
										else inputms[0][x].theta *= -(double).75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[5]*(data[0][0].time_spec/2.0+(double)1.0))
											* ((double) inputms[0][x].nsites+(double)inputms[0][x].nmhits)/(data[0][0].time_spec/2.0+(double)1.0);
									}
									else inputms[0][x].theta *= (double)-.75*(double)log((double)1.-(double)4./(double)3.*inputms[0][x].thetaml[5]) * ((double)inputms[0][x].nsites+(double)inputms[0][x].nmhits);
								}
								break;
						}							
						y = 1;
					}
					else {
						puts("\nSorry, ML theta estimation is not calculated. Do it from Main menu.\n");
						y = 0;
					}
					break;
			}
		}
		/*recombination*/
		y=0;
		do{
			#if COMMAND_LINE
			printf("\n\n     Indicate the level of recombination of 4Nr:\n\n");
			printf("     NOTE: The level of recombination can be later corrected by a factor in case using Non-standard neutral model\n");
  			
			printf(" 0 - 4Nr per nucleotide.\n");
			printf(" 1 - 4Nr per locus.\n");
			#endif
		
			do *kk = getchar();
			while(*kk<'0' || *kk>'1');
			
			switch(*kk) {
				case '0':
					xx = 0;
					do{
						#if COMMAND_LINE
						printf("\n\n     Values of 4Nr per nucleotide:\n\n");
						
						printf(" 0 - Same 4Nr value for all loci.\n");
						printf(" 1 - Different 4Nr values per loci.\n");
						#endif
					
						do *l = getchar();
						while(*l<'0' || *l>'1');
						
						switch(*l) {
							case '0':
								do {
									printf("\nIndicate the value of 4Nr/nt: ");
									scanf(" %lg",&value);
									if(value < (double)0.) puts("Error: the value must be postive or zero.");
								}while(value < (double)0.);
								
								for(x=0;x<data[0][0].n_loci;x++) {
									inputms[0][x].r = value *(double)inputms[0][x].nsites;
								}
								xx = 1;
								break;
							case '1':
								for(x=0;x<data[0][0].n_loci;x++) {
									do {
										printf("\nIndicate the value of 4Nr/nt for the locus %d: ",x);
										scanf(" %lg",&value);
										inputms[0][x].r = value;
										if(inputms[0][x].r < (double)0.) puts("Error: the value must be postive or zero.");
									}while(inputms[0][x].r < (double)0.);
								}
								xx = 1;
								break;
						}
					}while(xx==0);
					y = 1;
					break;
				case '1':
					for(x=0;x<data[0][0].n_loci;x++) {
						do {
							printf("\nIndicate the value of the 4Nr for the locus %d: ",x);
							scanf(" %lg",&value);
							if(value < (double)0.) puts("Error, the value must be positive or zero.");
							inputms[0][x].r = value;
						}while(inputms[0][x].r < (double)0.);
					}
					y = 1;
					break;
			}
		}while(y==0);

		for(x=0;x<data[0][0].n_loci;x++) if(inputms[0][x].r > 0.) break;
		if(x<data[0][0].n_loci) {
			puts("\n    Is Recombination in males supressed (e.g., Drosophila)? (y/n)");
			do *l = getchar();
			while(*l!='y' && *l!='n' && *l!='Y' && *l!='N');
			if(*l=='y' || *l=='Y') data[0][0].no_rec_males = 1;
			else data[0][0].no_rec_males = 0;
		}
        y = 0;
        do{
            #if COMMAND_LINE
            printf("\n\n     Indicate the ratio of transition/transversion per each locus:\n");
            printf("     (Only affects the variation when they are many multiple hits events (when assuming outgroup))\n\n");
            
            printf(" 0 - All ratios are equal.\n");
            printf(" 1 - Other values.\n");
            #endif
        
            do *k = getchar();
            while(*k<'0' || *k>'1');
            
            switch(*k) {
                case '0':
                    do {
                        puts("\nIndicate the ratio s/v for all loci (e.g., 0.5 is not biased): ");
                        scanf(" %lg",&inputms[0][0].ratio_sv);
                        if(inputms[0][0].ratio_sv <= (double)0.) puts("Error: the value must be positive.");
                    }while(inputms[0][0].ratio_sv <= (double)0.);
                    for(x=1;x<data[0][0].n_loci;x++) {
                        inputms[0][x].ratio_sv = inputms[0][0].ratio_sv;
                    }
                    y = 1;
                    break;
                case '1':
                    for(x=0;x<data[0][0].n_loci;x++) {
                        do {
                            printf("\nIndicate the value of the ratio s/v for the locus %d: ",x);
                            scanf(" %lg",&value);
                            if(value <= (double)0.) puts("Error: the value must be positive.");
                        }while(value <= (double)0.);
                        inputms[0][x].ratio_sv = value;
                    }
                    y = 1;
                    break;
            }
        }while(y==0);
        *dataobsequalsim = 0;
    }
    if(*dataindata == 0) *dataindata = 1;
	else if(*dataindata == 2) *dataindata = 3;
    return 1;
}

int open_menu_evolcond_sim(struct var **data,struct var2b **inputms,int *dataforsim/*,struct globvar *global*/,struct statistics *matrix,FILE *file_output)
{
    int x,y;
	int n;
    char k[1],l[1],m[1],kk[1];
    double sum,g;
    int valuei;
    double value,value2,ts;
    long int valueli;
	int nloci,remloci,numberloci;
	double factorth;
	
	while(1) {
		printf("\n\n");
		printf("**********************************************************\n");
		printf("**********************************************************\n");
		printf("* WARNING: Selective processes affect specific loci.     *\n");
		printf("* A demographic process should affect the entire genome. *\n");
		printf("*  The use of different demographic processes            *\n");
		printf("*  in different loci could not be justified.             *\n");
		printf("*  A selective process is modeled under                  *\n");
		printf("*  a stationary population background.                   *\n");
		printf("**********************************************************\n");
		printf("**********************************************************\n");
		printf("\n  Introduce the number of evolutionary processes \nconsidered in coalescent simulations:");
		printf("\n  (e.g., if one locus is under selection, we consider two events:\n selection in one locus and neutrality in the rest of loci).\n  ");
		do {
			scanf(" %d",&n);
			if(n <= 0 || n > data[0][0].n_loci) printf("Number of processes must between 1 and %d\n\n Introduce again the number of processes:\n",data[0][0].n_loci);
		} while(n <= 0 || n > data[0][0].n_loci);
		data[0][0].n_events = n;
		printf("\n\nNumber of evolutionary processes considered in coalescent simulations: %d",n);
		if(file_output) fprintf(file_output,"\n\nNumber of evolutionary processes considered in coalescent simulations: %d",n);
		data[0][0].neutrality = 0;
		data[0][0].subdivision = 0;
		data[0][0].changesizepop = 0;
		data[0][0].split_pop = 0;
		data[0][0].pselection = 0;
		
		remloci = data[0][0].n_loci;
		for(x=0;x<data[0][0].n_loci;x++) {
			inputms[0][x].definedmodel = 0;
			inputms[0][x].neutrality = 0;
			inputms[0][x].subdivision = 0;
			inputms[0][x].changesizepop = 0;
			inputms[0][x].split_pop = 0;
			inputms[0][x].pselection = 0;
		}
		
		while(n--) {
			#if COMMAND_LINE
			printf("\n\nEVOLUTIONARY PROCESS: %d\n",data[0][0].n_events - n);
			 
			printf("\n  Introduce parameters related to evolutionary conditions of simulations:\n\n");
			
			printf("   CHOOSE:\n");
			printf("    0 - Neutral panmictic model with constant population size.\n");
			printf("    1 - Growth/Decline/Bottleneck population size model (using a logistic equation).\n");
			printf("    2 - Subdivision (Wright island model).\n");
			printf("    3 - Population split and merge model.\n");
			printf("    4 - Strong positive selection in one or several loci.\n    ");
			#endif
		
			do *k = getchar();
			while(*k<'0' || *k>'4');
			
			if(*k=='0') data[0][0].neutrality += 1;
			if(*k=='1') data[0][0].changesizepop += 1;
			if(*k=='2') data[0][0].subdivision += 1;
			if(*k=='3') data[0][0].split_pop += 1;
			if(*k=='4') data[0][0].pselection += 1;

			if(n) {
				#if COMMAND_LINE
				printf("\n   Indicate the number of loci that follow this model\n\n");
				printf("   CHOOSE:\n");
				printf("    0 - All the remaining loci not yet assigned to any model.\n");
				printf("    1 - A group of loci or a single locus.\n");
				#endif
				
				do {
					*m = getchar();
					if(*m=='0' && n>0) {
						printf("\n    More processes are defined. Please indicate a group of loci.\n");
						*m = '1';
					}
				}while(*m<'0' || *m>'1');
			}
			else *m = 0;
			
			if(*m == '1' && n>0) {
				do{
					printf("\n    Indicate the number of loci: ");
					scanf(" %d",&nloci);
					if(nloci > remloci - n || nloci < 0) printf("\n    The number of loci must be between %d and %d.\n",0,remloci-n);
				}while(nloci > remloci - n || nloci < 0);
				for(y=0;y<nloci;y++) {
					do{
						printf("\n    Locus %d: Indicate the ID number of the locus [0,%d]: ",y,data[0][0].n_loci-1);
						scanf(" %d",&numberloci);
						printf("    %d:%s\n",numberloci,matrix[numberloci].gene);
						if(inputms[0][numberloci].definedmodel != 0 || numberloci < 0 || numberloci > data[0][0].n_loci-1)
							printf("\n    Sorry, this locus has a defined process or the number is out of range.");
					}while(inputms[0][numberloci].definedmodel != 0 || numberloci < 0 || numberloci > data[0][0].n_loci-1); 
					inputms[0][numberloci].definedmodel = data[0][0].n_events - n;
					if(*k=='0') inputms[0][numberloci].neutrality = data[0][0].neutrality;
					if(*k=='1') inputms[0][numberloci].changesizepop = data[0][0].changesizepop;
					if(*k=='2') inputms[0][numberloci].subdivision = data[0][0].subdivision;
					if(*k=='3') inputms[0][numberloci].split_pop = data[0][0].split_pop;
					if(*k=='4') inputms[0][numberloci].pselection = data[0][0].pselection;
				}
			}
			else {
				for(y=0;y<data[0][0].n_loci;y++) {
					if(inputms[0][y].definedmodel == 0) {
						inputms[0][y].definedmodel = data[0][0].n_events - n;
						if(*k=='0') inputms[0][y].neutrality = data[0][0].neutrality;
						if(*k=='1') inputms[0][y].changesizepop = data[0][0].changesizepop;
						if(*k=='2') inputms[0][y].subdivision = data[0][0].subdivision;
						if(*k=='3') inputms[0][y].split_pop = data[0][0].split_pop;
						if(*k=='4') inputms[0][y].pselection = data[0][0].pselection;
					}
				}
				printf("\n    All the remaining loci are included in the last process.\n");
			}
			/*if(*k != '0') {*//*correction of theta values by a factor in case non strictly neutral model...*/
				printf("\n     Do you want to use a factor correction of the level of variation for these loci (y/n)? ");
				do {
					*kk = getchar();
				}while(*kk !='y' && *kk !='Y' && *kk !='n' && *kk !='N');
				if(*kk =='Y' || *kk =='y') {
					do {
						printf("\n     Please indicate the factor of correction: ");
						scanf(" %lf",&factorth);
						if(factorth < (double)0) printf("\n     Sorry, the value must be positive.");
					}while(factorth < (double)0);

					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].theta = inputms[0][y].theta * factorth;
							inputms[0][y].r = inputms[0][y].r * factorth;
						}
					}
				}
			/*}*/
			
			switch(*k) {
				case '0':
					/*neutral panmictic*/
					/*
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							if((inputms[0][y].nrec = (double *)realloc(inputms[0][y].nrec,((unsigned)2)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].nrec[0] = (double)0;
							inputms[0][y].nrec[1] = (double)0;

							if((inputms[0][y].npast = (double *)realloc(inputms[0][y].npast,((unsigned)2)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].npast[0] = (double)0;
							inputms[0][y].npast[1] = (double)0;

							if((inputms[0][y].tpast = (double *)realloc(inputms[0][y].tpast,((unsigned)2)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].tpast[0] = (double)0;
							inputms[0][y].tpast[1] = (double)0;
							
							if((inputms[0][y].factor_pop = (double *)realloc(inputms[0][y].factor_pop,((unsigned)1)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].factor_pop[0] = (double)1.;
							
							if((inputms[0][y].freq = (double *)realloc(inputms[0][y].freq,((unsigned)1)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].freq[0] = (double)1;
						}
					}
					*/
					printf("\n   Neutral panmictic model with constant population size assigned.\n");
					break;
				case '1': 
					/*change in pop size*/
					printf("\n   All the events are considered here from the present to the past.");
					printf("\n   The population size is counted in relation to N initial (i.e., No = 1).\n");
					printf("\n   Note: Logistic equation is only used within events (from start to end). Between events the change is instantaneous.\n");
					do {
						printf("\n   Number of events? (e.g., 3. First event constant pop, second expansion, third reduction.) ");
						scanf(" %d",&valuei);
						if(valuei <= 0)
							puts("Error: the number of events must be one or more.\n");
					} while(valuei <= 0);
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].nintn = valuei;
						}
					}
					
					/*realloc nrec,npast,tpast*/
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							if((inputms[0][y].nrec = (double *)realloc(inputms[0][y].nrec,(2+inputms[0][y].nintn)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							memset(inputms[0][y].nrec,'\0',(2+inputms[0][y].nintn)*sizeof(double));

							if((inputms[0][y].npast = (double *)realloc(inputms[0][y].npast,(2+inputms[0][y].nintn)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							memset(inputms[0][y].npast,'\0',(2+inputms[0][y].nintn)*sizeof(double));

							if((inputms[0][y].tpast = (double *)realloc(inputms[0][y].tpast,(2+inputms[0][y].nintn)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}							
							memset(inputms[0][y].tpast,'\0',(2+inputms[0][y].nintn)*sizeof(double));

							for(x=0;x<inputms[0][y].nintn;x++) {
								if(x==0) {
									inputms[0][y].nrec[x] = (double)1.0;
									/*printf("\nThe relative value of N at the beginning of the first event is 1.0\n");*/
								}
							}
						}
					}
					printf("\n   The relative value of N at the beginning of the first event is 1.0\n");
					
					for(x=0;x<valuei;x++) {
						if(x>0) {
							do {
								printf("\n   Relative Population size N at the beggining of the event %d (less than 1 means expansion)? ",x);
								scanf(" %lg",&value);
								if(value <= (double)0.) 
									puts("Error: the value must be positive.\n");
							}while(value <= (double)0.);
							for(y=0;y<data[0][0].n_loci;y++) {
								if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
									inputms[0][y].nrec[x] = value;
								}
							}
						}
						do {
							printf("\n   Relative Population size N at the end of the event %d (less than 1 means expansion)? ",x);
							scanf(" %lg",&value);
							if(value <= (double)0.) 
								puts("Error: the value must be positive.\n");
						}while(value <= (double)0.);
						for(y=0;y<data[0][0].n_loci;y++) {
							if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
								inputms[0][y].npast[x] = value;
							}
						}
						
						do {
							printf("\n   Duration (in 4No generations) of the event %d? ",x);
							scanf(" %lg",&value);
							if(value <= (double)0.) 
								puts("Error: the value must be positive.\n");
						}while(value <= (double)0.);
						ts = value;
						for(y=0;y<data[0][0].n_loci;y++) {
							if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
								inputms[0][y].tpast[x] = value;
							}
						}

						if(x==0) {
							do {
								printf("\n   ts indicates that the curve will be logistic if ts is 0,");
								printf("\n   if the duration of ts is half of the first event, the curve will be exponential-like,");
								printf("\n   if ts is almost the same than the event, the curve wil be lineal-like.\n");
								printf("\n   Value of ts at the FIRST event? ");
								scanf(" %lg",&value);
								if(value < (double)0 || value > ts) 
									printf("Error: the value must be between 0 and %f.\n",ts);
							}while(value < (double)0 || value > ts);
							for(y=0;y<data[0][0].n_loci;y++) {
								if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
									inputms[0][y].ts = value;
								}
							}
						}
						/*
						for(y=0;y<data[0][0].n_loci;y++) {
							if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
								if((inputms[0][y].factor_pop = (double *)realloc(inputms[0][y].factor_pop,(1)*sizeof(double))) == 0) {
									puts(" NOT SUCCESSFUL.");
									if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
									return 0;
								}

								if((inputms[0][y].freq = (double *)realloc(inputms[0][y].freq,(1)*sizeof(double))) == 0) {
									puts(" NOT SUCCESSFUL.");
									if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
									return 0;
								}
							}
						}
						*/
					}
					printf("\n   Expansion/reduction model (using a logistic equation) assigned.\n");
					break;
				case '2':
					/*subdivision*/
					printf("\n   All the events are considered here from the present to the past.");
					printf("\n   The population size No is counted in relation to the studied deme at present time (i.e., No = 1).\n");
					printf("\n   This model assume ALL studied samples are collected from the same single deme!!.\n");
					do {
						printf("\n   Number of total demes? ");
						scanf(" %d",&valuei);
						if(valuei <= 1)
							puts("Error: the number of demes must be more than one.\n");
					} while(valuei <= 1);
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].npop = valuei;
							/*realloc factor_pop*/
							if((inputms[0][y].factor_pop = (double *) realloc(inputms[0][y].factor_pop,(inputms[0][y].npop)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							memset(inputms[0][y].factor_pop,'\0',(inputms[0][y].npop)*sizeof(double));
						}
					}
					
					do {
						printf("\n   Migration parameter (4mNo, where m is the migration rate) among demes? ");
						scanf(" %lg",&value);
						if(value <= (double)0.)
							puts("Error: the value must be positive.\n");
					} while(value <= (double)0.);
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].mig_rate = value;
						}
					}
					
					printf("\n   Equal population size for each deme (y/n)? ");
					do *l = getchar();
					while(*l!='y' && *l!='n' && *l!='Y' && *l!='N');
					if(*l == 'y' || *l == 'Y') {
						for(y=0;y<data[0][0].n_loci;y++) {
							if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
								for(x=0;x<valuei;x++) inputms[0][y].factor_pop[x] = (double)1.;
							}
						}
					}
					else {
						printf("\n   Relative Population size N for deme[0] is fixed to 1.0\n");
						for(y=0;y<data[0][0].n_loci;y++) {
							if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
								inputms[0][y].factor_pop[0] = (double)1.0;
							}
						}
						for(x=1;x<valuei;x++) {
							do {
								printf("\n   Relative Population size N for deme[%d]? ",x);
								scanf(" %lg",&value);
								if(value <= (double)0.) 
									puts("Error: the value must be positive.\n");
							}while(value <= (double)0.);
							for(y=0;y<data[0][0].n_loci;y++) {
								if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
									inputms[0][y].factor_pop[x] = value;
								}
							}
						}
					}
					/*
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							if((inputms[0][y].nrec = (double *)realloc(inputms[0][y].nrec,(2)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].nrec[0] = (double)0;
							inputms[0][y].nrec[1] = (double)0;

							if((inputms[0][y].npast = (double *)realloc(inputms[0][y].npast,(2)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].npast[0] = (double)0;
							inputms[0][y].npast[1] = (double)0;

							if((inputms[0][y].tpast = (double *)realloc(inputms[0][y].tpast,(2)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].tpast[0] = (double)0;
							inputms[0][y].tpast[1] = (double)0;

							if((inputms[0][y].freq = (double *)realloc(inputms[0][y].freq,(1)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].freq[0] = (double)1;
						}
					}
					*/
					printf("\n   Subdivision (Wright island) model assigned.\n");
					break;
				case '3':
					/*split pop*/
					printf("\n   All the events are considered here from the present to the past.");
					printf("\n   The population size No is counted in relation to the studied deme at present time (i.e., No = 1).\n");
					printf("\n   This model assume a single population that was divided in the past and joined again.\n");
					do {
						printf("\n   Time in 4N generations (present -> past) when the population was rejoined? ");
						scanf(" %lg",&value);
						if(value < (double)0.)
							puts("Error: the value must be positive or zero.\n");
					} while(value < (double)0.);
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].time_split = value;
						}
					}
					
					do {
						printf("\n   Time in 4N generations (present -> past) when the population was splitted? ");
						scanf(" %lg",&value2);
						if(value2 <= value)
							puts("Error: the value must be larger than the time of rejoining.\n");
					} while(value2 <= value);
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].time_scoal = value2;
						}
					}
					
					do {
						printf("\n   Number of total refuges during the split? ");
						scanf(" %d",&valuei);
						if(valuei <= 0)
							puts("Error: the number of refuges can not be negative or zero.\n");
					} while(valuei <= 0);
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].npop = valuei;
						}
					}
					
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							/*realloc factor_pop, freq*/
							if((inputms[0][y].factor_pop = (double *) realloc(inputms[0][y].factor_pop,(inputms[0][y].npop)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							memset(inputms[0][y].factor_pop,'\0',(inputms[0][y].npop)*sizeof(double));

							if((inputms[0][y].freq = (double *)realloc(inputms[0][y].freq,(inputms[0][y].npop)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							memset(inputms[0][y].freq,'\0',(inputms[0][y].npop)*sizeof(double));
						}
					}
					
					do {
						printf("\n   Migration parameter (4mNo, where m is the migration rate) among refuges? ");
						scanf(" %lg",&value);
						if(value < (double)0.)
							puts("Error: the value must be positive or zero.\n");
					} while(value < (double)0.);
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].mig_rate = value;
						}
					}

					printf("\n   Equal population size for each refuge than in present (y/n)? ");
					do *l = getchar();
					while(*l!='y' && *l!='n' && *l!='Y' && *l!='N');
					if(*l == 'y' || *l == 'Y') {
						for(y=0;y<data[0][0].n_loci;y++) {
							if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
								for(x=0;x<inputms[0][y].npop;x++) {
									inputms[0][y].factor_pop[x] = (double)1.;
								}
							}
						}
					}
					else {
						for(x=0;x<valuei;x++) {
							do {
								printf("\n   Relative Population size N for refuge[%d]? ",x);
								scanf(" %lg",&value);
								if(value <= (double)0.) 
									puts("Error: the value must be positive.\n");
							}while(value <= (double)0.);
							for(y=0;y<data[0][0].n_loci;y++) {
								if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
									inputms[0][y].factor_pop[x] = value;
								}
							}
						}
					}
					printf("\n   Equal contribution to the present population for each refuge (y/n)? ");
					do *l = getchar();
					while(*l!='y' && *l!='n' && *l!='Y' && *l!='N');
					if(*l == 'y' || *l == 'Y') {
						for(y=0;y<data[0][0].n_loci;y++) {
							if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
								for(x=0;x<inputms[0][y].npop;x++) {
									inputms[0][y].freq[x] = (double)1./(double)inputms[0][y].npop;
								}
							}
						}
					}
					else {
						do {
							sum = (double)0.;
							for(x=0;x<valuei;x++) {
								do {
									printf("\n   Contribution (to the present population) of the refuge[%d]? ",x);
									scanf(" %lg",&value);
									if(value < (double)0.) 
										puts("Error: the value must be positive or zero.\n");
								}while(value < (double)0.);
								for(y=0;y<data[0][0].n_loci;y++) {
									if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
										inputms[0][y].freq[x] = value;
									}
								}
								sum += value;
							}
							if(sum<=(double)0 && sum > (double)1.+(double)1E-05)
								puts("Sorry: the sum of the contribution frequencies must sum less than 1 and more than zero.\n");
						}while(sum<=(double)0 && sum > (double)1.+(double)1E-05);
					} 
					
					do {
						printf("\n   Relative Population size N for the ancestral population? ");
						scanf(" %lg",&value);
						if(value <= (double)0.) 
							puts("Error: the value must be positive.\n");
					}while(value <= (double)0.);
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].factor_anc = value;
						}
					}
					/*
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							if((inputms[0][y].nrec = (double *)realloc(inputms[0][y].nrec,(2)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].nrec[0] = (double)0;
							inputms[0][y].nrec[1] = (double)0;

							if((inputms[0][y].npast = (double *)realloc(inputms[0][y].npast,(2)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].npast[0] = (double)0;
							inputms[0][y].npast[1] = (double)0;

							if((inputms[0][y].tpast = (double *)realloc(inputms[0][y].tpast,(2)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].tpast[0] = (double)0;
							inputms[0][y].tpast[1] = (double)0;
						}
					}
					*/
					printf("\n   Population split and merge model assigned.\n");
					break;
				case '4': 
					/*pselection*/
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].pop_sel = (double)0.;
							inputms[0][y].sel_nt = 0;
							inputms[0][y].sinit = (double)0.;
						}
					}

					printf("\n   All the events are considered here from the present to the past.");
					do {
						printf("\n   Population size N (i.e., in number of individuals)? ");
						scanf(" %lg",&g);
						if(g <= (double)0.)
							puts("Error: the value must be positive.\n");
					} while(g <= (double)0.);
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].pop_size = (long int)g;
						}
					}
					
					printf("\n   Nucleotide position where the selection acts (relative to the assigned loci/us)? ");
					printf("\n   (position 1 is the extreme left nucleotide of the loci. Negative values are allowed)\n   ");
					scanf(" %ld",&valueli);
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].sel_nt = valueli - 1;
						}
					}
					do {
						printf("\n   Selective value (4Ns, where s is the selective coefficient) in the position %ld (relative to the assigned loci/us)? ",valueli+1);
						scanf(" %lg",&value);
						if(value <= (double)0.) 
							puts("Error: the value must be positive.\n");
					}while(value <= (double)0.);
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].pop_sel = value;
						}
					}					
					printf("\n   Time (in 4N generations) since the selective process ended (relative to the assigned loci/us)?");
					printf("\n   (Negative values indicate that the selctive process has not yet finished)\n   ");
					scanf(" %lg",&value);
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							inputms[0][y].sinit = value;
						}
					}
					/*
					for(y=0;y<data[0][0].n_loci;y++) {
						if(inputms[0][y].definedmodel == data[0][0].n_events - n) {
							if((inputms[0][y].nrec = (double *)realloc(inputms[0][y].nrec,(2)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].nrec[0] = (double)0;
							inputms[0][y].nrec[1] = (double)0;

							if((inputms[0][y].npast = (double *)realloc(inputms[0][y].npast,(2)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].npast[0] = (double)0;
							inputms[0][y].npast[1] = (double)0;

							if((inputms[0][y].tpast = (double *)realloc(inputms[0][y].tpast,(2)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].tpast[0] = (double)0;
							inputms[0][y].tpast[1] = (double)0;

							if((inputms[0][y].factor_pop = (double *)realloc(inputms[0][y].factor_pop,(1)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].factor_pop[0] = (double)1;

							if((inputms[0][y].freq = (double *)realloc(inputms[0][y].freq,(1)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL.");
								if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
								return 0;
							}
							inputms[0][y].freq[0] = (double)1;
						}
						printf("\n   Strong positive selection model in one or several loci assigned.\n");
					}
					*/
					break;				
			}
		}
		*dataforsim = 1;
		return 1;
	}
    return 1;
}


