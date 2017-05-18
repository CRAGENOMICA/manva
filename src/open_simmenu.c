/*
 *  open_simmenu.c
 *  MuLoNeTests
 *
 *  Created by sebas on Mon Feb 24 2003.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>

void open_simmenu(struct globvar **global,struct statistics **matrix,struct statmulo **matrixml,struct var **data,
    struct var2b **inputms,struct statistisim ***matrixsim, struct statistisimmuloc **matrixmlsim,
    struct horizontalstatsml **avgstatloci,struct LRTdist **lrdist/*,struct MLthetasim ***mthetasim*/,FILE *file_output)
{
    char k[1];

	void print_datasimulations(int,struct var *,struct var2b *,FILE *,int,struct statistics *);
    void menu_coalescent_sim(struct globvar **, struct statistics *, struct statmulo *,struct statistisim ***,
        struct statistisimmuloc **,struct horizontalstatsml **,struct var **,struct var2b **,struct LRTdist ** /*,struct MLthetasim ****/,FILE *);
    void open_display_sim_menu(struct statistics *,struct statmulo *,struct statistisim **,struct statistisimmuloc *,
        struct horizontalstatsml *,struct var *,struct var2b *,int *,int *,int *,int *,int,FILE *,int);
    void save_simmatrix(struct globvar *,struct statistics *,struct statmulo *,struct var *, struct var2b *,
        struct statistisim **, struct statistisimmuloc *, struct horizontalstatsml *,struct LRTdist * /*,struct MLthetasim ***/,FILE *);
    int open_filematrix(struct globvar **,struct statistics **,struct statmulo **,struct var **, struct var2b **,
        struct statistisim ***, struct statistisimmuloc **, struct horizontalstatsml **,struct LRTdist ** /*,struct MLthetasim ****/,FILE *);
        
    while(1) {
        #if COMMAND_LINE
        printf("\n\nMENU 3. Coalescent Monte Carlo simulations and Statistical inference:\n\n");
        
        printf("CHOOSE:\n");
        printf(" 0 - DO Coalescent Monte Carlo simulations.\n");
        printf(" 1 - Display SIMULATED results (plus observed) menu.\n");
        printf(" 2 - SAVE a file with all the information from observed and/or simulated data.\n");
        printf(" 3 - LOAD a file with all the information from observed and/or simulated data.\n");
        printf(" 4 - Back to Main menu.\n\n");
		/*
        if(file_output) {
            fprintf(file_output,"\n\n     MENU:\n     3. Coalescent Monte Carlo simulations and Statistical inference:\n\n");
            
            fprintf(file_output," 0 - DO Coalescent Monte Carlo simulations.\n");
            fprintf(file_output," 1 - Display SIMULATED results (plus observed) menu.\n");
            fprintf(file_output," 2 - SAVE a file with all the information from observed and/or simulated data.\n");
            fprintf(file_output," 3 - LOAD a file with all the information from observed and/or simulated data.\n");
            fprintf(file_output," 4 - Back to Main menu.\n\n");
        }
		*/
        #endif
    
        do *k = getchar();
        while(*k<'0' || *k>'4');
        
        if(file_output) fprintf(file_output,"OPTION CHOSEN: %c\n\n",*k);

        switch(*k) {
            case '0':
                menu_coalescent_sim(global,*matrix,*matrixml,matrixsim,matrixmlsim,avgstatloci,data,inputms,lrdist/*,mthetasim*/,file_output);
                break;
            case '1':
				if(global[0][0].montecarlo_sim) {
					print_datasimulations(global[0][0].dataindata,*data,*inputms,file_output,global[0][0].dataequalsim,*matrix);
					open_display_sim_menu(*matrix,*matrixml,*matrixsim,*matrixmlsim,*avgstatloci,
						*data,*inputms,&global[0][0].observed_data,&global[0][0].outgroup,&global[0][0].n_loci,
						&global[0][0].montecarlo_sim,global[0][0].onlymulo,file_output,global[0][0].dataequalsim);
				}
				else {
					printf("\n\n     Sorry. Coalescent simulations should be first performed.\n\n");
				}
                break;
            case '2':
                save_simmatrix(*global,*matrix,*matrixml,*data,*inputms,*matrixsim,*matrixmlsim,*avgstatloci,*lrdist/*,*mthetasim*/,file_output);
                break;
            case '3':
                if(open_filematrix(global,matrix,matrixml,data,inputms,matrixsim,matrixmlsim,avgstatloci,lrdist/*,mthetasim*/,file_output) == 0) {
					global[0][0].n_loci = 0;
					global[0][0].observed_data = 0;
					global[0][0].outgroup = -1;
					/*global[0][0].maxloc = 0;*/
					
					global[0][0].montecarlo_sim = 0;
					global[0][0].nlocimat = 0;
					global[0][0].nitermat = 0;
					global[0][0].dataforsim = 0;
					global[0][0].onlymulo = 0;
					
					global[0][0].dataindata = 0;
					global[0][0].dataequalsim = 0;

					global[0][0].ml0done = 0;
					global[0][0].nlocimatt = 0;
					global[0][0].nitermatt = 0;
					global[0][0].dataforsimt = 0;
					global[0][0].mltheta = 0;
                }
                break;
            case '4':
               	return;
                break;
        }
    }
    return;
}

void open_display_sim_menu(struct statistics *matrix,struct statmulo *matrixml,struct statistisim **matrixsim,struct statistisimmuloc *matrixmlsim, 
	struct horizontalstatsml *avgstatloci, struct var *data,struct var2b *inputms,int *observed_data,int *outgroup,int *n_loci,int *montecarlo_sim,
	int onlymulo, FILE *file_output, int dataobsequalsim)
{
    char k[1];
    void open_simdisplmenu(struct statistics *,struct statmulo *,struct statistisim **,struct statistisimmuloc *, struct horizontalstatsml *, 
		struct var *,struct var2b *,int *,int *,int * /*,int * */,int , FILE *,int,int );
    
    if(*montecarlo_sim) {    
        while(1) {
            #if COMMAND_LINE
            printf("\n\nMENU 3. Coalescent Monte Carlo simulations and Statistical inference:");
            printf("\n     3.1. Display coalescent simulation analysis menu:\n\n");
            
            printf("CHOOSE:\n");
            printf(" 0 - Display summary statistics.\n");
            printf(" 1 - Display neutrality tests.\n");
            printf(" 2 - Back to Main menu.\n");
			/*
            if(file_output) {
                fprintf(file_output,"\n\n     MENU:\n     3. Coalescent Monte Carlo simulations and Statistical inference:");
                fprintf(file_output,"\n     3.1. Display coalescent simulation analysis menu:\n\n");
                
                fprintf(file_output," 0 - Display summary statistics.\n");
                fprintf(file_output," 1 - Display neutrality tests.\n");
                fprintf(file_output," 2 - Back to Main menu.\n");
            }
			*/
            #endif
            
            do *k = getchar();
            while(*k<'0' || *k>'2');
            
            if(file_output) fprintf(file_output,"OPTION CHOSEN: %c\n\n",*k);

            switch(*k) {
                case '0':
                    open_simdisplmenu(matrix,matrixml,matrixsim,matrixmlsim,avgstatloci,data,inputms,
                        observed_data,outgroup,n_loci/*, montecarlo_sim*/,onlymulo,file_output,dataobsequalsim,0);
                    break;
                case '1':
                    open_simdisplmenu(matrix,matrixml,matrixsim,matrixmlsim, avgstatloci, data,inputms,
                        observed_data,outgroup,n_loci/*, montecarlo_sim*/,onlymulo,file_output,dataobsequalsim,1);
                    break;
                case '2':
                    return;
                    break;
            }
        }
    }
    else {
        printf("\n\n     coalescent simulations should be first performed.\n\n");
        printf("      Press 0.\n\n");
        do {
            *k = getchar();
        }while(*k!='0');
        return;
    }
    return;
}

void open_simdisplmenu(struct statistics *matrix,struct statmulo *matrixml,struct statistisim **matrixsim,struct statistisimmuloc *matrixmlsim,
	struct horizontalstatsml *avgstatloci,struct var *data,struct var2b *inputms,int *observed_data,int *outgroup,int *n_loci/*,int *montecarlo_sim*/,
	int onlymulo, FILE *file_output,int dataobsequalsim,int neuttest)
{
    char k[1];
    void open_simdisplmulomenu(struct statistics *,struct statmulo *,struct statistisim **,struct statistisimmuloc *, struct horizontalstatsml *,
		struct var *,struct var2b *,int *,int *,int * /*,int * */,int , FILE *,int,int,int);
    
    while(1) {
        #if COMMAND_LINE
        printf("\n\nMENU 3. Coalescent Monte Carlo simulations and Statistical inference:");
        printf("\n     3.1. Display coalescent simulation analysis menu:");
        if(!neuttest)
            printf("\n     3.1.0. Display statistics menu:\n\n");
        else
            printf("\n     3.1.1. Display neutrality tests menu:\n\n");
        
        printf("CHOOSE:\n");
        printf(" 0 - Display MULTILOCUS results.\n");
        printf(" 1 - Display results for each locus.\n");
        printf(" 2 - Back to previous menu.\n");
		/*
        if(file_output) {
            fprintf(file_output,"\n\n     MENU:\n     3. Coalescent Monte Carlo simulations and Statistical inference:");
            fprintf(file_output,"\n     3.1. Display coalescent simulation analysis menu:");
            if(!neuttest)
                fprintf(file_output,"\n     3.1.0. Display statistics menu:\n\n");
            else
                fprintf(file_output,"\n     3.1.1. Display neutrality tests menu:\n\n");
            
            fprintf(file_output," 0 - Display MULTILOCUS results.\n");
            fprintf(file_output," 1 - Display results for each locus.\n");
            fprintf(file_output," 2 - Back to previous menu.\n");
        }
		*/
        #endif
        
        do *k = getchar();
        while(*k<'0' || *k>'2');
        
        if(file_output) fprintf(file_output,"OPTION CHOSEN: %c\n\n",*k);

        switch(*k) {
            case '0':
                open_simdisplmulomenu(matrix,matrixml,matrixsim,matrixmlsim,avgstatloci,data,inputms,
                    observed_data, outgroup, n_loci,/* montecarlo_sim,*/onlymulo,file_output,dataobsequalsim,neuttest,1);
                break;
            case '1':
                if(onlymulo == 0)
                open_simdisplmulomenu(matrix,matrixml,matrixsim,matrixmlsim,avgstatloci,data,inputms,
                    observed_data,outgroup,n_loci,/* montecarlo_sim,*/onlymulo,file_output,dataobsequalsim,neuttest,0);
                else {
                    printf("\n Only multilocus analyses are available now.\n");
                    if(file_output) fputs("\n Only multilocus analyses are available now.\n",file_output);
                    return;
                }
                break;
            case '2':
                return;
                break;
        }
    }
    return;
}

void open_simdisplmulomenu(struct statistics *matrix,struct statmulo *matrixml,struct statistisim **matrixsim,struct statistisimmuloc *matrixmlsim,
	struct horizontalstatsml *avgstatloci,struct var *data,struct var2b *inputms,int *observed_data,int *outgroup,int *n_loci/*,int *montecarlo_sim*/,
	int onlymulo, FILE *file_output,int dataobsequalsim, int neuttest,int mulo)
{
    char k[1],l[1];
    void open_simdispldetgrid(struct var *,struct var2b *,struct statistisim **,struct statistisimmuloc *,struct horizontalstatsml *,int *,/*int *,*/
		int , FILE *,int,int,int,struct statistics *,int);
    void open_simdisplayprob(struct var *,struct var2b *,struct statistics *,struct statmulo *,struct statistisim **,struct statistisimmuloc *,
		struct horizontalstatsml *,int *,int * /*,int *,int * */,int , FILE *,int,int,int);
    void open_simdisplayhist(struct var *,struct var2b *,struct statistics *,struct statmulo *,struct statistisim **,struct statistisimmuloc *,
		/*struct horizontalstatsml *,int *, */int *,int * /*,int *,int */, FILE *,int,int,int);
    
    while(1) {
        #if COMMAND_LINE
        printf("\n\nMENU 3. Coalescent Monte Carlo simulations and Statistical inference:");
        printf("\n     3.1. Display coalescent simulation analysis menu:");
        if(!neuttest)
            printf("\n     3.1.0. Display statistics menu:");
        else
            printf("\n     3.1.1. Display neutrality tests menu:");
        if(mulo) {
            if(neuttest)
                printf("\n     3.1.1.0. Display multilocus neutrality tests menu:\n\n");
            else 
                printf("\n     3.1.0.0. Display multilocus statistics menu:\n\n");
        }
        else {
            if(neuttest)
                printf("\n     3.1.1.1. Display detailed neutrality tests menu:\n\n");
            else 
                printf("\n     3.1.0.1. Display detailed statistics menu:\n\n");
        }
        
        printf("CHOOSE:\n");
        printf(" 0 - Display the results in a table (WARNING: VERY LONG OUTPUT!).\n");
        printf(" 1 - Display the probability values (from observed data) or confident intervals.\n");
        printf(" 2 - Display histograms of simulated results.\n");
        printf(" 3 - Back to to previous menu.\n");
		/*
        if(file_output) {    
            fprintf(file_output,"\n\n     MENU:\n     3. Coalescent Monte Carlo simulations and Statistical inference:");
            fprintf(file_output,"\n     3.1. Display coalescent simulation analysis menu:");
            if(!neuttest)
                fprintf(file_output,"\n     3.1.0. Display statistics menu:");
            else
                fprintf(file_output,"\n     3.1.1. Display neutrality tests menu:");
            if(mulo) {
                if(neuttest)
                    fprintf(file_output,"\n     3.1.1.0. Display multilocus neutrality tests menu:\n\n");
                else 
                    fprintf(file_output,"\n     3.1.0.0. Display multilocus statistics menu:\n\n");
            }
            else {
                if(neuttest)
                    fprintf(file_output,"\n     3.1.1.1. Display detailed neutrality tests menu:\n\n");
                else 
                    fprintf(file_output,"\n     3.1.0.1. Display detailed statistics menu:\n\n");
            }
            
            fprintf(file_output,"\n 0 - Display the results in a table (WARNING: VERY LONG OUTPUT!).\n");
            fprintf(file_output," 1 - Display the probability values (from observed data) or confident intervals.\n");
            fprintf(file_output," 2 - Display histograms of simulated results.\n");
            fprintf(file_output," 3 - Back to to previous menu.\n");
        }
		*/
        #endif
        
        do *k = getchar();
        while(*k<'0' || *k>'3');
        
        if(file_output) fprintf(file_output,"\n\nOPTION CHOSEN: %c\n\n",*k);

        switch(*k) {
            case '0':
				/*
				printf("\nAre you sure (y/n)?\n");
				do *l = getchar();
				while(*l != 'y' && *l != 'n');
				*/
				*l='y';
				if(*l=='y') {
					printf(" This output is not displayed in the console.\n");
					printf(" The output will be send to a new output file (tab delimited).\n\n");
					open_simdispldetgrid(data,inputms,matrixsim,matrixmlsim,avgstatloci,outgroup/*,n_loci*/,
						onlymulo, file_output,neuttest,mulo,dataobsequalsim,matrix,*observed_data);
					if(file_output) printf(" ... Done.\n");
				}
                break;
            case '1':
                open_simdisplayprob(data,inputms,matrix,matrixml,matrixsim,matrixmlsim,avgstatloci,
                    observed_data, outgroup/*,n_loci,montecarlo_sim*/,onlymulo,file_output,dataobsequalsim,
                    neuttest,mulo);
                break;
            case '2':
                open_simdisplayhist(data,inputms,matrix,matrixml,matrixsim,matrixmlsim/*,avgstatloci,
                    observed_data*/,outgroup,n_loci/*, montecarlo_sim,onlymulo*/,file_output,dataobsequalsim,
                    neuttest,mulo);
                break;
            case '3':
                return;
                break;
        }
    }
    return;
}

