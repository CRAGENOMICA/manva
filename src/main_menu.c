/*
 *  main_menu.c
 *  MuLoNeTests
 *
 *  Created by sebas on Sat Feb 22 2003.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>

int open_mainmenu(struct globvar **global,struct statistics **matrix,struct statmulo **matrixml,struct var **datams,
    struct var2b **inputms,struct statistisim ***matrixsim,struct statistisimmuloc **matrixmlsim,
    struct horizontalstatsml **avgstatloci,struct LRTdist **lrdist/*,struct MLthetasim ***mthetasim*/,FILE **file_output)
{
    char k[1],l[1];
    void open_filemenu(struct globvar **,struct statistics **,struct statmulo **,struct var **, struct var2b **,
        struct statistisim ***, struct statistisimmuloc **, struct horizontalstatsml **,struct LRTdist ** /*,struct MLthetasim ****/,FILE **);
    void open_obsvmenu(struct statistics *,struct statmulo *,int *,int *,int *,char *,FILE *);
    void open_simmenu(struct globvar **,struct statistics **,struct statmulo **,struct var **, struct var2b **,
        struct statistisim ***, struct statistisimmuloc **, struct horizontalstatsml **,struct LRTdist ** /*,struct MLthetasim ****/,FILE *);
    void save_simmatrix(struct globvar *,struct statistics *,struct statmulo *,struct var *, struct var2b *,
        struct statistisim **, struct statistisimmuloc *, struct horizontalstatsml *,struct LRTdist * /*,struct MLthetasim ***/,FILE *);
    void open_helpmenu(FILE *);
	void open_paraml_menu(struct globvar **,struct statistics **,struct statmulo **,struct var **,
		struct var2b **,struct statistisim ***, struct statistisimmuloc **,
		struct horizontalstatsml **,struct LRTdist ** /*,struct MLthetasim ****/,FILE *);
	int open_exit();
    
    while(1) {
        #if COMMAND_LINE
        printf("\n\nMENU Main menu:\n\n");
        
        printf("CHOOSE:\n");
        printf(" 0 - Display FILE menu.\n");
        printf(" 1 - Display OBSERVED data menu.\n");
        printf(" 2 - ESTIMATION of variability by Maximum likelihood.\n");
        printf(" 3 - COALESCENT simulations and Statistical inference.\n");
        /*printf(" 4 - Display DOCUMENTATION and additional information.\n");*/
        printf(" 4 - EXIT of the program.\n\n");
		fflush(stdout);
		/*
		if(*file_output) { 
            fprintf(*file_output,"\n\n     Main menu:\n\n");
            
            fprintf(*file_output," 0 - Display FILE menu.\n");
            fprintf(*file_output," 1 - Display OBSERVED data menu.\n");
            fprintf(*file_output," 2 - ESTIMATION of variability by Maximum likelihood\n");
            fprintf(*file_output," 3 - COALESCENT simulations and Statistical inference.\n");
            fprintf(*file_output," 4 - Display DOCUMENTATION and additional information.\n");
            fprintf(*file_output," 5 - EXIT of the program.\n\n");
        }
		*/
        #endif
    
        do *k = getchar();
        while(*k<'0' || *k>'4');
        
        if(*file_output) fprintf(*file_output,"OPTION CHOSEN: %c\n\n",*k);

        switch(*k) {
            case '0':
                open_filemenu(global,matrix,matrixml,datams,inputms,matrixsim,matrixmlsim,avgstatloci,lrdist/*,mthetasim*/,file_output);
                break;
            case '1':
                open_obsvmenu(*matrix,*matrixml,&global[0][0].observed_data,&global[0][0].outgroup,
                    &global[0][0].n_loci,global[0][0].subset_positions,*file_output);
                break;
            case '2':
				#if MLtheta_open
				open_paraml_menu(global,matrix,matrixml,datams,inputms,matrixsim,matrixmlsim,avgstatloci,lrdist/*,mthetasim*/,*file_output);
				#else
				printf("\n Sorry. This section is not yet available.\n");
				printf("Press 0 to come back to the main Menu.\n");
				do *l = getchar();
				while(*l!='0' && *l!='9');
				if(*l == '9') open_paraml_menu(global,matrix,matrixml,datams,inputms,matrixsim,matrixmlsim,avgstatloci,lrdist/*,mthetasim*/,*file_output);
				#endif
                break;
            case '3':
                open_simmenu(global,matrix,matrixml,datams,inputms,matrixsim,matrixmlsim,avgstatloci,lrdist/*,mthetasim*/,*file_output);
                break;
            /*
			case '4':
                open_helpmenu(*file_output);
                break;
            */
			case '4':
                if(global[0][0].observed_data || global[0][0].dataforsim || global[0][0].montecarlo_sim) {
                    printf(" Do you want to keep all observed and simulated information in a MANVa file (y/n)? "); 
                    do *l = getchar();
                    while(*l !='y' && *l !='n' && *l !='Y' && *l !='N');
                    if(*l == 'y' || *l == 'Y') {
                        save_simmatrix(*global,*matrix,*matrixml,*datams,*inputms,*matrixsim,*matrixmlsim,*avgstatloci,*lrdist/*,*mthetasim*/,*file_output);
                    }
                }
                if((open_exit()) == 0)
					return 0;
				break;
        }
    }
}

