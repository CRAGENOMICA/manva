/*
 *  opensavematrix.c
 *  MuLoNeTests
 *
 *  Created by sebas on Sun Mar 09 2003.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void save_simmatrix(struct globvar *global,struct statistics *matrix,struct statmulo *matrixml,struct var *datams,
    struct var2b *inputms,struct statistisim **matrixsim, struct statistisimmuloc *matrixmlsim,
    struct horizontalstatsml *avgstatloci,struct LRTdist *lrdist/*,struct MLthetasim **mthetasim*/,FILE *file_output)
{
    FILE *extfile;
    char name_extfile[512]; /*511 characters...*/
    int x;
    
    printf("\n SAVE a MANVA file with extracted information:\n");

    if(global[0].observed_data || global[0].montecarlo_sim  || global[0].dataforsim || global[0].dataforsimt) {
        *name_extfile = 0;
        do {
            puts(" Please indicate the name of the file to save the extracted information: ");
            scanf(" %[^\n]",name_extfile);
        } while(*name_extfile == 0);
        name_extfile[511] = '\0';
        if(*name_extfile != 0) {
            if((extfile = fopen (name_extfile,"wb")) == 0) {
                puts("\n\nError in input/output");
                return;
            }
        }
		else {
			return;
		}
    }
    else {
        puts("\nSorry, There is no data in the memory.");
        if(file_output) fputs("\nSorry, There is no data in the memory.",file_output);
        return;
    }

    /*first the struct global*/
    /*global[0].maxloc = 0;*/
    if((fwrite(global,sizeof(struct globvar),1,extfile)) != 1) {
        puts("Error: it is not possible to save the data (0)");
        if(file_output) fputs("Error: it is not possible to save the data (0)",file_output);
		fclose(extfile);
        return;
    }
    /*then the observed data*/
    if(global[0].observed_data) {
        if((fwrite(matrix,sizeof(struct statistics),(unsigned)global[0].n_loci,extfile)) != (unsigned)global[0].n_loci) {
            puts("Error: it is not possible to save the data (1)");
            if(file_output) fputs("Error: it is not possible to save the data (1)",file_output);
			fclose(extfile);
            return;
        }
        if((fwrite(matrixml,sizeof(struct statmulo),1,extfile)) != 1) {
            puts("Error: it is not possible to save the data (2)");
            if(file_output) fputs("Error: it is not possible to save the data (2)",file_output);
			fclose(extfile);
            return;
        }
    }
    /*then the simulated data*/
    if(global[0].dataforsim == 1 || global[0].dataforsimt == 1) {
        if((fwrite(datams,sizeof(struct var),1,extfile)) != 1) {
            puts("Error: it is not possible to save the data (3)");
            if(file_output) fputs("Error: it is not possible to save the data (3)",file_output);
			fclose(extfile);
            return;
        }
        if((fwrite(inputms,sizeof(struct var2b),(unsigned)datams[0].n_loci,extfile)) != (unsigned)datams[0].n_loci) {
            puts("Error: it is not possible to save the data (4)");
            if(file_output) fputs("Error: it is not possible to save the data (4)",file_output);
			fclose(extfile);
            return;
        }
		for(x=0;x<datams[0].n_loci;x++) {
			if((fwrite(inputms[x].nrec,sizeof(double),(unsigned)(2+inputms[x].nintn),extfile)) != (unsigned)(2+inputms[x].nintn)) {
				puts("Error: it is not possible to save the data (5)");
				if(file_output) fputs("Error: it is not possible to save the data (5)",file_output);
				fclose(extfile);
				return;
			}
			if((fwrite(inputms[x].npast,sizeof(double),(unsigned)(2+inputms[x].nintn),extfile)) != (unsigned)(2+inputms[x].nintn)) {
				puts("Error: it is not possible to save the data (6)");
				if(file_output) fputs("Error: it is not possible to save the data (6)",file_output);
				fclose(extfile);
				return;
			}
			if((fwrite(inputms[x].tpast,sizeof(double),(unsigned)(2+inputms[x].nintn),extfile)) != (unsigned)(2+inputms[x].nintn)) {
				puts("Error: it is not possible to save the data (7)");
				if(file_output) fputs("Error: it is not possible to save the data (7)",file_output);
				fclose(extfile);
				return;
			}
			if(inputms[x].npop) {
				if((fwrite(inputms[x].factor_pop,sizeof(double),(unsigned)inputms[x].npop,extfile)) != (unsigned)inputms[x].npop) {
					puts("Error: it is not possible to save the data (8)");
					if(file_output) fputs("Error: it is not possible to save the data (8)",file_output);
					fclose(extfile);
					return;
				}
				if((fwrite(inputms[x].freq,sizeof(double),(unsigned)inputms[x].npop,extfile)) != (unsigned)inputms[x].npop) {
					puts("Error: it is not possible to save the data (9)");
					if(file_output) fputs("Error: it is not possible to save the data (9)",file_output);
					fclose(extfile);
					return;
				}
			}
		}
        if(global[0].montecarlo_sim == 1) {
            if((fwrite(matrixmlsim,sizeof(struct statistisimmuloc),(unsigned long)global[0].nitermat,extfile)) != (unsigned long)global[0].nitermat) {
                puts("Error: it is not possible to save the data (10)");
                if(file_output) fputs("Error: it is not possible to save the data (10)",file_output);
				fclose(extfile);
                return;
            }
            if((fwrite(avgstatloci,sizeof(struct horizontalstatsml),(unsigned)global[0].nlocimat,extfile)) != (unsigned)global[0].nlocimat) {
                puts("Error: it is not possible to save the data (11)");
                if(file_output) fputs("Error: it is not possible to save the data (11)",file_output);
				fclose(extfile);
                return;
            }
            if(global[0].onlymulo == 0) {
                for(x=0;x<global[0].nlocimat;x++) {
                    if((fwrite(matrixsim[x],sizeof(struct statistisim),(unsigned long)global[0].nitermat,extfile)) != (unsigned long)global[0].nitermat) {
                        puts("Error: it is not possible to save the data (12)");
                        if(file_output) fputs("Error: it is not possible to save the data (12)",file_output);
						fclose(extfile);
                        return;
                    }
                }
            }
        }
        if(global[0].ml0done == 1) {
            if((fwrite(lrdist,sizeof(struct LRTdist),(unsigned long)global[0].nitermatt,extfile)) != (unsigned long)global[0].nitermatt) {
                puts("Error: it is not possible to save the data (13)");
                if(file_output) fputs("Error: it is not possible to save the data (13)",file_output);
				fclose(extfile);
                return;
            }
			/*
			for(x=0;x<global[0].nlocimatt;x++) {
				if((fwrite(mthetasim[x],sizeof(struct MLthetasim),global[0].nitermatt,extfile)) != global[0].nitermatt) {
					puts("Error: it is not possible to save the data (14)");
					if(file_output) fputs("Error: it is not possible to save the data (14)",file_output);
					fclose(extfile);
					return;
				}
			}
			*/
		}
	}
	fflush(extfile);
	fclose(extfile);
    
    printf("\n\nMANVa format file:  %s. Saved. \n\n",name_extfile);
    if(file_output) fprintf(file_output,"\nMANVa format file:  %s. Saved. \n\n",name_extfile);
    
    return;
}

int open_filematrix(struct globvar **global,struct statistics **matrix,struct statmulo **matrixml,struct var **datams,
    struct var2b **inputms,struct statistisim ***matrixsim, struct statistisimmuloc **matrixmlsim,
    struct horizontalstatsml **avgstatloci,struct LRTdist **lrdist/*,struct MLthetasim ***mthetasim*/, 
	FILE *file_output)
{
    void print_datasimulations(int,struct var *,struct var2b *, FILE *,int,struct statistics *);
	void print_seqinfodata_ml0(struct globvar *,struct statistics * /*,struct statmulo * */,struct var *, struct var2b *,FILE *);
    void print_seqinfodata_sim(int,struct var *,struct var2b *, FILE *,int,struct statistics *);
    void print_evocond_sim(struct var *,struct var2b *, FILE * /*,int */,struct statistics *);
    void erase(struct globvar **,struct statistics ** /*,struct statmulo * */,struct statistisim ***,
        struct statistisimmuloc **, struct horizontalstatsml **,int *,struct var2b **,struct LRTdist **,struct var **,FILE *,int *);    
    FILE *extfile;
    char name_extfile[512]; /*511 characters.*/
    char *f;
    char k[1];
    int x;
    int ex=0,in=0;
	double *nrec,*npast,*tpast,*factor_pop,*freq;
	int n,m;
    
	char tripletsU[64][3] =
	{
		{"UUU"},/*0*/
		{"UUC"},/*1*/
		{"UUA"},/*2*/
		{"UUG"},/*3*/
		{"UCU"},/*4*/
		{"UCC"},/*5*/
		{"UCA"},/*6*/
		{"UCG"},/*7*/
		{"UAU"},/*8*/
		{"UAC"},/*9*/
		{"UAA"},/*10*/
		{"UAG"},/*11*/
		{"UGU"},/*12*/
		{"UGC"},/*13*/
		{"UGA"},/*14*/
		{"UGG"},/*15*/
		{"CUU"},/*16*/
		{"CUC"},/*17*/
		{"CUA"},/*18*/
		{"CUG"},/*19*/
		{"UCU"},/*20*/
		{"CCC"},/*21*/
		{"CCA"},/*22*/
		{"CCG"},/*23*/
		{"CAU"},/*24*/
		{"CAC"},/*25*/
		{"CAA"},/*26*/
		{"CAG"},/*27*/
		{"CGU"},/*28*/
		{"CGC"},/*29*/
		{"CGA"},/*30*/
		{"CGG"},/*31*/
		{"AUU"},/*32*/
		{"AUC"},/*33*/
		{"AUA"},/*34*/
		{"AUG"},/*35*/
		{"ACU"},/*36*/
		{"ACC"},/*37*/
		{"ACA"},/*38*/
		{"ACG"},/*39*/
		{"AAU"},/*40*/
		{"AAC"},/*41*/
		{"AAA"},/*42*/
		{"AAG"},/*43*/
		{"AGU"},/*44*/
		{"AGC"},/*45*/
		{"AGA"},/*46*/
		{"AGG"},/*47*/
		{"GUU"},/*48*/
		{"GUC"},/*49*/
		{"GUA"},/*50*/
		{"GUG"},/*51*/
		{"GCU"},/*52*/
		{"GCC"},/*53*/
		{"GCA"},/*54*/
		{"GCG"},/*55*/
		{"GAU"},/*56*/
		{"GAC"},/*57*/
		{"GAA"},/*58*/
		{"GAG"},/*59*/
		{"GGU"},/*60*/
		{"GGC"},/*61*/
		{"GGA"},/*62*/
		{"GGG"},/*63*/
	};
/*
	char tripletsN[64][3] =
	{
		{"111"},
		{"112"},
		{"114"},
		{"113"},
		{"121"},
		{"122"},
		{"124"},
		{"123"},
		{"141"},
		{"142"},
		{"144"},
		{"143"},
		{"131"},
		{"132"},
		{"134"},
		{"133"},
		{"211"},
		{"212"},
		{"214"},
		{"213"},
		{"221"},
		{"222"},
		{"224"},
		{"223"},
		{"241"},
		{"242"},
		{"244"},
		{"243"},
		{"231"},
		{"232"},
		{"234"},
		{"233"},
		{"411"},
		{"412"},
		{"414"},
		{"413"},
		{"421"},
		{"422"},
		{"424"},
		{"423"},
		{"441"},
		{"442"},
		{"444"},
		{"443"},
		{"431"},
		{"432"},
		{"434"},
		{"433"},
		{"311"},
		{"312"},
		{"314"},
		{"313"},
		{"321"},
		{"322"},
		{"324"},
		{"323"},
		{"341"},
		{"342"},
		{"344"},
		{"343"},
		{"331"},
		{"332"},
		{"334"},
		{"333"},
	};
*/	
    printf("\n\nLOAD a MANVA file with extracted information:\n\n");
    
    if(global[0][0].observed_data != 0 || global[0][0].montecarlo_sim != 0 || global[0][0].dataforsim != 0) {
        printf("\n\nThis command is going to erase all the information in the memory. Do you want to do that (y/n)? ");
        do *k = getchar();
        while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');
    }
    else *k = 'y';

    if(*k == 'n' || *k == 'N') return 1;
    else {
        erase(global,matrix/*,*matrixml*/,matrixsim,matrixmlsim,avgstatloci,&ex,inputms,lrdist,datams,file_output,&in);
    }

    if(file_output == 0) {
        puts("\n\n  You do not have an output file defined. It is recommended to send all results to an output file. ");
        puts("\n  Do you want to return to the File Menu and define an output file (y/n)? ");
        
        do *k = getchar();
        while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');        
        if(*k == 'y' || *k == 'Y') return 1;
    }

    do {
        puts("\n  Name of the file?");
        *name_extfile = 0;
        scanf(" %[^\n]",name_extfile);
    } while(*name_extfile == 0);
    name_extfile[511] = '\0';
    
    if((extfile = fopen (name_extfile,"rb")) == 0) {
        printf("\n  It is not possible to open the file %s.",name_extfile);
    }
    if(extfile) {
        if(!(f = (char *)malloc(BUFSIZ))) {
            puts("\nError: memory not reallocated. openfilematrix.c.1 \n");
            return 0;
        }
        setbuf(extfile,f);
        
        printf("\n Opening extract file %s:",name_extfile);
        if(file_output) fprintf(file_output,"\n Opening extract file %s:",name_extfile);

        /*first read the struct global*/
        if((fread(global[0],sizeof(struct globvar),1,extfile)) != 1) {
            puts(" NOT SUCCESSFUL (0).");
            if(file_output) fputs(" NOT SUCCESSFUL (0).",file_output);
			fclose(extfile);
            return 0;
        }
        /*then read the observed data*/
        if(global[0][0].observed_data && global[0][0].n_loci) {
            if((*matrix = (struct statistics *) realloc(*matrix,(global[0][0].n_loci)*sizeof(struct statistics))) == 0) {
                puts(" NOT SUCCESSFUL (1).");
                if(file_output) fputs(" NOT SUCCESSFUL (0).",file_output);
				fclose(extfile);
                return 0;
            }
			memset(*matrix,'\0',(global[0][0].n_loci)*sizeof(struct statistics));
            if((fread(*matrix,sizeof(struct statistics),(unsigned)global[0][0].n_loci,extfile)) != (unsigned)global[0][0].n_loci) {
                puts(" NOT SUCCESSFUL (2).");
                if(file_output) fputs(" NOT SUCCESSFUL (2).",file_output);
				fclose(extfile);
                return 0;
            }
            if((fread(*matrixml,sizeof(struct statmulo),1,extfile)) != 1) {
                puts(" NOT SUCCESSFUL (3).");
                if(file_output) fputs(" NOT SUCCESSFUL (3).",file_output);
				fclose(extfile);
                return 0;
            }
        }
        /*then read the simulated data*/
        if(global[0][0].dataindata > 0) {
			if(global[0][0].dataforsim == 1 || global[0][0].dataforsimt == 1) {
				if((fread(*datams,sizeof(struct var),1,extfile)) != 1) {
					puts(" NOT SUCCESSFUL (4).");
					if(file_output) fputs(" NOT SUCCESSFUL (4).",file_output);
					fclose(extfile);
					return 0;
				}
				/*memory locations*/
				nrec  = inputms[0][0].nrec;
				npast = inputms[0][0].npast;
				tpast = inputms[0][0].tpast;
				factor_pop = inputms[0][0].factor_pop;
				freq = inputms[0][0].freq;
				
				if(global[0][0].nlocimat) {
					if((*inputms = (struct var2b *) realloc(*inputms,(global[0][0].nlocimat)*sizeof(struct var2b))) == 0) {
						puts(" NOT SUCCESSFUL (5).");
						if(file_output) fputs(" NOT SUCCESSFUL (5).",file_output);
						fclose(extfile);
						return 0;
					}
					memset(*inputms,'\0',(global[0][0].nlocimat)*sizeof(struct var2b));
					if((fread(*inputms,sizeof(struct var2b),(unsigned long)global[0][0].nlocimat,extfile)) != (unsigned long)global[0][0].nlocimat) {
						puts(" NOT SUCCESSFUL (6).");
						if(file_output) fputs(" NOT SUCCESSFUL (6).",file_output);
						fclose(extfile);
						return 0;
					}
				}
				/*recover memory locations*/
				inputms[0][0].nrec = nrec;
				inputms[0][0].npast = npast;
				inputms[0][0].tpast = tpast;
				inputms[0][0].factor_pop = factor_pop;
				inputms[0][0].freq = freq;

				for(x=0;x<global[0][0].nlocimat;x++) {
					/*calloc/realloc and read*/
					if(x>0) {
						if((inputms[0][x].nrec = (double *)calloc((2+inputms[0][x].nintn),sizeof(double))) == 0) {
							puts(" NOT SUCCESSFUL (7).");
							if(file_output) fputs(" NOT SUCCESSFUL (7).",file_output);
							fclose(extfile);
							return 0;
						}
					}
					else {
						if((inputms[0][x].nrec = (double *)realloc(inputms[0][x].nrec,(2+inputms[0][x].nintn)*sizeof(double))) == 0) {
							puts(" NOT SUCCESSFUL (8).");
							if(file_output) fputs(" NOT SUCCESSFUL (8).",file_output);
							fclose(extfile);
							return 0;
						}
					}
					memset(inputms[0][x].nrec,'\0',(2+inputms[0][x].nintn)*sizeof(double));					
					if((fread(inputms[0][x].nrec,sizeof(double),(unsigned)(2+inputms[0][x].nintn),extfile)) != (unsigned)(2+inputms[0][0].nintn)) {
						puts("Error: it is not possible to read the data");
						if(file_output) fputs("Error: it is not possible to read the data",file_output);
						fclose(extfile);
						return 0;
					}
					
					if(x>0) {
						if((inputms[0][x].npast = (double *)calloc((2+inputms[0][x].nintn),sizeof(double))) == 0) {
							puts(" NOT SUCCESSFUL (9).");
							if(file_output) fputs(" NOT SUCCESSFUL (9).",file_output);
							fclose(extfile);
							return 0;
						}
					}
					else {
						if((inputms[0][x].npast = (double *)realloc(inputms[0][x].npast,(2+inputms[0][x].nintn)*sizeof(double))) == 0) {
							puts(" NOT SUCCESSFUL (10).");
							if(file_output) fputs(" NOT SUCCESSFUL (10).",file_output);
							fclose(extfile);
							return 0;
						}
					}
					memset(inputms[0][x].npast,'\0',(2+inputms[0][x].nintn)*sizeof(double));
					if((fread(inputms[0][x].npast,sizeof(double),(unsigned)(2+inputms[0][x].nintn),extfile)) != (unsigned)(2+inputms[0][x].nintn)) {
						puts("Error: it is not possible to read the data");
						if(file_output) fputs("Error: it is not possible to read the data",file_output);
						fclose(extfile);
						return 0;
					}
					
					if(x>0) {
						if((inputms[0][x].tpast = (double *)calloc((2+inputms[0][x].nintn),sizeof(double))) == 0) {
							puts(" NOT SUCCESSFUL (11).");
							if(file_output) fputs(" NOT SUCCESSFUL (11).",file_output);
							fclose(extfile);
							return 0;
						}
					}
					else {
						if((inputms[0][x].tpast = (double *)realloc(inputms[0][x].tpast,(2+inputms[0][x].nintn)*sizeof(double))) == 0) {
							puts(" NOT SUCCESSFUL (12).");
							if(file_output) fputs(" NOT SUCCESSFUL (12).",file_output);
							fclose(extfile);
							return 0;
						}
					}
					memset(inputms[0][x].tpast,'\0',(2+inputms[0][x].nintn)*sizeof(double));
					if((fread(inputms[0][x].tpast,sizeof(double),(unsigned)(2+inputms[0][x].nintn),extfile)) != (unsigned)(2+inputms[0][x].nintn)) {
						puts("Error: it is not possible to read the data");
						if(file_output) fputs("Error: it is not possible to read the data",file_output);
						fclose(extfile);
						return 0;
					}
					
					if(inputms[0][x].npop) {
						if(x>0) {
							if((inputms[0][x].factor_pop = (double *) calloc((inputms[0][x].npop),sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL (13).");
								if(file_output) fputs(" NOT SUCCESSFUL (13).",file_output);
								fclose(extfile);
								return 0;
							}
						}
						else {
							if((inputms[0][x].factor_pop = (double *) realloc(inputms[0][x].factor_pop,(inputms[0][x].npop)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL (14).");
								if(file_output) fputs(" NOT SUCCESSFUL (14).",file_output);
								fclose(extfile);
								return 0;
							}
						}
						memset(inputms[0][x].factor_pop,'\0',(inputms[0][x].npop)*sizeof(double));
						if((fread(inputms[0][x].factor_pop,sizeof(double),(unsigned)inputms[0][x].npop,extfile)) != (unsigned)inputms[0][x].npop) {
							puts("Error: it is not possible to read the data");
							if(file_output) fputs("Error: it is not possible to read the data",file_output);
							fclose(extfile);
							return 0;
						}
						
						if(x>0) {
							if((inputms[0][x].freq = (double *)calloc((inputms[0][x].npop),sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL (15).");
								if(file_output) fputs(" NOT SUCCESSFUL (15).",file_output);
								fclose(extfile);
								return 0;
							}
						}
						else {
							if((inputms[0][x].freq = (double *)realloc(inputms[0][x].freq,(inputms[0][x].npop)*sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL (16).");
								if(file_output) fputs(" NOT SUCCESSFUL (16).",file_output);
								fclose(extfile);
								return 0;
							}
						}
						memset(inputms[0][x].freq,'\0',(inputms[0][x].npop)*sizeof(double));
						if((fread(inputms[0][x].freq,sizeof(double),(unsigned)inputms[0][x].npop,extfile)) != (unsigned)inputms[0][x].npop) {
							puts("Error: it is not possible to read the data");
							if(file_output) fputs("Error: it is not possible to read the data",file_output);
							fclose(extfile);
							return 0;
						}
					}
					else {
						if(x>0) {
							if((inputms[0][x].factor_pop = (double *) calloc(1,sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL (17).");
								if(file_output) fputs(" NOT SUCCESSFUL (17).",file_output);
								fclose(extfile);
								return 0;
							}
						}
						else {
							if((inputms[0][x].factor_pop = (double *) realloc(inputms[0][x].factor_pop,(1*sizeof(double)))) == 0) {
								puts(" NOT SUCCESSFUL (18).");
								if(file_output) fputs(" NOT SUCCESSFUL (18).",file_output);
								fclose(extfile);
								return 0;
							}
						}
						inputms[0][x].factor_pop[0] = (double)1.0;
						if(x>0) {
							if((inputms[0][x].freq = (double *)calloc(1,sizeof(double))) == 0) {
								puts(" NOT SUCCESSFUL (19).");
								if(file_output) fputs(" NOT SUCCESSFUL (19).",file_output);
								fclose(extfile);
								return 0;
							}
						}
						else {
							if((inputms[0][x].freq = (double *)realloc(inputms[0][x].freq,(1*sizeof(double)))) == 0) {
								puts(" NOT SUCCESSFUL (20).");
								if(file_output) fputs(" NOT SUCCESSFUL (20).",file_output);
								fclose(extfile);
								return 0;
							}
						}
						inputms[0][x].freq[0] = (double)1.0;
					}
				}
				if(global[0][0].montecarlo_sim == 1 && global[0][0].nitermat && global[0][0].nlocimat) {
					if((*matrixmlsim = (struct statistisimmuloc *) realloc(*matrixmlsim,(global[0][0].nitermat)*sizeof(struct statistisimmuloc))) == 0) {
						puts(" NOT SUCCESSFUL (21).");
						if(file_output) fputs(" NOT SUCCESSFUL (21).",file_output);
						fclose(extfile);
						return 0;
					}
					memset(*matrixmlsim,'\0',(global[0][0].nitermat)*sizeof(struct statistisimmuloc));
					if((fread(*matrixmlsim,sizeof(struct statistisimmuloc),(unsigned long)global[0][0].nitermat,extfile)) != (unsigned long)global[0][0].nitermat) {
						puts(" NOT SUCCESSFUL (22).");
						if(file_output) fputs(" NOT SUCCESSFUL (22).",file_output);
						fclose(extfile);
						return 0;
					}
					
					if((*avgstatloci = (struct horizontalstatsml *) realloc(*avgstatloci,(global[0][0].nlocimat)*sizeof(struct horizontalstatsml))) == 0) {
						puts(" NOT SUCCESSFUL (23).");
						if(file_output) fputs(" NOT SUCCESSFUL (23).",file_output);
						fclose(extfile);
						return 0;
					}
					memset(*avgstatloci,'\0',(global[0][0].nlocimat)*sizeof(struct horizontalstatsml));
					if((fread(*avgstatloci,sizeof(struct horizontalstatsml),(unsigned)global[0][0].nlocimat,extfile)) != (unsigned)global[0][0].nlocimat) {
						puts(" NOT SUCCESSFUL (24).");
						if(file_output) fputs(" NOT SUCCESSFUL (24).",file_output);
						fclose(extfile);
						return 0;
					}
					if(global[0][0].onlymulo == 0) {
						if((*matrixsim = (struct statistisim **)realloc(*matrixsim,global[0][0].nlocimat*sizeof(struct statistisim *))) == 0) {
							puts(" NOT SUCCESSFUL (25).");
							if(file_output) fputs(" NOT SUCCESSFUL (25).",file_output);
							fclose(extfile);
							return 0;
						}
						/*if(global[0][0].maxloc < global[0][0].nlocimat) {
							for(x=0;x<global[0][0].maxloc;x++) {
								if((matrixsim[0][x] = (struct statistisim *)realloc(matrixsim[0][x],global[0][0].nitermat*sizeof(struct statistisim))) == 0) {
									puts(" NOT SUCCESSFUL.");
									if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
									fclose(extfile);
									return 0;
								}
								if((fread(matrixsim[x],sizeof(struct statistisim),global[0][0].nitermat,extfile)) != global[0][0].nitermat) {
									puts(" NOT SUCCESSFUL.");
									if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
									fclose(extfile);
									return 0;
								}
							}*/
							for(x=0/*global[0][0].maxloc*/;x<global[0][0].nlocimat;x++) {
								if((matrixsim[0][x] = (struct statistisim *)calloc(global[0][0].nitermat,sizeof(struct statistisim))) == 0) {
									puts(" NOT SUCCESSFUL (26).");
									if(file_output) fputs(" NOT SUCCESSFUL (26).",file_output);
									fclose(extfile);
									return 0;
								}
								if((fread(matrixsim[0][x],sizeof(struct statistisim),(unsigned long)global[0][0].nitermat,extfile)) != (unsigned long)global[0][0].nitermat) {
									puts(" NOT SUCCESSFUL (27).");
									if(file_output) fputs(" NOT SUCCESSFUL (27).",file_output);
									fclose(extfile);
									return 0;
								}
							}
							/*global[0][0].maxloc = global[0][0].nlocimat;*/
						/*}
						else {
							for(x=0;x<global[0][0].nlocimat;x++) {
								if((matrixsim[0][x] = (struct statistisim *)realloc(matrixsim[0][x],global[0][0].nitermat*sizeof(struct statistisim))) == 0) {
									puts(" NOT SUCCESSFUL.");
									if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
									fclose(extfile);
									return 0;
								}
							   if((fread(matrixsim[0][x],sizeof(struct statistisim),global[0][0].nitermat,extfile)) != global[0][0].nitermat) {
									puts(" NOT SUCCESSFUL.");
									if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
									fclose(extfile);
									return 0;
								}
							}
						}*/
					}
				}
			}
            if(global[0][0].dataforsimt == 1 && global[0][0].ml0done == 1 && global[0][0].nitermatt) {
                if((lrdist[0] = (struct LRTdist *)realloc(lrdist[0],(global[0][0].nitermatt)*sizeof(struct LRTdist))) == 0) {
                    puts(" NOT SUCCESSFUL (28).");
                    if(file_output) fputs(" NOT SUCCESSFUL (28).",file_output);
					fclose(extfile);
                    return 0;
                }
				memset(lrdist[0],'\0',(global[0][0].nitermatt)*sizeof(struct LRTdist));
                if((fread(lrdist[0],sizeof(struct LRTdist),(unsigned long)global[0][0].nitermatt,extfile)) != (unsigned long)global[0][0].nitermatt) {
                    puts(" NOT SUCCESSFUL (29).");
                    if(file_output) fputs(" NOT SUCCESSFUL (29).",file_output);
					fclose(extfile);
                    return 0;
                }
				/*
				if((*mthetasim = (struct MLthetasim **)realloc(*mthetasim,global[0][0].nlocimatt*sizeof(struct MLthetasim *))) == 0) {
					puts(" NOT SUCCESSFUL.");
					if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
					fclose(extfile);
					return 0;
				}
				if(global[0][0].maxloc < global[0][0].nlocimatt) {
					for(x=0;x<global[0][0].maxloc;x++) {
						if((mthetasim[0][x] = (struct MLthetasim *)realloc(mthetasim[0][x],global[0][0].nitermatt*sizeof(struct MLthetasim))) == 0) {
							puts(" NOT SUCCESSFUL.");
							if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
							fclose(extfile);
							return 0;
						}
						if((fread(mthetasim[x],sizeof(struct MLthetasim),global[0][0].nitermatt,extfile)) != global[0][0].nitermatt) {
							puts(" NOT SUCCESSFUL.");
							if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
							fclose(extfile);
							return 0;
						}
					}
					for(x=global[0][0].maxloc;x<global[0][0].nlocimatt;x++) {
						if((mthetasim[0][x] = (struct MLthetasim *)calloc(global[0][0].nitermatt,sizeof(struct MLthetasim))) == 0) {
							puts(" NOT SUCCESSFUL.");
							if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
							fclose(extfile);
							return 0;
						}
						if((fread(mthetasim[0][x],sizeof(struct MLthetasim),global[0][0].nitermatt,extfile)) != global[0][0].nitermatt) {
							puts(" NOT SUCCESSFUL.");
							if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
							fclose(extfile);
							return 0;
						}
					}
					global[0][0].maxloc = global[0][0].nlocimatt;
				}
				else {
					for(x=0;x<global[0][0].nlocimatt;x++) {
						if((mthetasim[0][x] = (struct MLthetasim *)realloc(mthetasim[0][x],global[0][0].nitermatt*sizeof(struct MLthetasim))) == 0) {
							puts(" NOT SUCCESSFUL.");
							if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
							fclose(extfile);
							return 0;
						}
						if((*mthetasim = (struct MLthetasim **)realloc(*mthetasim,global[0][0].nlocimatt*sizeof(struct MLthetasim *))) == 0) {
							puts(" NOT SUCCESSFUL.");
							if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
							fclose(extfile);
							return 0;
						}
						if((fread(mthetasim[0][x],sizeof(struct MLthetasim),global[0][0].nitermatt,extfile)) != global[0][0].nitermatt) {
							puts(" NOT SUCCESSFUL.");
							if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
							fclose(extfile);
							return 0;
						}
					}
                }
				*/
			}
        }
        printf("\n\nSUCCESSFULL LOADING.\n\n");
        if(file_output) fputs("\n\nSUCCESSFULL LOADING.\n\n",file_output);
    }   
    
    /*observed data*/
    if(global[0][0].observed_data) {
        printf("  Number of loci: %d\n",global[0][0].n_loci);
        printf("  Outgroup lines: %s\n",global[0][0].name_outgroup);
        printf("  Excluded lines: %s\n\n",global[0][0].name_excluded);
        printf("  Included lines: %s\n\n",global[0][0].name_ingroups);
        printf("\n");
		if(global[0][0].gfffiles == 1) {
			printf("  Subset of positions included: %s\n\n",global[0][0].subset_positions);
			if(global[0][0].ifgencode == 1) {
				printf("  Last Genetic code included: %s\n\n",global[0][0].code_name);
				for(n=0;n<64;n += 16) {
					for(m=n;m<n+4;m += 1) {
						printf("  %c%c%c %c   %c%c%c %c   %c%c%c %c   %c%c%c %c\n",
							tripletsU[m+0][0],tripletsU[m+0][1],tripletsU[m+0][2],global[0][0].genetic_code[m+0],
							tripletsU[m+4][0],tripletsU[m+4][1],tripletsU[m+4][2],global[0][0].genetic_code[m+4],
							tripletsU[m+8][0],tripletsU[m+8][1],tripletsU[m+8][2],global[0][0].genetic_code[m+8],
							tripletsU[m+12][0],tripletsU[m+12][1],tripletsU[m+12][2],global[0][0].genetic_code[m+12]);
					}
				}
				printf("\n\n");
			}
		}
		if(file_output) {
            fprintf(file_output," Number of loci: %d\n",global[0][0].n_loci);
            fprintf(file_output," Outgroup lines: %s\n",global[0][0].name_outgroup);
            fprintf(file_output," Excluded lines: %s\n\n",global[0][0].name_excluded);
            fprintf(file_output," Included lines: %s\n\n",global[0][0].name_ingroups);
            fputs("\n",file_output);
			if(global[0][0].gfffiles == 1) {
				fprintf(file_output," Subset of positions included: %s\n\n",global[0][0].subset_positions);
				if(global[0][0].ifgencode == 1) {
					fprintf(file_output," Last Genetic code included: %s\n\n",global[0][0].code_name);
					for(n=0;n<64;n += 16) {
						for(m=n;m<n+4;m += 1) {
							fprintf(file_output,"  %c%c%c %c   %c%c%c %c   %c%c%c %c   %c%c%c %c\n",
								tripletsU[m+0][0],tripletsU[m+0][1],tripletsU[m+0][2],global[0][0].genetic_code[m+0],
								tripletsU[m+4][0],tripletsU[m+4][1],tripletsU[m+4][2],global[0][0].genetic_code[m+4],
								tripletsU[m+8][0],tripletsU[m+8][1],tripletsU[m+8][2],global[0][0].genetic_code[m+8],
								tripletsU[m+12][0],tripletsU[m+12][1],tripletsU[m+12][2],global[0][0].genetic_code[m+12]);
						}
					}
					fprintf(file_output,"\n\n");
				}
			}
        }
        
        for(x=0;x<global[0][0].n_loci;x++) {
            printf(" %s\tnsamples: %d\tnoutgroups: %d\n",matrix[0][x].gene,matrix[0][x].nsamples,matrix[0][x].noutgroups);
            if(file_output) 
                fprintf(file_output," %s\tnsamples: %d\tnoutgroups: %d\n",
                    matrix[0][x].gene,matrix[0][x].nsamples,matrix[0][x].noutgroups);
        }
    }
    /*estimated parameters by ML*/
    if(global[0][0].ml0done == 1) {
        printf("\nMultilocus theta estimation by ML is included.\n");
        if(file_output) fputs("\nMultilocus theta estimation by ML is included.\n",file_output);
        
        if(global[0][0].dataequalsim == 0) {
            printf("\nAttention. Parameter estimation is not corresponding to the observed data.\n");
            if(file_output) 
                fputs("\nAttention.Parameter estimation is not corresponding to the observed data.\n",
                    file_output);
        }
		print_seqinfodata_ml0(*global,*matrix/*,*matrixml*/,*datams,*inputms,file_output);
    }
    else {
        if(global[0][0].dataforsim == 1) {
            printf("\nMultilocus theta estimation results is NOT included.");
            if(file_output) fputs("\nMultilocus theta estimation results is NOT included.",file_output);
            printf("\nBUT data ready to run is included.");
            if(file_output) fputs("\nBUT data ready to run is included.",file_output);
			print_seqinfodata_ml0(*global,*matrix/*,*matrixml */,*datams,*inputms,file_output);
        }
    }
	/*simulated data*/
    if(global[0][0].montecarlo_sim == 1) {
        printf("\nMonte Carlo Coalescent simulations results are also included.\n");
        if(file_output) fputs("\nMonte Carlo Coalescent simulations results are also included.\n",file_output);
        
        if(global[0][0].dataequalsim == 0) {
            printf("\nAttention. Monte Carlo simulations are not corresponding to the observed data.\n");
            if(file_output) 
                fputs("\nAttention. Monte Carlo simulations are not corresponding to the observed data.\n",
                    file_output);
        }
        print_datasimulations(global[0][0].dataindata,*datams,*inputms,file_output,global[0][0].dataequalsim,*matrix);
    }
    else {
        if(global[0][0].dataforsim == 1) {
            printf("\nMonte Carlo Coalescent simulations results are NOT included.");
            if(file_output) fputs("\nMonte Carlo Coalescent simulations results are NOT included.",file_output);
            printf("\nBUT Monte Carlo Coalescent simulations parameters are included.");
            if(file_output) fputs("\nBUT Monte Carlo Coalescent simulations parameters are included.",file_output);
            print_seqinfodata_sim(global[0][0].dataindata,*datams,*inputms,file_output,global[0][0].dataequalsim,*matrix);
            print_evocond_sim(*datams,*inputms,file_output/*,global[0][0].dataequalsim*/,*matrix);
        }
    }
    
    return 1;
}
