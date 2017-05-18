/*
 *  file_menu.c
 *  MuLoNeTests
 *
 *  Created by sebas on Sat Feb 22 2003.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

void open_filemenu(struct globvar **global,struct statistics **matrix,struct statmulo **matrixml,struct var **datams,
    struct var2b **inputms,struct statistisim ***matrixsim, struct statistisimmuloc **matrixmlsim,
    struct horizontalstatsml **avgstatloci,struct LRTdist **lrdist/*,struct MLthetasim ***mthetasim*/,FILE **file_output)
{
    static int once = 0;
    #if COMMAND_LINE
    static char *name_folderinput; /*256 characters..*/
    static char *name_folderinputgff; /*256 characters..*/
    static char *namefile_folder;
    static char *executation;
    FILE *filelist;
    char *f;
    #endif
    char k[1],l[1]/*,kk[1]*/;
    int n;
    static char *name_fileinput;
    static char *name_fileinputfolder;
    static char *name_fileinputgff;
    static char *name_fileinputfoldergff;
    static char *name_fileinputGFF;
    static char *name_fileinputfolderGFF;
    static int excludelines = -1;
    static int includelines = -1;
        
    int open_files1(/*struct statistics **,struct statmulo **,int *,*/int * /*,int * */,int *,FILE *,char **,char ** /*,double * */,int *, char **);
    void open_files2(struct statistics **,struct statmulo **,int *,int *,int * /*,int * */,int *,
		FILE *, char **, char **,char **,char **,int ,int *,int ,struct lianinput **,double ,
		char **,char **,char **,char **,int ,char * /*,int ,char * */,char *,int *,char **);
    void open_outfile(/*struct statistics *,struct statmulo *,int *,int *,int *,*/FILE **);
    void close_outfile(/*struct statistics *,struct statmulo *,int *,int *,int *,*/FILE **);
    void erase(struct globvar **,struct statistics ** /*,struct statmulo * */,struct statistisim ***,struct statistisimmuloc **,
        struct horizontalstatsml **,int *,struct var2b **,struct LRTdist ** /*,struct MLthetasim ****/,struct var **,FILE *,int *);
    void save_simmatrix(struct globvar *,struct statistics *,struct statmulo *,struct var *, struct var2b *,
        struct statistisim **, struct statistisimmuloc *, struct horizontalstatsml *,struct LRTdist * /*,struct MLthetasim ***/,FILE *);
    int open_filematrix(struct globvar **,struct statistics **,struct statmulo **,struct var **, struct var2b **,
        struct statistisim ***, struct statistisimmuloc **, struct horizontalstatsml **,struct LRTdist ** /*,struct MLthetasim ****/,FILE *);
    int dohka_obs(struct statistics **,struct statmulo **,int *, int *);
    
    int *observed_data;
    int *outgroup;
    int *n_loci;
    int *montecarlo_sim;
    char *name_outgroup;
    char *name_excluded;
    char *name_ingroups;

    static char *lian_file;
	static int lian_in = -1;
	static struct lianinput *liandata;
	FILE *lianp;
	static double factor_chrn = (double)0.0;
	int x,y;
	
	static char *pathgff;/*256 characters*/
	char *ff,*ff1;
	
	void gff_function(int *,int *,char *,int *,char *,char *,FILE **);

    if(once == 0) { /* only initialize once. */
        once = 1;
        if ((name_fileinput = (char *)calloc(256,sizeof(char ))) == 0) {
            puts("\nError: memory not reallocated. filemenu.1 \n");
            exit(1);
        }
        if ((name_fileinputfolder = (char *)calloc(256,sizeof(char *))) == 0) {
            puts("\nError: memory not reallocated. filemenu.1 \n");
            exit(1);
        }
        if ((name_fileinputgff = (char *)calloc(256,sizeof(char ))) == 0) {
            puts("\nError: memory not reallocated. filemenu.1 \n");
            exit(1);
        }
        if ((name_fileinputfoldergff = (char *)calloc(256,sizeof(char *))) == 0) {
            puts("\nError: memory not reallocated. filemenu.1 \n");
            exit(1);
        }
        if ((name_fileinputGFF = (char *)calloc(256,sizeof(char ))) == 0) {
            puts("\nError: memory not reallocated. filemenu.1 \n");
            exit(1);
        }
        if ((name_fileinputfolderGFF = (char *)calloc(256,sizeof(char *))) == 0) {
            puts("\nError: memory not reallocated. filemenu.1 \n");
            exit(1);
        }
        if ((lian_file = (char *)calloc(256,sizeof(char ))) == 0) {
            puts("\nError: memory not reallocated. filemenu.17 \n");
            exit(1);
        }
        if ((liandata = (struct lianinput *)calloc(1,sizeof(struct lianinput))) == 0) {
            puts("\nError: memory not reallocated. filemenu.17 \n");
            exit(1);
        }
        if ((liandata[0].namesam = (char **)calloc(1,sizeof(char *))) == 0) {
            puts("\nError: memory not reallocated. filemenu.17 \n");
            exit(1);
        }
        if ((liandata[0].samhap = (int **)calloc(1,sizeof(int *))) == 0) {
            puts("\nError: memory not reallocated. filemenu.17 \n");
            exit(1);
        }
        if ((pathgff = (char *)calloc(256,sizeof(char ))) == 0) {
            puts("\nError: memory not reallocated. filemenu.1 \n");
            exit(1);
        }
        #if COMMAND_LINE
        if ((name_folderinput = (char *)calloc(256,sizeof(char ))) == 0) {
            puts("\nError: memory not reallocated. filemenu.1 \n");
            exit(1);
        }
        if ((name_folderinputgff = (char *)calloc(256,sizeof(char ))) == 0) {
            puts("\nError: memory not reallocated. filemenu.1 \n");
            exit(1);
        }
        if ((namefile_folder = (char *)calloc(256,sizeof(char ))) == 0) {
            puts("\nError: memory not reallocated. filemenu.1 \n");
            exit(1);
        }
        if ((executation = (char *)calloc(300,sizeof(char ))) == 0) {
            puts("\nError: memory not reallocated. filemenu.1 \n");
            exit(1);
        }
        #endif
		liandata[0].nsam = 0;
		liandata[0].nloci = 0;
    }
    while(1) {
        if(*file_output) fflush(*file_output);

        observed_data = &global[0][0].observed_data;
        outgroup = &global[0][0].outgroup;
        n_loci = &global[0][0].n_loci;
        montecarlo_sim = &global[0][0].montecarlo_sim;
        name_outgroup = global[0][0].name_outgroup;
        name_excluded = global[0][0].name_excluded;
        name_ingroups = global[0][0].name_ingroups;

        *name_fileinput = 0;
        n = 1;
        #if COMMAND_LINE
        n = 0;
        printf("\n\nMENU 0. File menu:\n\n");
        printf("CHOOSE:\n");
        printf(" 0 - CREATE an OUTPUT file to send the results.\n");
        printf(" 1 - CLOSE current OUTPUT file.\n");
        printf(" 2 - OPEN a set of INPUT data files (nbrf/pir or FASTA format) included in a folder.\n");
        printf(" 3 - INCLUDE single INPUT files separately.\n");
        printf(" 4 - CLOSE all included INPUT loci and ERASE information from the memory.\n");
        printf(" 5 - SAVE a file with all the information from observed and/or simulated data.\n");
        printf(" 6 - LOAD a file with all the information from observed and/or simulated data.\n");
        printf(" 7 - Display Main menu.\n\n");
		
        if(*file_output) {
            fprintf(*file_output,"\n\n     0. File menu:\n\n");
            fprintf(*file_output," 0 - CREATE an OUTPUT file to send the results.\n");
            fprintf(*file_output," 1 - CLOSE current OUTPUT file.\n");
            fprintf(*file_output," 2 - OPEN a set of INPUT data files (nbrf/pir or FASTA format) included in a folder.\n");
            fprintf(*file_output," 3 - INCLUDE single INPUT files separately.\n");
            fprintf(*file_output," 4 - CLOSE all included INPUT loci and ERASE information from the memory.\n");
            fprintf(*file_output," 5 - SAVE a file with all the information from observed and/or simulated data.\n");
            fprintf(*file_output," 6 - LOAD a file with all the information from observed and/or simulated data.\n");
            fprintf(*file_output," 7 - Display Main menu.\n\n");
        }
        #endif
        
        do *k = getc(stdin);
        while(*k<'0' || *k>'7');
        
        if(*file_output) fprintf(*file_output,"OPTION CHOSEN: %c\n\n",*k);
        
        switch(*k) {
            case '0':
                open_outfile(/* *matrix,*matrixml,observed_data,outgroup,n_loci,*/file_output);
                break;
            case '1':
                close_outfile(/* *matrix,*matrixml,observed_data,outgroup,n_loci,*/file_output);
                break;
            case '2':
            #if COMMAND_LINE
				/*
				if(global[0][0].observed_data != 0 || global[0][0].montecarlo_sim != 0 || global[0][0].dataforsim != 0) {
					printf("\n\nThis command is going to erase all the information in the memory. Do you want to do that (y/n)? ");
					do *kk = getchar();
					while(*kk!='y' && *kk!='n' && *kk!='Y' && *kk!='N');
				}
				else *kk = 'y';

				if(*kk == 'n' || *kk == 'N') break;
				else {
					global[0][0].n_loci = 0;
					global[0][0].observed_data = 0;
					global[0][0].montecarlo_sim = 0;
					global[0][0].dataforsim = 0;
					global[0][0].outgroup = -1;

					erase(global,matrix,*matrixml,matrixsim,matrixmlsim,avgstatloci,&excludelines,inputms,mlparam,datams,*file_output);
				}
				*/
                /* do 'ls' to a file named 'listfilesfolder.mlt'. Then read one by one from that folder.*/
				if(*file_output == 0) {
					printf("\n\n  You do not have an output file defined.\n It is recommended to send all results to an output file. ");
					printf("\n  Do you want to return to the File Menu and define an output file (y/n)? ");
					
					do *k = getchar();
					while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');        
					if(*k == 'y' || *k == 'Y') break;
				}

                if(global[0][0].observed_data == 0) {
					if(lian_in == -1) {
						printf("\n");
						printf("******************************************************************\n");
						printf("*          MANVa assume independence among loci.                 *\n");
						printf("* In case each locus have the same individuals in the sample,    *\n");
						printf("* it is possible to analyze the LINKAGE DISEQUILIBRIUM among loci*\n");
						printf("*          using LIAN 3.5 (Haubold and Hudson, 2000)             *\n");
						printf("*         (http://adenine.biz.fh-weihenstephan.de/lian/)         *\n");
						printf("******************************************************************\n");

						printf("MANVa can generate an input file for LIAN program \nONLY when the samples (excluding outgroups)"); 
						printf("\n have the SAME NAME across all loci.");
						printf("\n\n  Do you want to generate an input file compatible to LIAN (y/n)? ");
						
						do *k = getchar();
						while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');
						
						if(*k == 'y' || *k == 'Y') {
							lian_in = 1;
						}
						else if(*k == 'n' || *k == 'N') lian_in = 0;
					}
					if(!(open_files1(/*matrix,matrixml,observed_data,*/outgroup/*,n_loci*/,&excludelines,
					*file_output,&name_outgroup,&name_excluded/*,&factor_chrn*/,&includelines,&name_ingroups)))
						break;
				}
				else {
					printf("\n  Adding more loci to the current data:\n");
					global[0][0].dataequalsim = 0;
				}
				printf("\n\n  Correction factor for the number of chromosome copies \n  (i.e., autosomal is 1.0, X chromosome 0.75, Y chromosome 0.25): ");
				scanf(" %lf",&factor_chrn);
				
                #if UNIX_CL
                printf("\n  Path to the folder with the alignments (FASTA/NBRF formats)? ");
                /*printf("\n(Characters allowed are numbers, letters and '_#@$%%&*()?=+-\\/.')\n");*/
                *name_folderinput = 0;
                scanf(" %[^\n]",name_folderinput);
                strncat(name_folderinput,"/\0",256);
                /**/
                for(n=2;n<256;n++) 
                    if(name_folderinput[n] == '\0') {
                        if(name_folderinput[n-1] == '/' && name_folderinput[n-2] == '/') {
                            name_folderinput[n-1] = '\0';
							break;
						}
                        if(name_folderinput[n-1] == '/' && name_folderinput[n-2] == ' ') {
                            name_folderinput[n-2] = '/';
                            name_folderinput[n-1] = '\0';
							break;
						}
                        if(name_folderinput[n-1] == '/' && name_folderinput[n-2] == '/' && name_folderinput[n-3] == ' ') {
                            name_folderinput[n-3] = '/';
                            name_folderinput[n-2] = '\0';
							break;
						}
                    }
                /**/
                #elif DOS_CL
                printf("\n  Path to the folder with the alignments (FASTA/NBRF formats)? ");
                /*printf("\n(Characters allowed are numbers, letters and '_#@$%%&*()?=+-\\/.')\n");*/
                *name_folderinput = 0;
                scanf(" %[^\n]",name_folderinput);
				if(name_folderinput[0] == 34){
					memmove(name_folderinput,name_folderinput+1,256);
					n = (int)strcspn(name_folderinput,"\"\0");
					name_folderinput[n] = '\\';
					name_folderinput[n+1] = '\0';
				}
				else strncat(name_folderinput,"\\",256);
                /**/
				if(name_folderinput[0] == '/' || name_folderinput[0] == '\\') {
					x=0;
					while(name_folderinput[x] != '\0' && x < 254) x++;
					name_folderinput[x+1] = '\0';
					for(n=x;n>0;n--) {
						name_folderinput[n] = name_folderinput[n-1];					
					}
					name_folderinput[0] = '.';
				}
                for(n=0;n<256;n++) {
                    if(name_folderinput[n] == '/') {
                        name_folderinput[n] = '\\';
                    }
					if(name_folderinput[n] == '\0') break;
				}
                for(n=2;n<256;n++) {
                    if(name_folderinput[n] == '\0') {
                        if(name_folderinput[n-1] == '\\' && name_folderinput[n-2] == '\\') 
                            name_folderinput[n-1] = '\0';
                        break;
                    }
				}
                /**/
                #endif
                executation[0] = 0;
                #if UNIX_CL
                strncat(executation,"/bin/ls ",20);
                #elif DOS_CL
                strncat(executation,"dir /A:-D /B ",20);
                #endif
				#if DOS_CL 
				strncat(executation,"\"",256);
				#endif
				strncat(executation,name_folderinput,256);
				#if DOS_CL
				strncat(executation,"\"",256);
				#endif
				strncat(executation," > ",3);
                strncat(executation,"listfilesfolder.mlt",20);
                if(system(executation)) {
					printf("\n Error in path file or writing the file 'listfilesfolder.mlt' in the directory. ");
					if(*n_loci == 0) {
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
						excludelines = -1;
						includelines = -1;

						liandata[0].nsam = 0;
						liandata[0].nloci = 0;
						lian_in = -1;
						
						global[0][0].gfffiles = 0;
						global[0][0].subset_positions[0] = '\0';
						global[0][0].ifgencode = 0;
					}
					break;
				}
				
				/*ask for gff files and include information*/
				gff_function(&global[0][0].gfffiles,&global[0][0].observed_data,global[0][0].subset_positions,&global[0][0].ifgencode,global[0][0].code_name,global[0][0].genetic_code,file_output);
				/*ask for path and keep it (pathgff)*/
				#if GFF_ACTIVE && UNIX_CL
				if(global[0][0].gfffiles == 1) {
					printf("\n  Path to the folder with .gff/.GFF files? ");
					printf("\n  (The files in the folder must have the same name\n  than the alignent files but extension .gff or .GFF) ");
					*name_folderinputgff = 0;
					scanf(" %[^\n]",name_folderinputgff);
					strncat(name_folderinputgff,"/",256);
					/**/
					for(n=2;n<256;n++) {
						if(name_folderinputgff[n] == '\0') {
							if(name_folderinputgff[n-1] == '/' && name_folderinputgff[n-2] == '/') 
								name_folderinputgff[n-1] = '\0';
							break;
						}
					}
				}
                #elif GFF_ACTIVE && DOS_CL
				if(global[0][0].gfffiles == 1) {
					printf("\n  Path to the folder with .gff/.GFF files? ");
					printf("\n  (The files in the folder must have the same name\n  than the alignent files but extension .gff or .GFF) ");
					*name_folderinputgff = 0;
					scanf(" %[^\n]",name_folderinputgff);
					if(name_folderinputgff[0] == 34){
						memmove(name_folderinputgff,name_folderinputgff+1,256);
						n = (int)strcspn(name_folderinputgff,"\"\0");
						name_folderinputgff[n] = '\\';
						name_folderinputgff[n+1] = '\0';
					}
					else strncat(name_folderinputgff,"\\",256);
					if(name_folderinputgff[0] == '/' || name_folderinputgff[0] == '\\') {
						x=0;
						while(name_folderinputgff[x] != '\0' && x < 254) x++;
						name_folderinputgff[x+1] = '\0';
						for(n=x;n>0;n--) {
							name_folderinputgff[n] = name_folderinputgff[n-1];					
						}
						name_folderinputgff[0] = '.';
					}
					for(n=0;n<256;n++) {
						if(name_folderinputgff[n] == '/') {
							name_folderinputgff[n] = '\\';
						}
						if(name_folderinputgff[n] == '\0') break;
					}
					for(n=2;n<256;n++) {
						if(name_folderinputgff[n] == '\0') {
							if(name_folderinputgff[n-1] == '\\' && name_folderinputgff[n-2] == '\\') 
								name_folderinputgff[n-1] = '\0';
							break;
						}
					}
				}
                #endif
				
                /*open and read every file.*/
                if((filelist = fopen ("listfilesfolder.mlt","r")) == 0) {
                    puts("\n  It is not possible to open the file listfilesfolder.mlt");
                    if(*file_output) fputs("\n  It is not possible to open the file listfilesfolder.mlt",*file_output);
                    return;
                }
                if(!(f = (char *)malloc(BUFSIZ))) {
                    puts("\n  Error: memory not reallocated. file_menu.1 \n");
					if(*n_loci == 0) {
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
						excludelines = -1;
						includelines = -1;

						liandata[0].nsam = 0;
						liandata[0].nloci = 0;
						lian_in = -1;
						
						global[0][0].gfffiles = 0;
						global[0][0].subset_positions[0] = '\0';
						global[0][0].ifgencode = 0;
					}
                    return;
                }
                setbuf(filelist,f);

                name_fileinputfolder[0] = 0;
                name_fileinput[0] = 0;

				printf("\n In case reading codons:\n");
				printf("   GU: Gap/uncertainty at one position or codon, not analyzed.\n");
				printf("   MM: Multiple mutations at one position or codon, not analyzed.\n");
				printf("   SC: Stop codon, not analyzed.\n");
				if(*file_output) {
					fprintf(*file_output,"\n  In case reading codons:\n");
					fprintf(*file_output,"   GU: Gap/uncertainty at one position or codon, not analyzed.\n");
					fprintf(*file_output,"   MM: Multiple mutations at one position or codon, not analyzed.\n");
					fprintf(*file_output,"   SC: Stop codon, not analyzed.\n");
				}
				while(1) {
                    if(!(fscanf(filelist,"%s",namefile_folder))) {
                        puts("\n  It is not possible to read the file listfilesfolder.mlt");
                        if(*file_output) 
                            fputs("\n  It is not possible to read the file listfilesfolder.mlt",*file_output);
                        free(f);
                        fclose(filelist);
						if(*n_loci == 0) {
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
							excludelines = -1;
							includelines = -1;

							liandata[0].nsam = 0;
							liandata[0].nloci = 0;
							lian_in = -1;
							
							global[0][0].gfffiles = 0;
							global[0][0].subset_positions[0] = '\0';
							global[0][0].ifgencode = 0;
						}
                        return;
                    }
                    if(feof(filelist)) break;
                    /*realloc *matrix*/
                    if((*matrix = (struct statistics *) realloc(*matrix,((*n_loci)+1)*sizeof(struct statistics))) == 0) {
                        puts("  Matrix not reallocated! NOT SUCCESSFUL.");
                        if(*file_output) fputs("  Matrix not reallocated! NOT SUCCESSFUL.",*file_output);
						if(*n_loci == 0) {
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
							excludelines = -1;
							includelines = -1;

							liandata[0].nsam = 0;
							liandata[0].nloci = 0;
							lian_in = -1;
							
							global[0][0].gfffiles = 0;
							global[0][0].subset_positions[0] = '\0';
							global[0][0].ifgencode = 0;
						}
                        return;
                    }
					memset((*matrix)+(*n_loci),'\0',sizeof(struct statistics));
                    
					strncat(name_fileinputfolder,name_folderinput,256);
					strncat(name_fileinputfolder,namefile_folder,256);
                    strncat(name_fileinput,namefile_folder,256);
					
					if(global[0][0].gfffiles == 1) {
						name_fileinputfoldergff[0] = 0;
						name_fileinputgff[0] = 0;
						name_fileinputfolderGFF[0] = 0;
						name_fileinputGFF[0] = 0;
						
						strncat(name_fileinputgff,namefile_folder,256);
						ff = strstr(name_fileinputgff,".");
						ff1 = ff;
						while(ff1 != NULL) {/*search the last '.'*/
							if((ff1 = strstr(ff+1,".")) != NULL)
								ff = ff1;
						}
						if(ff) {
							*ff = '\0';
							strcat(ff,".gff");
						}
						else strncat(name_fileinputgff,".gff",256);
						strncat(name_fileinputfoldergff,name_folderinputgff,256);
						strncat(name_fileinputfoldergff,name_fileinputgff,256);
						strncat(name_fileinputGFF,namefile_folder,256);
						ff = strstr(name_fileinputGFF,".");
						if(ff) {
							*ff = '\0';
							strcat(ff,".GFF");
						}
						else strncat(name_fileinputGFF,".GFF",256);
						strncat(name_fileinputfolderGFF,name_folderinputgff,256);
						strncat(name_fileinputfolderGFF,name_fileinputGFF,256);
					}
					
                    open_files2(matrix,matrixml,observed_data,outgroup,n_loci/*,montecarlo_sim*/,&excludelines,
                                *file_output,&name_outgroup,&name_excluded,&name_fileinputfolder,&name_fileinput,1,
                                &global[0][0].dataforsim,lian_in,&liandata,factor_chrn,
								&name_fileinputfoldergff,&name_fileinputgff,&name_fileinputfolderGFF,&name_fileinputGFF,
								global[0][0].gfffiles,global[0][0].subset_positions/*,global[0][0].ifgencode,global[0][0].code_name*/,
								global[0][0].genetic_code,&includelines,&name_ingroups);
                    name_fileinputfolder[0] = 0;
                    name_fileinput[0] = 0;
                }
                free(f);
                fclose(filelist);

                remove("listfilesfolder.mlt");
                
                printf("\n\n  Total number of loci introduced: %d.",*n_loci);
                if(*file_output) fprintf(*file_output,"\n  Total number of loci introduced: %d.",*n_loci);
    
                if(*n_loci == 0) {
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
                    excludelines = -1;
					includelines = -1;

					liandata[0].nsam = 0;
					liandata[0].nloci = 0;
					lian_in = -1;
					
					global[0][0].gfffiles = 0;
					global[0][0].subset_positions[0] = '\0';
					global[0][0].ifgencode = 0;
                }
                else {
					/*global[0][0].dataindata = 0;*/
				}
                /*
				    printf("\n\n Do you want to save all the loaded information in a MANVa file (y/n)? ");
                    do {
                        *l = getc(stdin);
                    }while(*l !='y' && *l!='Y' && *l !='n' && *l!='N');
                    if(*l == 'n' || *l == 'n') break;
                    if(*outgroup == 1 && *n_loci > 1) {

                        if(!(dohka_obs(matrix,matrixml,outgroup,n_loci))) {
                            puts("\n Error: HKA could not be calculated, sorry.");
                            if(*file_output) fputs("\n Error: HKA could not be calculated, sorry.",*file_output);
                        }
                    }
                    save_simmatrix(*global,*matrix,*matrixml,*datams,*inputms,*matrixsim,*matrixmlsim,*avgstatloci,*lrdist,*file_output);
                }
				*/
                break;
            #endif
            case '3':
				if(*file_output == 0) {
					printf("\n\n  You do not have an output file defined. \n  It is recommended to send all results to an output file. ");
					printf("\n  Do you want to return to the File Menu and define an output file (y/n)? ");
					
					do *k = getchar();
					while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');        
					if(*k == 'y' || *k == 'Y') break;
				}

                if(global[0][0].observed_data == 0) {
					if(lian_in == -1) {
						printf("\n");
						printf("******************************************************************\n");
						printf("*          MANVa assume independence among loci.                 *\n");
						printf("* In case each locus have the same individuals in the sample,    *\n");
						printf("* it is possible to analyze the LINKAGE DISEQUILIBRIUM among loci*\n");
						printf("*          using LIAN 3.5 (Haubold and Hudson, 2000)             *\n");
						printf("*         (http://adenine.biz.fh-weihenstephan.de/lian/)         *\n");
						printf("******************************************************************\n");

						printf("MANVa can generate an input file for LIAN program\n ONLY when the samples (excluding outgroups)"); 
						printf("\n have the SAME NAME across all loci.");
						printf("\n\n  Do you want to generate an input file compatible to LIAN (y/n)? ");
						
						do *k = getchar();
						while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');
						
						if(*k == 'y' || *k == 'Y') {
							lian_in = 1;
						}
						else if(*k == 'n' || *k == 'N') lian_in = 0;
					}
					if(!(open_files1(/*matrix,matrixml,observed_data,*/outgroup/*,n_loci*/,&excludelines,
					*file_output,&name_outgroup,&name_excluded/*,&factor_chrn*/,&includelines,&name_ingroups)))
						break;
				}
				else {
					printf("\n  Adding more loci to the current data:\n");
					global[0][0].dataequalsim = 0;
				}
				
				printf("\n\n  Correction factor for the number of chromosome copies\n (i.e., autosomal is 1.0, X chromosome 0.75, Y chromosome 0.25): ");
				scanf(" %lf",&factor_chrn);
				
				/*ask for gff files and include information*/
				gff_function(&global[0][0].gfffiles,&global[0][0].observed_data,global[0][0].subset_positions,&global[0][0].ifgencode,global[0][0].code_name,global[0][0].genetic_code,file_output);

                if((*matrix = (struct statistics *) realloc(*matrix,((*n_loci)+1)*sizeof(struct statistics))) == 0) {
                    puts("  Error: Matrix not reallocated. NOT SUCCESSFUL.");
                    if(*file_output) fputs("  Error: Matrix not reallocated. NOT SUCCESSFUL.",*file_output);
                    return;
                }
				memset((*matrix)+(*n_loci),'\0',sizeof(struct statistics));
                name_fileinputfolder[0] = 0;
                name_fileinput[0] = 0;

				printf("\n In case reading codons:\n");
				printf("   SC: Stop codon, not analyzed.\n");
				printf("   GU: Gap/uncertainty at codon, not analyzed.\n");
				printf("   MM: Multiple mutations at codon, not analyzed.\n");
				if(*file_output) {
					fprintf(*file_output,"\n  In case reading codons:\n");
					fprintf(*file_output,"   SC: Stop codon, not analyzed.\n");
					fprintf(*file_output,"   GU: Gap/uncertainty at codon, not analyzed.\n");
					fprintf(*file_output,"   MM: Multiple mutations at codon, not analyzed.\n");
				}

                open_files2(matrix,matrixml,observed_data,outgroup,n_loci/*,montecarlo_sim*/,&excludelines,
                    *file_output,&name_outgroup,&name_excluded,&name_fileinputfolder,&name_fileinput,0,
                    &global[0][0].dataforsim,lian_in,&liandata,factor_chrn,
					&name_fileinputfoldergff,&name_fileinputgff,&name_fileinputfolderGFF,&name_fileinputGFF,
					global[0][0].gfffiles,global[0][0].subset_positions/*,global[0][0].ifgencode,global[0][0].code_name*/,
					global[0][0].genetic_code,&includelines,&name_ingroups);

                if(*n_loci == 0) {
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

					liandata[0].nsam = 0;
					liandata[0].nloci = 0;
					lian_in = -1;

					global[0][0].gfffiles = 0;
					global[0][0].subset_positions[0] = '\0';
					global[0][0].ifgencode = 0;
                }
                else {
					/*global[0][0].dataindata = 0;*/
				}
				/*
                    printf("\n\n Do you want to save all the loaded information in a MANVa file (y/n)? ");
                    do {
                        *l = getc(stdin);
                    }while(*l !='y' && *l!='Y' && *l !='n' && *l!='N');
                    if(*l == 'n' || *l == 'n') break;
                    if(*outgroup == 1 && *n_loci > 1) {
                        if(!(dohka_obs(matrix,matrixml,outgroup,n_loci))) {
                            puts("\n Error: HKA could not be calculated, sorry.");
                            if(*file_output) fputs("\n Error: HKA could not be calculated, sorry.",*file_output);
                        }
                    }
                    save_simmatrix(*global,*matrix,*matrixml,*datams,*inputms,*matrixsim,*matrixmlsim,*avgstatloci,*lrdist,*file_output);
                }
				*/
                break;
            case '4':
				if(lian_in == 1) {
					printf("\n  LIAN file generation aborted.");
					if(*file_output) fprintf(*file_output,"\n LIAN file generation aborted.");
					
					for(x=0;x<liandata[0].nsam;x++) {
						free(liandata[0].namesam[x]);
						free(liandata[0].samhap[x]);
					}
					if((liandata[0].namesam = (char **)realloc(liandata[0].namesam,1*sizeof(char *))) == 0) {
						puts("  Matrix not reallocated! NOT SUCCESSFUL.");
						if(*file_output) fputs("  Matrix not reallocated! NOT SUCCESSFUL.",*file_output);
						return;
					}
					if((liandata[0].samhap = (int **)realloc(liandata[0].samhap,1*sizeof(int *))) == 0) {
						puts("  Matrix not reallocated! NOT SUCCESSFUL.");
						if(*file_output) fputs("  Matrix not reallocated! NOT SUCCESSFUL.",*file_output);
						return;
					}
				}
				liandata[0].nsam = 0;
				liandata[0].nloci = 0;
				lian_in = -1;

                erase(global,matrix/*,*matrixml*/,matrixsim,matrixmlsim,avgstatloci,&excludelines,inputms,lrdist/*,mthetasim*/,datams,*file_output,&includelines);
                break;
            case '5':
                /* before (re)calculate hka.*/
                if(*outgroup == 1 && *n_loci > 1) {
                    if(!(dohka_obs(matrix,matrixml,outgroup,n_loci))) {
                        puts("\n  Error: HKA could not be calculated, sorry.");
                        if(*file_output) fputs("\n Error: HKA could not be calculated, sorry.",*file_output);
                    }
                }
				/*Save the file LIAN in case the user define it*/
				if(lian_in == 1 && liandata[0].nloci > 1 && liandata[0].nsam > 2) {
					printf("\n\n  Please indicate the name of LIAN input file: ");
					do{
						lian_file[0] = 0;
						scanf(" %[^\n]",lian_file);
					} while(lian_file[0] == 0);
					printf("\n  Saving LIAN file.\n\n");
					if(*file_output) fprintf(*file_output,"\n  Saving LIAN file.\n\n");
					
					if((lianp = fopen (lian_file,"w")) == 0) {
						puts("\n  Error in input/output. LIAN file not saved.");
						if(*file_output) fputs("\n  Error in input/output. LIAN file not saved.",*file_output);
					}
					else {
						fprintf(lianp,"%d %d\n",liandata[0].nsam,liandata[0].nloci);
						for(x=0;x<liandata[0].nsam;x++) {
							fprintf(lianp,"%s",liandata[0].namesam[x]);
							for(y=0;y<liandata[0].nloci;y++) {
								fprintf(lianp," %d",liandata[0].samhap[x][y]);
							}
						}
						fprintf(lianp,"# \n# %d samples in %d loci. \n",liandata[0].nsam,liandata[0].nloci);
						fprintf(lianp,"# \n# Names of the loci included: \n");
						if(global[0][0].subset_positions[0]) fprintf(lianp,"# Subset of positions considered:  %s",global[0][0].subset_positions);
						for(y=0;y<liandata[0].nloci;y++) {
							fprintf(lianp,"# %d:%-20s\n",y,matrix[0][y].gene);
						}
						fprintf(lianp,"# \n# Input file generated with ");
						fprintf(lianp,MANVa);
						fclose(lianp);
					}
				}
				else {
					if(lian_in == 1 && (liandata[0].nloci < 2 || liandata[0].nsam < 3)) {
						printf("\n  Sorry. Not enough information to generate a LIAN file.\n\n");
						if(*file_output) fprintf(*file_output,"\n S orry. Not enough information to generate a LIAN file.\n\n");
					}
				}
				if(lian_in == 1) {
					for(x=0;x<liandata[0].nsam;x++) {
						free(liandata[0].namesam[x]);
						free(liandata[0].samhap[x]);
					}
					if((liandata[0].namesam = (char **)realloc(liandata[0].namesam,1*sizeof(char *))) == 0) {
						puts("  Matrix not reallocated! NOT SUCCESSFUL.");
						if(*file_output) fputs("  Matrix not reallocated! NOT SUCCESSFUL.",*file_output);
						return;
					}
					if((liandata[0].samhap = (int **)realloc(liandata[0].samhap,1*sizeof(int *))) == 0) {
						puts("  Matrix not reallocated! NOT SUCCESSFUL.");
						if(*file_output) fputs("  Matrix not reallocated! NOT SUCCESSFUL.",*file_output);
						return;
					}
				} 
				liandata[0].nsam = 0;
				liandata[0].nloci = 0;
				lian_in = -1;

                save_simmatrix(*global,*matrix,*matrixml,*datams,*inputms,*matrixsim,*matrixmlsim,*avgstatloci,*lrdist/*,*mthetasim*/,*file_output);
                break;
            case '6':
				if(lian_in == 1) {
					printf("\n  LIAN file generation aborted.");
					if(*file_output) fprintf(*file_output,"\n  LIAN file generation aborted.");
					
					for(x=0;x<liandata[0].nsam;x++) {
						free(liandata[0].namesam[x]);
						free(liandata[0].samhap[x]);
					}
					if((liandata[0].namesam = (char **)realloc(liandata[0].namesam,1*sizeof(char *))) == 0) {
						puts("  Matrix not reallocated! NOT SUCCESSFUL.");
						if(*file_output) fputs("  Matrix not reallocated! NOT SUCCESSFUL.",*file_output);
						return;
					}
					if((liandata[0].samhap = (int **)realloc(liandata[0].samhap,1*sizeof(int *))) == 0) {
						puts("  Matrix not reallocated! NOT SUCCESSFUL.");
						if(*file_output) fputs("  Matrix not reallocated! NOT SUCCESSFUL.",*file_output);
						return;
					}

					lian_in = -1;
				}
                if(open_filematrix(global,matrix,matrixml,datams,inputms,matrixsim,matrixmlsim,avgstatloci,lrdist/*,mthetasim*/,*file_output) == 0) {
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

					global[0][0].gfffiles = 0;
					global[0][0].subset_positions[0] = '\0';
					global[0][0].ifgencode = 0;
                }
				liandata[0].nsam = 0;
				liandata[0].nloci = 0;
				lian_in = -1;
                break;
            case '7':
                /* go to main menu, but before that, calculate hka.*/
				while (*n_loci > 0) {
					printf("\n\n  Do you want to save all the loaded information in a MANVa file (y/n)? ");
					do {
						*l = getc(stdin);
					}while(*l !='y' && *l!='Y' && *l !='n' && *l!='N');
					if(*outgroup == 1 && *n_loci > 1) {
						if(!(dohka_obs(matrix,matrixml,outgroup,n_loci))) {
							puts("\n  Error: HKA could not be calculated, sorry.");
							if(*file_output) fputs("\n  Error: HKA could not be calculated, sorry.",*file_output);
						}
					}
					if(*l == 'N' || *l == 'n') break;
					else save_simmatrix(*global,*matrix,*matrixml,*datams,*inputms,*matrixsim,*matrixmlsim,*avgstatloci,*lrdist/*,mthetasim*/,*file_output);
					break;
				}
				/*Save the file LIAN in case the user define it*/
				if(lian_in == 1 && liandata[0].nloci > 1 && liandata[0].nsam > 2) {
					printf("\n\n  Please indicate the name of LIAN input file: ");
					do{
						lian_file[0] = 0;
						scanf(" %[^\n]",lian_file);
					} while(lian_file[0] == 0);

					printf("\n  Saving LIAN file.\n\n");
					if(*file_output) fprintf(*file_output,"\n  Saving LIAN file.\n\n");
					
					if((lianp = fopen (lian_file,"w")) == 0) {
						puts("\n  Error in input/output. LIAN file not saved.");
						if(*file_output) fputs("\n  Error in input/output. LIAN file not saved.",*file_output);
					}
					else {
						/*fprintf(lianp,"%d %d\n",liandata[0].nsam,liandata[0].nloci);*/
						for(x=0;x<liandata[0].nsam;x++) {
							fprintf(lianp,"%s",liandata[0].namesam[x]);
							for(y=0;y<liandata[0].nloci;y++) {
								fprintf(lianp," %d",liandata[0].samhap[x][y]);
							}
							fprintf(lianp,"\n");
						}
						fprintf(lianp,"# \n# %d shared samples in %d loci. \n",liandata[0].nsam,liandata[0].nloci);
						fprintf(lianp,"# \n# Names of the loci included: \n");
						if(global[0][0].subset_positions[0]) fprintf(lianp,"# Subset of positions considered:  %s",global[0][0].subset_positions);
						for(y=0;y<liandata[0].nloci;y++) {
							fprintf(lianp,"# %d:%-20s\n",y,matrix[0][y].gene);
						}
						fprintf(lianp,"# \n# Input file generated with ");
						fprintf(lianp,MANVa);
						fclose(lianp);
					}
				}
				else {
					if(lian_in != 1) {
						liandata[0].nsam = 0;
						liandata[0].nloci = 0;
						lian_in = -1;
						return;
					}
					else {
						if(liandata[0].nloci < 2 || liandata[0].nsam < 3) {
							printf("\n  Sorry. Not enough information to generate a LIAN file.\n\n");
							if(*file_output) fprintf(*file_output,"\n  Sorry. Not enough information to generate a LIAN file.\n\n");
						}
					}
				}
				if(lian_in == 1) {
					for(x=0;x<liandata[0].nsam;x++) {
						free(liandata[0].namesam[x]);
						free(liandata[0].samhap[x]);
					}
					if((liandata[0].namesam = (char **)realloc(liandata[0].namesam,1*sizeof(char *))) == 0) {
						puts("  Matrix not reallocated! NOT SUCCESSFUL.");
						if(*file_output) fputs("  Matrix not reallocated! NOT SUCCESSFUL.",*file_output);
						return;
					}
					if((liandata[0].samhap = (int **)realloc(liandata[0].samhap,1*sizeof(int *))) == 0) {
						puts("  Matrix not reallocated! NOT SUCCESSFUL.");
						if(*file_output) fputs("  Matrix not reallocated! NOT SUCCESSFUL.",*file_output);
						return;
					}
				} 
				liandata[0].nsam = 0;
				liandata[0].nloci = 0;
				lian_in = -1;
				
                return;
                break;
        }
    }
}

void erase(struct globvar **global,struct statistics **matrix/*,struct statmulo *matrixml*/,struct statistisim ***matrixsim, struct statistisimmuloc **matrixmlsim,
           struct horizontalstatsml **avgstatloci,int *excludelines,struct var2b **inputms,struct LRTdist **lrdist/*,struct MLthetasim ***mthetasim*/,
		   struct var **datams,FILE *file_output,int *includelines)
{
    char k[1];
    long int i;
    int erase = 0;
    
    if(global[0][0].observed_data != 0 || global[0][0].montecarlo_sim != 0 || global[0][0].dataforsim != 0 || global[0][0].dataforsimt != 0) {
        printf("\n\n  Erase menu:\n");    
        printf("\n  Are you sure you want to erase all the information\n from the memory of the program (y/n)? ");
        do *k = getchar();
        while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');
        erase = 1;
    }
    else *k = 'y';
    
    if(*k == 'y' || *k == 'Y') {
        /*reallloc *matrix to 1*/
        if((*matrix = (struct statistics *) realloc(*matrix,1*sizeof(struct statistics))) == 0) {
            puts(" Matrix not reallocated! NOT SUCCESSFUL.");
            if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
            return;
        }
		memset(*matrix,'\0',sizeof(struct statistics));
		if(global[0][0].dataindata > 0) {
			inputms[0][0].nintn = 0;
			if((inputms[0][0].nrec = (double *)realloc(inputms[0][0].nrec,2*sizeof(double))) == 0) {
				puts(" Value not reallocated! NOT SUCCESSFUL.");
				if(file_output) fputs(" Value not reallocated! NOT SUCCESSFUL.",file_output);
				return;
			}
			memset(inputms[0][0].nrec,'\0',2*sizeof(double));
			if((inputms[0][0].npast = (double *)realloc(inputms[0][0].npast,2*sizeof(double))) == 0) {
				puts(" Value not reallocated! NOT SUCCESSFUL.");
				if(file_output) fputs(" Value not reallocated! NOT SUCCESSFUL.",file_output);
				return;
			}
			memset(inputms[0][0].npast,'\0',2*sizeof(double));
			if((inputms[0][0].tpast = (double *)realloc(inputms[0][0].tpast,2*sizeof(double))) == 0) {
				puts(" Value not reallocated! NOT SUCCESSFUL.");
				if(file_output) fputs(" Value not reallocated! NOT SUCCESSFUL.",file_output);
				return;
			}
			memset(inputms[0][0].tpast,'\0',2*sizeof(double));
			if((inputms[0][0].freq = (double *)realloc(inputms[0][0].freq,1*sizeof(double))) == 0) {
				puts(" Value not reallocated! NOT SUCCESSFUL.");
				if(file_output) fputs(" Value not reallocated! NOT SUCCESSFUL.",file_output);
				return;
			}
			memset(inputms[0][0].freq,'\0',sizeof(double));
			if((inputms[0][0].factor_pop = (double *)realloc(inputms[0][0].factor_pop,1*sizeof(double))) == 0) {
				puts(" Value not reallocated! NOT SUCCESSFUL.");
				if(file_output) fputs(" Value not reallocated! NOT SUCCESSFUL.",file_output);
				return;
			}
			memset(inputms[0][0].factor_pop,'\0',sizeof(double));
			
			inputms[0][0].nrec[0] = (double)1;
			inputms[0][0].npast[0] = (double)1;
			inputms[0][0].tpast[0] = (double)1;
			inputms[0][0].freq[0] = (double)1;
			inputms[0][0].factor_pop[0] = (double)1;
			
			for(i=1;i<global[0][0].nlocimat;i++) {
				free(inputms[0][i].nrec);
				free(inputms[0][i].npast);
				free(inputms[0][i].tpast);
				free(inputms[0][i].factor_pop);
				free(inputms[0][i].freq);
			}
			
			if((*inputms = (struct var2b *)realloc(*inputms,1*sizeof(struct var2b))) == 0) {
				puts(" Matrix not reallocated! NOT SUCCESSFUL.");
				if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
				return;
			}
			memset(*inputms,'\0',sizeof(struct var2b));
			memset(*datams,'\0',sizeof(struct var));
			
			if(global[0][0].montecarlo_sim) {
				if((*matrixmlsim = (struct statistisimmuloc *)realloc(*matrixmlsim,1*sizeof(struct statistisimmuloc))) == 0) {
					puts(" Matrix not reallocated! NOT SUCCESSFUL.");
					if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
					return;
				}
				memset(*matrixmlsim,'\0',sizeof(struct statistisimmuloc));
				if((*avgstatloci = (struct horizontalstatsml *)realloc(*avgstatloci,1*sizeof(struct horizontalstatsml))) == 0) {
					puts(" Matrix not reallocated! NOT SUCCESSFUL.");
					if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
					return;
				}
				memset(*avgstatloci,'\0',sizeof(struct horizontalstatsml));
				if(global[0][0].onlymulo == 0) {
					for(i=0;i<global[0][0].nlocimat;i++) 
						free(matrixsim[0][i]);
					if((*matrixsim = (struct statistisim **)realloc(*matrixsim,1*sizeof(struct statistisim *))) == 0) {
						puts(" Matrix not reallocated! NOT SUCCESSFUL.");
						if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
						return;
					}
				}
			}
			if(global[0][0].dataforsimt == 1) {
				if(global[0][0].ml0done == 1) {
					if((*lrdist = (struct LRTdist *)realloc(*lrdist,1*sizeof(struct LRTdist))) == 0) {
						puts(" Matrix not reallocated! NOT SUCCESSFUL.");
						if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
						return;
					}
					memset(*lrdist,'\0',sizeof(struct LRTdist));
					/*
					for(i=0;i<global[0][0].nlocimatt;i++) 
						free(mthetasim[0][i]);
					if((*mthetasim = (struct MLthetasim **)realloc(*mthetasim,1*sizeof(struct MLthetasim *))) == 0) {
						puts(" Matrix not reallocated! NOT SUCCESSFUL.");
						if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
						return;
					}
					*/
				}
			}
		}
        global[0][0].n_loci = 0;
        global[0][0].observed_data = 0;
        global[0][0].outgroup = -1;
        /*global[0][0].maxloc = 0;*/
        *excludelines = -1;
		*includelines = -1;
        
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
		
		global[0][0].gfffiles = 0;
		global[0][0].ifgencode = 0;
		
        if(erase == 1) {
            printf("\n\nData in the memory has been erased.\n\n");
            if(file_output) 
				fprintf(file_output,"\n\nData in the memory has been erased.\n\n");
        }
        return;
    }
    else if(*k == 'n' || *k == 'N') return;
        else exit(1);
	return;
}

void erase_results(struct globvar **global,struct statistisim ***matrixsim, struct statistisimmuloc **matrixmlsim,struct horizontalstatsml **avgstatloci/*,struct var2b **inputms,struct LRTdist **lrdist*/,struct var **data,FILE *file_output)
{
    long int i;
    
    data[0][0].n_iter = (long unsigned)0;
    data[0][0].seed1 = (long int)0;
    data[0][0].tlimit = (double)0;
	
	if(global[0][0].dataindata > 0) {
		/*
		inputms[0][0].nintn = 0;
		if((inputms[0][0].nrec = (double *)realloc(inputms[0][0].nrec,2*sizeof(double))) == 0) {
			puts(" Value not reallocated! NOT SUCCESSFUL.");
			if(file_output) fputs(" Value not reallocated! NOT SUCCESSFUL.",file_output);
			return;
		}
		memset(inputms[0][0].nrec,'\0',2*sizeof(double));
		if((inputms[0][0].npast = (double *)realloc(inputms[0][0].npast,2*sizeof(double))) == 0) {
			puts(" Value not reallocated! NOT SUCCESSFUL.");
			if(file_output) fputs(" Value not reallocated! NOT SUCCESSFUL.",file_output);
			return;
		}
		memset(inputms[0][0].npast,'\0',2*sizeof(double));
		if((inputms[0][0].tpast = (double *)realloc(inputms[0][0].tpast,2*sizeof(double))) == 0) {
			puts(" Value not reallocated! NOT SUCCESSFUL.");
			if(file_output) fputs(" Value not reallocated! NOT SUCCESSFUL.",file_output);
			return;
		}
		memset(inputms[0][0].tpast,'\0',2*sizeof(double));
		if((inputms[0][0].freq = (double *)realloc(inputms[0][0].freq,1*sizeof(double))) == 0) {
			puts(" Value not reallocated! NOT SUCCESSFUL.");
			if(file_output) fputs(" Value not reallocated! NOT SUCCESSFUL.",file_output);
			return;
		}
		memset(inputms[0][0].freq,'\0',sizeof(double));
		if((inputms[0][0].factor_pop = (double *)realloc(inputms[0][0].factor_pop,1*sizeof(double))) == 0) {
			puts(" Value not reallocated! NOT SUCCESSFUL.");
			if(file_output) fputs(" Value not reallocated! NOT SUCCESSFUL.",file_output);
			return;
		}
		memset(inputms[0][0].factor_pop,'\0',sizeof(double));
		
		inputms[0][0].nrec[0] = (double)1;
		inputms[0][0].npast[0] = (double)1;
		inputms[0][0].tpast[0] = (double)1;
		inputms[0][0].freq[0] = (double)1;
		inputms[0][0].factor_pop[0] = (double)1;
		
		for(i=1;i<global[0][0].nlocimat;i++) {
			free(inputms[0][i].nrec);
			free(inputms[0][i].npast);
			free(inputms[0][i].tpast);
			free(inputms[0][i].factor_pop);
			free(inputms[0][i].freq);
		}
		
		if((*inputms = (struct var2b *)realloc(*inputms,1*sizeof(struct var2b))) == 0) {
			puts(" Matrix not reallocated! NOT SUCCESSFUL.");
			if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
			return;
		}
		memset(*inputms,'\0',sizeof(struct var2b));
		*/
		if(global[0][0].montecarlo_sim) {
			if((*matrixmlsim = (struct statistisimmuloc *)realloc(*matrixmlsim,1*sizeof(struct statistisimmuloc))) == 0) {
				puts(" Matrix not reallocated! NOT SUCCESSFUL.");
				if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
				return;
			}
			memset(*matrixmlsim,'\0',sizeof(struct statistisimmuloc));
			if((*avgstatloci = (struct horizontalstatsml *)realloc(*avgstatloci,1*sizeof(struct horizontalstatsml))) == 0) {
				puts(" Matrix not reallocated! NOT SUCCESSFUL.");
				if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
				return;
			}
			memset(*avgstatloci,'\0',sizeof(struct horizontalstatsml));
			if(global[0][0].onlymulo == 0) {
				for(i=0;i<global[0][0].nlocimat;i++) 
					free(matrixsim[0][i]);
				if((*matrixsim = (struct statistisim **)realloc(*matrixsim,1*sizeof(struct statistisim *))) == 0) {
					puts(" Matrix not reallocated! NOT SUCCESSFUL.");
					if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
					return;
				}
			}
		}
		/*
		if(global[0][0].dataforsimt == 1) {
			if(global[0][0].ml0done == 1) {
				if((*lrdist = (struct LRTdist *)realloc(*lrdist,1*sizeof(struct LRTdist))) == 0) {
					puts(" Matrix not reallocated! NOT SUCCESSFUL.");
					if(file_output) fputs(" Matrix not reallocated! NOT SUCCESSFUL.",file_output);
					return;
				}
				memset(*lrdist,'\0',sizeof(struct LRTdist));
			}
		}
		*/
	}
	global[0][0].montecarlo_sim = 0;
	/*global[0][0].nlocimat = 0;*/
	global[0][0].nitermat = 0;
	/*global[0][0].dataforsim = 0;*/
	global[0][0].onlymulo = 0;
	
	/*global[0][0].dataindata = 0;*/
	/*global[0][0].dataequalsim = 0;*/

	/*global[0][0].ml0done = 0;*/
	/*global[0][0].nlocimatt = 0;*/
	/*global[0][0].nitermatt = 0;*/
	/*global[0][0].dataforsimt = 0;*/
	/*global[0][0].mltheta = 0;*/
					
	printf("\n\nResults in the memory have been erased.\n\n");
	if(file_output) 
		fprintf(file_output,"\n\nResults in the memory have been erased.\n\n");

	return;
}

void open_outfile(/*struct statistics *matrix,struct statmulo *matrixml,int *observed_data,int *outgroup,int *n_loci,*/FILE **file_output)
{
    char k[1];
    char name_outfile[512];	/* be careful with the limits... 511 characters */
    time_t now;
    struct tm *date;
    char s[80];
    
    time(&now);
    date = localtime(&now);
    strftime(s,80,"%c",date);
    
    if(*file_output == 0) {
        *name_outfile = 0;
        printf("\n\n  Please indicate the output file: ");
        scanf(" %[^\n]",name_outfile);
        if(*name_outfile != 0) {
            if((*file_output = fopen (name_outfile,"w+")) == 0) {
                puts("\n\n  Error in input/output");
                return;
            }
        }
        fputs(MANVa,*file_output);
        if(*file_output) fprintf(*file_output,"MANVa output file:  %s \n\n",name_outfile);
        printf("\n\n  The output file has been assigned.");
       	return;
    }
    else {
        printf("\n\n  Do you want to replace the output file (y/n)? ");
        
        do *k = getchar();
        while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');
        
        if(*k == 'y' || *k == 'Y') {
            *name_outfile = 0;
            printf("\n\n  Please indicate the name of the new output file: ");
            scanf(" %[^\n]",name_outfile);
            if(*name_outfile != 0) {
                fclose(*file_output);
                if((*file_output = fopen (name_outfile,"w+")) == 0) {
                    puts("\n\n  Error in input/output");
                    return;
                }
            }
            fputs(MANVa,*file_output);
            if(*file_output) fprintf(*file_output,"MANVa (Multilocus Analysis of Nucleotide Variation) output file:  %s \n\n",s);
            printf("\n\n  A new output file has been assigned.");
            return;
        }
        else if(*k == 'n' || *k == 'N') {
            printf("\n\n  Do you want to close the output file (y/n)? ");
            do *k = getchar();
            while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');
            if(*k == 'y' || *k == 'Y') {
                fclose(*file_output);
                file_output = 0;
                printf("\n\n  Output file closed.");
            }
            return;
        }
    }
}

void close_outfile(/*struct statistics *matrix,struct statmulo *matrixml,int *observed_data,int *outgroup,int *n_loci,*/FILE **file_output)
{
    char k[1];
    
    if(*file_output != 0) {
        printf("\n\n  Do you want to close the output file (y/n)? ");
        do *k = getchar();
        while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');
        if(*k == 'y' || *k == 'Y') {
            fclose(*file_output);
            *file_output = 0;
            printf("\n\n  Output file closed.");
        }
    }
    return;
}

int open_files1(/*struct statistics **matrix,struct statmulo **matrixml,int *observed_data,*/int *outgroup/*,int *n_loci*/,int *excludelines,FILE *file_output, char **name_outgroup, char **name_excluded/*, double *factor_chrn*/,int *includelines, char **name_ingroups)
{
    char k[1];
    char stemp[512];
    
    /* alert message in case output file is not yet defined */
	/*
    if(file_output == 0) {
        printf("\n\n  You do not have an output file defined. It is recommended to send all results to an output file. ");
        printf("\n  Do you want to return to the File Menu and define an output file (y/n)? ");
        
        do *k = getchar();
        while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');        
        if(*k == 'y' || *k == 'Y') return 0;
        else if(*k == 'n' || *k == 'N');
            else return 1;
    }
	*/
    /* Open files from the current folder. If n_loci > 0, then add to the matrix */

    printf("\n\n  Open files:");
	
	if(*outgroup == -1) {
        printf("\n\n  Are you including an outgroup (y/n)? ");
        
        do *k = getchar();
        while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');
        
        if(*k == 'y' || *k == 'Y') {
            printf("\n\n  Please indicate the word (or part of the word) to identify the OUTGROUP lines\n  in ALL your files (max. 50 characters): ");
            do{
                stemp[0] = 0;
                scanf(" %[^\n]",stemp);
            } while(stemp[0] == 0);
            memcpy(name_outgroup[0],stemp,50);
            *outgroup = 1;
            printf("\n  Lines with the characters '%s' will be outgroup.",name_outgroup[0]);
            if(file_output) fprintf(file_output,"\n  Lines with the characters '%s' will be outgroup.",name_outgroup[0]);
        }
        else if(*k == 'n' || *k == 'N') *outgroup = 0;
    }
	else {
		printf("\n  Lines with the characters '%s' will be outgroup.",name_outgroup[0]);
		if(file_output) fprintf(file_output,"\n  Lines with the characters '%s' will be outgroup.",name_outgroup[0]);
	}
    if(*excludelines == -1 && *includelines == -1) {
        printf("\n\n  Ingroup lines. Would you like to:\n    1 - Exclude specific lines \n    2 - Include specific lines \n    3 - Include all\n\n    ");
        
        do *k = getchar();
        while(*k<'1' || *k>'3');
        
        if(*k == '1') {
            printf("\n\n  Please indicate the word (or part of the word) to identify the EXCLUDED lines\n  in ALL your files (max. 50 characters): ");
            do{
                stemp[0] = 0;
                scanf(" %[^\n]",stemp);
            } while(stemp[0] == 0);
            memcpy(name_excluded[0],stemp,50);
            if(*outgroup == 1) {
                if((strstr(name_outgroup[0],name_excluded[0]) != 0) || (strstr(name_excluded[0],name_outgroup[0]) != 0)) {
                    printf("\n  otgroup name and excluded name are coincidents, sorry.");
                    if(file_output) fputs("\n  otgroup name and excluded name are coincidents, sorry.",file_output);
                    *outgroup = -1;
                    *excludelines = -1;
                    *includelines = -1;
                    return 0;
                }
            }
            *excludelines = 1;
            printf("\n  Lines with the characters '%s' will be excluded.\n",name_excluded[0]);
            if(file_output) 
                fprintf(file_output,"\n  Lines with the characters '%s' will be excluded.\n",name_excluded[0]);
			*includelines = -1;
        }
        if(*k == '2') {
            printf("\n\n  Please indicate the word (or part of the word) to identify the INGROUP lines\n  in ALL your files (max. 50 characters): ");
            do{
                stemp[0] = 0;
                scanf(" %[^\n]",stemp);
            } while(stemp[0] == 0);
            memcpy(name_ingroups[0],stemp,50);
            if(*outgroup == 1) {
                if((strstr(name_outgroup[0],name_ingroups[0]) != 0) || (strstr(name_ingroups[0],name_outgroup[0]) != 0)) {
                    printf("\n  otgroup name and ingroups name are coincidents, sorry.");
                    if(file_output) fputs("\n  otgroup name and ingroups name are coincidents, sorry.",file_output);
                    *outgroup = -1;
                    *excludelines = -1;
                    *includelines = -1;
                    return 0;
                }
            }
            *includelines = 1;
            printf("\n  Lines with the characters '%s' will be the ingroup.\n",name_ingroups[0]);
            if(file_output) 
                fprintf(file_output,"\n  Lines with the characters '%s' will be the ingroup.\n",name_ingroups[0]);
			*excludelines = -1;
        }
    }
	else {
		if(*excludelines != -1) {
			printf("\n  Lines with the characters '%s' will be excluded.\n",name_excluded[0]);
			if(file_output) 
				fprintf(file_output,"\n  Lines with the characters '%s' will be excluded.\n",name_excluded[0]);
		}
		else {
			printf("\n  Lines with the characters '%s' will be ingroup.\n",name_ingroups[0]);
			if(file_output) 
				fprintf(file_output,"\n  Lines with the characters '%s' will be ingroup.\n",name_ingroups[0]);
		}
	}
    return 1;
}
    
void open_files2(struct statistics **matrix,struct statmulo **matrixml,int *observed_data,int *outgroup,int *n_loci/*,int *montecarlo_sim*/,int *excludelines,
	FILE *file_output, char **name_outgroup, char **name_excluded,char **name_fileinputfolder,char **name_fileinput,int once,int *dataforsim,int lian_in,struct lianinput **liandata,double factor_chrn,
	char **name_fileinputfoldergff,char **name_fileinputgff,char **name_fileinputfolderGFF,char **name_fileinputGFF,int gfffiles,char *subset_positions/*,int ifgencode,char *codename*/,
	char *genetic_code,int *includelines,char **name_ingroups)
{    
    FILE *file_input;
    int get_obsdata(struct statistics **,struct statmulo ** /*,int * */,int *,int *,FILE *,FILE *,char *,int,char *,int,struct lianinput **,double,char **,char **,char ** /*,char ** */,int,
		char * /*,int,char * */,char *,int,char *);
    char *sl,*vect;
	int lenc,x;
	#if DOS_CL
	int n;
	#endif

	if(*n_loci == 32767) {
		printf("\n  Sorry, no more loci are allowed (maximum is 32.767 loci).");
		return;
	}
	if(once == 0) {
		printf("\n  Name of the data file (NBRF/pir or FASTA format)? \n");
		scanf(" %[^\n]",*name_fileinputfolder);
		lenc = (int)strlen(*name_fileinputfolder);
		vect = *name_fileinputfolder;
		#if UNIX_CL
		while((sl = memchr(vect,'/',lenc)) != 0) {
			vect = sl+1;
			lenc = (int)strlen(vect+1);
		}
		#elif DOS_CL
		if(name_fileinputfolder[0][0] == '/' || name_fileinputfolder[0][0] == '\\') {
			x=0;
			while(name_fileinputfolder[0][x] != '\0' && x < 254) x++;
			name_fileinputfolder[0][x+1] = '\0';
			for(n=x;n>0;n--) {
				name_fileinputfolder[0][n] = name_fileinputfolder[0][n-1];					
			}
			name_fileinputfolder[0][0] = '.';
		}
		for(n=0;n<256;n++) {
			if(name_fileinputfolder[0][n] == '/') {
				name_fileinputfolder[0][n] = '\\';
			}
			if(name_fileinputfolder[0][n] == '\0') break;
		}
		for(n=2;n<256;n++) {
			if(name_fileinputfolder[0][n] == '\0') {
				if(name_fileinputfolder[0][n-1] == '\\' && name_fileinputfolder[0][n-2] == '\\') 
					name_fileinputfolder[0][n-1] = '\0';
				break;
			}
		}
		while((sl = memchr(vect,'\\',lenc)) != 0) {
			vect = sl+1;
			lenc = (int)strlen(vect+1);
		}
		#endif
		for(x=0;x<=lenc;x++) 
			name_fileinput[0][x] = vect[x];
		name_fileinput[0][lenc+1] = '\0';
		
		if(gfffiles == 1) {
			printf("\n  Path and name of the .gff file (GFF format)? \n");
			scanf(" %[^\n]",*name_fileinputfoldergff);
			lenc = (int)strlen(*name_fileinputfoldergff);
			vect = *name_fileinputfoldergff;
			#if UNIX_CL
			while((sl = memchr(vect,'/',lenc)) != 0) {
				vect = sl+1;
				lenc = (int)strlen(vect+1);
			}
			#elif DOS_CL
			if(name_fileinputfoldergff[0][0] == '/' || name_fileinputfoldergff[0][0] == '\\') {
				x=0;
				while(name_fileinputfoldergff[0][x] != '\0' && x < 254) x++;
				name_fileinputfolder[0][x+1] = '\0';
				for(n=x;n>0;n--) {
					name_fileinputfoldergff[0][n] = name_fileinputfoldergff[0][n-1];					
				}
				name_fileinputfoldergff[0][0] = '.';
			}
			for(n=0;n<256;n++) {
				if(name_fileinputfoldergff[0][n] == '/') {
					name_fileinputfoldergff[0][n] = '\\';
				}
				if(name_fileinputfoldergff[0][n] == '\0') break;
			}
			for(n=2;n<256;n++) {
				if(name_fileinputfoldergff[0][n] == '\0') {
					if(name_fileinputfoldergff[0][n-1] == '\\' && name_fileinputfoldergff[0][n-2] == '\\') 
						name_fileinputfoldergff[0][n-1] = '\0';
					break;
				}
			}
			while((sl = memchr(vect,'\\',lenc)) != 0) {
				vect = sl+1;
				lenc = (int)strlen(vect+1);
			}
			#endif
			name_fileinputfolderGFF[0][0] = 0;
			strncat(*name_fileinputfolderGFF,*name_fileinputfoldergff,256);
			for(x=0;x<=lenc;x++) {
				name_fileinputgff[0][x] = vect[x];
				name_fileinputGFF[0][x] = vect[x];
			}
			name_fileinputgff[0][lenc+1] = '\0';
			name_fileinputGFF[0][lenc+1] = '\0';
		}
	}
	
	if((file_input = fopen (*name_fileinputfolder,"r")) == 0) {
		/*printf("\n  It is not possible to open the file %s.",*name_fileinputfolder);*/
	}
	if(file_input) {		
		printf("\n Opening file %s:",*name_fileinputfolder);
		if(file_output) fprintf(file_output,"\n Opening file %s:",*name_fileinputfolder);
		if(get_obsdata(matrix,matrixml/*,observed_data*/,outgroup,n_loci,file_output,file_input,name_outgroup[0],
				*excludelines,name_excluded[0],lian_in,liandata,factor_chrn,
				name_fileinputfoldergff,name_fileinputgff,name_fileinputfolderGFF/*,name_fileinputGFF*/,
				gfffiles,subset_positions/*,ifgencode,codename*/,genetic_code,*includelines,name_ingroups[0])) {
			memcpy(matrix[0][*n_loci].gene,*name_fileinput,20);
			matrix[0][*n_loci].gene[20] = '\0';
			matrix[0][*n_loci].theta_mlS = (double)-10000;
			*observed_data = 1;
			/**montecarlo_sim = 0;*/
			*dataforsim = 0;
			*n_loci += 1;
			printf(" SUCCESSFUL.");
			if(file_output) fputs(" SUCCESSFUL.",file_output);
		}
		else {
			printf(" NOT SUCCESSFUL.");
			if(file_output) fputs(" NOT SUCCESSFUL.",file_output);
		}
		fclose(file_input);
	}
    
    return;
}
void gff_function(int *gfffiles,int *observed_data,char *subset_positions,int *ifgencode,char *code_name,char *genetic_code,FILE **file_output)
{
	char k[1];
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
	if(*gfffiles == 0) {
		/*ask if gff files included*/
		printf("\n\n  Do you want to use a subset of positions from data (y/n)? ");
		printf("\n  (Warning: this option needs a GFF-format file associated to each alignment) ");		
		do *k = getchar();
		while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');		
		if(*k == 'y' || *k == 'Y') {
			*gfffiles = 1;
		}
		else if(*k == 'n' || *k == 'N') *gfffiles = -1;
	
		if(*gfffiles == 1) {
			if(*observed_data == 0) {
				/*if yes, ask for kind of subset, menu*/
				printf("\n\n  Please chose the subset of positions:\n\n ");
				printf("\n  0 - Coding. ");
				printf("\n  1 - NonCoding. ");
				printf("\n  2 - Synonymous. ");
				printf("\n  3 - NonSynonymous. ");
				printf("\n  4 - Silent (Synonymous + NonCoding). ");
				printf("\n  5 - Others.");
				printf("\n  6 - NO SUBSET (All positions).");
				printf("\n\n  OPTION: ");
				do *k = getc(stdin);
				while(*k<'0' || *k>'6');		
				switch(*k) {
					case '0':
						strncpy(subset_positions,"coding\0",20);
						break;
					case '1':
						strncpy(subset_positions,"noncoding\0",20);
						break;
					case '2':
						strncpy(subset_positions,"synonymous\0",20);
						break;
					case '3':
						strncpy(subset_positions,"nonsynonymous\0",20);
						break;
					case '4':
						strncpy(subset_positions,"silent\0",20);
						break;
					case '5':
						printf("\n\n  Please specify the name of the field to be analyzed\n (less than 20 characters):\n");
						scanf(" %[^\n]",subset_positions);
						break;
					case '6':
						*gfffiles = -1;
						return;
						break;
				}
				if(*file_output) fprintf(*file_output,"\n\n  Subset of positions to be analyzed: %s\n\n",subset_positions);
			}
		}
	}
	if(*gfffiles == 1) {
		if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0) {
			/*ask for genetic code, menu*/
			*ifgencode = 1;
			printf("\n\n  Please chose the genetic code:\n\n ");
			printf("\n  0 - Nuclear Universal. ");
			printf("\n  1 - mtDNA Drosophila. ");
			printf("\n  2 - mtDNA Mammals. ");
			printf("\n  3 - Other. ");
			printf("\n\n  OPTION: ");
			
			do *k = getc(stdin);
			while(*k<'0' || *k>'3');		
			switch(*k) {
				case '0':
					genetic_code[0] = 'F';
					genetic_code[1] = 'F';
					genetic_code[2] = 'L';
					genetic_code[3] = 'L';
					genetic_code[4] = 'S';
					genetic_code[5] = 'S';
					genetic_code[6] = 'S';
					genetic_code[7] = 'S';
					genetic_code[8] = 'Y';
					genetic_code[9] = 'Y';
					genetic_code[10] = '*';
					genetic_code[11] = '*';
					genetic_code[12] = 'C';
					genetic_code[13] = 'C';
					genetic_code[14] = '*';
					genetic_code[15] = 'W';
					genetic_code[16] = 'L';
					genetic_code[17] = 'L';
					genetic_code[18] = 'L';
					genetic_code[19] = 'L';
					genetic_code[20] = 'P';
					genetic_code[21] = 'P';
					genetic_code[22] = 'P';
					genetic_code[23] = 'P';
					genetic_code[24] = 'H';
					genetic_code[25] = 'H';
					genetic_code[26] = 'Q';
					genetic_code[27] = 'Q';
					genetic_code[28] = 'R';
					genetic_code[29] = 'R';
					genetic_code[30] = 'R';
					genetic_code[31] = 'R';
					genetic_code[32] = 'I';
					genetic_code[33] = 'I';
					genetic_code[34] = 'I';
					genetic_code[35] = 'M';
					genetic_code[36] = 'T';
					genetic_code[37] = 'T';
					genetic_code[38] = 'T';
					genetic_code[39] = 'T';
					genetic_code[40] = 'N';
					genetic_code[41] = 'N';
					genetic_code[42] = 'K';
					genetic_code[43] = 'K';
					genetic_code[44] = 'S';
					genetic_code[45] = 'S';
					genetic_code[46] = 'R';
					genetic_code[47] = 'R';
					genetic_code[48] = 'V';
					genetic_code[49] = 'V';
					genetic_code[50] = 'V';
					genetic_code[51] = 'V';
					genetic_code[52] = 'A';
					genetic_code[53] = 'A';
					genetic_code[54] = 'A';
					genetic_code[55] = 'A';
					genetic_code[56] = 'D';
					genetic_code[57] = 'D';
					genetic_code[58] = 'E';
					genetic_code[59] = 'E';
					genetic_code[60] = 'G';
					genetic_code[61] = 'G';
					genetic_code[62] = 'G';
					genetic_code[63] = 'G';
					
					strncpy(code_name,"Nuclear Universal\0",50);
					printf("\n  Nuclear Universal code:\n\n");
					if(*file_output) fprintf(*file_output,"\n  Nuclear Universal code:\n\n");
					break;
				case '1':
					genetic_code[0] = 'F';
					genetic_code[1] = 'F';
					genetic_code[2] = 'L';
					genetic_code[3] = 'L';
					genetic_code[4] = 'S';
					genetic_code[5] = 'S';
					genetic_code[6] = 'S';
					genetic_code[7] = 'S';
					genetic_code[8] = 'Y';
					genetic_code[9] = 'Y';
					genetic_code[10] = '*';
					genetic_code[11] = '*';
					genetic_code[12] = 'C';
					genetic_code[13] = 'C';
					genetic_code[14] = 'W';
					genetic_code[15] = 'W';
					genetic_code[16] = 'L';
					genetic_code[17] = 'L';
					genetic_code[18] = 'L';
					genetic_code[19] = 'L';
					genetic_code[20] = 'P';
					genetic_code[21] = 'P';
					genetic_code[22] = 'P';
					genetic_code[23] = 'P';
					genetic_code[24] = 'H';
					genetic_code[25] = 'H';
					genetic_code[26] = 'Q';
					genetic_code[27] = 'Q';
					genetic_code[28] = 'R';
					genetic_code[29] = 'R';
					genetic_code[30] = 'R';
					genetic_code[31] = 'R';
					genetic_code[32] = 'I';
					genetic_code[33] = 'I';
					genetic_code[34] = 'M';
					genetic_code[35] = 'M';
					genetic_code[36] = 'T';
					genetic_code[37] = 'T';
					genetic_code[38] = 'T';
					genetic_code[39] = 'T';
					genetic_code[40] = 'N';
					genetic_code[41] = 'N';
					genetic_code[42] = 'K';
					genetic_code[43] = 'K';
					genetic_code[44] = 'S';
					genetic_code[45] = 'S';
					genetic_code[46] = 'S';
					genetic_code[47] = 'S';
					genetic_code[48] = 'V';
					genetic_code[49] = 'V';
					genetic_code[50] = 'V';
					genetic_code[51] = 'V';
					genetic_code[52] = 'A';
					genetic_code[53] = 'A';
					genetic_code[54] = 'A';
					genetic_code[55] = 'A';
					genetic_code[56] = 'D';
					genetic_code[57] = 'D';
					genetic_code[58] = 'E';
					genetic_code[59] = 'E';
					genetic_code[60] = 'G';
					genetic_code[61] = 'G';
					genetic_code[62] = 'G';
					genetic_code[63] = 'G';

					strncpy(code_name,"mtDNA Drosophila\0",50);
					printf("\n  mtDNA Drosophila code:\n\n");
					if(*file_output) fprintf(*file_output,"\n  mtDNA Drosophila code:\n\n");
					break;
				case '2':
					genetic_code[0] = 'F';
					genetic_code[1] = 'F';
					genetic_code[2] = 'L';
					genetic_code[3] = 'L';
					genetic_code[4] = 'S';
					genetic_code[5] = 'S';
					genetic_code[6] = 'S';
					genetic_code[7] = 'S';
					genetic_code[8] = 'Y';
					genetic_code[9] = 'Y';
					genetic_code[10] = '*';
					genetic_code[11] = '*';
					genetic_code[12] = 'C';
					genetic_code[13] = 'C';
					genetic_code[14] = 'W';
					genetic_code[15] = 'W';
					genetic_code[16] = 'L';
					genetic_code[17] = 'L';
					genetic_code[18] = 'L';
					genetic_code[19] = 'L';
					genetic_code[20] = 'P';
					genetic_code[21] = 'P';
					genetic_code[22] = 'P';
					genetic_code[23] = 'P';
					genetic_code[24] = 'H';
					genetic_code[25] = 'H';
					genetic_code[26] = 'Q';
					genetic_code[27] = 'Q';
					genetic_code[28] = 'R';
					genetic_code[29] = 'R';
					genetic_code[30] = 'R';
					genetic_code[31] = 'R';
					genetic_code[32] = 'I';
					genetic_code[33] = 'I';
					genetic_code[34] = 'M';
					genetic_code[35] = 'M';
					genetic_code[36] = 'T';
					genetic_code[37] = 'T';
					genetic_code[38] = 'T';
					genetic_code[39] = 'T';
					genetic_code[40] = 'N';
					genetic_code[41] = 'N';
					genetic_code[42] = 'K';
					genetic_code[43] = 'K';
					genetic_code[44] = 'S';
					genetic_code[45] = 'S';
					genetic_code[46] = '*';
					genetic_code[47] = '*';
					genetic_code[48] = 'V';
					genetic_code[49] = 'V';
					genetic_code[50] = 'V';
					genetic_code[51] = 'V';
					genetic_code[52] = 'A';
					genetic_code[53] = 'A';
					genetic_code[54] = 'A';
					genetic_code[55] = 'A';
					genetic_code[56] = 'D';
					genetic_code[57] = 'D';
					genetic_code[58] = 'E';
					genetic_code[59] = 'E';
					genetic_code[60] = 'G';
					genetic_code[61] = 'G';
					genetic_code[62] = 'G';
					genetic_code[63] = 'G';
					
					strncpy(code_name,"mtDNA Mammals\0",50);
					printf("\n  mtDNA Mammals code:\n\n");
					if(*file_output) fprintf(*file_output,"\n  mtDNA Mammals code:\n\n");
					break;
				case '3':
					printf("\n\n  Please specify the corresponding aminacid for each codon:\n  (one letter (UPPERCASE) aa code. Stop codon is '*')\n");
					for(n=0;n<64;n++) {
						printf("%c%c%c : ",tripletsU[n][0],tripletsU[n][1],tripletsU[n][2]);								
						do {
							*k = getchar();
							if(*k == '*') break;
						} while(*k<'A' || *k>'Z');
						genetic_code[n] = *k;
					}
					printf("\n\n  Name of the genetic code? ");
					scanf(" %19s",code_name);
					break;
			}
			for(n=0;n<64;n += 16) {
				for(m=n;m<n+4;m += 1) {
					printf("  %c%c%c %c   %c%c%c %c   %c%c%c %c   %c%c%c %c\n",
						tripletsU[m+0][0],tripletsU[m+0][1],tripletsU[m+0][2],genetic_code[m+0],
						tripletsU[m+4][0],tripletsU[m+4][1],tripletsU[m+4][2],genetic_code[m+4],
						tripletsU[m+8][0],tripletsU[m+8][1],tripletsU[m+8][2],genetic_code[m+8],
						tripletsU[m+12][0],tripletsU[m+12][1],tripletsU[m+12][2],genetic_code[m+12]);
					if(*file_output) {
						fprintf(*file_output,"  %c%c%c %c   %c%c%c %c   %c%c%c %c   %c%c%c %c\n",
							tripletsU[m+0][0],tripletsU[m+0][1],tripletsU[m+0][2],genetic_code[m+0],
							tripletsU[m+4][0],tripletsU[m+4][1],tripletsU[m+4][2],genetic_code[m+4],
							tripletsU[m+8][0],tripletsU[m+8][1],tripletsU[m+8][2],genetic_code[m+8],
							tripletsU[m+12][0],tripletsU[m+12][1],tripletsU[m+12][2],genetic_code[m+12]);
					}
				}
			}
		}
		else *ifgencode = 0;
	}			
	return;
}

