/*
 *  get_obsdata.c
 *  MuLoNeTests
 *
 *  Created by sebas on Sun Feb 23 2003.
 *
 */

#include "MuLoNeTests.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

int get_obsdata(struct statistics **matrix,struct statmulo **matrixml/*,int *observed_data*/,int *outgroup,int *n_loci,FILE *file_output,FILE *file_input, char *name_outgroup,
	int excludelines,char *name_excluded,int lian_in,struct lianinput **liandata,double factor_chrn,
	char **name_fileinputfoldergff,char **name_fileinputgff,char **name_fileinputfolderGFF/*,char **name_fileinputGFF*/,int gfffiles,char *subset_positions/*,int ifgencode,char *codename*/,
	char *genetic_code,int includelines,char *name_ingroups)
{    
    char *f;
    long int count,xx;
    int c;
    int n_sam;
    long int n_sit;
    int nseq;
    int maxsam;

    static char **names = 0; /* limit of the name of the lines to 50 characters. be careful ... */
    static char *DNA_matr = 0;
    int n_samp;
    long int n_site;
    
    int n_excl = 0;
    
    int x;
    int var_char(FILE * /*,char * */,long int *,int *,int *,long int *,int *,int *,char ***,char ** /*,int * */,
        long int *,int,char *,int *,int,char *,char *,int);    
    int get_obsstats(struct statistics **,struct statmulo ** /*,int * */,int *,int *,FILE * /*,FILE * */,char *,
        int, long int,char *,char **,int,int,struct lianinput **,double,double *,double *);

	static double *matrix_sizepos = 0; /*size of each position, 0 not include, >0 include, Syn/Nsyn positions are allowed*/
	static double *matrix_segrpos = 0; /*always 1 except in those biallelic positions that are not desired (e.g., choose syn but we do not want nsyn)*/
	static long int maxsites = 0;
	/*include **name_fileinputfoldergff,**name_fileinputgff,*subset_positions,ifgencode,*codename,*genetic_code,*matrix_sizepos,n_samp,n_site,*DNA_matr*/
	int use_gff(char **,char **,char **,/*char **, */char *, /*int,char *, */char *,double *,int,long int,char *,double *,FILE *);
	
    count = 0;
    c = 0;
    n_sam = 0;
    n_sit = 0;
    nseq  = 0;
    maxsam= 128;
    n_samp= 0;
    n_site= 0;

    if(names == 0) { /* only initialize once. Check */
        if ((names = (char **)calloc(128,sizeof(char *))) == 0) {
            puts("\nError: memory not reallocated. get_obsdata.1 \n");
            return(0);
        }
        for(x=0;x<128;x++) {
            if ((names[x] = (char *)calloc(50,sizeof(char))) == 0) {
                puts("\nError: memory not reallocated. get_obsdata.2 \n");
                return(0);
            }
        }    
        if ((DNA_matr = (char *)calloc(10000,sizeof(char))) == 0) {
            puts("\nError: memory not reallocated. get_obsdata.3 \n");
            return(0);
        }
    }
    
    if(!(f = (char *)malloc(BUFSIZ))) {
        puts("\nError: memory not reallocated. get_obsdata.4 \n");
        return(0);
    }
    setbuf(file_input,f);

    c = fgetc(file_input);
    while (c != 0 && c != -1) {
        while(c == 32 || c == 9 || c == 13 || c == 10 || c == '*') c = fgetc(file_input);
        n_sit = 0;
        if(!(var_char(file_input/*,f*/,&count,&c,&n_sam,&n_sit,&nseq,&maxsam,&names,&DNA_matr/*,&n_samp*/,
            &n_site,excludelines,name_excluded,&n_excl,includelines,name_ingroups,name_outgroup,*outgroup)))
            return 0;
        /*if(c !='>') c = fgetc(file_input);*//*???*/
    }
    n_samp = n_sam;
    free(f);
    
    if(n_sam == 0 || n_site == 0) return(0);
    else {
        if(n_sam > 32167) {
            puts("Error: too much samples. Only 32167 samples per loci are allowed.\n");
            return(0);
        }
		/*init matrix_sizepos*/
		if(matrix_sizepos == 0) {
			if((matrix_sizepos = (double *)malloc(n_site*sizeof(double))) == 0) {
				puts("Error: memory not reallocated. get_obsstat.2"); 
				return(0);
			}
			if((matrix_segrpos = (double *)malloc(n_site*sizeof(double))) == 0) {
				puts("Error: memory not reallocated. get_obsstat.2"); 
				return(0);
			}
			maxsites = n_site;
		}
		else{
			if(n_site > maxsites) {
				if((matrix_sizepos = (double *)realloc(matrix_sizepos,n_site*sizeof(double))) == 0) {
					puts("Error: memory not reallocated. get_obsstat.2b"); 
					return(0);
				}
				if((matrix_segrpos = (double *)realloc(matrix_segrpos,n_site*sizeof(double))) == 0) {
					puts("Error: memory not reallocated. get_obsstat.2b"); 
					return(0);
				}
				maxsites = n_site;
			}
		}
		for(xx=0;xx<n_site;xx++) {
			matrix_sizepos[xx] = (double)1;
			matrix_segrpos[xx] = (double)1;
		}
		/*here include a function to filter positions: to read gff files (if necessary)*/
		if(gfffiles == 1) {
			/*include **name_fileinputfoldergff,**name_fileinputgff,*subset_positions,ifgencode,*codename,*genetic_code,*matrix_sizepos,n_samp,n_site,*DNA_matr*/
			/*the function read the gff file and cut the DNA_matr, also gives the number of positions in matrix_sizepos and count the total in n_site*/
			/*only modify values in matrix_sizepos*/
			if(use_gff(name_fileinputfoldergff,name_fileinputgff,name_fileinputfolderGFF,/*name_fileinputGFF,*/subset_positions,/*ifgencode,codename,*/
				genetic_code,matrix_sizepos,n_samp,n_site,DNA_matr,matrix_segrpos,file_output) == 0) {
				/*if error realloc DNA_matr*/
				if((DNA_matr = realloc(DNA_matr,((long int)1*sizeof(char)))) == 0) {
					puts("Error: realloc error varchar.01\n");
					return(0);
				}    
				return(0);
			}
		}
        /*function to analyze all data*/
        if(get_obsstats(matrix,matrixml/*,observed_data*/,outgroup,n_loci,file_output/*,file_input*/,name_outgroup,
            n_samp,n_site,DNA_matr,names,n_excl,lian_in,liandata,factor_chrn,matrix_sizepos,matrix_segrpos)) {
            if((DNA_matr = realloc(DNA_matr,((long int)1*sizeof(char)))) == 0) {
                puts("Error: realloc error varchar.1\n");
                return(0);
            }    
            return(1);
        }
        else
			return(0);
    }
    return(1);
}

int var_char(FILE *file_input/*,char *f*/,long int *count,int *c,int *n_sam,long int *n_sit,int *nseq,int *maxsam,char ***names,char **DNA_matr/*,int *n_samp*/,
	long int *n_site,int excludelines,char *name_excluded,int *n_excl,int includelines,char *name_ingroups,char *name_outgroup,int outgroup)
{
    int  aa = 0;
    int  bb = 0;
    long int  dd;
    double ee;
    char *strin;
    long int t;
    
    int assigna(FILE * /*,char * */,int *,int * /*,int * */,int *,char ***);

    aa = assigna(file_input/*,f*/,c,nseq/*,n_sam*/,maxsam,names);
	if(aa == 1) {
		if(outgroup > 0) {
			if(((strin = strstr(names[0][*nseq-1],name_outgroup)) == 0)) {
				if(excludelines > 0) {
					if(((strin = strstr(names[0][*nseq-1],name_excluded)) != 0)) {
						*nseq -= 1;
						*n_excl += 1;
						*c = fgetc(file_input);
						while(*c != '*' && *c != -1 && *c != 0 && *c != '>')
							*c = fgetc(file_input);
						aa = 0;
					}
				}
				if(includelines > 0) {
					if(((strin = strstr(names[0][*nseq-1],name_ingroups)) == 0)) {
						*nseq -= 1;
						*n_excl += 1;
						*c = fgetc(file_input);
						while(*c != '*' && *c != -1 && *c != 0 && *c != '>')
							*c = fgetc(file_input);
						aa = 0;
					}
				}
			}
		}
		else {
			if(excludelines > 0) {
				if(((strin = strstr(names[0][*nseq-1],name_excluded)) != 0)) {
					*nseq -= 1;
					*n_excl += 1;
					*c = fgetc(file_input);
					while(*c != '*' && *c != -1 && *c != 0 && *c != '>')
						*c = fgetc(file_input);
					aa = 0;
				}
			}
			if(includelines > 0) {
				if(((strin = strstr(names[0][*nseq-1],name_ingroups)) == 0)) {
					*nseq -= 1;
					*n_excl += 1;
					*c = fgetc(file_input);
					while(*c != '*' && *c != -1 && *c != 0 && *c != '>')
						*c = fgetc(file_input);
					aa = 0;
				}
			}
		}
	}
    if(aa == 1) {
        while(bb == 0) {
            dd = (long int)floor((double)*count/(double)10000);
            ee = (double)*count/(double)10000;
            
            if(dd == ee) {
                if((t=(((long int)dd+(long int)1)*(long int)10000)) > 2147483647) {
                    puts("Error: too much positions.\n");
                    return(0);
                }
                if((*DNA_matr = realloc(*DNA_matr,((long int)dd+(long int)1)*(long int)10000*sizeof(char))) == 0) {
                    puts("Error: realloc error varchar.1\n");
                    return(0);
                }    
            }
            switch(*c/* = fgetc(file_input)*/) {
                case 'T':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '1';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 't':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '1';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'U':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '1';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'u':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '1';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'C':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '2';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'c':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '2';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'G':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '3';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'g':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '3';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'A':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '4';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'a':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '4';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 0:
                    bb = 1; /*in FASTA case*/
                    break;
                case -1:
                    bb = 1; /*in FASTA case*/
                    break;
                case 10:
                    break;
                case 13:
                    break;
                case 32:
                    break;
                case '*': /* in NBRF case*/
                    bb = 1;
                    break;
                case '>': /* in FASTA case*/
                    bb = 1;
                    break;
                case '-':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '5';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case '?':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'N':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'n':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                /*degenerated are converted to N*/
                case 'W':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'w':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'M':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'm':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'R':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'r':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'Y':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'y':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'K':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'k':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'S':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 's':
                    DNA_matr[0][(((unsigned long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                default:
                    puts("Unexpected value in file");
                    printf("%d",*c);
                    return(0);
                    break;
            }
            if(*c !='>'/*|| *c != 0 || *c != -1*/) *c = fgetc(file_input);
        }
        *n_sam += 1;
        if(*n_site == 0) *n_site = *n_sit;
        else if(*n_site != *n_sit) {
            puts("The number of sites are not equal in all lines in the alignment.");
            return(0);
        }
    }
    return 1;
}

int assigna(FILE *file_input/*,char *f */,int *c,int *nseq /*,int *n_sam */,int *maxsam,char ***names)
{
    int N_VAR = 2;
    char var_file[2][50]  =
    {
        { ">"},
        { ">DL;"},
    };

    int i_;
    int j;
    int nn;
    int x;
    int c0; 
    
    j = 0;
    for(i_=0;i_<N_VAR;i_++) {
        while((var_file[i_][j]) == *c && (var_file[i_][j]) != '\0' && c != '\0') {
            *c = fgetc(file_input);
            j++;
        }
    }
    if(j<4 && j>0) i_= 1;/*FASTA*/
    else if(j==4) i_= 2;/*NBRF*/
        else {
            i_=0;/*NO ACCEPTED*/
            if(*c != 0 && *c != -1) *c = fgetc(file_input);
        }
        
    if((i_ == 2  && *c != 0 && *c != -1)) { /* NBRF files */
        nn = 0;
        while(*c != '\0' && *c != 13 && *c != 10 && nn < 50-2) {
            if(*c != '\t' && *c != 32) {
				names[0][*nseq][nn] = (char)*c;
				nn++;
			}
			*c = fgetc(file_input);
        }
        names[0][*nseq][nn] = '\0';
        *nseq += 1;
        if(*nseq == *maxsam) {
            *maxsam += 32;
            if(*maxsam > 32767) {
                puts("\n Sorry, no more samples are allowed.");
                    return 0;
            }
            if ((*names = (char **)realloc(*names,*maxsam*sizeof(char *))) == 0) {
                puts("\nError: memory not reallocated. assigna.1 \n");
                return(0);
            }
            for(x=*nseq;x<*maxsam;x++) {
                if ((names[0][x] = (char *)calloc(50,sizeof(char))) == 0) {
                    puts("\nError: memory not reallocated. assigna.2 \n");
                    return(0);
                }
            }
        }
        /*use unix or macos or dos format. begin*/
        while(*c != '\0' && *c != 13 && *c != 10 && *c != -1 && *c != 0)
            *c = fgetc(file_input);

        c0 = *c;
        *c = fgetc(file_input);
        if(c0 == 13 && *c == 10) *c = fgetc(file_input);
        while(*c != 10 && *c != 13 && *c != -1 && c != 0) 
            *c = fgetc(file_input);
        if(*c == -1 || *c == 0) {
            puts("\n Unexpected end of file");
            return(0);
        }
        c0 = *c;
        *c = fgetc(file_input);
        if(c0 == 13 && *c == 10) *c = fgetc(file_input);
        if(*c == -1 || *c == 0) {
            puts("\n Unexpected end of file");
            return(0);
        }
        /*use unix or macos or dos format. end*/
        
        return(1);
    }
    else {	
        if(i_ == 1 && *c != 0 && *c != -1) {	/* FASTA files */
            nn = 0;
            while(*c != '\0' && *c != 13 && *c != 10 && nn < 50-2) {
				if(*c != '\t' && *c != 32) {
					names[0][*nseq][nn] = (char)*c;
					nn++;
				}
				*c = fgetc(file_input);
            }
            names[0][*nseq][nn] = '\0';
            *nseq += 1;
            if(*nseq == *maxsam) {
                *maxsam += 32;
                if(*maxsam > 32767) {
                    puts("\n Sorry, no more samples are allowed.");
                        return 0;
                }
                if ((*names = (char **)realloc(*names,*maxsam*sizeof(char *))) == 0) {
                    puts("\nError: memory not reallocated. assigna.1 \n");
                    return(0);
                }
                for(x=*nseq;x<*maxsam;x++) {
                    if ((names[0][x] = (char *)calloc(50,sizeof(char))) == 0) {
                        puts("\nError: memory not reallocated. assigna.2 \n");
                        return(0);
                    }
                }
            }
            /*use unix or macos or dos format. begin*/
            while(*c != '\0' && *c != 13 && *c != 10 && *c != -1 && *c != 0)
                *c = fgetc(file_input);
    
            c0 = *c;
            *c = fgetc(file_input);
            if(c0 == 13 && *c == 10) *c = fgetc(file_input);
            if(*c == -1 || *c == 0) {
                puts("\n Unexpected end of file");
                return 0;
            }
            /*use unix or macos or dos format. end*/
            return(1);
        }
        else return(0);
    }
    return 0;
}
