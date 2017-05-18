/*
 *  allocationsim.c
 *  MuLoNeTests
 *
 *  Created by sonsins on Tue Mar 18 2003.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int reallocvectors(int *onlymulo/*,int *maxloc*/,struct statistisim ***matrixsim,struct statistisimmuloc **matrixmlsim,struct horizontalstatsml **avgstatloci,long int nitermat, int nlocimat, FILE *file_output)
{
    int x,y;
	*onlymulo = 0;
    if(!(nlocimat == 0) && !(nitermat == 0)) {
        /*matrixsimml*/
        if((*matrixmlsim = (struct statistisimmuloc *) realloc(*matrixmlsim,(long int)(nitermat)*sizeof(struct statistisimmuloc))) == 0) {
            puts("\nError in allocation. realloc_simmatrix.3 \n");
            return 0;
        }
		memset(*matrixmlsim,'\0',(long int)(nitermat)*sizeof(struct statistisimmuloc));
        /*avgstatloci*/
        if((*avgstatloci = (struct horizontalstatsml *) realloc(*avgstatloci,(long int)(nlocimat)*sizeof(struct horizontalstatsml))) == 0) {
            puts("\nError in allocation. realloc_vectors.4 \n");
            return 0;
        }
		memset(*avgstatloci,'\0',(long int)(nlocimat)*sizeof(struct horizontalstatsml));
        /*matrixsim*/
        if((matrixsim[0] = (struct statistisim **)realloc(matrixsim[0],(long int)(nlocimat)*sizeof(struct statistisim *))) == 0) {
            puts("\n memory not reallocated for keeping all results. Only multilocus analysis available, sorry.");
            if(file_output) 
                fputs("\n memory not reallocated for keeping all results. Only multilocus analysis available, sorry. ",
                    file_output);
            *onlymulo = 1;
            return 0;
        }
        /*if(*maxloc < nlocimat) {*/
            /*
			for(x=0;x<*maxloc;x++) {
                if((matrixsim[0][x] = (struct statistisim *) realloc(matrixsim[0][x],nitermat*sizeof(struct statistisim))) == 0) {
                    puts("\n memory not reallocated for keeping all results. Only multilocus analysis available, sorry.");
                    if(file_output) 
                        fputs("\n memory not reallocated for keeping all results. Only multilocus analysis available, sorry. ",
                            file_output);
                    *onlymulo = 1;
                    return 0;
                }
            }
			*/
            for(x=0/**maxloc*/;x<nlocimat;x++) {
                if((matrixsim[0][x] = (struct statistisim *) calloc((unsigned)nitermat,sizeof(struct statistisim))) == 0) {
                    puts("\n memory not reallocated for keeping all results. Only multilocus analysis available, sorry.");
                    if(file_output) 
                        fputs("\n memory not reallocated for keeping all results. Only multilocus analysis available, sorry. ",
                            file_output);
                    *onlymulo = 1;
					/*free occupied memory*/
					for(y=0;y<x;y++) {
						free(matrixsim[0][y]);
					}
                    return 0;
                }
            }
        /*}*/
        /**maxloc = nlocimat;*/
    }
    else return 1;

    return 0;
}

int freeallocvectors(int *onlymulo/*,int *maxloc*/,struct statistisim ***matrixsim/*,long int nitermat*/, int nlocimat,int dataindata,int *dataforsim/*, FILE *file_output*/)
{
	int x;
	
	if(*dataforsim == 1) {
		if(dataindata == 1 || dataindata == 3) {
			if(*onlymulo==0) {
				for(x=0;x<nlocimat;x++) {
					free(matrixsim[0][x]);
				}
			}
		}
	}
	*dataforsim = 0;
	return 0;
}

int reallocvectorsML(/*int *maxloc,*/struct LRTdist **lrdist/*,struct MLthetasim ***mthetasim*/,long int nitermatt, int nlocimatt/*, FILE *file_output*/)
{
    /*int x;*/
    
    if(!(nlocimatt == 0) && !(nitermatt == 0)) {
        /*lrdist*/
        if((*lrdist = (struct LRTdist *)realloc(*lrdist,(long int)(nitermatt)*sizeof(struct LRTdist))) == 0) {
            puts("\nError in allocation. realloc_simmatrix.3 \n");
            return 0;
        }
		memset(*lrdist,'\0',(long int)(nitermatt)*sizeof(struct LRTdist));
        /*mthetasim*/
        /*
		if((*mthetasim = (struct MLthetasim **)realloc(*mthetasim,(long int)(nlocimatt)*sizeof(struct MLthetasim *))) == 0) {
            puts("\nError in allocation. p6\n");
            if(file_output) 
                fputs("\nError in allocation. p6\n",
                    file_output);
            return 1;
        }
        if(*maxloc < nlocimatt) {
            for(x=0;x<*maxloc;x++) {
                if((mthetasim[0][x] = (struct MLthetasim *)realloc(mthetasim[0][x],nitermatt*sizeof(struct MLthetasim))) == 0) {
                    puts("\nError in allocation. realloc_simmatrix.4 \n");
                    if(file_output) 
                        fputs("\nError in allocation. realloc_simmatrix.4 \n",
                            file_output);
                    return 1;
                }
            }
            for(x=*maxloc;x<nlocimatt;x++) {
                if((mthetasim[0][x] = (struct MLthetasim *)calloc(nitermatt,sizeof(struct MLthetasim))) == 0) {
                    puts("\nError in allocation. realloc_simmatrix.5 \n");
                    if(file_output) 
                        fputs("\nError in allocation. realloc_simmatrix.5 \n",
                            file_output);
                    return 1;
                }
            }
        }
		*/
        /**maxloc = nlocimatt;*/
    }
    else return 1;

	return 0;
}


