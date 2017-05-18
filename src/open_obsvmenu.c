/*
 *  open_obsvmenu.c
 *  MuLoNeTests
 *
 *  Created by sebas on Mon Feb 24 2003.
 *
 */

#include "MuLoNeTests.h"
#include <stdio.h>
#include <stdlib.h>

void open_obsvmenu(struct statistics *matrix,struct statmulo *matrixml,int *observed_data,int *outgroup,int *n_loci,char *subset_positions,FILE *file_output)
{
    char k[1];
    void open_displstmenu(struct statistics *,struct statmulo * /*,int * */,int *,int *,char *,FILE *);
    void open_displntmenu(struct statistics *,struct statmulo * /*,int * */,int *,int *,char *,FILE *);
    
    if(*observed_data) {    
        while(1) {
            #if COMMAND_LINE
            
            printf("\n\nMENU 1. Observed analysis menu:\n\n");
            
            printf("CHOOSE:\n");
            printf(" 0 - Display summary statistics.\n");
            printf(" 1 - Display neutrality tests.\n");
            printf(" 2 - Back to Main menu.\n");
			/*
            if(file_output) {
                
                fprintf(file_output,"\n\n     MENU:\n     1. Observed analysis menu:\n\n");
                
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
                    open_displstmenu(matrix,matrixml/*,observed_data */,outgroup,n_loci,subset_positions,file_output);
                    break;
                case '1':
                    open_displntmenu(matrix,matrixml/*,observed_data */,outgroup,n_loci,subset_positions,file_output);
                    break;
                case '2':
                    return;
                    break;
            }
        }
    }
    else {
        printf("\n\n     Input data must be first introduced.\n\n");
        printf("      Press 0.\n\n");
        do {
            *k = getchar();
        }while(*k!='0');
        return;
    }
    return;
}

void open_displstmenu(struct statistics *matrix,struct statmulo *matrixml/*,int *observed_data */,int *outgroup,int *n_loci,char *subset_positions,FILE *file_output)
{
    char k[1];
    void open_displstmulo(struct statmulo * /*,int * */,int * /*,int * */,char *,FILE *);
    void open_displstdetmenu(struct statistics * /*,int * */,int *,int *,char *,FILE *);
    
    while(1) {
        #if COMMAND_LINE
        if(file_output) fflush(file_output);
        
        printf("\n\nMENU 1. Observed analysis menu:");
        printf("\n     1.0. Display statistics menu:\n\n");
        
        printf("CHOOSE:\n");
        printf(" 0 - Display summary MULTILOCUS statistics.\n");
        printf(" 1 - Display all detailed statistics.\n");
        printf(" 2 - Back to Observed analysis menu.\n");
		/*
        if(file_output) {
            
            fprintf(file_output,"\n\n     MENU:\n     1. Observed analysis menu:");
            fprintf(file_output,"\n     1.0. Display statistics menu:\n\n");
            
            fprintf(file_output," 0 - Display summary MULTILOCUS statistics.\n");
            fprintf(file_output," 1 - Display all detailed statistics.\n");
            fprintf(file_output," 2 - Back to Observed analysis menu.\n");
        }
		*/
        #endif
        
        do *k = getchar();
        while(*k<'0' || *k>'2');
        
        if(file_output) fprintf(file_output,"OPTION CHOSEN: %c\n\n",*k);

        switch(*k) {
            case '0':
                open_displstmulo(matrixml/*,observed_data */,outgroup/*,n_loci*/,subset_positions,file_output);
                break;
            case '1':
                open_displstdetmenu(matrix/*,observed_data */,outgroup,n_loci,subset_positions,file_output);
                break;
            case '2':
                return;
                break;
        }
    }
    return;
}

void open_displntmenu(struct statistics *matrix,struct statmulo *matrixml/*,int *observed_data */,int *outgroup,int *n_loci,char *subset_positions,FILE *file_output)
{
    char k[1];
    void open_displntmulo(struct statistics *,struct statmulo * /*,int * */,int *,int *,char *,FILE *);
    void open_displntdetmenu(struct statistics * /*,int * */,int *,int *,char *,FILE *);
    
    while(1) {
        if(file_output) fflush(file_output);
        #if COMMAND_LINE
        
        printf("\n\nMENU 1. Observed analysis menu:");
        printf("\n     1.1. Display neutrality tests menu:\n\n");
        
        printf("CHOOSE:\n");
        printf(" 0 - Display MULTILOCUS Neutrality tests.\n");
        printf(" 1 - Display all detailed neutrality tests.\n");
        printf(" 2 - Back to Observed analysis menu.\n");
		/*
        if(file_output) {
            
            fprintf(file_output,"\n\n     MENU:\n     1. Observed analysis menu:");
            fprintf(file_output,"\n     1.1. Display neutrality tests menu:\n\n");
            
            fprintf(file_output," 0 - Display MULTILOCUS Neutrality tests.\n");
            fprintf(file_output," 1 - Display all detailed neutrality tests.\n");
            fprintf(file_output," 2 - Back to Observed analysis menu.\n");
        }
		*/
        #endif
        
        do *k = getchar();
        while(*k<'0' || *k>'2');
        
        if(file_output) fprintf(file_output,"OPTION CHOSEN: %c\n\n",*k);

        switch(*k) {
            case '0':
                open_displntmulo(matrix,matrixml/*,observed_data*/,outgroup,n_loci,subset_positions,file_output);
                break;
            case '1':
                open_displntdetmenu(matrix/*,observed_data */,outgroup,n_loci,subset_positions,file_output);
                break;
            case '2':
                return;
                break;
        }
    }
    return;
}

void open_displstdetmenu(struct statistics *matrix/*,int *observed_data */,int *outgroup,int *n_loci,char *subset_positions,FILE *file_output)
{
    char k[1];
    void open_displstdetgrid(struct statistics * /*,int * */,int *,int *,char *,FILE *);
    void open_displayhistobs(struct statistics * /*,int * */,int *,int * /*,char * */,FILE *,int);
    
    while(1) {
        if(file_output) fflush(file_output);
        #if COMMAND_LINE
        
        printf("\n\nMENU 1. Observed analysis menu:");
        printf("\n     1.0. Display statistics menu:");
        printf("\n     1.0.1. Display detailed statistics menu:\n\n");
        
        printf("CHOOSE:\n");
        printf(" 0 - Display a table with the statistics for each locus.\n");
        printf(" 1 - Display histograms with the statistics.\n");
        printf(" 2 - Back to Display statistics menu.\n");
		/*
        if(file_output) {
            
            fprintf(file_output,"\n\n     MENU:\n     1. Observed analysis menu:");
            fprintf(file_output,"\n     1.0. Display statistics menu:");
            fprintf(file_output,"\n     1.0.1. Display detailed statistics menu:\n\n");
            
            fprintf(file_output," 0 - Display a table with the statistics for each locus.\n");
            fprintf(file_output," 1 - Display histograms with the statistics.\n");
            fprintf(file_output," 2 - Back to Display statistics menu.\n");
        }
		*/
        #endif
        
        do *k = getchar();
        while(*k<'0' || *k>'2');
        
        if(file_output) fprintf(file_output,"\nOPTION CHOSEN: %c\n\n",*k);

        switch(*k) {
            case '0':
                open_displstdetgrid(matrix,/*observed_data,*/outgroup,n_loci,subset_positions,file_output);
                break;
            case '1':
                open_displayhistobs(matrix,/*observed_data,*/outgroup,n_loci/*,subset_positions*/,file_output,0);
                break;
            case '2':
                return;
                break;
        }
    }
    return;
}

void open_displntdetmenu(struct statistics *matrix/*,int *observed_data */,int *outgroup,int *n_loci,char *subset_positions,FILE *file_output)
{
    char k[1];
    void open_displntdetgrid(struct statistics * /*,int * */,int *,int *,char *,FILE *);
    void open_displayhistobs(struct statistics * /*,int * */,int *,int *,/* char *, */FILE *,int);
    
    while(1) {
        #if COMMAND_LINE
        if(file_output) fflush(file_output);
        
        printf("\n\nMENU 1. Observed analysis menu:");
        printf("\n     1.1. Display neutrality tests menu:");
        printf("\n     1.1.1. Display detailed neutrality tests menu:\n\n");
        
        printf("CHOOSE:\n");
        printf(" 0 - Display a table with the tests of neutrality for each locus.\n");
        printf(" 1 - Display histograms with tests of neutrality.\n");
        printf(" 2 - Back to Display Neutrality tests menu.\n");
		/*
        if(file_output) {
            
            fprintf(file_output,"\n\n     MENU:\n     1. Observed analysis menu:");
            fprintf(file_output,"\n     1.1. Display neutrality tests menu:");
            fprintf(file_output,"\n     1.1.1. Display detailed neutrality tests menu:\n\n");
            
            fprintf(file_output," 0 - Display a table with the tests of neutrality for each locus.\n");
            fprintf(file_output," 1 - Display histograms with tests of neutrality.\n");
            fprintf(file_output," 2 - Back to Display Neutrality tests menu.\n");
        }
		*/
        #endif
        
        do *k = getchar();
        while(*k<'0' || *k>'2');
        
        if(file_output) fprintf(file_output,"OPTION CHOSEN: %c\n\n",*k);

        switch(*k) {
            case '0':
                open_displntdetgrid(matrix/*,observed_data */,outgroup,n_loci,subset_positions,file_output);
                break;
            case '1':
                open_displayhistobs(matrix/*,observed_data */,outgroup,n_loci,/*subset_positions,*/file_output,1);
                break;
            case '2':
                return;
                break;
        }
    }
    return;
}

