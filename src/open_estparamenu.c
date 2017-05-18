/*
 *  open_estparamenu.c
 *  ManVa
 *
 *  Created by sebas on September 26th 2005.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void open_paraml_menu(struct globvar **global,struct statistics **matrix,struct statmulo **matrixml,struct var **data,
    struct var2b **inputms,struct statistisim ***matrixsim, struct statistisimmuloc **matrixmlsim,
    struct horizontalstatsml **avgstatloci,struct LRTdist **lrdist/*,struct MLthetasim ***mthetasim*/,FILE *file_output)
{
    char k[1];

    void menu_paraml0_sim(struct globvar **,struct statistics **,struct statmulo **,struct var **, struct var2b **,
        struct statistisim ***, struct statistisimmuloc **, struct horizontalstatsml **,struct LRTdist ** /*,struct MLthetasim ****/,FILE *);    
	void print_display_sim(struct globvar *,struct statistics * /*,struct statmulo * */,struct var *, struct var2b * /*,struct LRTdist * *//*,struct MLthetasim ***/,FILE *);

	void save_simmatrix(struct globvar *,struct statistics *,struct statmulo *,struct var *, struct var2b *,
        struct statistisim **, struct statistisimmuloc *, struct horizontalstatsml *,struct LRTdist * /*,struct MLthetasim ***/,FILE *);
    int open_filematrix(struct globvar **,struct statistics **,struct statmulo **,struct var **, struct var2b **,
        struct statistisim ***, struct statistisimmuloc **, struct horizontalstatsml **,struct LRTdist ** /*,struct MLthetasim ****/,FILE *);
        
    while(1) {
        #if COMMAND_LINE
        printf("\n\nMENU 2. Estimation of variability by Maximum likelihood:\n\n");
        
        printf("CHOOSE:\n");
        printf(" 0 - Estimation of variation under NEUTRAL PANMICTIC model with NULL recombination.\n");
        printf(" 1 - DISPLAY levels of variation under NEUTRAL PANMICTIC model with NULL recombination.\n");
        printf(" 2 - SAVE a file with all the information from observed and/or simulated data.\n");
        printf(" 3 - LOAD a file with all the information from observed and/or simulated data.\n");
        printf(" 4 - Back to Main menu.\n\n");
		/*
        if(file_output) {
            fprintf(file_output,"\n\n     MENU:\n     2. Estimation of varibility by Maximum likelihood:\n\n");            
            fprintf(file_output," 0 - Estimation of variation under NEUTRAL PANMICTIC model with NULL recombination.\n");
            fprintf(file_output," 1 - DISPLAY levels of variation under NEUTRAL PANMICTIC model with NULL recombination.\n");
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
				menu_paraml0_sim(global,matrix,matrixml,data,inputms,matrixsim,matrixmlsim,avgstatloci,lrdist/*,mthetasim*/,file_output);
                break;
            case '1':
				if(global[0][0].ml0done == 1) {
					print_display_sim(*global,*matrix/*,*matrixml*/,*data,*inputms/*,*lrdist *//*,*mthetasim*/,file_output);
				}
				else {
					printf("\n\n   Sorry. Estimation must be first performed.");
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

void print_display_sim(struct globvar *global,struct statistics *matrix/*,struct statmulo *matrixml*/,struct var *data,struct var2b *inputms/*,
	struct LRTdist *lrdist *//*,struct MLthetasim **mthetasim*/,FILE *file_output)
{
	int j;
	void print_seqinfodata_ml0(struct globvar *,struct statistics * /*,struct statmulo * */,struct var *, struct var2b *,FILE *);
	
	double AIC1,AIC2,AIC3,AIC4,AICn,AICmin,sumwAIC;
	int bestAIC;
		
	if(file_output) {
		fprintf(file_output,"\n\n     MENU:\n     2. Estimation of parameters by Maximum likelihood:\n\n");            
		fprintf(file_output,"\n     2.1 Display levels of variation under NEUTRAL PANMICTIC model with NULL recombination.\n");
	}

	print_seqinfodata_ml0(global,matrix/*,matrixml */,data,inputms,file_output);

	/*display the results from ML: */
	printf("\n\n  2. ESTIMATION BY MAXIMUM LIKELIHOOD OF THETA GIVEN THE NUMBER OF SEGREGATING SITES.");
	printf("\n  2.1. DISPLAY levels of variation under NEUTRAL PANMICTIC model with NULL recombination.");
	printf("\n     Calculation using the equation 9.5 from Tavare Theor. Pop. Biol. 26: 119-164.\n");
	printf("\n     (Precission until 5e-12 for each locus).\n");

    printf("\n\nMaximum likelihood estimate considering one theta.\n");
    printf("\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM1/region\n");
    for(j=0;j<data[0].n_loci;j++) {
        if(global[0].dataequalsim) printf("%d:%s",j,matrix[j].gene);
		else printf("%d",j);
		printf("\t%g",inputms[j].factor_chrn);
		printf("\t%d",inputms[j].nsam);
		printf("\t%ld",inputms[j].nsites);
		printf("\t%d",inputms[j].S);
		printf("\t%g\n",inputms[j].thetaml[1]*(double)inputms[j].nsites);
	}
    printf("\ntheta/nt\n %g",data[0].theta1_ml[0]);
    printf("\nML\n %g\n",data[0].ml1);
	
    if(data[0].n_loci > 2) {
        printf("\nMaximum likelihood estimate considering two groups of independent loci");
        
        printf("\nlocus_group1: (%d)\n",data[0].nloci2m[0]);
        printf("\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM2/region\n");
        j = 0;
		do {
			if(inputms[j].thetaml[2] == data[0].theta2_ml[0]) {
				if(global[0].dataequalsim) printf("%d:%s",j,matrix[j].gene);
				else printf("%d",j);
				printf("\t%g",inputms[j].factor_chrn);
				printf("\t%d",inputms[j].nsam);
				printf("\t%ld",inputms[j].nsites);
				printf("\t%d",inputms[j].S);
				printf("\t%g\n",inputms[j].thetaml[2]*(double)inputms[j].nsites);
			}
			j++;
		} while(j<data[0].n_loci);
        printf("\nlocus_group2: (%d)\n",data[0].nloci2m[1]);
		printf("\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM2/region\n");
		j = 0;
        do {
			if(inputms[j].thetaml[2] == data[0].theta2_ml[1]) {
				if(global[0].dataequalsim) printf("%d:%s",j,matrix[j].gene);
				else printf("%d",j);
				printf("\t%g",inputms[j].factor_chrn);
				printf("\t%d",inputms[j].nsam);
				printf("\t%ld",inputms[j].nsites);
				printf("\t%d",inputms[j].S);
				printf("\t%g\n",inputms[j].thetaml[2]*(double)inputms[j].nsites);
				j++;
			}
			j++;
		} while(j<data[0].n_loci);
        printf("\ntheta1ml/nt\ttheta2ml/nt\n%g\t%g\n",data[0].theta2_ml[0],data[0].theta2_ml[1]);
		printf("\nmlgroup1\tmlgroup2\tML\n%g\t%g\t%g\n",data[0].ml2[0],data[0].ml2[1],data[0].ml2[0]+data[0].ml2[1]);
    }

    if(data[0].n_loci > 3) {
        printf("\nMaximum likelihood estimate considering three groups of independent loci");
        printf("\nlocus_group1: (%d)\n",data[0].nloci3m[0]);
		printf("\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM3/region\n");
		j = 0;
        do {
			if(inputms[j].thetaml[3] == data[0].theta3_ml[0]) {
				if(global[0].dataequalsim) printf("%d:%s",j,matrix[j].gene);
				else printf("%d",j);
				printf("\t%g",inputms[j].factor_chrn);
				printf("\t%d",inputms[j].nsam);
				printf("\t%ld",inputms[j].nsites);
				printf("\t%d",inputms[j].S);
				printf("\t%g\n",inputms[j].thetaml[3]*(double)inputms[j].nsites);
			}
			j++;
		} while(j<data[0].n_loci);
        printf("\nlocus_group2: (%d)\n",data[0].nloci3m[1]);
        printf("\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM3/region\n");
		j=0;
        do {
			if(inputms[j].thetaml[3] == data[0].theta3_ml[1]) {
				if(global[0].dataequalsim) printf("%d:%s",j,matrix[j].gene);
				else printf("%d",j);
				printf("\t%g",inputms[j].factor_chrn);
				printf("\t%d",inputms[j].nsam);
				printf("\t%ld",inputms[j].nsites);
				printf("\t%d",inputms[j].S);
				printf("\t%g\n",inputms[j].thetaml[3]*(double)inputms[j].nsites);
			}
			j++;
		} while(j<data[0].n_loci);
        printf("\nlocus_group3: (%d)\n",data[0].nloci3m[2]);
        printf("\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM3/region\n");
		j = 0;
		do {
			if(inputms[j].thetaml[3] == data[0].theta3_ml[2] && j<data[0].n_loci) {
				if(global[0].dataequalsim) printf("%d:%s",j,matrix[j].gene);
				else printf("%d",j);
				printf("\t%g",inputms[j].factor_chrn);
				printf("\t%d",inputms[j].nsam);
				printf("\t%ld",inputms[j].nsites);
				printf("\t%d",inputms[j].S);
				printf("\t%g\n",inputms[j].thetaml[3]*(double)inputms[j].nsites);
			}
			j++;
		} while(j<data[0].n_loci);
        printf("\ntheta1ml/nt\ttheta2ml/nt\ttheta3ml/nt\n%g\t%g\t%g\n",data[0].theta3_ml[0],data[0].theta3_ml[1],data[0].theta3_ml[2]);
        printf("\nmlgroup1\tmlgroup2\tmlgroup3\tML\n%g\t%g\t%g\t%g\n",data[0].ml3[0],data[0].ml3[1],data[0].ml3[2],data[0].ml3[0]+data[0].ml3[1]+data[0].ml3[2]);
    }

    if(data[0].n_loci > 4) {    
        printf("\nMaximum likelihood estimate considering four groups of independent loci");
        printf("\nlocus_group1: (%d)\n",data[0].nloci4m[0]);
        printf("\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM4/region\n");
		j = 0;
        do {
			if(inputms[j].thetaml[4] == data[0].theta4_ml[0] && j<data[0].n_loci) {
				if(global[0].dataequalsim) printf("%d:%s",j,matrix[j].gene);
				else printf("%d",j);
				printf("\t%g",inputms[j].factor_chrn);
				printf("\t%d",inputms[j].nsam);
				printf("\t%ld",inputms[j].nsites);
				printf("\t%d",inputms[j].S);
				printf("\t%g\n",inputms[j].thetaml[4]*(double)inputms[j].nsites);
			}
			j++;
		} while(j<data[0].n_loci);
        printf("\nlocus_group2: (%d)\n",data[0].nloci4m[1]);
        printf("\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM4/region\n");
		j=0;
        do {
			if(inputms[j].thetaml[4] == data[0].theta4_ml[1] && j<data[0].n_loci) {
				if(global[0].dataequalsim) printf("%d:%s",j,matrix[j].gene);
				else printf("%d",j);
				printf("\t%g",inputms[j].factor_chrn);
				printf("\t%d",inputms[j].nsam);
				printf("\t%ld",inputms[j].nsites);
				printf("\t%d",inputms[j].S);
				printf("\t%g\n",inputms[j].thetaml[4]*(double)inputms[j].nsites);
			}
			j++;
		} while(j<data[0].n_loci);
        printf("\nlocus_group3: (%d)\n",data[0].nloci4m[2]);
        printf("\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM4/region\n");
		j=0;
        do {
			if(inputms[j].thetaml[4] == data[0].theta4_ml[2] && j<data[0].n_loci) {
				if(global[0].dataequalsim) printf("%d:%s",j,matrix[j].gene);
				else printf("%d",j);
				printf("\t%g",inputms[j].factor_chrn);
				printf("\t%d",inputms[j].nsam);
				printf("\t%ld",inputms[j].nsites);
				printf("\t%d",inputms[j].S);
				printf("\t%g\n",inputms[j].thetaml[4]*(double)inputms[j].nsites);
			}
			j++;
		} while(j<data[0].n_loci);
        printf("\nlocus_group4: (%d)\n",data[0].nloci4m[3]);
        printf("\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM4/region\n");
		j = 0;
        do {
			if(inputms[j].thetaml[4] == data[0].theta4_ml[3] && j<data[0].n_loci) {
				if(global[0].dataequalsim) printf("%d:%s",j,matrix[j].gene);
				else printf("%d",j);
				printf("\t%g",inputms[j].factor_chrn);
				printf("\t%d",inputms[j].nsam);
				printf("\t%ld",inputms[j].nsites);
				printf("\t%d",inputms[j].S);
				printf("\t%g\n",inputms[j].thetaml[4]*(double)inputms[j].nsites);
			}
			j++;
		} while(j<data[0].n_loci);
        printf("\ntheta1ml/nt\ttheta2ml/nt\ttheta3ml/nt\ttheta4ml/nt\n%g\t%g\t%g\t%g\n",
			data[0].theta4_ml[0],data[0].theta4_ml[1],data[0].theta4_ml[2],data[0].theta4_ml[3]);
        printf("\nmlgroup1\tmlgroup2\tmlgroup3\tmlgroup4\tML\n%g\t%g\t%g\t%g\t%g\n",
			data[0].ml4[0],data[0].ml4[1],data[0].ml4[2],data[0].ml4[3],data[0].ml4[0]+data[0].ml4[1]+data[0].ml4[2]+data[0].ml4[3]);
    }

    printf("\nMaximum likelihood estimate considering each locus independently");
    printf("\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaMl/region\tthetaMl/nt\tml\n");
    for(j=0;j<data[0].n_loci;j++) {
        if(global[0].dataequalsim) printf("%d:%s",j,matrix[j].gene);
		else printf("%d",j);
			printf("\t%g",inputms[j].factor_chrn);
		printf("\t%d",inputms[j].nsam);
		printf("\t%ld",inputms[j].nsites);
		printf("\t%d",inputms[j].S);
		printf("\t%g",inputms[j].thetaml[5]*(double)inputms[j].nsites);
		printf("\t%g",inputms[j].thetaml[5]);
		printf("\t%g\n",inputms[j].mltheta);
	}
    printf("ML:\t\t\t\t\t\t\t%g",data[0].mli_theta);
    puts("\n");

	puts("\nSUMMARY\n");
    puts("\nMaximum likelihood theta values estimated using different models:\n");
    printf("\nM1: %g (%d loci)",data[0].theta1_ml[0],data[0].n_loci);
    if(data[0].n_loci > 2) 
        printf("\nM2: %g (%d loci), %g (%d loci)",
			data[0].theta2_ml[0],data[0].nloci2m[0],data[0].theta2_ml[1],data[0].n_loci-data[0].nloci2m[0]);
    if(data[0].n_loci > 3) 
        printf("\nM3: %g (%d loci), %g (%d loci), %g (%d loci)",
			data[0].theta3_ml[0],data[0].nloci3m[0],data[0].theta3_ml[1],data[0].nloci3m[1],data[0].theta3_ml[2],data[0].n_loci-data[0].nloci3m[0]-data[0].nloci3m[1]);
    if(data[0].n_loci > 4) {
        printf("\nM4: %g (%d loci), %g (%d loci), %g (%d loci), %g (%d loci)",
			data[0].theta4_ml[0],data[0].nloci4m[0],data[0].theta4_ml[1],data[0].nloci4m[1],data[0].theta4_ml[2],data[0].nloci4m[2],data[0].theta4_ml[3],data[0].n_loci-data[0].nloci4m[0]-data[0].nloci4m[1]-data[0].nloci4m[2]);
    }
    puts("\n");

    puts("\nModels:\n");
    printf(" M1: One theta value for all loci.\tML = %g\n",data[0].ml1);
    if(data[0].n_loci > 2) {
        printf(" M2: Two independent theta values.\tML = %g\n",data[0].ml2[2]);
        if(data[0].n_loci > 3) {
            printf(" M3: Three independent theta values.\tML = %g\n",data[0].ml3[3]);
            if(data[0].n_loci > 4) {
                printf(" M4: Four independent theta values.\tML = %g\n",data[0].ml4[4]);
            }
        }
    }
    if(data[0].n_loci > 1) printf(" M%d: %d independent theta values.\tML = %g\n",data[0].n_loci,data[0].n_loci,data[0].mli_theta);
	
	/******** AKAIKE INFORMATION CRITERION (AIC) **********/
	/*calculate AIC: the number of free parameters are not clear for me... I put the highest number...*/
	AIC1 = -(double)2. * (data[0].ml1 + (double)1.);
	if(data[0].n_loci > 2) AIC2 = -(double)2. * (data[0].ml2[2] + (double)3.);
	if(data[0].n_loci > 3) AIC3 = -(double)2. * (data[0].ml3[3] + (double)5.);
	if(data[0].n_loci > 4) AIC4 = -(double)2. * (data[0].ml4[4] + (double)7.);
	if(data[0].n_loci > 1) AICn = -(double)2. * (data[0].mli_theta + (double)data[0].n_loci);
	/*calculate minimum AIC value*/
	AICmin = AIC1; bestAIC = 1;
	if(data[0].n_loci > 2) if(AIC2 < AICmin) {AICmin = AIC2; bestAIC = 2;}
	if(data[0].n_loci > 3) if(AIC3 < AICmin) {AICmin = AIC3; bestAIC = 3;}
	if(data[0].n_loci > 4) if(AIC4 < AICmin) {AICmin = AIC4; bestAIC = 4;}
	if(data[0].n_loci > 1) if(AICn < AICmin) {AICmin = AICn; bestAIC = data[0].n_loci;}
	/*calculate sumwr*/
	sumwAIC = (double)exp(-0.5*(AIC1-AICmin));
	if(data[0].n_loci > 2) sumwAIC += (double)exp(-0.5*(AIC2-AICmin));
	if(data[0].n_loci > 3) sumwAIC += (double)exp(-0.5*(AIC3-AICmin));
	if(data[0].n_loci > 4) sumwAIC += (double)exp(-0.5*(AIC4-AICmin));
	if(data[0].n_loci > 1) sumwAIC += (double)exp(-0.5*(AICn-AICmin));	
	/*print*/
	/*
	printf("\nAKAIKE INFORMATION CRITERION (AIC):\n");
	printf("Warning: These results are orientative. Statistical inference can obtained by Monte Carlo simulations\n");
	printf("\nModel\tlnL\tK\tAIC\tdelta\tweight\n");
	printf("\nM1\t%g\t%d\t%g\t%g\t%g\n",data[0].ml1,1,AIC1,AIC1-AICmin,(double)exp(-0.5*(AIC1-AICmin))/sumwAIC);
	if(data[0].n_loci > 2) printf("\nM2\t%g\t%d\t%g\t%g\t%g\n",data[0].ml2[2],3,AIC2,AIC2-AICmin,(double)exp(-0.5*(AIC2-AICmin))/sumwAIC);
	if(data[0].n_loci > 3) printf("\nM3\t%g\t%d\t%g\t%g\t%g\n",data[0].ml3[3],5,AIC3,AIC3-AICmin,(double)exp(-0.5*(AIC3-AICmin))/sumwAIC);
	if(data[0].n_loci > 4) printf("\nM4\t%g\t%d\t%g\t%g\t%g\n",data[0].ml4[4],7,AIC4,AIC4-AICmin,(double)exp(-0.5*(AIC4-AICmin))/sumwAIC);
	if(data[0].n_loci > 1) printf("\nM%d\t%g\t%d\t%g\t%g\t%g\n",data[0].n_loci,data[0].mli_theta,data[0].n_loci,AICn,AICn-AICmin,(double)exp(-0.5*(AICn-AICmin))/sumwAIC);
	
	if(bestAIC == 1) printf("\nThe selected model is M1\n");
	else if(bestAIC == 2) printf("\nThe selected model is M2\n");
		else if(bestAIC == 3) printf("\nThe selected model is M3\n");
			else if(bestAIC == 4) printf("\nThe selected model is M4\n");
				else if(bestAIC == data[0].n_loci) printf("\nThe selected model is M%d\n",data[0].n_loci);
	printf("Supported models should be those with delta values lower than 2.\n");
	printf("Not supported models should be those with delta values larger than 10.\n");
	*/
	/***********************************************/

    if(data[0].n_loci > 1) {
		printf("\nSTATISTICAL INFERENCE OBTAINED BY MONTE CARLO SIMULATIONS");
		printf("\n\nLikelihood ratio test (LRT) probabilities obtained with %ld coalescent simulations:\n",global[0].nitermatt);
        if(data[0].n_loci >= 2) {
            if(data[0].p1 <= (double)1.)
                printf("\nLR(M1 vs M2) = %g. Probability = %g\n",data[0].LR12,data[0].p1);
            else
                printf("\nLR(M1 vs M2) = %g. Probability not calculated\n",data[0].LR12);
		}
	}
    if(data[0].n_loci > 2) {
        if(data[0].n_loci >= 3) {
            if(data[0].p2 <= (double)1.)
                printf("LR(M2 vs M3) = %g. Probability = %g\n",data[0].LR23,data[0].p2);
            else
                printf("LR(M2 vs M3) = %g. Probability not calculated\n",data[0].LR23);
        }
   }
    if(data[0].n_loci > 3) {
        if(data[0].n_loci >= 4) {
            if(data[0].p3 <= (double)1.)
                printf("LR(M3 vs M4) = %g. Probability = %g\n",data[0].LR34,data[0].p3);
            else
                printf("LR(M3 vs M4) = %g. Probability not calculated\n",data[0].LR34);
        }
    }
    if(data[0].n_loci > 4) {
        if(data[0].p4 <= (double)1.)
            printf("LR(M4 vs M%d) = %g. Probability = %g\n",data[0].n_loci,data[0].LR4l,data[0].p4);
        else
            printf("LR(M4 vs M%d) = %g. Probability not calculated\n",data[0].n_loci,data[0].LR4l);
    }
	
	/*THE BEST MODEL*/
	if(data[0].bestmodel == 1) printf("\nTHE BEST MODEL is M%d. (A single theta/nt per all loci).\n\n",data[0].bestmodel);
	else {
		if(data[0].bestmodel > 1 && data[0].bestmodel < 5) 
			printf("\nTHE BEST MODEL is M%d. (%d different theta/nt distributed in the loci).\n\n",data[0].bestmodel,data[0].bestmodel);
		if(data[0].bestmodel == 5) printf("\nTHE BEST MODEL is M%d. (Each loci has a different theta/nt).\n\n",data[0].n_loci);
	}
	

	/*PRINT TO THE FILE*/
	if(file_output) {
		fprintf(file_output,"\n\n  2. ESTIMATION BY MAXIMUM LIKELIHOOD OF THETA GIVEN THE NUMBER OF SEGREGATING SITES.");
		fprintf(file_output,"\n  2.1. DISPLAY levels of variation under NEUTRAL PANMICTIC model with NULL recombination.");
		fprintf(file_output,"\n     Calculation using the equation 9.5 from Tavare Theor. Pop. Biol. 26: 119-164.\n");
		fprintf(file_output,"\n     (Precission until 5e-12 for each locus).\n");

		fprintf(file_output,"\n\nMaximum likelihood estimate considering one theta.\n");
		fprintf(file_output,"\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM1/region\n");
		for(j=0;j<data[0].n_loci;j++) {
			if(global[0].dataequalsim) fprintf(file_output,"%d:%s",j,matrix[j].gene);
			else fprintf(file_output,"%d",j);
			fprintf(file_output,"\t%g",inputms[j].factor_chrn);
			fprintf(file_output,"\t%d",inputms[j].nsam);
			fprintf(file_output,"\t%ld",inputms[j].nsites);
			fprintf(file_output,"\t%d",inputms[j].S);
			fprintf(file_output,"\t%g\n",inputms[j].thetaml[1]*(double)inputms[j].nsites);
		}
		fprintf(file_output,"\ntheta/nt\n %g",data[0].theta1_ml[0]);
		fprintf(file_output,"\nML\n %g\n",data[0].ml1);
		
		if(data[0].n_loci > 2) {
			fprintf(file_output,"\nMaximum likelihood estimate considering two groups of independent loci");
			
			fprintf(file_output,"\nlocus_group1: (%d)\n",data[0].nloci2m[0]);
			fprintf(file_output,"\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM2/region\n");
			j = 0;
			do {
				if(inputms[j].thetaml[2] == data[0].theta2_ml[0] && j<data[0].n_loci) {
					if(global[0].dataequalsim) fprintf(file_output,"%d:%s",j,matrix[j].gene);
					else fprintf(file_output,"%d",j);
					fprintf(file_output,"\t%g",inputms[j].factor_chrn);
					fprintf(file_output,"\t%d",inputms[j].nsam);
					fprintf(file_output,"\t%ld",inputms[j].nsites);
					fprintf(file_output,"\t%d",inputms[j].S);
					fprintf(file_output,"\t%g\n",inputms[j].thetaml[2]*(double)inputms[j].nsites);
				}
				j++;
			} while(j<data[0].n_loci);
			fprintf(file_output,"\nlocus_group2: (%d)\n",data[0].nloci2m[1]);
			fprintf(file_output,"\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM2/region\n");
			j = 0;
			do {
				if(inputms[j].thetaml[2] == data[0].theta2_ml[1] && j<data[0].n_loci) {
					if(global[0].dataequalsim) fprintf(file_output,"%d:%s",j,matrix[j].gene);
					else fprintf(file_output,"%d",j);
					fprintf(file_output,"\t%g",inputms[j].factor_chrn);
					fprintf(file_output,"\t%d",inputms[j].nsam);
					fprintf(file_output,"\t%ld",inputms[j].nsites);
					fprintf(file_output,"\t%d",inputms[j].S);
					fprintf(file_output,"\t%g\n",inputms[j].thetaml[2]*(double)inputms[j].nsites);
				}
				j++;
			} while(j<data[0].n_loci);
			fprintf(file_output,"\ntheta1ml/nt\ttheta2ml/nt\n%g\t%g\n",data[0].theta2_ml[0],data[0].theta2_ml[1]);
			fprintf(file_output,"\nmlgroup1\tmlgroup2\tML\n%g\t%g\t%g\n",data[0].ml2[0],data[0].ml2[1],data[0].ml2[2]);
		}

		if(data[0].n_loci > 3) {
			fprintf(file_output,"\nMaximum likelihood estimate considering three groups of independent loci");
			fprintf(file_output,"\nlocus_group1: (%d)\n",data[0].nloci3m[0]);
			fprintf(file_output,"\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM3/region\n");
			j = 0;
			do {
				if(inputms[j].thetaml[3] == data[0].theta3_ml[0] && j<data[0].n_loci) {
					if(global[0].dataequalsim) fprintf(file_output,"%d:%s",j,matrix[j].gene);
					else fprintf(file_output,"%d",j);
					fprintf(file_output,"\t%g",inputms[j].factor_chrn);
					fprintf(file_output,"\t%d",inputms[j].nsam);
					fprintf(file_output,"\t%ld",inputms[j].nsites);
					fprintf(file_output,"\t%d",inputms[j].S);
					fprintf(file_output,"\t%g\n",inputms[j].thetaml[3]*(double)inputms[j].nsites);
				}
				j++;
			} while(j<data[0].n_loci);
			fprintf(file_output,"\nlocus_group2: (%d)\n",data[0].nloci3m[1]);
			fprintf(file_output,"\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM3/region\n");
			j=0;
			do {
				if(inputms[j].thetaml[3] == data[0].theta3_ml[1] && j<data[0].n_loci) {
					if(global[0].dataequalsim) fprintf(file_output,"%d:%s",j,matrix[j].gene);
					else fprintf(file_output,"%d",j);
					fprintf(file_output,"\t%g",inputms[j].factor_chrn);
					fprintf(file_output,"\t%d",inputms[j].nsam);
					fprintf(file_output,"\t%ld",inputms[j].nsites);
					fprintf(file_output,"\t%d",inputms[j].S);
					fprintf(file_output,"\t%g\n",inputms[j].thetaml[3]*(double)inputms[j].nsites);
				}
				j++;
			} while(j<data[0].n_loci);
			fprintf(file_output,"\nlocus_group3: (%d)\n",data[0].nloci3m[2]);
			fprintf(file_output,"\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM3/region\n");
			j = 0;
			do {
				if(inputms[j].thetaml[3] == data[0].theta3_ml[2] && j<data[0].n_loci) {
					if(global[0].dataequalsim) fprintf(file_output,"%d:%s",j,matrix[j].gene);
					else fprintf(file_output,"%d",j);
					fprintf(file_output,"\t%g",inputms[j].factor_chrn);
					fprintf(file_output,"\t%d",inputms[j].nsam);
					fprintf(file_output,"\t%ld",inputms[j].nsites);
					fprintf(file_output,"\t%d",inputms[j].S);
					fprintf(file_output,"\t%g\n",inputms[j].thetaml[3]*(double)inputms[j].nsites);
				}
				j++;
			} while(j<data[0].n_loci);
			fprintf(file_output,"\ntheta1ml/nt\ttheta2ml/nt\ttheta3ml/nt\n%g\t%g\t%g\n",data[0].theta3_ml[0],data[0].theta3_ml[1],data[0].theta3_ml[2]);
			fprintf(file_output,"\nmlgroup1\tmlgroup2\tmlgroup3\tML\n%g\t%g\t%g\t%g\n",data[0].ml3[0],data[0].ml3[1],data[0].ml3[2],data[0].ml3[3]);
		}

		if(data[0].n_loci > 4) {    
			fprintf(file_output,"\nMaximum likelihood estimate considering four groups of independent loci");
			fprintf(file_output,"\nlocus_group1: (%d)\n",data[0].nloci4m[0]);
			fprintf(file_output,"\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM4/region\n");
			j = 0;
			do {
				if(inputms[j].thetaml[4] == data[0].theta4_ml[0] && j<data[0].n_loci) {
					if(global[0].dataequalsim) fprintf(file_output,"%d:%s",j,matrix[j].gene);
					else fprintf(file_output,"%d",j);
					fprintf(file_output,"\t%g",inputms[j].factor_chrn);
					fprintf(file_output,"\t%d",inputms[j].nsam);
					fprintf(file_output,"\t%ld",inputms[j].nsites);
					fprintf(file_output,"\t%d",inputms[j].S);
					fprintf(file_output,"\t%g\n",inputms[j].thetaml[4]*(double)inputms[j].nsites);
				}
				j++;
			} while(j<data[0].n_loci);
			fprintf(file_output,"\nlocus_group2: (%d)\n",data[0].nloci4m[1]);
			fprintf(file_output,"\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM4/region\n");
			j=0;
			do {
				if(inputms[j].thetaml[4] == data[0].theta4_ml[1] && j<data[0].n_loci) {
					if(global[0].dataequalsim) fprintf(file_output,"%d:%s",j,matrix[j].gene);
					else fprintf(file_output,"%d",j);
					fprintf(file_output,"\t%g",inputms[j].factor_chrn);
					fprintf(file_output,"\t%d",inputms[j].nsam);
					fprintf(file_output,"\t%ld",inputms[j].nsites);
					fprintf(file_output,"\t%d",inputms[j].S);
					fprintf(file_output,"\t%g\n",inputms[j].thetaml[4]*(double)inputms[j].nsites);
				}
				j++;
			} while(j<data[0].n_loci);
			fprintf(file_output,"\nlocus_group3: (%d)\n",data[0].nloci4m[2]);
			fprintf(file_output,"\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM4/region\n");
			j=0;
			do {
				if(inputms[j].thetaml[4] == data[0].theta4_ml[2] && j<data[0].n_loci) {
					if(global[0].dataequalsim) fprintf(file_output,"%d:%s",j,matrix[j].gene);
					else fprintf(file_output,"%d",j);
					fprintf(file_output,"\t%g",inputms[j].factor_chrn);
					fprintf(file_output,"\t%d",inputms[j].nsam);
					fprintf(file_output,"\t%ld",inputms[j].nsites);
					fprintf(file_output,"\t%d",inputms[j].S);
					fprintf(file_output,"\t%g\n",inputms[j].thetaml[4]*(double)inputms[j].nsites);
				}
				j++;
			} while(j<data[0].n_loci);
			fprintf(file_output,"\nlocus_group4: (%d)\n",data[0].nloci4m[3]);
			fprintf(file_output,"\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaM4/region\n");
			j = 0;
			do {
				if(inputms[j].thetaml[4] == data[0].theta4_ml[3] && j<data[0].n_loci) {
					if(global[0].dataequalsim) fprintf(file_output,"%d:%s",j,matrix[j].gene);
					else fprintf(file_output,"%d",j);
					fprintf(file_output,"\t%g",inputms[j].factor_chrn);
					fprintf(file_output,"\t%d",inputms[j].nsam);
					fprintf(file_output,"\t%ld",inputms[j].nsites);
					fprintf(file_output,"\t%d",inputms[j].S);
					fprintf(file_output,"\t%g\n",inputms[j].thetaml[4]*(double)inputms[j].nsites);
				}
				j++;
			} while(j<data[0].n_loci);
			fprintf(file_output,"\ntheta1ml/nt\ttheta2ml/nt\ttheta3ml/nt\ttheta4ml/nt\n%g\t%g\t%g\t%g\n",
				data[0].theta4_ml[0],data[0].theta4_ml[1],data[0].theta4_ml[2],data[0].theta4_ml[3]);
			fprintf(file_output,"\nmlgroup1\tmlgroup2\tmlgroup3\tmlgroup4\tML\n%g\t%g\t%g\t%g\t%g\n",
				data[0].ml4[0],data[0].ml4[1],data[0].ml4[2],data[0].ml4[3],data[0].ml4[4]);
		}

		fprintf(file_output,"\nMaximum likelihood estimate considering each locus independently");
		fprintf(file_output,"\nname_loci\tfactor_chrn\t#samples\t#sites\tS\tthetaMl/region\tthetaMl/nt\tml\n");
		for(j=0;j<data[0].n_loci;j++) {
			if(global[0].dataequalsim) fprintf(file_output,"%d:%s",j,matrix[j].gene);
			else fprintf(file_output,"%d",j);
			fprintf(file_output,"\t%g",inputms[j].factor_chrn);
			fprintf(file_output,"\t%d",inputms[j].nsam);
			fprintf(file_output,"\t%ld",inputms[j].nsites);
			fprintf(file_output,"\t%d",inputms[j].S);
			fprintf(file_output,"\t%g",inputms[j].thetaml[5]*(double)inputms[j].nsites);
			fprintf(file_output,"\t%g",inputms[j].thetaml[5]);
			fprintf(file_output,"\t%g\n",inputms[j].mltheta);
		}
		fprintf(file_output,"ML:\t\t\t\t\t\t\t%g",data[0].mli_theta);
		fprintf(file_output,"\n");

		fputs("\nSUMMARY\n",file_output);
		fputs("\nMaximum likelihood theta values estimated using different models:\n",file_output);
		fprintf(file_output,"\nM1: %g (%d loci)",data[0].theta1_ml[0],data[0].n_loci);
		if(data[0].n_loci > 2) 
			fprintf(file_output,"\nM2: %g (%d loci), %g (%d loci)",
				data[0].theta2_ml[0],data[0].nloci2m[0],data[0].theta2_ml[1],data[0].n_loci-data[0].nloci2m[0]);
		if(data[0].n_loci > 3) 
			fprintf(file_output,"\nM3: %g (%d loci), %g (%d loci), %g (%d loci)",
				data[0].theta3_ml[0],data[0].nloci3m[0],data[0].theta3_ml[1],data[0].nloci3m[1],data[0].theta3_ml[2],data[0].n_loci-data[0].nloci3m[0]-data[0].nloci3m[1]);
		if(data[0].n_loci > 4) {
			fprintf(file_output,"\nM4: %g (%d loci), %g (%d loci), %g (%d loci), %g (%d loci)",
				data[0].theta4_ml[0],data[0].nloci4m[0],data[0].theta4_ml[1],data[0].nloci4m[1],data[0].theta4_ml[2],data[0].nloci4m[2],data[0].theta4_ml[3],data[0].n_loci-data[0].nloci4m[0]-data[0].nloci4m[1]-data[0].nloci4m[2]);
		}
		fprintf(file_output,"\n");

		fputs("\nModels:\n",file_output);
		fprintf(file_output," M1: One theta value for all loci.\tML = %g\n",data[0].ml1);
		if(data[0].n_loci > 2) {
			fprintf(file_output," M2: Two independent theta values.\tML = %g\n",data[0].ml2[2]);
			if(data[0].n_loci > 3) {
				fprintf(file_output," M3: Three independent theta values.\tML = %g\n",data[0].ml3[3]);
				if(data[0].n_loci > 4) {
					fprintf(file_output," M4: Four independent theta values.\tML = %g\n",data[0].ml4[4]);
				}
			}
		}
		if(data[0].n_loci > 1) fprintf(file_output," M%d: %d independent theta values.\tML = %g\n",data[0].n_loci,data[0].n_loci,data[0].mli_theta);
		
		/******** AKAIKE INFORMATION CRITERION (AIC) **********/
		/*calculate AIC: the number of free parameters are not clear for me... I put the highest number...*/
		AIC1 = -(double)2. * (data[0].ml1 + (double)1.);
		if(data[0].n_loci > 2) AIC2 = -(double)2. * (data[0].ml2[2] + (double)3.);
		if(data[0].n_loci > 3) AIC3 = -(double)2. * (data[0].ml3[3] + (double)5.);
		if(data[0].n_loci > 4) AIC4 = -(double)2. * (data[0].ml4[4] + (double)7.);
		if(data[0].n_loci > 1) AICn = -(double)2. * (data[0].mli_theta + (double)data[0].n_loci);
		/*calculate minimum AIC value*/
		AICmin = AIC1; bestAIC = 1;
		if(data[0].n_loci > 2) if(AIC2 < AICmin) {AICmin = AIC2; bestAIC = 2;}
		if(data[0].n_loci > 3) if(AIC3 < AICmin) {AICmin = AIC3; bestAIC = 3;}
		if(data[0].n_loci > 4) if(AIC4 < AICmin) {AICmin = AIC4; bestAIC = 4;}
		if(data[0].n_loci > 1) if(AICn < AICmin) {AICmin = AICn; bestAIC = data[0].n_loci;}
		/*calculate sumwr*/
		sumwAIC = (double)exp(-0.5*(AIC1-AICmin));
		if(data[0].n_loci > 2) sumwAIC += (double)exp(-0.5*(AIC2-AICmin));
		if(data[0].n_loci > 3) sumwAIC += (double)exp(-0.5*(AIC3-AICmin));
		if(data[0].n_loci > 4) sumwAIC += (double)exp(-0.5*(AIC4-AICmin));
		if(data[0].n_loci > 1) sumwAIC += (double)exp(-0.5*(AICn-AICmin));	
		/*print*/
		/*
		printf("\nAKAIKE INFORMATION CRITERION (AIC):\n");
		printf("Warning: These results are orientative. Statistical inference can obtained by Monte Carlo simulations\n");
		printf("\nModel\tlnL\tK\tAIC\tdelta\tweight\n");
		fputs("\nModel\tlnL\tK\tAIC\tdelta\tweight\n",file_output);
		fprintf(file_output,"\nM1\t%g\t%d\t%g\t%g\t%g\n",data[0].ml1,1,AIC1,AIC1-AICmin,(double)exp(-0.5*(AIC1-AICmin))/sumwAIC);
		if(data[0].n_loci > 2) fprintf(file_output,"\nM2\t%g\t%d\t%g\t%g\t%g\n",data[0].ml2[2],3,AIC2,AIC2-AICmin,(double)exp(-0.5*(AIC2-AICmin))/sumwAIC);
		if(data[0].n_loci > 3) fprintf(file_output,"\nM3\t%g\t%d\t%g\t%g\t%g\n",data[0].ml3[3],5,AIC3,AIC3-AICmin,(double)exp(-0.5*(AIC3-AICmin))/sumwAIC);
		if(data[0].n_loci > 4) fprintf(file_output,"\nM4\t%g\t%d\t%g\t%g\t%g\n",data[0].ml4[4],7,AIC4,AIC4-AICmin,(double)exp(-0.5*(AIC4-AICmin))/sumwAIC);
		if(data[0].n_loci > 1) fprintf(file_output,"\nM%d\t%g\t%d\t%g\t%g\t%g\n",data[0].n_loci,data[0].mli_theta,data[0].n_loci,AICn,AICn-AICmin,(double)exp(-0.5*(AICn-AICmin))/sumwAIC);
		
		if(bestAIC == 1) printf("\nThe selected model is M1\n");
		else if(bestAIC == 2) printf("\nThe selected model is M2\n");
			else if(bestAIC == 3) printf("\nThe selected model is M3\n");
				else if(bestAIC == 4) printf("\nThe selected model is M4\n");
					else if(bestAIC == data[0].n_loci) printf("\nThe selected model is M%d\n",data[0].n_loci);
		printf("Supported models should be those with delta values lower than 2.\n");
		printf("Not supported models should be those with delta values larger than 10.\n");
		*/
		/***********************************************/

		if(data[0].n_loci > 1) {
			fprintf(file_output,"\nSTATISTICAL INFERENCE OBTAINED BY MONTE CARLO SIMULATIONS");
			fprintf(file_output,"\n\nLikelihood ratio test (LRT) probabilities obtained with %ld coalescent simulations:\n",global[0].nitermatt);
			if(data[0].n_loci >= 2) {
				if(data[0].p1 <= (double)1.)
					fprintf(file_output,"\nLR(M1 vs M2) = %g. Probability = %g\n",data[0].LR12,data[0].p1);
				else
					fprintf(file_output,"\nLR(M1 vs M2) = %g. Probability not calculated\n",data[0].LR12);
			}
		}
		if(data[0].n_loci > 2) {
			if(data[0].n_loci >= 3) {
				if(data[0].p2 <= (double)1.)
					fprintf(file_output,"LR(M2 vs M3) = %g. Probability = %g\n",data[0].LR23,data[0].p2);
				else
					fprintf(file_output,"LR(M2 vs M3) = %g. Probability not calculated\n",data[0].LR23);
			}
		}
		if(data[0].n_loci > 3) {
			if(data[0].n_loci >= 4) {
				if(data[0].p3 <= (double)1.)
					fprintf(file_output,"LR(M3 vs M4) = %g. Probability = %g\n",data[0].LR34,data[0].p3);
				else
					fprintf(file_output,"LR(M3 vs M4) = %g. Probability not calculated\n",data[0].LR34);
			}
		}
		if(data[0].n_loci > 4) {
			if(data[0].p4 <= (double)1.)
				fprintf(file_output,"LR(M4 vs M%d) = %g. Probability = %g\n",data[0].n_loci,data[0].LR4l,data[0].p4);
			else
				fprintf(file_output,"LR(M4 vs M%d) = %g. Probability not calculated\n",data[0].n_loci,data[0].LR4l);
		}

		/*THE BEST MODEL*/
		if(data[0].bestmodel == 1) fprintf(file_output,"\nTHE BEST MODEL is M%d. (A single theta/nt per all loci).\n\n",data[0].bestmodel);
		else {
			if(data[0].bestmodel > 1 && data[0].bestmodel < 5) 
				fprintf(file_output,"\nTHE BEST MODEL is M%d. (%d different theta/nt distributed in the loci).\n\n",data[0].bestmodel,data[0].bestmodel);
			if(data[0].bestmodel == 5) fprintf(file_output,"\nTHE BEST MODEL is M%d. (Each loci has a different theta/nt).\n\n",data[0].n_loci);
		}
	}
	
	return;
}

void menu_paraml0_sim(struct globvar **global,struct statistics **matrix,struct statmulo **matrixml,struct var **data,
    struct var2b **inputms,struct statistisim ***matrixsim, struct statistisimmuloc **matrixmlsim,
    struct horizontalstatsml **avgstatloci,struct LRTdist **lrdist/*,struct MLthetasim ***mthetasim*/,FILE *file_output)
{
    char k[1],l[1];
    long int valueld;
	double valuef;
    long int recom;

    int open_menu_seqinfdataml0_sim(struct globvar **,struct statistics ** /*,struct statmulo ** */,struct var **,struct var2b **,FILE *);
	void print_seqinfodata_ml0(struct globvar *,struct statistics * /*,struct statmulo * */,struct var *, struct var2b *,FILE *);
	void save_simmatrix(struct globvar *,struct statistics *,struct statmulo *,struct var *, struct var2b *,
        struct statistisim **, struct statistisimmuloc *, struct horizontalstatsml *,struct LRTdist * /*,struct MLthetasim ***/,FILE *);
	int thetamlS(struct globvar **,struct var **,struct var2b **,struct statistics *,struct LRTdist ** /*,struct MLthetasim ****/,FILE *);
	int reallocvectorsML(/*int *,*/struct LRTdist ** /*,struct MLthetasim ****/,long int, int/*, FILE * */);
	/*void erase_results(struct globvar **,struct statistisim ***, struct statistisimmuloc **,struct horizontalstatsml ** *//*,struct var2b **,struct LRTdist **,*//*struct var **data,FILE *);*/
	
    if(global[0][0].ml0done == 1) {
        printf(" Do you want to overwrite the results from memory (y/n)? \n");
        do *k = getchar();
        while(*k!='y' && *k!='n' && *k!='Y' && *k!='N');
        
        if(*k == 'y' || *k == 'Y') {
			/**/
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
			global[0][0].ml0done = 0;
			/*global[0][0].nlocimatt = 0;*/
			global[0][0].nitermatt = 0;
			global[0][0].mltheta = 0;
			/*
			erase_results(global,matrixsim,matrixmlsim,avgstatloci,lrdist,data,file_output);
			*/
			printf("\n\nResults in the memory have been erased.\n\n");
			if(file_output) 
				fprintf(file_output,"\n\nResults in the memory have been erased.\n\n");
        }
        else return;
    }    
    while(1) {
        #if COMMAND_LINE
        printf("\n\nMENU 2. Estimation of parameters by Maximum likelihood:");
        printf("\n     2.0. Estimation of levels of variation under neutral panmictic model with null recombination:\n\n");
        
        printf("CHOOSE:\n");
        printf(" 0 - Introduce data to estimate parameters.\n");
        printf(" 1 - Run program to estimate parameters.\n");
        printf(" 2 - Back to Main menu.\n\n");
		/*
        if(file_output) {
            fprintf(file_output,"\n\n     MENU:\n     2. Estimation of parameters by Maximum likelihood:");
            fprintf(file_output,"\n     2.0. Estimation of levels of variation under neutral panmictic model with null recombination:\n\n");
            
            fprintf(file_output," 0 - Introduce data to estimate parameters.\n");
            fprintf(file_output," 1 - Run program to estimate parameters.\n");
            fprintf(file_output," 2 - Back to Main menu.\n\n");
        }
		*/
        #endif
    
        do *k = getchar();
        while(*k<'0' || *k>'2');
        
        if(file_output) fprintf(file_output,"OPTION CHOSEN: %c\n\n",*k);

		switch(*k) {
            case '0':
                /*menu for sequence information data*/
                if(global[0][0].dataindata > 0) {
                    print_seqinfodata_ml0(*global,*matrix/*,*matrixml */,*data,*inputms,file_output);
					printf("Do you want to use the current data values (data for loci, i.e., nsamples, etc. introduced before) (y/n)? ");
                    do *l = getchar();
                    while(*l!='y' && *l!='n' && *l!='Y' && *l!='N');
                    if(*l == 'n' || *l == 'n') {
						global[0][0].dataindata = 0;
						if(open_menu_seqinfdataml0_sim(global,matrix,/*matrixml, */data,inputms,file_output) == 0) {
                            return;
                        }
                        print_seqinfodata_ml0(*global,*matrix/*,*matrixml */,*data,*inputms,file_output);
                    }
					else {
						if(global[0][0].dataindata == 1) {
							if(open_menu_seqinfdataml0_sim(global,matrix/*,matrixml */,data,inputms,file_output) == 0) {
								return;
							}
							print_seqinfodata_ml0(*global,*matrix/*,*matrixml */,*data,*inputms,file_output);
						}
					}
                }
                else {
					if(open_menu_seqinfdataml0_sim(global,matrix/*,matrixml*/,data,inputms,file_output) == 0) {
                        return;
                    }
                    print_seqinfodata_ml0(*global,*matrix/*,*matrixml */,*data,*inputms,file_output);
                }
				printf("\n Input data for theta estimation assigned.\n");
				if(file_output) fprintf(file_output,"\n Input data for theta estimation assigned.\n");
                break;
            case '1':
                if(global[0][0].dataforsimt == 0) {
                    printf("Estimation of parameters only run when all data is included.\n");
                    if(file_output)
                        fputs("Estimation of parameters only run when all data is included.\n",
                            file_output);
                    return;
                }
                /*ask for parameters: interval of theta/nt,sections,n_iter,seed1,*/
                do {
                    printf("\nMinimum value of theta/nt? ");
                    scanf(" %lf",&valuef);
                    data[0][0].thetamin_bp = valuef;
                    if(data[0][0].thetamin_bp <= (double)0 && data[0][0].thetamin_bp >= (double)1) puts("Error: the value must be between 0 and 1 (both excluded).");
                }while(data[0][0].thetamin_bp <= (double)0);
                do {
                    printf("\nMaximum value of theta/nt? ");
                    scanf(" %lf",&valuef);
                    data[0][0].thetamax_bp = valuef;
                    if(data[0][0].thetamax_bp <= data[0][0].thetamin_bp && data[0][0].thetamax_bp >= (double)1) puts("Error: the value must be between minimum theta and 1 (both excluded).");
                }while(data[0][0].thetamax_bp <= data[0][0].thetamin_bp && data[0][0].thetamax_bp >= (double)1);

                do {
                    printf("\nNumber of sections of the interval 'thetamax - thetamin'? ");
                    scanf(" %ld",&valueld);
                    data[0][0].steps = valueld;
                    if(data[0][0].steps <= (double)0 && data[0][0].steps >= 1000000) puts("Error: the value must be between 1 and 1000000.");
                }while(data[0][0].steps <= (double)0 && data[0][0].steps >= 1000000);

                do {
                    if(data[0][0].n_loci <= 10) recom = 1000;
                    else if(data[0][0].n_loci <= 100) recom = 1000;
                    else if(data[0][0].n_loci <= 1000) recom = 100;
                    else recom = 100;
                    printf("\nWARNING: This method needs coalescent simulations to calculate the probability for Likelihood Ratios between models.");
                    printf("\nHow many iterations per loci (recommended <= %ld)? ",recom);
                    scanf(" %ld",&valueld);
                    data[0][0].n_iter2 = valueld;
                    if(data[0][0].n_iter2 <= 0) puts("Error: the value must be positive.");
                }while(data[0][0].n_iter2 <= 0);
                do {
                    printf("\nSeed for the pseudorandom numbers (positive value)? ");
                    scanf(" %ld",&valueld);
                    data[0][0].seed2 = valueld;
                    if(data[0][0].seed2 <= 0) puts("Error: the value must be positive.");
                }while(data[0][0].seed2 <= 0);
				
				/*reallocate vector and matrix*/
                if(reallocvectorsML(/*&(global[0][0].maxloc),*/lrdist/*,mthetasim*/,data[0][0].n_iter2,data[0][0].n_loci/*,file_output*/) == 1) {
					return;
				}
				
                /*link to thetamlS PROGRAM*/                
                if(thetamlS(global,data,inputms,*matrix,lrdist/*,mthetasim*/,file_output) == 1) {
                    puts("Error: Estimation of parameters could not be calculated.\n");
                    if(file_output) fputs("Error: Estimation of parameters could not be calculated.\n",file_output);
                    return;
                }
                /*successful simulations*/
                global[0][0].nitermatt = data[0][0].n_iter2;
                global[0][0].nlocimatt = global[0][0].nlocimat = data[0][0].n_loci;
                global[0][0].ml0done = 1;
                
                printf("\n\n Do you want to save this information in a MANVa file (y/n)? ");
                do {
                    *l = getc(stdin);
                }while(*l !='y' && *l !='Y' && *l !='n' && *l!='N');
                if(*l == 'y' || *l == 'Y') {
                    save_simmatrix(*global,*matrix,*matrixml,*data,*inputms,*matrixsim,*matrixmlsim,*avgstatloci,*lrdist/*,*mthetasim*/,file_output);
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

void print_seqinfodata_ml0(struct globvar *global,struct statistics *matrix/*,struct statmulo *matrixml */,struct var *data, struct var2b *inputms,FILE *file_output) 
{
	/*show all values to start estimation*/
    int x;
    
    printf("\n\n Parameters for Maximum likelihood estimation of multilocus theta/s:");
    printf("\n Parameters related to sequence information data:\n");
    printf("\n Number of loci: %d\n",data[0].n_loci);
    
    if(global[0].dataindata == 2 || global[0].dataindata == 3) 
		printf("\nnloci\tfactor_chrn\tnsamples\tnsites\tSegregating_sites\n");
    else printf("\nnloci\tfactor_chrn\tnsamples\tnsites\n");
	for(x=0;x<data[0].n_loci;x++) {
        if(global[0].dataequalsim == 1)
            printf("%d:%s\t",x,matrix[x].gene);
        else
             printf("locus[%d]\t",x);
        printf("%g\t",inputms[x].factor_chrn);
        printf("%d\t",inputms[x].nsam);
        printf("%ld\t",inputms[x].nsites);
        if(global[0].dataindata == 2 || global[0].dataindata == 3) {
			printf("%d",inputms[x].S);
		}
		printf("\n");
    }
    printf("\n\n");

    if(file_output) {
        fputs("\n\n Parameters for Maximum likelihood estimation of multilocus theta/s:",file_output);
        fputs("\n Parameters related to sequence information data:\n",file_output);
        fprintf(file_output,"\n Number of loci: %d\n",data[0].n_loci);
        if(global[0].dataindata == 2 || global[0].dataindata == 3) 
			fputs("\nnloci\tfactor_chrn\tnsamples\tnsites\tSegregating_sites\n",file_output);
		else fputs("\nnloci\tfactor_chrn\tnsamples\tnsites\n",file_output);
        for(x=0;x<data[0].n_loci;x++) {
            if(global[0].dataequalsim == 1)
                fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
            else
                fprintf(file_output,"locus[%d]\t",x);
            fprintf(file_output,"%g\t",inputms[x].factor_chrn);
            fprintf(file_output,"%d\t",inputms[x].nsam);
            fprintf(file_output,"%ld\t",inputms[x].nsites);
            if(global[0].dataindata == 2 || global[0].dataindata == 3) {
				fprintf(file_output,"%d\t",inputms[x].S);
			}
			fprintf(file_output,"\n");
        }
        fputs("\n\n",file_output);
    }
    return;
}

void print_seqinfoobsdata_ml0(struct globvar *global,struct statistics *matrix/*,struct statmulo *matrixml,struct var *data, struct var2b *inputms*/,FILE *file_output) 
{
	/*show all values to start estimation*/
    int x;
    
    printf("\n Parameters related to observed sequence information data:\n");
    printf("\n Number of loci: %d\n",global[0].n_loci);
    
	printf("\nnloci\tfactor_chrn\tnsamples\tnsites\tSegregating_sites\n");
	for(x=0;x<global[0].n_loci;x++) {
		printf("%d:%s\t",x,matrix[x].gene);
        printf("%g\t",matrix[x].factor_chrn);
        printf("%d\t",matrix[x].nsamples);
        printf("%.2f\t",matrix[x].nsites);
		printf("%d",matrix[x].biallsites);
		printf("\n");
    }
    printf("\n\n");

    if(file_output) {
        fputs("\n Parameters related to observed sequence information data:\n",file_output);
        fprintf(file_output,"\n Number of loci: %d\n",global[0].n_loci);
		fputs("\nnloci\tfactor_chrn\tnsamples\tnsites\tSegregating_sites\n",file_output);
        for(x=0;x<global[0].n_loci;x++) {
                fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
            fprintf(file_output,"%g\t",matrix[x].factor_chrn);
            fprintf(file_output,"%d\t",matrix[x].nsamples);
            fprintf(file_output,"%.2f\t",matrix[x].nsites);
			fprintf(file_output,"%d\t",matrix[x].biallsites);
			fprintf(file_output,"\n");
        }
        fputs("\n\n",file_output);
    }
    return;
}

int open_menu_seqinfdataml0_sim(struct globvar **global,struct statistics **matrix/*,struct statmulo **matrixml*/,struct var **data,struct var2b **inputms, FILE *file_output)
{
    char ll[1];
	char mm[1];
    int valuei;
    long int valueli;
	double valuef;
	int x;
	void print_seqinfoobsdata_ml0(struct globvar *,struct statistics * /*,struct statmulo *,struct var *, struct var2b * */,FILE *);
	void print_seqinfodata_ml0(struct globvar *,struct statistics * /*,struct statmulo * */,struct var *, struct var2b *,FILE *);
    
	*mm = *ll = '0';

	if(global[0][0].dataindata &&  global[0][0].dataequalsim == 0) {
		print_seqinfodata_ml0(*global,*matrix/*,*matrixml*/,*data,*inputms,file_output);
		printf("Do you want to use these current values to estimate multilocus theta (y/n)? ");
		do *mm = getchar();
		while(*mm!='y' && *mm!='n' && *mm!='Y' && *mm!='N');
	}
	else *mm = 'n';

	if((global[0][0].observed_data && (*mm == 'n' || *mm == 'N')) || 
	   (global[0][0].observed_data && (*mm == 'y' || *mm == 'Y') && global[0][0].dataindata == 1 && global[0][0].dataequalsim)) {
		print_seqinfoobsdata_ml0(*global,*matrix/*,*matrixml,*data,*inputms */,file_output);
		printf("Do you want to use the values of the observed data to estimate multilocus theta (y/n)? ");
		do *ll = getchar();
		while(*ll!='y' && *ll!='n' && *ll!='Y' && *ll!='N');
	}
	else {
		if((*mm == 'y' || *mm == 'Y') && global[0][0].dataindata > 1) *ll = '\0';
		else *ll = 'n';
	}
		
	if(*ll == 'n' || *ll == 'N') {
		/*no use current data*/
		x = 0;
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
		if((inputms[0][x].freq = (double *) realloc(inputms[0][x].freq,sizeof(double))) == 0) {
			puts("\nError: memory not reallocated. main.12 \n");
			return(0);
		}
		memset(inputms[0][x].freq,'\0',sizeof(double));
		if((inputms[0][x].factor_pop = (double *) realloc(inputms[0][x].factor_pop,sizeof(double))) == 0) {
			puts("\nError: memory not reallocated. main.8 \n");
			return(0);
		}
		memset(inputms[0][x].factor_pop,'\0',sizeof(double));
		
		inputms[0][x].nintn = 0;
		inputms[0][x].npop = 0;

		for(x=1;x<data[0][0].n_loci;x++) {
			free(inputms[0][x].factor_pop);
			free(inputms[0][x].nrec);
			free(inputms[0][x].npast);
			free(inputms[0][x].tpast);
			free(inputms[0][x].freq);
		}
				
		global[0][0].montecarlo_sim = 0;
		global[0][0].nlocimat = 0;
		global[0][0].nitermat = 0;
		global[0][0].dataforsim = 0;
		global[0][0].dataindata = 0;
		global[0][0].dataequalsim = 0;
		global[0][0].ml0done = 0;
		global[0][0].nlocimatt = 0;
		global[0][0].nitermatt = 0;
		global[0][0].dataforsimt = 0;
		global[0][0].mltheta = 0;		
	}
	
	if(*ll == 'Y' || *ll == 'y') {
		/*use observed data, assumed dataequalsim=1 if dataindata>1*/
		if(global[0][0].dataindata == 0) {
			data[0][0].n_loci = global[0][0].n_loci;
			if((inputms[0] = (struct var2b *)realloc(inputms[0],(data[0][0].n_loci)*sizeof(struct var2b))) == 0) {
				puts("\nError: memory not reallocated. seqinfodata.4 \n");
				return(0);
			}			
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
			}

            for(x=0;x<data[0][0].n_loci;x++) {
                inputms[0][x].factor_chrn = matrix[0][x].factor_chrn;
                inputms[0][x].nsam = matrix[0][x].nsamples;
                inputms[0][x].nsites = (long int)matrix[0][x].nsites;
				inputms[0][x].S = matrix[0][x].biallsites;
            }
			global[0][0].dataindata = 2;
		}
		if(global[0][0].dataindata == 1) {
            for(x=0;x<data[0][0].n_loci;x++) {
				inputms[0][x].S = matrix[0][x].biallsites;
            }
			global[0][0].dataindata = 3;
		}
		global[0][0].dataequalsim = 1;
		global[0][0].nlocimatt = global[0][0].nlocimat = global[0][0].n_loci;
		global[0][0].dataforsimt = 1;
	}
	if(*ll == 'N' || *ll == 'n') {
		/*no observed data*/
		if(global[0][0].dataindata == 0) {
			do {
				puts("\nIndicate the number of loci:\n");
				scanf(" %d",&valuei);
				data[0][0].n_loci = valuei;
				if(data[0][0].n_loci <= 0) puts("Error. The value must be positive.");
			}while(data[0][0].n_loci <= 0);

			if((inputms[0] = (struct var2b *)realloc(inputms[0],data[0][0].n_loci*sizeof(struct var2b))) == 0) {
				puts("\nError: memory not reallocated. main.4 \n");
				exit(0);
			}

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
			}

			puts("\nIndicate the correction factor for chromosome population size (i.e., autosomal is 1.0, X chromosome 0.75, Y chromosome 0.25):\n");
			for(x=0;x<data[0][0].n_loci;x++) {
				do{
					printf("\n   Indicate the correction factor for chromosome population size for the locus %d: ",x);
					scanf(" %lg",&valuef);
					inputms[0][x].factor_chrn = valuef;
					if(inputms[0][x].factor_chrn <= (double)0) puts("Error. The value must be higher than 0.");
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
			puts("\nIndicate the number of sites:\n");
			for(x=0;x<data[0][0].n_loci;x++) {
				do {
					printf("\n   Indicate the number of sites for the locus %d: ",x);
					scanf(" %ld",&valueli);
					inputms[0][x].nsites = valueli;
					if(inputms[0][x].nsites <= 0) puts("Error. The value must be positive.");
				}while(inputms[0][x].nsites <= 0);
			}
			puts("\nIndicate the number of segregating sites per locus:\n");
			for(x=0;x<data[0][0].n_loci;x++) {
				do {
					printf("\n   Indicate the number of segregating sites for the locus %d: ",x);
					scanf(" %d",&valuei);
					inputms[0][x].S = valuei;
					if(inputms[0][x].S < 0) puts("Error. The value must be positive or zero.");
				}while(inputms[0][x].S < 0);
			}
			for(x=0;x<data[0][0].n_loci;x++) {
				inputms[0][x].npop = 1;
				inputms[0][x].nintn = 0;
			}
			global[0][0].dataindata = 2;
		}
		if(global[0][0].dataindata == 1) {
			puts("\nIndicate the number of segregating sites per locus:\n");
			for(x=0;x<data[0][0].n_loci;x++) {
				do {
					printf("\n   Indicate the number of segregating sites for the locus %d: ",x);
					scanf(" %d",&valuei);
					inputms[0][x].S = valuei;
					if(inputms[0][x].S < 0) puts("Error. The value must be positive or zero.");
				}while(inputms[0][x].S < 0);
			}
			global[0][0].dataindata = 3;
		}
		global[0][0].dataequalsim = 0;
		global[0][0].nlocimatt = global[0][0].nlocimat = data[0][0].n_loci;
		global[0][0].dataforsimt = 1;
	}
	
    return 1;
}
