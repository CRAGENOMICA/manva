/*
 *  print_parsim.c
 *  MuLoNeTests
 *
 *  Created by sebas on Sat Apr 05 2003.
 *
 */
 
#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>

void print_datasimulations(int dataindata,struct var *data,struct var2b *inputms, FILE *file_output,int dataobsequalsim,struct statistics *matrix)
{
    int x,y,z;

    /*print ALL parameters of simulations (var & var2)*/
    printf("\n\n Parameters for coalescent simulations:\n");
	if(dataindata == 1 || dataindata == 3) {
		printf("\n Number of iterations per loci: %ld",data[0].n_iter);
		printf("\n Seed: %ld",data[0].seed1);
	}
    printf("\n Number of loci: %d\n",data[0].n_loci);
    if(dataindata == 1 || dataindata == 3) {
		if(data[0].time_spec != -1.) {
			printf(" Time of divergence to the ancestor in 2N gen.: %g\n",data[0].time_spec);
			printf(" Multiple hits are considered in simulations, but not counted in the analysis.\n");
		}
		if(data[0].tlimit < 9999)
			printf("\n The length of the sample tree is limited to a length of 4N x %g generations.\n",data[0].tlimit);
    }
    printf("\n Note that the simulations will include the number of sites plus positions with mhits.\n");

	if(data[0].no_rec_males == 1) puts("\n Recombination supressed in males\n");
    if(dataindata == 1 || dataindata == 3) printf("\n\nnloci\tfactor_chrn\tnsamples\tnsites\tnmhits\ttheta(4Nu)_locus\trec(4Nr)_locus\ttrans/transv\n");
	else printf("\n\nnloci\tfactor_chrn\tnsamples\tnsites\tnmhits\n");
    for(x=0;x<data[0].n_loci;x++) {
        if(dataobsequalsim == 1)
            printf("%d:%s\t",x,matrix[x].gene);
        else
             printf("\tlocus[%d]\t",x);
        printf("%g\t",inputms[x].factor_chrn);
        printf("%d\t",inputms[x].nsam);
        printf("%ld\t",inputms[x].nsites);
        printf("%ld\t",inputms[x].nmhits);
        if(dataindata == 1 || dataindata == 3) {
			printf("%g\t",inputms[x].theta * (double)1/inputms[x].factor_chrn);
			printf("%g\t",inputms[x].r * (double)1/inputms[x].factor_chrn);
			printf("%g\t",inputms[x].ratio_sv);
		}
		printf("\n");
    }
	if(dataindata == 1 || dataindata == 3) {
		printf("\n\nNumber of evolutive processes ocurred in different loci: %d\n",data[0].n_events);
		if(data[0].n_events > 1) {
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
		}
		if(data[0].neutrality) {
			if(data[0].n_events==1) printf("\n\n Neutral panmictic model process in Coalescent simulations.");
			else {
				if(data[0].n_events==2)  printf("\n\n Neutral panmictic model process except in loci indicated below.\n");
				else {
					printf("\n\n Neutral panmictic model process in loci indicated below:\n");
					for(x=0;x<data[0].n_loci;x++) {
						if(inputms[x].neutrality) printf("\n   %d:%s",x,matrix[x].gene);
					}
				}
			}
		}
		if(data[0].changesizepop) {
			printf("\n\n Change in population size. Logistic equation within events:");		
			for(z=1;z<=data[0].changesizepop;z++) {	
				if(data[0].changesizepop > 1) printf("\n\n\tProcess %d\n",z);
				x=0;
				while(inputms[x].changesizepop != z) x++;
				printf("\n\tNumber of events: %d\n",inputms[x].nintn);
				printf("\n\t\tInitial_N\tFinal_N\tDuration(in 4N gen)\n");
				for(y=0;y<inputms[x].nintn;y++) {
					printf("\t\t%g\t",inputms[x].nrec[y]);
					printf("%g\t",inputms[x].npast[y]);
					printf("%g\n",inputms[x].tpast[y]);
				}
				printf("\n\t\tts: %g\n",inputms[x].ts);
				if(data[0].n_events==1) printf("\n\n All loci follow this model.");
				else {
					printf("\n\t\tLoci following these parameters:");
					for(x=0;x<data[0].n_loci;x++) {
						if(inputms[x].changesizepop == z) printf("\n\t\t\t%d:%s",x,matrix[x].gene);
					}
				}
			}
		}
		if(data[0].subdivision) {
			printf("\n\n Subdivision. Equal migration and extinction rates among demes:");
			printf("\n    The level of variation and recombination indicated is within deme.\n");
			printf("\n    In order to do correct comparisons using the island isotropic migration model,");
			printf("\n    the level of variation and recombination within deme is divided by the number of demes when simulations are performed.\n");
			printf("\n    All samples are in the same deme.");
			for(z=1;z<=data[0].subdivision;z++) {	
				if(data[0].subdivision > 1) printf("\n\n\tProcess %d\n",z);
				x=0;
				while(inputms[x].subdivision != z) x++;
				printf("\n    Total number of demes: %d",inputms[x].npop);
				printf("\n    Migration among demes (4Nm): %g",inputms[x].mig_rate);
				printf("\n    Population size N (in relation to initial N) of every population:\n\t");
				for(y=0;y<inputms[x].npop;y++) 
					printf("%g\t",inputms[x].factor_pop[y]);
				if(data[0].n_events==1) printf("\n\n All loci follow this model.");
				else {
					printf("\n    Loci following these parameters:");
					for(x=0;x<data[0].n_loci;x++) {
						if(inputms[x].subdivision == z) printf("\n    %d:%s",x,matrix[x].gene);
					}
				}
			}
		}
		if(data[0].split_pop) {
			printf("\n\n Refugia model. Split populattion in several isolated demes and rejoin:");
			for(z=1;z<=data[0].split_pop;z++) {	
				if(data[0].split_pop > 1) printf("\n\n\tProcess %d\n",z);
				x=0;
				while(inputms[x].split_pop != z) x++;
				printf("\n\n    Split population in the past in %d refuges:",inputms[x].npop);
				printf("\n    Duration (in 4N generations) until the split: %g",inputms[x].time_split);
				printf("\n    Duration (in 4N generations) from the split until rejoining in a single population: %g",inputms[x].time_scoal);
				printf("\n    Migration among refuges (4Nm): %g",inputms[x].mig_rate);
				printf("\n    Population size N (in relation to initial N) of every population in the splitted phase:\n\t");
				for(y=0;y<inputms[x].npop;y++) 
					printf("%g\t",inputms[x].factor_pop[y]);
				printf("\n     Population size N (in relation to the initial N) of the ancestral population: %g",
					inputms[x].factor_anc);
				printf("\n     Average contribution of every population (in the splitted phase) to the present population:\n\t");
				for(y=0;y<inputms[x].npop;y++) 
					printf("%g\t",inputms[x].freq[y]);
				if(data[0].n_events==1) printf("\n\n All loci follow this model.");
				else {
					printf("\n    Loci following these parameters:");
					for(x=0;x<data[0].n_loci;x++) {
						if(inputms[x].split_pop == z) printf("\n    %d:%s",x,matrix[x].gene);
					}
				}
			}

		}
		if(data[0].pselection) {
			printf("\n\n Positive selection process");
			for(z=1;z<=data[0].pselection;z++) {	
				if(data[0].pselection > 1) printf("\n\n\tProcess %d\n",z);
				x=0;
				while(inputms[x].pselection != z) x++;
				printf("\n    Population size (in number of individuals N): %ld",inputms[x].pop_size);
				printf("\n    Selection coefficient 4Ns: %g",inputms[x].pop_sel);
				printf("\n    Nucleotide position selected: %ld",inputms[x].sel_nt+1);
				printf("\n    Time (in 4N generations) since the selective process ended: %g\n",inputms[x].sinit);
				if(data[0].n_events==1) printf("\n\n All loci follow this model.");
				else {
					printf("\n    Loci following these parameters:");
					for(y=0;y<data[0].n_loci;y++) {
						if(inputms[y].pselection == z) printf("\n    %d:%s",y,matrix[y].gene);
					}
				}
			}
		}
		printf("\n\n");
    }
	
    if(file_output) {
        fputs("\n\n Parameters used in coalescent simulations:\n",file_output);
        if(dataindata == 1 || dataindata == 3) {
			fprintf(file_output,"\n Number of iterations per loci: %ld",data[0].n_iter);
			fprintf(file_output,"\n Seed: %ld",data[0].seed1);
		}
        fprintf(file_output,"\n Number of loci: %d\n",data[0].n_loci);
        if(dataindata == 1 || dataindata == 3) {
			if(data[0].time_spec != -1.) {
				fprintf(file_output," Time of divergence to the ancestor in 2N gen.: %g\n",data[0].time_spec);
				fputs(" Multiple hits are considered in simulations, but not counted in the analysis.\n",file_output);
			}
			if(data[0].tlimit < 9999)
				fprintf(file_output,
					" The length of the sample tree is limited to a length of 4x%g N generations.\n",data[0].tlimit);
        }
		
		fputs("\n Note that the simulations will include the number of sites plus positions with mhits.\n",file_output);

		if(data[0].no_rec_males == 1) fputs("\n Recombination supressed in males\n",file_output);
        if(dataindata == 1 || dataindata == 3) fputs("\nnloci\tfactor_chrn\tnsamples\tnsites\tnmhits\ttheta(4Nu)_locus\trec(4Nr)_locus\ttrans/transv\n",file_output);
        else fputs("\nnloci\tfactor_chrn\tnsamples\tnsites\tnmhits\n",file_output);
		for(x=0;x<data[0].n_loci;x++) {
            if(dataobsequalsim == 1)
                fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
            else
                fprintf(file_output,"\tlocus[%d]\t",x);
            fprintf(file_output,"%g\t",inputms[x].factor_chrn);
            fprintf(file_output,"%d\t",inputms[x].nsam);
            fprintf(file_output,"%ld\t",inputms[x].nsites);
            fprintf(file_output,"%ld\t",inputms[x].nmhits);
            if(dataindata == 1 || dataindata == 3) {
				fprintf(file_output,"%g\t",inputms[x].theta * (double)1/inputms[x].factor_chrn);
				fprintf(file_output,"%g\t",inputms[x].r * (double)1/inputms[x].factor_chrn);
				fprintf(file_output,"%g\t",inputms[x].ratio_sv);
			}
			fprintf(file_output,"\n");
        }
		if(dataindata == 1 || dataindata == 3) {
			fprintf(file_output,"\n\nNumber of evolutive processes ocurred in different loci: %d\n\n",data[0].n_events);
			if(data[0].n_events > 1) {
				fprintf(file_output,"**********************************************************\n");
				fprintf(file_output,"**********************************************************\n");
				fprintf(file_output,"* WARNING: Selective processes affect specific loci.     *\n");
				fprintf(file_output,"* A demographic process should affect the entire genome. *\n");
				fprintf(file_output,"*  The use of different demographic processes            *\n");
				fprintf(file_output,"*  in different loci could not be justified.             *\n");
				fprintf(file_output,"*  A selective process is modeled under                  *\n");
				fprintf(file_output,"*  a stationary population background.                   *\n");
				fprintf(file_output,"**********************************************************\n");
				fprintf(file_output,"**********************************************************\n");
			}
			if(data[0].neutrality) {
				if(data[0].n_events==1) fprintf(file_output,"\n\n Neutral panmictic model process in Coalescent simulations.");
				else {
					if(data[0].n_events==2)  fprintf(file_output,"\n\n Neutral panmictic model process except in loci indicated below.\n");
					else {
						fprintf(file_output,"\n\n Neutral panmictic model process in loci indicated below:\n");
						for(x=0;x<data[0].n_loci;x++) {
							if(inputms[x].neutrality) fprintf(file_output,"\n   %d:%s",x,matrix[x].gene);
						}
					}
				}
			}
			if(data[0].changesizepop) {
				fprintf(file_output,"\n\n Change in population size. Logistic equation within events:");		
				for(z=1;z<=data[0].changesizepop;z++) {	
					if(data[0].changesizepop > 1) fprintf(file_output,"\n\n\tProcess %d\n",z);
					x=0;
					while(inputms[x].changesizepop != z) x++;
					fprintf(file_output,"\n\n\tNumber of events: %d\n",inputms[x].nintn);
					fprintf(file_output,"\n\tInitial_N\tFinal_N\tDuration(in 4N gen)\n");
					for(y=0;y<inputms[x].nintn;y++) {
						fprintf(file_output,"\t%g\t",inputms[x].nrec[y]);
						fprintf(file_output,"%g\t",inputms[x].npast[y]);
						fprintf(file_output,"%g\n",inputms[x].tpast[y]);
					}
					fprintf(file_output,"\n\t\tts: %g\n",inputms[x].ts);
					if(data[0].n_events==1) fprintf(file_output,"\n\n All loci follow this model.");
					else {
						fprintf(file_output,"\n\tLoci following these parameters:");
						for(x=0;x<data[0].n_loci;x++) {
							if(inputms[x].changesizepop == z) fprintf(file_output,"\n\t\t%d:%s",x,matrix[x].gene);
						}
					}
				}
			}
			if(data[0].subdivision) {
				fprintf(file_output,"\n\n Subdivision. Equal migration and extinction rates among demes:");
				fprintf(file_output,"\n    The level of variation and recombination indicated is within deme.\n");
				fprintf(file_output,"\n    In order to do correct comparisons using the island isotropic migration model,");
				fprintf(file_output,"\n    the level of variation and recombination within deme is divided by the number of demes when simulations are performed.\n");
				fprintf(file_output,"\n    All samples are in the same deme.");
				for(z=1;z<=data[0].subdivision;z++) {	
					if(data[0].subdivision > 1) fprintf(file_output,"\n\n\tProcess %d\n",z);
					x=0;
					while(inputms[x].subdivision != z) x++;
					fprintf(file_output,"\n    Total number of demes: %d",inputms[x].npop);
					fprintf(file_output,"\n    Migration among demes (4Nm): %g",inputms[x].mig_rate);
					fprintf(file_output,"\n    Population size N (in relation to initial N) of every population:\n\t");
					for(y=0;y<inputms[x].npop;y++) 
						fprintf(file_output,"%g\t",inputms[x].factor_pop[y]);
					if(data[0].n_events==1) fprintf(file_output,"\n\n All loci follow this model.");
					else {
					fprintf(file_output,"\n\tLoci following these parameters:");
						for(x=0;x<data[0].n_loci;x++) {
							if(inputms[x].subdivision == z) fprintf(file_output,"\n\t%d:%s",x,matrix[x].gene);
						}
					}
				}
			}
			if(data[0].split_pop) {
				fprintf(file_output,"\n\n Refugia model. Split populattion in several isolated demes and rejoin:");
				for(z=1;z<=data[0].split_pop;z++) {	
					if(data[0].split_pop > 1) fprintf(file_output,"\n\n\tProcess %d\n",z);
					x=0;
					while(inputms[x].split_pop != z) x++;
					fprintf(file_output,"\n\n Split population in the past in %d refuges:",inputms[x].npop);
					fprintf(file_output,"\n   Duration (in 4N generations) until the split: %g",inputms[x].time_split);
					fprintf(file_output,"\n   Duration (in 4N generations) from the split until rejoining in a single population: %g",inputms[x].time_scoal);
					fprintf(file_output,"\n   Migration among refuges (4Nm): %g",inputms[x].mig_rate);
					fprintf(file_output,"\n   Population size N (in relation to initial N) of every population in the splitted phase:\n\t");
					for(y=0;y<inputms[x].npop;y++) 
						fprintf(file_output,"%g\t",inputms[x].factor_pop[y]);
					fprintf(file_output,"\n     Population size N (in relation to the initial N) of the ancestral population: %g",
						inputms[x].factor_anc);
					fprintf(file_output,"\n    Average contribution of every population (in the splitted phase) to the present population:\n\t");
					for(y=0;y<inputms[x].npop;y++) 
						fprintf(file_output,"%g\t",inputms[x].freq[y]);
					if(data[0].n_events==1) fprintf(file_output,"\n\n All loci follow this model.");
					else {
					fprintf(file_output,"\n\tLoci following these parameters:");
						for(x=0;x<data[0].n_loci;x++) {
							if(inputms[x].split_pop == z) fprintf(file_output,"\n\t%d:%s",x,matrix[x].gene);
						}
					}
				}

			}
			if(data[0].pselection) {
				fprintf(file_output,"\n\n Positive selection process");
				for(z=1;z<=data[0].pselection;z++) {	
					if(data[0].pselection > 1) fprintf(file_output,"\n\n\tProcess %d\n",z);
					x=0;
					while(inputms[x].pselection != z) x++;
					fprintf(file_output,"   Population size (in number of individuals N): %ld",inputms[x].pop_size);
					fprintf(file_output,"\n    Selection coefficient 4Ns: %g",inputms[x].pop_sel);
					fprintf(file_output,"\n    Nucleotide position selected: %ld",inputms[x].sel_nt+1);
					fprintf(file_output,"\n    Time (in 4N generations) since the selective process ended: %g\n",inputms[x].sinit);
					if(data[0].n_events==1) fprintf(file_output,"\n\n All loci follow this model.");
					else {
					fprintf(file_output,"\n\tLoci following these parameters:");
						for(y=0;y<data[0].n_loci;y++) {
							if(inputms[y].pselection == z) fprintf(file_output,"\n\t%d:%s",y,matrix[y].gene);
						}
					}
				}
			}
			fputs("\n\n",file_output);
		}
    }
}

void print_seqinfodata_sim(int dataindata,struct var *data,struct var2b *inputms, FILE *file_output,int dataobsequalsim,struct statistics *matrix)
{
    int x;
    
    printf("\n\n Parameters for coalescent simulations:");
    printf("\n Parameters related to sequence information data:\n");
    printf("\n Number of loci: %d\n",data[0].n_loci);
    
    printf("\n Note that the simulations will include the number of sites plus positions with mhits.\n");

	if((dataindata == 1 || dataindata == 3) && data[0].no_rec_males == 1) puts("\n Recombination supressed in males\n");
    if(dataindata == 1 || dataindata == 3) printf("\nnloci\tfactor_chrn\tnsamples\tnsites\tnmhits\ttheta(4Nu)_locus\trec(4Nr)_locus\ttrans/transv\n");
    else printf("\nnloci\tfactor_chrn\tnsamples\tnsites\tnmhits\n");
	for(x=0;x<data[0].n_loci;x++) {
        if(dataobsequalsim == 1)
            printf("%d:%s\t",x,matrix[x].gene);
        else
             printf("locus[%d]\t",x);
        printf("%g\t",(double)inputms[x].factor_chrn);
        printf("%d\t",inputms[x].nsam);
        printf("%ld\t",inputms[x].nsites);
        printf("%ld\t",inputms[x].nmhits);
        if(dataindata == 1 || dataindata == 3) {
			printf("%g\t",inputms[x].theta * (double)1/inputms[x].factor_chrn);
			printf("%g\t",inputms[x].r * (double)1/inputms[x].factor_chrn);
			printf("%g\t",inputms[x].ratio_sv);
		}
		printf("\n");
    }
    printf("\n\n");

    if(file_output) {
        fputs("\n\n Parameters for coalescent simulations:",file_output);
        fputs("\n Parameters related to sequence information data:\n",file_output);
        fprintf(file_output,"\n Number of loci: %d\n",data[0].n_loci);

		fputs("\n Note that the simulations will include the number of sites plus positions with mhits.\n",file_output);

		if((dataindata == 1 || dataindata == 3) && data[0].no_rec_males == 1) fputs("\n Recombination supressed in males\n",file_output);
        if(dataindata == 1 || dataindata == 3) fputs("\nnloci\tfactor_chrn\tnsamples\tnsites\tnmhits\ttheta(4Nu)_locus\trec(4Nr)_locus\ttrans/transv\n",file_output);
		else fputs("\nnloci\tfactor_chrn\tnsamples\tnsites\tnmhits\n",file_output);
        for(x=0;x<data[0].n_loci;x++) {
            if(dataobsequalsim == 1)
                fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
            else
                fprintf(file_output,"locus[%d]\t",x);
            fprintf(file_output,"%g\t",(double)inputms[x].factor_chrn);
            fprintf(file_output,"%d\t",inputms[x].nsam);
            fprintf(file_output,"%ld\t",inputms[x].nsites);
            fprintf(file_output,"%ld\t",inputms[x].nmhits);
            if(dataindata == 1 || dataindata == 3) {
				fprintf(file_output,"%g\t",inputms[x].theta * (double)1/inputms[x].factor_chrn);
				fprintf(file_output,"%g\t",inputms[x].r * (double)1/inputms[x].factor_chrn);
				fprintf(file_output,"%g\t",inputms[x].ratio_sv);
			}
			fprintf(file_output,"\n");
        }
        fputs("\n\n",file_output);
    }
    return;
}

void print_seqinfodataobs_sim(/*int dataindata,struct var *data,struct var2b *inputms,*/ FILE *file_output/*,int dataobsequalsim*/,struct statistics *matrix, struct globvar *global)
{
    int x;
    
    printf("\n Parameters related to observed sequence information data:\n");
    printf("\n Number of loci: %d\n",global[0].n_loci);
    
	printf("\nnloci\tfactor_chrn\tnsamples\tnsites\tnmhits\n");
	for(x=0;x<global[0].n_loci;x++) {
        printf("%d:%s\t",x,matrix[x].gene);
        printf("%g\t",matrix[x].factor_chrn);
        printf("%d\t",matrix[x].nsamples);
        printf("%.2f\t",matrix[x].nsites);
		printf("%d\t",matrix[x].nmhits);
		printf("\n");
    }
    printf("\n\n");

    if(file_output) {
        fputs("\n Parameters related to observed sequence information data:\n",file_output);
        fprintf(file_output,"\n Number of loci: %d\n",global[0].n_loci);
        fputs("\nnloci\tfactor_chrn\tnsamples\tnsites\tnmhits\n",file_output);
        for(x=0;x<global[0].n_loci;x++) {
            fprintf(file_output,"%d:%s\t",x,matrix[x].gene);
            fprintf(file_output,"%g\t",matrix[x].factor_chrn);
            fprintf(file_output,"%d\t",matrix[x].nsamples);
            fprintf(file_output,"%.2f\t",matrix[x].nsites);
            fprintf(file_output,"%d\t",matrix[x].nmhits);
			fprintf(file_output,"\n");
        }
        fputs("\n\n",file_output);
    }
    return;
}

void print_evocond_sim(struct var *data,struct var2b *inputms, FILE *file_output/*,int dataobsequalsim*/,struct statistics *matrix)
{
    int x,y,z;
    
    printf("\n\n Parameters for coalescent simulations:");
    printf("\n Parameters related to evolutionary conditions of simulations:\n");
    printf("\n\nNumber of evolutive processes ocurred in different loci: %d\n\n",data[0].n_events);
	if(data[0].n_events > 1) {
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
	}
    if(data[0].neutrality) {
        if(data[0].n_events==1) printf("\n\n Neutral panmictic model process in Coalescent simulations.");
		else {
			if(data[0].n_events==2)  printf("\n\n Neutral panmictic model process for all loci except in loci indicated below:\n");
			else {
				printf("\n\n Neutral panmictic model process in loci indicated below:\n");
				for(x=0;x<data[0].n_loci;x++) {
					if(inputms[x].neutrality) printf("\n   %d:%s",x,matrix[x].gene);
				}
			}
		}
    }
    if(data[0].changesizepop) {
        printf("\n\n Change in population size. Logistic equation within events:");		
		for(z=1;z<=data[0].changesizepop;z++) {	
			if(data[0].changesizepop > 1) printf("\n\n\tProcess %d\n",z);
			x=0;
			while(inputms[x].changesizepop != z) x++;
			printf("\n\tNumber of events: %d\n",inputms[x].nintn);
			printf("\n\t\tInitial_N\tFinal_N\tDuration(in 4N gen)\n");
			for(y=0;y<inputms[x].nintn;y++) {
				printf("\t\t%g\t",inputms[x].nrec[y]);
				printf("%g\t",inputms[x].npast[y]);
				printf("%g\n",inputms[x].tpast[y]);
			}
			printf("\n\t\tts: %g\n",inputms[x].ts);
			if(data[0].n_events==1) printf("\n\n All loci follow this model.");
			else {
				printf("\n\t\tLoci following these parameters:");
				for(x=0;x<data[0].n_loci;x++) {
					if(inputms[x].changesizepop == z) printf("\n\t\t%d:%s",x,matrix[x].gene);
				}
			}
		}
    }
    if(data[0].subdivision) {
        printf("\n\n Subdivision. Equal migration and extinction rates among demes:");
        printf("\n    The level of variation and recombination indicated is within deme.\n");
        printf("\n    In order to do correct comparisons using the island isotropic migration model,");
        printf("\n    the level of variation and recombination within deme is divided by the number of demes when simulations are performed.\n");
        printf("\n    All samples are in the same deme.");
		for(z=1;z<=data[0].subdivision;z++) {	
			if(data[0].subdivision > 1) printf("\n\tProcess %d\n",z);
			x=0;
			while(inputms[x].subdivision != z) x++;
			printf("\n    Total number of demes: %d",inputms[x].npop);
			printf("\n    Migration among demes (4Nm): %g",inputms[x].mig_rate);
			printf("\n    Population size N (in relation to initial N) of every population:\n\t");
			for(y=0;y<inputms[x].npop;y++) 
				printf("%g\t",inputms[x].factor_pop[y]);
			if(data[0].n_events==1) printf("\n\n All loci follow this model.");
			else {
				printf("\n    Loci following these parameters:");
				for(x=0;x<data[0].n_loci;x++) {
					if(inputms[x].subdivision == z) printf("\n    %d:%s",x,matrix[x].gene);
				}
			}
		}
    }
    if(data[0].split_pop) {
        printf("\n\n Refugia model. Split populattion in several isolated demes and rejoin:");
		for(z=1;z<=data[0].split_pop;z++) {	
			if(data[0].split_pop > 1) printf("\n\n\tProcess %d\n",z);
			x=0;
			while(inputms[x].split_pop != z) x++;
			printf("\n\n    Split population in the past in %d refuges:",inputms[x].npop);
			printf("\n    Duration (in 4N generations) until the split: %g",inputms[x].time_split);
			printf("\n    Duration (in 4N generations) from the split until rejoining in a single population: %g",inputms[x].time_scoal);
			printf("\n    Migration among refuges (4Nm): %g",inputms[x].mig_rate);
			printf("\n    Population size N (in relation to initial N) of every population in the splitted phase:\n\t");
			for(y=0;y<inputms[x].npop;y++) 
				printf("%g\t",inputms[x].factor_pop[y]);
			printf("\n     Population size N (in relation to the initial N) of the ancestral population: %g",
				inputms[x].factor_anc);
			printf("\n     Average contribution of every population (in the splitted phase) to the present population:\n\t");
			for(y=0;y<inputms[x].npop;y++) 
				printf("%g\t",inputms[x].freq[y]);
			if(data[0].n_events==1) printf("\n\n All loci follow this model.");
			else {
				printf("\n     Loci following these parameters:");
				for(x=0;x<data[0].n_loci;x++) {
					if(inputms[x].split_pop == z) printf("\n    %d:%s",x,matrix[x].gene);
				}
			}
		}

    }
    if(data[0].pselection) {
        printf("\n\n Positive selection process");
		for(z=1;z<=data[0].pselection;z++) {	
			if(data[0].pselection > 1) printf("\n\n\tProcess %d\n",z);
			x=0;
			while(inputms[x].pselection != z) x++;
			printf("\n    Population size (in number of individuals N): %ld",inputms[x].pop_size);
			printf("\n    Selection coefficient 4Ns: %g",inputms[x].pop_sel);
			printf("\n    Nucleotide position selected: %ld",inputms[x].sel_nt+1);
			printf("\n    Time (in 4N generations) since the selective process ended: %g\n",inputms[x].sinit);
			if(data[0].n_events==1) printf("\n\n All loci follow this model.");
			else {
				printf("\n     Loci following these parameters:");
				for(y=0;y<data[0].n_loci;y++) {
					if(inputms[y].pselection == z) printf("\n    %d:%s",y,matrix[y].gene);
				}
			}
		}
    }
    printf("\n\n");
    
    if(file_output) {
        fputs("\n\n Parameters for coalescent simulations:",file_output);
        fputs("\n Parameters related to evolutionary conditions of simulations:\n",file_output);
		fprintf(file_output,"\n\nNumber of evolutive processes ocurred in different loci: %d\n",data[0].n_events);
		if(data[0].n_events > 1) {
			fprintf(file_output,"**********************************************************\n");
			fprintf(file_output,"**********************************************************\n");
			fprintf(file_output,"* WARNING: Selective processes affect specific loci.     *\n");
			fprintf(file_output,"* A demographic process should affect the entire genome. *\n");
			fprintf(file_output,"*  The use of different demographic processes            *\n");
			fprintf(file_output,"*  in different loci could not be justified.             *\n");
			fprintf(file_output,"*  A selective process is modeled under                  *\n");
			fprintf(file_output,"*  a stationary population background.                   *\n");
			fprintf(file_output,"**********************************************************\n");
			fprintf(file_output,"**********************************************************\n");
		}
		if(data[0].neutrality) {
			if(data[0].n_events==1) fprintf(file_output,"\n\n Neutral panmictic model process in Coalescent simulations.");
			else {
				if(data[0].n_events==2)  fprintf(file_output,"\n\n Neutral panmictic model process except in loci indicated below.\n");
				else {
					fprintf(file_output,"\n\n Neutral panmictic model process in loci indicated below:\n");
					for(x=0;x<data[0].n_loci;x++) {
						if(inputms[x].neutrality) fprintf(file_output,"\n   %d:%s",x,matrix[x].gene);
					}
				}
			}
		}
		if(data[0].changesizepop) {
			fprintf(file_output,"\n\n Change in population size. Logistic equation within events:");		
			for(z=1;z<=data[0].changesizepop;z++) {	
				if(data[0].changesizepop > 1) fprintf(file_output,"\n\n\tProcess %d\n",z);
				x=0;
				while(inputms[x].changesizepop != z) x++;
				fprintf(file_output,"\n\n\tNumber of events: %d\n",inputms[x].nintn);
				fprintf(file_output,"\n\tInitial_N\tFinal_N\tDuration(in 4N gen)\n");
				for(y=0;y<inputms[x].nintn;y++) {
					fprintf(file_output,"\t%g\t",inputms[x].nrec[y]);
					fprintf(file_output,"%g\t",inputms[x].npast[y]);
					fprintf(file_output,"%g\n",inputms[x].tpast[y]);
				}
				fprintf(file_output,"\n\t\tts: %g\n",inputms[x].ts);
				if(data[0].n_events==1) fprintf(file_output,"\n\n All loci follow this model.");
				else {
					fprintf(file_output,"\n\t\tLoci following these parameters:");
					for(x=0;x<data[0].n_loci;x++) {
						if(inputms[x].changesizepop == z) fprintf(file_output,"\n\t\t%d:%s",x,matrix[x].gene);
					}
				}
			}
		}
		if(data[0].subdivision) {
			fprintf(file_output,"\n\n Subdivision. Equal migration and extinction rates among demes:");
			fprintf(file_output,"\n    The level of variation and recombination indicated is within deme.\n");
			fprintf(file_output,"\n    In order to do correct comparisons using the island isotropic migration model,");
			fprintf(file_output,"\n    the level of variation ad recombination within deme is divided by the number of demes when simulations are performed.\n");
			fprintf(file_output,"\n    All samples are in the same deme.");
			for(z=1;z<=data[0].subdivision;z++) {	
				if(data[0].subdivision > 1) fprintf(file_output,"\n\n\tProcess %d\n",z);
				x=0;
				while(inputms[x].subdivision != z) x++;
				fprintf(file_output,"\n    Total number of demes: %d",inputms[x].npop);
				fprintf(file_output,"\n    Migration among demes (4Nm): %g",inputms[x].mig_rate);
				fprintf(file_output,"\n    Population size N (in relation to initial N) of every population:\n\t");
				for(y=0;y<inputms[x].npop;y++) 
					fprintf(file_output,"%g\t",inputms[x].factor_pop[y]);
				if(data[0].n_events==1) fprintf(file_output,"\n\n All loci follow this model.");
				else {
					fprintf(file_output,"\n\t\tLoci following these parameters:");
					for(x=0;x<data[0].n_loci;x++) {
						if(inputms[x].subdivision == z) fprintf(file_output,"\n\t%d:%s",x,matrix[x].gene);
					}
				}
			}
		}
		if(data[0].split_pop) {
			fprintf(file_output,"\n\n Refugia model. Split populattion in several isolated demes and rejoin:");
			for(z=1;z<=data[0].split_pop;z++) {	
				if(data[0].split_pop > 1) fprintf(file_output,"\n\n\tProcess %d\n",z);
				x=0;
				while(inputms[x].split_pop != z) x++;
				fprintf(file_output,"\n\n Split population in the past in %d refuges:",inputms[x].npop);
				fprintf(file_output,"\n   Duration (in 4N generations) until the split: %g",inputms[x].time_split);
				fprintf(file_output,"\n   Duration (in 4N generations) from the split until rejoining in a single population: %g",inputms[x].time_scoal);
				fprintf(file_output,"\n   Migration among refuges (4Nm): %g",inputms[x].mig_rate);
				fprintf(file_output,"\n   Population size N (in relation to initial N) of every population in the splitted phase:\n\t");
				for(y=0;y<inputms[x].npop;y++) 
					fprintf(file_output,"%g\t",inputms[x].factor_pop[y]);
				fprintf(file_output,"\n     Population size N (in relation to the initial N) of the ancestral population: %g",
					inputms[x].factor_anc);
				fprintf(file_output,"\n    Average contribution of every population (in the splitted phase) to the present population:\n\t");
				for(y=0;y<inputms[x].npop;y++) 
					fprintf(file_output,"%g\t",inputms[x].freq[y]);
				if(data[0].n_events==1) fprintf(file_output,"\n\n All loci follow this model.");
				else {
					fprintf(file_output,"\n\t\tLoci following these parameters:");
					for(x=0;x<data[0].n_loci;x++) {
						if(inputms[x].split_pop == z) fprintf(file_output,"\n\t%d:%s",x,matrix[x].gene);
					}
				}
			}

		}
		if(data[0].pselection) {
			fprintf(file_output,"\n\n Positive selection process");
			for(z=1;z<=data[0].pselection;z++) {	
				if(data[0].pselection > 1) fprintf(file_output,"\n\n\tProcess %d\n",z);
				x=0;
				while(inputms[x].pselection != z) x++;
				fprintf(file_output,"   Population size (in number of individuals N): %ld",inputms[x].pop_size);
				fprintf(file_output,"\n    Selection coefficient 4Ns: %g",inputms[x].pop_sel);
				fprintf(file_output,"\n    Nucleotide position selected: %ld",inputms[x].sel_nt+1);
				fprintf(file_output,"\n    Time (in 4N generations) since the selective process ended: %g\n",inputms[x].sinit);
				if(data[0].n_events==1) fprintf(file_output,"\n\n All loci follow this model.");
				else {
					fprintf(file_output,"\n\t\tLoci following these parameters:");
					for(y=0;y<data[0].n_loci;y++) {
						if(inputms[y].pselection == z) fprintf(file_output,"\n\t%d:%s",y,matrix[y].gene);
					}
				}
			}
		}
        fputs("\n\n",file_output);
    }
	return;
}

