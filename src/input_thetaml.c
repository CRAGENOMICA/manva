/*
 *  input_thetaml.c
 *  thetamlS
 *
 *  Created by Sebastian E. Ramos-Onsins on April 09 2004.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>

struct parameters_thetaml {
    int nloci;
    long int *numberloci;
    int *nsamA;
    double *factor_chrn;
    int *SA;
    double *length;
    double thetamin_user;
    double thetamax_user;
    int steps_user;
    long int ncoalsim;
    long int seed;
};

int input_datathetasml(struct var **datams, struct var2b **inputms,struct parameters_thetaml **data)
{
    int i;
    
	/*Introduce data from datams and inputms to data (thetamlS)*/
	data[0][0].nloci = datams[0][0].n_loci;
	if(!(data[0][0].numberloci = (long int *)malloc((data[0][0].nloci)*sizeof(long int)))) {
		puts("\nError: Not enough memory.\n");
		return 1;
	}
	if(!(data[0][0].nsamA = (int *)malloc((data[0][0].nloci)*sizeof(int)))) {
		puts("\nError: Not enough memory.\n");
		return 1;
	}
	if(!(data[0][0].factor_chrn = (double *)malloc((data[0][0].nloci)*sizeof(double)))) {
		puts("\nError: Not enough memory.\n");
		return 1;
	}
	if(!(data[0][0].SA = (int *)malloc((data[0][0].nloci)*sizeof(int)))) {
		puts("\nError: Not enough memory.\n");
		return 1;
	}
	if(!(data[0][0].length = (double *)malloc((data[0][0].nloci)*sizeof(double)))) {
		puts("\nError: Not enough memory.\n");
		return 1;
	}
	
	for(i=0;i<data[0][0].nloci;i++) {
		data[0][0].numberloci[i] = (long int)i;
		data[0][0].factor_chrn[i] = inputms[0][i].factor_chrn;
		data[0][0].nsamA[i] = inputms[0][i].nsam;
		data[0][0].SA[i] = inputms[0][i].S;
		data[0][0].length[i] = (double)inputms[0][i].nsites;
	}
	data[0][0].thetamin_user = datams[0][0].thetamin_bp;
	data[0][0].thetamax_user = datams[0][0].thetamax_bp;
	data[0][0].steps_user = datams[0][0].steps;
	data[0][0].seed = datams[0][0].seed2;
	data[0][0].ncoalsim = datams[0][0].n_iter2;
		
    return 0;
}

