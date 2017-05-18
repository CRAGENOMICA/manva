/*
 *  open_helpmenu.c
 *  MuLoNeTests
 *
 *  Created by sebas on Mon Feb 24 2003.
 *
 */

#include "MuLoNeTests.h"
#include <stdio.h>
#include <stdlib.h>

void open_helpmenu(FILE *file_output)
{
    /*
    Introduction
        tests based on the original idea of J. Hey.
        independent loci
        preferred non coding regions
        idea is obtaining tests that use a number of different independent genealogies.
        SNPs or nucleotide sequences are allowed.
    Architecture of the program
        code
        interface
        operative systems
    menu of the program.
        file
        observed data analiysis
        coalescent simulations
        help
    input data.
        DNA sequences or SNPs.
        format.
            FASTA
            NBRF/PIR
        limitations.
    analysis.
        assumptions.
            independence of the loci.
            gaps eliminated from the entire alignment.
            no mhits considered. eliminated from all the set of data (including outgroup).
            outgroup: distance, be careful.
            For S=0, included in HKA, not in the other tests.
            Shared polymorpisms included in HKA, not in the other tests.
        statistics.
            outgroup or not.
            based on seg sites freq.
                waterson
                tajima
                fuli
                faywu
            haplotypes.
                nhapl
            linkage.
                b
                q
                za
            divergence
                shared
                fixed
        neutrality tests.
            one locus.
                outgroup or not.
                calculation or not depending n and S
                based on seg sites freq.
                    tajimad
                    R2
                    fulid
                    fulif
                    faywuh
                        be careful, it is not weighted.
                haplotypes.
                    fsfu (rozas program)
                linkage.
                    bwall
                    qwall
                    zarozas
                divergence
                    hka.
                        estimation of parameters. Secant approach.
                        Contrasted by using X2 and simulations. 
                        Jukes and Cantor correction.
            multilocus.
                weighting problem: ideal: use the same number of samples and same variability for each locus.
                average and variance from calculated loci. other moments.
        non-parametric tests
            wilcoxon's tests to detect differentiation. 
    simulations.
        coalescent simulations.
            parameters.
                Using parameters from observed data.
                introduced by hand.
            recommended number of iterations.
            neutral equilibrium model.
            recombination within and free recombination among.
        display options.
            all iterations.
            confidence intreval and prob. obs data.
            histograms
        advanced options.
            changing pop size. increase/decrease.
            subdivision under island model.
            split populations (refugia).
            selction in one or several loci.
    Histograms.
        observed data.
        simulated data.
        combination observed/simulated.
        Kolmogorov-Smirnov test.
    authors.
        adresses.
        contact.
        reference.
    references.
    */
	if(file_output) fprintf(file_output,"\n");

    return;
}

