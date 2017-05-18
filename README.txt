MANVa: Multilocus Analysis of Nucleotide Variation package

Sebastian E. Ramos-Onsins* and Thomas Mitchell-Olds**
Previous address: Max-Planck Intitute of Chemical Ecology
Present address: *Centre of research in Agricultural Genomics; ***Duke University
 
This package contains the application MANVa


MANVA:

Current Released Version: MANVa: v0.9892beta (20120515)

The application program MANVa is designed to:

(i) Analyze empirical data from a single population and an outgroup species (if included). It is designed for analyzing a large number of independent loci.
(ii) It is able to read GFF files for each alignment. For example, silent, synonymous and non-synonymous positions can be analyzed separately.
(iii) Perform multilocus coalescent simulations conditioned on the population mutation rate theta (4Nm , where N is the effective population size and m is the mutational rate).
(iv) NOT YET AVAILABLE: Calculate the levels of variation by maximum likelihood for the multilocus set. Probabilities are obtained by coalescent simulation.
(v) Calculate probabilities and confidence intervals for the observed data when comparing with simulated data.
(vi) Histogrames for the observed and simulated data are also obtained.

This program is based on a previous version of Hudson’s coalescent program ms (Hudson, 2002) and modified for the above purposes purposes.

HOW TO INTRODUCE EMPIRICAL DATA:
Aligned sequences for a population and an outgroup should be in nbrf/pir or multiFASTA format.

COMPILE:
For PCwindows, change in MuLoNeTests.h the definition to:

#define DOS_CL 1
#define UNIX_CL 0

Compile with gcc:

gcc *.c -lm -o manva -Wall -pedantic

RUN:
Simply execute the program. It is not necessary an input file.

DOCUMENTATION:
In next versions.
