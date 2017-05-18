/*
 *  calc_multiloci.c
 *  MuLoNeTests
 *
 *  Created by sebas on Thu Feb 27 2003.
 *
 */

#include "MuLoNeTests.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int calcobsmultiloci(struct statistics *matrix,struct statmulo **matrixml,int *outgroup, int *n_loci/*, FILE *file_output*/)
{
    int x;
    double jc;

    if((x = *n_loci) == 0) {
        matrixml[0][0].nloci = 0;
        matrixml[0][0].Snsites = (double)0;
        matrixml[0][0].Sbiallsites = 0;
        matrixml[0][0].Sbiallsitesn = 0;
        matrixml[0][0].Snmhits = 0;
        matrixml[0][0].Sshared = 0;
        matrixml[0][0].Sfixed = 0;
        matrixml[0][0].Stransitions = 0;
        matrixml[0][0].Stransversions = 0;
        matrixml[0][0].Snhapl = 0;
        matrixml[0][0].Snhaplsam = (double)0;
        matrixml[0][0].Shapldiv = (double)0;
        matrixml[0][0].Sewtest = (double)0;
        matrixml[0][0].Sza = (double)0.;
        matrixml[0][0].Sb = 0;
        matrixml[0][0].Sq = 0;
        matrixml[0][0].SR2 = (double)0.;
        matrixml[0][0].SRvpi = (double)0.;
        matrixml[0][0].SRm = (long int)0.;
        matrixml[0][0].Stheta_wat = (double)0.;
        matrixml[0][0].Stheta_watn = (double)0.;
        matrixml[0][0].Stheta_taj = (double)0.;
        matrixml[0][0].Stheta_tajn = (double)0.;
        matrixml[0][0].Stheta_fw = (double)0.;
        matrixml[0][0].Stheta_fuli = (double)0.;
        matrixml[0][0].Stheta_fulin = (double)0.;
        matrixml[0][0].Stheta_L = (double)0.;
        matrixml[0][0].Sndivergence = (double)0.;
        matrixml[0][0].Stheta_wat_nut = (double)0.;
        matrixml[0][0].Stheta_watn_nut = (double)0.;
        matrixml[0][0].Stheta_taj_nut = (double)0.;
        matrixml[0][0].Stheta_tajn_nut = (double)0.;
        matrixml[0][0].Stheta_fw_nut = (double)0.;
        matrixml[0][0].Stheta_fuli_nut = (double)0.;
        matrixml[0][0].Stheta_fulin_nut = (double)0.;
        matrixml[0][0].Stheta_L_nut = (double)0.;
        matrixml[0][0].Sndivergence = (double)0.;
        matrixml[0][0].Sndivergencejc = (double)0.;
        matrixml[0][0].Sndivergence_nut = (double)0.;
        matrixml[0][0].Sndivergencejc_nut = (double)0.;
        matrixml[0][0].StajimaD = (double)0.;
        matrixml[0][0].SfuliD = (double)0.;
        matrixml[0][0].SfuliF = (double)0.;
        matrixml[0][0].SfuliDn = (double)0.;
        matrixml[0][0].SfuliFn = (double)0.;
        matrixml[0][0].SfuFs = (double)0.;
        matrixml[0][0].SfaywuH = (double)0.;
        matrixml[0][0].SfaywuHo = (double)0.;
        matrixml[0][0].SzengE = (double)0.;
        matrixml[0][0].SrZA = (double)0.;
        matrixml[0][0].SwB = (double)0.;
        matrixml[0][0].SwQ = (double)0.;
    
        matrixml[0][0].S2biallsites = 0;
        matrixml[0][0].S2biallsitesn = 0;
        matrixml[0][0].S2shared = 0;
        matrixml[0][0].S2fixed = 0;
        matrixml[0][0].S2nhapl = 0;
        matrixml[0][0].S2nhaplsam = (double)0;
        matrixml[0][0].S2hapldiv = (double)0;
		matrixml[0][0].S2ewtest = (double)0;
        matrixml[0][0].S2za = (double)0.;
        matrixml[0][0].S2b = 0;
        matrixml[0][0].S2q = 0;
        matrixml[0][0].S2R2 = (double)0.;
        matrixml[0][0].S2Rvpi = (double)0.;
        matrixml[0][0].S2Rm = (long int)0.;
        matrixml[0][0].S2theta_wat = (double)0.;
        matrixml[0][0].S2theta_watn = (double)0.;
        matrixml[0][0].S2theta_taj = (double)0.;
        matrixml[0][0].S2theta_tajn = (double)0.;
        matrixml[0][0].S2theta_fw = (double)0.;
        matrixml[0][0].S2theta_fuli = (double)0.;
        matrixml[0][0].S2theta_fulin = (double)0.;
        matrixml[0][0].S2theta_L = (double)0.;
        matrixml[0][0].S2ndivergence = (double)0.;
        matrixml[0][0].S2ndivergencejc = (double)0.;
        matrixml[0][0].S2theta_wat_nut = (double)0.;
        matrixml[0][0].S2theta_watn_nut = (double)0.;
        matrixml[0][0].S2theta_taj_nut = (double)0.;
        matrixml[0][0].S2theta_tajn_nut = (double)0.;
        matrixml[0][0].S2theta_fw_nut = (double)0.;
        matrixml[0][0].S2theta_fuli_nut = (double)0.;
        matrixml[0][0].S2theta_fulin_nut = (double)0.;
        matrixml[0][0].S2theta_L_nut = (double)0.;
        matrixml[0][0].S2ndivergence_nut = (double)0.;
        matrixml[0][0].S2ndivergencejc_nut = (double)0.;
        matrixml[0][0].S2tajimaD = (double)0.;
        matrixml[0][0].S2fuliD = (double)0.;
        matrixml[0][0].S2fuliF = (double)0.;
        matrixml[0][0].S2fuliDn = (double)0.;
        matrixml[0][0].S2fuliFn = (double)0.;
        matrixml[0][0].S2fuFs = (double)0.;
        matrixml[0][0].S2faywuH = (double)0.;
        matrixml[0][0].S2faywuHo = (double)0.;
        matrixml[0][0].S2zengE = (double)0.;
        matrixml[0][0].S2rZA = (double)0.;
        matrixml[0][0].S2wB = (double)0.;
        matrixml[0][0].S2wQ = (double)0.;
    
        matrixml[0][0].nltajD = 0;
        matrixml[0][0].nlflDn = 0;
        matrixml[0][0].nlflF  = 0;
        matrixml[0][0].nlflFn = 0;
        matrixml[0][0].nlFs = 0;
        matrixml[0][0].nlH = 0;
        matrixml[0][0].nlHo = 0;
        matrixml[0][0].nlE = 0;
        matrixml[0][0].nlZ = 0;
        matrixml[0][0].nlB = 0;
        matrixml[0][0].nlQ = 0;
        matrixml[0][0].nlR2 = 0;
    }       
    
    matrixml[0][0].nloci += 1;
    matrixml[0][0].Snsites += matrix[x].nsites;
    matrixml[0][0].Snmhits += matrix[x].nmhits;
    matrixml[0][0].Sbiallsites += matrix[x].biallsites;
    matrixml[0][0].Sbiallsitesn += matrix[x].biallsitesn;
    matrixml[0][0].Sshared += matrix[x].shared;
    matrixml[0][0].Sfixed += matrix[x].fixed;
    matrixml[0][0].Stransitions += matrix[x].transitions;
    matrixml[0][0].Stransversions += matrix[x].transversions;
    matrixml[0][0].Snhapl += matrix[x].nhapl;
	matrixml[0][0].Snhaplsam += matrix[x].nhaplsam;
	matrixml[0][0].Shapldiv += matrix[x].hapldiv;
	matrixml[0][0].Sewtest += matrix[x].ewtest;
    matrixml[0][0].Sza += matrix[x].za;
    matrixml[0][0].Sb += matrix[x].b;
    matrixml[0][0].Sq += matrix[x].q;
	matrixml[0][0].SRvpi += matrix[x].Rvpi;
	matrixml[0][0].SRm += matrix[x].Rm;
    matrixml[0][0].Stheta_wat += matrix[x].theta_wat*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_watn += matrix[x].theta_watn*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_taj += matrix[x].theta_taj*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_tajn += matrix[x].theta_tajn*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_fw += matrix[x].theta_fw*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_L += matrix[x].theta_L*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_fuli += matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_fulin += matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Sndivergence += matrix[x].ndivergence;
    if(matrix[x].ndivergence/(double)matrix[x].nsites < (double)0.75 && matrixml[0][0].Sndivergencejc != (double)-10000) {
        if(matrixml[0][0].Sndivergencejc != (double) -10000) {
            jc = (double)-0.75 * (double)log((double)1. - (double)4./(double)3. * matrix[x].ndivergence/(double)matrix[x].nsites) * (double)matrix[x].nsites;
            matrixml[0][0].Sndivergencejc += jc;
        }
    }
    else matrixml[0][0].Sndivergencejc = (double)-10000;

    matrixml[0][0].Stheta_wat_nut += matrix[x].theta_wat/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_watn_nut += matrix[x].theta_watn/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_taj_nut += matrix[x].theta_taj/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_tajn_nut += matrix[x].theta_tajn/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_fw_nut += matrix[x].theta_fw/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_L_nut += matrix[x].theta_L/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_fuli_nut += matrix[x].theta_fuli/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Stheta_fulin_nut += matrix[x].theta_fulin/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].Sndivergence_nut += matrix[x].ndivergence/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn);
    if(matrix[x].ndivergence/(double)matrix[x].nsites < (double)0.75 && matrixml[0][0].Sndivergencejc_nut != (double)-10000) {
        if(matrixml[0][0].Sndivergencejc_nut != (double) -10000) {
            jc = (double)-0.75 * (double)log((double)1. - (double)4./(double)3. * matrix[x].ndivergence/(double)matrix[x].nsites);
            matrixml[0][0].Sndivergencejc_nut += jc;
        }
    }
    else matrixml[0][0].Sndivergencejc_nut = (double)-10000;

    matrixml[0][0].S2biallsites += matrix[x].biallsites * (matrix[x].biallsites);
    matrixml[0][0].S2biallsitesn += matrix[x].biallsitesn * (matrix[x].biallsitesn);
    matrixml[0][0].S2shared += matrix[x].shared * matrix[x].shared;
    matrixml[0][0].S2fixed += matrix[x].fixed * matrix[x].fixed;
    matrixml[0][0].S2nhapl += matrix[x].nhapl * matrix[x].nhapl;
	matrixml[0][0].S2nhaplsam += matrix[x].nhaplsam * matrix[x].nhaplsam;
	matrixml[0][0].S2hapldiv += matrix[x].hapldiv * matrix[x].hapldiv;
	matrixml[0][0].S2ewtest += matrix[x].ewtest * matrix[x].ewtest;
    matrixml[0][0].S2za += matrix[x].za * matrix[x].za;
    matrixml[0][0].S2b += matrix[x].b * matrix[x].b;
    matrixml[0][0].S2q += matrix[x].q * matrix[x].q;
	matrixml[0][0].S2Rvpi += matrix[x].Rvpi * matrix[x].Rvpi;
	matrixml[0][0].S2Rm += matrix[x].Rm * matrix[x].Rm;
    matrixml[0][0].S2theta_wat += matrix[x].theta_wat*((double)1/matrix[x].factor_chrn) * matrix[x].theta_wat*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_watn += matrix[x].theta_watn*((double)1/matrix[x].factor_chrn) * matrix[x].theta_watn*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_taj += matrix[x].theta_taj*((double)1/matrix[x].factor_chrn) * matrix[x].theta_taj*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_tajn += matrix[x].theta_tajn*((double)1/matrix[x].factor_chrn) * matrix[x].theta_tajn*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_fw += matrix[x].theta_fw*((double)1/matrix[x].factor_chrn) * matrix[x].theta_fw*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_L += matrix[x].theta_L*((double)1/matrix[x].factor_chrn) * matrix[x].theta_L*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_fuli += matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn) * matrix[x].theta_fuli*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_fulin += matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn) * matrix[x].theta_fulin*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2ndivergence += matrix[x].ndivergence * matrix[x].ndivergence;
    if(matrix[x].ndivergence/(double)matrix[x].nsites < (double)0.75 && matrixml[0][0].S2ndivergencejc != (double)-10000) {
        if(matrixml[0][0].S2ndivergencejc != (double) -10000) {
            jc = (double)-0.75 * (double)log((double)1. - (double)4./(double)3. * matrix[x].ndivergence/(double)matrix[x].nsites) * (double)matrix[x].nsites;
            matrixml[0][0].S2ndivergencejc += jc * jc;
            matrixml[0][0].nldivjc += 1;
        }
    }
    else matrixml[0][0].S2ndivergencejc = (double)-10000;
    matrixml[0][0].S2theta_wat_nut += matrix[x].theta_wat/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn) * matrix[x].theta_wat/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_watn_nut += matrix[x].theta_watn/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) * matrix[x].theta_watn/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_taj_nut += matrix[x].theta_taj/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn) * matrix[x].theta_taj/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_tajn_nut += matrix[x].theta_tajn/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) * matrix[x].theta_tajn/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_fw_nut += matrix[x].theta_fw/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn) * matrix[x].theta_fw/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_L_nut += matrix[x].theta_L/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn) * matrix[x].theta_L/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_fuli_nut += matrix[x].theta_fuli/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn) * matrix[x].theta_fuli/(double)(matrix[x].nsites - matrix[x].shared)*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2theta_fulin_nut += matrix[x].theta_fulin/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn) * matrix[x].theta_fulin/(double)matrix[x].nsites*((double)1/matrix[x].factor_chrn);
    matrixml[0][0].S2ndivergence_nut += matrix[x].ndivergence/(double)matrix[x].nsites * matrix[x].ndivergence/(double)matrix[x].nsites;
    if(matrix[x].ndivergence/(double)matrix[x].nsites < (double)0.75 && matrixml[0][0].S2ndivergencejc_nut != (double)-10000) {
        if(matrixml[0][0].S2ndivergencejc_nut != (double) -10000) {
            jc = (double)-0.75 * (double)log((double)1. - (double)4./(double)3. * matrix[x].ndivergence/(double)matrix[x].nsites);
            matrixml[0][0].S2ndivergencejc_nut += jc * jc;
        }
    }
    else matrixml[0][0].S2ndivergencejc_nut = (double)-10000;

    
    if(*outgroup) {
        if(matrix[x].fuliD >-9999) {
            matrixml[0][0].SfuliD += matrix[x].fuliD;
            matrixml[0][0].S2fuliD += matrix[x].fuliD * matrix[x].fuliD;
            matrixml[0][0].nlflD += 1;
        }
        if(matrix[x].faywuH != (double) -10000) {
            matrixml[0][0].SfaywuH += matrix[x].faywuH;
            matrixml[0][0].S2faywuH += matrix[x].faywuH * matrix[x].faywuH;
            matrixml[0][0].nlH += 1;
        }
        if(matrix[x].faywuHo != (double) -10000) {
            matrixml[0][0].SfaywuHo += matrix[x].faywuHo;
            matrixml[0][0].S2faywuHo += matrix[x].faywuHo * matrix[x].faywuHo;
            matrixml[0][0].nlHo += 1;
        }
        if(matrix[x].zengE != (double) -10000) {
            matrixml[0][0].SzengE += matrix[x].zengE;
            matrixml[0][0].S2zengE += matrix[x].zengE * matrix[x].zengE;
            matrixml[0][0].nlE += 1;
        }
        if(matrix[x].fuliF >-9999) {
            matrixml[0][0].SfuliF += matrix[x].fuliF;
            matrixml[0][0].S2fuliF += matrix[x].fuliF * matrix[x].fuliF;
            matrixml[0][0].nlflF += 1;
        }
    }
    if(matrix[x].tajimaD >-9999) {
        matrixml[0][0].StajimaD += matrix[x].tajimaD;
        matrixml[0][0].S2tajimaD += matrix[x].tajimaD * matrix[x].tajimaD;
        matrixml[0][0].nltajD += 1;
    }
    if(matrix[x].R2 >-9999) {
        matrixml[0][0].SR2 += matrix[x].R2;
        matrixml[0][0].S2R2 += matrix[x].R2 * matrix[x].R2;
        matrixml[0][0].nlR2 += 1;
    }
    if(matrix[x].fuliDn >-9999) {
        matrixml[0][0].SfuliDn += matrix[x].fuliDn;
        matrixml[0][0].S2fuliDn += matrix[x].fuliDn * matrix[x].fuliDn;
        matrixml[0][0].nlflDn += 1;
    }
    if(matrix[x].fuliFn >-9999) {
        matrixml[0][0].SfuliFn += matrix[x].fuliFn;
        matrixml[0][0].S2fuliFn += matrix[x].fuliFn * matrix[x].fuliFn;
        matrixml[0][0].nlflFn += 1;
    }
    if(matrix[x].fuFs != (double) -10000) {
        matrixml[0][0].SfuFs += matrix[x].fuFs;
        matrixml[0][0].S2fuFs += matrix[x].fuFs * matrix[x].fuFs;
        matrixml[0][0].nlFs += 1;
    }
    if(matrix[x].rZA >-9999) {
        matrixml[0][0].SrZA += matrix[x].rZA;
        matrixml[0][0].S2rZA += matrix[x].rZA * matrix[x].rZA;
        matrixml[0][0].nlZ += 1;
    }
    if(matrix[x].wB >-9999) {
        matrixml[0][0].SwB += matrix[x].wB;
        matrixml[0][0].S2wB += matrix[x].wB * matrix[x].wB;
        matrixml[0][0].nlB += 1;
    }
    if(matrix[x].wQ >-9999) {
        matrixml[0][0].SwQ += matrix[x].wQ;
        matrixml[0][0].S2wQ += matrix[x].wQ * matrix[x].wQ;
        matrixml[0][0].nlQ += 1;
    }
    
    return 1;
}
