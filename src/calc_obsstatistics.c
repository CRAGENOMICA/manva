/*
 *  calc_obsstatistics.c
 *  MuLoNeTests
 *  BIOMETRY FUNCTIONS INSIDE
 *  Created by sebas on Tue Feb 25 2003.
 *
 */

#include "MuLoNeTests.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int calc_statistics(struct statistics **matrix,int *outgroup,int *n_loci/*,FILE *file_output*/,int n_samp, double algsites,
	char *matrix_pol,int *matrix_freq,int nnnsam, int nnnoutg, int bial_sites,char **names,char *name_outgroup,int lian_in,struct lianinput **liandata,double factor_chrn,int mhits)
{
    int x,a,b;
    int chapl(int,char *,int *,int,int, int,char **,char *,int,struct lianinput **,int *);
	double testHap(int, int *);
    int clinkage(struct statistics **,char *,int *,int *,/*int *,*/int,int,int,int);
    int cthetas(struct statistics **,char *,int *,int *,int *,int,int,int,int);
	int recombinationpi(struct statistics **,char *, int *,int * /*, int * */,int, int,int,int);
	int Min_rec_obs(int,struct statistics **,char *,int *,int, int,int,int);
	double EWtest(int, int *);

    long int n_comb,dif;
	int *fhapl;
	
	if((fhapl = (int *)malloc(nnnsam*sizeof(int))) == 0) {
		puts("calloc error fhapl.1");
		return -1;
    }
    
    matrix[0][*n_loci].nsamples = nnnsam;
    matrix[0][*n_loci].noutgroups = nnnoutg;
    matrix[0][*n_loci].nsites = algsites;
    matrix[0][*n_loci].factor_chrn = factor_chrn;
    matrix[0][*n_loci].nmhits = mhits;
    matrix[0][*n_loci].biallsites = 0;
    matrix[0][*n_loci].biallsitesn = 0;
    matrix[0][*n_loci].shared = 0;
    matrix[0][*n_loci].fixed = 0;
    
    matrix[0][*n_loci].nhapl = 0;
	matrix[0][*n_loci].nhaplsam = (double)0;
	matrix[0][*n_loci].hapldiv = (double)0;
	matrix[0][*n_loci].ewtest = (double)0;
	matrix[0][*n_loci].Rvpi = (double)0;
	matrix[0][*n_loci].Rm = 0;
    matrix[0][*n_loci].za = (double)0.;
    matrix[0][*n_loci].b = (int)0;
    matrix[0][*n_loci].q = (int)0;
    matrix[0][*n_loci].theta_wat = (double)0.;
    matrix[0][*n_loci].theta_watn = (double)0.;
    matrix[0][*n_loci].theta_taj = (double)0.;
    matrix[0][*n_loci].theta_tajn = (double)0.;
    matrix[0][*n_loci].theta_fw = (double)0.;
    matrix[0][*n_loci].theta_fuli = (double)0.;
    matrix[0][*n_loci].theta_fulin = (double)0.;
    matrix[0][*n_loci].theta_L = (double)0.;
    matrix[0][*n_loci].ndivergence = (double)0.;
    
    for(x=0;x<bial_sites;x++) {
        if(*(matrix_freq+x) > 0 && *(matrix_freq+x) < nnnsam) matrix[0][*n_loci].biallsitesn += 1;
    	if(*(matrix_freq+x) >= 0 && *(matrix_freq+x) < nnnsam) matrix[0][*n_loci].biallsites += 1;
        if(*(matrix_freq+x) == 0) matrix[0][*n_loci].shared += 1;
        if(*(matrix_freq+x) == nnnsam) matrix[0][*n_loci].fixed += 1;
    }
    if((matrix[0][*n_loci].nhapl = chapl(n_samp,matrix_pol,matrix_freq,nnnoutg,nnnsam,bial_sites,names,name_outgroup,lian_in,liandata,fhapl)) < 0) return 0;
    matrix[0][*n_loci].nhaplsam = (double)matrix[0][*n_loci].nhapl/(double)matrix[0][*n_loci].nsamples;
	matrix[0][*n_loci].hapldiv = testHap(matrix[0][*n_loci].nsamples,fhapl);
	matrix[0][*n_loci].ewtest = EWtest(matrix[0][*n_loci].nsamples,fhapl);
	if(clinkage(matrix,matrix_pol,matrix_freq,n_loci/*,outgroup*/,n_samp,nnnsam,nnnoutg,bial_sites) < 0) 
		return 0;
    if(cthetas(matrix,matrix_pol,matrix_freq,n_loci,outgroup,n_samp,nnnsam,nnnoutg,bial_sites) < 0) 
		return 0;
    if(*outgroup) {
        n_comb = (long)nnnsam * (long)nnnoutg;
        for(x=0;x<bial_sites;x++) {
            if((dif = *(matrix_freq+x)) > 0) { /*polymorphic or fixed sample*/
                matrix[0][*n_loci].ndivergence += ((double)dif * (double)nnnoutg)/(double)n_comb;
            }
            else if(dif < 0) {/*polymorphic outgroup*/
                    matrix[0][*n_loci].ndivergence += ((double)-dif * (double)nnnsam)/(double)n_comb;
                }
                else {/*shared polymorphism*/
                    for(a=0;a<nnnoutg;a++)
                        for(b=nnnoutg;b<n_samp;b++)
                            if(*(matrix_pol+((a*(n_samp))+x)) != *(matrix_pol+((b*(n_samp))+x)))
                                matrix[0][*n_loci].ndivergence += (double)1./(double)n_comb;
                }
        }
    }
    else matrix[0][*n_loci].ndivergence = -1;
    if(recombinationpi(matrix,matrix_pol,matrix_freq,n_loci/*,outgroup*/,n_samp,nnnsam,nnnoutg,bial_sites) < 0) 
		return 0;
	matrix[0][*n_loci].Rm = Min_rec_obs(0,matrix,matrix_pol,matrix_freq,n_samp,nnnsam,nnnoutg,bial_sites);

	free(fhapl);
    return 1;
}

int chapl(int n_samp,char *matrix_pol,int *matrix_freq,int nnnoutg,int nnnsam,int bial_sites,char **names,char *name_outgroup,int lian_in,struct lianinput **liandata, int *fhapl)
{
    /*calcular n_hapl*/
    int *haplotype = 0;
    char *hapl = 0;
    int nhapl,nhapl2;
    int a,b,y0;
	
	/*for lian*/
	int aa,bb,c;
	char *strin;
	
	/*main*/
    if((haplotype = (int *)malloc(nnnsam*sizeof(int))) == 0) {
        puts("calloc error chapl.1");
        return -1;
    }
    if((hapl = (char *)calloc(nnnsam*bial_sites,sizeof(char))) == 0) {
        puts("calloc error chapl.2");
        return -1;
    }
	for(a=0;a<nnnsam;a++) fhapl[a] = (int)0;

    for(a=0;a<bial_sites;a++)
        if(*(matrix_freq+a) > 0) /*avoid shared. Only polymorphisms in the sampled species*/
            for(y0=nnnoutg;y0<n_samp;y0++)
                hapl[(y0-nnnoutg)*bial_sites+a] = *(matrix_pol+((a*(n_samp))+y0));                    

	/*calculating freq of haplotypes and number of haplotypes from hapl (goes froom 0 to nnnsam)*/
	nhapl = 0;                
	for(a=0;a<nnnsam;a++) haplotype[a] = 1;            
	for(a=0;a<nnnsam-1;a++) {
		if(haplotype[a]) {
			nhapl += 1;
			for(b=a+1;b<nnnsam;b++) {
				if(haplotype[b]) {
					if(memcmp(hapl+a*bial_sites,hapl+b*bial_sites,bial_sites) == 0) { 
						haplotype[a] += 1;
						haplotype[b] = 0;
					}
				}
			}
		}
	}
	if(haplotype[a]) nhapl += 1;	
	for(a=0;a<nnnsam;a++) fhapl[a] = haplotype[a]; /*calcular freq. haplotips*/
	
	
	/*include here the lian struct data*/
	if(lian_in == 1) {
		/*Giving a number to each haplotype*/
		nhapl2 = 0;                
		for(a=0;a<nnnsam;a++) haplotype[a] = 0;
		for(a=0;a<nnnsam-1;a++) {
			if(haplotype[a]==0) {
				nhapl2 += 1;
				haplotype[a] = nhapl2;
				for(b=a+1;b<nnnsam;b++) {
					if(haplotype[b]==0) {
						if(memcmp(hapl+a*bial_sites,hapl+b*bial_sites,bial_sites) == 0) { 
							haplotype[a] = nhapl2;
							haplotype[b] = nhapl2;
						}
					}
				}
			}
		}
		/*if liandata is new (first locus)*/
		if(liandata[0][0].nloci == 0) {
			liandata[0][0].nsam = nnnsam;
			/*realloc liandata for a new sample*/
			if((liandata[0][0].namesam = (char **)realloc(liandata[0][0].namesam,(liandata[0][0].nsam)*sizeof(char *))) == 0) {
				puts(" Matrix not reallocated! NOT SUCCESSFUL.");
				return -10000;
			}
			if((liandata[0][0].samhap = (int **)realloc(liandata[0][0].samhap,(liandata[0][0].nsam)*sizeof(int *))) == 0) {
				puts(" Matrix not reallocated! NOT SUCCESSFUL.");
				return -10000;
			}
			for(a=0,c=0;a<nnnsam;c++) {
				if(nnnoutg > 0) {
					if((strin = strstr(names[c],name_outgroup)) == 0) {
						if((liandata[0][0].namesam[a] = (char *)calloc(128,sizeof(char ))) == 0) {
							puts(" Matrix not reallocated! NOT SUCCESSFUL.");
							return -10000;
						}
						if((liandata[0][0].samhap[a] = (int *)calloc(1,sizeof(int ))) == 0) {
							puts(" Matrix not reallocated! NOT SUCCESSFUL.");
							return -10000;
						}
						for(b=0;b<128;b++)
							liandata[0][0].namesam[a][b] = names[c][b];
						liandata[0][0].namesam[a][127] ='\0';					
						liandata[0][0].samhap[a][liandata[0][0].nloci] = haplotype[a];/*matrix_pol[][]->hapl[]->haplotype[] sequential*/
						a++;
					}
				}
				else {
					if((liandata[0][0].namesam[a] = (char *)calloc(128,sizeof(char ))) == 0) {
						puts(" Matrix not reallocated! NOT SUCCESSFUL.");
						return -10000;
					}
					if((liandata[0][0].samhap[a] = (int *)calloc(1,sizeof(int ))) == 0) {
						puts(" Matrix not reallocated! NOT SUCCESSFUL.");
						return -10000;
					}
					for(b=0;b<128;b++)
						liandata[0][0].namesam[a][b] = names[c][b];
					liandata[0][0].namesam[a][127] ='\0';					
					liandata[0][0].samhap[a][liandata[0][0].nloci] = haplotype[a];/*matrix_pol[][]->hapl[]->haplotype[] sequential*/
					a++;
				}
			}
		}		
		else {
			/*realloc liandata for each new loci*/
			for(a=0;a<liandata[0][0].nsam;a++) {
				if((liandata[0][0].samhap[a] = (int *)realloc(liandata[0][0].samhap[a],(liandata[0][0].nloci+1) * sizeof(int ))) == 0) {
					puts(" Matrix not reallocated! NOT SUCCESSFUL.");
					return -10000;
				}
				liandata[0][0].samhap[a][liandata[0][0].nloci] = -1;
			}
			/*take only those samples shared with other loci*/
			for(a=0;a<liandata[0][0].nsam;a++) {
				for(b=0,c=0;c<n_samp;c++) {
					if(nnnoutg > 0) {
						if((strin = strstr(names[c],name_outgroup)) == 0) {
							if((strin = strstr(liandata[0][0].namesam[a],names[c])) != 0) {
								liandata[0][0].samhap[a][liandata[0][0].nloci] = haplotype[b];/*matrix_pol[][]->hapl[]->haplotype[] sequential*/
								break;
							}
							b++;
						}
					}
					else {
						if((strin = strstr(liandata[0][0].namesam[a],names[c])) != 0) {
							liandata[0][0].samhap[a][liandata[0][0].nloci] = haplotype[b];/*matrix_pol[][]->hapl[]->haplotype[] sequential*/
							break;
						}
						b++;
					}
				}
			}
			/*eliminate the samples that are not shared*/
			for(a=0;a<liandata[0][0].nsam;a++) {
				if(liandata[0][0].samhap[a][liandata[0][0].nloci] == -1) {
					for(b=a;b<liandata[0][0].nsam-1;b++) {
						for(aa=0;aa<=127;aa++) 
							liandata[0][0].namesam[b][aa] = liandata[0][0].namesam[b+1][aa];
						for(bb=0;bb<=liandata[0][0].nloci;bb++) 
							liandata[0][0].samhap[b][bb] = liandata[0][0].samhap[b+1][bb];
					}
					free(liandata[0][0].namesam[liandata[0][0].nsam-1]);
					free(liandata[0][0].samhap[liandata[0][0].nsam-1]);
					if((liandata[0][0].namesam = (char **)realloc(liandata[0][0].namesam,(liandata[0][0].nsam-1)*sizeof(char *))) == 0) {
						puts(" Matrix not reallocated! NOT SUCCESSFUL.");
						return -10000;
					}
					if((liandata[0][0].samhap = (int **)realloc(liandata[0][0].samhap,(liandata[0][0].nsam-1)*sizeof(int *))) == 0) {
						puts(" Matrix not reallocated! NOT SUCCESSFUL.");
						return -10000;
					}
					liandata[0][0].nsam -= 1;
					a--;
				}
			}
		}
		liandata[0][0].nloci += 1;
	}
	
	free(haplotype);
    free(hapl);
    
    return nhapl;
}

int clinkage(struct statistics **matrix,char *matrix_pol, int *matrix_freq,int *n_loci/*,int *outgroup*/,int n_samp, int nnnsam,int nnnoutg,int bial_sites)
{
    int a,b,c,d,j,k,i;
    int A = 0;
	int B = 0;
    int val10,val20,val21;
    
    double za;
    int za1,za2,valza1,valza2;
    double val11,AA,BB,CC;
	int **veca;
	/*long int comb;*/
    
    if((veca = (int **)calloc(bial_sites,sizeof(int *))) == 0) {
        puts("calloc error veca.1");
        return -1;
    }
	for(d=0;d<bial_sites;d++) {
		if((veca[d] = (int *)calloc(n_samp,sizeof(int))) == 0) {
			puts("calloc error veca.2");
			return -1;
		}
	}
	
    a = 0;
    za = 0;
    for(j=0;j<bial_sites-1;) {
        k = j;
        while(k+1 < bial_sites) { /*calculate k*//*not shared counted...*/
            if(*(matrix_freq+k) >/*=*/ 0 && *(matrix_freq+k) < nnnsam) break; /*fixed and pol. outgroup not counted.*/
            else k++;
        }
        j = k+1;
        while(j < bial_sites) { /*calculate j*//*not shared counted...*/
            if(*(matrix_freq+j) >/*=*/ 0 && *(matrix_freq+j) < nnnsam) break;/*fixed and pol. outgroup not counted*/
            else j++;
        }
        if(j < bial_sites) {
            /*b and q*/
            val20 = val21 = -1;
            b = 0;
            /*za*/
            if(*(matrix_freq+k) > 0 && *(matrix_freq+j) > 0) { /*shared not allowed in za*/
                za1 = za2 = 48 + 1;
                valza1 = *(matrix_freq+k);
                valza2 = *(matrix_freq+j);
            }

            for(i=nnnoutg;i<n_samp;i++) {
                /*b and q*/
                val10 = (*(matrix_pol+((k*(n_samp))+i)) - 48)*4 + (*(matrix_pol+((j*(n_samp))+i)) - 48);
                if(val20 == -1) val20 = val10;
                else{
                    if(val20 != val10) {
                        if(val21 == -1) val21 = val10;
                        else  if(val21 != val10) {
                            b = 1;/*incongruent*/
                            break;
                        }
                    }
                }
            }
            /* za*/
            if(*(matrix_freq+k) > 0 && *(matrix_freq+j) > 0) { /*shared not allowed in za*/
                val11 = 0;
                for(i=nnnoutg;i<n_samp;i++)
                    if((*(matrix_pol+((j*(n_samp))+i))) == za1 && (*(matrix_pol+((k*(n_samp))+i))) == za2) val11++;
                /* R2 = (p1q1 - (p1*q1))"2 / (p1*(1-p1) * (q1*(1-q1))) */
                AA = (double)valza1/(double)nnnsam;
                BB = (double)valza2/(double)nnnsam;
                CC = (double)val11/(double)nnnsam;
                /*if(vala != ((double)nsam-1.0) && valb != ((double)nsam-1.0)) if only informative */
                za += ((CC - (AA*BB)) * (CC - (AA*BB))) / (AA*((double)1.0 - AA)*BB*((double)1.0 - BB));
            }
            /* b and q*/
            if(!b) {/*congruent*/
                B += 1;
                if(!A) {
					for(i=nnnoutg;i<n_samp;i++)
						veca[A][i] = *(matrix_pol+((j*(n_samp))+i));
					A += 1;
				}
				else {
					a = 0;
					for(c=0;c<A;c++) {
						d = 0;
						for(i=nnnoutg;i<n_samp;i++) {
							if(veca[c][i] == *(matrix_pol+((j*(n_samp))+i))) {
								d += 1;
							}
						}
						if(d == n_samp - nnnoutg || d == 0) {
							a = 1;
							break;
						}
					}
					if(!a) {
						for(i=nnnoutg;i<n_samp;i++)
							veca[A][i] = *(matrix_pol+((j*(n_samp))+i));
						A += 1;
					}
				}
            }
        }
    }

    matrix[0][*n_loci].q = B+A;
    matrix[0][*n_loci].b = B;
    matrix[0][*n_loci].za = za;
	
	for(d=0;d<bial_sites;d++) free(veca[d]);
	free(veca);
	
	/*ZnS*//*
    za = 0;
	comb = 0;
    for(j=0,k=0;k<bial_sites-1;) {
        while(k+1 < bial_sites) {
            if(*(matrix_freq+k) > 0 && *(matrix_freq+k) < nnnsam) break;
            else k++;
        }
        j = k+1;
		while(j < bial_sites) {
			while(j < bial_sites) {
				if(*(matrix_freq+j) > 0 && *(matrix_freq+j) < nnnsam) break;
				else j++;
			}
			if(j < bial_sites) {
				if(*(matrix_freq+k) > 0 && *(matrix_freq+j) > 0) {
					za1 = za2 = 48 + 1;
					valza1 = *(matrix_freq+k);
					valza2 = *(matrix_freq+j);
					val11 = 0;
					for(i=nnnoutg;i<n_samp;i++)
						if((*(matrix_pol+((j*(n_samp))+i))) == za1 && (*(matrix_pol+((k*(n_samp))+i))) == za2) val11++;
					AA = (double)valza1/(double)nnnsam;
					BB = (double)valza2/(double)nnnsam;
					CC = (double)val11/(double)nnnsam;
					//if(valza1 != ((double)nsam-1.0) && valza2 != ((double)nsam-1.0)) if only informative
					za += ((CC - (AA*BB)) * (CC - (AA*BB))) / (AA*((double)1.0 - AA)*BB*((double)1.0 - BB));
					comb++;
				}
			}
			j++;
		}
		k++;
    }
	if(comb) matrix[0][*n_loci].za = za/(double)comb;
	else matrix[0][*n_loci].za = -10000;
	*//**/

    return 1;
}

/*Haplotype tests from Depaulis et al.*/
double testHap(int Nsample, int *Freqhap)
{/*Depaulis statistics*/
    int i;
    double H;
    
	if(Nsample < 2) return -10000.;
	
    H = 0;
    for(i=0;i<Nsample;i++)
        H += Freqhap[i] * Freqhap[i];
    H = (double)1 - H/((double)Nsample*(double)Nsample);
    /*and weigthed: haplotype diversity*/
    H = H*(double)Nsample/(double)(Nsample-1);
    
    return H;
}

/*Ewens-Watterson test*/
double EWtest(int Nsample, int *Freqhap)
{
    int i;
    double H;
    
	if(Nsample < 2) return -10000.;

    H = 0;
    for(i=0;i<Nsample;i++)
        H += Freqhap[i] * Freqhap[i];
    H = (double)H/((double)Nsample*(double)Nsample);
    
    return H;
}

int cthetas(struct statistics **matrix,char *matrix_pol,int *matrix_freq,int *n_loci,int *outgroup,int n_samp, int nnnsam,int nnnoutg,int bial_sites)
{
    int x,dif,a,b;
    double an;
    double invcomb = (double)2.0/((double)nnnsam*(nnnsam-(double)1.));
		
    an = (double)0.;
    for(x=1;x<nnnsam;x++) an += (double)1./(double)x;
    
    matrix[0][*n_loci].theta_wat  = (double)matrix[0][*n_loci].biallsites/an;
    matrix[0][*n_loci].theta_watn = (double)matrix[0][*n_loci].biallsitesn/an;
    
    matrix[0][*n_loci].theta_taj   = (double)0.;
    matrix[0][*n_loci].theta_tajn  = (double)0.;
    matrix[0][*n_loci].theta_fuli  = (double)0.;
    matrix[0][*n_loci].theta_fulin = (double)0.;
    matrix[0][*n_loci].theta_fw = (double)0.;
    matrix[0][*n_loci].theta_L = (double)0.;

    for(x=0;x<bial_sites;x++) {
        if(*(matrix_freq+x) > 0 && *(matrix_freq+x) < nnnsam) {
            dif = *(matrix_freq+x);
            matrix[0][*n_loci].theta_tajn += (double)dif*(nnnsam-dif);
            matrix[0][*n_loci].theta_taj += (double)dif*(nnnsam-dif);
            if(*(matrix_freq+x) == 1 || *(matrix_freq+x) == nnnsam-1) matrix[0][*n_loci].theta_fulin += 1;
            if(*outgroup) {
                matrix[0][*n_loci].theta_fw += dif*dif;
                matrix[0][*n_loci].theta_L += dif;
                if(*(matrix_freq+x) == 1) matrix[0][*n_loci].theta_fuli += 1;
            }
        }
    	if(*(matrix_freq+x) == 0) {/*shared*/
            b = *(matrix_pol+((x*(n_samp))+nnnsam-1));
            dif  = 0;
            for(a=nnnoutg;a<nnnsam-1;a++) {
                /*for(b=a;b<nnnsam;b++) {*/
                /*    if(*(matrix_pol+((x*(n_samp))+a)) != *(matrix_pol+((x*(n_samp))+b)))*/
                /*        matrix[0][*n_loci].theta_tajn += 1;*/
                /*}*/
                if(*(matrix_pol+((x*(n_samp))+a)) != b) dif += 1;
            }
            matrix[0][*n_loci].theta_tajn += (double)dif*(nnnsam-dif);
        }
    }
    matrix[0][*n_loci].theta_taj  *= invcomb;
    matrix[0][*n_loci].theta_tajn *= invcomb;
    matrix[0][*n_loci].theta_fulin *= (double)(nnnsam - 1)/(double)nnnsam;
    if(*outgroup) {
		matrix[0][*n_loci].theta_fw *= invcomb;
		matrix[0][*n_loci].theta_L *= (double)1/((double)(nnnsam - 1));		
	}
    return 1;
}

int recombinationpi(struct statistics **matrix,char *matrix_pol, int *matrix_freq,int *n_loci/*,int *outgroup*/,int n_samp, int nnnsam,int nnnoutg,int bial_sites)
{
	/*(from J. rozas code) From HUDSON. Genet. Res. 1987 50:245-250*/
	double km,Sk2,hij,Shij,Shj2,thetam,gcn;
	int i,j,x,dif,kij;
	double pi;
	double C;
	double x1,x2,hx2,xacc;
	int zbracC(double *,double *,double,int);
	double zriddrC(double,double,double,double,int);
	
	if(matrix[0][*n_loci].biallsitesn == 0) {/*option: I decide to give 0 when no bial_sites are*/
		matrix[0][*n_loci].Rvpi = (double)0;
		return 1;
	}
	
	/*calculating km, Shij and Shj2*/
	pi = hij = Shij = Shj2 = (double)0;
    for(x=0;x<bial_sites;x++) {
        if(*(matrix_freq+x) > 0 && *(matrix_freq+x) < nnnsam) {
            dif = *(matrix_freq+x);
            pi += (double)dif*((double)nnnsam-(double)dif);
			hij = (double)1 - (dif*dif + ((double)nnnsam - dif)*((double)nnnsam - dif))/((double)nnnsam*(double)nnnsam);
			Shij += hij;
			Shj2 += (hij * hij);
		}
        if(*(matrix_freq+x) == 0) {/*shared*/
			dif = 0;
			for(i=nnnoutg;i<n_samp-1;i++) 
				if(*(matrix_pol+((x*(n_samp))+i)) != *(matrix_pol+((x*(n_samp))+n_samp-1))) dif++;
            pi += (double)dif*((double)nnnsam-(double)dif);
			hij = (double)1 - (dif*dif + ((double)nnnsam - dif)*((double)nnnsam - dif))/((double)nnnsam*(double)nnnsam);
			Shij += hij;
			Shj2 += (hij * hij);
		}
	}
	km = (double)2*pi/((double)nnnsam*(double)nnnsam);/*note that is not theta_taj*/
	
	/*calculating Sk2*/
	Sk2 = (double)0;
	for(i=nnnoutg;i<n_samp-1;i++) {
		for(j=i+1;j<n_samp;j++) {
			kij = 0;
			for(x=0;x<bial_sites;x++) 
				if(*(matrix_pol+((x*(n_samp))+i)) != *(matrix_pol+((x*(n_samp))+j))) kij += 1;
			Sk2 += (((double)kij - km) * ((double)kij - km));
		}
	}
	/*include the diagonal and the lower matrix*/
	Sk2 *= ((double)2);
	Sk2 += ((double)nnnsam * (km * km));
	Sk2 /= ((double)nnnsam * (double)nnnsam);
	
	/*calculating theta*/
	thetam = Shij * ((double)nnnsam / ((double)nnnsam - (double)1));
	
	/*calculating gCn*/
	gcn = (Sk2 - Shij + Shj2)/(thetam*thetam);
	
	/*Bisection method to find the solution*/
	/*calculate the range*/
	x1 = (double)0.00001;
	x2 = hx2 = (double)matrix[0][*n_loci].biallsitesn;
	
	/*if zbracC = 0, that means that funtionC is not crossing 0*/
	/*
	 if(zbracC(&x1,&x2,gcn,nnnsam) == 0) {
		C = (double)-1.11e30;
		return -1;
	}
	*/
	/**/
	while(zbracC(&x1,&x2,gcn,nnnsam) == 0) {
		if(x2 > 1e10) {
			matrix[0][*n_loci].Rvpi = 1e10;
			return 1;
		}
		x2 = hx2 * (double)10;
	}
	/**/
	
	/*estimate the value of C*/
	xacc = (double)1e-6*x2; /* accuracy of five numbers*/
	C = (double)zriddrC(x1,x2,xacc,gcn,nnnsam);
	if(C == -1.11e30) return -1; /*na*/
	
	matrix[0][*n_loci].Rvpi = C;
	return 1;
}

double functiongCn0(double C,double gCn,int nsam)
{
	double I1,I2,fC,r97,P1,P2,P3,P4;
	
	if(C==(double)0) return -10000;
	r97 = (double)sqrt((double)97);
	I2 = (double)1/r97 * (double)log((((double)2*C + (double)13 - r97) * ((double)13 + r97)) / (((double)2*C + (double)13 + r97) * ((double)13 - r97)));
	I1 = (double)0.5 * (double)log((C*C + (double)13*C + (double)18)/(double)18) - (double)13/(double)2 * I2;
	
	P1 = (-C+(C-(double)1)*I1 + (double)2*((double)7*C+(double)9)*I2);
	P2 = (C*C/(double)2 + C + ((double)5-C)*I1 - (double)18*(C+(double)1)*I2)/(double)nsam;
	P3 = (-(C*C)/(double)2 + (double)2*C - (double)2*(C+(double)9)*I1 - (double)4*((double)2*C+(double)9)*I2)/((double)nsam*(double)nsam);
	P4 = ((double)1/(double)nsam)*((double)-2*C + (double)2*(C+(double)7)*I1 + (double)12*(C+(double)3)*I2)/((double)nsam*(double)nsam);
	fC = ((double)2/(C*C)) * (P1 + P2 + P3 + P4);
	/*
	fC = ((double)2/(C*C)) * ((-C+(C-(double)1)*I1 + (double)2*((double)7*C+(double)9)*I2)
	                           +(C*C/(double)2 + C + ((double)5-C)*I1 - (double)18*(C+(double)1)*I2)/(double)nsam
							   +(-(C*C)/(double)2 + (double)2*C - (double)2*(C+(double)9)*I1 - (double)4*((double)2*C+(double)9)*I2)/((double)nsam*(double)nsam)
							   +((double)1/(double)nsam)*((double)-2*C + (double)2*(C+(double)7)*I1 + (double)12*(C+(double)3)*I2)/((double)nsam*(double)nsam));
	*/
	return gCn - fC;
}

int zbracC(double *x1,double *x2,double gCn,int nsam)
{
	/* Based on Numerical Recipes in C. Press et al. 1992. 
    We need *x1 and *x2 be the range where a root is within them. We expand geometrically the range until finding
	(one a positive and one a negative value). If not, return 0.
	*/

    double f1,f2;
    double functiongCn0(double,double,int);
    int k=60;
    
    if(*x1 == *x2) return 0;

    f1 = functiongCn0(*x1,gCn,nsam);
    f2 = functiongCn0(*x2,gCn,nsam);
	
	if(f1*f2 < (double)0) return 1;

    while(k--) {
        if(fabs(f1) < fabs(f2)) {
            *x1 += (double)1.5 * (*x1 - *x2);
            f1 = functiongCn0(*x1,gCn,nsam);
        }
        else {
            *x2 += (double)1.5 * (*x2 - *x1);
            f2 = functiongCn0(*x2,gCn,nsam);
        }
        if(f1*f2 < (double)0) return 1;
    }
    return 0;
}

double zriddrC(double xlow,double xhigh,double xacc,double gCn,int nsam)
{
	/* Based on Numerical Recipes in C. Press et al. 1992., p. 358 an on
	Ridders, 1979, IEEE Transactions on Circuits and systems, Vol. Cas-26, No. 11, pp. 979-980.
	*/
    int k=60;
	double flow,fhigh;
	double f1,f2,f3,f4;
	double x1,x2,x3,x4;
	double den,num,nsign;
    double functiongCn0(double,double,int);

    flow  = functiongCn0(xlow,gCn,nsam);
    fhigh = functiongCn0(xhigh,gCn,nsam);

	if(flow  == (double)0) return xlow;
	if(fhigh == (double)0) return xhigh;
	if(flow*fhigh > (double)0) 
		return (double)-1e32;
	
	x1 = xlow;
	x2 = xhigh;
	f1 = flow;
	f2 = fhigh;
		
	while(k--) {
		x3 = (x1+x2)/(double)2;
		f3 = functiongCn0(x3,gCn,nsam);
		if(f1 - f2 < (double)0) nsign = (double)-1;
		else nsign = (double)1;
		num = (x3-x1) * f3 * nsign;
		den = (double)sqrt((double)f3*(double)f3 - (double)f1*(double)f2);
		if(den <= xacc && -den <= xacc) return x3;
		x4 = x3 + num/den;
		f4 = functiongCn0(x4,gCn,nsam);
		if(f4 <= xacc && -f4 <= xacc) return x4;
		if(f3*f4<(double)0) {
			x1 = x3;
			f1 = f3;
			x2 = x4;
			f2 = f4;
		}
		else {
			if(f1*f4<(double)0) {
				x2 = x4;
				f2 = f4;
			}
			else {
				if(f2*f4<(double)0) {
					x1 = x4;
					f1 = f4;
				}
			}
		}
		if(fabs(x1-x2) <= xacc) return x1;
	}	
	return (double)-1e32;
}

/*Wall's program for calculating minimum recombination events*/
int Min_rec_obs(int x,struct statistics **matrix,char *matrix_pol, int *matrix_freq,int n_samp, int nnnsam,int nnnoutg,int bial_sites)
{  /* Calculate min # rec. events */
	int a, b, c, e, gtest, flag = 0;
	/*int h;*/
	int t11,t12,t21,t22;
  
	if (bial_sites < 2 || x >= (bial_sites-1)) return (0);
	
	for (a=x+1; a<bial_sites; ++a) {
        while(a < bial_sites) { /*calculate a*//*not shared counted...*/
            if(*(matrix_freq+a) >/*=*/ 0 && *(matrix_freq+a) < nnnsam) break;/*fixed and pol. outgroup not counted*/
            else a++;
        }
		if(a < bial_sites) {
			for (b=x; b<a; ++b) {
				while(b < a) { /*calculate b*//*not shared counted...*/
					if(*(matrix_freq+b) >/*=*/ 0 && *(matrix_freq+b) < nnnsam) break;/*fixed and pol. outgroup not counted*/
					else b++;
				}
				if(b < a) {
					for(e=nnnoutg;e<n_samp;e++) {
						if(e==nnnoutg) {
							t21 = *(matrix_pol+((b*(n_samp))+e));
							t11 = *(matrix_pol+((a*(n_samp))+e));
						}
						else {
							if(*(matrix_pol+((b*(n_samp))+e)) != t21) t22 = *(matrix_pol+((b*(n_samp))+e));
							if(*(matrix_pol+((a*(n_samp))+e)) != t11) t12 = *(matrix_pol+((a*(n_samp))+e));
						}
					}
					
					gtest = 0;
					for(e=nnnoutg;e<n_samp;e++) {
						if (*(matrix_pol+((b*(n_samp))+e)) == t21 && *(matrix_pol+((a*(n_samp))+e)) == t11) {
							++gtest;
						break;
						}
					}
					for(e=nnnoutg;e<n_samp;e++) {
						if (*(matrix_pol+((b*(n_samp))+e)) == t21 && *(matrix_pol+((a*(n_samp))+e)) == t12) {
							++gtest;
							break;
						}
					}
					for(e=nnnoutg;e<n_samp;e++) {
						if (*(matrix_pol+((b*(n_samp))+e)) == t22 && *(matrix_pol+((a*(n_samp))+e)) == t11) {
							++gtest;
							break;
						}
					}
					for(e=nnnoutg;e<n_samp;e++) {
						if (*(matrix_pol+((b*(n_samp))+e)) == t22 && *(matrix_pol+((a*(n_samp))+e)) == t12) {
							++gtest;
							break;
						}
					}
					if (gtest == 4) {
						flag = 1;
						break;
					}
				}
			}
		}
		if (flag == 1) break;
	}
	if (a >= bial_sites) return (0);
	else {
		c = Min_rec_obs(a,matrix,matrix_pol,matrix_freq,n_samp,nnnsam,nnnoutg,bial_sites);
		return (1+c);
	}
	return 0;
}
