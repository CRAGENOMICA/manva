/*
 *  get_obsdatastats.c
 *  MuLoNeTests
 *
 *  Created by sebas on Mon Feb 24 2003.
 *
 */

#include "MuLoNeTests.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int get_obsstats(struct statistics **matrix,struct statmulo **matrixml/*,int *observed_data */,int *outgroup,int *n_loci,FILE *file_output/*,FILE *file_input*/, char *name_outgroup,
	int n_samp, long int n_site, char *DNA_matr, char **names,int n_excl,int lian_in,struct lianinput **liandata,double factor_chrn, double *matrix_sizepos,double *matrix_segrpos)
{
    int calc_statistics(struct statistics **,int *,int * /*,FILE * */,int,double,char *,int *,int,int,int,char **,char *,int,struct lianinput **,double,int);
    int calc_neutests(struct statistics **,int *,int *, int * /*, FILE * */);
    int calcobsmultiloci(struct statistics *,struct statmulo **,int *, int * /*, FILE * */);
    
    static char *matrix_pol = 0;
    static int *matrix_freq = 0;
    static int *unic = 0;
    static int *matrix_pos = 0;
    /*static int *mhitbp = 0;*/
    /*static long int maxnsite = 128;*/
    static int maxbialsites = 256;
    static int maxnsamp = 32;
    static int *nnsam = 0;
    static int nnnsam;
    static int *nnoutg = 0;
    static int nnnoutg;

    double algsites;
    int bial_sites = 0;
    int mhits = 0;

    int m_0,k_0,m_1,k_1;
    int v,w,x,y,z;

    int v0,v1,y0,y1,b,c,d;
    char *strin;
    double _sites;
    
    if(matrix_pol == 0) {
        /* Two pointers indicating the samples that are current samples and the outgroup samples */
        if((nnsam = (int *) calloc(maxnsamp,sizeof(int))) == 0) {
            puts("Error: memory not reallocated. get_obsstats.1"); 
            return(0);
        }
        if((nnoutg = (int *) calloc(maxnsamp,sizeof(int))) == 0) {
            puts("Error: memory not reallocated. get_obsstats.2"); 
            return(0);
        }
        /* matrix of polymorphisms: only 0 and 1 */
        if((matrix_pol = (char *) calloc (maxnsamp*maxbialsites, sizeof(char))) == 0) {
            puts("Error: memory not reallocated. get_obsstat.3"); 
            return(0);
        }
        /* indicates the position and the frequency */
        if((matrix_pos = (int *) calloc (maxbialsites, sizeof(int))) == 0) {
            puts("Error: memory not reallocated. get_obsstat.4"); 
            return(0);
        }
        if((matrix_freq = (int *) calloc (maxbialsites, sizeof(int))) == 0) {
            puts("Error: memory not reallocated. get_obsstat.5"); 
            return(0);
        }
        if((unic = (int *) calloc (maxnsamp, sizeof(int))) == 0) {
            puts("Error: memory not reallocated. get_obsstat.6"); 
            return(0);
        }
        /* indicates the position of the mhits */
		/*
        if((mhitbp = (int *) calloc (maxnsite, sizeof(int))) == 0) {
            puts("Error: memory not reallocated. get_obsstat.6"); 
            return(0);
        }
		*/
    }
    if(n_samp > maxnsamp) {
        /* Reallocation in case the value be larger than specified */
        if((matrix_pol = (char *) realloc (matrix_pol,(n_samp*maxbialsites)*sizeof(char))) == 0) {
            puts("Error: memory not reallocated. get_obsstat.7b"); 
            return(0);
        }
        if((nnsam = (int *) realloc(nnsam,n_samp*sizeof(int))) == 0) {
            puts("Error: memory not reallocated. get_obsstats.9"); 
            return(0);
        }
        if((nnoutg = (int *) realloc(nnoutg,n_samp*sizeof(int))) == 0) {
            puts("Error: memory not reallocated. get_obsstats.10"); 
            return(0);
        }
        if((unic = (int *) realloc(unic,n_samp*sizeof(int))) == 0) {
            puts("Error: memory not reallocated. get_obsstats.11"); 
            return(0);
        }
        maxnsamp = n_samp;
    }
	/*
    if(n_site > maxnsite) {
        if((mhitbp = (int *) realloc (mhitbp,n_site*sizeof(int))) == 0) {
            puts("Error: memory not reallocated. get_obsstat.8"); 
            return(0);
        }
        maxnsite = n_site;
    }
	*/    
    algsites = (double)n_site; /* number of aligned positions, excluding gaps and mhits */
    nnnsam = nnnoutg = 0;
    
    /* calculate number of samples in outgroup and in the current sample */
    if(*outgroup) {
        /* compare names with outgroup names */
        for(x=0;x<n_samp;x++) {
            if((strin = strstr(names[x],name_outgroup)) != 0) {
                nnoutg[nnnoutg] = x;
                nnnoutg += 1;
            }
            else {
                nnsam[nnnsam] = x;
                nnnsam += 1;                
            }
        }
        if(nnnsam < 2 || nnnoutg == 0) {
            printf(" n_samples: %d, n_outgroups: %d, n_excluded %d.",nnnsam,nnnoutg,n_excl);
            if(file_output) 
                fprintf(file_output," n_samples: %d, n_outgroups: %d, n_excluded %d.",nnnsam,nnnoutg,n_excl);
            if(nnnsam < 2) {
                printf(" NOT ENOUGH SAMPLES.");
                if(file_output) fprintf(file_output," NOT ENOUGH SAMPLES.");
            }
            if(nnnoutg == 0) {
                printf(" NO OUTGROUP AVAILABLE.");
                if(file_output) fprintf(file_output," NO OUTGROUP AVAILABLE.");
            }
            return(0);
        }
    }
    else {
        /* if no outgroup, then all sequences are samples */
        if(n_samp < 2) {
            printf(" n_samples: %d, n_outgroups: %d, n_excluded %d.",nnnsam,nnnoutg,n_excl);
            if(file_output) 
                fprintf(file_output," n_samples: %d, n_outgroups: %d, n_excluded %d.",nnnsam,nnnoutg,n_excl);
            return(0);
            if(nnnsam < 2) {
                printf(" NOT ENOUGH SAMPLES.");
                if(file_output) fprintf(file_output," NOT ENOUGH SAMPLES.");
            }
        }
        for(x=0;x<n_samp;x++) nnsam[x] = x;
        nnnsam = n_samp;
    }
    /*init unics*/
    for(y=0;y<n_samp;y++) {
        unic[y] = 0;
    }
    /* find positions that are biallelic, excluding the rest. */
    /* IMPORTANT: multiple hits are eliminated from analysis */
    _sites = (double)0;
    matrix[0][*n_loci].transitions = 0;
	matrix[0][*n_loci].transversions = 0;
	
	for(x=0;x<n_site;x++) {
		/*eliminate those positions that are not included in analysis*/
		if(matrix_sizepos[x] == (double)0) {
			algsites -= (double)1;
			continue;
		}
		/*look for biallelic*/
        m_0 = k_0 = m_1 = k_1 = y = z = v = 0;
        do {
            m_0 = *(DNA_matr+(((unsigned long)n_site*(unsigned long)0)+(unsigned long)x));
            w   = *(DNA_matr+(((unsigned long)n_site*(unsigned long)y)+(unsigned long)x));            
            if(w > 48 + 4) { /* the position y (also y == 0) is not a, g, t or c */
                z = 1;
                v = 0;
                algsites -=matrix_sizepos[x];
            }
            else { /* the position y is a, g, t or c */
                if(w != m_0 && w != k_0) { /* if position 0 and y are different, and also different of the k_0 */
                    if(k_0 == 0) { /* if k_0 is 0 (initialized), the k_0 will be w */
                        k_0 = w;
                        k_1++; /* counting the numbers of w different from m_0 and not > 4 and not mhits */
                    }
                    else { /* if k_0 has already a value */
                        if((w <= 48 + 4) && (m_0 <= 48 + 4))  
                            v = 1; /* if w and m_0 are different and are also different from k_0 = mhit */
                    }
                }
                else { /* in case w=m_0 or w=k_0 */
                    if(w == m_0) m_1++; /* if w equal to m_0 count m_1 */
                    else k_1++; /* if w is different from m_0 */
                }
            }
            y++;
        } while(z==0 && y < n_samp);
        
        /*count transition and transversions*/
		if(k_0 != 0 && k_0 != m_0 && v == 0 && z == 0) {
			if((k_0 == '1' && m_0 == '2') || (k_0 == '2' && m_0 == '1') || 
			   (k_0 == '2' && m_0 == '3') || (k_0 == '3' && m_0 == '2')) 
					matrix[0][*n_loci].transitions += 1;
			else matrix[0][*n_loci].transversions += 1;
		}
		
		/* mhit position */
        if(v == 1) {
            /*mhitbp[mhits] = x+1;*/
            mhits++;
            algsites -= matrix_sizepos[x];
            z = 1;
        }
        if(z==0 && v == 0) {
			_sites += matrix_sizepos[x];
			algsites -= (double)1 - matrix_sizepos[x];
		}
        /* do the matrix of biallelic positions: FIRST the outgroup/s, SECOND the current samples */
        if(z==0 && m_1<n_samp && matrix_segrpos[x]/*to eliminate biallelic syn/nsyn not desired*/) {
			if(*outgroup) {
                /* with outgroup: 0 is the outgroup, 1 new mutations (be careful with muts in outgroup) */
                y0 = nnnoutg;
                b = 0;
                do { /* look for polymorphisms in the outgroup */
                    y0 --;
                    y1 = nnoutg[y0];
                    if(*(DNA_matr+(((unsigned long)n_site*(unsigned long)y1)+(unsigned long)x)) == k_0) b++;
				} while(y0);
                if(b != 0 && b != nnnoutg) { /* in case the outgroup be polymorphic */
                    v0 = nnnsam;
                    c = 0;
                    do { /* look for polymorphisms in the sample not outgroup */
                        v0 --;
                        v1 = nnsam[v0];
                        if(*(DNA_matr+(((unsigned long)n_site*(unsigned long)y1)+(unsigned long)x)) == k_0) c++;
                    }while(v0);
                    if(c == 0 || c == nnnsam) { /* in case the sample is not polymorphic */
                        for(y=nnnoutg;y<n_samp;y++) *(matrix_pol+((bial_sites*(n_samp))+y)) = '0'; /*monomorphic = 0*/
                        if(c == 0) k_0 = m_0; /* indicate what nt. is 0 */
                        d = 0;
                        for(y=0;y<nnnoutg;y++) {
                            if(*(DNA_matr+(((unsigned long)n_site*(unsigned long)nnoutg[y])+(unsigned long)x)) == k_0)
                                *(matrix_pol+((bial_sites*(n_samp))+y)) = '0';
                            else {
                                *(matrix_pol+((bial_sites*(n_samp))+y)) = '1';
                                d -= 1; 
                            }
                        }
                        *(matrix_pos+bial_sites) = x+1;
                        *(matrix_freq+bial_sites) = d; /* negative value for polymorphism in the outgroup */
                    }
                    else { /* shared polymorphisms */
                        if(m_1 > k_1)  {
                            k_0 = m_0;
                            k_1 = m_1;
                        }
                        for(y=0;y<n_samp;y++) { /* shared, we don't know what is the ancestor (0 or 1) */
                            if(*(DNA_matr+(((unsigned long)n_site*(unsigned long)y)+(unsigned long)x)) == k_0)
                                *(matrix_pol+((bial_sites*(n_samp))+y)) = '0';
                            else {
                                *(matrix_pol+((bial_sites*(n_samp))+y)) = '1';
                            }
                        }
                        *(matrix_pos+bial_sites) = x+1;
                        *(matrix_freq+bial_sites) = 0; /* 0 for shared. special case. */
                        /*calculate unics: with shared we do not count*/
                        /*
                        d = 0;
                        for(y=nnnoutg;y<n_samp;y++) {
                            if(*(DNA_matr+(((unsigned long)n_site*(unsigned long)y)+(unsigned long)x)) != k_0)
                                d += 1;
                        }
                        for(y=nnnoutg;y<n_samp;y++) {
                            if((d==1 && *(matrix_pol+((bial_sites*(n_samp))+y)) == '1') ||
                            (d==(n_samp-nnnoutg)-1 && *(matrix_pol+((bial_sites*(n_samp))+y)) == '0'))
                                unic[y-nnnoutg] += 1;
                        }
                        */
                    }    
                }
                else { /* in case the outgroup be monomorphic */
                    for(y=0;y<nnnoutg;y++) *(matrix_pol+((bial_sites*(n_samp))+y)) = '0'; /*monomorphic: outg = 0*/
                    if(b == 0) k_0 = m_0; /* indicate what nt. is 0 */
                    d = 0;
                    for(y=nnnoutg;y<n_samp;y++) {
                        if(*(DNA_matr+(((unsigned long)n_site*(unsigned long)nnsam[y-nnnoutg])+(unsigned long)x)) == k_0)
                            *(matrix_pol+((bial_sites*(n_samp))+y)) = '0';
                        else {
                            *(matrix_pol+((bial_sites*(n_samp))+y)) = '1';
                            d += 1;
                        }
                    }
                    *(matrix_pos+bial_sites) = x+1;
                    *(matrix_freq+bial_sites) = d;
                    /*calculate unics*/
                    for(y=nnnoutg;y<n_samp;y++) {
                        if((d==1 && *(matrix_pol+((bial_sites*(n_samp))+y)) == '1') ||
                           (d==(n_samp-nnnoutg)-1 && *(matrix_pol+((bial_sites*(n_samp))+y)) == '0'))
                            unic[y-nnnoutg] += 1;
                    }
                }          
            }
            else { /* not outgroup */
                /* k_0 will be the large frequency */
                if(m_1 > k_1)  {
                    k_0 = m_0;
                    k_1 = m_1;
                }
                d = 0;
                for(y=0;y<n_samp;y++) { /* non-outgroup: 0 indicates higher frequency */
                    if(*(DNA_matr+(((unsigned long)n_site*(unsigned long)y)+(unsigned long)x)) == k_0)
                        *(matrix_pol+((bial_sites*(n_samp))+y)) = '0';
                    else {
                        *(matrix_pol+((bial_sites*(n_samp))+y)) = '1';
                        d += 1;
                    }
                }
                *(matrix_pos+bial_sites) = x+1;
                *(matrix_freq+bial_sites) = d;
                /*calculate unics*/
                for(y=0;y<n_samp;y++) {
                    if((d==1 && *(matrix_pol+((bial_sites*(n_samp))+y)) == '1') ||
                       (d==n_samp-1 && *(matrix_pol+((bial_sites*(n_samp))+y)) == '0')) 
                        unic[y] += 1;
                }
            }
            /* one more position */
            bial_sites++;
            /* reallocations */
            if(bial_sites == maxbialsites) {
                if(maxbialsites == 32767) {
                    printf("\n Sorry, it is only accepted a maximum of 32767 biallelic sites per loci. It has been cut at position %d.",x+1);
                    if(file_output) fprintf(file_output,"\n Sorry, it is only accepted a maximum of 32767 biallelic sites per loci. It has been cut at position %d.",x+1);
                    break;
                }
                else if(maxbialsites > 32767 - 128) maxbialsites = 32767;
                    else maxbialsites += 128;
                if((matrix_pol = realloc(matrix_pol,(maxnsamp)*(maxbialsites)*sizeof(char))) == 0) {
                    puts("Error: memory not reallocated. get_obsstat.11"); 
                    return(0);
                }
                if((matrix_pos = realloc(matrix_pos,(maxbialsites)*sizeof(int))) == 0) {
                    puts("Error: memory not reallocated. get_obsstat.12");
                    return(0);
                }
                if((matrix_freq = realloc(matrix_freq,(maxbialsites)*sizeof(int))) == 0) {
                    puts("Error: memory not reallocated. get_obsstat.13");
                    return(0);
                }
            }
        }
    } 
	/*
	if(_sites != algsites) {
		printf("CHECK BUG in get_obsatastats.c");
	}
	*/
    printf(" n_samples: %d, n_outgroups: %d, n_excluded: %d, valid sites: %.2f.",nnnsam,nnnoutg,n_excl,_sites);
    if(file_output) fprintf(file_output," n_samples: %d, n_outgroups: %d, n_excluded: %d, valid sites: %.2f.",
        nnnsam,nnnoutg,n_excl,_sites);
    if(_sites == (double)0) {
        printf(" Not valid sites available in this file:  ");
        if(file_output) fprintf(file_output," Not valid sites available in this file: ");
        return 0;
    }
    /* calculate all necessary statistics and keep in matrix */
    if(!(calc_statistics(matrix,outgroup,n_loci/*,file_output*/,n_samp,algsites,
         matrix_pol,matrix_freq,nnnsam,nnnoutg,bial_sites,names,name_outgroup,lian_in,liandata,factor_chrn,mhits))) {
        return 0;
    }
    if(!(calc_neutests(matrix,unic,outgroup,n_loci/*,file_output*/))) {        
        return 0;
    }
    if(!(calcobsmultiloci(*matrix,matrixml,outgroup,n_loci/*,file_output*/))) {
        return 0;
    }
    return 1;
}
