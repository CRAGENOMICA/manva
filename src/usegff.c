
#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char tripletsU[64][3] =
{
	{"UUU"},/*0*/
	{"UUC"},/*1*/
	{"UUA"},/*2*/
	{"UUG"},/*3*/
	{"UCU"},/*4*/
	{"UCC"},/*5*/
	{"UCA"},/*6*/
	{"UCG"},/*7*/
	{"UAU"},/*8*/
	{"UAC"},/*9*/
	{"UAA"},/*10*/
	{"UAG"},/*11*/
	{"UGU"},/*12*/
	{"UGC"},/*13*/
	{"UGA"},/*14*/
	{"UGG"},/*15*/
	{"CUU"},/*16*/
	{"CUC"},/*17*/
	{"CUA"},/*18*/
	{"CUG"},/*19*/
	{"UCU"},/*20*/
	{"CCC"},/*21*/
	{"CCA"},/*22*/
	{"CCG"},/*23*/
	{"CAU"},/*24*/
	{"CAC"},/*25*/
	{"CAA"},/*26*/
	{"CAG"},/*27*/
	{"CGU"},/*28*/
	{"CGC"},/*29*/
	{"CGA"},/*30*/
	{"CGG"},/*31*/
	{"AUU"},/*32*/
	{"AUC"},/*33*/
	{"AUA"},/*34*/
	{"AUG"},/*35*/
	{"ACU"},/*36*/
	{"ACC"},/*37*/
	{"ACA"},/*38*/
	{"ACG"},/*39*/
	{"AAU"},/*40*/
	{"AAC"},/*41*/
	{"AAA"},/*42*/
	{"AAG"},/*43*/
	{"AGU"},/*44*/
	{"AGC"},/*45*/
	{"AGA"},/*46*/
	{"AGG"},/*47*/
	{"GUU"},/*48*/
	{"GUC"},/*49*/
	{"GUA"},/*50*/
	{"GUG"},/*51*/
	{"GCU"},/*52*/
	{"GCC"},/*53*/
	{"GCA"},/*54*/
	{"GCG"},/*55*/
	{"GAU"},/*56*/
	{"GAC"},/*57*/
	{"GAA"},/*58*/
	{"GAG"},/*59*/
	{"GGU"},/*60*/
	{"GGC"},/*61*/
	{"GGA"},/*62*/
	{"GGG"},/*63*/
};

char tripletsN[64][3] =
{
	{"111"},
	{"112"},
	{"114"},
	{"113"},
	{"121"},
	{"122"},
	{"124"},
	{"123"},
	{"141"},
	{"142"},
	{"144"},
	{"143"},
	{"131"},
	{"132"},
	{"134"},
	{"133"},
	{"211"},
	{"212"},
	{"214"},
	{"213"},
	{"221"},
	{"222"},
	{"224"},
	{"223"},
	{"241"},
	{"242"},
	{"244"},
	{"243"},
	{"231"},
	{"232"},
	{"234"},
	{"233"},
	{"411"},
	{"412"},
	{"414"},
	{"413"},
	{"421"},
	{"422"},
	{"424"},
	{"423"},
	{"441"},
	{"442"},
	{"444"},
	{"443"},
	{"431"},
	{"432"},
	{"434"},
	{"433"},
	{"311"},
	{"312"},
	{"314"},
	{"313"},
	{"321"},
	{"322"},
	{"324"},
	{"323"},
	{"341"},
	{"342"},
	{"344"},
	{"343"},
	{"331"},
	{"332"},
	{"334"},
	{"333"},
};


int use_gff(char **name_fileinputfoldergff,char **name_fileinputgff,char **name_fileinputfolderGFF/*,char **name_fileinputGFF*/,char *subset_positions,
		/* int ifgencode,char *codename, */char *genetic_code,double *matrix_sizepos,int n_samp,long int n_site,char *DNA_matr,double *matrix_segrpos,
		FILE *file_output)
{
	FILE *file_gff;
	char *row,*f,cstrand[1],cframe[2],aaseq[1],aaput[1];
	int i,j,n,nn,m,q,nrows,stop,countpath;
	char fields[9][256];
	struct valuesgff 
	{
		char filename[256];
		char feature[256];
		char strand[1];
		long int start;
		long int end;
		char frame[1];
		char seqname[256];
	}*fieldsgff;
	char *seqid/*, *fileid*/;
	double *cmat,*cmatnc,*cmatsil;
	long int ii,ii2,k,startframe,endframe,endtrp;
	char *cod3n,*cod3put;
	int tripletnsamp(char *,char *,char,double *,int,long int,long int,long int,FILE *);
	long int start,end;
	double cod3f,cod3ft[3];
	
	int gg;
	char *pst;
	
	/*read gff file: read name_fileinputfoldergff, if it does not work, read name_fileinputfolderGFF*/
	/*name_fileinputgff/name_fileinputGFF is only for displaying the name on the screen.*/
	/*fields to read in GFF-format: seqname, noread(source), feature, start, end, noread(score), strand, frame, noread(rest)*/
	if((file_gff = fopen (*name_fileinputfoldergff,"r")) == 0) {
		if((file_gff = fopen (*name_fileinputfolderGFF,"r")) == 0) {
			printf("\n  It is not possible to open the file %s",*name_fileinputfoldergff);
			return 0; /*error*/
		}
	}
	if(file_gff) {
		/*init*/
		if(!(f = (char *)malloc(BUFSIZ))) {
			puts("\nError: memory not reallocated. use_gff.1 \n");
			return 0; /*error*/
		}
		if(!(row = (char *)malloc(1024*sizeof(char)))) {
			puts("\nError: memory not reallocated. use_gff.2 \n");
			return 0; /*error*/
		}
		if(!(fieldsgff = (struct valuesgff *)calloc(1,sizeof(struct valuesgff)))) {
			puts("\nError: memory not reallocated. use_gff.3 \n");
			return 0; /*error*/
		}
		setbuf(file_gff,f);
		/*read rows*/
		nrows = 0;
		while(1) {
			fgets(row, 1024*sizeof(char), file_gff);
			if(feof(file_gff)) break;
			/*i=0;*/
			/*while((row[i] = fgetc(file_gff)) != 0 && row[i] != 10 && row[i] != 13 && feof(file_gff)!= 1) {*/
			/*	i++;*/
			/*	if(i >= 1024) while((j = fgetc(file_gff)) != 0  && j != 10 && j != 13 && (feof(file_gff)) != 1); */
			/*}*/
			/*row[i] = 0;*/
			
			i=0;
			while(row[i] == 32 || row[i] == '\t') {
				if(row[i] == 10 || row[i] == 13 || row[i] == 0) break;
				i++;
				if(i >= 1024) break;
			}
			if(row[i] == '#') continue;
			if(i >= 1024) continue;
						
			/*include fields in variables*/
			j = k = 0;
			while(row[i] != 10 && row[i] != 13 && row[i] != 0) {
				while(row[i] == 32 || row[i] == '\t') {
					if(row[i] == 10 || row[i] == 13 || row[i] == 0) break;
					i++;
					if(i >= 1024) break;
				}
				k=0;
				while(row[i] != 32 && row[i] != '\t' && row[i] != 10 && row[i] != 13 && row[i] != 0) {
					fields[j][k] = row[i];
					k++;
					if(k >= 256) break;
					i++;
					if(i >= 1024) break;
				}
				fields[j][k] = '\0';
				j++;
				if(j==9) {
					break;
				}
			}
			/*add fields in struct*/
			for(n=0;n<j;n++) {
				switch(n) {
					case 0:
						fieldsgff[nrows].filename[0] = '\0';
						strncat(fieldsgff[nrows].filename,fields[n],256*sizeof(char));
						break;
					case 1:
						break;
					case 2:
						fieldsgff[nrows].feature[0] = '\0';
						strncat(fieldsgff[nrows].feature,fields[n],256*sizeof(char));
						break;
					case 3:
						fieldsgff[nrows].start = atol(fields[n]);
						break;
					case 4:
						fieldsgff[nrows].end = atol(fields[n]);
						break;
					case 5:
						break;
					case 6:
						fieldsgff[nrows].strand[0] = fields[n][0];
						break;
					/*case 7:
						fieldsgff[nrows].frame[0] = fields[n][0];
						break;*/
					case 7:
						if(fields[n][0] == '1') {
							fieldsgff[nrows].frame[0] = '2';
						}
						else {
							if(fields[n][0] == '2') {fieldsgff[nrows].frame[0] = '1';
							}
							else {
								fieldsgff[nrows].frame[0] = fields[n][0];
							}
						}
						break;
					case 8:
						fieldsgff[nrows].seqname[0] = '\0';
						pst = strstr(fields[n],"Parent=");
						if(pst) {
							while(*pst != '=') pst++;
							pst++;
							gg = 0;
							while(*pst != ';' && *pst != '\0') {
								fieldsgff[nrows].seqname[gg] = *pst;
								pst++;
								gg++;
								if(gg==255) break;
									
							}
							fieldsgff[nrows].seqname[gg] = '\0';
						}
						break;
				}
			}
			/*minimum five fields: seqname,whatever(dot),feature,start,end.*/
			if(j < 9) continue;
			/*filtering ...*/
			if(fieldsgff[nrows].start < (long int)1)  {
				printf("GFF file error: start is lower than 1. Row not analyzed. ");
				if(file_output) fprintf(file_output,"GFF file error: start is larger than end. Row not analyzed. ");
				continue;
			}
			if(fieldsgff[nrows].end > (long int)n_site)  {
				printf("GFF file error: end is larger than number of total sites. Row not analyzed. ");
				if(file_output) fprintf(file_output,"GFF file error: start is larger than end. Row not analyzed. ");
				continue;
			}
			if(fieldsgff[nrows].start > fieldsgff[nrows].end) {
				printf("GFF file error: start is larger than end. Row not analyzed. ");
				if(file_output) fprintf(file_output,"GFF file error: start is larger than end. Row not analyzed. ");
				continue;
			}
			/*internally starts from 0 to n_site-1*/
			fieldsgff[nrows].start -= (long int)1;
			fieldsgff[nrows].end -= (long int)1;
			/*we assume by default starting from strand '+' and frame 0*/
			/**/
			if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0) {
				if(strcmp(fieldsgff[nrows].feature,"CDS") == 0) {
					if(fieldsgff[nrows].strand[0] != '\0' && fieldsgff[nrows].strand[0] != '+' && 
					   fieldsgff[nrows].strand[0] != '-'  && fieldsgff[nrows].strand[0] != '.') {
							fieldsgff[nrows].strand[0] = 32;
					}
					if(fieldsgff[nrows].frame[0] != '\0' && fieldsgff[nrows].frame[0] != '1' && 
					   fieldsgff[nrows].frame[0] != '2'  && fieldsgff[nrows].frame[0] != '0' && 
					   fieldsgff[nrows].strand[0] != '.') {
							fieldsgff[nrows].frame[0] = 32;
					}
				}
			}
			/**/
			/*if row accepted*/
			nrows += 1;
			if(!(fieldsgff = (struct valuesgff *)realloc(fieldsgff,(nrows+1)*sizeof(struct valuesgff)))) {
				puts("\nError: memory not reallocated. use_gff.3 \n");
				free(row);
				free(f);
				free(fieldsgff);
				fclose(file_gff);
				return 0; /*error*/
			}		
			if(feof(file_gff)) break;
		}
		free(row);
		free(f);
		fclose(file_gff);
		
		for(n=0;n<9;n++) {
			switch(n) {
				case 0:
					fieldsgff[nrows].filename[0] = '\0';
					strncat(fieldsgff[nrows].filename,*name_fileinputgff,256*sizeof(char));
					break;
				case 1:
					break;
				case 2:
					fieldsgff[nrows].feature[0] = '\0';
					strncat(fieldsgff[nrows].feature,"all",256*sizeof(char));
					break;
				case 3:
					fieldsgff[nrows].start = (long int)0;
					break;
				case 4:
					fieldsgff[nrows].end = (long int)(n_site-1);
					break;
				case 5:
					break;
				case 6:
					fieldsgff[nrows].strand[0] = 32;
					break;
				case 7:
					fieldsgff[nrows].frame[0] = 32;
					break;
				case 8:
					fieldsgff[nrows].seqname[0] = '\0';
					break;
			}
		}
		nrows += 1;
		
		/*count positions for each seqname, if exons from different genes are overlapped, 
			the positions counted will NOT be correctly assigned!!*/

		/*set to zero all values in matrix_sizepos*/
		if(strcmp(subset_positions,"noncoding") == 0) for(ii=0;ii<n_site;ii++) matrix_sizepos[ii] = (double)1;
		else for(ii=0;ii<n_site;ii++) matrix_sizepos[ii] = (double)0;
		/*init a current matrix_sizepos for each seqname*/
		if((cmat = (double *)calloc((unsigned long)n_site,sizeof(double))) == 0) {
			puts("\nError: memory not reallocated. use_gff.4 \n");
			return 0; /*error*/
		}
		if((cmatnc = (double *)calloc((unsigned long)n_site,sizeof(double))) == 0) {
			puts("\nError: memory not reallocated. use_gff.4b \n");
			return 0; /*error*/
		}
		if((cmatsil = (double *)calloc((unsigned long)n_site,sizeof(double))) == 0) {
			puts("\nError: memory not reallocated. use_gff.4c \n");
			return 0; /*error*/
		}
		for(ii=0;ii<n_site;ii++) cmatsil[ii] = (double)1;
		/*init a 3*n_samp matrix for syn/nsyn triplets*/
		if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0) {
			if((cod3n = (char *)calloc((unsigned long)(3*(n_samp)),sizeof(char))) == 0) {
				puts("\nError: memory not reallocated. use_gff.5 \n");
				return 0; /*error*/
			}
			if((cod3put = (char *)calloc((unsigned long)3,sizeof(char))) == 0) {
				puts("\nError: memory not reallocated. use_gff.5 \n");
				return 0; /*error*/
			}
		}
		/*check that CDS regions are not overlapped*/
		if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0) {
			for(n=0;n<nrows;n++) {
				for(m=0;m<nrows;m++) {
					if(m==n) continue;
					/*if(strcmp(fieldsgff[m].filename,fieldsgff[n].filename) == 0) {*/
						if(strcmp(fieldsgff[m].feature,"CDS") == 0 && strcmp(fieldsgff[n].feature,"CDS") == 0) {
							if(fieldsgff[m].start > fieldsgff[n].start && fieldsgff[m].start < fieldsgff[n].end) {
								printf(" Overlapping CDS regions. ");
								if(file_output) fprintf(file_output," Overlapping CDS regions. ");
								return 0;
							}
							if(fieldsgff[m].end > fieldsgff[n].start && fieldsgff[m].end < fieldsgff[n].end) {
								printf(" Overlapping CDS regions. ");
								if(file_output) fprintf(file_output," Overlapping CDS regions. ");
								return 0;
							}
						}
					/*}*/
				}
			}
		}
		
		/*take the selected regions*/
		for(n=0;n<nrows;n++) {
			/*take seqname of each row and check if there are coincidences before, if yes, not count*/
			seqid = fieldsgff[n].seqname;
			/*fileid = fieldsgff[n].filename;*/
			i=n-1;
			k=0;
			while(i>=0) {
				if(strcmp(seqid,fieldsgff[i].seqname) == 0 /*&& strcmp(fileid,fieldsgff[i].filename) == 0*/) {
					k = 1;
					break;
				}
				i--;
			}
			if(k) continue;
			
			/*include in current matrix_sizepos all rows with the same seqname and the indicated feature*/				
			
			if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0) {
				cstrand[0] = fieldsgff[n].strand[0];
				/*frame is only read for the first CDS in the gene (at -/+ strand)*/
				cframe[0] = cframe[1] = 'N';
				startframe = endframe = -1;
			
				for(ii=0;ii<n_site;ii++) cmat[ii] = (double)0;
				if(strcmp(subset_positions,"silent") == 0) for(ii=0;ii<n_site;ii++) cmatnc[ii] = (double)1;
				
				for(j=n;j<nrows;j++) {
					if(strcmp(fieldsgff[j].feature,"CDS") == 0) {
						if(strcmp(seqid,fieldsgff[j].seqname) == 0) {
							for(ii=fieldsgff[j].start;ii<=fieldsgff[j].end;ii++) cmat[ii] = (double)1;
							if(cstrand[0] == 32 || cstrand[0] == '.') cstrand[0] = fieldsgff[j].strand[0];
							else {
								if((cstrand[0] == '+' && fieldsgff[j].strand[0] == '-') || (cstrand[0] == '-' && fieldsgff[j].strand[0] == '+')) {
									cstrand[0] = '*';
									break;
								}
							}
							if(startframe == -1) {/*init*/
								cframe[0] = fieldsgff[j].frame[0];
								cframe[1] = fieldsgff[j].frame[0];
								startframe = fieldsgff[j].start;
								endframe = fieldsgff[j].end;
							}
							if(startframe > fieldsgff[j].start) {
								startframe = fieldsgff[j].start;
								cframe[0] = fieldsgff[j].frame[0];
							}
							if(endframe < fieldsgff[j].end) {
								endframe = fieldsgff[j].end;
								cframe[1] = fieldsgff[j].frame[0];
							}
						}
					}
					if(strcmp(subset_positions,"silent") == 0) {
						if(strcmp(fieldsgff[j].feature,"CDS") == 0) {
							if(strcmp(seqid,fieldsgff[j].seqname) == 0)
								for(ii=fieldsgff[j].start;ii<=fieldsgff[j].end;ii++) cmatnc[ii] = (double)0;
						}
					}
				}
				if(cstrand[0] == 32) cstrand[0] = '+';
				if(cframe[0] == 32) cframe[0] = '0';
				if(cframe[1] == 32) cframe[1] = '0';
			}
			else {
				if(strcmp(subset_positions,"noncoding") == 0) {
					for(ii=0;ii<n_site;ii++) cmat[ii] = (double)1;
					for(j=n;j<nrows;j++) {
						if(strcmp(fieldsgff[j].feature,"CDS") == 0) {
							if(strcmp(seqid,fieldsgff[j].seqname) == 0)
								for(ii=fieldsgff[j].start;ii<=fieldsgff[j].end;ii++) cmat[ii] = (double)0;
						}
					}
				}
				else {
					if(strcmp(subset_positions,"coding") == 0) {
						for(ii=0;ii<n_site;ii++) cmat[ii] = (double)0;
						for(j=n;j<nrows;j++) {
							if(strcmp(fieldsgff[j].feature,"CDS") == 0) {
								if(strcmp(seqid,fieldsgff[j].seqname) == 0)
									for(ii=fieldsgff[j].start;ii<=fieldsgff[j].end;ii++) cmat[ii] = (double)1;
							}
						}
					}
					else {
						for(ii=0;ii<n_site;ii++) cmat[ii] = (double)0;
						for(j=n;j<nrows;j++) {
							if(strcmp(fieldsgff[j].feature,subset_positions) == 0) {
								if(strcmp(seqid,fieldsgff[j].seqname) == 0)
									for(ii=fieldsgff[j].start;ii<=fieldsgff[j].end;ii++) cmat[ii] = (double)1;
							}
						}
					}
				}
			}
			/*count positions in matrix_sizepos*/
			/*first look at subset_positions if syn/nsyn/silent*/
			if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0) {
				/*search syn/nsyn postions: REMEMBER, LOOK AT 'CDS' FEATURE IN GFF FORMAT!*/
				/*look at strand and frame fields*/
				if(cstrand[0] == '*') {
					printf(" Error in GFF-format: annotation in %s is not considered. ",seqid);
					if(file_output) fprintf(file_output," Error in GFF-format: annotation in is %s not considered. ",seqid);
					continue;
				}
				else {
					/*check frame*/
					if(cstrand[0] == '-') {
						start = (long int)n_site-1;
						end = 0;
						k = -1;
						endtrp = startframe;
						startframe = endframe;
						cframe[0] = cframe[1];
					}
					else {
						start = 0;
						end = (long int)n_site-1;
						endtrp = endframe;
						k = 1;
					}
					/*start looking at coding region sequences*/
					for(ii=start;ii*k <= end*k;ii += k) {							
						while(ii*k <= end*k && cmat[ii] == (double)0) ii += k;
						if(ii*k > end*k) break;
						if(ii == startframe) {
							if(cframe[0] == '1') {
								cmat[ii] = (double)0;
								do{
									ii += k;
								}while(ii*k <= end*k && cmat[ii] == (double)0);
								cmat[ii] = (double)0;
								continue;
							}
							if(cframe[0] == '2') {
								cmat[ii] = (double)0;
								continue;
							}
						}
						/*function to read matrix with 3 * n_samp char, reverse-complementary if cstrand is '-', return 0 if more than 1 mutation in triplet or gaps/uncertainty*/
						if(tripletnsamp(cod3n,DNA_matr,cstrand[0],cmat,n_samp,(long int)n_site,end,ii,file_output) == 0) {
							cmat[ii] = (double)0;
							do{
								ii += k;
							}while(ii*k <= end*k && cmat[ii] == (double)0);
							if(ii*k > end*k) break;
							cmat[ii] = (double)0;
							do{
								ii += k;
							}while(ii*k <= end*k && cmat[ii] == (double)0);
							if(ii*k > end*k) break;
							cmat[ii] = (double)0;
							continue;
						}
						/*count positions in cmat*/
						stop = 0;
						cod3ft[0] = cod3ft[1] = cod3ft[2] = (double)0;
						for(j=0;j<n_samp;j++) {
							for(q=0;q<64;q++) {
								if(memcmp(tripletsN[q],cod3n+3*j,3) == 0) {
									aaseq[0] = genetic_code[q];
									break;
								}
							}
							if(aaseq[0] == '*') {
								ii2 = ii;
								do{
									ii2 += k;
								}while(cmat[ii2] == (double)0);
								do{
									ii2 += k;
								}while(cmat[ii2] == (double)0);
								if(endtrp != ii2) {
									printf(" SC %ld.",ii+1);
									if(file_output) fprintf(file_output," SC %ld.",ii+1);
								}
								stop = 1;
								break;
							}
							for(i=0;i<3;i++) {
								cod3f = (double)0;
								countpath = 0;
								for(nn=1;nn<=4;nn++) {
									cod3put[0] = cod3n[3*j+0];
									cod3put[1] = cod3n[3*j+1];
									cod3put[2] = cod3n[3*j+2];
									switch(nn) {
										case 1:
											if(cod3put[i]=='1') break;
											else cod3put[i] = '1';
											break;
										case 2:
											if(cod3put[i]=='2') break;
											else cod3put[i] = '2';
											break;
										case 3:
											if(cod3put[i]=='3') break;
											else cod3put[i] = '3';
											break;
										case 4:
											if(cod3put[i]=='4') break;
											else cod3put[i] = '4';
											break;
									}
									if(memcmp(cod3put,cod3n+3*j,3) == 0) continue;
									for(q=0;q<64;q++) {
										if(memcmp(tripletsN[q],cod3put,3) == 0) {
											aaput[0] = genetic_code[q];
											break;
										}
									}
									if(aaput[0] == '*') continue;
									else countpath += 1;
																		
									if(aaseq[0] == aaput[0]) {
										if(!(strcmp(subset_positions,"nonsynonymous") == 0)) cod3f += (double)1; /*syn*/
									}else 
										if(strcmp(subset_positions,"nonsynonymous") == 0) cod3f += (double)1; /*nsyn*/
								}
								if(countpath) cod3ft[i] += cod3f/(double)countpath;
							}
						}
						if(stop == 1) {
							cmat[ii] = (double)0;
							do{
								ii += k;
							}while(cmat[ii] == (double)0);
							if(ii*k > end*k) break;
							cmat[ii] = (double)0;
							do{
								ii += k;
							}while(cmat[ii] == (double)0);
							if(ii*k > end*k) break;
							cmat[ii] = (double)0;
							continue;
						}
						/*count biallelic sites in matrix_segrpos: REMEMBER overlapping coding regions are not well calculated!*/
						for(j=1;j<n_samp;j++) if(memcmp(cod3n+3*0,cod3n+3*j,3) != 0) break;
						if(j < n_samp) {
							for(i=0;i<3;i++) if(memcmp(cod3n+3*0+i,cod3n+3*j+i,1) != 0) break;
							for(q=0;q<64;q++) {
								if(memcmp(tripletsN[q],cod3n+3*0,3) == 0) {
									aaseq[0] = genetic_code[q];
									break;
								}
							}
							for(q=0;q<64;q++) {
								if(memcmp(tripletsN[q],cod3n+3*j,3) == 0) {
									aaput[0] = genetic_code[q];
									break;
								}
							}
							if(aaseq[0] == aaput[0]) {
								if((strcmp(subset_positions,"nonsynonymous") == 0)) matrix_segrpos[ii+i*k] = (double)0; /*reject syn*/
							}
							else
								if(!(strcmp(subset_positions,"nonsynonymous") == 0)) matrix_segrpos[ii+i*k] = (double)0; /*reject nsyn*/
						}
						/*weight cmat positions*/
						cmat[ii] = cod3ft[0]/((double)n_samp);
						do{
							ii += k;
						}while(cmat[ii] == (double)0);
						cmat[ii] = cod3ft[1]/((double)n_samp);
						do{
							ii += k;
						}while(cmat[ii] == (double)0);
						cmat[ii] = cod3ft[2]/((double)n_samp);
					}
					/*include cmat in matrix_sizepos*/
					for(ii=0;ii<(long int)n_site;ii++) {
						if(matrix_sizepos[ii] < cmat[ii]) matrix_sizepos[ii] = cmat[ii];
					}
					/*noncoding for silent calculated separatedly*/
					if(strcmp(subset_positions,"silent") == 0) {
						for(ii=0;ii<(long int)n_site;ii++) {
							if(cmatsil[ii] > cmatnc[ii]) cmatsil[ii] = cmatnc[ii];
						}
					}
				}
			}
			else {
				/*include cmat in matrix_sizepos for other options*/
				if(strcmp(subset_positions,"noncoding") == 0) {
					for(ii=0;ii<(long int)n_site;ii++) {
						if(matrix_sizepos[ii] > cmat[ii]) matrix_sizepos[ii] = cmat[ii];
					}
				}
				else {
					for(ii=0;ii<(long int)n_site;ii++) {
						if(matrix_sizepos[ii] < cmat[ii]) matrix_sizepos[ii] = cmat[ii];
					}
				}
			}
		}
		/*in silent: add noncoding to syn positions*/
		if(strcmp(subset_positions,"silent") == 0) {
			for(ii=0;ii<(long int)n_site;ii++) {
				if(cmatsil[ii] == (double)1) matrix_sizepos[ii] = cmatsil[ii];
			}
		}

		free(cmat);
		free(cmatnc);
		free(cmatsil);
		if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0) {
			free(cod3n);
			free(cod3put);
		}
	}
	free(fieldsgff);

	return 1; /*ok*/
}

/*function to read matrix with 3 * n_samp char, reverse-complementary if cstrand is '-', return 0 if more than 1 mutation in triplet or gaps/uncertainty*/
int tripletnsamp(char *cod3n,char *DNA_matr,char strand,double *cmat,int n_samp,long int n_site,long int end,long int ii,FILE *file_output)
{
	int i,j,k,x,z;
	long int ii2;
	char triplet1[3],triplet2[3];
	if(strand == '-') k = -1;
	else k = 1;
	
	/*make matrix with triplets*/
	ii2 = ii;
	for(i=0;i<3;i++) {
		for(j=0;j<n_samp;j++) {
			cod3n[3*j+i] = DNA_matr[n_site*j+ii2];
			/*read complementary when strand is '-'*/
			if(k==-1) {
				if(cod3n[3*j+i] == '1') cod3n[3*j+i] = '4';
				else {
					if(cod3n[3*j+i] == '2') cod3n[3*j+i] = '3';
					else {
						if(cod3n[3*j+i] == '3') cod3n[3*j+i] = '2';
						else if(cod3n[3*j+i] == '4') cod3n[3*j+i] = '1';
					}
				}
			}
		}
		if(i < 2) {
			do {
				ii2 += k;
			}while(ii2*k <= end*k && cmat[ii2] == (double)0);
		}
		if(ii2*k > end*k) 
			return 0;
	}
	
	/*check if gaps/uncertainty*/
	for(i=0;i<3*n_samp;i++) {
		if(!(cod3n[i] == '1' || cod3n[i] == '2' || cod3n[i] == '3' || cod3n[i] == '4' )) {
			printf(" GU %ld.",ii+1);
			if(file_output) fprintf(file_output," GU %ld.",ii+1);
			return 0;
		}
	}
	/*check if more than 1 mutation*/
	triplet1[0] = cod3n[3*0+0];
	triplet1[1] = cod3n[3*0+1];
	triplet1[2] = cod3n[3*0+2];
	triplet2[0] = '\0';
	for(i=1;i<n_samp;i++) {
		if(memcmp(triplet1,cod3n+(3*i),3*sizeof(char)) != 0) {
			if(triplet2[0] == '\0') {
				triplet2[0] = cod3n[3*i+0];
				triplet2[1] = cod3n[3*i+1];
				triplet2[2] = cod3n[3*i+2];
				if(memcmp(triplet2,triplet1,3*sizeof(char)) != 0) {
					z=0; 
					for(x=0;x<3;x++) {
						if(triplet1[x]!=triplet2[x]) 
							z++;
					}
					if(z>1) {
						if(file_output) fprintf(file_output," MM %ld.",ii+1);
						return 0;
					}
				}
			}
			else {
				if(memcmp(triplet2,cod3n+(3*i),3*sizeof(char)) != 0) {
					printf(" MM %ld.",ii+1);
					if(file_output) fprintf(file_output," MM %ld.",ii+1);
					return 0;
				}
			}
		}
	}
	
	return 1;
}



