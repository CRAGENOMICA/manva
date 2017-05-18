/*
 *  histograms.c
 *  MuLoNeTests
 *
 *  Created by sonsins on Fri Mar 07 2003.

 *
 */

#include "MuLoNeTests.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct sorted_kolmsm {
	double sortedv;
	int pop;
};
struct kms_test {
	double f_a;
	double f_b;
	double abs_d;
};

int start_histogram(double *vectorobs,double *vectorsim,int totalobs,long int totalsim,char *nameobs,char *namesim,FILE *file_output)
{
    char k[1],l[1];
    
    int nsegments;
    int nblocks;
    int limitx;
    double minx;
    double maxx;
    
    int do_histogram(double *,double *,int,long int,char *,char *,int *,int *,int *,double *,double *,FILE *);
    
    nsegments = 30;
    nblocks = 30;
    limitx = 0;
    minx = (double)0.;
    maxx = (double)0.;
    
    #if COMMAND_LINE    
    if(namesim == 0)
        printf("\n\n     Histogram for %s.\n\n",nameobs);
    else if(nameobs == 0)
            printf("\n\n     Histogram for %s.\n\n",namesim);
        else printf("\n\n     Histogram for %s - %s.\n\n",nameobs,namesim);
    
    printf("CHOOSE:\n");
    printf(" 0 - Do Histogram with automatic options.\n");
    printf(" 1 - Advanced options.\n");
    printf(" 2 - Go to the next histogram.\n");
    printf(" 3 - Back to previous menu.\n");
	/*
    if(file_output) {
        fprintf(file_output,"\n 0 - Do Histogram with automatic options.\n");
        fprintf(file_output," 1 - Advanced options.\n");
        fprintf(file_output," 2 - Go to the next histogram.\n");
        fprintf(file_output," 3 - Back to previous menu.\n");
    }
	*/
    #endif
    
    do *k = getchar();
    while(*k<'0' || *k>'3');
    
    if(file_output) fprintf(file_output,"\nCHOSEN: %c\n",*k);

    switch(*k) {
        case '0':
            if(!(do_histogram(vectorobs,vectorsim,totalobs,totalsim,
                nameobs,namesim,&nsegments,&nblocks,&limitx,&minx,&maxx,file_output)))
                return 0;
            break;
        case '1':
            /*Ask all options:*/
            do {
                printf("\n  Please indicate the number of segments you want to divide the data (30 is reccomended):");
                scanf(" %d",&nsegments);
            }while(nsegments < 1);
            do {
                printf("\n  Please indicate the number of partitions in the Y-axis (30 is reccomended):");
                scanf(" %d",&nblocks);
            }while(nblocks < 1);
            printf("\n\n  Do you want to indicate the maximum and minimum limits of the X-axis (y/n)? ");
            do *l = getchar();
            while(*l!='y' && *l!='n' && *l!='Y' && *l!='N');            
            if(*l == 'y' || *l == 'Y') {
                limitx = 1;
                printf("\n  Please indicate the minimum value in the X-axis:");
                scanf(" %lg",&minx);
                do {
                    printf("\n  Please indicate the maximum value in the X-axis:");
                    scanf(" %lg",&maxx);
                    if(maxx <= minx) printf("\n The maximum value must be higher than the minimum value.");
                } while(maxx <= minx);
            }
            else limitx = 0;            
            if(!(do_histogram(vectorobs,vectorsim,totalobs,totalsim,nameobs,namesim,
                &nsegments,&nblocks,&limitx,&minx,&maxx,file_output)))
                return 0;
            break;
        case '2':
            break;
        case '3':
            return 2;
            break;
    }
    return 1;
}

int do_histogram(double *vectorobs,double *vectorsim,int totalobs,long int totalsim,char *nameobs,char *namesim,int *nsegments,int *nblocks,int *limitx,double *minx,double *maxx,FILE *file_output)
{
    char **plot;
    long int x,y,z;
    int j,k,l,m;
    int outsiderxmin = 0;
    int outsiderxmax = 0;
    double intervalx,intervalyobs,intervalysim;
    int limitplotx,limitploty;
    long int *segment;
	struct sorted_kolmsm *sorted_tot;
	struct kms_test *kstest;
	void sort_struct(struct sorted_kolmsm *, long int, long int);
	int nsections;
    double Dks,Dalfa,alfa;
	double add;

    if(!((vectorobs && totalobs) || (vectorsim && totalsim) || (*nsegments) || (*nblocks))) return 0;
    if((*limitx) == 0) {
        /* Assignate the max and minimum limits.*/
        if(totalobs && vectorobs) {
            (*minx) = (*maxx) = vectorobs[0];
            for(x=1;x<(long int)totalobs;x++) {
                if(vectorobs[x] > (*maxx)) (*maxx) = vectorobs[x];
                if(vectorobs[x] < (*minx)) (*minx) = vectorobs[x];
            }
        }
        if(totalsim && vectorsim) {
            if(!(totalobs && vectorobs)) {
                (*minx) = (*maxx) = vectorsim[0];
                y = 1;
            }
            else y = 0;
            for(x=y;x<(long int)totalsim;x++) {
                if(vectorsim[x] > (*maxx)) (*maxx) = vectorsim[x];
                if(vectorsim[x] < (*minx)) (*minx) = vectorsim[x];
            }
        }
    }
    else {
        /* check for outsiders from the assigned limits.*/
        if(totalobs && vectorobs) {
            x = 0;
            while((x < totalobs) && (outsiderxmax == 0 || outsiderxmin == 0)) {
                if(vectorobs[x] > (*maxx)) outsiderxmax = 1;
                if(vectorobs[x] < (*minx)) outsiderxmin = 1;
                x++;
            }
        }
        if(totalsim && vectorsim) {
            x = 0;
            while((x < totalsim) && (outsiderxmax == 0 || outsiderxmin == 0)) {
                if(vectorsim[x] > (*maxx)) outsiderxmax = 1;
                if(vectorsim[x] < (*minx)) outsiderxmin = 1;
                x++;
            }
        }
    }
	add = ((*maxx)-(*minx))/(double)10000;
	(*maxx) += add;
	(*minx) -= add;
    /*in case all numbers have the same value*/
    if(*minx == *maxx) {
        *minx -= 3.;
        *maxx += 3.;
    }
    
    /*Fill the segments.*/
    intervalx = ((*maxx)-(*minx))/(double)(*nsegments);
    if((segment = (long int *) calloc(2*(*nsegments)+4,sizeof(long int))) == 0) {
        printf("\nError: memory not reallocated. do_histogram.1 \n");
        return (0);
    }
    if(totalobs && vectorobs) {
        for(x=0;x<(long int)totalobs;x++) {
            if(outsiderxmin) {
                if(vectorobs[x] < (*minx)) {
                    segment[0] += (long int)1;
                    continue;
                }
            }
            for(z=1,y=2;y<(long int)(2*(*nsegments)+2);y += 2,z++) {            
                if(vectorobs[x] <= (double)z*intervalx + (*minx)) {
                    segment[y] += (long int)1;
                    break;
                }
            }
            if(outsiderxmax) {
                if(vectorobs[x] > (*maxx)) {
                    segment[2*(*nsegments)+2] += (long int)1;
                    continue;
                }
            }
        }
    }
    if(totalsim && vectorsim) {
        for(x=0;x<(long int)totalsim;x++) {
            if(outsiderxmin) {
                if(vectorsim[x] < (*minx)) {
                    segment[1] += (long int)1;
                    continue;
                }
            }
            for(z=1,y=3;y<(long int)(2*(*nsegments)+2);y += 2,z++) {            
                if(vectorsim[x] <= (double)z*intervalx + (*minx)) {
                    segment[y] += (long int)1;
                    break;
                }
            }
            if(outsiderxmax) {
                if(vectorsim[x] > (*maxx)) {
                    segment[2*(*nsegments)+3] += (long int)1;
                    continue;
                }
            }
        }
    }

    /*count the intervals y.*/
    y = 0;
    for(x=0;x<(long int)(2*(*nsegments)+4);x += 2) {
        if(y < segment[x]) y = segment[x];
    }
    intervalyobs = /*ceil*/((double)y/(double)(*nblocks));
    y = 0;
    for(x=1;x<(long int)(2*(*nsegments)+4);x += 2) {
        if(y < segment[x]) y = segment[x];
    }
    intervalysim = /*ceil*/((double)y/(double)(*nblocks));
    
    /*Do the plot.*/
    /*axis.*/
    limitplotx = 4 + 2*outsiderxmin + 3*(*nsegments) + 2*outsiderxmax + 4;
    limitploty = 1 + (*nblocks);
    
    if((plot = (char **) malloc((unsigned long)limitploty*sizeof(char *))) == 0) {
        printf("\nError: memory not reallocated. do_histogram.2 \n");
        return (0);
    }
    for(j=0;j<limitploty;j++) {
        if((plot[j] = (char *) malloc((unsigned long)limitplotx*sizeof(char))) == 0) {
            printf("\nError: memory not reallocated. do_histogram.3 \n");
            return (0);
        }
        memset(plot[j],32,(unsigned long)limitplotx);
    }
    
    j = 0;
    plot[limitploty-1][j++] = '|';
    plot[limitploty-1][j++] = ' ';
    plot[limitploty-1][j++] = '_';
    plot[limitploty-1][j++] = '_';
    if(outsiderxmin) plot[limitploty-1][j++] = '/';
    if(outsiderxmin) plot[limitploty-1][j++] = '/';

    k = j;
    for(j=k;j<k+(*nsegments)*3;) {
        plot[limitploty-1][j++] = ' ';
        plot[limitploty-1][j++] = '_';
        plot[limitploty-1][j++] = '_';
    }

    if(outsiderxmax) plot[limitploty-1][j++] = '/';
    if(outsiderxmax) plot[limitploty-1][j++] = '/';
    plot[limitploty-1][j++] = ' ';
    plot[limitploty-1][j++] = '_';
    plot[limitploty-1][j++] = '_';

    plot[limitploty-1][j++] = '\0';
    for(l=0;l<limitploty-1;l++) {
        plot[l][0] = '|';
        plot[l][limitplotx-1] = '\0';
    }
    plot[limitploty-1][limitplotx-1] = '\0';
    
    /*histogram. observed is character x (before 165), simulated 79.*/

    j = 0;
    m = 2;
    k = (int)floor(((double)segment[j] / intervalyobs));
    for(l=1;l<k+1;l++) plot[limitploty-1-l][m] = 'x'/*165*/;
    j++;
    m++;
    k = (int)floor(((double)segment[j] / intervalysim));
    for(l=1;l<k+1;l++) plot[limitploty-1-l][m] = 79;
    j++;
    if(outsiderxmin) m += 4;
    else m +=2;

    for(j=2;j<2*(*nsegments)+2;) {
        k = (int)floor((int)((double)segment[j] / intervalyobs));
        for(l=1;l<k+1;l++) plot[limitploty-1-l][m] = 'x'/*165*/;
        j++;
        m++;
        k = (int)floor(((double)segment[j] / intervalysim));
        for(l=1;l<k+1;l++) plot[limitploty-1-l][m] = 79;
        j++;
        m += 2;
    }

    if(outsiderxmax) m += 2;
    k = (int)floor(((double)segment[j] / intervalyobs));
    for(l=1;l<k+1;l++) plot[limitploty-1-l][m] = 'x'/*165*/;
    j++;
    m++;
    k = (int)floor(((double)segment[j] / intervalysim));
    for(l=1;l<k+1;l++) plot[limitploty-1-l][m] = 79;
    
    /*Print plot. Titles, legends and plot. then table/s.*/
    printf("\n\n HISTOGRAM showing the distribution of the values: ");
    if(file_output) fputs("\n\n HISTOGRAM showing the distribution of the values: ",file_output);
    
    /*names.*/
    if(totalobs && vectorobs) {
        printf("%s",nameobs);
        if(file_output) fprintf(file_output,"%s",nameobs);    
        if(totalsim && vectorsim) {
            printf(" and ");
            if(file_output) fputs(" and ",file_output);
        }
    }
    if(totalsim && vectorsim) {
        printf("%s",namesim);
        if(file_output) fprintf(file_output,"%s",namesim);
    }
    printf("\n\n");
    if(file_output) fputs("\n\n",file_output);
    printf(" Only included polymorphisms that are not shared with the outgroup lines (if included).\n\n");
    if(file_output) fputs(" Only included polymorphisms that are not shared with the outgroup lines (if included).\n\n",file_output);
    
    /*plot.*/
    for(j=0;j<limitploty;j++) printf("%s\n",plot[j]);
    if(file_output) for(j=0;j<limitploty;j++) fprintf(file_output,"%s\n",plot[j]);
    
    /*legends: min,max X. interval X. interval Yobs and Ysim. outsiders.*/
    
    printf("  Axis Y: Frequency for each interval. Axis-X: Value of the statistic.");
    if(file_output) fprintf(file_output,"  Axis Y: Frequency for each interval. Axis-X: Value of the plotted statistic.");
    printf("\n\n");
    if(file_output) fputs("\n\n",file_output);
    if(outsiderxmin) {
        printf("\n Number of values on the histogram lower than %.3g (minimum limit) are indicated on the left of the plot before '//'.",(*minx));
        if(file_output) fprintf(file_output,"\n Number of values on the histogram lower than %.3g (minimum limit) are indicated on the left of the plot before '//'.",(*minx));
    }
    if(outsiderxmax) {
        printf("\n Number of values on the histogram higher than %.3g (maximum limit) are indicated on the right of the plot after '//'.",(*maxx));
        if(file_output) fprintf(file_output,"\n Number of values on the histogram higher than %.3g (maximum limit) are indicated on the right of the plot after '//'.",(*maxx));
    }
    if(outsiderxmin == 0 && outsiderxmax == 0)
        printf("\n Histogram showing a column of data each interval of %.3g and starting from %.3G as a minimum and %.3g as a maximum.",intervalx,(*minx)-intervalx,(*maxx)+intervalx);
    else printf("\n Histogram showing a column of data each interval of %.3g.",intervalx);
    if(totalobs && vectorobs)
        printf("\n Each block in observed values indicates %.2f cases with values in the specific interval.",intervalyobs);
    if(totalsim && vectorsim)
        printf("\n Each block in simulated values indicates %.2f cases with values in the specific interval.\n",intervalysim);
    if(file_output) {
        if(outsiderxmin == 0 && outsiderxmax == 0)
            fprintf(file_output,"\n Histogram showing a column of data each interval of %.3g and starting from %.3g as a minimum and %.3g as a maximum.",intervalx,(*minx)-intervalx,(*maxx)+intervalx);
        else fprintf(file_output,"\n Histogram showing a column of data each interval of %.3g.",intervalx);
        if(totalobs && vectorobs)
            fprintf(file_output,"\n Each block in observed values indicates %.2f cases with values in the specific interval.",intervalyobs);
        if(totalsim && vectorsim)
            fprintf(file_output,"\n Each block in simulated values indicates %.2f cases with values in the specific interval.\n",intervalysim);
    }
    
    /*legend obs and sim. 165 and 79.*/
    if(totalobs && vectorobs) {
        printf("\n Observed values are represented with %c.",'x'/*165*/);
        if(file_output) fprintf(file_output,"\n Observed values are represented with %c.",'x'/*165*/);
    }
    if(totalsim && vectorsim) {
        printf("\n Expected values are represented with %c.\n",79);
        if(file_output) fprintf(file_output,"\n Simulated values are represented with %c.\n",79);
    }

    /*vectors.(Table)*/
    printf(" \n In the Table is represented the intermediate value of the interval\n");
    if(file_output) fprintf(file_output," \n In the Table is represented the intermediate value of the interval\n");

    printf(" \n Table. (same order than in the histogram)\n\n");
    if(file_output) fprintf(file_output," \n Table. (same order than in the histogram)\n\n");
    
    printf(" Interval:\t");
    if(file_output) fputs(" Interval:\t",file_output);
    if(outsiderxmin) {
        printf("<%g\t",(*minx));
        if(file_output) fprintf(file_output,"<%g\t",(*minx));
    }
    else {
        printf("%g\t",(*minx)-intervalx/2.);
        if(file_output) fprintf(file_output,"%g\t",(*minx)-intervalx/2.);
    }
    for(j=0;j<(*nsegments);j ++) {
        printf("%g\t",(*minx)+((double)j*intervalx)+intervalx/2.);
        if(file_output) fprintf(file_output,"%g\t",(*minx)+((double)j*intervalx)+intervalx/2.);
    }
    if(outsiderxmax) {
        printf(">%g\t",(*maxx));
        if(file_output) fprintf(file_output,">%g\t",(*maxx));
    }
    else {
        printf("%g\t",(*maxx)+intervalx/2.);
        if(file_output) fprintf(file_output,"%g\t",(*maxx)+intervalx/2.);
    }
        
    if(totalobs && vectorobs) {
        printf("\n");
        if(file_output) fputs("\n",file_output);
       printf(" Observed:\t");
        if(file_output) fputs(" Observed:\t",file_output);
        for(j=0;j<2*(*nsegments)+4;j += 2) {
            printf("%ld\t",segment[j]);
            if(file_output) fprintf(file_output,"%ld\t",segment[j]);
        }
    }
    if(totalsim && vectorsim) {
        printf("\n");
        if(file_output) fputs("\n",file_output);
        printf(" Expected:\t");
        if(file_output) fputs(" Expected:\t",file_output);
        for(j=1;j<2*(*nsegments)+4;j += 2) {
            printf("%ld\t",segment[j]);
            if(file_output) fprintf(file_output,"%ld\t",segment[j]);
        }
    }
    printf("\n");
    if(file_output) fputs("\n",file_output);

    if(totalobs > 25 && totalsim > 40 /*&& *nsegments==0*/)/*NOT YET COMPLETELY ACTIVE*/
    {
        /*Do test Kolmogorov-Smirnov.*/
		if((sorted_tot = (struct sorted_kolmsm *) calloc(totalobs+totalsim,sizeof(struct sorted_kolmsm))) == 0) {
			printf("\nError: memory not reallocated. do_histogram.10 \n");
			return (0);
		}
		/*introduce data*/
		for(x=0;x<totalobs;x++) {
			sorted_tot[x].sortedv = (double)vectorobs[x];
			sorted_tot[x].pop = (int)1;
		}
		for(x=totalobs;x<totalsim+totalobs;x++) {
			sorted_tot[x].sortedv = (double)vectorsim[x-totalobs];
			sorted_tot[x].pop = (int)2;
		}
		/*sort struct*/
		sort_struct(sorted_tot,0,totalsim+totalobs-1);
		/*number of sections*/
		nsections = *nsegments;
		if(outsiderxmin)
			nsections += 1;
		if(outsiderxmax)
			nsections += 1;
		/*init ks_test struct*/
		if((kstest = (struct kms_test *) calloc(nsections,sizeof(struct kms_test))) == 0) {
			printf("\nError: memory not reallocated. do_histogram.10 \n");
			return (0);
		}
		/*calculate f_a,f_b and abs_d*/
		x=y=z=0;
		for(j=0;j<nsections;j++) {
            if(outsiderxmin && j == 0) {
                while(sorted_tot[x].sortedv < (*minx)) {
					if(sorted_tot[x].pop == 1) y++;
					if(sorted_tot[x].pop == 2) z++;
					x++;
				}
			}
            else {
				if(outsiderxmax && j == nsections-1) {
					while(sorted_tot[x].sortedv > (*maxx)) {
						if(sorted_tot[x].pop == 1) y++;
						if(sorted_tot[x].pop == 2) z++;
						x++;
					}
				}
				else {
					while(x < totalobs+totalsim && sorted_tot[x].sortedv <= (*minx + (double)intervalx * (double)(j+1))) {
						if(sorted_tot[x].pop == 1) y++;
						if(sorted_tot[x].pop == 2) z++;
						x++;
					}
				}
			}
			kstest[j].f_a = (double)y/(double)totalobs;
			kstest[j].f_b = (double)z/(double)totalsim;
			kstest[j].abs_d = (double)fabs(kstest[j].f_a - kstest[j].f_b);
		}
		/* obtain D*/
		Dks = kstest[0].abs_d;
		for(j=1;j<nsections;j++) {
			if(Dks < kstest[j].abs_d) 
				Dks = kstest[j].abs_d;
		}
		/*critical values*/
		alfa = (double)0.01;
		Dalfa = (double)sqrt(-(double)log((double)0.5*alfa)/(double)2.);
        Dalfa *= (double)sqrt((double)(totalobs+totalsim)/(double)(totalobs*totalsim));
		if(Dks > Dalfa) {
			/*printf prob < 0.01*/
			printf("\n Kolmogorov-Smirnov Test: D = %g; D(P = 0.01) = %g; P < 0.01 SIGNIFICANT ",Dks,Dalfa);
			if(file_output) fprintf(file_output,"\n Kolmogorov-Smirnov Test: D = %g; D(P = 0.01) = %g; P < 0.01 SIGNIFICANT ",Dks, Dalfa);
		}
		else {
			alfa = (double)0.05;
			Dalfa = (double)sqrt(-(double)log((double)0.5*alfa)/(double)2.);
			Dalfa *= (double)sqrt((double)(totalobs+totalsim)/(double)(totalobs*totalsim));
			if(Dks > Dalfa) {
				/*printf prob < 0.05*/
				printf("\n Kolmogorov-Smirnov Test: D = %g; D(P = 0.05) = %g; P < 0.05 SIGNIFICANT ",Dks,Dalfa);
				if(file_output) fprintf(file_output,"\n Kolmogorov-Smirnov Test: D = %g; D(P = 0.05) = %g; P < 0.05 SIGNIFICANT ",Dks, Dalfa);
			}
			else {
				/*printf prob Non-significant*/
				printf("\n Kolmogorov-Smirnov Test: D = %g; D(P = 0.05) = %g; NS",Dks,Dalfa);
				if(file_output) fprintf(file_output,"\n Kolmogorov-Smirnov Test: D = %g; D(P = 0.05) = %g; NS",Dks,Dalfa);
			}
		}
		if(outsiderxmin || outsiderxmax) {
			printf("\n Results of this test migth be strongly affected by the tail/s of the histogram distribution!!\n");
			if(file_output) fprintf(file_output,"\n Results of this test migth be strongly affected by the tail/s of the histogram distribution!!\n");
		}
		free(sorted_tot);
		free(kstest);
    }

    free(segment);
    for(j=0;j<limitploty;j++) free(plot[j]);
    free(plot);
    return 1;
}

/*From reference manual in C: Herbert Schildt; 2000*/
void sort_struct(struct sorted_kolmsm item[], long int left, long int right)
{
	register long int i,j;
	double x;
	struct sorted_kolmsm temp;
	
	i = left;
	j = right;
	x = item[(left+right)/2].sortedv;
	
	do {
		while((item[i].sortedv < x) && (j < right)) i++;
		while((x < item[j].sortedv) && (j > left )) j--;
		if(i<=j) {
			temp = item[i];
			item[i] = item[j];
			item[j] = temp;
			i++;
			j--;
		}
	} while(i<=j);	
	if(left < j) sort_struct(item,left,j);
	if(i < right) sort_struct(item,i,right);
}

