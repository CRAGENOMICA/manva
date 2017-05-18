/*From Biometry*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*Sort abs(data) in qsort function*/ /*CHECK IF fabs WORKS!*/
int compabs(const void *i, const void *j)
{
    if(fabs(*(double *)i) < fabs(*(double *)j)) return -1;
    if(fabs(*(double *)i) > fabs(*(double *)j)) return  1;
    
    return 0;
}

/*calculate value Ts in Wilcoxon's Signed Rank Test (p. 443 in Biometry)*/
long int wilcoxonrst_ts(double *vector1,double *vector2,int size)
{
    double *vector12;
    int x;
    int compabs(const void *,const void *);
    long int sum,less;
    
    if((vector12 = (double *) malloc(size * sizeof(double))) == 0) {
        puts("Error in allocation memory. wilcoxonrst_ts.1");
        return -1;
    }
    
    for(x=0;x<size;x++) vector12[x] = vector1[x] - vector2[x];
    qsort(vector12,size,sizeof(double),compabs);
    sum = less = 0;
    for(x=0;x<size;x++) {
        if(vector12[x] > 0) sum += (long int)x;
        else if(vector12[x] < 0) less += (long int)x;
    }
    free(vector12);
    
    if(sum < less) return sum;
    else return less; 
    
    return 0;
}

double wilcoxontstonormal(int n,long int Ts)
{
    double ts;
    if(n > 50) {
        ts = ((double)Ts - ((double)n*((double)n+(double)1.0))/(double)4.0)/(double)sqrt(((double)n*((double)n+(double)0.5)*((double)n+(double)1.0))/(double)12.0);
        return ts;
    }
    
    return -1.;
}

