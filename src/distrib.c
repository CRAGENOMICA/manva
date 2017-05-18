/*NUMERICAL RECIPES IN C*/

#include <math.h>
#define PI 3.14159265359

#include <stdio.h>
#include <stdlib.h>

double bnldev(double pp, int n) 
{
	/*Based on Numerical Recipes in C, Press et al. 1992 and on
	Fishman 1979 J. American Statistical Association Vol. 74, No. 366, pp 418-423*/
	
	double ran1(void);	
	double p,np;
	int N;
	int nn;
	double r;
	double A,B,C,D,V,s;
	double m,mu;	
	static double *f=0;
	double poisso(double);
	static int max = 200;
	
	if(f == 0) {
		if((f=(double *)calloc(max+1,sizeof(double))) == 0) {
			printf("Error allocation memory");
			return -10000;
		}
		f[1] = (double)0;
		for(N=1;N<max;N++)
			f[N+1] = f[N] + (double)log((double)N);
	}
	if(n > max) {
		if((f=(double *)realloc(f,n*sizeof(double))) == 0) {
			printf("Error allocation memory");
			return -10000;
		}
		for(N=max;N<n;N++)
			f[N+1] = f[N] + (double)log((double)N);
		max = n;
	}
	
	if(pp > 0.5) p = (double)1.-pp;
	else p = pp;
	
	np = n * p;
	
	if(n==0) {
		puts("Error bindist");
		return (double)-10000.;
	}
	if(p==(double)0) {
		if(pp > 0.5) return (double)n;
		return (double)0;
	}
	
	if(n < 20) {
		/*Bernouilli Method*/
		nn = n;
		N=0;
		while(nn--) 
			if(ran1()<p) N++;		
	}
	else {
		if(np < (double)10) {
			/*Rejection Method: BI Algorithm*/
			s = (double)1- p;
			A = (double)1;
			B = p/s;
			C = ((double)n+(double)1)*B;
			D = A;
			N = 0;
			V = ran1()/(double)pow(s,(double)n);
			while(V > A) {
				N++;
				D *= (C/(double)N - B);
				A += D;
				if(N > n) break;
			}
		}
		else {
			/*Poisson method: BP Algorithm*/
			mu = n - (double)floor((double)(n*((double)1 - p)));
			if(n*((double)1-p) - (double)floor((double)(n*((double)1-p))) > p)
				mu = p*((double)floor((double)(n*((double)1-p))) + (double)1) / ((double)1-p);
			r = ((double)1/p - (double)1) * mu;
			s = (double)log((double)r);
			m = (double)floor((double)(r));
			do {
				do {
					N = (int)poisso(mu);
				}while((int)N > n);
				V = -(double)log((double)ran1());
			}while(V < (m-(double)(n - N))*s - f[(int)m+1] + f[(int)(n-N)+1]);
		}
	}
	if(pp > 0.5) N = n - N;
	return (double)N;
}

/*MODIFICATION TO LARGE N NUMBERS WITH N*P >= 1*/
double freqlongbnldev(double pp, long int n) 
{
	/*Based on Numerical Recipes in C, Press et al. 1992 and on
	Fishman 1979 J. American Statistical Association Vol. 74, No. 366, pp 418-423*/
	
	double ran1(void);	
	double p,np;
	long int N;
	long int nn;
	double A,B,C,D,V,s;
		
	if(pp > 0.5) p = (double)1.-pp;
	else p = pp;
	
	np = n * p;
	
	if(n==0) {
		puts("Error bindist");
		return (double)-10000.;
	}
	if(p==(double)0) {
		if(pp > 0.5) return (double)n;
		return (double)0;
	}
	
	if(n < 20) {
		/*Bernouilli Method*/
		nn = n;
		N=0;
		while(nn--) 
			if(ran1()<p) N++;		
	}
	else {
		/*Rejection Method: BI Algorithm*/
		s = (double)1- p;
		A = (double)1;
		B = p/s;
		C = ((double)n+(double)1)*B;
		D = A;
		N = 0;
		V = ran1()/(double)pow(s,(double)n);
		while(V > A) {
			N++;
			D *= (C/(double)N - B);
			A += D;
			if(N > n) break;
		}
	}
	if(pp > 0.5) N = n - N;
	return (double)N;
}

double factln_(unsigned long x)
{
	/*Based on Numerical Recipes in C, Press et al. 1992 and on
	Lanczos 1964 J. SIAM Numer. Anal. Ser. B, Vol. 1 pp 86-96.*/
	
	/*G(n+1) = n!*/
	/*do the log(n!)*/
	double gammln(double);
	static double factlog[120];
		
	/*
	if(x < 0) { 
		puts("Error factln");
		return (double)-10000.;
	}
	*/
	if(x == 0) return 0.;
	
	if(x < 120) {
		if(factlog[x] == (double)0) {
			factlog[x] = gammln((double)x+(double)1.0);
			return factlog[x];
		}
		else return factlog[x];
	}
	return (gammln((double)x+(double)1.0));
}

double poisso(double lambda) 
{
	/*Based on Atkinson 1979 Appl. Statist. 28: No. 1, pp, 29-35.*/
	double ran1(void);	
	double r,s;
	int N;
	
	double factln_(unsigned long);
	double alfa,beta,k,X;
	double rand1,rand2;
	static double c = (double)0.6;
	
	if(lambda < (double)0) {
		puts("Error in poisso()");
		return (double)-10000.;
	}
	
	if(lambda == (double)0) return (double)0;
	if(lambda <= (double)20) {
		r = (double)exp(-(double)lambda);
		N = (int)0;
		s = (double)1;
		do {
			s *= ran1();
			if(s >= r) N += (int)1;
			else break;
		}while(1);
	}
	else {
		beta = (double)PI * (double)1./(double)sqrt((double)3*(double)lambda);
		alfa = beta * (double)lambda;
		k = (double)log((double)c) - (double)lambda - (double)log((double)beta);
		do{
			rand1 = ran1();
			X = (alfa-(double)log(((double)1-(double)rand1)/(double)rand1))/beta;
			if(X >= (double)-0.5) {
				N = (int)(X + (double)0.5);
				rand2 = ran1();
				if(alfa - beta*X +
				  (double)log((double)(rand2/(((double)1+(double)exp((double)(alfa-beta*X)))*((double)1+(double)exp((double)(alfa-beta*X)))))) 
				  <= k + (double)N*(double)log((double)lambda) - (double)factln_((unsigned long)N)) 
					break;
			}
		}while(1);
	}
	return (double)N;
}



