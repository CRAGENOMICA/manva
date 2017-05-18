/*
 *  variances.c
 *  MuLoNeTests
 *
 *  Created by sonsins on Fri Mar 21 2003.
 *
 */

#include "MuLoNeTests.h"
#include "mhmlspmlnt.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double variances(double sumx2,double sumx,int sample)
{
    double var;
    double n = (double)sample;
    
    if(n < 3) return (double)-10000;
    
    var = (sumx2 - sumx*sumx/(double)n)/((double)n-(double)1.);
    
    return var;
}

double varianceit(double sumx2,double sumx,long int sample)
{
    double var;
    double n = (double)sample;
    
    if(n < 3) return (double)-10000;

    var = (sumx2 - sumx*sumx/(double)n)/((double)n-(double)1.);
    
    return var;
}

