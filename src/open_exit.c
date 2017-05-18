/*
 *  open_exit.c
 *  MuLoNeTests
 *
 *  Created by sebas on Mon Feb 24 2003.
 *
 */

#include "MuLoNeTests.h"
#include <stdio.h>
#include <stdlib.h>

int open_exit()
{
    char l;
    
    puts(" Are you sure you want to exit the program (y/n)? ");
    
    do l = getchar();
    while(l !='y' && l !='n' && l !='Y' && l !='N');
    
    if(l == 'y' || l == 'Y') return 0;
    else if(l == 'n' || l == 'N') return 1;
        
    return 1;
}

