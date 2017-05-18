/*
 *  comblocisorted.c
 *  thetamlS
 *
 *  Created by sebas on Tue May 04 2004.
 *
 */

/*M2*/
int findml2sorted(int numloci2m1,int *sortedlociml,int *nloci2ml,int nloci,int steps,double *thetaM2,double *mlM2,double **ml)
{
    int x,y;
    double th1=0;
	double th2=0;
    double currml1=0;
	double currml2=0;
	double ml1=0;
	double ml2=0;
    
    for(x=0;x<steps;x++) {
        currml1 = currml2 = (double)0.;
        for(y=0;y<numloci2m1;y++) currml1 += ml[sortedlociml[y]][x];
        for(y=numloci2m1;y<nloci;y++) currml2 += ml[sortedlociml[y]][x];
        if(currml1 > ml1 || x==0) {
            ml1 = currml1;
            th1 = (double)x;
        }
        if(currml2 > ml2 || x==0) {
            ml2 = currml2;
            th2 = (double)x;
        }
    }
    if((ml1 + ml2 > mlM2[0] + mlM2[1]) || nloci2ml[0] == -1) {
        mlM2[0] = ml1;
        mlM2[1] = ml2;
        thetaM2[0] = th1;
        thetaM2[1] = th2;
        *nloci2ml = numloci2m1;
    }
    return 0;
}
int findml2coalsorted(int numloci2m1,int *sortedlociml,int nloci,int steps,double *mlM2,double **ml)
{
    int x,y;
    double currml1=0;
	double currml2=0;
	double ml1=0;
	double ml2=0;
    
    for(x=0;x<steps;x++) {
        currml1 = currml2 = (double)0.;
        for(y=0;y<numloci2m1;y++) currml1 += ml[sortedlociml[y]][x];
        for(y=numloci2m1;y<nloci;y++) currml2 += ml[sortedlociml[y]][x];
        if(currml1 > ml1 || x==0) {
            ml1 = currml1;
        }
        if(currml2 > ml2 || x==0) {
            ml2 = currml2;
        }
    }
    if((ml1 + ml2 > mlM2[0] + mlM2[1])) {
        mlM2[0] = ml1;
        mlM2[1] = ml2;
    }
    return 0;
}
/*M3*/
int findml3sorted(int numloci3m1,int numloci3m2,int *sortedlociml,int *nloci3ml1,int *nloci3ml2,int nloci,int steps,double *thetaM3,double *mlM3,double **ml)
{
    int x,y;
    double th1=0;
	double th2=0;
	double th3=0;
    double currml1=0;
	double currml2=0;
	double currml3=0;
	double ml1=0;
	double ml2=0;
	double ml3=0;
    
    for(x=0;x<steps;x++) {
        currml1 = currml2 = currml3 = (double)0.;
        for(y=0;y<numloci3m1;y++) currml1 += ml[sortedlociml[y]][x];
        for(y=numloci3m1;y<numloci3m1+numloci3m2;y++) currml2 += ml[sortedlociml[y]][x];
        for(y=numloci3m1+numloci3m2;y<nloci;y++) currml3 += ml[sortedlociml[y]][x];
        if(currml1 > ml1 || x==0) {
            ml1 = currml1;
            th1 = (double)x;
        }
        if(currml2 > ml2 || x==0) {
            ml2 = currml2;
            th2 = (double)x;
        }
        if(currml3 > ml3 || x==0) {
            ml3 = currml3;
            th3 = (double)x;
        }
    }
    if((ml1 + ml2 + ml3 > mlM3[0] + mlM3[1] + mlM3[2]) || nloci3ml1[0] == -1){
        mlM3[0] = ml1;
        mlM3[1] = ml2;
        mlM3[2] = ml3;
        thetaM3[0] = th1;
        thetaM3[1] = th2;
        thetaM3[2] = th3;
        *nloci3ml1 = numloci3m1;
        *nloci3ml2 = numloci3m2;
    }
    return 0;
}
int findml3coalsorted(int numloci3m1,int numloci3m2,int *sortedlociml,int nloci,int steps,double *mlM3,double **ml)
{
    int x,y;
    double currml1=0;
	double currml2=0;
	double currml3=0;
	double ml1=0;
	double ml2=0;
	double ml3=0;
    
    for(x=0;x<steps;x++) {
        currml1 = currml2 = currml3 = (double)0.;
        for(y=0;y<numloci3m1;y++) currml1 += ml[sortedlociml[y]][x];
        for(y=numloci3m1;y<numloci3m1+numloci3m2;y++) currml2 += ml[sortedlociml[y]][x];
        for(y=numloci3m1+numloci3m2;y<nloci;y++) currml3 += ml[sortedlociml[y]][x];
        if(currml1 > ml1 || x==0) {
            ml1 = currml1;
        }
        if(currml2 > ml2 || x==0) {
            ml2 = currml2;
        }
        if(currml3 > ml3 || x==0) {
            ml3 = currml3;
        }
    }
    if((ml1 + ml2 + ml3 > mlM3[0] + mlM3[1] + mlM3[2])){
        mlM3[0] = ml1;
        mlM3[1] = ml2;
        mlM3[2] = ml3;
    }
    return 0;
}
int findml3sortedmh(int numloci3m1,int numloci3m2,int *sortedlociml,int *relnloci3ml1,int *relnloci3ml2,int *nloci3ml1,int *nloci3ml2,int nloci,int steps,double *relthetaM3,double *relmlM3,double *thetaM3,double *mlM3,double **ml,double logranloci)
{
    int x,y;
    double th1=0;
	double th2=0;
	double th3=0;
    double currml1=0;
	double currml2=0;
	double currml3=0;
	double ml1=0;
	double ml2=0;
	double ml3=0;
    
    for(x=0;x<steps;x++) {
        currml1 = currml2 = currml3 = (double)0.;
        for(y=0;y<numloci3m1;y++) currml1 += ml[sortedlociml[y]][x];
        for(y=numloci3m1;y<numloci3m1+numloci3m2;y++) currml2 += ml[sortedlociml[y]][x];
        for(y=numloci3m1+numloci3m2;y<nloci;y++) currml3 += ml[sortedlociml[y]][x];
        if(currml1 > ml1 || x==0) {
            ml1 = currml1;
            th1 = (double)x;
        }
        if(currml2 > ml2 || x==0) {
            ml2 = currml2;
            th2 = (double)x;
        }
        if(currml3 > ml3 || x==0) {
            ml3 = currml3;
            th3 = (double)x;
        }
    }
    /*acceptation of ml using MH*/
    if(((ml1 + ml2 + ml3) - (relmlM3[0] + relmlM3[1] + relmlM3[2]) >= logranloci) || nloci3ml1[0] == -1){
        relmlM3[0] = ml1;
        relmlM3[1] = ml2;
        relmlM3[2] = ml3;
        relthetaM3[0] = th1;
        relthetaM3[1] = th2;
        relthetaM3[2] = th3;
        *relnloci3ml1 = numloci3m1;
        *relnloci3ml2 = numloci3m2;
    }
    /*find best ml*/
    if((ml1 + ml2 + ml3 > mlM3[0] + mlM3[1] + mlM3[2]) || nloci3ml1[0] == -1){
        mlM3[0] = ml1;
        mlM3[1] = ml2;
        mlM3[2] = ml3;
        thetaM3[0] = th1;
        thetaM3[1] = th2;
        thetaM3[2] = th3;
        *nloci3ml1 = numloci3m1;
        *nloci3ml2 = numloci3m2;
    }
    return 0;
}
int findml3coalsortedmh(int numloci3m1,int numloci3m2,int *sortedlociml,int nloci,int steps,double *relmlM3,double *mlM3,double **ml,double logranloci)
{
    int x,y;
    double currml1=0;
	double currml2=0;
	double currml3=0;
	double ml1=0;
	double ml2=0;
	double ml3=0;
    
    for(x=0;x<steps;x++) {
        currml1 = currml2 = currml3 = (double)0.;
        for(y=0;y<numloci3m1;y++) currml1 += ml[sortedlociml[y]][x];
        for(y=numloci3m1;y<numloci3m1+numloci3m2;y++) currml2 += ml[sortedlociml[y]][x];
        for(y=numloci3m1+numloci3m2;y<nloci;y++) currml3 += ml[sortedlociml[y]][x];
        if(currml1 > ml1 || x==0) {
            ml1 = currml1;
        }
        if(currml2 > ml2 || x==0) {
            ml2 = currml2;
        }
        if(currml3 > ml3 || x==0) {
            ml3 = currml3;
        }
    }
    /*acceptation of ml using MH*/
    if(((ml1 + ml2 + ml3) - (relmlM3[0] + relmlM3[1] + relmlM3[2]) >= logranloci)){
        relmlM3[0] = ml1;
        relmlM3[1] = ml2;
        relmlM3[2] = ml3;
    }
    /*find best ml*/
    if((ml1 + ml2 + ml3 > mlM3[0] + mlM3[1] + mlM3[2])){
        mlM3[0] = ml1;
        mlM3[1] = ml2;
        mlM3[2] = ml3;
    }
    return 0;
}
/*M4*/
int findml4sorted(int numloci4m1,int numloci4m2,int numloci4m3,int *sortedlociml,int *nloci4ml1,int *nloci4ml2,int *nloci4ml3,int nloci,int steps,double *thetaM4,double *mlM4,double **ml)
{
    int x,y;
    double th1=0;
	double th2=0;
	double th3=0;
	double th4=0;
    double currml1=0;
	double currml2=0;
	double currml3=0;
	double currml4=0;
	double ml1=0;
	double ml2=0;
	double ml3=0;
	double ml4=0;
    
    for(x=0;x<steps;x++) {
        currml1 = currml2 = currml3 = currml4 = (double)0.;
        for(y=0;y<numloci4m1;y++) currml1 += ml[sortedlociml[y]][x];
        for(y=numloci4m1;y<numloci4m1+numloci4m2;y++) currml2 += ml[sortedlociml[y]][x];
        for(y=numloci4m1+numloci4m2;y<numloci4m1+numloci4m2+numloci4m3;y++) currml3 += ml[sortedlociml[y]][x];
        for(y=numloci4m1+numloci4m2+numloci4m3;y<nloci;y++) currml4 += ml[sortedlociml[y]][x];
        if(currml1 > ml1 || x==0) {
            ml1 = currml1;
            th1 = (double)x;
        }
        if(currml2 > ml2 || x==0) {
            ml2 = currml2;
            th2 = (double)x;
        }
        if(currml3 > ml3 || x==0) {
            ml3 = currml3;
            th3 = (double)x;
        }
        if(currml4 > ml4 || x==0) {
            ml4 = currml4;
            th4 = (double)x;
        }
    }
    if((ml1 + ml2 + ml3 + ml4 > mlM4[0] + mlM4[1] + mlM4[2] + mlM4[3]) || nloci4ml1[0] == -1){
        mlM4[0] = ml1;
        mlM4[1] = ml2;
        mlM4[2] = ml3;
        mlM4[3] = ml4;
        thetaM4[0] = th1;
        thetaM4[1] = th2;
        thetaM4[2] = th3;
        thetaM4[3] = th4;
        *nloci4ml1 = numloci4m1;
        *nloci4ml2 = numloci4m2;
        *nloci4ml3 = numloci4m3;
    }
    return 0;
}
int findml4coalsorted(int numloci4m1,int numloci4m2,int numloci4m3,int *sortedlociml,int nloci,int steps,double *mlM4,double **ml)
{
    int x,y;
    double currml1=0;
	double currml2=0;
	double currml3=0;
	double currml4=0;
	double ml1=0;
	double ml2=0;
	double ml3=0;
	double ml4=0;
    
    for(x=0;x<steps;x++) {
        currml1 = currml2 = currml3 = currml4 = (double)0.;
        for(y=0;y<numloci4m1;y++) currml1 += ml[sortedlociml[y]][x];
        for(y=numloci4m1;y<numloci4m1+numloci4m2;y++) currml2 += ml[sortedlociml[y]][x];
        for(y=numloci4m1+numloci4m2;y<numloci4m1+numloci4m2+numloci4m3;y++) currml3 += ml[sortedlociml[y]][x];
        for(y=numloci4m1+numloci4m2+numloci4m3;y<nloci;y++) currml4 += ml[sortedlociml[y]][x];
        if(currml1 > ml1 || x==0) {
            ml1 = currml1;
        }
        if(currml2 > ml2 || x==0) {
            ml2 = currml2;
        }
        if(currml3 > ml3 || x==0) {
            ml3 = currml3;
        }
        if(currml4 > ml4 || x==0) {
            ml4 = currml4;
        }
    }
    if((ml1 + ml2 + ml3 + ml4 > mlM4[0] + mlM4[1] + mlM4[2] + mlM4[3])){
        mlM4[0] = ml1;
        mlM4[1] = ml2;
        mlM4[2] = ml3;
        mlM4[3] = ml4;
    }
    return 0;
}
int findml4sortedmh(int numloci4m1,int numloci4m2,int numloci4m3,int *sortedlociml,int *relnloci4ml1,int *relnloci4ml2,int *relnloci4ml3,int *nloci4ml1,int *nloci4ml2,int *nloci4ml3,int nloci,int steps,double *relthetaM4,double *relmlM4,double *thetaM4,double *mlM4,double **ml,double logranloci)
{
    int x,y;
    double th1=0;
	double th2=0;
	double th3=0;
	double th4=0;
    double currml1=0;
	double currml2=0;
	double currml3=0;
	double currml4=0;
	double ml1=0;
	double ml2=0;
	double ml3=0;
	double ml4=0;
    
    for(x=0;x<steps;x++) {
        currml1 = currml2 = currml3 = currml4 = (double)0.;
        for(y=0;y<numloci4m1;y++) currml1 += ml[sortedlociml[y]][x];
        for(y=numloci4m1;y<numloci4m1+numloci4m2;y++) currml2 += ml[sortedlociml[y]][x];
        for(y=numloci4m1+numloci4m2;y<numloci4m1+numloci4m2+numloci4m3;y++) currml3 += ml[sortedlociml[y]][x];
        for(y=numloci4m1+numloci4m2+numloci4m3;y<nloci;y++) currml4 += ml[sortedlociml[y]][x];
        if(currml1 > ml1 || x==0) {
            ml1 = currml1;
            th1 = (double)x;
        }
        if(currml2 > ml2 || x==0) {
            ml2 = currml2;
            th2 = (double)x;
        }
        if(currml3 > ml3 || x==0) {
            ml3 = currml3;
            th3 = (double)x;
        }
        if(currml4 > ml4 || x==0) {
            ml4 = currml4;
            th4 = (double)x;
        }
    }
    /*acceptation of ml using MH*/
    if((((ml1 + ml2 + ml3 + ml4) - (relmlM4[0] + relmlM4[1] + relmlM4[2] + relmlM4[3])) >= logranloci)
        || nloci4ml1[0] == -1){
        relmlM4[0] = ml1;
        relmlM4[1] = ml2;
        relmlM4[2] = ml3;
        relmlM4[3] = ml4;
        relthetaM4[0] = th1;
        relthetaM4[1] = th2;
        relthetaM4[2] = th3;
        relthetaM4[3] = th4;
        *relnloci4ml1 = numloci4m1;
        *relnloci4ml2 = numloci4m2;
        *relnloci4ml3 = numloci4m3;
    }
    /*find best ml*/
    if((ml1 + ml2 + ml3 + ml4 > mlM4[0] + mlM4[1] + mlM4[2] + mlM4[3]) || nloci4ml1[0] == -1){
        mlM4[0] = ml1;
        mlM4[1] = ml2;
        mlM4[2] = ml3;
        mlM4[3] = ml4;
        thetaM4[0] = th1;
        thetaM4[1] = th2;
        thetaM4[2] = th3;
        thetaM4[3] = th4;
        *nloci4ml1 = numloci4m1;
        *nloci4ml2 = numloci4m2;
        *nloci4ml3 = numloci4m3;
    }
    return 0;
}
int findml4coalsortedmh(int numloci4m1,int numloci4m2,int numloci4m3,int *sortedlociml,int nloci,int steps, double *relmlM4,double *mlM4,double **ml,double logranloci)
{
    int x,y;
    double currml1=0;
	double currml2=0;
	double currml3=0;
	double currml4=0;
	double ml1=0;
	double ml2=0;
	double ml3=0;
	double ml4=0;
    
    for(x=0;x<steps;x++) {
        currml1 = currml2 = currml3 = currml4 = (double)0.;
        for(y=0;y<numloci4m1;y++) currml1 += ml[sortedlociml[y]][x];
        for(y=numloci4m1;y<numloci4m1+numloci4m2;y++) currml2 += ml[sortedlociml[y]][x];
        for(y=numloci4m1+numloci4m2;y<numloci4m1+numloci4m2+numloci4m3;y++) currml3 += ml[sortedlociml[y]][x];
        for(y=numloci4m1+numloci4m2+numloci4m3;y<nloci;y++) currml4 += ml[sortedlociml[y]][x];
        if(currml1 > ml1 || x==0) {
            ml1 = currml1;
        }
        if(currml2 > ml2 || x==0) {
            ml2 = currml2;
        }
        if(currml3 > ml3 || x==0) {
            ml3 = currml3;
        }
        if(currml4 > ml4 || x==0) {
            ml4 = currml4;
        }
    }
    /*acceptation of ml using MH*/
    if((((ml1 + ml2 + ml3 + ml4) - (relmlM4[0] + relmlM4[1] + relmlM4[2] + relmlM4[3])) >= logranloci)){
        relmlM4[0] = ml1;
        relmlM4[1] = ml2;
        relmlM4[2] = ml3;
        relmlM4[3] = ml4;
    }
    /*find best ml*/
    if((ml1 + ml2 + ml3 + ml4 > mlM4[0] + mlM4[1] + mlM4[2] + mlM4[3])){
        mlM4[0] = ml1;
        mlM4[1] = ml2;
        mlM4[2] = ml3;
        mlM4[3] = ml4;
    }
    return 0;
}
