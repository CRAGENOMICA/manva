/*Hudson streec file modified*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SEGINC 1000
#define SELECTION_RECW_ALLOWED	1
#define SELECTION_CLASSIC 0
#define FSTFROMNODIVISION 0

long int nchrom;
long int begs;
unsigned long nsegs;
unsigned long nlinks;
static long int *nnodes = NULL;
double t, cleft,pc,lnpc;

long int total_nts,sel_nts_glob,new_chrom;	/*in case selection*/
int ifsel_glob;					/*in case selection*/

static unsigned long seglimit = SEGINC;
static unsigned long maxchr;

/*Cada chrom té un nombre de segments, els quals es troben a una població. La direccio dels segments es indicada*/
/*Cada chrom té els seus propis segments. S'indica l'inici i final de cada segment, així com el nombre del segment amb el que continua aquest chrom*/
struct seg {
    unsigned long beg;
    unsigned long end;
    int desc;
};
struct chromo {
    unsigned long nseg;
    long int pop;
    struct seg *pseg;
};
static struct chromo *chrom = NULL;
/*Els arbres estan fets en nodes (contenen el temps de coalescencia i el nombre del node amb el que conecta per dalt)*/
/*tambe tenim els segments-length que contenen l'inici nt(eg 10-500, es 10) del segment, la direccio al nodes (arbre) al qual esta conectat, i el nombre del següent segment (eg seria el segment que conte des del 501). Hi ha un segl per cada fragment produit per recombinacio que te un arbre associat.*/
struct node {
    int abv;
    int ndes;
    double time;
} *ptree1,*ptree2;
struct segl {
    long int beg;
    struct node *ptree;
    long int next;
};
static struct segl *seglst = NULL;

/* El resultat és la matriu de segments-length, que té associat tots els arbres per a cada segment. */
/* Les matrius chrom i pseg només són neccesàries per fer la coalescència, però després ja no són necessàries. */
struct segl *segtre_mig(long int npop,int nsam,int *inconfig,unsigned long nsites,double r,double f,double track_len,
    double mig_rate,unsigned long *pnsegs/*,unsigned long iteration*/,
    double *factor,int ifselection, 
    double pop_sel, double sinit,double pop_size,long int sel_nt,int *all_sel,int *selnsam,
    int nintn,double *nrec,double *npast, double *tpast,
    int split_pop, double time_split, double time_scoal, double factor_anc, double *freq,
    double tlimit,int iflogistic,double Tts,double factor_chrn)
{
    int j,dec,c1,c2,ind,rchrom,intn;
    int migrant,*config;
    long int pop,source_pop;/*modificat a long int*/
    
    double trm,tcoal,ttemp,rft,clefta,alphag;
    double prec,cin,prect,mig,ran,coal_prob,prob,rdum; 
    int re(int),ca(int,unsigned long,int,int);
    void pick2_chrom(long int,int *,int *,int *);/*modificat a long int*/
    int cinr(int,unsigned long),cleftr(int); 

    double ran1(void);
    double bnldev(double,int);
    
    double xdt,eps,eps2n,tf,ts,tot_rec,coal_prob0,coal_prob1,tcoal0,tcoal1,freqend;	/*for selection*/
    double poisso(double xm); /*finding the time xdt reach to eps (before the deterministic selective event starts)*/
	/*/double freqpa;*/ /*for selective module*/
	/*long int timeg;*/ /*for selective module*/
	double nref,Ts;
	double xacc,x1,x2;
	int zbracn(double(*)(double,double,double,double,double,double,double,double,double,double),double *,double *,double,double,double,double,double,double,double,double,double);
	double zriddrn(double(*)(double,double,double,double,double,double,double,double,double,double),double,double,double,double,double,double,double,double,double,double,double,double);
	double functiont_freqp_sel(double,double,double,double,double,double,double,double,double,double);
	double functiont_freqq_nsel(double,double,double,double,double,double,double,double,double,double);
	double functiont_logistic(double,double,double,double,double,double,double,double,double,double);
	
	double tpop;
	int coalpop;
	int a,b;
	double migran,migsum,migppop;

	#if SELECTION_CLASSIC
	double dt,rec_prob;
	#endif
	
	double sumfreq,pip,spip;			/*refugia*/
    
    #if FSTFROMNODIVISION == 1
		int cf; /*per Fst*/
	#endif

    sel_nts_glob = sel_nt;			/*for selection*/
    ifsel_glob = ifselection;			/*for selection*/

    
    if(chrom == NULL) {
        maxchr = nsam + 50;	/* + 50 ... */
        if((chrom = (struct chromo *)malloc((unsigned)(maxchr*sizeof(struct chromo)))) == NULL)
            perror("malloc error. segtre_mig.1");
    }
    if(nnodes == NULL) {
        if((nnodes = (long int *)malloc((unsigned)(seglimit*sizeof(long int)))) == NULL)
            perror("malloc error. segtre_mig.2");
    }
    if(seglst == NULL) {
        if((seglst = (struct segl *)malloc((unsigned)(seglimit*sizeof(struct segl)))) == NULL)
            perror("malloc error. segtre_mig.3");
    }
    /*nou el 9.5.03*/
    if((int)maxchr < nsam) {
        maxchr = (unsigned long)nsam + 50;
        if((chrom = (struct chromo *)realloc(chrom,(maxchr*sizeof(struct chromo)))) == NULL)
            perror("malloc error. segtre_mig.1");
    }
    
    config = (int *)malloc((unsigned)((npop+1)*sizeof(int)));	/* nou vector de mida de poblacions */
    if(config==NULL) perror("malloc error segtre_mig.4");
    for(pop=0;pop<npop;pop++) config[pop] = inconfig[pop];	/* inicialitzar amb valors del input */
    for(pop=ind=0;pop<npop;pop++)	/* matriu chrom, indica els individus de cada població i els segments/ind */
        for(j=0;j<inconfig[pop];j++,ind++) {
            chrom[ind].nseg = 1;
            if(!(chrom[ind].pseg = (struct seg *)malloc((unsigned)sizeof(struct seg))))
                perror("malloc error segtre_mig.5");
            (chrom[ind].pseg)->beg = 0;
            (chrom[ind].pseg)->end = nsites - 1;
            (chrom[ind].pseg)->desc = ind;
            chrom[ind].pop = pop;
        }
    seglst[0].beg = 0;
    if(!(seglst[0].ptree = (struct node *)calloc((unsigned)(2*nsam),sizeof(struct node))))
        perror("calloc error segtre_mig.6");
    nnodes[0] = nsam -1;	/* nombre de fragments per a cada segment .. */
    nchrom = nsam; 		/* nombre de fragments que hi han a l'arbre .. */
    nlinks = ((unsigned long)(nsam))*(nsites-1);
    
    total_nts = nlinks; 		/*included for selection*/
    if(ifselection == 0) sel_nt = 0;	/*included for selection*/
   
    nsegs = 1;
    t = (double)0;
	r /= (double)(nsites-1);
	trm = 0.;

    if(ifselection == 0 || (ifselection==1 && pop_sel < (double)10)) {
        if(split_pop == 0) {
        /********************* ROUTINE FROM HUDSON modified *************************/
            if(f > 0.0) pc = (track_len - (double)1)/track_len;	/* tot això conversió gènica */
            else pc = (double)1;
            lnpc = (double)log(pc);
            cleft = nsam * ((double)1 - pow(pc,(double)(nsites-1)));
            rft = r*f*track_len;
            intn = 1;
            alphag = (double)2*(double)10/(tpast[intn]-tpast[intn-1]);/*logistic*/
			if(iflogistic) {
				Ts = (double)Tts;
				nref = ((double)nrec[1] + ((double)npast[1]-(double)nrec[1])/((double)1 + (double)exp(-(double)alphag*(Ts-((double)tpast[0]+(double)tpast[1])/(double)2))));
			}
            /******************************** Main loop *********************************/
            while(nchrom > 1) {                 
                /*modification to increase the spped with high R, if t > tlimit r = 0.;*/
                if(t > tlimit)
                    r = (double)0; 
                                
                prec = (double)nlinks * r;	/* prob. recombinació. Cada posicio entre 2nt te una prob. de recombinar */
                cin  = (double)nlinks * r * f;/* prob. cin event: INACTIVAT a les funcions ca,rec,xover ...*/
                clefta = cleft * rft;/* prob. cleft event: INACTIVAT a les funcions ca,rec,xover ...*/
                prect = prec + cin + clefta;/* prob. recombinació + cin + cleft */
                 
                coal_prob = (double)0;
                for(pop=0;pop<npop;pop++)
                    coal_prob += ((double)config[pop])*(config[pop]-(double)1)*factor_chrn/factor[pop]; /* arbres en funció de 4No */
                    
                mig = nchrom * mig_rate;	/* prob. migració. Cada ind. té una probabilitat de migrar */
                
                if(prect+mig > (double)0) {
                    while((rdum = ran1()) == 1.0);/* de fet, mai pot ser 1, es [0.1) en ran1 ... */
                    trm = -log((double)1-rdum)/(prect+mig);/* temps a la recombinació, (conversió), migració o extincio*/
                }/* no es troben relacionats amb No, per tant no afecta canvi demogràfic */
                
				coalpop = 0;
				tcoal = tpop = trm + 999999.;
				for(pop=0;pop<npop;pop++) {
					coal_prob = ((double)config[pop])*(config[pop]-(double)1)*factor_chrn/factor[pop]; /* arbres en funció de 4No in a diploid sp*/
					if(coal_prob > (double)0) {
						while((rdum = (double)ran1()) == 1.0);
						/* ponderar canvi demogràfic: ACTIVAT. totes les poblacions canvien a la vegada*//*5.5.03*/
						if(iflogistic==0) {
							/*instantaneous*/
							tpop = -(double)log((double)1.0 - (double)rdum)/coal_prob;
							tpop *= npast[intn];
						}
						else {
							/*logistic*/
							tpop = -(double)log((double)1.0-(double)rdum)/coal_prob;
							/*calculate the range*/
							if(intn > 1) Ts = (double)0;
							x1 = (double)0.;
							x2 = (double)0.5;
							if(zbracn(functiont_logistic,&x1,&x2,(double)tpop,(double)t,(double)tpast[intn-1],(double)tpast[intn],(double)alphag,(double)nrec[intn],(double)npast[intn],nref,Ts) == 1) {
								xacc = (double)1e-6*x2; 
								tpop = (double)zriddrn(functiont_logistic,x1,x2,xacc,(double)tpop,(double)t,(double)tpast[intn-1],(double)tpast[intn],(double)alphag,(double)nrec[intn],(double)npast[intn],nref,Ts);
							}
							/*printf("%f\n",tcoal);*/
						}
					}
					if(tpop < tcoal) {
						tcoal = tpop;
						coalpop = pop;
					}
				}
/*
				if(coal_prob > (double)0) {
                    while((rdum = ran1()) == 1.0);
					if(iflogistic==0) {
						tcoal = -(double)log(1.0 - rdum)/coal_prob;
						tcoal *= npast[intn];
					}
					else {
						tcoal = -log(1.0-rdum)/coal_prob;
						if(intn > 1) Ts = (double)0;
						x1 = (double)0.;
						x2 = (double)0.5;
						if(zbracn(functiont_logistic,&x1,&x2,(double)tcoal,(double)t,(double)tpast[intn-1],(double)tpast[intn],(double)alphag,(double)nrec[intn],(double)npast[intn],nref,Ts) == 1) {
							xacc = (double)1e-6*x2; 
							tcoal = (double)zriddrn(functiont_logistic,x1,x2,xacc,(double)tcoal,(double)t,(double)tpast[intn-1],(double)tpast[intn],(double)alphag,(double)nrec[intn],(double)npast[intn],nref,Ts);
						}
					}
                }
                else tcoal = trm + 999999.;
*/                
                if((prect + mig) > (double)0) ttemp = ((trm < tcoal) ? trm : tcoal);
                else ttemp = tcoal;
                /*en cas canvi No. Si ttemp > tpast aleshores tornem a calcular probs. i t=ttemp : Ara ACTIVAT*//*5.5.03*/
                if((intn < nintn) && (t + ttemp > tpast[intn])) {
                    t = tpast[intn++];
                    alphag = (double)2*(double)10/(tpast[intn] - tpast[intn-1]);/*logistic*/
                }
                else {/*no changes in size of pop*/
                    t += ttemp;
                    if(((prect+mig ) > (double)0) && (trm < tcoal)) {
                        if((ran = ran1()) < (prec/(prect + mig))) {
                            /* recombination */
                            rchrom = re(nsam);	/* rchrom és l'individu a on té lloc la recombinació */		
                            config[chrom[rchrom].pop] += 1;  /* aumenta 1 la població on pertany l'individu */
                        }
                        else {
                            /**/
                            if(ran  < ((prec + clefta)/(prect + mig))) {
                                /* cleft event*/
                                rchrom = cleftr(nsam);
                                config[chrom[rchrom].pop] += 1;
                            }
                            else {
                                if(ran < (prect/(prect + mig))) {
                                    /* cin event*/
                                    rchrom = cinr(nsam,nsites);
                                    if(rchrom >= 0) config[chrom[rchrom].pop] += 1;
                                }
                                else /*if (ran < ((prect + mig)/(prect + mig)))*/
                                    {
                                    /* migration event */
									migran = mig*ran1();
									migsum = 0.;
									a = b = 0;
									migppop = config[a] * mig;
									while(migsum + migppop < migran) {
										migsum += migppop;
										b++;
										if(b==npop) {
											b=0;
											a++;
											/*if(a==npop) break;*//*error*/
										}
										migppop = config[a] * mig;
									}
									migppop = mig;
									for(j=0;j<nchrom;j++) {
										if(chrom[j].pop == a) {
											if(migsum + migppop < migran) migsum += migppop;
											else break;
										}
									}
									migrant = j;
									source_pop = b;
									/*end new code*/
									config[chrom[migrant].pop] -= 1;
									config[source_pop] += 1;
									chrom[migrant].pop = source_pop;
								}
                            }
                        }
                    }
                    else /*if(((prect+mig) > 0.0 && tcoal <= trm) || (prect+mig) == 0.0)*/ {
                        /* coalescent event */
                        /* pick the two, c1, c2 */
                        /*
						ran = ran1();
                        prob = 0.0;	
                        for(pop=0;pop<npop;pop++) {
                            prob += (((double)config[pop])*(config[pop]-(double)1)/factor[pop])*factor_chrn/coal_prob;
                            if(ran < prob) break;
                        }
                        pick2_chrom(pop,config,&c1,&c2);	
						*/
                        pick2_chrom(coalpop,config,&c1,&c2);	
                        dec = ca(nsam,nsites,c1,c2);	/* dec és el nombre de fragments a restar */
                        config[coalpop] -= dec;			/* si hi ha MRCA aleshores és més d'un */
                    }/*coal event*/
                }/*event*/
            }/*en cas chrom > 1*/
        }
        else {
            /*********************** REFUGIA **************************/
            sumfreq = 0.;
            for(pop=0;pop<npop;pop++) sumfreq += freq[pop];
            for(pop=0;pop<npop;pop++) freq[pop] /= sumfreq;
            /* Main loop */                    
            while(nchrom > 1) {

                /*modification to increase the spped with high R, if t > tlimit r = 0.;*/
                if(t > tlimit) 
                    r = 0.; 

                prec = prect = (double)nlinks * r;
                
                if((t >= time_split) && (t < time_scoal)) {
                    coal_prob = (double)0;
                    for(pop=0;pop<npop;pop++)
                        coal_prob += ((double)config[pop])*(config[pop]-1.0)*factor_chrn/factor[pop];
                }
                else {
                    if(t < time_split)
                        coal_prob = ((double)config[0])*((double)config[0]-1.0)*factor_chrn;
                    else /*t >= time_scoal*/
                        coal_prob = ((double)config[0])*((double)config[0]-1.0)*factor_chrn/factor_anc;
                }
                if(coal_prob > 0.0) {
                    while((rdum = ran1()) == 1.0);
                    tcoal = -log(1.0 - rdum)/coal_prob;
                }
                else tcoal = trm + 999999.;
               
			    if(t<time_split || t>=time_scoal) mig = 0.;
				else mig = nchrom * mig_rate;	/* prob. migració. Cada ind. té una probabilitat de migrar */
				
				if(prect+mig > (double)0) {
                    while((rdum = ran1()) == 1.0);	
                    trm = -log(1.0-rdum)/(prect+mig);	
                }
                
                if((prect) > (double)0) ttemp = ((trm < tcoal) ? trm : tcoal);
                else ttemp = tcoal;
                if((t < time_split) && (t+ttemp >= time_split)) {
                    t = time_split;
                    /* From one to npop pops */
                    for(j=0;j<nchrom;j++) {
                        if(ran1() <= (1.0/*sumfreq*/ - freq[0])) {
                            pip = ran1()*(1.0/*sumfreq*/ - freq[0]);
                            spip = (double)0;
                            for(pop=1;pop<npop;pop++)
                                if((spip += freq[pop]) > pip) break;
                            chrom[j].pop = pop;
                            config[pop] += 1;
                            config[0] -= 1;
                        }
                    }
                }
                else {
                    if((t >= time_split) && (t < time_scoal) && (t+ttemp >= time_scoal)) {
                        t = time_scoal;
                        /* From npop to one pop */
                        for(j=0;j<nchrom;j++) chrom[j].pop = 0;
                        config[0] = nchrom;
                        for(pop=1;pop<npop;pop++) config[pop] = 0;
                        /*printf("config[0]=%d\tconfig[1]=%d\n",config[0],config[1]);*/
                    }
                    else {
                        t += ttemp;
                                                    
                        if(((prect+mig) > (double)0) && (trm < tcoal)) {
							if((ran = ran1()) < (prect/(prect + mig))) {
								/* recombination */
								rchrom = re(nsam);		
								config[chrom[rchrom].pop] += 1;
							}
							else {
								/* migration event */
								migrant = (int)((double)nchrom*ran1());
								while((source_pop = (int)(npop*ran1())) == chrom[migrant].pop);
								config[chrom[migrant].pop] -= 1;
								config[source_pop] += 1;
								chrom[migrant].pop = source_pop;
							}
                        }
                        else {
                            /* coalescent event */
                            /* pick the two, c1, c2 */
                            while((ran = ran1()) == 1.0);
                            prob = (double)0;
                            if((t >= time_split/* + ttemp*/) && (t/*  + 1e-10 *//*error*/< time_scoal/* + ttemp */)) {
                                for(pop=0;pop<npop;pop++) {
                                    prob += (((double)config[pop])*(config[pop]-1.)*factor_chrn/factor[pop])/coal_prob;
                                    if(ran < prob) break;
                                }
                            }
                            else pop = 0;
                            pick2_chrom(pop,config,&c1,&c2);
                            dec = ca(nsam,nsites,c1,c2);
                            config[pop] -= dec;
                        }
                    }
                }
            }
        }
    }
    else {/******************************SELECTION MODULE*********************************************/
        /*EQUATION 3A and 3B IN STEPHAN 1992.*/
        /*selected pop is 0. Non-selected pop is 1*/
		eps = /**/1./(pop_sel);/**//*100./(2.*(double)pop_size);*//*/1./(2.*(double)pop_size);*//*5./(pop_sel/2.);*//*frequency of selective allele when selective phase finish*/
		/*look for the generation when the determinsitc event starts (ts). Binomial, but using a poisson given p << 0.1 and n >> 30*/
		/**/
		eps2n = eps * (double)2* (double)pop_size;
		/*
		do{
			timeg = 0;
			freqpa = (double)1;
			while(freqpa > (double)0 && freqpa < eps2n) {
				freqpa = poisso((double)freqpa);
				timeg++;
			}
		}while(freqpa < eps2n);
		ts = (double)timeg/((double)4*(double)pop_size);
		*/
		/**//*we consider no non-selected allele is recombined to our sample from 1 to 1-eps generations.*/
		/**/ts = (double)0;/**/
		freqend = eps;/*in case selective mutation is a new mutation, if not, the user would define freqend: not done yet*/
		tf = (pop_sel*(ts)+(double)log((1.-freqend-eps+eps*freqend)/(eps*freqend)))/(pop_sel);/*time starting/ending selective phase*//*the strength of selection is 2Nes because affects the individual*/
		*all_sel = config[0];
        if(sinit < (double)0) {/*selection has not fixed yet the selective allele */
			xdt = (double)1 - eps/(eps+((double)1-eps)*(double)exp((double)(-pop_sel * ((double)0 -sinit-ts))));/*freq sel. allele in N individuals*/
			if(-sinit < ts) {
				*all_sel = (int)config[0];/*all alleles have the favorable mutation*/
				for(ind=0;ind<*all_sel;ind++) selnsam[ind] = ind;/*fav lines in a vector (will be useful in mig+sel)*/
				for(ind=*all_sel;ind<config[0];ind++) chrom[ind].pop = 1; /*chrom with unfav allele*/
				config[0] = nchrom;
				config[1] = 0;
			} 
			else {
				if(-sinit >= tf) {
					*all_sel = (int)0;/*any allele has the favorable mutation*/
					for(ind=0;ind<*all_sel;ind++) selnsam[ind] = ind;/*fav lines in a vector (will be useful in mig+sel)*/
					for(ind=*all_sel;ind<config[0];ind++) chrom[ind].pop = 1; /*chrom with unfav allele*/
					config[0] = 0;
					config[1] = nchrom;
				}
				else {
					*all_sel = (int)bnldev((double)xdt,config[0]);/*how many alleles have the favorable mutation*/
					for(ind=0;ind<*all_sel;ind++) selnsam[ind] = ind;/*fav lines in a vector (will be useful in mig+sel)*/
					for(ind=*all_sel;ind<config[0];ind++) chrom[ind].pop = 1; /*chrom with unfav allele*/
					config[1] = config[0] - *all_sel;
					config[0] = *all_sel;
				}
			}
        }
        else xdt = 1.;
        /* selective position is at 'sel_nt' starting from the left of the studied sequence (it can be neg)*/
        /* recombination includes all the region until the selected position */
        /* recombination is only performed when nt from the studied region are involved. */
        /* For the rest we only assign the group given xdt and r=r/(nsites-1). sites are the sites in the region.*/
		total_nts = ((((long int)nsites-1 > (long int)sel_nt-1) ? (long int)nsites-1 : (long int)sel_nt-1) - (((long int)sel_nt < 0) ? (long int)sel_nt : 0));
        /*it is necessary to know the total positions that are able to recombine (until the selected pos.)*/
        r *= (double)(nsites-1)/(double)total_nts; 
        /*in selection module, R=4Nr for the fragment studied + all until the selective position*/
        total_nts *= nchrom; /*like nlinks but for all the region*/
        /*************************************** Main loop for selection *****************************************/
		#if SELECTION_CLASSIC
        while(nchrom > 1) {
            /*modification to increase the spped with high R, if t > tlimit r = 0.;*/
            if(t > tlimit) 
                r = (double)0; 

            prect = (double)nlinks * r;        /* 'tot_rec' is r in all region.*/
            tot_rec = (double)total_nts * r;   /*like prect but for all the region*/      

            /************ non-selective phase ****************************/
            if((t < sinit + ts) || (t >= sinit + tf)) { 
                coal_prob = ((double)nchrom)*(nchrom-(double)1)*factor_chrn;
                if(prect > (double)0) {
                    while((rdum = ran1()) == (double)1);	
                    trm = -log(1.0-rdum)/(prect);	
                }
                while((rdum = ran1()) == 1.0);
                tcoal = -log(1.0 - rdum)/coal_prob;
                if((prect) > (double)0) ttemp = ((trm < tcoal) ? trm : tcoal);
                else ttemp = tcoal;
                t += ttemp;
				if((t < sinit + ts && t - ttemp < sinit + ts) || (t >= sinit + tf && t - ttemp >= sinit + tf)) {                                            
					if(((prect) > (double)0) && (trm < tcoal)) {
						#if SELECTION_RECW_ALLOWED
						/* recombination */
						rchrom = re(nsam);		
						config[chrom[rchrom].pop] += 1;
						#endif
					}
					else {
						/* coalescent event */
						if(config[0] == 0) pop = 1;
						else pop = 0;
						/* pick the two, c1, c2 */
						pick2_chrom(pop,config,&c1,&c2);	
						dec = ca(nsam,nsites,c1,c2);		
						config[pop] -= dec;
					}
				}
				else {
					if((t >= sinit + ts) && (t - ttemp < sinit + ts)) {
						/*starting selective process. We assume no recombination between selected-deselected occured before*/
						t = sinit + ts;
						xdt = 1. - eps;
					}
				}
            }
            /************************* selective phase ***************************/
            else {
				/**/dt = (tf-ts)/((tf-ts)*(4.*(double)pop_size));/**//*do each generation individually in discrete steps*/
				/*dt = 1./(100.*pop_sel/2.);*//*braverman*/
				xdt = (double)1 - eps/(eps+((double)1-eps)*(double)exp(-pop_sel * (t+dt-sinit-ts)));
				coal_prob = dt * factor_chrn * (((double)config[0])*(config[0]-(double)1)/(           xdt) + 
							                  ((double)config[1])*(config[1]-(double)1)/((double)1- xdt));
				rec_prob  = dt * tot_rec;
				while(0.1 < (coal_prob+rec_prob)) {/*in case two or more events can happen at the same time...*/
					dt =  dt * 0.1;
					xdt = (double)1 - eps/(eps+((double)1-eps)*(double)exp(-pop_sel * (t+dt-sinit-ts)));
					coal_prob = dt*factor_chrn * (((double)config[0])*(config[0]-(double)1)/((double)xdt) + 
									              ((double)config[1])*(config[1]-(double)1)/((double)1- xdt));
					rec_prob  = dt*tot_rec;
				}
				ttemp = dt;
				if(t+ttemp < sinit+tf) {/*in case xdt <= eps and config[0] > 0, we join ind pop 0 to 1*/
                    t += ttemp;
                    /*the program would be faster if we keep the product of coal_prob+rec_prob until having an event*/
                    if((rdum = ran1()) < (coal_prob + rec_prob)) { /*then we have one event, otherwise nothing*/
                        if((rdum = ran1()) >= coal_prob/(coal_prob + rec_prob)) {/*rec or coal?*/
                            /* recombination */
                            while((rdum = ran1()) == 1.0);
                            if(rdum < (double)nlinks/(double)total_nts) { /*recombination is only effective within nlinks*/
                                #if SELECTION_RECW_ALLOWED
                                rchrom = re(nsam);
                                while((rdum = ran1()) == 1.0);
								rdum = freqend+rdum*((double)1-(double)eps-freqend);/*rdum goes from freqend to 1-eps*/
                                if(rdum < xdt) chrom[new_chrom].pop = 0;
                                else chrom[new_chrom].pop = 1;
                                config[chrom[new_chrom].pop] += 1;
                                #endif
                            }
                            else {/*recombination in one chrom outside studied region.We do changes only some times.Slow*/
                                while((rdum = ran1()) == 1.0);
                                new_chrom = rdum * nchrom;	/*find the chrom*/
                                config[chrom[new_chrom].pop] -= 1;/*we check the pop. We erase that from the former pop*/
                                while((rdum = ran1()) == 1.0);
								rdum = freqend+rdum*((double)1-(double)eps-freqend);/*rdum goes from freqend to 1-eps*/
                                if(rdum < xdt) chrom[new_chrom].pop = 0;
                                else chrom[new_chrom].pop = 1;
                                config[chrom[new_chrom].pop] += 1;
                                /*we assign the new pop. many times is the same than before*/
                            }
                        }
                        else { /*coalescent*/
                            while((ran = ran1()) == 1.0);
                            prob = (ttemp*((double)config[0])*(config[0]-(double)1)*factor_chrn/xdt)/coal_prob;
                            if(ran < prob) pop = 0;
                            else pop = 1; 
                            /* pick the two, c1, c2 */
                            pick2_chrom(pop,config,&c1,&c2);	
                            dec = ca(nsam,nsites,c1,c2);	
                            config[pop] -= dec;
                        }
                    }
                }
                else {/*End selective phase. we cut at 1-eps.*/
					t = sinit + tf;
					while(config[0] > 1) {
						/*coalescent*/
						pick2_chrom(pop=0,config,&c1,&c2);	
						dec = ca(nsam,nsites,c1,c2);	
						config[pop] -= dec;
						t += 1E-17;
					}
					for(new_chrom=0;new_chrom<nchrom;new_chrom++) {
						if(chrom[new_chrom].pop == 0) {
							config[0] -= 1;
							chrom[new_chrom].pop = 1;
							config[1] += 1;
							break;
						}
					}
                }
            }
        }
		/************ IN CASE USING SELECTIVE APPROACH BY INTEGRATING THE DETERMINISTIC EQUATION ***********/
		#else
		while(nchrom > 1) {
            /*modification to increase the speed with high R, if t > tlimit, r = 0.;*/
            if(t > tlimit) 
                r = (double)0; 

            prect = (double)nlinks * (double)r;        /* 'tot_rec' is r in all studied region.*/
            tot_rec = (double)total_nts * (double)r;   /*like prect but for all the region (in case selection is out from studied region is bigger)*/      

			/*recombination*/
			if(t < sinit + ts || t >= sinit + tf) 
				tot_rec = prect;/*in case no selection, we do not check recombinants in the non-studied region*/
			if(tot_rec > (double)0) {
				while((rdum = ran1()) == 1.0);	
				trm = -log((double)1-rdum)/tot_rec;
			}
			/*selected alleles*/
			if(t < sinit + tf) 
				coal_prob0 = ((double)config[0])*(config[0]-(double)1)*factor_chrn;
			else coal_prob0 = (double)0;
			if(coal_prob0 > (double)0) {
				while((rdum = ran1()) == 1.0);				
				tcoal0 = -log((double)1 - rdum)/coal_prob0;
				/*integrating xdt in function of t, we calculate the time of coalescent directly*/
				if((t >= sinit + ts) && (t < sinit + tf)) {
					/*calculate the range*/
					x1 = 0.;
					x2 = (double)sinit + (double)tf;
					if(zbracn(functiont_freqp_sel,&x1,&x2,(double)tcoal0,(double)t,(double)ts,(double)sinit,(double)eps,(double)pop_sel,(double)0,(double)0,(double)0) == 1) {
						/*estimate the value of T*/
						xacc = (double)1e-6*x2; /* accuracy*/
						tcoal0 = (double)zriddrn(functiont_freqp_sel,x1,x2,xacc,(double)tcoal0,(double)t,(double)ts,(double)sinit,(double)eps,(double)pop_sel,(double)0,(double)0,(double)0);
					}
				}
			}
			else tcoal0 = trm + 99999.;
			/*non-selected alleles*/
			if(t >= sinit + ts)
				coal_prob1 = ((double)config[1])*(config[1]-(double)1)*factor_chrn;
			else coal_prob1 = (double)0;
			if(coal_prob1 > (double)0) {
				while((rdum = ran1()) == 1.0);				
				tcoal1 = -log((double)1 - rdum)/coal_prob1;
				/*integrating 1-xdt in function of t, we calculate the time of coalescent directly*/
				if((t >= sinit + ts) && (t < sinit + tf)) {
					/*calculate the range*/
					x1 = 0.;
					x2 = (double)sinit + (double)tf;
					if(zbracn(functiont_freqq_nsel,&x1,&x2,(double)tcoal1,(double)t,(double)ts,(double)sinit,(double)eps,pop_sel,0,0,0) == 1) {
						/*estimate the value of T*/
						xacc = (double)1e-6*x2; /* accuracy*/
						tcoal1 = (double)zriddrn(functiont_freqq_nsel,x1,x2,xacc,(double)tcoal1,(double)t,(double)ts,(double)sinit,(double)eps,pop_sel,0,0,0);
					}
				}
			}
			else tcoal1 = trm + 99999.;
			
			/*decision*/
			tcoal = ((tcoal0 < tcoal1) ? tcoal0 : tcoal1); 
			if((tot_rec) > (double)0) ttemp = ((trm < tcoal) ? trm : tcoal);
			else ttemp = tcoal;
			
			if(t+ttemp >= sinit + ts && t < sinit + ts) {
				/*starting selective process. We assume no recombination between selected-deselected occured before*/
				ttemp = sinit + ts - t + 1E-17/*precission error*/;
			}
			if((t+ttemp >= sinit + tf && t < sinit + tf && t >= sinit + ts) || (coal_prob0 == (double)0 && coal_prob1 == (double)0 && tot_rec == (double)0)) {
				/*finishing selective process at 1-eps...*/
				ttemp = sinit + tf - t + 1E-17/*precission error*/;
				while(config[0] > 1) {
					/*coalescent*/
					pick2_chrom(pop=0,config,&c1,&c2);	
					dec = ca(nsam,nsites,c1,c2);	
					config[pop] -= dec;
					t += 1E-17;
				}
				for(new_chrom=0;new_chrom<nchrom;new_chrom++) {
					if(chrom[new_chrom].pop == 0) {
						config[0] -= 1;
						chrom[new_chrom].pop = 1;
						config[1] += 1;
						break;
					}
				}
			}
			
			/*modify the frequency of selective alleles: xdt changes from ts to tf*/
			xdt = (double)1 - eps/(eps+((double)1-eps)*(double)exp(-pop_sel * (t+ttemp-sinit-ts)));
			t += ttemp;  
			
			if(ttemp == trm){
				/* recombination */
				while((rdum = ran1()) == 1.0);
				if(rdum < (double)nlinks/(double)total_nts || (t < sinit + ts || t >= sinit + tf)) { 
					/*recombination is only effective within nlinks*/
					#if SELECTION_RECW_ALLOWED
					rchrom = re(nsam);
					while((rdum = ran1()) == 1.0);
					rdum = freqend+rdum*((double)1-(double)eps-freqend);/*rdum goes from freqend to 1-eps*/
					if(rdum < xdt) chrom[new_chrom].pop = 0;
					else chrom[new_chrom].pop = 1;
					config[chrom[new_chrom].pop] += 1;
					#endif
				}
				else {/*recombination in one chrom outside studied region.We do changes only some times.Slow*/
					while((rdum = ran1()) == 1.0);
					new_chrom = (long int)((double)rdum * (double)nchrom);	/*find the chrom*/
					config[chrom[new_chrom].pop] -= 1;/*we check the pop. We erase that from the former pop*/
					while((rdum = ran1()) == 1.0);
					rdum = freqend+rdum*((double)1-(double)eps-freqend);/*rdum goes from freqend to 1-eps*/
					if(rdum < xdt) chrom[new_chrom].pop = 0;
					else chrom[new_chrom].pop = 1;
					config[chrom[new_chrom].pop] += 1;
					/*we assign the new pop. many times is the same than before*/
				}
			}
			if(ttemp == tcoal) { /*coalescent*/
				pop = ((tcoal0 < tcoal1) ? 0 : 1);
				/* pick the two, c1, c2 */
				pick2_chrom(pop,config,&c1,&c2);	
				dec = ca(nsam,nsites,c1,c2);	
				config[pop] -= dec;
			}
        }
		#endif
    }    
	/*printf("\n%f\t%f",-log(ran1())/2.,t);*/
	/************************ End of making tree *********************************/
    *pnsegs = nsegs;
    
    free(config);
    return(seglst);
}
/*************************************************** Hudson routine *******************************************/
/* recombination */
int re(int nsam)
{
    struct seg *pseg;
    unsigned long el,spot,is;
    int lsg,lsgm1,ic;
	double ran;
    
    double ran1(void);
    void xover(int,int,unsigned long);
    
    /* First generate a random x-over spot, then locate it as to chrom and seg */
    ran = ran1();
	spot = (unsigned long)floor((nlinks * ran + (double)1)); 	/* suma 1, per evitar el zero al càlcul */
    /* get chromosome number (ic) */
    for(ic=0;ic<nchrom;ic++) { 	/* Busca els valors MÀXIM i MÍNIM de l'individu escollit. (Què nt. segment i individu) */
        lsg = chrom[ic].nseg;	/* el nombre de segments a l'individu ic */
        lsgm1 = lsg - 1;	/* de fet lsg-1 conté la informació de l'ultim segment */
        pseg = chrom[ic].pseg;	/* punter al primer segment de l'individu ic */
        el = ((pseg+lsgm1)->end) - (pseg->beg); /* mida del segments (max-min) a l'individu ic*/
        if(spot <= el) break; 	/* anem restant 'el' fins trobar l'individu */
        spot -= el;
    }
    is = pseg->beg + spot - 1; 	/* posició dins l'individu ic */
    xover(nsam,ic,is);
    return(ic);
}

void xover(int nsam,int ic,unsigned long is)
{
    struct seg *pseg,*pseg2;
    long int i,lsg,lsgm1,newsg,jseg,k,in;
    double len;
    double ran();
    
    pseg = chrom[ic].pseg;	/* punter al primer segment de l'individu ic */
    lsg  = chrom[ic].nseg;	/* nombre de segments de ic */
    len  = (pseg + lsg-1)->end - pseg->beg;	/* max-min de l'individu ic */

    cleft -= (double)1 - pow(pc,len);    /* per conversió */
    /* get segment number (jseg)*/
    for(jseg=0;is >= (pseg+jseg)->end;jseg++);	/* Busca el segment a on es troba la recombinació */
    if(is >= (pseg+jseg)->beg) in = 1;		/* mira si la recombinació es troba entre els segments o dins el segment */
    else in = 0;				/* in=0 entre segments, separació d'individus, pero no es fa un nou arbre */
    newsg = lsg - jseg; 			/* Això indica el desplaçament dels fragments de l'individu */
    
    /* copy LAST part of chrom to nchrom */
    nchrom++; 					/* fem un individu més, variable externa, val per tot el fitxer */
    if((unsigned long int)nchrom >= (unsigned long int)maxchr) {
        maxchr += 50;				/* afegeix 50 independent individus cada vegada que hem d'ampliar */
        if(!(chrom = (struct chromo *)realloc(chrom,(unsigned long)(maxchr*sizeof(struct chromo)))))
            perror("realloc error. xover.1");
    }
    if(!(pseg2 = chrom[nchrom-1].pseg = (struct seg *)calloc((unsigned)newsg,sizeof(struct seg))))/* pseg2 apunta a chrom */
        perror("calloc error. xover.2");	/* pseg2 apunta a chrom[nchrom-1].pseg. Són newsg segments nous */
    chrom[nchrom-1].nseg = newsg; 		/* el nou individu te newsg segments, els de la dreta */	
    chrom[nchrom-1].pop = chrom[ic].pop;	/* la mateixa població de la qual prové, es clar *//*pero no en cas de seleccio*/
    pseg2->end = (pseg+jseg)->end;		/* el primer segment pseg2->end apunta a final del segment de jseg, ok */
    
    if(in) {					/* només al cas que haguem de fer nous arbres */
        pseg2->beg = is + 1;			/* pseg2->beg és a on hi ha hagut la rec. + 1 */
        (pseg+jseg)->end = is;			/* aleshores (pseg+jseg) acaba a is, el punt de rec., clar */
    }
    else pseg2->beg = (pseg+jseg)->beg;		/* si no es fan nous arbres, l'inici de pseg2 és a l'inici del segment */
    
    pseg2->desc = (pseg+jseg)->desc;		/* el nombre del desc és el mateix que el de pseg+jseg. */
    for(k=1;k<newsg;k++) {
        (pseg2+k)->beg = (pseg+jseg+k)->beg;	/* creant tots els segments per la dreta del nou individu */
        (pseg2+k)->end = (pseg+jseg+k)->end;
        (pseg2+k)->desc = (pseg+jseg+k)->desc;	/* el descendent de seg indica l'individu de tree. IMPORTANT ! */
    }
    lsg = chrom[ic].nseg = lsg - newsg + in;	/* el nombre de segments de ic és jseg més in (1 o 0) */
    lsgm1 = lsg - 1;				/* l'últim fragment és lsgm1 */
    nlinks -= pseg2->beg - (pseg+lsgm1)->end;	
    /* posicions a recombinar: restem el principi d'un segment amb el final de l'altre: normalment resta només un nt. pero 
    els max i min canvien depenent dels segments, aixi que es poden restar molts més nt. si cau entre segments */
    
    if(ifsel_glob) { /*in case SELECTION*/
        total_nts -= (long int)pseg2->beg - (long int)(pseg+lsgm1)->end;/*restem entre zones recombinants*/
        if((long int)sel_nts_glob > (long int)(pseg+lsgm1)->end) total_nts += sel_nts_glob - (long int)(pseg+lsgm1)->end;
        /*sumem si sel_nts dreta del chrom esquerra*/
        if((long int)sel_nts_glob < (long int)pseg2->beg) total_nts += (long int)pseg2->beg - sel_nts_glob;		
        /*sumem si sel_nts esquerra del chrom dreta*/
        if(sel_nts_glob < (long int)is) 
            new_chrom = nchrom-1;/*quin es el chrom separat de sel_nts? assignar a 'new_chrom'*/
        else 
            new_chrom = ic;
    }
    
    len = (pseg+lsgm1)->end - pseg->beg;	/* conversió */
    cleft += 1.0 - pow(pc,len);    		/* conversió */
    len = (pseg2 + newsg - 1)->end - pseg2->beg;/* llargada del nou individu */
    cleft += 1.0 - pow(pc,len);   		/* conversió */
    if(!(chrom[ic].pseg = (struct seg *)realloc(chrom[ic].pseg,(unsigned long)(lsg*sizeof(struct seg)))))
        perror("realloc error. xover.3");	/* només es deixen lsg(=jseg+in) segments a ic */
    if(in) {
        begs = (long int)pseg2->beg;	/* inici del primer segment */
        for(i=0,k=0;(k < (long int)nsegs-1) && (begs > seglst[seglst[i].next].beg-1);i= seglst[i].next, k++);
        /* arribem a on es troba begs, i tenim el nombre del segment (i)-> o arribem al final, l'ultim segment */
        if(begs != seglst[i].beg) {	/* en cas que no hi hagi hagut recombinació al mateix lloc */
            /* new tree */
            if(nsegs >= seglimit) {	
                seglimit += SEGINC;
                if(!(nnodes = (long int *)realloc(nnodes,(unsigned)(sizeof(long int)*seglimit))))
                    perror("realloc error. xover.4");
                if(!(seglst = (struct segl *)realloc(seglst,(unsigned)(sizeof(struct segl)*seglimit))))
                    perror("realloc error. xover.5");
            }
            seglst[nsegs].next = seglst[i].next;	/* Crear un segment entre els altres segments */
            seglst[i].next = nsegs;			/* MOLT BO ! */
            seglst[nsegs].beg = begs;
            if(!(seglst[nsegs].ptree = (struct node *)calloc((unsigned)(2*nsam),sizeof(struct node))))
                perror("calloc error. xover.6");	/* crear els nodes de l'arbre del nou segment */
            nnodes[nsegs] = nnodes[i];			/* el nombre de nodes de nsegs és el mateix que i */
            ptree1 = seglst[i].ptree;			/* punter a l'arbre de i */
            ptree2 = seglst[nsegs].ptree;		/* punter a l'arbre de nsegs */
            nsegs++;					/* nsegs és un segment més gran */
            for(k=0;k<=nnodes[i];k++) {			/* donem els mateixos valors a i que a nsegs */
                (ptree2+k)->abv = (ptree1+k)->abv;
                (ptree2+k)->time = (ptree1+k)->time;
            }
        }
    }
}
/* coalescent functions */
void pick2_chrom(long int pop,int *config,int *pc1,int *pc2)
{
    int c1,c2,cs,cb,i,count;
    void pick2(int,int *,int *);
    double ran1(void);
    
    pick2(config[pop],&c1,&c2);	/* trobar els dos individus de la població pop que tindran coalescència, nombre c1 i c2 */
    cs = (c1 > c2) ? c2 : c1;	/* Ara buscar quins individus pertanyen a pop, i trobar dins pop quin es c1 i c2 */
    cb = (c1 > c2) ? c1 : c2;	/* cs és el mínim i cb el màxim entre c1 i c2 */
    i = count = 0;
    for(;;) {
        while(chrom[i].pop != pop) i++;
        if(count == cs) break;	/* anem buscant individus de pop, quan trobem, count++ fins cs, i ja tenim l'individu */
        count++;
        i++;
    }
    *pc1 = i;
    i++;
    count++;
    for(;;) {
        while(chrom[i].pop != pop) i++;
        if(count == cb) break;	/* igual per cb */
        count++;
        i++;
    }
    *pc2 = i;
}
void pick2(int n, int *i,int *j)
{
    double ran1(void);
   
    *i = (int)floor(n*ran1());
    while((*j = (int)floor(n * ran1())) == *i); /* dos valors aleatoris entre 0 i n-1 */
}
/* coalescent */
int ca(int nsam,unsigned long nsites,int c1, int c2)
{
    int yes1,yes2,seg1,seg2;
	long int seg;
    unsigned long start,end;
    long int tseg,desc,k;
    struct seg *pseg;
    struct node *ptree;
    int isseg(unsigned long,int,int *);
    unsigned long links(int);
    unsigned long calc_total_nts(int);
     
    seg1=0;	/* valor de la primera posició del segment actiu de l'individu c1 */
    seg2=0;	/* valor de la primera posició del segment actiu de l'individu c2 */
    
    if(!(pseg=(struct seg *)calloc((unsigned)nsegs,sizeof(struct seg))))	/* vector de segments pel nou node */
        perror("calloc error. ca.1");
    tseg = -1;						/* nombre de segments del nou node */
    
    for(seg=0,k=0;k<(long int)nsegs;seg=seglst[seg].next,k++) {	/* mirem tots els segments */
        start = seglst[seg].beg;			/* 1a posició del segment que mirem */
        yes1  = isseg(start,c1,&seg1);			/* yes1=1 si el segment es troba dins c1 */
        yes2  = isseg(start,c2,&seg2);			/* yes2=1 si el segment es troba dins c2 */
        if(yes1 || yes2) {				/* si un dels dos té el segment */
            tseg++;					/* sumem el segment */
            (pseg+tseg)->beg = seglst[seg].beg;		/* l'inici del segment pel nou node és l'inici del segment actiu */
            end = (unsigned long)(k < (long int)nsegs-1 ? (long int)seglst[seglst[seg].next].beg-1 : (long int)nsites-1);
            (pseg+tseg)->end = end;			/* i el final és el final del segment */
            
            if(yes1 && yes2) {				/* si els dos tenen el node hi ha COALESCÈNCIA */
                nnodes[seg]++;				/* per aquell segment sumem el nombre de nodes de l'arbre */
                if(nnodes[seg] >= (2*nsam-2)) tseg--;	/* en cas sigui el MRCA, aleshores restem el segment */
                else (pseg+tseg)->desc = nnodes[seg];	/*si no,el nombre del node del segment de l'individu és nnodes[seg]*/
                ptree = seglst[seg].ptree;		/* punter a l'arbre de seglst[seg] */
                desc = (chrom[c1].pseg+seg1)->desc;	/* el node de c1 és desc */
                (ptree+desc)->abv = nnodes[seg];	/* el node de desc apunta a nnodes[seg], el nou node */
                desc = (chrom[c2].pseg + seg2)->desc;	/* el node de c2 és desc */
                (ptree+desc)->abv = nnodes[seg];	/* el node de desc apunta a nnodes[seg], el nou node */
                (ptree+nnodes[seg])->time = (double)t;		/* per últim, indicar el temps de la coalescència */
            }
            else (pseg+tseg)->desc = (yes1 ? (chrom[c1].pseg + seg1)->desc : (chrom[c2].pseg + seg2)->desc);
        }	/* en cas només un dels individus té el segment, indicar el desc de l'individu que el té */
    }

    nlinks -= links(c1);	/* mida de posicions de recombinació. restar els individus que tenen coal i sumar el nou */
    if(ifsel_glob) total_nts -= (long int)calc_total_nts(c1); /*in case SELECTION*/
    cleft  -= (double)1 - pow(pc,(double)links(c1));	/* conversió */
    free(chrom[c1].pseg);	/* els segments de c1 ara ja no interessen */
    if(tseg < 0) {		/* en cas el nou node sigui MRCA per TOTS els segments que contenien els descendents */
        free(pseg);				/* en aquest cas tampoc interessa res de pseg */
        chrom[c1].pseg = chrom[nchrom-1].pseg;	/* assignem les direccions de c1 a l'últim cromosoma, c1 queda lliure */
        chrom[c1].nseg = chrom[nchrom-1].nseg;
        chrom[c1].pop  = chrom[nchrom-1].pop;
        if(c2==nchrom-1) c2 = c1;		/* c2 també queda liure, si c2 és l'últim, c2=c1 */
        nchrom--;				/* un individu menys! */
    }
    else {			/* en cas no hi hagi MRCA per tots segments */
        if(!(pseg = (struct seg *)realloc(pseg,(unsigned)((tseg+1)*sizeof(struct seg)))))
            perror("realloc error. ca.2");
        chrom[c1].pseg = pseg;			/* ara c1 passa a tenir els valors de pseg */
        chrom[c1].nseg = tseg + 1;
        nlinks += links(c1);
        if(ifsel_glob) total_nts += (long int)calc_total_nts(c1); /*in case SELECTION*/
        cleft += 1.0 - pow(pc,(double)links(c1)); /* conversió */
    }

    nlinks -= links(c2);	
    if(ifsel_glob) total_nts -= (long int)calc_total_nts(c2); 	/*in case SELECTION*/
    cleft  -= 1.0 - pow(pc,(double)links(c2)); /* conversió */
    free(chrom[c2].pseg);			/* eliminem c2 i apuntem a l'últim individu */
    chrom[c2].pseg = chrom[nchrom-1].pseg;
    chrom[c2].nseg = chrom[nchrom-1].nseg;
    chrom[c2].pop  = chrom[nchrom-1].pop;
    nchrom--;					/* un individu menys! */
    if(tseg<0) return(2); /* decrease of nchrom is two */
    else return(1);
}

/* isseg: does chromosome c contain the segment on seglst which starts at start? *psg is the segment of chrom[c] at wich one is to begin looking. */
int isseg(unsigned long start, int c, int *psg) /**psg és inicialment 0, però varia quan anem avançant pels segments */
{
    int ns;
    struct seg *pseg;
    
    ns = chrom[c].nseg;		/* nombre de segments de l'individu c */
    pseg = chrom[c].pseg;	/* punter al primer segment de l'individu c */
    
    /*Sylvain diu que es incorrecte, aleshores: for(;((*psg) < ns) && ((pseg+(*psg))->beg <= start);  ++(*psg)) */
    /*for(;((pseg+(*psg))->beg <= start) && ((*psg) < ns); ++(*psg))*//*des del segment psg fins que sigui més gran de start*/
    for(;((*psg) < ns) && ((pseg+(*psg))->beg <= start);  ++(*psg)) /*Sylvain modification*/
        if((pseg+(*psg))->end >= start) return(1);	/* psg ja és 1+ perque ++(*psg), i no (*psg)++ ?? NO POT SER */
        /* en cas el final del segment sigui més gran o igual start, tenim el segment buscat a l'individu c */
        /* com c pot tenir segments més grans, start ha d'estar entre o = begin i end, aleshores el segment és inclós */
    return(0);
}
unsigned long calc_total_nts(int c)	/*IN CASE SELECTION*/
{
    unsigned long links(int);
    long int len;
    int ns;
    len = (long int)links(c);					/*afegim regio estudiada*/
    ns = chrom[c].nseg - 1;
    if((long int)sel_nts_glob > (long int)(chrom[c].pseg + ns)->end) 
        len += sel_nts_glob - (long int)(chrom[c].pseg + ns)->end;/*sumem si sel_nts dreta del chrom esquerra*/
    if((long int)sel_nts_glob < (long int)(chrom[c].pseg)->beg) 
        len += (long int)(chrom[c].pseg)->beg - sel_nts_glob;	/*sumem si sel_nts esquerra del chrom dreta*/
    return (unsigned long)len;
}

unsigned long links(int c)
{
    int ns;
    ns = chrom[c].nseg - 1;
    return((chrom[c].pseg + ns)->end - (chrom[c].pseg)->beg);	/* max - min de l'individu c */
}

/* No tinc clar encara si pot treballar amb seleccio.... */
int cleftr(int nsam)
{
    struct seg *pseg ;
    int ic;
    double x, sum;
    unsigned long len,is;
    void xover(int, int, unsigned long);
    unsigned long links(int);
    double ran1(void);

    while((x = (double)cleft*ran1()) == 0.0) ;
    sum = 0.0;
    ic = -1;

    while(sum < x) {
        sum += (double)1 - (double)pow((double)pc,links(++ic)); /*trobem l'individu ic*/
    }
    pseg = chrom[ic].pseg; 	/*pseg apunta a l'individu ic*/
    len = links(ic);		/*mirem la llargada del max al min de ic*/
    is = pseg->beg + (unsigned long)floor(1.0 + (double)log(1.0 - (1.0- pow( pc, len))*ran1())/lnpc)-1; /*localitzem el punt de rec a is*/
    xover(nsam,ic,is);		/*recombinacio*/
    return(ic);			/*a ic*/
}
/* No tinc clar encara si pot treballar amb seleccio.... */
int cinr(int nsam, unsigned long nsites)
{
    struct seg *pseg ;
    int lsg, lsgm1,ic;
    unsigned long spot,len,el,is,endic;
    
    void xover(int, int, unsigned long); 
    int ca(int,unsigned long,int,int) ;
    double ran1(void);

    spot = (unsigned long)floor(nlinks*ran1() + (double)1); /*el punt escollit a tota la regio*/
    /* get chromosome number (ic) */
    for(ic=0;ic<nchrom;ic++) { 	/* Busca els valors MÀXIM i MÍNIM de l'individu escollit. (Què nt. segment i individu) */
        lsg = chrom[ic].nseg; /* el nombre de segments a l'individu ic */
        lsgm1 = lsg - 1;/* de fet lsg-1 conté la informació de l'ultim segment */
        pseg = chrom[ic].pseg;/* punter al primer segment de l'individu ic */
        el = ((pseg+lsgm1)->end) - (pseg->beg);/* mida dels segments (max-min) a l'individu ic*/
        if(spot <= el) break;/* anem restant 'el' fins trobar l'individu */
        spot -= el ;
    }
    is = pseg->beg + spot -1;			/* posició dins l'individu ic */
    endic = (pseg+lsgm1)->end; 			/*posicio final de l'individu ic*/
    xover(nsam,ic,is);

    len = (unsigned long)floor(1.0 + (double)log(ran1())/lnpc);	/*llargada del evente de conversio...*/
    if(is+len >= endic) return(ic);  		/*si es mes llarg, es igual que rec, acabem*/
    if(is+len < (chrom[nchrom-1].pseg)->beg){	/*si es mes curt que l'inici del nou chrom*/
        ca(nsam,nsites,ic,nchrom-1);			/*aleshores eliminem el nou chrom, coalescencia*/
        return(-1);					/*i no ha passat res (...?)*/
    }
    xover(nsam,nchrom-1,is+len);		/*tornem a recombinar el fragment a llargada is+len*/
    ca(nsam,nsites,ic,nchrom-1);		/*... i fem coalescencia a l'extrem. conversio finalitzada*/
    return(ic);					/*i tenim un mes*/
}

/*MORE FUNCTIONS*/

double functiont_freqp_sel(double x,double tcoal,double t,double ts,double sinit,double eps,double pop_sel,double no1,double no2,double no3)
/*include no1 and no2 to have all functions with the same number of variables*/
{
	double fdT;

	no1=no2=no3=(double)0;
	fdT  = ((eps*(double)exp(pop_sel*(t+       x-sinit-(ts)))/pop_sel)+(t+       x)*((double)1-eps))/((double)1-eps);
	fdT -= ((eps*(double)exp(pop_sel*(t+(double)0-sinit-(ts)))/pop_sel)+(t+(double)0)*((double)1-eps))/((double)1-eps);
	fdT -= tcoal;
	
	return fdT;
}

double functiont_freqq_nsel(double x,double tcoal,double t,double ts,double sinit,double eps,double pop_sel,double no1,double no2,double no3)
/*include no1 and no2 to have all functions with the same number of variables*/
{
	double fdT;
	
	no1=no2=no3=(double)0;
	fdT  = (((eps-(double)1)*(double)exp(-pop_sel*(t+       x-sinit-(ts)))/pop_sel)+(t+       x)*eps)/eps;
	fdT -= (((eps-(double)1)*(double)exp(-pop_sel*(t+(double)0-sinit-(ts)))/pop_sel)+(t+(double)0)*eps)/eps;
	fdT -= tcoal;
	
	return fdT;
}

double functiont_logistic(double x,double tcoal,double t,double tpast0,double tpast1,double alphag,double nrec1,double npast,double nrecf,double Ts)
{
	double fdT,expw0,expw1,logo0,logo1;
	
	if(fabs(nrec1-npast) < (double)1E-08) {
		fdT = nrecf * (((t+       x)/npast)) - nrecf * (((t+(double)0)/npast)) - tcoal;
	}
	else {
		expw0 = -alphag*((t+Ts+       x)-(tpast0+tpast1)/(double)2);
		expw1 = -alphag*((t+Ts+(double)0)-(tpast0+tpast1)/(double)2);
		if(expw0 > (double)60) {
			logo0 = expw0;
			logo1 = expw1;
		}
		else {
			logo0 = (double)log(npast+(double)exp(expw0)*nrec1);
			logo1 = (double)log(npast+(double)exp(expw1)*nrec1);
		}
	
		fdT  =  (nrecf * (((t+       x)/npast) + ((nrec1-npast)*logo0)/(alphag*npast*nrec1)));
		fdT -=  (nrecf * (((t+(double)0)/npast) + ((nrec1-npast)*logo1)/(alphag*npast*nrec1)));
		fdT -=  tcoal;
	}
		
	return fdT;
}

int zbracn(double (*func)(double,double,double,double,double,double,double,double,double,double),double *x1,double *x2,double a0,double a1,double a2,double a3,double a4,double a5,double a6,double a7,double a8)
{
	/* Based on Numerical Recipes in C. Press et al. 1992. 
    We need *x1 and *x2 be the range where a root is within them. We expand geometrically the range until finding
	(one a positive and one a negative value). If not, return 0.
	*/

    double f1,f2;
    int k=60;
    
    if(*x1 == *x2) return 0;

    f1 = (*func)(*x1,a0,a1,a2,a3,a4,a5,a6,a7,a8);
    f2 = (*func)(*x2,a0,a1,a2,a3,a4,a5,a6,a7,a8);
	
	if(f1*f2 < (double)0) return 1;

    while(k--) {
        if(fabs(f1) < fabs(f2)) {
            *x1 += (double)1.5 * (*x1 - *x2);
            f1 = (*func)(*x1,a0,a1,a2,a3,a4,a5,a6,a7,a8);
        }
        else {
            *x2 += (double)1.5 * (*x2 - *x1);
            f2 = (*func)(*x2,a0,a1,a2,a3,a4,a5,a6,a7,a8);
        }
        if(f1*f2 < (double)0) return 1;
    }
    return 0;
}

double zriddrn(double (*func)(double,double,double,double,double,double,double,double,double,double),double xlow,double xhigh,double xacc,double a0,double a1,double a2,double a3,double a4,double a5,double a6,double a7, double a8)
{
	/* Based on Numerical Recipes in C. Press et al. 1992., p. 358 an on
	Ridders, 1979, IEEE Transactions on Circuits and systems, Vol. Cas-26, No. 11, pp. 979-980.
	*/
    int k=60;
	double flow,fhigh;
	double f1,f2,f3,f4;
	double x1,x2,x3,x4;
	double den,num,nsign;
	

    flow  = (*func)(xlow ,a0,a1,a2,a3,a4,a5,a6,a7,a8);
    fhigh = (*func)(xhigh,a0,a1,a2,a3,a4,a5,a6,a7,a8);

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
		f3 = (*func)(x3,a0,a1,a2,a3,a4,a5,a6,a7,a8);
		if(f1 - f2 < (double)0) nsign = (double)-1;
		else nsign = (double)1;
		num = (x3-x1) * f3 * nsign;
		den = (double)sqrt((double)f3*(double)f3 - (double)f1*(double)f2);
		if(den <= xacc && -den <= xacc) return x3;
		x4 = x3 + num/den;
		f4 = (*func)(x4,a0,a1,a2,a3,a4,a5,a6,a7,a8);
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


