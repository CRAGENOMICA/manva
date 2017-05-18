/*LINK TO HUDSON'S MS PROGRAM*/


#include "mhmlspmlnt.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int mainms (struct var **data,struct var2b **inputms,struct statistisim ***matrixsim,struct statistisimmuloc **matrixmlsim,struct horizontalstatsml **avgstatloci,int onlymulo,FILE *file_output)
{	    
    int ms(struct var2 *,struct statistisim **,struct statistisimmuloc **,struct horizontalstatsml **,int,FILE *);
    int dohka_sim(struct statistisim ***,struct statistisimmuloc **,double *, int *,long int);
    int compare_(const void *,const void *);
    void init_seed1(long int);
    
    int v,x,z;
    long int i;
    static int y=0;
    static double *vecsort;
	double valmed;
    
    static struct var2 *inputp = 0;

    /*seed*/
    init_seed1(data[0][0].seed1);

    /*Build var2*/
    if(inputp == 0) {
        /*calloc inputp*/
        if((inputp = (struct var2 *) calloc(1,sizeof(struct var2))) == 0) {
            puts("\nError: memory not reallocated for running simulations. menu_coal_sim.0 \n");
            if(file_output)
                fputs("\nError: memory not reallocated for running simulations. menu_coal_sim.0 \n",file_output);
            return 0;
        }
        if((inputp[0].config = (int *)calloc(1,sizeof(int))) == 0) {
            puts("Error: memory not reallocated. main_sm2.7");
            return 0;
        }
        if((inputp[0].factor_pop = (double *)calloc(1,sizeof(double))) == 0) {
            puts("Error: memory not reallocated. main_sm2.8");
            return 0;
        }
        if((inputp[0].nrec = (double *)calloc(3,sizeof(double))) == 0) {
            puts("Error: memory not reallocated. main_sm2.9");
            return 0;
        }
        if((inputp[0].npast = (double *)calloc(3,sizeof(double))) == 0) {
            puts("Error: memory not reallocated. main_sm2.10");
            return 0;
        }
        if((inputp[0].tpast = (double *)calloc(3,sizeof(double))) == 0) {
            puts("Error: memory not reallocated. main_sm2.11");
            return 0;
        }
        if((inputp[0].freq = (double *)calloc(1,sizeof(double))) == 0) {
            puts("Error: memory not reallocated. main_sm2.12");
            return 0;
        }
    }

    /*init matrixmlsim*/
    for(i=0;i<(long int)data[0][0].n_iter;i++) {
        matrixmlsim[0][i].Snsites = 0;
        matrixmlsim[0][i].Sbiallsites = 0;
        matrixmlsim[0][i].S2biallsites = 0;
        matrixmlsim[0][i].Sfixed = 0;
        matrixmlsim[0][i].S2fixed = 0;
        matrixmlsim[0][i].Snhapl = 0;
        matrixmlsim[0][i].S2nhapl = 0;
        matrixmlsim[0][i].Snhaplsam = 0;
        matrixmlsim[0][i].S2nhaplsam = 0;
        matrixmlsim[0][i].Shapldiv = (double)0;
        matrixmlsim[0][i].S2hapldiv = (double)0;
        matrixmlsim[0][i].Sewtest = (double)0;
        matrixmlsim[0][i].S2ewtest = (double)0;
        matrixmlsim[0][i].Sza = (double)0;
        matrixmlsim[0][i].S2za = (double)0;
        matrixmlsim[0][i].Sb = 0;
        matrixmlsim[0][i].S2b = 0;
        matrixmlsim[0][i].Sq = 0;
        matrixmlsim[0][i].S2q = 0;
        matrixmlsim[0][i].SRvpi = (double)0;
        matrixmlsim[0][i].S2Rvpi = (double)0;
        matrixmlsim[0][i].SRm = (long int)0;
        matrixmlsim[0][i].S2Rm = (long int)0;
    
        matrixmlsim[0][i].Stheta_wat = (double)0;
        matrixmlsim[0][i].S2theta_wat = (double)0;
        matrixmlsim[0][i].Stheta_wat_nut = (double)0;
        matrixmlsim[0][i].S2theta_wat_nut = (double)0;
        matrixmlsim[0][i].Stheta_taj = (double)0;
        matrixmlsim[0][i].S2theta_taj = (double)0;
        matrixmlsim[0][i].Stheta_taj_nut = (double)0;
        matrixmlsim[0][i].S2theta_taj_nut = (double)0;
        matrixmlsim[0][i].Stheta_fw = (double)0;
        matrixmlsim[0][i].S2theta_fw = (double)0;
        matrixmlsim[0][i].Stheta_fw_nut = (double)0;
        matrixmlsim[0][i].S2theta_fw_nut = (double)0;
        matrixmlsim[0][i].Stheta_fuli = (double)0;
        matrixmlsim[0][i].S2theta_fuli = (double)0;
        matrixmlsim[0][i].Stheta_L = (double)0;
        matrixmlsim[0][i].S2theta_L = (double)0;
        matrixmlsim[0][i].Stheta_fuli_nut = (double)0;
        matrixmlsim[0][i].S2theta_fuli_nut = (double)0; 
        matrixmlsim[0][i].Stheta_fulin = (double)0;
        matrixmlsim[0][i].S2theta_fulin = (double)0;
        matrixmlsim[0][i].Stheta_fulin_nut = (double)0;
        matrixmlsim[0][i].S2theta_fulin_nut = (double)0; 
        matrixmlsim[0][i].Stheta_L_nut = (double)0;
        matrixmlsim[0][i].S2theta_L_nut = (double)0; 
        matrixmlsim[0][i].Sndivergence = (double)0;
        matrixmlsim[0][i].S2ndivergence = (double)0;
        matrixmlsim[0][i].Sndivergence_nut = (double)0;
        matrixmlsim[0][i].S2ndivergence_nut = (double)0; 
        matrixmlsim[0][i].Sndivergence = (double)0;
        matrixmlsim[0][i].S2ndivergence = (double)0;
        matrixmlsim[0][i].Sndivergence_nut = (double)0;
        matrixmlsim[0][i].S2ndivergence_nut = (double)0; 
    
        matrixmlsim[0][i].StajimaD = (double)0;
        matrixmlsim[0][i].S2tajimaD = (double)0;
        matrixmlsim[0][i].nltajD = (int)0;
        matrixmlsim[0][i].SfuliD = (double)0;
        matrixmlsim[0][i].S2fuliD = (double)0;
        matrixmlsim[0][i].nlflD = (int)0;
        matrixmlsim[0][i].SfuliDn = (double)0;
        matrixmlsim[0][i].S2fuliDn = (double)0;
        matrixmlsim[0][i].nlflDn = (int)0;
        matrixmlsim[0][i].SfuliF = (double)0;
        matrixmlsim[0][i].S2fuliF = (double)0;
        matrixmlsim[0][i].nlflF = (int)0;
        matrixmlsim[0][i].SfuliFn = (double)0;
        matrixmlsim[0][i].S2fuliFn = (double)0;
        matrixmlsim[0][i].nlflFn = (int)0;
        matrixmlsim[0][i].SfuFs = (double)0;
        matrixmlsim[0][i].S2fuFs = (double)0;
        matrixmlsim[0][i].nlFs = (int)0;
        matrixmlsim[0][i].SfaywuH = (double)0;
        matrixmlsim[0][i].S2faywuH = (double)0;
        matrixmlsim[0][i].SfaywuHo = (double)0;
        matrixmlsim[0][i].S2faywuHo = (double)0;
        matrixmlsim[0][i].SzengE = (double)0;
        matrixmlsim[0][i].S2zengE = (double)0;
        matrixmlsim[0][i].nlH = (int)0;
        matrixmlsim[0][i].nlHo = (int)0;
        matrixmlsim[0][i].nlE = (int)0;
        matrixmlsim[0][i].SrZA = (double)0;
        matrixmlsim[0][i].S2rZA = (double)0;
        matrixmlsim[0][i].nlZ = (int)0;
        matrixmlsim[0][i].SwB = (double)0;
        matrixmlsim[0][i].S2wB = (double)0;
        matrixmlsim[0][i].nlB = (int)0;
        matrixmlsim[0][i].SwQ = (double)0;
        matrixmlsim[0][i].S2wQ = (double)0;
        matrixmlsim[0][i].nlQ = (int)0;
        matrixmlsim[0][i].SR2 = (double)0;
        matrixmlsim[0][i].S2R2 = (double)0;
        matrixmlsim[0][i].nlR2 = (int)0;
    }

    printf("\n Each dot indicate that approximately 2%% of the simulation is done.");
    printf("\n         1    2    3    4    5    6    7    8    9  100%%");
    printf("\n RUN ");
    fflush(stdout);
    if(file_output) {
        fprintf(file_output,"\n Each dot indicate that approximately 2%% of the simulation is done.");
        fprintf(file_output,"\n         1    2    3    4    5    6    7    8    9  100%%");
        fprintf(file_output,"\n RUN ");
        fflush(file_output);
    }
    /* COALESCENT SIMULATIONS FOR EACH INDEPENDENT LOCUS */
    for(x=0;x<data[0][0].n_loci;x++) {    
        inputp[0].howmany = data[0][0].n_iter;
        inputp[0].nsam = inputms[0][x].nsam;
        inputp[0].nsites = inputms[0][x].nsites + inputms[0][x].nmhits;
        inputp[0].factor_chrn = (double)1/inputms[0][x].factor_chrn; /*REMEMBER: We must do 1/factor_chrn !*/
		/*recombination*/
		inputp[0].r = inputms[0][x].r * (double)1/inputms[0][x].factor_chrn;
		if(inputms[0][x].factor_chrn == (double)1) 
			if(data[0][0].no_rec_males == 1) inputp[0].r = inputp[0].r * (double)0.5;
		if(inputms[0][x].factor_chrn == (double)0.75) 
			inputp[0].r = inputp[0].r * (double)2/(double)3;
		if(inputms[0][x].factor_chrn <= (double)0.5) 
			inputp[0].r = (double)0;
        /**/
		if(inputms[0][x].subdivision == 1) {
            inputp[0].theta = inputms[0][x].theta/(double)inputms[0][x].npop * (double)1/inputms[0][x].factor_chrn;
            inputp[0].r = inputms[0][x].r/(double)inputms[0][x].npop;
            inputp[0].T_out = data[0][0].time_spec * (double)inputms[0][x].npop;
        }
        else {
            inputp[0].theta = inputms[0][x].theta * (double)1/inputms[0][x].factor_chrn;
            inputp[0].T_out = data[0][0].time_spec;
        }
		inputp[0].T_out /= 2.0; /*time in 4N generations*/
        inputp[0].nloci = x;
        inputp[0].tloci = data[0][0].n_loci;
        if(data[0][0].time_spec >= (double)0) inputp[0].mhits = 1;
        else inputp[0].mhits = 0;
        inputp[0].ratio_sv = inputms[0][x].ratio_sv;
        
        if(inputms[0][x].subdivision == 0 && inputms[0][x].split_pop == 0) {
            inputp[0].npop = 1;
            inputp[0].migrate = (double)0;
            if((inputp[0].config = (int *)realloc(inputp[0].config,sizeof(int))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.7");
                return 0;
            }
            if((inputp[0].factor_pop = (double *)realloc(inputp[0].factor_pop,sizeof(double))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.8");
                return 0;
            }
            inputp[0].config[0] = inputp[0].nsam;
            inputp[0].factor_pop[0] = (double)1.; 
        }
        else {
            inputp[0].npop = inputms[0][x].npop;
            inputp[0].migrate = inputms[0][x].mig_rate;
            if((inputp[0].config = (int *)realloc(inputp[0].config,inputms[0][x].npop*sizeof(int))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.7");
                return 0;
            }
            if((inputp[0].factor_pop = (double *)
                realloc(inputp[0].factor_pop,inputms[0][x].npop*sizeof(double))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.8");
                return 0;
            }
            for(v=0;v<inputms[0][x].npop;v++) {
                if(v==0) {
                    inputp[0].config[v] = inputp[0].nsam;
                    inputp[0].factor_pop[v] = (double)1.0;
                }
                else {
                    inputp[0].config[v] = 0;
                    inputp[0].factor_pop[v] = inputms[0][x].factor_pop[v];
                }
            }
        }
        
        if(inputms[0][x].changesizepop == 0) {
            inputp[0].nintn = 0;
            if((inputp[0].nrec = (double *)realloc(inputp[0].nrec,2*sizeof(double))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.7");
                return 0;
            }
            if((inputp[0].npast = (double *)realloc(inputp[0].npast,2*sizeof(double))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.7");
                return 0;
            }
            if((inputp[0].tpast = (double *)realloc(inputp[0].tpast,2*sizeof(double))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.7");
                return 0;
            }
            inputp[0].nrec[0] = (double)1.;
            inputp[0].npast[0] = (double)1.;
            inputp[0].tpast[0] = (double)0;
            inputp[0].nrec[1] = (double)1.;
            inputp[0].npast[1] = (double)1.;
            inputp[0].tpast[1] = (double)1.;
            inputp[0].iflogistic = 1;
            inputp[0].ts = (double)0;
        }
        else {
            inputp[0].nintn = inputms[0][x].nintn;
			if((inputp[0].nrec = (double *)realloc(inputp[0].nrec,(2+inputms[0][x].nintn)*sizeof(double))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.7");
                return 0;
            }
            if((inputp[0].npast = (double *)realloc(inputp[0].npast,(2+inputms[0][x].nintn)*sizeof(double))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.7");
                return 0;
            }
            if((inputp[0].tpast = (double *)realloc(inputp[0].tpast,(2+inputms[0][x].nintn)*sizeof(double))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.7");
                return 0;
            }
            inputp[0].iflogistic = 1;
			inputp[0].ts = inputms[0][0].ts;

            inputp[0].nrec[0] = (double)1.;
            inputp[0].npast[0] = (double)1.;
            inputp[0].tpast[0] = (double)0;
            inputp[0].nrec[1] = (double)1.;
            inputp[0].npast[1] = (double)1.;
            inputp[0].tpast[1] = (double)1.;

            for(v=1;v<=inputms[0][x].nintn;v++) {
                inputp[0].nrec[v] = inputms[0][x].nrec[v-1];
                inputp[0].npast[v] = inputms[0][x].npast[v-1];
                inputp[0].tpast[v] = inputms[0][x].tpast[v-1] + inputp[0].tpast[v-1];
            }
        }
        
        if(inputms[0][x].split_pop == 0) {
            inputp[0].split_pop = 0;
            inputp[0].time_split = (double)0;
            inputp[0].time_scoal = (double)0;
            inputp[0].factor_anc = (double)1.;
            if((inputp[0].freq = (double *)realloc(inputp[0].freq,sizeof(double))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.7");
                return 0;
            }
            inputp[0].freq[0] = 1;
        }
        else {
            inputp[0].split_pop = inputms[0][x].split_pop;
            inputp[0].time_split = inputms[0][x].time_split;
            inputp[0].time_scoal = inputms[0][x].time_scoal;
            inputp[0].factor_anc = inputms[0][x].factor_anc;
            if((inputp[0].freq = (double *)realloc(inputp[0].freq,inputms[0][x].npop*sizeof(double))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.7");
                return 0;
            }
            for(v=0;v<inputms[0][x].npop;v++) {
                inputp[0].freq[v] = inputms[0][x].freq[v];
            }
        }
        if(inputms[0][x].pselection == 0) {
            inputp[0].ifselection = 0;
            inputp[0].pop_size = 0;
            inputp[0].pop_sel = (double)0;
            inputp[0].sel_nt = 0;
            inputp[0].sinit = (double)0;
        }
        else {
            if((inputp[0].config = (int *)realloc(inputp[0].config,2*sizeof(int))) == 0) {
                puts("Error: memory not reallocated. menu_coal_sim.7");
                return 0;
            }
            inputp[0].npop = 2;
            inputp[0].ifselection = 1;
            inputp[0].pop_size = inputms[0][x].pop_size;
            inputp[0].pop_sel = inputms[0][x].pop_sel;
            inputp[0].sel_nt = (double)inputms[0][x].sel_nt;
            inputp[0].sinit = inputms[0][x].sinit;
        }
        inputp[0].tlimit = data[0][0].tlimit;

        if(onlymulo == 1) {
            /*init avgstatloci*/
            avgstatloci[0][x].biallsites = 0;
            avgstatloci[0][x].fixed = 0;
            avgstatloci[0][x].nhapl = 0;
            avgstatloci[0][x].nhaplsam = 0;
            avgstatloci[0][x].hapldiv = (double)0;
            avgstatloci[0][x].ewtest = (double)0;
            avgstatloci[0][x].za = (double)0;
            avgstatloci[0][x].b = 0;
            avgstatloci[0][x].q = 0;
            avgstatloci[0][x].Rvpi = (double)0;
            avgstatloci[0][x].Rm = (double)0;
            avgstatloci[0][x].theta_wat = (double)0;
            avgstatloci[0][x].theta_taj = (double)0;
            avgstatloci[0][x].theta_fw = (double)0;
            avgstatloci[0][x].theta_fuli = (double)0;
            avgstatloci[0][x].theta_fulin = (double)0;
            avgstatloci[0][x].theta_L = (double)0;
            avgstatloci[0][x].ndivergence = (double)0;
            avgstatloci[0][x].theta_wat_nut = (double)0;
            avgstatloci[0][x].theta_taj_nut = (double)0;
            avgstatloci[0][x].theta_fw_nut = (double)0;
            avgstatloci[0][x].theta_fuli_nut = (double)0;
            avgstatloci[0][x].theta_fulin_nut = (double)0;
            avgstatloci[0][x].theta_L_nut = (double)0;
            avgstatloci[0][x].ndivergence_nut = (double)0;
            avgstatloci[0][x].tajimaD = (double)0;
            avgstatloci[0][x].fuliD = (double)0;
            avgstatloci[0][x].fuliDn = (double)0;
            avgstatloci[0][x].fuliF = (double)0;
            avgstatloci[0][x].fuliFn = (double)0;
            avgstatloci[0][x].fuFs = (double)0;
            avgstatloci[0][x].faywuH = (double)0;
            avgstatloci[0][x].faywuHo = (double)0;
            avgstatloci[0][x].zengE = (double)0;
            avgstatloci[0][x].rZA = (double)0;
            avgstatloci[0][x].wB = (double)0;
            avgstatloci[0][x].wQ = (double)0;
            avgstatloci[0][x].R2 = (double)0;
        }
        /***************** Hudson's coalescent (some few modifications) *****************/
        if(ms(inputp,*matrixsim,matrixmlsim,avgstatloci,onlymulo,file_output)) return 0;
    }
        
    printf(" Simulation finished.\n");
    fflush(stdout);
    if(file_output) {
        fprintf(file_output," Simulation finished\n");
        fflush(file_output);
    }

    if(onlymulo == 0) {
        /*calculate hka for each iteration once simulations have finished*/
        if(data[0][0].time_spec > (double)0 && data[0][0].n_loci > 1) {
			printf("\n Calculating HKA test...");
			fflush(stdout);
             for(i=0;i<(long int)data[0][0].n_iter;i++) {
                if((dohka_sim(matrixsim,matrixmlsim,&data[0][0].time_spec,&data[0][0].n_loci,i)) == 0) {
                    /*We decide to eliminate those cases we can not count*/
                    /**/
					for(v=0;v<data[0][0].n_loci;v++) {
                        matrixsim[0][v][i].hka = (double)-10000;
                        matrixsim[0][v][i].hka_theta = (double)-10000;
                    }
                    matrixmlsim[0][i].Shka = (double)-10000;
                    matrixmlsim[0][i].hka_T = (double)-10000;
					/**/
                    /*
                    puts("Error: hka could not be calculated.\n");
                    if(file_output) fputs("Error: hka could not be calculated.\n",file_output);
                    return 0;
                    */
                }
            }
            printf("\n HKA calculated succesfully.");
            if(file_output) fputs("\n HKA calculated succesfully.",file_output);
        }
    }

    /*in case onlymulo == 0 then we calculate the median for the statistics*/
    if(onlymulo == 0) {
        printf("\n Calculating median...\n");
		fflush(stdout);
        if(y==0) {
            if((vecsort = (double *)calloc(1,sizeof(double))) == 0) {
                puts("Error in alloc. main_sm");
                return 0;
            }
            y=1;
        }
        if((vecsort = (double *)realloc(vecsort,inputp[0].howmany*sizeof(double))) == 0) {
            puts("Error in alloc. main_sm");
            return 0;
        }
        for(x=0;x<inputp[0].tloci;x++) {
            for(z=0;z<37;z++){
                /*introduce all values in the vector*/
                switch(z) {
                    case 0:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].biallsites;
                        break;
                    case 1:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].fixed;
                        break;
                    case 2:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].nhapl;
                        break;
                    case 3:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].za;
                        break;
                    case 4:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].b;
                        break;
                    case 5:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].q;
                        break;
                    case 6:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].theta_wat * ((double)1/inputms[0][x].factor_chrn);
                        break;
                    case 7:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].theta_taj * ((double)1/inputms[0][x].factor_chrn);
                        break;
                    case 8:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].theta_fw * ((double)1/inputms[0][x].factor_chrn);
                        break;
                    case 9:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].theta_fuli * ((double)1/inputms[0][x].factor_chrn);
                        break;
                    case 10:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].theta_fulin * ((double)1/inputms[0][x].factor_chrn);
                        break;
                    case 11:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].ndivergence;
                        break;
                    case 12:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].tajimaD;
                        break;
                    case 13:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].fuliD;
                        break;
                    case 14:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].fuliDn;
                        break;
                    case 15:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].fuliF;
                        break;
                    case 16:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].fuliFn;
                        break;
                    case 17:
                        for(i=0;i<(long int)inputp[0].howmany;i++) {
                            if((double)matrixsim[0][x][i].fuFs == (double)-10000)
								vecsort[i] = (double)-1e6;
							else vecsort[i] = (double)matrixsim[0][x][i].fuFs;							
						}
                        break;
                    case 18:
                        for(i=0;i<(long int)inputp[0].howmany;i++) {
                            if((double)matrixsim[0][x][i].faywuH == (double)-10000)
								vecsort[i] = (double)-1e6;
                            else vecsort[i] = (double)matrixsim[0][x][i].faywuH;
						}
                        break;
                    case 19:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].rZA;
                        break;
                    case 20:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].wB;
                        break;
                    case 21:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].wQ;
                        break;
                    case 22:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].R2;
                        break;
                    case 23:
                        for(i=0;i<(long int)inputp[0].howmany;i++) {
                            vecsort[i] = (double)matrixsim[0][x][i].Rvpi;
							if(vecsort[i] != (double)-10000) vecsort[i] *= ((double)1/inputms[0][x].factor_chrn);
						}
                        break;
                    case 24:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].nhaplsam;
                        break;
                    case 25:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].hapldiv;
                        break;
                    case 26:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].ewtest;
                        break;
                    case 27:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].theta_L;
                        break;
                    case 28:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].faywuHo;
                        break;
                    case 29:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].zengE;
                        break;
					case 30:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].theta_wat_nut * ((double)1/inputms[0][x].factor_chrn);
                        break;
					case 31:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].theta_taj_nut * ((double)1/inputms[0][x].factor_chrn);
                        break;
					case 32:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].theta_fw_nut * ((double)1/inputms[0][x].factor_chrn);
                        break;
					case 33:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].theta_fuli_nut * ((double)1/inputms[0][x].factor_chrn);
                        break;
					case 34:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].theta_fulin_nut * ((double)1/inputms[0][x].factor_chrn);
                        break;
					case 35:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].ndivergence_nut;
                        break;
					case 36:
                        for(i=0;i<(long int)inputp[0].howmany;i++)
                            vecsort[i] = (double)matrixsim[0][x][i].theta_L_nut;
                        break;
                }
                if(inputp[0].howmany > 1)
                    qsort(vecsort,inputp[0].howmany,sizeof(double),compare_);/*all values are sorted*/
                /*look for all values with -10000. Substract and then find the median for the rest of values*/
                i=0;
                while(i<(long int)inputp[0].howmany && (vecsort[i] == (double)-10000 || vecsort[i] == (double)-1e6)) i++;
                /*assign the value to avgstatloci*/
                if(i>=(long int)inputp[0].howmany) {
					valmed = -10000;
				}
				else {
					if((double)(inputp[0].howmany-i)/(double)2 == (long int)((double)(inputp[0].howmany-i)/(double)2)) {
						valmed = (vecsort[(long int)floor(((double)(inputp[0].howmany-i)/(double)2.+(double)i)) - 1] + 
								  vecsort[(long int)floor(((double)(inputp[0].howmany-i)/(double)2.+(double)i))    ])/(double)2;
					}
					else {
						valmed = vecsort[(long int)floor(((double)(inputp[0].howmany-i)/(double)2.+(double)i))];
					}
				}
				switch(z) {
					case 0:
						avgstatloci[0][x].biallsites = valmed;
						break;
					case 1:
						avgstatloci[0][x].fixed = valmed;
						break;
					case 2:
						avgstatloci[0][x].nhapl = valmed;
						break;
					case 3:
						avgstatloci[0][x].za = valmed;
						break;
					case 4:
						avgstatloci[0][x].b = valmed;
						break;
					case 5:
						avgstatloci[0][x].q = valmed;
						break;
					case 6:
						avgstatloci[0][x].theta_wat = valmed;
						break;
					case 7:
						avgstatloci[0][x].theta_taj = valmed;
						break;
					case 8:
						avgstatloci[0][x].theta_fw = valmed;
						break;
					case 9:
						avgstatloci[0][x].theta_fuli = valmed;
						break;
					case 10:
						avgstatloci[0][x].theta_fulin = valmed;
						break;
					case 11:
						avgstatloci[0][x].ndivergence = valmed;
						break;
					case 12:
						avgstatloci[0][x].tajimaD = valmed;
						break;
					case 13:
						avgstatloci[0][x].fuliD = valmed;
						break;
					case 14:
						avgstatloci[0][x].fuliDn = valmed;
						break;
					case 15:
						avgstatloci[0][x].fuliF = valmed;
						break;
					case 16:
						avgstatloci[0][x].fuliFn = valmed;
						break;
					case 17:
						avgstatloci[0][x].fuFs = valmed;
						break;
					case 18:
						avgstatloci[0][x].faywuH = valmed;
						break;
					case 19:
						avgstatloci[0][x].rZA = valmed;
						break;
					case 20:
						avgstatloci[0][x].wB = valmed;
						break;
					case 21:
						avgstatloci[0][x].wQ = valmed;
						break;
					case 22:
						avgstatloci[0][x].R2 = valmed;
						break;
					case 23:
						avgstatloci[0][x].Rvpi = valmed;
						break;
					case 24:
						avgstatloci[0][x].nhaplsam = valmed;
						break;
					case 25:
						avgstatloci[0][x].hapldiv = valmed;
						break;
					case 26:
						avgstatloci[0][x].ewtest = valmed;
						break;
					case 27:
						avgstatloci[0][x].theta_L = valmed;
						break;
					case 28:
						avgstatloci[0][x].faywuHo = valmed;
						break;
					case 29:
						avgstatloci[0][x].zengE = valmed;
						break;
					case 30:
						avgstatloci[0][x].theta_wat_nut = valmed;
						break;
					case 31:
						avgstatloci[0][x].theta_taj_nut = valmed;
						break;
					case 32:
						avgstatloci[0][x].theta_fw_nut = valmed;
						break;
					case 33:
						avgstatloci[0][x].theta_fuli_nut = valmed;
						break;
					case 34:
						avgstatloci[0][x].theta_fulin_nut = valmed;
						break;
					case 35:
						avgstatloci[0][x].ndivergence_nut = valmed;
						break;
					case 36:
						avgstatloci[0][x].theta_L_nut = valmed;
						break;
				}
            }
        }
        /**/
        if((vecsort = (double *)realloc(vecsort,sizeof(double))) == 0) {
            puts("Error in alloc. main_sm");
            return 0;
        }
        /**/
        printf(" median calculated succesfully.\n");
        if(file_output) fputs("\n median calculated succesfully.\n",file_output);
    }
    return 1;
}

/*compare two double numbers in a long int list*/
int compare_(const void *i,const void *j) 
{
    if(*(double *)i < *(double *)j) return -1;
    if(*(double *)i > *(double *)j) return  1;
    return 0;
}

