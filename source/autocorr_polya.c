/*
 * gcc -O3 -o summary summary.c -lm
 */

#include <stdio.h>
#include <float.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <gsl/gsl_statistics_double.h>

#include "static_quark_potential.h"

double *auto_covariance(double **data, double tot_mean, int *n_dataset, int n_replicas) {
    double *cov = NULL;
    int n_tot = 0;
    cov = malloc(n_dataset[0] * sizeof(double));
    for (int replica = 0; replica < n_replicas; replica++) 
        n_tot += n_dataset[replica];
    for (int time = 0; time < n_dataset[0]; time++) {
        cov[time] = 0.;
        for (int replica = 0; replica < n_replicas; replica++) {
            for (int i = 0; i < n_dataset[replica] - time; i++) {
                cov[time] += (data[replica][i] - tot_mean) * (data[replica][i + time] - tot_mean);
            }
        }
        cov[time] /= (double)(n_tot - n_replicas * time);
    }
    return cov;
}

double gamma_error(
    double **data, 
    double tot_mean, 
    int *n_dataset, 
    int n_replicas, 
    double *std_dev, 
    double *t_corr, 
    double *t_corr_err
) {
	double integ_corr_time = 0.5, corr_time, g = 1, err;
	double *cov = auto_covariance(data, tot_mean, n_dataset, n_replicas);
    int window, n_tot = 0;
    for (int replica = 0; replica < n_replicas; replica++) 
        n_tot += n_dataset[replica];
	
    if (std_dev) 
        *std_dev = sqrt(n_tot * cov[0] / (n_tot - 1));

	for(window = 1; (window < n_dataset[0]) && (g > 0); window++){
		integ_corr_time += cov[window] / cov[0];
		corr_time = 1 / log((2 * integ_corr_time + 1) / (2 * integ_corr_time - 1));
		g = exp(-((double)window / corr_time)) - corr_time / sqrt((double)(window * n_tot));
	}
    integ_corr_time *= 1. + (2. * (double)window + 1.) / (double)n_tot;
	if (t_corr) 
        *t_corr = integ_corr_time;
    if (t_corr_err) 
        *t_corr_err = 4. * ((double)window + 0.5 - integ_corr_time) * integ_corr_time * integ_corr_time / (double)n_tot;


	err = sqrt(2. * integ_corr_time * cov[0] / (double)n_tot);
	free(cov);
    return err;
}

int analyze_datasets(int arg_offset, char **argv, int n_replicas) {
    char name[1000], s_temp[1000];
    double tot_mean, std_dev, t_corr, err, t_corr_err, beta;
    double **data = NULL, *means = NULL, *n_dataset_fl = NULL;
    int *n_dataset = NULL, L_t;
    data = malloc(n_replicas * sizeof(double *));
    means = malloc(n_replicas * sizeof(double));
    n_dataset = malloc(n_replicas * sizeof(int));
    n_dataset_fl = malloc(n_replicas * sizeof(double));
    if ((data == NULL) || (means == NULL) || (n_dataset == NULL) || (n_dataset_fl == NULL)) {
        printf("Error: Failed allocating memory for data! Exiting...\n");
        return 1;
    }
    for (int replica = 0; replica < n_replicas; replica++) {
        PAR par;
        FILE *input = fopen(argv[arg_offset + replica], "rb");
        if(!input){
            printf("Could not open file %s. Exiting...\n", argv[arg_offset + replica]);
            return 1;
        }

        fread(&par, sizeof(PAR), 1, input);
        fclose(input);

        n_dataset[replica] = par.n_configs;
        n_dataset_fl[replica] = (double)n_dataset[replica];
        if (replica == 0) {
            beta = par.beta;
            L_t = par.L_t;
        } else {
            if ((par.beta != beta) || (par.L_t != L_t)) {
                printf("Error: simulation Parameters aren't the same! Exiting...");
                return 1;
            }
        }

        strcpy(s_temp, argv[arg_offset + replica]);
        s_temp[strlen(s_temp) - 4] = '\0';
        sprintf(name, "polyakov/%s_polyakov.bin", s_temp);
        input = NULL;
        input = fopen(name, "rb");
        if (!input) {
            printf("Could not open polyakov file %s. Exiting...\n", argv[arg_offset + replica]);
            return 1;
        }

        data[replica] = malloc(n_dataset[replica] * sizeof(double));

        for(int i = 0; i < par.n_configs; i++){
            fread(data[replica] + i, sizeof(double), 1, input);
        }
        fclose(input);

        means[replica] = gsl_stats_mean(data[replica], 1, par.n_configs);
    }
    
    tot_mean = gsl_stats_wmean(n_dataset_fl, 1, means, 1, n_replicas);

    for (int replica = 0; replica < n_replicas; replica++) 
        err = gamma_error(data, tot_mean, n_dataset, n_replicas, &std_dev, &t_corr, &t_corr_err);
//        err = error_auto_weight_simple(data[0], tot_mean, n_dataset[0], &std_dev, &t_corr);
    
    printf("%f\t%d\t%e\t%e\t%e\t%e\n", beta, L_t, tot_mean, err, t_corr, t_corr_err);
    return 0;
}


int main(int argc, char **argv) {
    int n_replicas;
    if(argc < 2){
		printf("Usage: ./autocorr {INPUTFILE_1} {...}\n");
		return 1;
	}
    printf("beta\tL_t\tmean\terror\tcorr_time\n");
    for (int arg = 1; arg < argc;) {
        n_replicas = atoi(argv[arg]);
        if (analyze_datasets(arg + 1, argv, n_replicas)) 
            exit(EXIT_FAILURE);
        arg += n_replicas + 1;
    }
	
    return 0;
}
