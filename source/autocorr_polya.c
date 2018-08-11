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

double auto_covariance(double **data, double tot_mean, int *n_dataset, int n_replicas, int time, int n_tot) {
    double cov = 0.;
    for (int replica = 0; replica < n_replicas; replica++) {
        for (int i = 0; i < n_dataset[replica] - time; i++) {
            cov += (data[replica][i] - tot_mean) * (data[replica][i + time] - tot_mean);
        }
    }
    cov /= (double)(n_tot - n_replicas * time);
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
    int window, n_tot = 0;
    for (int replica = 0; replica < n_replicas; replica++) n_tot += n_dataset[replica];
	double cov_0 = auto_covariance(data, tot_mean, n_dataset, n_replicas, 0, n_tot);
    if (std_dev) *std_dev = sqrt(n_tot * cov_0 / (n_tot - 1));

	for(window = 1; (window < n_dataset[0]) && (g > 0); window++){
        double cov_w = auto_covariance(data, tot_mean, n_dataset, n_replicas, window, n_tot);
		integ_corr_time += cov_w / cov_0;
		corr_time = 1 / log((2 * integ_corr_time + 1) / (2 * integ_corr_time - 1));
		g = exp(-((double)window / corr_time)) - corr_time / sqrt((double)(window * n_tot));
	}
    integ_corr_time *= 1. + (2. * (double)window + 1.) / (double)n_tot;
	if (t_corr) *t_corr = integ_corr_time;
    if (t_corr_err) *t_corr_err = 4. * ((double)window + 0.5 - integ_corr_time) * integ_corr_time * integ_corr_time / (double)n_tot;


	err = sqrt(2. * integ_corr_time * cov_0 / (double)n_tot);
    return err;
}

int analyze_datasets(int arg_offset, char **argv, int n_replicas) {
    char name[1000], s_temp[1000];
    double tot_mean, 
           tot_mean_abs, 
           tot_mean2, 
           tot_mean4, 
           std_dev, 
           std_dev_abs, 
           std_dev2, 
           std_dev4, 
           t_corr, 
           t_corr_abs, 
           t_corr2, 
           t_corr4, 
           err, 
           err_abs, 
           err2, 
           err4, 
           t_corr_err, 
           t_corr_err_abs, 
           t_corr_err2, 
           t_corr_err4, 
           beta, 
           cov_2_4_tot;
    double **data = NULL, 
           **data_abs = NULL, 
           **data2 = NULL, 
           **data4 = NULL, 
           *means = NULL, 
           *means_abs = NULL, 
           *means2 = NULL, 
           *means4 = NULL, 
           *n_dataset_fl = NULL, 
           *cov_2_4 = NULL;
    int *n_dataset = NULL, L_t, L;
    data = malloc(n_replicas * sizeof(double *));
    data_abs = malloc(n_replicas * sizeof(double *));
    data2 = malloc(n_replicas * sizeof(double *));
    data4 = malloc(n_replicas * sizeof(double *));
    means = malloc(n_replicas * sizeof(double));
    means_abs = malloc(n_replicas * sizeof(double));
    means2 = malloc(n_replicas * sizeof(double));
    means4 = malloc(n_replicas * sizeof(double));
    n_dataset = malloc(n_replicas * sizeof(int));
    n_dataset_fl = malloc(n_replicas * sizeof(double));
    cov_2_4 = malloc(n_replicas * sizeof(double));
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
            L = par.L;
        } else {
            if ((par.beta != beta) || (par.L_t != L_t) || (par.L != L)) {
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
        data_abs[replica] = malloc(n_dataset[replica] * sizeof(double));
        data2[replica] = malloc(n_dataset[replica] * sizeof(double));
        data4[replica] = malloc(n_dataset[replica] * sizeof(double));

        for(int i = 0; i < par.n_configs; i++){
            fread(data[replica] + i, sizeof(double), 1, input);
        }
        for (int i = 0; i < par.n_configs; i++) {
            data_abs[replica][i] = fabs(data[replica][i]);
            data2[replica][i] = data[replica][i] * data[replica][i];
            data4[replica][i] = data2[replica][i] * data2[replica][i];
        }
        fclose(input);

        means[replica] = gsl_stats_mean(data[replica], 1, par.n_configs);
        means_abs[replica] = gsl_stats_mean(data_abs[replica], 1, n_dataset[replica]);
        means2[replica] = gsl_stats_mean(data2[replica], 1, n_dataset[replica]);
        means4[replica] = gsl_stats_mean(data4[replica], 1, n_dataset[replica]);
        cov_2_4[replica] = gsl_stats_covariance_m(data2[replica], 1, data4[replica], 1, n_dataset[replica], means2[replica], means4[replica]);
    }
    
    tot_mean = gsl_stats_wmean(n_dataset_fl, 1, means, 1, n_replicas);
    tot_mean_abs = gsl_stats_wmean(n_dataset_fl, 1, means_abs, 1, n_replicas);
    tot_mean2 = gsl_stats_wmean(n_dataset_fl, 1, means2, 1, n_replicas);
    tot_mean4 = gsl_stats_wmean(n_dataset_fl, 1, means4, 1, n_replicas);
    cov_2_4_tot = gsl_stats_wmean(n_dataset_fl, 1, cov_2_4, 1, n_replicas);

    for (int replica = 0; replica < n_replicas; replica++) {
        err = gamma_error(data, tot_mean, n_dataset, n_replicas, &std_dev, &t_corr, &t_corr_err);
        err_abs = gamma_error(data_abs, tot_mean_abs, n_dataset, n_replicas, &std_dev_abs, &t_corr_abs, &t_corr_err_abs);
        err2 = gamma_error(data2, tot_mean2, n_dataset, n_replicas, &std_dev2, &t_corr2, &t_corr_err2);
        err4 = gamma_error(data4, tot_mean4, n_dataset, n_replicas, &std_dev4, &t_corr4, &t_corr_err4);
    }
    free(n_dataset);
    free(n_dataset_fl);
    free(cov_2_4);
    for (int replica = 0; replica < n_replicas; replica++) { 
        free(data[replica]);
        free(data_abs[replica]);
        free(data2[replica]);
        free(data4[replica]);
    }
    free(data);
    free(data_abs);
    free(data2);
    free(data4);
    free(means);
    free(means_abs);
    free(means2);
    free(means4);
    printf("%g %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
            beta, 
            L_t, 
            L, 
            tot_mean, 
            err, 
            t_corr, 
            std_dev,
            tot_mean_abs, 
            err_abs, 
            t_corr_abs, 
            std_dev_abs, 
            tot_mean2, 
            err2, 
            t_corr2, 
            std_dev2, 
            tot_mean4, 
            err4, 
            t_corr4, 
            std_dev4, 
            cov_2_4_tot,
			t_corr_err, 
			t_corr_err_abs, 
			t_corr_err2, 
			t_corr_err4
    );
    return 0;
}


int main(int argc, char **argv) {
    int n_replicas;
    if(argc < 2){
		printf("Usage: ./autocorr {INPUTFILE_1} {...}\n");
		return 1;
	}
    for (int arg = 1; arg < argc;) {
        n_replicas = atoi(argv[arg]);
        if (analyze_datasets(arg + 1, argv, n_replicas)) 
            exit(EXIT_FAILURE);
        arg += n_replicas + 1;
    }
	
    return 0;
}
