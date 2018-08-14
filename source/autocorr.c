#include <stdio.h>
#include <float.h>
#include <math.h>
#include <complex.h>
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
    if (t_corr_err) *t_corr_err = sqrt(4. * ((double)window + 0.5 - integ_corr_time) * integ_corr_time * integ_corr_time / (double)n_tot);


	err = sqrt(2. * integ_corr_time * cov_0 / (double)n_tot);
    return err;
}

int analyze_datasets(int arg_offset, char **argv, int n_replicas) {
    double **tot_mean, 
           **std_dev, 
           **t_corr, 
           **err, 
           **t_corr_err, 
           beta; 
    double ****data = NULL, 
           ***means = NULL, 
           *n_dataset_fl = NULL; 
    int *n_dataset = NULL, L_t, L;

    n_dataset = malloc(n_replicas * sizeof(int));
    n_dataset_fl = malloc(n_replicas * sizeof(double));
    
    for (int replica = 0; replica < n_replicas; replica++) {
        PAR par;
        FILE *input = fopen(argv[arg_offset + replica], "rb");
        if(!input){
            printf("Could not open file %s. Exiting...\n", argv[arg_offset + replica]);
            return 1;
        }

        fread(&par, sizeof(PAR), 1, input);

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

        data = malloc((L_t - 1) * sizeof(double ***));
        
        if (replica == 0) {
            for(int i = 0; i < L_t - 1; i++) {
                data[i] = malloc((L - 1) * sizeof(double **));
                for(int j = 0; j < L - 1; j++) {
                    data[i][j] = malloc(n_replicas * sizeof(double *));
                }
            }
        }
        for(int i = 0; i < L_t - 1; i++) {
            for(int j = 0; j < L - 1; j++) {
                data[i][j][replica] = malloc(n_dataset[replica] * sizeof(double));        
            }
        }
		for (int k = 0; k < n_dataset[replica]; k++) {
			for (int i = 0; i < L_t - 1; i++) {
				for (int j = 0; j < L - 1; j++) { 
                    fread(data[i][j][replica] + k, sizeof(double), 1, input);
				}
			}
        }
        fclose(input);
    }
        
    means = malloc((L_t - 1) * sizeof(double **));
    for (int i = 0; i < L_t - 1; i++) {
        means[i] = malloc((L - 1) * sizeof(double *));
        for (int j = 0; j < L - 1; j++) {
            means[i][j] = malloc(n_replicas * sizeof(double));
        }
    }

	if ((data == NULL) || (means == NULL) || (n_dataset == NULL) || (n_dataset_fl == NULL)) {
        printf("Error: Failed allocating memory for data! Exiting...\n");
        return 1;
    }
    for (int replica = 0; replica < n_replicas; replica++) {

        for (int i = 0; i < L_t - 1; i++) {
            for (int j = 0; j < L - 1; j++) {
                means[i][j][replica] = gsl_stats_mean(data[i][j][replica], 1, n_dataset[replica]);
            }
        }
    }
    tot_mean = malloc((L_t - 1) * sizeof(double *));
    for (int i = 0; i < L_t - 1; i++) {
        tot_mean[i] = malloc((L - 1) * sizeof(double));
        for (int j = 0; j < L - 1; j++) {
            tot_mean[i][j] = gsl_stats_wmean(n_dataset_fl, 1, means[i][j], 1, n_replicas);
        }
    }

    std_dev = malloc((L_t - 1) * sizeof(double *));
    t_corr = malloc((L_t - 1) * sizeof(double *));
    t_corr_err = malloc((L_t - 1) * sizeof(double *));
    err = malloc((L_t - 1) * sizeof(double *));
    for(int i = 0; i < L_t - 1; i++){
        std_dev[i] = malloc((L - 1) * sizeof(double));
        t_corr[i] = malloc((L - 1) * sizeof(double));
        t_corr_err[i] = malloc((L - 1) * sizeof(double));
        err[i] = malloc((L - 1) * sizeof(double));
        for (int j = 0; j < L - 1; j++) {
            err[i][j] = gamma_error(data[i][j], tot_mean[i][j], n_dataset, n_replicas, std_dev[i] + j, t_corr[i] + j, t_corr_err[i] + j);
			printf("%d,%d\t%g\t%g\n", i + 1, j + 1, t_corr[i][j], t_corr_err[i][j]);
        }
    }
    
    free(n_dataset);
    free(n_dataset_fl);
    for (int i = 0; i < L_t - 1; i++) {
        for (int j = 0; j < L - 1; j++) {
            for (int replica = 0; replica < n_replicas; replica++) { 
                free(data[i][j][replica]);
            }
            free(data[i][j]);
            free(means[i][j]);
        }
        free(data[i]);
        free(means[i]);
    }
    free(data);
    free(means);

	FILE *output;
	char filename[1000];
	char mkdir[1000];
    sprintf(filename, "L%d_Lt%d_beta%f", L, L_t, beta);
	snprintf(mkdir, 100, "mkdir %s_pltdata", filename);
	system(mkdir);
	sprintf(mkdir, "%s", filename);
	for (int r = 0; r < L - 1; r++) {
		snprintf(filename, 100, "%s_pltdata/r%d.csv", mkdir, r + 1);
		output = fopen(filename, "w");
		if (output == NULL) {
			printf("Failed opening file %s for writing plot data. Exiting", filename);
		}
		for(int t = 0; t < L_t - 1; t++){
			fprintf(
                output, 
                "%d\t%g\t%g\n", 
                t + 1, 
                tot_mean[t][r], 
                err[t][r]
            );
		}
		fclose(output);
	}
	for (int i = 0; i < L_t - 1; i++) {
		free(tot_mean[i]);
		free(err[i]);
		free(t_corr[i]);
		free(t_corr_err[i]);
		free(std_dev[i]);
	}
	free(tot_mean);
	free(err);
	free(t_corr);
	free(t_corr_err);
	free(std_dev);
	
    return 0;
}

int main(int argc, char **argv) {
    int n_replicas;
    if(argc < 2) {
		printf("Usage: ./autocorr {n_files1} {file_1...} {...}\n");
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
