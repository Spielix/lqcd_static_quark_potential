/*
 * gcc -O3 -o summary summary.c -lm
 */

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <strings.h>

#include "static_quark_potential.h"


double arith_mittel(double *x, unsigned n){
	unsigned i;
	double sum=0;
	for(i = 0; i < n; i++) sum += x[i];
	return sum/n;
}

double stabw(double *x, unsigned n, double mu){
	unsigned i;
	double sum=0, summand;
	for(i = 0; i < n; i++){
		summand = x[i]-mu;
		sum += summand*summand;
	}
	return sqrt(sum/(n-1));
}

unsigned hoch2(unsigned n){
	unsigned potenz=1, x=2;
	while(n){
		if(n%2)	potenz *= x;
		x *= x;
		n /= 2;
	}
	return potenz;
}

	/*very fast transform: no time waisted with memory allocation.
	needs input, output and dummy arrays of lengths 2^r*/
double complex *vfft(double *f, double complex *h, double complex *g, unsigned n, int hin){
	unsigned r = (unsigned)(log2(n)*(1+DBL_EPSILON)), m=n/2, l=1;
	unsigned i, k, j;
	unsigned a, b;
	double complex w, wSchritt;
	const double wurzelN = 1/sqrt(n);
	const double complex phase=-hin*M_PI*I;

	for(j = 0; j < n; j++) g[j] = f[j]*wurzelN;
	
	for(i = 0; i < r; i++){
		w = cexp(phase/m);
		for(k = 0; k < l; k++){
			wSchritt = 1;
			for(j = 0; j < m; j++){
				a = 2*k*m + j;
				b = a+m;
				g[a] += g[b];
				g[b] = wSchritt*(g[a] - 2*g[b]);
				wSchritt *= w;
			}
		}
		m /= 2;
		l *= 2; 
	}

		/*sortieren*/
	j = 0;
	for(i = 0; i < n-1; i++){
		h[i] = g[j];
		for(m = n/2; m <= j; m /= 2){
			j -= m;
		}
		j += m;
	}
	h[n-1] = g[n-1];

	return h;
}

	/*very fast transform: no time waisted with memory allocation.
	needs input, output and dummy arrays of lengths 2^r
	input and output may be identical, dummy array may not*/
double complex *cvfft(double complex *f, double complex *h, double complex *g, unsigned n, int hin){
	unsigned r = (unsigned)(log2(n)*(1+DBL_EPSILON)), m=n/2, l=1;
	unsigned i, k, j;
	unsigned a, b;
	double complex w, wSchritt;
	const double wurzelN = 1/sqrt(n);
	const double complex phase=-hin*M_PI*I;

	for(j = 0; j < n; j++) g[j] = f[j]*wurzelN;
	
	for(i = 0; i < r; i++){
		w = cexp(phase/m);
		for(k = 0; k < l; k++){
			wSchritt = 1;
			for(j = 0; j < m; j++){
				a = 2*k*m + j;
				b = a+m;
				g[a] += g[b];
				g[b] = wSchritt*(g[a] - 2*g[b]);
				wSchritt *= w;
			}
		}
		m /= 2;
		l *= 2; 
	}

		/*sortieren*/
	j = 0;
	for(i = 0; i < n-1; i++){
		h[i] = g[j];
		for(m = n/2; m <= j; m /= 2){
			j -= m;
		}
		j += m;
	}
	h[n-1] = g[n-1];

	return h;
}

	/*input- and output-arrays of lengths n
	and four dummy arrays of lengths m=2^log(2n-1)*/
double complex *chirp_z_vfft(double *f, double complex *h, double complex *a, double complex *b, double complex *bft, double complex *dummy, unsigned n, unsigned m, int hin){
	unsigned k, test=hoch2((unsigned)(log2(n)*(1+DBL_EPSILON)));
	const double norm = sqrt((double)m/n);
	const double complex phase = hin*M_PI*I/n;

	if(n == test) return vfft(f, h, dummy, n, hin);

	for(k = 0; k < n; k++){
		b[k] = cexp(phase*k*k);
		a[k] = f[k]*conj(b[k]);
	}
	for(; k <= m-n; k++){
		a[k] = 0;
		b[k] = 0;
	}
	for(; k < m; k++){
		a[k] = 0;
		b[k] = b[m-k];
	}

	cvfft(a, a, dummy, m, 1);
	cvfft(b, bft, dummy, m, 1);
	for(k = 0; k < m; k++) a[k] *= bft[k];
	cvfft(a, a, dummy, m, -1);
	for(k = 0; k < n; k++) h[k] = a[k]*conj(b[k])*norm;

	return h;
}

double *fast_auto_cov(double *x, double mu, unsigned n){
	double complex *cov;
	double complex *ab, *bft, *dummy;
	double *abs_sq;
	double a, b;
	const double wurzelN=1/sqrt(n);
	unsigned i, m;

	m = hoch2((unsigned)(log2(2*n-1)*(1-DBL_EPSILON)+1));
	abs_sq = malloc(n*sizeof(double));
	cov = malloc(m*sizeof(double complex));
	ab = malloc(m*sizeof(double complex));
	bft = malloc(m*sizeof(double complex));
	dummy = malloc(m*sizeof(double complex));

	for(i = 0; i < n; i++) abs_sq[i] = x[i]-mu;
	chirp_z_vfft(abs_sq, cov, cov, ab, bft, dummy, n, m, -1);
	/*cimagprint_array(cov, n);*/
	for(i = 0; i < n; i++){
		a = creal(cov[i]);
		b = cimag(cov[i]);
		abs_sq[i] = a*a+b*b;
	}
	chirp_z_vfft(abs_sq, cov, cov, ab, bft, dummy, n, m, 1);
	for(i = 0; i < n; i++) abs_sq[i] = creal(cov[i])*wurzelN;
	free(cov);
	free(ab);
	free(bft);
	free(dummy);

	return abs_sq;
}

double error_auto_weight_simple(double *x, double mu, unsigned n, double *std_dev, double *t_corr){
	unsigned i;
	double cov_int;
	double *cov=fast_auto_cov(x, mu, n);

	if(std_dev) *std_dev = sqrt(cov[0]*n/(n-1));

	cov_int = cov[0]*0.5;
	for(i = 1; i < n && cov[i] > 0; i++){
		cov_int += cov[i];
	}
	if(t_corr) *t_corr = cov_int/cov[0];
	free(cov);

	return sqrt(2*cov_int/(n-1));
}

double error_auto_weight(double *x, double mu, unsigned n, double *std_dev, double *t_corr){
	unsigned i;
	double t_int=0.5, t, g=1;
	double cov_0_inv;
	double *cov=fast_auto_cov(x, mu, n);

	cov_0_inv = 1/cov[0];
	if(std_dev) *std_dev = sqrt(n/(cov_0_inv*(n-1)));

	for(i = 1; i < n && g > 0; i++){
		t_int += cov[i]*cov_0_inv;
		t = 1/log((2*t_int+1)/(2*t_int-1));
		g = exp(-(i/t))-t/sqrt(i*n);
	}
	if(t_corr) *t_corr = t_int;
	free(cov);

	return sqrt(2*t_int/(cov_0_inv*(n-1)));
}

int main(int argc, char **argv){
	if(argc != 2){
		printf("Usage: ./autocorr {INPUTFILE}");
		return 1;
	}
	FILE *input = fopen(argv[1],"rb");
	if(!input){
		printf("Could not open file\n");
		return 1;
	}
	PAR par;
	fread(&par, sizeof(PAR), 1, input);
	double ***data = malloc((par.L)*sizeof(double**));
	for(int i = 0; i < par.L-1; i++) data[i] = malloc(par.L*sizeof(double*));
	for(int i = 0; i < par.L-1; i++){
		for(int j = i; j < par.L-1; j++){
			data[i][j] = malloc(par.n_configs*sizeof(double));
		}
	}
	for(int i = 0; i < par.n_configs; i++){
		for(int j = 0; j < par.L-1; j++){
			for(int k = j; k < par.L-1; k++){
				fread(data[j][k]+i, sizeof(double), 1, input);
			}
		}
	}
	fclose(input);
	double **mu = malloc((par.L)*sizeof(double*));
	double **std_dev = malloc((par.L)*sizeof(double*));
	double **t_corr = malloc((par.L)*sizeof(double*));
	double **err = malloc((par.L)*sizeof(double*));
	for(int i = 0; i < par.L; i++){
		mu[i] = malloc(par.L*sizeof(double));
		std_dev[i] = malloc(par.L*sizeof(double));
		t_corr[i] = malloc(par.L*sizeof(double));
		err[i] = malloc(par.L*sizeof(double));
	}
	for(int i = 0; i < par.L-1; i++) {
		for(int j = i; j < par.L-1; j++){
			mu[i][j] = arith_mittel(data[i][j], par.n_configs);
		}
	}
	for(int i = 0; i < par.L-1; i++) {
		for(int j = i; j < par.L-1; j++) {
			err[i][j] = error_auto_weight(data[i][j], mu[i][j], par.n_configs, std_dev[i]+j, t_corr[i]+j);
		}
	}
	double *wt = malloc(par.L*sizeof(double));
	double *wt_err = malloc(par.L*sizeof(double));
	FILE *output;
	char filename[100];
	char mkdir[100];
	for(int r = 0; r < par.L-1; r++){
		for(int t = r-1; t >= 0; t--){
			wt[t] = mu[t][r];
			wt_err[t] = err[t][r];
		}
		for(int t = r; t < par.L; t++){
			wt[t] = mu[r][t];
			wt_err[t] = err[r][t];
		}
		snprintf(mkdir, 100, "mkdir %s_pltdata", argv[1]);
		system(mkdir);
		snprintf(filename, 100, "%s_pltdata/%d.txt", argv[1], r);
		output = fopen(filename, "w");
		for(int t = 0; t < par.L-2; t++){
			fprintf(output, "%d\t%g\t%g\n", t, log(wt[t]/wt[t+1]), sqrt((wt_err[t]/wt[t])*(wt_err[t]/wt[t])+(wt_err[t+1]/wt[t+1])*(wt_err[t+1]/wt[t+1])));
		}
		fclose(output);
	}
		
	return 0;
}