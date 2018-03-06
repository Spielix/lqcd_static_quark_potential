#include "static_quark_potential.h"

int sign(double x);

#ifdef MAIN__RAND_SU2_C__
int main(){
    PAR *par = NULL;
    /* allocate memory for simulation parameters */
     par = malloc(sizeof(PAR));
     if (par == NULL) {
        printf("Error: Failed allocating memory for simulation parameters. Exiting...\n");
        exit(EXIT_FAILURE);
     }

    
    /* initialize simulation parameters with statndard values */
    par->L = 0;
    par->seed = 0;
    par->n_configs = 10;
    par->gen_type = gsl_rng_ranlxs0;
    par->n_su2 = 20;
    par->eps = .24;
    par->n_hits = 10;
    par->n_therm = 500;
    par->n_corr = 50;
    par->beta = 0.;
    par->ran_gen = gsl_rng_alloc(par->gen_type);
    
    gsl_matrix_complex **su2 = malloc(20*sizeof(gsl_matrix_complex*));
    for (int i = 0; i < 20; i++){
		su2[i] = gsl_matrix_complex_calloc(2,2);
	}
	init_su2(par, su2);
	
	gsl_matrix_complex *mat_tmp = gsl_matrix_complex_calloc(2,2);
	gsl_matrix_complex_set(mat_tmp, 0, 0, gsl_complex_rect(4, 0));
	
	printf("%i\n", check_su2(mat_tmp, mat_tmp, 1e-5));
	/*
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			temp = gsl_matrix_complex_get(su2[5], i, j);
			printf("%g + %gi\t", GSL_REAL(temp), GSL_IMAG(temp));
		}
		printf("\n");
	}*/
	return 0;
 }
 #endif

int init_su2(PAR *par, gsl_matrix_complex **su2){
	gsl_matrix_complex **s = malloc(4*sizeof(gsl_matrix_complex*));
	gsl_complex I = gsl_complex_rect(0,1);
	gsl_complex UNIT = gsl_complex_rect(1,0);
	/*init Pauli Matrices*/
	for(int i = 0; i < 4; i++){
		s[i] = gsl_matrix_complex_calloc(2,2);
	}
	gsl_matrix_complex_set(s[0], 0, 0, UNIT);
	gsl_matrix_complex_set(s[0], 1, 1, UNIT);
	gsl_matrix_complex_set(s[1], 0, 1, UNIT);
	gsl_matrix_complex_set(s[1], 1, 0, UNIT);
	gsl_matrix_complex_set(s[2], 0, 1, gsl_complex_mul_real(I,-1));
	gsl_matrix_complex_set(s[2], 1, 0, I);
	gsl_matrix_complex_set(s[3], 0, 0, gsl_complex_mul_real(UNIT, -1));
	gsl_matrix_complex_set(s[3], 1, 1, UNIT);
	
	
	
	double r[4];
	double r_abs;
	
	gsl_matrix_complex **temp = malloc(4*sizeof(gsl_matrix_complex*));
	for(int i = 0; i < 4; i++){
		temp[i] = gsl_matrix_complex_calloc(2,2);
	}
	
	for(int m = 0; m < par->n_su2/2; m++){
		
		for(int i = 0; i < 4; i++){
			
			r[i] = gsl_rng_uniform_pos(par->ran_gen) - 0.5;
		}
		r_abs = sqrt(r[1]*r[1]+r[2]*r[2]+r[3]*r[3]);
		
		for(int i = 0; i < 4; i++){
			gsl_matrix_complex_memcpy(temp[i], s[i]);
		}
		gsl_matrix_complex_scale(temp[0], gsl_complex_rect(sign(r[0])*sqrt(1-par->eps*par->eps), 0));
		for(int i = 1; i < 4; i++){
			gsl_matrix_complex_scale(temp[i], gsl_complex_rect(0, par->eps*r[i]/r_abs));
		}
		for(int i = 0; i < 4; i++){
			gsl_matrix_complex_add(su2[m], temp[i]);
		}
	}
	for(int m = par->n_su2/2; m < par->n_su2; m++){
		psl_matrix_complex_dagger_memcpy(su2[m], su2[m-par->n_su2/2]);
	}
	for(int i = 0; i < 4; i++){
		gsl_matrix_complex_free(s[i]);
		gsl_matrix_complex_free(temp[i]);
	}
	free(s);
	free(temp);
	return 0;
}


int sign(double x) 
{ 
    if (x==0) 
        return 0; 
    else 
        return (x>0) ? 1 : -1; 
}

int check_su2(gsl_matrix_complex *matrix, gsl_matrix_complex *dagger, double epsilon){
	gsl_matrix_complex *unity = gsl_matrix_complex_calloc(2,2);
	gsl_matrix_complex_set(unity, 0, 0, gsl_complex_rect(1,0));
	gsl_matrix_complex_set(unity, 1, 1, gsl_complex_rect(1,0));
	gsl_matrix_complex *check_unitary = gsl_matrix_complex_calloc(2,2);
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), matrix, dagger, gsl_complex_rect(0,0), check_unitary);
	gsl_matrix_complex_sub(unity, check_unitary);
	int ret = 0;
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			if(gsl_complex_abs(gsl_matrix_complex_get(unity, i, j)) > epsilon){
				ret++;
				break;
			}
		}
	}
	double tmp = gsl_complex_abs(gsl_complex_sub(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,0), gsl_matrix_complex_get(matrix,1,1)), gsl_complex_mul(gsl_matrix_complex_get(matrix,1,0), gsl_matrix_complex_get(matrix,0,1))));
	if(tmp > 1+epsilon) ret++;
	return ret;
}
