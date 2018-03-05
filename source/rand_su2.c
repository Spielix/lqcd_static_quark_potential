#include "static_quark_potential.h"

int sign(double x);


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
	
	for(int m = 0; m < par->n_su3/2; m++){
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
