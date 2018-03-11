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
    par->eps = 0;
    par->n_hits = 10;
    par->n_therm = 500;
    par->n_corr = 50;
    par->beta = 0.;
    par->tadpole = 0.;
    par->ran_gen = gsl_rng_alloc(par->gen_type);
    
    gsl_matrix_complex **su2 = malloc(20*sizeof(gsl_matrix_complex*));
    for (int i = 0; i < 20; i++){
		su2[i] = gsl_matrix_complex_calloc(2,2);
	}
	init_su2(par, su2);
	
	gsl_complex temp;
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			temp = gsl_matrix_complex_get(su2[5], i, j);
			printf("%g + %gi\t", GSL_REAL(temp), GSL_IMAG(temp));
		}
		printf("\n");
	}
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			temp = gsl_matrix_complex_get(su2[2], i, j);
			printf("%g + %gi\t", GSL_REAL(temp), GSL_IMAG(temp));
		}
		printf("\n");
	}
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			temp = gsl_matrix_complex_get(su2[8], i, j);
			printf("%g + %gi\t", GSL_REAL(temp), GSL_IMAG(temp));
		}
		printf("\n");
	}
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
		gsl_matrix_complex_scale(temp[0], gsl_complex_rect(/*sign(r[0])* */sqrt(1-par->eps*par->eps), 0));
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
	gsl_matrix_complex_free(unity);
	gsl_matrix_complex_free(check_unitary);
	return ret;
}

/* calculate the position of a specific link in the array */
static inline int ind(int i, int j, int k, int l, int dir, int le) {
    return (((i * le + j) * le + k) * le + l) * 8 + dir;
}

/* calculate the position of a specific link in the array, while applying periodic boundary conditions to all
 * coordinates */
static inline int periodic_ind(int i, int j, int k, int l, int dir, int le_0, int le) {
    return (((((i + le_0) % le_0) * le + ((j + le) % le)) * le + ((k + le) % le)) * le + ((l + le) % le)) * 8 + dir;
}

static inline int msl_mat_mul(gsl_matrix_complex *A, gsl_matrix_complex *B, gsl_matrix_complex *C){
	return gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), A, B, gsl_complex_rect(0,0), C);
}

/* calculate the position of a specific link in the gauge array */
static inline int gauge_ind(int i, int j, int k, int l, int dagger, int le) {
    return (((le * i + j) * le + k) * le + l) * 2 + dagger;
}

/* calculate the position of a specific link in the gauge array, while applying periodic boundary conditions to all
 * coordinates */
static inline int gauge_periodic_ind(int i, int j, int k, int l, int dagger, int le_0, int le) {
    return (((le * ((i  + le_0) % le_0) + ((j + le) % le)) * le + ((k + le) % le)) * le + ((l + le) % le)) * 2 + dagger;
}

void init_gauge(PAR *par, gsl_matrix_complex **gauge, gsl_matrix_complex **su2) {
    /* initialize all independent links with random SU(2) matrices */
    for (int i = 0; i < par->L_t; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    gsl_matrix_complex_memcpy(
                        gauge[gauge_ind(i, j, k, l, 0, par->L)], 
                        su2[(int)(gsl_rng_uniform(par->ran_gen) * par->n_su2)]
                    );
                }
            }
        }
    }

    /* generate the daggered matrices */
    for (int i = 0; i < par->L_t; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    psl_matrix_complex_dagger_memcpy(
                        gauge[gauge_ind(i, j, k, l, 1, par->L)], 
                        gauge[gauge_ind(i, j, k, l, 0, par->L)]
                    );
                }
            }
        }
    }
}

double gauge_inv(PAR *par, gsl_matrix_complex **lattice){
	gsl_matrix_complex **gauge = malloc((par->L*par->L*par->L*par->L*8)*sizeof(gsl_matrix_complex*)), 
            **gauged_lattice = malloc((par->L*par->L*par->L*par->L*8)*sizeof(gsl_matrix_complex*));
	for(int i = 0; i < par->L*par->L*par->L*par->L*8; i++) 
        gauge[i] = gsl_matrix_complex_calloc(2,2);
	for(int i = 0; i < par->L*par->L*par->L*par->L*8; i++) {
		gauged_lattice[i] = gsl_matrix_complex_calloc(2,2);
		gsl_matrix_complex_memcpy(gauged_lattice[i], lattice[i]);
	}
	gsl_matrix_complex **temp = malloc(par->n_su2*sizeof(gsl_matrix_complex*));
	for(int i = 0; i < par->n_su2; i++){
		temp[i] = gsl_matrix_complex_calloc(2,2);
	}

	init_su2(par, temp);
	init_gauge(par, gauge, temp);
	for(int a = 0; a < par->L_t; a++){
		for(int b = 0; b < par->L; b++){
			for(int c = 0; c < par->L; c++){
				for(int d = 0; d < par->L; d++){
					for(int dir = 0; dir < 4; dir++){
						msl_mat_mul(
							gauge[gauge_ind(a, b, c, d, 0, par->L)],
							gauged_lattice[ind(a, b, c, d, dir, par->L)],
							par->m_workspace);
						msl_mat_mul(
							par->m_workspace,
							gauge[gauge_periodic_ind(
								a + (dir == 0),
								b + (dir == 1),
								c + (dir == 2),
								d + (dir == 3), 
								1, 
                                par->L_t, 
                                par->L
                            )],
							gauged_lattice[ind(a, b, c, d, dir, par->L)]);
					}
				}
			}
		}
	}
	for (int i = 0; i < par->L_t; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    for (int dir_dagger = 4; dir_dagger < 8; dir_dagger++) {
                        psl_matrix_complex_dagger_memcpy(
                            gauged_lattice[ind(i, j, k, l, dir_dagger, par->L)], 
                            gauged_lattice[periodic_ind(
                                i - (dir_dagger == 4), 
                                j - (dir_dagger == 5), 
                                k - (dir_dagger == 6), 
                                l - (dir_dagger == 7), 
                                dir_dagger - 4, 
                                par->L_t, 
                                par->L
                            )]
                        );
                    }
                }
            }
        }
	}
	double action[2] = {0};
	measure_action_r(par, lattice, action);
	measure_action_r(par, gauged_lattice, action+1);
    for(int i = 0; i < par->n_su2; i++) 
		gsl_matrix_complex_free(temp[i]);
    for (int i = 0; i < par->L * par->L * par->L * par->L * 2; i++) 
        gsl_matrix_complex_free(gauge[i]);
    for (int i = 0; i < par->L * par->L * par->L * par->L * 8; i++) 
        gsl_matrix_complex_free(gauged_lattice[i]);
	free(temp);
	free(gauge);
	free(gauged_lattice);
	return fabs(action[0]-action[1]);
}

/* applicates a random local gauge transformation onto the lattice. all operators should be invariant under this
 * transformation */
int gauge_transform_lattice(PAR *par, gsl_matrix_complex **lattice){
	gsl_matrix_complex **gauge = malloc((par->L_t*par->L*par->L*par->L*2)*sizeof(gsl_matrix_complex*)); 
    if (gauge == NULL) {
        printf("Error: Failed allocating memory for lattice of gauge matrices. Exiting...\n");
        return 1;
    }
	gsl_matrix_complex **su2_gauge = malloc(par->n_su2*sizeof(gsl_matrix_complex*));
    if (su2_gauge == NULL) {
        printf("Error: Failed allocating memory for array of SU(2) matrices for the gauge transformation. Exiting...\n");
        free(gauge);
        return 1;
    }
	for(int i = 0; i < par->L_t*par->L*par->L*par->L*2; i++) {
        gauge[i] = gsl_matrix_complex_calloc(2,2);
        if (gauge[i] == NULL) {
            printf("Error: Failed allocating memory for gauge matrices. Exiting...\n");
            for (int j = 0; j < i; j++) 
                gsl_matrix_complex_free(gauge[j]);
            free(gauge);
            free(su2_gauge);
            return 1;
        }
    }
	for(int i = 0; i < par->n_su2; i++){
		su2_gauge[i] = gsl_matrix_complex_calloc(2,2);
        if (su2_gauge[i] == NULL) {
            printf("Error: Failed allocating memory for SU(2) matrices for the gauge transformation. Exiting...\n");
            for (int j = 0; j < i; j++) 
                gsl_matrix_complex_free(su2_gauge[j]);
            for (int j = 0; j < par->L_t * par->L * par->L * par->L * 2; j++) 
                gsl_matrix_complex_free(gauge[j]);
            free(su2_gauge);
            free(gauge);
            return 1;
        }
	}

	init_su2(par, su2_gauge);
	init_gauge(par, gauge, su2_gauge);
	for(int a = 0; a < par->L_t; a++){
		for(int b = 0; b < par->L; b++){
			for(int c = 0; c < par->L; c++){
				for(int d = 0; d < par->L; d++){
					for(int dir = 0; dir < 4; dir++){
						msl_mat_mul(
							gauge[gauge_ind(a, b, c, d, 0, par->L)],
							lattice[ind(a, b, c, d, dir, par->L)],
							par->m_workspace);
						msl_mat_mul(
							par->m_workspace,
							gauge[gauge_periodic_ind(
								a + (dir == 0),
								b + (dir == 1),
								c + (dir == 2),
								d + (dir == 3), 
								1, 
                                par->L_t, 
                                par->L
                            )],
							lattice[ind(a, b, c, d, dir, par->L)]);
					}
				}
			}
		}
	}
	for (int i = 0; i < par->L_t; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    for (int dir_dagger = 4; dir_dagger < 8; dir_dagger++) {
                        psl_matrix_complex_dagger_memcpy(
                            lattice[ind(i, j, k, l, dir_dagger, par->L)], 
                            lattice[periodic_ind(
                                i - (dir_dagger == 4), 
                                j - (dir_dagger == 5), 
                                k - (dir_dagger == 6), 
                                l - (dir_dagger == 7), 
                                dir_dagger - 4, 
                                par->L_t, 
                                par->L
                            )]
                        );
                    }
                }
            }
        }
	}
	
    for(int i = 0; i < par->n_su2; i++) 
		gsl_matrix_complex_free(su2_gauge[i]);
    for (int i = 0; i < par->L_t * par->L * par->L * par->L * 2; i++) 
        gsl_matrix_complex_free(gauge[i]);
	free(su2_gauge);
	free(gauge);
	
    return 0;
}



