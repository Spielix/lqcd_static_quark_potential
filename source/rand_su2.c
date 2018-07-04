#include "static_quark_potential.h"

int sign(double x);

/* #ifdef MAIN__RAND_SU2_C__
int main(){
    PAR *par = NULL;
   */ /* allocate memory for simulation parameters */
  /*   par = malloc(sizeof(PAR));
     if (par == NULL) {
        printf("Error: Failed allocating memory for simulation parameters. Exiting...\n");
        exit(EXIT_FAILURE);
     }
*/
    
    /* initialize simulation parameters with standard values */
  /*  par->L = 0;
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
    
    double *su2 = calloc(20 * 8, sizeof(double));
	init_su2(par, su2);
	
	gsl_complex temp;
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			temp = gsl_matrix_complex_get(gsl_matrix_complex_view_array(su2 + 8 * 5, 2, 2), i, j);
			printf("%g + %gi\t", GSL_REAL(temp), GSL_IMAG(temp));
		}
		printf("\n");
	}
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			temp = gsl_matrix_complex_get(gsl_matrix_complex_view_array(su2 + 8 * 2, 2, 2), i, j);
			printf("%g + %gi\t", GSL_REAL(temp), GSL_IMAG(temp));
		}
		printf("\n");
	}
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			temp = gsl_matrix_complex_get(gsl_matrix_complex_view_array(su2 + 8 * 8, 2, 2), i, j);
			printf("%g + %gi\t", GSL_REAL(temp), GSL_IMAG(temp));
		}
		printf("\n");
	}
	return 0;
 }
 #endif */

void init_su2(PAR *par, double *su2){
	double pauli[32] = {0};
	gsl_complex 
        z_i = gsl_complex_rect(0., 1.), 
	    z_one = gsl_complex_rect(1., 0.);
	
    /* initialize Pauli matrices */
    gsl_matrix_complex_view m_pauli_view = gsl_matrix_complex_view_array(pauli, 2, 2);
	gsl_matrix_complex_set(&m_pauli_view.matrix, 0, 0, z_one);
	gsl_matrix_complex_set(&m_pauli_view.matrix, 1, 1, z_one);

    m_pauli_view = gsl_matrix_complex_view_array(pauli + 8, 2, 2);
	gsl_matrix_complex_set(&m_pauli_view.matrix, 0, 1, z_one);
	gsl_matrix_complex_set(&m_pauli_view.matrix, 1, 0, z_one);

    m_pauli_view = gsl_matrix_complex_view_array(pauli + 16, 2, 2);
	gsl_matrix_complex_set(&m_pauli_view.matrix, 0, 1, gsl_complex_mul_real(z_i, -1.));
	gsl_matrix_complex_set(&m_pauli_view.matrix, 1, 0, z_i);

    m_pauli_view = gsl_matrix_complex_view_array(pauli + 24, 2, 2);
	gsl_matrix_complex_set(&m_pauli_view.matrix, 0, 0, gsl_complex_mul_real(z_one, -1.));
	gsl_matrix_complex_set(&m_pauli_view.matrix, 1, 1, z_one);
	
	double r[4] = {0};
	double r_abs;
	
	double m_temp_arr[32] = {0};
	gsl_matrix_complex_view m_temp_view, m_su2_view, m_su2_dagg_view;

	for (int m = 0; m < par->n_su2 / 2; m++) {
		
		for (int i = 0; i < 4; i++) 
			r[i] = gsl_rng_uniform_pos(par->ran_gen) - 0.5;
		
        r_abs = sqrt(r[1] * r[1] + r[2] * r[2] + r[3] * r[3]);
		
		for (int i = 0; i < 4; i++) {
            m_temp_view = gsl_matrix_complex_view_array(m_temp_arr + 8 * i, 2, 2);
            m_pauli_view = gsl_matrix_complex_view_array(pauli + 8 * i, 2, 2);
			gsl_matrix_complex_memcpy(
                &m_temp_view.matrix, 
                &m_pauli_view.matrix 
            );
		}
        m_temp_view = gsl_matrix_complex_view_array(m_temp_arr, 2, 2);
		gsl_matrix_complex_scale(
            &m_temp_view.matrix, 
            gsl_complex_rect( /* sign(r[0]) * */ sqrt(1. - par->eps * par->eps), 0.)
        );
		for (int i = 1; i < 4; i++) {
            m_temp_view = gsl_matrix_complex_view_array(m_temp_arr + 8 * i, 2, 2);
			gsl_matrix_complex_scale(
                &m_temp_view.matrix, 
                gsl_complex_rect(0., par->eps * r[i] / r_abs)
            );
		}
        m_su2_view = gsl_matrix_complex_view_array(su2 + 8 * m, 2, 2);
		for (int i = 0; i < 4; i++) {
            m_temp_view = gsl_matrix_complex_view_array(m_temp_arr + 8 * i, 2, 2);
			gsl_matrix_complex_add(
                &m_su2_view.matrix, 
                &m_temp_view.matrix
            );
		}
	}
	for (int m = par->n_su2 / 2; m < par->n_su2; m++) {
        m_su2_dagg_view = gsl_matrix_complex_view_array(su2 + 8 * m, 2, 2);
        m_su2_view = gsl_matrix_complex_view_array(su2 + 8 * m - par->n_su2 / 2, 2, 2);
		psl_matrix_complex_dagger_memcpy(
            &m_su2_dagg_view.matrix, 
            &m_su2_view.matrix
        );
	}
}


int sign(double x) 
{ 
    if (x==0) 
        return 0; 
    else 
        return (x>0) ? 1 : -1; 
}

int check_su2(gsl_matrix_complex *matrix, gsl_matrix_complex *dagger, double epsilon){
	double 
        m_unity[8] = {0},
        m_check_unitary[8] = {0};
    gsl_matrix_complex_view 
        m_unity_view = gsl_matrix_complex_view_array(m_unity, 2, 2),
        m_check_unitary_view = gsl_matrix_complex_view_array(m_check_unitary, 2, 2);

	gsl_matrix_complex_set_identity(&m_unity_view.matrix);

	gsl_blas_zgemm(
        CblasNoTrans, 
        CblasNoTrans, 
        gsl_complex_rect(1, 0), 
        matrix, 
        dagger, 
        gsl_complex_rect(0, 0), 
        &m_check_unitary_view.matrix
    );
	gsl_matrix_complex_sub(&m_unity_view.matrix, &m_check_unitary_view.matrix);
	int ret = 0;
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			if (gsl_complex_abs(gsl_matrix_complex_get(&m_unity_view.matrix, i, j)) > epsilon){
				ret++;
				break;
			}
		}
	}
	double tmp = gsl_complex_abs(
        gsl_complex_sub(
            gsl_complex_mul(
                gsl_matrix_complex_get(
                    matrix, 
                    0, 
                    0
                ), 
                gsl_matrix_complex_get(
                    matrix, 
                    1, 
                    1
                )
            ), 
            gsl_complex_mul(
                gsl_matrix_complex_get(
                    matrix, 
                    1, 
                    0
                ), 
                gsl_matrix_complex_get(
                    matrix, 
                    0, 
                    1
                )
            )
        )
    );
	if(tmp > 1+epsilon) ret++;
	return ret;
}

/* calculate the position of a specific link in the array */
static inline int ind(int i, int j, int k, int l, int dir, int le) {
    return (((i * le + j) * le + k) * le + l) * 4 + dir;
}

/* calculate the position of a specific link in the array, while applying periodic boundary conditions to all
 * coordinates */
static inline int periodic_ind(int i, int j, int k, int l, int dir, int le_0, int le) {
    return (((((i + le_0) % le_0) * 
        le + ((j + le) % le)) * 
        le + ((k + le) % le)) * 
        le + ((l + le) % le)) * 
        4 + dir;
}

static inline int msl_mat_mul(int dag_1, int dag_2, gsl_matrix_complex *A, gsl_matrix_complex *B, gsl_matrix_complex *C){
	return gsl_blas_zgemm(
        dag_1 ? CblasConjTrans : CblasNoTrans, 
        dag_2 ? CblasConjTrans : CblasNoTrans, 
        gsl_complex_rect(1, 0), 
        A, 
        B, 
        gsl_complex_rect(0, 0), 
        C
    );
}

/* calculate the position of a specific link in the gauge array */
static inline int gauge_ind(int i, int j, int k, int l, int le) {
    return ((le * i + j) * le + k) * le + l;
}

/* calculate the position of a specific link in the gauge array, 
 * while applying periodic boundary conditions to all coordinates */
static inline int gauge_periodic_ind(int i, int j, int k, int l, int le_0, int le) {
    return ((le * ((i  + le_0) % le_0) + 
        ((j + le) % le)) * le + 
        ((k + le) % le)) * le + 
        ((l + le) % le);
}

void init_gauge(PAR *par, double *gauge, double *su2) {
    /* initialize all independent links with random SU(2) matrices */
    gsl_matrix_complex_view m_gauge_view, m_su2_view;

    for (int i = 0; i < par->L_t; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    
                    m_gauge_view = gsl_matrix_complex_view_array(
                        gauge + 8 * gauge_ind(i, j, k, l, par->L), 
                        2, 
                        2
                    );
                    m_su2_view = gsl_matrix_complex_view_array(
                        su2 + 8 * (int)(gsl_rng_uniform(par->ran_gen) * par->n_su2), 
                        2, 
                        2
                    );
                    gsl_matrix_complex_memcpy(
                        &m_gauge_view.matrix, 
                        &m_su2_view.matrix 
                    );
                }
            }
        }
    }
}

/* double gauge_inv(PAR *par, gsl_complex *lattice){
	gsl_complex 
        *gauge = calloc(par->L_t * par->L * par->L * par->L * 8 * 8, sizeof(double)), 
        *gauged_lattice = calloc(par->L_t * par->L * par->L * par->L * 8 * 8, sizeof(double));
	for(int i = 0; i < par->L_t * par->L * par->L * par->L * 8; i++) {
		gsl_matrix_complex_memcpy(
            gsl_matrix_complex_view_array(gauged_lattice + 8 * i, 2, 2), 
            gsl_matrix_complex_view_array(lattice + 8 * i, 2, 2)
        );
	}
	gsl_complex *temp = calloc(par->n_su2 * 8, sizeof(double));

	init_su2(par, temp);
	init_gauge(par, gauge, temp);

	for(int a = 0; a < par->L_t; a++){
		for(int b = 0; b < par->L; b++){
			for(int c = 0; c < par->L; c++){
				for(int d = 0; d < par->L; d++){

					for(int dir = 0; dir < 4; dir++){
						
                        msl_mat_mul(
                            0, 
                            0, 
							gsl_matrix_complex_view_array(
                                gauge + 8 * gauge_ind(a, b, c, d, par->L), 
                                2, 
                                2
                            ),
							gsl_matrix_complex_view_array(
                                gauged_lattice + 8 * ind(a, b, c, d, dir, par->L), 
                                2, 
                                2
                            ),
							par->m_workspace);
						msl_mat_mul(
                            0, 
                            1, 
							par->m_workspace,
							gsl_matrix_complex_view_array(
                                gauge + 8 * gauge_periodic_ind(
                                    a + (dir == 0),
                                    b + (dir == 1),
                                    c + (dir == 2),
                                    d + (dir == 3), 
                                    par->L_t, 
                                    par->L
                                ), 
                                2, 
                                2
                            ),
							gsl_matrix_complex_view_array(
                                gauged_lattice + 8 * ind(a, b, c, d, dir, par->L), 
                                2, 
                                2
                            )
                        );
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
                            gsl_matrix_complex_view_array(
                                gauged_lattice + 8 * ind(i, j, k, l, dir_dagger, par->L), 
                                2, 
                                2
                            ), 
                            gsl_matrix_complex_view_array(
                                gauged_lattice + 8 * periodic_ind(
                                    i - (dir_dagger == 4), 
                                    j - (dir_dagger == 5), 
                                    k - (dir_dagger == 6), 
                                    l - (dir_dagger == 7), 
                                    dir_dagger - 4, 
                                    par->L_t, 
                                    par->L
                                ), 
                                2, 
                                2
                            )
                        );
                    }
                }
            }
        }
	}
	double action[2] = {0};
	measure_action_r(par, lattice, action);
	measure_action_r(par, gauged_lattice, action+1);
	free(temp);
	free(gauge);
	free(gauged_lattice);
	return fabs(action[0]-action[1]);
} */

/* Applies a random local gauge transformation onto the lattice. All operators should be invariant under this
 * transformation */
int gauge_transform_lattice(PAR *par, double *lattice){
	double *gauge = calloc(par->L_t * par->L * par->L * par->L * 8, sizeof(double)); 
    if (gauge == NULL) {
        printf("Error: Failed allocating memory for lattice of gauge matrices. Exiting...\n");
        return 1;
    }
	double *su2_gauge = calloc(par->n_su2 * 8, sizeof(double));
    if (su2_gauge == NULL) {
        printf("Error: Failed allocating memory for array of SU(2) matrices for the gauge transformation. Exiting...\n");
        free(gauge);
        return 1;
    }
    gsl_matrix_complex_view m_latt_view, m_gauge_view;

	init_su2(par, su2_gauge);
	init_gauge(par, gauge, su2_gauge);
	for(int a = 0; a < par->L_t; a++){
		for(int b = 0; b < par->L; b++){
			for(int c = 0; c < par->L; c++){
				for(int d = 0; d < par->L; d++){
					
                    for(int dir = 0; dir < 4; dir++){
						
                        m_latt_view = gsl_matrix_complex_view_array(
                            lattice + 8 * ind(a, b, c, d, dir, par->L), 
                            2, 
                            2
                        );
                        m_gauge_view = gsl_matrix_complex_view_array(
                            gauge + 8 * gauge_ind(a, b, c, d, par->L), 
                            2, 
                            2
                        );
                        msl_mat_mul(
                            0, 
                            0, 
							&m_gauge_view.matrix,
							&m_latt_view.matrix,
							par->m_workspace
                        );
                        m_gauge_view = gsl_matrix_complex_view_array(
                            gauge + 8 * gauge_periodic_ind(
                                a + (dir == 0),
                                b + (dir == 1),
                                c + (dir == 2),
                                d + (dir == 3), 
                                par->L_t, 
                                par->L
                            ), 
                            2, 
                            2
                        );
						msl_mat_mul(
                            0, 
                            1, 
							par->m_workspace,
							&m_gauge_view.matrix,
							&m_latt_view.matrix
                        );
					}
				}
			}
		}
	}
	free(su2_gauge);
	free(gauge);
	
    return 0;
}

