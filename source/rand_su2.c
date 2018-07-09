#include "static_quark_potential.h"

int sign(double x);

void init_su2(const PAR *par, double *su2){
    double norm,
    x_0 = sqrt(1. - par->eps * par->eps);

    for (int i = 0; i < par->n_su2 / 2; i++) {
        norm = 0.;
        su2[4 * i] = x_0;
        for (int j = 1; i < 4; i++) {
            su2[4 * i + j] = gsl_rng_uniform_pos(par->ran_gen) - 0.5;
            norm += su2[4 * i + j] * su2[4 * i + j];
        }
        norm = sqrt(norm);
        for (int j = 1; j < 4; j++) 
            su2[4 * i + j] *= par->eps / norm;
        psl_su2_dagger_memcpy(su2 + par->n_su2 / 2 + 4 * i, su2 + 4 * i);
    }	
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

void init_gauge(const PAR *par, double *gauge, const double *su2) {
    int rand_ind;

    /* initialize all independent links with random SU(2) matrices */
    for (int i = 0; i < par->L_t; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    
                    rand_ind = (int)(gsl_rng_uniform(par->ran_gen) * par->n_su2);
                    
                    for (int m = 0; m < 4; m++) {
                        gauge[4 * gauge_ind(i, j, k, l, par->L) + m] = su2[4 * rand_ind + m];
                    }
                }
            }
        }
    }
}

/* Applies a random local gauge transformation onto the lattice. All operators should be invariant under this
 * transformation */
int gauge_transform_lattice(const PAR *par, double *lattice){
    double m_temp[4] = {0};
	double *gauge = calloc(par->L_t * par->L * par->L * par->L * 4, sizeof(double)); 
    if (gauge == NULL) {
        printf("Error: Failed allocating memory for lattice of gauge matrices. Exiting...\n");
        return 1;
    }
	double *su2_gauge = calloc(par->n_su2 * 4, sizeof(double));
    if (su2_gauge == NULL) {
        printf("Error: Failed allocating memory for array of SU(2) matrices for the gauge transformation. Exiting...\n");
        free(gauge);
        return 1;
    }

	init_su2(par, su2_gauge);
	init_gauge(par, gauge, su2_gauge);
	for(int a = 0; a < par->L_t; a++){
		for(int b = 0; b < par->L; b++){
			for(int c = 0; c < par->L; c++){
				for(int d = 0; d < par->L; d++){
					
                    for(int dir = 0; dir < 4; dir++){
						
                        psl_su2_product_notrans_notrans(
							gauge + 4 * gauge_ind(a, b, c, d, par->L), 
							lattice + 4 * ind(a, b, c, d, dir, par->L), 
							m_temp
                        );
                        
						psl_su2_product_notrans_conjtrans(
							m_temp,
							gauge + 4 * gauge_periodic_ind(
                                a + ((dir == 0) ? 1 : 0),
                                b + ((dir == 1) ? 1 : 0),
                                c + ((dir == 2) ? 1 : 0),
                                d + ((dir == 3) ? 1 : 0), 
                                par->L_t, 
                                par->L
                            ),
							lattice + 4 * ind(a, b, c, d, dir, par->L)
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

