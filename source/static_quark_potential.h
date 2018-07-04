#ifndef STATIC_QUARK_POTENTIAL_H__
#define STATIC_QUARK_POTENTIAL_H__


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
/* #include <unistd.h>
#include <fcntl.h> */
#include <time.h>
#include <math.h>

#include <gsl/gsl_rng.h>    /* random number generation */
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_blas.h>   /* Basic Linear Algebra Subprograms */
#include <gsl/gsl_linalg.h>     /* Only for debugging (determinants) */

typedef struct PAR {
    int n_configs;
    int L;
    int L_t;
    long int seed;
    gsl_rng_type *gen_type;
    gsl_rng *ran_gen;
    int n_su2;
    double eps;
    int n_hits;
    int n_therm;
    int n_corr;
    double beta;
    double tadpole;
    gsl_matrix_complex *m_workspace;
    double *temp_lat;
    int *temp_lat_filled;
} PAR;


void init_su2(PAR *par, double *su2);
void psl_matrix_complex_dagger(gsl_matrix_complex *m);
void psl_matrix_complex_dagger_memcpy(gsl_matrix_complex *dest, gsl_matrix_complex *src);
void psl_matrix_complex_product_3_add(
    PAR *par, 
    int dag_1, 
    int dag_2, 
    int dag_3, 
    const gsl_matrix_complex *matrix_1, 
    const gsl_matrix_complex *matrix_2, 
    const gsl_matrix_complex *matrix_3, 
    gsl_matrix_complex *m_sum
);
void psl_matrix_complex_product_4(
    PAR *par,
    int dag_1, 
    int dag_2, 
    int dag_3, 
    int dag_4, 
    const gsl_matrix_complex *matrix_1, 
    const gsl_matrix_complex *matrix_2, 
    const gsl_matrix_complex *matrix_3, 
    const gsl_matrix_complex *matrix_4,
    gsl_matrix_complex *m_result
);
double plaquette_op(
    PAR *par,
    double *lattice, 
    int i_start, 
    int j_start, 
    int k_start, 
    int l_start, 
    int dir_1, 
    int dir_2
);
void psl_matrix_complex_unitarize(gsl_matrix_complex *matrix);
int measure_polyakov(PAR *par, double *result, double *lattice, char *file_name);
int measure(PAR *par, gsl_matrix *results, double *lattice, char *file_name);
void measure_action_l(PAR *par, double *lattice, double *action);
void measure_action_r(PAR *par, double *lattice, double *action);
void unitarize_lattice(PAR *par, double *lattice);
void update_lattice(PAR *par, double *lattice, double *su3, double *acceptance);
void init_lattice(PAR *par, double *lattice, double *su3);
int simulate(PAR *par, double *lattice);
int read_args(PAR *par, char *arg);
int check_su2(gsl_matrix_complex *matrix, gsl_matrix_complex *dagger, double epsilon);
int gauge_transform_lattice(PAR *par, double *lattice);

/* calculate the product of a given matrix m_product with n_matrices links along the direction dir from a
 * starting point and save it under m_product*/
void lattice_line_product(
    const PAR *par, 
    double *lattice, 
    int i_start, 
    int j_start, 
    int k_start, 
    int l_start, 
    int dir, 
    int n_matrices, 
    gsl_matrix_complex *m_product
);

/* calculates a planar Wilson-loop */
void lattice_loop_rect(
    const PAR *par, 
    double *lattice, 
    int i_start, 
    int j_start, 
    int k_start, 
    int l_start, 
    int dir_1, 
    int dir_2, 
    int L_1, 
    int L_2, 
    double *result
);

#endif
