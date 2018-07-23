#ifndef STATIC_QUARK_POTENTIAL_H__
#define STATIC_QUARK_POTENTIAL_H__


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
/* #include <unistd.h>
#include <fcntl.h> */
#include <time.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>    /* random number generation */
//#include <gsl/gsl_complex.h>
//#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_double.h>
//#include <gsl/gsl_matrix_complex_double.h>
//#include <gsl/gsl_vector_complex_double.h>
//#include <gsl/gsl_blas.h>   /* Basic Linear Algebra Subprograms */
//#include <gsl/gsl_linalg.h>     /* Only for debugging (determinants) */

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
    double *temp_lat;
    int *temp_lat_filled;
} PAR;

void psl_su2_memcpy(double *dest, const double *src);
void psl_su2_dagger(double *m);
void psl_su2_dagger_memcpy(double *dest, const double *src);
void psl_su2_product_notrans_notrans(const double *m_1, const double *m_2, double *m_3);
void psl_su2_product_conjtrans_conjtrans(const double *m_1, const double *m_2, double *m_3);
void psl_su2_product_notrans_conjtrans(const double *m_1, const double *m_2, double *m_3);
void psl_su2_product_conjtrans_notrans(const double *m_1, const double *m_2, double *m_3);

void init_su2(const PAR *par, double *su2);
void psl_su2_product_3_add(
    int dag_1, 
    int dag_2, 
    int dag_3, 
    const double *m_1, 
    const double *m_2, 
    const double *m_3, 
    double *m_sum
);
void psl_su2_product_4(
    int dag_1, 
    int dag_2, 
    int dag_3, 
    int dag_4, 
    const double *m_1, 
    const double *m_2, 
    const double *m_3, 
    const double *m_4,
    double *m_result
);
double plaquette_op(
    const PAR *par,
    const double *lattice, 
    int i_start, 
    int j_start, 
    int k_start, 
    int l_start, 
    int dir_1, 
    int dir_2
);
void psl_su2_unitarize(double *matrix);
int measure_polyakov(const PAR *par, double *result, const double *lattice, const char *file_name);
int measure(const PAR *par, gsl_matrix *results, const double *lattice, const char *file_name);
void measure_action_l(const PAR *par, const double *lattice, double *action);
void measure_action_r(const PAR *par, const double *lattice, double *action);
void unitarize_lattice(const PAR *par, double *lattice);
void metropolis_update_lattice(const PAR *par, double *lattice, const double *su3, double *acceptance);
void overrelaxation_update_lattice(const PAR *par, double *lattice, const double *su3);
void init_lattice(const PAR *par, double *lattice, const double *su3);
int simulate(PAR *par, double *lattice);
int read_args(PAR *par, char *arg);
int gauge_transform_lattice(const PAR *par, double *lattice);

/* calculate the product of a given matrix m_product with n_matrices links along the direction dir from a
 * starting point and save it under m_product*/
void lattice_line_product(
    const PAR *par, 
    const double *lattice, 
    int i_start, 
    int j_start, 
    int k_start, 
    int l_start, 
    int dir, 
    int n_matrices, 
    double *m_product
);

void wilson_line_product(
    const PAR *par, 
    const double *lattice, 
    int i_start, 
    int j_start, 
    int k_start, 
    int l_start, 
    int dir, 
    int n_matrices, 
    double *m_product
);

/* calculates a planar Wilson-loop */
void lattice_loop_rect(
    const PAR *par, 
    const double *lattice, 
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
