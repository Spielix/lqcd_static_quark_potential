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
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_blas.h>   /* Basic Linear Algebra Subprograms */
#include <gsl/gsl_linalg.h>     /* Only for debugging (determinants) */

typedef struct PAR {
    int n_configs;
    int L;
    long int seed;
    gsl_rng_type *gen_type;
    gsl_rng *ran_gen;
    int n_su2;
    double eps;
    int n_hits;
    int n_therm;
    int n_corr;
    double beta;
    gsl_matrix_complex *m_workspace;
} PAR;


int init_su2(PAR *par, gsl_matrix_complex **su2);
void psl_matrix_complex_dagger(gsl_matrix_complex *m);
void psl_matrix_complex_dagger_memcpy(gsl_matrix_complex *dest, gsl_matrix_complex *src);
void psl_matrix_complex_product_3_add(PAR *par,const gsl_matrix_complex *matrix_1, const gsl_matrix_complex *matrix_2, const gsl_matrix_complex *matrix_3, gsl_matrix_complex *m_sum);
int psl_matrix_complex_unitarize(gsl_matrix_complex *matrix);
int measure(PAR *par, double *results, gsl_matrix_complex **lattice);
int measure_aa_a2a(PAR *par, double *results, gsl_matrix_complex **lattice);
int measure_action_l(PAR *par, gsl_matrix_complex **lattice, double *action);
int measure_action_r(PAR *par, gsl_matrix_complex **lattice, double *action);
int unitarize_lattice(PAR *par, gsl_matrix_complex **lattice);
int update_lattice(PAR *par, gsl_matrix_complex **lattice, gsl_matrix_complex **su3, double * acceptance);
int init_su3(PAR *par, gsl_matrix_complex **su3);
void init_lattice(PAR *par, gsl_matrix_complex **lattice, gsl_matrix_complex **su3);
int simulate(PAR *par, gsl_matrix_complex **lattice);
int read_args(PAR *par, char *arg);
int check_su2(gsl_matrix_complex *matrix, gsl_matrix_complex *dagger, double epsilon);
double gauge_inv(PAR *par, gsl_matrix_complex **lattice);
int gauge_transform_lattice(PAR *par, gsl_matrix_complex **lattice);

#endif
