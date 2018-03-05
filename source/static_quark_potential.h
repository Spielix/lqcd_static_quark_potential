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
    int n_su3;
    double eps;
    int n_hits;
    int n_therm;
    int n_corr;
    double beta;
    gsl_matrix_complex *m_temp_1;
} PAR;


int init_su2(PAR *par, gsl_matrix_complex **su2);

#endif
