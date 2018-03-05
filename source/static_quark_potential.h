#ifndef STATIC_QUARK_POTENTIAL_H__
#define STATIC_QUARK_POTENTIAL_H__

#include <gsl/gsl_rng.h>

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
} PAR;

#endif
