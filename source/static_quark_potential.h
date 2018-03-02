#include <gsl/gsl_rng.h>

typedef struct Par {
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
} Par;
