#include <gsl/gsl_rng.h>

typedef struct Par {
    int nconfigs;
    int L;
    long int seed;
    gsl_rng_type *gen_type;
    gsl_rng *ran_gen;
    int num_su3;
    double eps;
} Par;
