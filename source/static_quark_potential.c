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

#include "static_quark_potential.h"


void measure(Par *par, double *results, int *spin) {
}


void print_results(Par *par, double *results) {
}


int init_su3(Par *par, gsl_matrix_complex **su3) {
    int n, i, j;
    gsl_vector_complex *vec[3], *tempv;
    gsl_complex tempz;

    /* allocating memory for three vectors which will temporarily store the column vectors of the matrices for
     * unitarization */
    for (i = 0; i < 3; i++) {
        vec[i] = gsl_vector_complex_alloc(3);
    }
    tempv = gsl_vector_complex_alloc(3);
    if ((vec[0] == NULL) || (vec[1] == NULL) || (vec[2] == NULL) || (tempv == NULL)) {
        printf("Error: Allocating memory for vector in init_su3 failed. Exiting...\n");
        return 1;
    }

    for (n = 0; n < par->num_su3 / 2; n++) {
        /* initialize random hermitian matrix */
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                gsl_matrix_complex_set(
                    su3[n], 
                    i, 
                    j, 
                    gsl_complex_rect((i == j), par->eps * (gsl_rng_uniform_pos(par->ran_gen) * 2. - 1.))
                );
            }
        }
        /* unitarize the matrix (Gram-Schmidt) */
        for (i = 0; i < 3; i++)
            gsl_matrix_complex_get_col(vec[i], su3[n], i);

        gsl_blas_zdscal(1. / gsl_blas_dznrm2(vec[0]), vec[0]);     /* normalization */

        gsl_blas_zdotu(vec[0], vec[1], &tempz);
        gsl_vector_complex_memcpy(tempv, vec[0]);
        gsl_vector_complex_scale(tempv, tempz);
        gsl_vector_complex_sub(vec[1], tempv);
        gsl_blas_zdscal(1. / gsl_blas_dznrm2(vec[1]), vec[1]);     /* normalization */

        gsl_blas_zdotu(vec[0], vec[2], &tempz);
        gsl_vector_complex_memcpy(tempv, vec[0]);
        gsl_vector_complex_scale(tempv, tempz);
        gsl_vector_complex_sub(vec[2], tempv);
        gsl_blas_zdotu(vec[1], vec[2], &tempz);
        gsl_vector_complex_memcpy(tempv, vec[1]);
        gsl_vector_complex_scale(tempv, tempz);
        gsl_vector_complex_sub(vec[2], tempv);
        gsl_blas_zdscal(1. / gsl_blas_dznrm2(vec[2]), vec[2]);     /* normalization */

        for (i = 0; i < 3; i++)
            gsl_matrix_complex_set_col(su3[n], i, vec[i]);
    }

    /* create inverse matrices by transposing and taking the complex conjugate of each element */
    for (n = par->num_su3 / 2; n < par->num_su3; n++) {
        gsl_matrix_complex_transpose_memcpy(su3[n], su3[n - par->num_su3 / 2]);
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                gsl_matrix_complex_set(
                    su3[n], 
                    i, 
                    j, 
                    gsl_complex_conjugate(gsl_matrix_complex_get(su3[n], i, j))
                );
            }
        }
    }

    for (i = 0; i < 3; i++)
        gsl_vector_complex_free(vec[i]);
    gsl_vector_complex_free(tempv);

    /* DEBUG: 
    for (n = 0; n < par->num_su3 / 2; n++) {
        
    } */

    return 0;
}


int sim(Par *par, int *spin) {
    double results [2];
    int i, succ;
    
    gsl_matrix_complex **su3;

    for (i = 1; i < par->num_su3; i++) {
        su3[i] = gsl_matrix_complex_alloc(3, 3);
        if (su3[i] == NULL) {
            printf("Error: failed allocating memory for su3 matrices. Exiting...\n");
            return 1;
        }
    }
    
    succ = !init_su3(par, su3);

    for (i = 0; i < par->num_su3; i++)
        gsl_matrix_complex_free(su3[i]);

    return !succ;
}


int read_args(Par *par, char *arg) {
    static int *spin = NULL;
    int success = 1;
    char *s;
    
    /* Running the simulation with already given parameters */
    if (!strcmp(arg, "run")) {
        if (!par->L) {
            printf("Give system size L!\n");
            return 1;
        }
        
        /* Initialization of randum number generator */
        par->ran_gen = gsl_rng_alloc(par->gen_type);
        if (par->ran_gen == NULL) {
            if (spin) free(spin);
            printf("Error: Failed allocating memory for the random number generator! Exiting...\n");
            return 1;
        }
        if (par->seed) gsl_rng_set(par->ran_gen, par->seed);
        else gsl_rng_set(par->ran_gen, (long)time(NULL));
        
        success = !sim(par, spin);
        
        gsl_rng_free(par->ran_gen);
        if (spin) free(spin);

        return !success;
    }

    /* interpretation of simulation parameters */
    s = strchr(arg, '=');

    if (!s) {
        fprintf(stderr, "Command '%s' not recognized, expected format: '<name>=<value>'. Exiting...\n", arg);
        return 1;
    }

    *s++ = '\0';

    if (!strcmp(arg, "L")) {
        int L2, tempL;
   
        /* checking if the given lattice size is a power of 2 */
        tempL = strtod(s, NULL);
        if ((tempL > 0) && !(tempL & (tempL - 1))) par->L = tempL;
        else {
            fprintf(stderr, "L has to be a power of 2. Exiting...\n");
            return 1;
        }

        L2 = par->L * par->L;
        
        spin = realloc(spin, L2 * sizeof(int));
        if (spin == NULL) {
            printf("Error: Failed allocating memory for ising_lattice! Exiting...\n");
            return 1;
        }
        return 0;
    }

    if (!strcmp(arg, "seed")) {
        par->seed = strtol(s, NULL, 0);
        return 0;
    }
    
    if (!strcmp(arg, "nconfigs")) {
        par->nconfigs = strtol(s, NULL, 0);
        return 0;
    }
    
    if (!strcmp(arg, "num_su3")) {
        par->num_su3 = strtol(s, NULL, 0);
        return 0;
    }
    
    if (!strcmp(arg, "eps")) {
        par->eps = strtod(s, NULL);
        return 0;
    }

/*    if (!strcmp(arg, "gen_type")) {
        par->gen_type = strtol(s, NULL, 0);
        return 0;
    } */

    fprintf(stderr, "No such variable name: '%s'. Exiting...\n", arg);
    return 1;
}


int main(int argc, char *argv[])
{
    int iarg;
    Par *par = malloc(sizeof(Par));

    par->L = 0;
    par->seed = 0;
    par->nconfigs = 1;
    par->gen_type = gsl_rng_ranlxs0;
    par->num_su3 = 100;
    par->eps = 1.;
  
    if (argc == 1) {
        printf("Usage: %s L=16 nconfigs=100 run\n", argv[0]);
        printf("Optional arguments (with defaults) L=%d seed=%ld", par->L, par->seed);
        free(par);
        exit(EXIT_SUCCESS);
    }
    
    /* read_args interprets the arguments given to the program and starts it when "run" appears */
    for (iarg = 1; iarg < argc; iarg++)
        if (read_args(par, argv[iarg])) {
            free(par);
            exit(EXIT_FAILURE);
        }

    free(par);
}
