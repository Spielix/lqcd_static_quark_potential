#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h> /* For RANLXS random number generation */

#include "static_quark_potential.h"


void measure(Par *par, double *results, int *spin) {
}


void print_results(Par *par, double *results) {
}


void init_spin(Par *par, int *spin) {
    int i, L2 = par->L * par->L;

    /* Setting random Spin values on the lattice */
    for (i = 0; i < L2; i++) {
        /* Getting a random number from gsl_rng */
        spin[i] = (int)(2. * gsl_rng_uniform(par->ran_gen));
        if (!spin[i]) spin[i] = -1;
    }
}


int sim(Par *par, int *spin) {
    double results [2];
    int i;

    for (i = 0; i < par->nconfigs; i++) {
        init_spin(par, spin);
        measure(par, results, spin);
        print_results(par, results);
    }
    return 0;
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
