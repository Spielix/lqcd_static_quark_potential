#include "static_quark_potential.h"

// #define NO_BIG_LOOPS

int main(int argc, char *argv[]) {
    PAR par;
    FILE *data_file;

    data_file = fopen(argv[1], "r");
    if (data_file == NULL) {
        printf("Failed opening data file. Exiting...\n");
        exit(EXIT_FAILURE);
    }

    if (fread(&par, sizeof(PAR), 1, data_file) != 1) {
        printf("Failed reading parameters from file. Exiting...\n");
        fclose(data_file);
        exit(EXIT_FAILURE);
    }

    printf("The simulation parameters where:\nL_t = %d, L = %d\nbeta = %g, epsilon = %g\nseed = %ld\nn_samples = %d, n_therm = %d, n_corr = %d, n_su2 = %d, n_hits = %d\n", 
        par.L_t, par.L, par.beta, par.eps, par.seed, par.n_configs, par.n_therm, par.n_corr, par.n_su2, par.n_hits);

    fclose(data_file);

    exit(EXIT_SUCCESS);
}
