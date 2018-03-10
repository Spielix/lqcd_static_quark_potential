#include "static_quark_potential.h"

int main(int argc, char *argv[]) {
    PAR par_old, par_new;
    FILE *data_file_old, *data_file_new;
    int counter = 0;
    char s_temp[1000];

    data_file_new = fopen(argv[1], "w");
    if (data_file_new == NULL) {
        printf("Failed creating/opening new data file. Exiting...\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 2; i < argc; i++) {
        data_file_old = fopen(argv[i], "r");
        if (data_file_old == NULL) {
            printf("Failed opening old data file named %s. Exiting...\n", argv[i]);
            fclose(data_file_new);
            if (remove(argv[1])) 
                printf("Warning: Failed removing unfinished data file.\n");
            exit(EXIT_FAILURE);
        }

        if (fread(&par_old, sizeof(PAR), 1, data_file_old) != 1) {
            printf("Failed reading parameters from file %s. Exiting...\n", argv[i]);
            fclose(data_file_old);
            fclose(data_file_new);
            if (remove(argv[1])) 
                printf("Warning: Failed removing unfinished data file.\n");
            exit(EXIT_FAILURE);
        }
        if (i == 2) {
            par_new = par_old;
            if (fwrite(&par_new, sizeof(PAR), 1, data_file_new) != 1) {
                printf("Failed writing parameter struct to new data file. Exiting...\n");\
                fclose(data_file_old);
                fclose(data_file_new);
                if (remove(argv[1])) 
                    printf("Warning: Failed removing unfinished data file.\n");
                exit(EXIT_FAILURE);
            }
        } else {
            if ((par_old.L != par_new.L) || (par_old.L_t != par_new.L_t)) {
                printf("Data files with differing lattice size detected. Exiting...\n");
                fclose(data_file_old);
                fclose(data_file_new);
                if (remove(argv[1])) 
                    printf("Warning: Failed removing unfinished data file.\n");
                exit(EXIT_FAILURE);
            }
            if (par_old.seed == par_new.seed) {
                printf("Data files with same seed detected. Exiting...\n");
                fclose(data_file_old);
                fclose(data_file_new);
                if (remove(argv[1])) 
                    printf("Warning: Failed removing unfinished data file.\n");
                exit(EXIT_FAILURE);
            }
            if (
                (par_old.beta != par_new.beta) || 
                (par_old.tadpole != par_new.tadpole) || 
                (par_old.eps != par_new.eps)
            ) {
                printf("Data files with differing beta / tadpole / eps detected. Exiting...\n");
                fclose(data_file_old);
                fclose(data_file_new);
                if (remove(argv[1])) 
                    printf("Warning: Failed removing unfinished data file.\n");
                exit(EXIT_FAILURE);
            }
            if ((par_old.n_corr != par_new.n_corr) || (par_old.n_therm != par_new.n_therm)) 
                printf("Warning: Data files with differing n_corr / n_therm detected. Continuing...\n");
        }
        while (!feof(data_file_old)) {
            double data;
            if (fread(&data, sizeof(double), 1, data_file_old) != 1) {
                printf("Failed reading data from data file. Exiting...\n");
                fclose(data_file_old);
                fclose(data_file_new);
                if (remove(argv[1])) 
                    printf("Warning: Failed removing unfinished data file.\n");
                exit(EXIT_FAILURE);
            }
            if (fwrite(&data, sizeof(double), 1, data_file_new) != 1) {
                printf("Failed writing data to file. Exiting...\n");
                fclose(data_file_old);
                fclose(data_file_new);
                if (remove(argv[1])) 
                    printf("Warning: Failed removing unfinished data file.\n");
                exit(EXIT_FAILURE);
            }
            counter++;
        }
        fclose(data_file_old);
    }

    par_new.n_configs = counter;
    rewind(data_file_new);

    if (fwrite(&par_new, sizeof(PAR), 1, data_file_new) != 1) {
        printf("Failed writing  new n_configs to data file. Exiting...\n");
        fclose(data_file_new);
        if (remove(argv[1])) 
            printf("Warning: Failed removing unfinished data file.\n");
        exit(EXIT_FAILURE);
    }

    fclose(data_file_new);

    for (int i = 2; i < argc; i++) {
        sprintf(s_temp, "%s_went_into_%s", argv[i], argv[1]);
        if (rename(argv[i], s_temp)) {
            printf("Warning: Failed renaming file %s.\n", argv[i]);
        }
    }

    exit(EXIT_SUCCESS);
}
