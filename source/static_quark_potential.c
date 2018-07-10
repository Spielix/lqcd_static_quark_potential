#include "static_quark_potential.h"
// #define TADP
// #define GAUGE
#define POLYA
// #define POLYA_GAUGE
// #define DEBUG

static inline int ind(int i, int j, int k, int l, int dir, int le);
static inline int periodic_ind(int i, int j, int k, int l, int dir, int le_0, int le);
static inline int temp_ind(int i, int j, int k, int l, int dir, int ll, int le_max, int le);
static inline int temp_periodic_ind(int i, int j, int k, int l, int dir, int ll, int le_max, int le_0, int le);

int main(int argc, char *argv[])
{
    PAR *par = NULL;

    /* allocate memory for simulation parameters */
     par = malloc(sizeof(PAR));
     if (par == NULL) {
        printf("Error: Failed allocating memory for simulation parameters. Exiting...\n");
        exit(EXIT_FAILURE);
     }

    
    /* initialize simulation parameters with standard values */
    par->L = 0;
    par->L_t = 0;
    par->seed = 0;
    par->n_configs = 10;
    par->gen_type = (gsl_rng_type *)gsl_rng_mt19937;
    par->n_su2 = 10000;
    par->eps = 0.4;
    par->n_hits = 10;
    par->n_therm = 500;
    par->n_corr = 50;
    par->beta = 2.3;
    par->tadpole = 1.;

    if (argc == 1) {
        printf("Usage: %s L=16 nconfigs=100 run\n", argv[0]);
        printf("Optional arguments (with defaults) L=%d seed=%ld", par->L, par->seed);
        free(par);
        exit(EXIT_SUCCESS);
    }
    
    /* read_args interprets the arguments given to the program and starts it when "run" appears */
    for (int iarg = 1; iarg < argc; iarg++)
        if (read_args(par, argv[iarg])) {
            free(par);
            exit(EXIT_FAILURE);
        }

    free(par);
    exit(EXIT_SUCCESS);
}

int read_args(PAR *par, char *arg) {
    static double *lattice = NULL;
    int success = 1;
    char *s;
    
    if (!strcmp(arg, "run")) {
        if ((!par->L) || (!par->L_t)) {
            printf("Give system size L and L_t!\n");
            return 1;
        }
        if (par->beta == 0.) {
            printf("Give beta!\n");
            return 1;
        }
        /* preparing for the simulation */
        /* initialization of random number generator */
        par->ran_gen = gsl_rng_alloc(par->gen_type);
        if (par->ran_gen == NULL) {
            if (lattice) 
                free(lattice);
            printf("Error: Failed allocating memory for the random number generator! Exiting...\n");
            return 1;
        }
        if (!(par->seed)) 
            par->seed = (long)time(NULL);
        gsl_rng_set(par->ran_gen, par->seed);
        
        /* running the simulation with already given parameters */
        success = !simulate(par, lattice);
        
        gsl_rng_free(par->ran_gen);
        
        free(lattice);
        free(par->temp_lat);
        free(par->temp_lat_filled);

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
        unsigned long num_of_links, 
            temp_lat_size;
   
        /* checking if the given lattice size is a power of 2 */
        par->L = strtod(s, NULL);

        /* if ((tempL > 0) && !(tempL & (tempL - 1))) par->L = tempL;
        else {
            fprintf(stderr, "L has to be a power of 2. Exiting...\n");
            return 1;
        } */
        
        /* allocate memory for the lattice (link matrices) */
        if (par->L_t) {
            num_of_links = (unsigned long)(par->L_t * par->L * par->L * par->L) * (unsigned long)4;
            if (par->L_t > par->L) 
                temp_lat_size = num_of_links * (unsigned long)par->L_t;
            else 
                temp_lat_size = num_of_links * (unsigned long)par->L;
            
            lattice = malloc(num_of_links * 4 * sizeof(double));
            par->temp_lat = malloc(temp_lat_size * 4 * sizeof(double));
            par->temp_lat_filled = malloc(temp_lat_size * sizeof(int));

            if ((lattice == NULL) || (par->temp_lat == NULL) || (par->temp_lat_filled == NULL)) {
                printf("Error: Failed allocating memory for the lattice and/or temporary lattice %p %p %p! Exiting...\n", (void *)lattice, (void *)par->temp_lat, (void *)par->temp_lat_filled);
                return 1;
            }
        }
        return 0;
    }

    if (!strcmp(arg, "L_t")) {
        unsigned long num_of_links, 
            temp_lat_size;
   
        /* checking if the given lattice size is a power of 2 */
        par->L_t = strtod(s, NULL);

        /* if ((tempL > 0) && !(tempL & (tempL - 1))) par->L = tempL;
        else {
            fprintf(stderr, "L has to be a power of 2. Exiting...\n");
            return 1;
        } */
        
        /* allocate memory for the lattice (link matrices) */
        if (par->L) {
            num_of_links = (unsigned long)(par->L_t * par->L * par->L * par->L) * (unsigned long)4;
            if (par->L_t > par->L) 
                temp_lat_size = num_of_links * (unsigned long)par->L_t;
            else 
                temp_lat_size = num_of_links * (unsigned long)par->L;
            
            lattice = malloc(num_of_links * 4 * sizeof(double));
            par->temp_lat = malloc(temp_lat_size * 4 * sizeof(double));
            par->temp_lat_filled = malloc(temp_lat_size * sizeof(int));

            if ((lattice == NULL) || (par->temp_lat == NULL) || (par->temp_lat_filled == NULL)) {
                printf(
                    "Error: Failed allocating memory for the lattice and/or temporary lattice %ld %ld %ld! Exiting...\n", 
                    (long int)lattice, 
                    (long int)par->temp_lat, 
                    (long int)par->temp_lat_filled
                );
                return 1;
            }
        }
        return 0;
    }

    if (!strcmp(arg, "seed")) {
        par->seed = strtol(s, NULL, 0);
        return 0;
    }
    
    if (!strcmp(arg, "n_configs")) {
        par->n_configs = strtol(s, NULL, 0);
        return 0;
    }
    
    if (!strcmp(arg, "n_su2")) {
        par->n_su2 = strtol(s, NULL, 0);
        if (par->n_su2 % 2) 
            par->n_su2 = par->n_su2 + 1;
        return 0;
    }
    
    if (!strcmp(arg, "eps")) {
        par->eps = strtod(s, NULL);
        return 0;
    }

    if (!strcmp(arg, "beta")) {
        par->beta = strtod(s, NULL);
        return 0;
    }

    if (!strcmp(arg, "n_hits")) {
        par->n_hits = strtol(s, NULL, 0);
        return 0;
    }

    if (!strcmp(arg, "n_therm")) {
        par->n_therm = strtol(s, NULL, 0);
        return 0;
    }

    if (!strcmp(arg, "n_corr")) {
        par->n_corr = strtol(s, NULL, 0);
        return 0;
    }

    if (!strcmp(arg, "tadpole")) {
        par->tadpole = strtod(s, NULL);
        return 0;
    }

/*    if (!strcmp(arg, "gen_type")) {
        par->gen_type = strtol(s, NULL, 0);
        return 0;
    } */

    fprintf(stderr, "No such variable name: '%s'. Exiting...\n", arg);
    return 1;
}

int simulate(PAR *par, double *lattice) {
    double acceptance = 0.;
    double *su2 = NULL;
#ifdef POLYA
    double polya_temp_res, polya_res = 0;
#endif
#ifdef TADP
    double tadpole_result[100];
#endif
   
    /* prepare data file */
    char file_name[100];
    FILE *data_file;
    sprintf(file_name, "data/L%d_Lt%d_beta%3.2f_u%3.2f_seed%ld.bin", par->L, par->L_t, par->beta, par->tadpole, par->seed);
    data_file = fopen(file_name, "w");
    if (data_file == NULL) {
        printf("Error: Failed opening data file for preparation. Exiting...\n");
        return 1;
    }
    if (fwrite(par, sizeof(PAR), 1, data_file) == -1) {
        printf("Failed writing parameters to file. Exiting...\n");
        fclose(data_file);
        return 1;
    }
    fclose(data_file);

    /* allocate memory for a specified number of random SU(2) matrices */
    su2 = calloc(par->n_su2 * 4, sizeof(double));
    if (su2 == NULL) {
        printf("Error: Failed allocating memory for array of SU(2)-matrix addresses. Exiting...\n");
        return 1;
    }
#ifndef TADP
#ifndef POLYA
    double *results = NULL;
    gsl_matrix_view m_result_view;
    /* allocate memory for a specified number of result-matrices */
    results = malloc(par->n_configs * (par->L_t - 1) * (par->L - 1) * sizeof(double));
    if (results == NULL) {
        printf("Error: Failed allocating memory for array of results. Exiting...\n");
        free(su2);
        return 1;
    }
#endif
#endif
    
    printf(
        "Starting Lattice-QCD simulation with parameters: L=%d, seed=%ld, eps=%4.3f, beta=%4.3f, tadpole=%4.3f\n", 
        par->L, 
        par->seed,
        par->eps,
        par->beta, 
        par->tadpole
    );

    /* generate the SU(2) matrices */
    printf("Generating %d random SU(2)-matrices...\n", par->n_su2);
    init_su2(par, su2);
    //for (int jj = 0; jj < 4 * par->n_su2; jj++) 
    //    printf("%d\t%g\n", jj % 4, su2[jj]);

    /* DEBUG: 
    for (int i = 80; i < 88; i++) {
        printf("%f\t", su2[i]);
    }
    printf("\n");
    
    gsl_complex z_temp;
    const gsl_complex z_zero = gsl_complex_rect(0., 0.), 
                    z_one = gsl_complex_rect(1., 0.);
    gsl_matrix_complex *m_temp = NULL, 
                       *m_temp_0 = NULL;
    gsl_permutation *p = NULL;
    
    m_temp = gsl_matrix_complex_alloc(2, 2);
    m_temp_0 = gsl_matrix_complex_alloc(2, 2);
    p = gsl_permutation_alloc(2);
    for (int i = 0; i < par->n_su2 / 2; i++) {
        int signum;
        gsl_complex det;
        
        gsl_matrix_complex_memcpy(m_temp, su2[i]);
        gsl_linalg_complex_LU_decomp(m_temp, p, &signum);
        det = gsl_linalg_complex_LU_det(m_temp, signum);
        printf("DEBUG: det = %4.3f + %4.3fi\n", GSL_REAL(det), GSL_IMAG(det));
        printf("DEBUG: abs(det) = %4.3f\n", gsl_complex_abs(det));
    }
    for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 2; k++) {
            z_temp = gsl_matrix_complex_get(su2[5], j, k);
            printf("%3.2f+%3.2fi\t", GSL_REAL(z_temp), GSL_IMAG(z_temp));
        }
        printf("\n");
    }
    for (int i = 0; i < par->n_su2 / 2; i++) {
        gsl_blas_zgemm(
            CblasNoTrans, 
            CblasNoTrans, 
            z_one, 
            su2[i], 
            su2[i + par->n_su2 / 2], 
            z_zero,
            m_temp
        );
        printf("DEBUG:\n");
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                z_temp = gsl_matrix_complex_get(m_temp, j, k);
                printf("%3.2f+%3.2fi\t", GSL_REAL(z_temp), GSL_IMAG(z_temp));
            }
            printf("\n");
        }
    }
    
    gsl_matrix_complex_memcpy(m_temp, su2[(int)(par->n_su2 * gsl_rng_uniform(par->ran_gen))]);
    for (int i = 0; i < par->n_hits; i++) {
        gsl_blas_zgemm(
            CblasNoTrans, 
            CblasNoTrans, 
            z_one, 
            su2[(int)(par->n_su2 * gsl_rng_uniform(par->ran_gen))],
            m_temp, 
            z_zero, 
            m_temp_0
        );
        gsl_matrix_complex_memcpy(m_temp, m_temp_0);
    }
    z_temp = z_zero;
    for (int i = 0; i < 2; i++) 
        z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp_0, i, i));
    printf("DEBUG: Trace = %e\n", GSL_REAL(z_temp)); 
    gsl_matrix_complex_free(m_temp);
    gsl_matrix_complex_free(m_temp_0);
    gsl_permutation_free(p);
     DEBUG_END */

    /* initialize all links in the lattice */
    printf("Initializing %d x %d^3 lattice...\n", par->L_t, par->L);
    init_lattice(par, lattice, su2);
    
    /* thermalize the lattice */
    printf("Thermalizing (%d sweeps)...\n", par->n_therm);
    for (int i = 0; i < par->n_therm; i++) {
        update_lattice(par, lattice, su2, &acceptance);
        if (i % 1000 == 0) 
            unitarize_lattice(par, lattice);
        /* if (check_su2(lattice[ind(3, 3, 3, 3, 0, par->L)], lattice[ind(4, 3, 3, 3, 4)], 1e-10)) {
            printf("DEBUG; %d\n", i);
            return 0;
        } */
    }
    acceptance /= (double)par->n_therm * 
        (double)par->L_t * pow((double)par->L, 3.) * 
        4. * (double)par->n_hits;
    printf("Thermalization-acceptance: %3.2f\n", acceptance);
    acceptance = 0.;

    /* take measurements after every n_corr-th update-sweep */
    printf(
        "Taking measurements of %d configurations(n_corr=%d, n_hits=%d):\n",
        par->n_configs, 
        par->n_corr, 
        par->n_hits
    );

    printf("L_r x L_t\n");
    for (int i = 0, index = 0; i < par->L - 1; i++) {
        for (int j = 0; j < par->L_t - 1; j++) {
            printf("%dx%d\t\t", i + 1, j + 1);
            index++;
            if (index == 12) 
                break;
        }
        if (index == 12) 
            break;
    }
    printf("\n");
    /* measure_tadpole_alt(par, lattice, tadpole_result);
    printf("%12e\n", tadpole_result[0]);
    tadpole_result[0] = 0; */
    
    for (int i = 0; i < par->n_configs; i++) {
        
        if ((i + 1) % 10 == 0) 
            init_su2(par, su2);
        
        for (int j = 0; j < par->n_corr; j++) {
            update_lattice(par, lattice, su2, &acceptance);

        }
        
        /* DEBUG: measure the action of the lattice with different direction in each Wilson-loop
        double action_1, action_2;
        if (measure_action_r(par, lattice, &action_1)) {
            for (int j = 0; j < par->n_su2; j++) 
                gsl_matrix_complex_free(su2[j]);
            free(su2);
            return 1;
        }
        if (measure_action_l(par, lattice, &action_2)) {
            for (int j = 0; j < par->n_su2; j++) 
                gsl_matrix_complex_free(su2[j]);
            free(su2);
            return 1;
        }
        printf("DEBUG: action_l = %7.2g\nDEBUG: action_r = %7.2g\n", action_1, action_2);
        DEBUG END */

        unitarize_lattice(par, lattice);
#ifndef TADP      
#ifndef POLYA
        m_result_view = gsl_matrix_view_array(
            results + (par->L_t - 1) * (par->L - 1) * i, 
            par->L_t - 1, 
            par->L - 1
        );
        if (measure(
            par, 
            &m_result_view.matrix, 
            lattice, 
            file_name
        )) {
            free(su2);
            free(results);
            return 1;
        }
        for (int k = 0, index = 0; k < par->L - 1; k++) {
            for (int l = 0; l < par->L_t - 1; l++) {
                printf("%7.3e\t", gsl_matrix_get(
                    &m_result_view.matrix,
                    l, 
                    k
                ));
                index++;
                if (index == 12) 
                    break;
            }
            if (index == 12) 
                break;
        }
        printf("\n");

#endif
#endif
#ifdef POLYA
        if (measure_polyakov(par, &polya_temp_res, lattice, file_name)) {
            free(su2);
            return 1;
        }
        printf("%e\n", polya_temp_res);
        polya_res += polya_temp_res;
#endif
#ifdef TADP
        /* measure tadpole and print out the results */
        measure_tadpole_alt(par, lattice, tadpole_result + i);
        printf("%12e\n", tadpole_result[i]);
#endif
#ifdef GAUGE
        /* printf("%g\n", gauge_inv(par, lattice)); */

        /* DEBUG: gauge transform the lattice and look at the effect on measurements  
        double action;
        measure_action_r(par, lattice, &action);
        printf("action = %7.2e\n", action); */

        if (gauge_transform_lattice(par, lattice)) {
            free(su2);
            free(results);
            return 1;
        }
        if (measure(
            par, 
            &m_result_view.matrix, 
            lattice, 
            "plz_rm_this_gauge_stuff"
        )) {
            free(su2);
            free(results);
            return 1;
        }
        for (int k = 0, index = 0; k < par->L - 1; k++) {
            for (int l = 0; l < par->L_t - 1; l++) {
                printf("%7.3e\t", gsl_matrix_get(
                    &m_result_view.matrix, 
                    l, 
                    k
                ));
                index++;
                if (index == 12) 
                    break;
            }
            if (index == 12) 
                break;
        }
        printf("g\n");
#endif
#ifdef POLYA_GAUGE
        if (gauge_transform_lattice(par, lattice)) {
            free(su2);
            return 1;
        }
        if (measure_polyakov(par, &polya_temp_res, lattice, "plz_rm_this_gauge_stuff")) {
            free(su2);
            return 1;
        }
        printf("%e\tg\n", polya_temp_res);
#endif
    }

#ifdef TADP
    double tadpole_avg = 0.;
    for (int i = 0; i < par->n_configs; i++) {
        tadpole_avg += tadpole_result[i];
    }
    tadpole_avg /= par->n_configs;
    printf("Result: <P> = %10.5e\n", tadpole_avg);
#endif
#ifdef POLYA
    polya_res /= (double)par->n_configs;
    printf("\nAverage:\t%f\n", polya_res);
#endif
#ifndef TADP
#ifndef POLYA
    printf("\nAverages:\n");
    gsl_matrix_view m_temp_result_view;
    m_result_view = gsl_matrix_view_array(results, par->L_t - 1, par->L - 1);
    for (int i = 1; i < par->n_configs; i++) 
        m_temp_result_view = gsl_matrix_view_array(results + 4 * i, par->L_t - 1, par->L - 1);
        gsl_matrix_add(
            &m_result_view.matrix, 
            &m_temp_result_view.matrix
        );
    gsl_matrix_scale(
        &m_result_view.matrix, 
        1. / (double)(par->n_configs)
    );
    for (int k = 0, index = 0; k < par->L - 1; k++) {
        for (int l = 0; l < par->L_t - 1; l++) {
            printf(
                "%7.3e\t", 
                gsl_matrix_get(&m_result_view.matrix, l, k)
            );
            index++;
            if (index == 15) 
                break;
        }
        if (index == 15) 
            break;
    }
    printf("\n");
#endif
#endif

    acceptance /= 
        (double)par->n_configs * 
        (double)par->n_corr * 
        (double)par->L_t * 
        pow((double)par->L, 3.) * 
        4. * 
        (double)par->n_hits;
    printf("Acceptance: %3.2f\n", acceptance);
    
    free(su2);
#ifndef TADP
#ifndef POLYA
    free(results);
#endif
#endif
    return 0;
}

/* calculate the position of a specific link in the array */
static inline int ind(int i, int j, int k, int l, int dir, int le) {
    return (((i * le + j) * le + k) * le + l) * 4 + dir;
}

/* calculate the position of a specific link in the array, while applying periodic boundary conditions to all
 * coordinates */
static inline int periodic_ind(int i, int j, int k, int l, int dir, int le_0, int le) {
    return (((((i + le_0) % le_0) * 
        le + ((j + le) % le)) * 
        le + ((k + le) % le)) * 
        le + ((l + le) % le)) * 
        4 + dir;
}

/* calculate the position of a specific link in the temp_lat array */
static inline int temp_ind(int i, int j, int k, int l, int dir, int ll, int le_max, int le) {
    return ((((i * le + j) * le + k) * le + l) * 4 + dir) * le_max + ll;
}

/* calculate the position of a specific link in the temp_lat array, 
 * while applying periodic boundary conditions to all coordinates */
static inline int temp_periodic_ind(
    int i, 
    int j, 
    int k, 
    int l, 
    int dir, 
    int ll, 
    int le_max, 
    int le_0, 
    int le
) {
    return ((((((i + le_0) % le_0) * 
        le + ((j + le) % le)) * 
        le + ((k + le) % le)) * 
        le + ((l + le) % le)) * 
        4 + dir) *
        le_max + ll;
}

void psl_su2_memcpy(double *dest, const double *src) {
    for (int i = 0; i < 4; i++) 
        dest[i] = src[i];
}

void psl_su2_dagger(double *m) {
    for (int i = 1; i < 4; i++) 
        m[i] = -m[i];
}

void psl_su2_dagger_memcpy(double *dest, const double *src) {
    dest[0] = src[0];
    for (int i = 1; i < 4; i++) 
        dest[i] = -src[i];
}
/*
void psl_su2_product_test(int dag_1, int dag_2, double *m_1, double *m_2, double *m_3) {
    int signum_1 = (dag_1 ? -1 : 1);
    int signum_2 = (dag_2 ? -1 : 1);
    int signum = signum_1 * signum_2;
    m_3[0] = m_1[0] * m_2[0] - signum * (m_1[1] * m_2[1] + m_1[2] * m_2[2] + m_1[3] * m_2[3]);
    m_3[1] = signum_2 * m_1[0] * m_2[1] + signum_1 * m_1[1] * m_2[0] + signum * (m_1[2] * m_2[3] - m_1[3] * m_2[2]);
    m_3[2] = signum_2 * m_1[0] * m_2[2] + signum_1 * m_1[2] * m_2[0] + signum * (m_1[3] * m_2[1] - m_1[1] * m_2[3]);
    m_3[3] = signum_2 * m_1[0] * m_2[3] + signum_1 * m_1[3] * m_2[0] + signum * (m_1[1] * m_2[2] - m_1[2] * m_2[1]);
}
*/

void psl_su2_product_notrans_notrans(const double *m_1, const double *m_2, double *m_3) {
    m_3[0] = m_1[0] * m_2[0] - m_1[1] * m_2[1] - m_1[2] * m_2[2] - m_1[3] * m_2[3];
    m_3[1] = m_1[0] * m_2[1] + m_1[1] * m_2[0] + m_1[2] * m_2[3] - m_1[3] * m_2[2];
    m_3[2] = m_1[0] * m_2[2] + m_1[2] * m_2[0] + m_1[3] * m_2[1] - m_1[1] * m_2[3];
    m_3[3] = m_1[0] * m_2[3] + m_1[3] * m_2[0] + m_1[1] * m_2[2] - m_1[2] * m_2[1];
}

void psl_su2_product_conjtrans_conjtrans(const double *m_1, const double *m_2, double *m_3) {
    m_3[0] = m_1[0] * m_2[0] - m_1[1] * m_2[1] - m_1[2] * m_2[2] - m_1[3] * m_2[3];
    m_3[1] = -m_1[0] * m_2[1] - m_1[1] * m_2[0] + m_1[2] * m_2[3] - m_1[3] * m_2[2];
    m_3[2] = -m_1[0] * m_2[2] - m_1[2] * m_2[0] + m_1[3] * m_2[1] - m_1[1] * m_2[3];
    m_3[3] = -m_1[0] * m_2[3] - m_1[3] * m_2[0] + m_1[1] * m_2[2] - m_1[2] * m_2[1];
}

void psl_su2_product_notrans_conjtrans(const double *m_1, const double *m_2, double *m_3) {
    m_3[0] = m_1[0] * m_2[0] + m_1[1] * m_2[1] + m_1[2] * m_2[2] + m_1[3] * m_2[3];
    m_3[1] = -m_1[0] * m_2[1] + m_1[1] * m_2[0] - m_1[2] * m_2[3] + m_1[3] * m_2[2];
    m_3[2] = -m_1[0] * m_2[2] + m_1[2] * m_2[0] - m_1[3] * m_2[1] + m_1[1] * m_2[3];
    m_3[3] = -m_1[0] * m_2[3] + m_1[3] * m_2[0] - m_1[1] * m_2[2] + m_1[2] * m_2[1];
}

void psl_su2_product_conjtrans_notrans(const double *m_1, const double *m_2, double *m_3) {
    m_3[0] = m_1[0] * m_2[0] + m_1[1] * m_2[1] + m_1[2] * m_2[2] + m_1[3] * m_2[3];
    m_3[1] = m_1[0] * m_2[1] - m_1[1] * m_2[0] - m_1[2] * m_2[3] + m_1[3] * m_2[2];
    m_3[2] = m_1[0] * m_2[2] - m_1[2] * m_2[0] - m_1[3] * m_2[1] + m_1[1] * m_2[3];
    m_3[3] = m_1[0] * m_2[3] - m_1[3] * m_2[0] - m_1[1] * m_2[2] + m_1[2] * m_2[1];
}

/* this function calculates the matrix-product of three SU(2)-matrices and -->ADDs<-- the result onto the m_sum
 * matrix */
void psl_su2_product_3_add(
    int dag_1,
    int dag_2, 
    int dag_3, 
    const double *m_1, 
    const double *m_2, 
    const double *m_3, 
    double *m_sum
) {
    double m_temp[8] = {0};

#ifdef DEBUG
    if (((dag_1 != 0) && (dag_1 != 1)) || ((dag_2 != 0) && (dag_2 != 1)) || ((dag_3 != 0) && (dag_3 != 1))) 
        printf("Warning: dagger-flag not set to zero or one in psl_su2_product_3_add. Continuing...\n");
#endif

    if ((dag_1 == 0) && (dag_2 == 0))
        psl_su2_product_notrans_notrans(m_1, m_2, m_temp);
    else if ((dag_1 == 1) && (dag_2 == 1))
        psl_su2_product_conjtrans_conjtrans(m_1, m_2, m_temp);
    else if ((dag_1 == 0) && (dag_2 == 1))
        psl_su2_product_notrans_conjtrans(m_1, m_2, m_temp);
    else 
        psl_su2_product_conjtrans_notrans(m_1, m_2, m_temp);

    if (dag_3 == 0) 
        psl_su2_product_notrans_notrans(m_temp, m_3, m_temp + 4);
    else 
        psl_su2_product_notrans_conjtrans(m_temp, m_3, m_temp + 4);

    for (int i = 0; i < 4; i++) 
        m_sum[i] += m_temp[4 + i];
}

/* calculates the product of 4 SU(2)-matrices */
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
) {
    double m_temp[4] = {0};

#ifdef DEBUG
    if (((dag_1 != 0) && (dag_1 != 1)) || 
        ((dag_2 != 0) && (dag_2 != 1)) || 
        ((dag_3 != 0) && (dag_3 != 1)) || 
        ((dag_4 != 0) && (dag_4 != 1))
    ) 
        printf("Warning: dagger-flag not set to zero or one in psl_su2_product_3_add. Continuing...\n");
#endif

    if ((dag_1 == 0) && (dag_2 == 0))
        psl_su2_product_notrans_notrans(m_1, m_2, m_result);
    else if ((dag_1 == 1) && (dag_2 == 1))
        psl_su2_product_conjtrans_conjtrans(m_1, m_2, m_result);
    else if ((dag_1 == 0) && (dag_2 == 1))
        psl_su2_product_notrans_conjtrans(m_1, m_2, m_result);
    else 
        psl_su2_product_conjtrans_notrans(m_1, m_2, m_result);

    if (dag_3 == 0) 
        psl_su2_product_notrans_notrans(m_result, m_3, m_temp);
    else 
        psl_su2_product_notrans_conjtrans(m_result, m_3, m_temp);
    
   if (dag_4 == 0) 
        psl_su2_product_notrans_notrans(m_temp, m_4, m_result);
    else 
        psl_su2_product_notrans_conjtrans(m_temp, m_4, m_result);
}


/* calculates a plaquette */
double plaquette_op(
    const PAR *par,
    const double *lattice, 
    int i_start, 
    int j_start, 
    int k_start, 
    int l_start, 
    int dir_1, 
    int dir_2
) {
    double m_temp[8] = {0};
    
    psl_su2_product_notrans_notrans(
        lattice + 4 * ind(i_start, j_start, k_start, l_start, dir_1, par->L), 
        lattice + 4 * periodic_ind(
            i_start + ((dir_1 == 0) ? 1 : 0), 
            j_start + ((dir_1 == 1) ? 1 : 0), 
            k_start + ((dir_1 == 2) ? 1 : 0), 
            l_start + ((dir_1 == 3) ? 1 : 0),
            dir_2, 
            par->L_t, 
            par->L
        ), 
        m_temp
    );
    psl_su2_product_notrans_conjtrans(
        m_temp, 
        lattice + 4 * periodic_ind(
            i_start + ((dir_2 == 0) ? 1 : 0), 
            j_start + ((dir_2 == 1) ? 1 : 0), 
            k_start + ((dir_2 == 2) ? 1 : 0), 
            l_start + ((dir_2 == 3) ? 1 : 0), 
            dir_1, 
            par->L_t, 
            par->L
        ), 
        m_temp + 4
    );
    psl_su2_product_notrans_conjtrans(
        m_temp + 4, 
        lattice + 4 * ind(
            i_start, 
            j_start, 
            k_start, 
            l_start, 
            dir_2, 
            par->L
        ), 
        m_temp
    );
    /* return half of real part of trace of m_temp */
    return m_temp[0];
}

/* calculates the product of 6 matrices 
void psl_matrix_complex_product_6(
    PAR *par,
    const gsl_matrix_complex *matrix_1, 
    const gsl_matrix_complex *matrix_2, 
    const gsl_matrix_complex *matrix_3, 
    const gsl_matrix_complex *matrix_4,
    const gsl_matrix_complex *matrix_5, 
    const gsl_matrix_complex *matrix_6,
    gsl_matrix_complex *m_result
) {
    const gsl_complex z_zero = gsl_complex_rect(0., 0.), 
                    z_one = gsl_complex_rect(1., 0.);

    psl_matrix_complex_product_4(
        par,
        matrix_1, 
        matrix_2, 
        matrix_3, 
        matrix_4, 
        m_result
    );
    gsl_blas_zgemm(
        CblasNoTrans, 
        CblasNoTrans, 
        z_one, 
        m_result, 
        matrix_5, 
        z_zero, 
        par->m_workspace
    );
    gsl_blas_zgemm(
        CblasNoTrans, 
        CblasNoTrans, 
        z_one, 
        par->m_workspace, 
        matrix_6, 
        z_zero, 
        m_result
    );
} */

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
) {
#ifdef DEBUG
    if ((dir_1 >= 4) || (dir_2 >= 4)) 
        printf("Warning: lattice_loop_rect was used with dir_1 or dir_2 >= 4. Something probably went wrong!\n");
#endif
    double m_temp[4] = {0};
    m_temp[0] = 1;  /* set identity */

    wilson_line_product(
        par, 
        lattice, 
        i_start, 
        j_start, 
        k_start, 
        l_start, 
        dir_1, 
        L_1, 
        m_temp
    );
    wilson_line_product(
        par, 
        lattice, 
        i_start + ((dir_1 == 0) ? L_1 : 0), 
        j_start + ((dir_1 == 1) ? L_1 : 0), 
        k_start + ((dir_1 == 2) ? L_1 : 0), 
        l_start + ((dir_1 == 3) ? L_1 : 0), 
        dir_2, 
        L_2, 
        m_temp
    );
    wilson_line_product(
        par, 
        lattice, 
        i_start + ((dir_2 == 0) ? L_2 : 0), 
        j_start + ((dir_2 == 1) ? L_2 : 0), 
        k_start + ((dir_2 == 2) ? L_2 : 0), 
        l_start + ((dir_2 == 3) ? L_2 : 0), 
        dir_1 + 4, 
        L_1, 
        m_temp
    );
    wilson_line_product(
        par, 
        lattice, 
        i_start, 
        j_start, 
        k_start, 
        l_start, 
        dir_2 + 4, 
        L_2, 
        m_temp
    );

    *result = m_temp[0];
}

/* calculate the product of a given matrix m_product with n_matrices links along the direction dir from a
 * starting point and save it under m_product */
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
) {
    int sign, 
        dir_abs, 
        dag;
    double m_temp[4] = {0};
#ifdef DEBUG
    if ((dir < 0) || (dir >= 8)) 
        printf("Warning: lattice_line_product was used with dir not in the right range. Something probably went wrong.\n");
#endif
    if (dir < 4) { 
        sign = 1;
        dir_abs = dir;
        dag = 0;
    } else {
        sign = -1;
        dir_abs = dir - 4;
        dag = 1;
    } 

#ifdef DEBUG
    if ((dag != 0) && (dag != 1)) 
        printf("Warning: conjugate transpose flag isn't set to 0 or 1 in lattice_line_product. Continuing...\n");
#endif
    if (dag == 0) 
        for (int x = 0; x < n_matrices; x++) {
            psl_su2_product_notrans_notrans(
                m_product, 
                lattice + 4 * periodic_ind(
                    i_start + ((dir_abs == 0) ? sign * x - dag : 0), 
                    j_start + ((dir_abs == 1) ? sign * x - dag : 0), 
                    k_start + ((dir_abs == 2) ? sign * x - dag : 0), 
                    l_start + ((dir_abs == 3) ? sign * x - dag : 0), 
                    dir_abs, 
                    par->L_t, 
                    par->L
                ), 
                m_temp
            );
            psl_su2_memcpy(m_product, m_temp);
        }
    else 
        for (int x = 0; x < n_matrices; x++) {
            psl_su2_product_notrans_conjtrans(
                m_product, 
                lattice + 4 * periodic_ind(
                    i_start + ((dir_abs == 0) ? sign * x - dag : 0), 
                    j_start + ((dir_abs == 1) ? sign * x - dag : 0), 
                    k_start + ((dir_abs == 2) ? sign * x - dag : 0), 
                    l_start + ((dir_abs == 3) ? sign * x - dag : 0), 
                    dir_abs, 
                    par->L_t, 
                    par->L
                ), 
                m_temp
            );
            psl_su2_memcpy(m_product, m_temp);
        }
}

/* calculate the product of a given matrix m_product with n_matrices links along the direction dir from a
 * starting point and save it under m_product. --> Only works in Wilson-measurement loop! <-- */
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
) {
    int dir_abs,
        dag, 
        temp_index,
        last_temp_index,
        L_max;
    double m_temp[4] = {0};
#ifdef DEBUG
    if ((dir < 0) || (dir >= 8)) 
        printf("Warning: lattice_line_product was used with dir not in the right range. Something probably went wrong.\n");
#endif
    /* get some more practical parameters */
    if (dir < 4) { 
        dir_abs = dir;
        dag = 0;
    } else {
        dir_abs = dir - 4;
        dag = 1;
    }
    L_max = (par->L_t > par->L) ? par->L_t : par->L;

    /* generate index of result matrix (not daggered) */
    temp_index = temp_periodic_ind(
        i_start, 
        j_start, 
        k_start, 
        l_start, 
        dir_abs, 
        n_matrices - 1, 
        L_max, 
        par->L_t, 
        par->L
    ); 

    /* Look up if this matrix was calculated already in this measurement */
    if (!(par->temp_lat_filled[temp_index])) {
        if (n_matrices == 1) {  /* temp_lattice isn't filled at all at this starting point and direction */
            /* copy first (n = 1) matrix from lattice */
            for (int i = 0; i < 4; i++) 
                par->temp_lat[4 * temp_index + i] = lattice[4 * periodic_ind(
                    i_start, 
                    j_start, 
                    k_start, 
                    l_start, 
                    dir_abs, 
                    par->L_t, 
                    par->L
                ) + i];

            par->temp_lat_filled[temp_index] = 1;
        } else { 
            /* generate the index of the next shorter line product at the same position and direction  */
            last_temp_index = temp_periodic_ind(
                i_start, 
                j_start, 
                k_start, 
                l_start, 
                dir_abs, 
                n_matrices - 2, 
                L_max, 
                par->L_t,
                par->L
            );
#ifdef DEBUG
            if (par->temp_lat_filled[last_temp_index] == 0) 
                printf("Warning: lattice_line_product used a matrix of temp_lat that wasn't filled already in this measurement.\n");
#endif
            /* now multiply the missing matrix at the end of the line product onto the old line product and 
             * save it in the temp_lattice */
            psl_su2_product_notrans_notrans(
                par->temp_lat + 4 * last_temp_index, 
                lattice + 4 * periodic_ind(
                    i_start + ((dir_abs == 0) ? (n_matrices - 1) : 0), 
                    j_start + ((dir_abs == 1) ? (n_matrices - 1) : 0), 
                    k_start + ((dir_abs == 2) ? (n_matrices - 1) : 0), 
                    l_start + ((dir_abs == 3) ? (n_matrices - 1) : 0), 
                    dir_abs, 
                    par->L_t, 
                    par->L
                ), 
                par->temp_lat + 4 * temp_index
            );
            par->temp_lat_filled[temp_index] = 1;
        }
    }
    /* multiply m_product with the new line product and save it in m_product */
#ifdef DEBUG
        if ((dag != 0) && (dag != 1)) 
            printf("Warning: conjugate transpose flag isn't set to 0 or 1 in wilson_line_product. Continuing...\n");
#endif

    if (dag == 0) 
        psl_su2_product_notrans_notrans(
            m_product, 
            par->temp_lat + 4 * temp_index, 
            m_temp
        );
    else 
        psl_su2_product_notrans_conjtrans(
            m_product, 
            par->temp_lat + 4 * temp_index, 
            m_temp
        );
    for (int i = 0; i < 4; i++) 
        m_product[i] = m_temp[i];
}


void psl_su2_unitarize(double *matrix) {
    double norm = 0;
    for (int i = 0; i < 4; i++) 
        norm += matrix[i] * matrix[i];
    norm = sqrt(norm);
    for (int i = 0; i < 4; i++) 
        matrix[i] /= norm;
}

int measure_polyakov(const PAR *par, double *result, const double *lattice, const char *file_name) {
    FILE *data_file;
    char polyakov_file_name[1000];
    double m_temp[4] = {0};

    *result = 0.;

    strcpy(polyakov_file_name, file_name);
    polyakov_file_name[strlen(polyakov_file_name) - 4] = '\0';
    strcat(polyakov_file_name, "_polyakov.bin");

    data_file = fopen(polyakov_file_name, "ab");
    if (data_file == NULL) {
        printf("Error: Failed opening file to write Polyakov-loop data. Exiting...\n");
        return 1;
    }

    for (int j = 0; j < par->L; j++) {
        for (int k = 0; k < par->L; k++) {
            for (int l = 0; l < par->L; l++) {
                /* set m_temp to be a identity matrix */
                m_temp[0] = 1;
                for (int i = 1; i < 4; i++)
                    m_temp[i] = 0;
                lattice_line_product(
                    par, 
                    lattice, 
                    0, 
                    j, 
                    k, 
                    l, 
                    0, 
                    par->L_t, 
                    m_temp
                );
                /* take half of the real part of the trace of m_temp */
                *result += m_temp[0];
            }
        }
    }

    /* normalize result and apply tadpole correction */
    *result  /= (double)(par->L * par->L * par->L) * gsl_pow_int(par->tadpole, par->L_t);

    /* write result to file */
    if (fwrite(result, sizeof(double), 1, data_file) != 1) {
        printf("Error: Failed writing Polyakov-loop data to file. Exiting...\n");
        fclose(data_file);
        return 1;
    }

    fclose(data_file);

    return 0;
}

int measure(const PAR *par, gsl_matrix *results, const double *lattice, const char *file_name) {
    double result;
    int L_max; 
    unsigned long temp_lat_size;
    FILE *data_file;

    if (par->L_t > par->L) 
        L_max = par->L_t;
    else 
        L_max = par->L;

    temp_lat_size = (unsigned long)(par->L_t * par->L * par->L * par->L) * (unsigned long)(4 * L_max);
    
    data_file = fopen(file_name, "ab");
    if (data_file == NULL) {
        printf("Error: Failed opening data-file to write results of Wilson-loops. Exiting...\n");
        return 1;
    }

    gsl_matrix_set_zero(results);

    for (unsigned long i = 0; i < temp_lat_size; i++) 
        par->temp_lat_filled[i] = 0;

    for (int L_0 = 1; L_0 <= par->L_t - 1; L_0++) {
        for (int L_i = 1; L_i <= par->L - 1; L_i++) {
            
            for (int dir = 1; dir < 4; dir++) {

                for (int i = 0; i < par->L_t; i++) {
                    for (int j = 0; j < par->L; j++) {
                        for (int k = 0; k < par->L; k++) {
                            for (int l = 0; l < par->L; l++) {

                                lattice_loop_rect(
                                    par, 
                                    lattice, 
                                    i, 
                                    j, 
                                    k, 
                                    l, 
                                    0, 
                                    dir, 
                                    L_0, 
                                    L_i, 
                                    &result
                                );
                                gsl_matrix_set(
                                    results, 
                                    L_0 - 1, 
                                    L_i - 1, 
                                    gsl_matrix_get(results, L_0 - 1, L_i - 1) + result
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    /* normalize results and apply tadpole correction */
    for (int i = 0; i < par->L_t - 1; i++) {
        for (int j = 0; j < par->L - 1; j++) {
            gsl_matrix_set(
                results, 
                i, 
                j, 
                gsl_matrix_get(results, i, j) / (double)(par->L_t * par->L * par->L * par->L * 3) / gsl_pow_int(par->tadpole, 2 * (i + j + 2))
            );
        }
    }

    /* write results to file */
    for (int i = 0; i < par->L_t - 1; i++) {
        for (int j = 0; j < par->L - 1; j++) {
            result = gsl_matrix_get(results, i, j);
            if (fwrite(&result, sizeof(double), 1, data_file) != 1) {
                printf("Error: Failed writing results of Wilson-loops to data file. Exiting...\n");
                fclose(data_file);
                return 1;
            }
        }
    }

    fclose(data_file);

    return 0;
}

void measure_action_r(const PAR *par, const double *lattice, double *action) {
    double m_temp[4] = {0}; 

    *action = 0.;

    for (int i = 0; i < par->L_t; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {

                    for (int dir_1 = 0; dir_1 < 4; dir_1++) {
                        for (int dir_2 = dir_1 + 1; dir_2 < 4; dir_2++) {

                            /* calculate a x a Wilson-loop */
                            psl_su2_product_4(
                                0, 
                                0, 
                                1, 
                                1, 
                                lattice + 4 * ind(i, j, k, l, dir_1, par->L),
                                lattice + 4 * periodic_ind(
                                    i + (dir_1 == 0), 
                                    j + (dir_1 == 1), 
                                    k + (dir_1 == 2), 
                                    l + (dir_1 == 3), 
                                    dir_2,
                                    par->L_t, 
                                    par->L
                                ),
                                lattice + 4 * periodic_ind(
                                    i + (dir_2 == 0), 
                                    j + (dir_2 == 1), 
                                    k + (dir_2 == 2), 
                                    l + (dir_2 == 3), 
                                    dir_1,
                                    par->L_t, 
                                    par->L
                                ),
                                lattice + 4 * ind(i, j, k, l, dir_2, par->L),
                                m_temp
                            );
                            /* take half of the real part of the trace of m_temp */
                            *action += m_temp[0];
                        }
                    }
                }
            }
        }
    }
    /* factor -beta from Wilson action */
    *action *= -par->beta / (double)(par->tadpole * par->tadpole * par->tadpole * par->tadpole); 
}

void measure_action_l(const PAR *par, const double *lattice, double *action) {
    double m_temp[4] = {0}; 

    *action = 0.;

    for (int i = 0; i < par->L_t; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {

                    for (int dir_2 = 0; dir_2 < 4; dir_2++) {
                        for (int dir_1 = dir_2 + 1; dir_1 < 4; dir_1++) {

                            /* calculate a x a Wilson-loop */
                            psl_su2_product_4(
                                0, 
                                0, 
                                1, 
                                1, 
                                lattice + 4 * ind(i, j, k, l, dir_1, par->L),
                                lattice + 4 * periodic_ind(
                                    i + (dir_1 == 0), 
                                    j + (dir_1 == 1), 
                                    k + (dir_1 == 2), 
                                    l + (dir_1 == 3), 
                                    dir_2,
                                    par->L_t, 
                                    par->L
                                ),
                                lattice + 4 * periodic_ind(
                                    i + (dir_2 == 0), 
                                    j + (dir_2 == 1), 
                                    k + (dir_2 == 2), 
                                    l + (dir_2 == 3), 
                                    dir_1,
                                    par->L_t, 
                                    par->L
                                ),
                                lattice + 4 * ind(i, j, k, l, dir_2, par->L),
                                m_temp
                            );
                            /* take half of the real part of the trace of m_temp */
                            *action += m_temp[0];
                        }
                    }
                }
            }
        }
    }
    /* factor -beta from Wilson action */
    *action *= -par->beta / (double)(par->tadpole * par->tadpole * par->tadpole * par->tadpole); 
}

void unitarize_lattice(const PAR *par, double *lattice) {

    for (int i = 0; i < par->L_t; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {

                    for (int dir = 0; dir < 4; dir++) {
                        psl_su2_unitarize(lattice + 4 * ind(i, j, k, l, dir, par->L));
                    }
                }
            }
        }
    }
}

void update_lattice(const PAR *par, double *lattice, const double *su2, double *acceptance) {
    double Delta_S;
    double m_workspace_arr[16] = {0}; 
    const int
        latt_proposal_index_offset = 0,
        latt_diff_index_offset = 4, 
        loops_index_offset = 8,
        staples_index_offset =  12;

    //for (int i = 0; i < par->L * par->L * par->L * par->L_t * 4 * 4; i++) 
    //    printf("%d\t%g\n", i % 4, lattice[i]);

    /* Loop through the whole lattice (all independent links) */
    for (int i = 0; i < par->L_t; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {

                    for (int dir_1 = 0; dir_1 < 4; dir_1++) {
                        /* calculate the sum over all staples 
                         * (parts of plaquette operators which contain the link that is to be updated) 
                         * which is needed to calculate the differences in action later */
                        for (int ii = 0; ii < 4; ii++) 
                            m_workspace_arr[staples_index_offset + ii] = 0;
                        
                        for (int dir_2 = 0; dir_2 < 4; dir_2++) {
                            
                            if (dir_1 == dir_2) 
                                continue;
                            
                            psl_su2_product_3_add(
                                0, 
                                1, 
                                1, 
                                lattice + 4 * periodic_ind(
                                    i + ((dir_1 == 0) ? 1 : 0), 
                                    j + ((dir_1 == 1) ? 1 : 0),
                                    k + ((dir_1 == 2) ? 1 : 0),
                                    l + ((dir_1 == 3) ? 1 : 0),
                                    dir_2,
                                    par->L_t, 
                                    par->L
                                ),
                                lattice + 4 * periodic_ind(
                                    i + ((dir_2 == 0) ? 1 : 0),
                                    j + ((dir_2 == 1) ? 1 : 0),
                                    k + ((dir_2 == 2) ? 1 : 0),
                                    l + ((dir_2 == 3) ? 1 : 0),
                                    dir_1,
                                    par->L_t, 
                                    par->L
                                ),
                                lattice + 4 * ind(i, j, k, l, dir_2, par->L),
                                m_workspace_arr + staples_index_offset
                            );

                            psl_su2_product_3_add(
                                1, 
                                1, 
                                0, 
                                lattice + 4 * periodic_ind(
                                    i + ((dir_1 == 0) ? 1 : 0) + ((dir_2 == 0) ? -1 : 0), 
                                    j + ((dir_1 == 1) ? 1 : 0) + ((dir_2 == 1) ? -1 : 0),
                                    k + ((dir_1 == 2) ? 1 : 0) + ((dir_2 == 2) ? -1 : 0),
                                    l + ((dir_1 == 3) ? 1 : 0) + ((dir_2 == 3) ? -1 : 0),
                                    dir_2,
                                    par->L_t, 
                                    par->L
                                ),
                                lattice + 4 * periodic_ind(
                                    i + ((dir_2 == 0) ? -1 : 0),
                                    j + ((dir_2 == 1) ? -1 : 0),
                                    k + ((dir_2 == 2) ? -1 : 0),
                                    l + ((dir_2 == 3) ? -1 : 0),
                                    dir_1,
                                    par->L_t, 
                                    par->L
                                ),
                                lattice + 4 * periodic_ind(
                                    i + ((dir_2 == 0) ? -1 : 0), 
                                    j + ((dir_2 == 1) ? -1 : 0),
                                    k + ((dir_2 == 2) ? -1 : 0),
                                    l + ((dir_2 == 3) ? -1 : 0),
                                    dir_2,
                                    par->L_t, 
                                    par->L
                                ),
                                m_workspace_arr + staples_index_offset
                            );
                            //for (int jj = 0; jj < 4; jj++) 
                            //    printf("DEBUG: %d\t%g\n", jj, m_workspace_arr[staples_index_offset + jj]);
                        }
                        /* do the specified number of MC-updates of this link */
                        for (int n = 0; n < par->n_hits; n++) {
                            /* multiply with a random SU(2)-matrix and save in a temporary matrix, until we
                             * know if we will use it */
                            psl_su2_product_notrans_notrans(
                                su2 + 4 * (int)(par->n_su2 * gsl_rng_uniform(par->ran_gen)), 
                                lattice + 4 * ind(i, j, k, l, dir_1, par->L),
                                m_workspace_arr + latt_proposal_index_offset
                            );
                            
                            /* measure the difference in action the new link operator does for the lattice */
                            for (int ii = 0; ii < 4; ii++) 
                                m_workspace_arr[latt_diff_index_offset + ii] = 
                                    m_workspace_arr[latt_proposal_index_offset + ii] - 
                                    lattice[4 * ind(i, j, k, l, dir_1, par->L) + ii];

                            psl_su2_product_notrans_notrans(
                                m_workspace_arr + latt_diff_index_offset,
                                m_workspace_arr + staples_index_offset,
                                m_workspace_arr + loops_index_offset
                            );

                            Delta_S = -par->beta * 
                                m_workspace_arr[loops_index_offset] / 
                                (double)(par->tadpole * par->tadpole * par->tadpole * par->tadpole);
                            
                            /* accept/reject step */
                            if (Delta_S <= 0.) {
                                /* save updated matrix */
                                for (int ii = 0; ii < 4; ii++) 
                                    lattice[4 * ind(i, j, k, l, dir_1, par->L) + ii] = 
                                        m_workspace_arr[latt_proposal_index_offset + ii];

                                *acceptance += 1.;
                                                                
                            } else {
                                if (gsl_rng_uniform(par->ran_gen) < exp(-Delta_S)) {
                                    /* save updated matrix */
                                    for (int ii = 0; ii < 4; ii++) 
                                        lattice[4 * ind(i, j, k, l, dir_1, par->L) + ii] = 
                                            m_workspace_arr[latt_proposal_index_offset + ii];
                                    
                                    *acceptance += 1.;
                                                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void init_lattice(const PAR *par, double *lattice, const double *su2) {
    int su2_index;
    /* initialize all independent links with random SU(3) matrices */
    for (int i = 0; i < par->L_t; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {

                    for (int dir = 0; dir < 4; dir++) {

                        su2_index = (int)(gsl_rng_uniform(par->ran_gen) * par->n_su2);
                        
                        for (int ii = 0; ii < 4; ii++) 
                            lattice[4 * ind(i, j, k, l, dir, par->L) + ii] = 
                            su2[4 * su2_index + ii];
                    }
                }
            }
        }
    }
}

