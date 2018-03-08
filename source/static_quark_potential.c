#include "static_quark_potential.h"


static inline int ind(int i, int j, int k, int l, int dir, int le);
static inline int periodic_ind(int i, int j, int k, int l, int dir, int le);


#ifndef MAIN_STATIC__C
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
    par->seed = 0;
    par->n_configs = 10;
    par->gen_type = (gsl_rng_type *)gsl_rng_ranlxs0;
    par->n_su2 = 100;
    par->eps = 1.;
    par->n_hits = 10;
    par->n_therm = 500;
    par->n_corr = 50;
    par->beta = 0.;
    par->tadpole = 0.;

    par->m_workspace = gsl_matrix_complex_calloc(2, 2);
    if (par->m_workspace == NULL) {
        printf("Error: Failed allocating memory for matrix-workspace. Exiting...\n");
        free(par);
        exit(EXIT_FAILURE);
    }
  
    if (argc == 1) {
        printf("Usage: %s L=16 nconfigs=100 run\n", argv[0]);
        printf("Optional arguments (with defaults) L=%d seed=%ld", par->L, par->seed);
        gsl_matrix_complex_free(par->m_workspace);
        free(par);
        exit(EXIT_SUCCESS);
    }
    
    /* read_args interprets the arguments given to the program and starts it when "run" appears */
    for (int iarg = 1; iarg < argc; iarg++)
        if (read_args(par, argv[iarg])) {
            gsl_matrix_complex_free(par->m_workspace);
            free(par);
            exit(EXIT_FAILURE);
        }

    gsl_matrix_complex_free(par->m_workspace);
    free(par);
    exit(EXIT_SUCCESS);
}
#endif

int read_args(PAR *par, char *arg) {
    static gsl_matrix_complex **lattice = NULL;
    int success = 1;
    char *s;
    
    if (!strcmp(arg, "run")) {
        if (!par->L) {
            printf("Give system size L!\n");
            return 1;
        }
        if (par->beta == 0.) {
            printf("Give beta!\n");
            return 1;
        }
        if (par->eps == 1.) {
            printf("Warning: Standard value of epsilon is 1.\n");
        }
        /* preparing for the simulation */
        int num_of_links = par->L * par->L * par->L * par->L * 8;
        
        /* initialization of random number generator */
        par->ran_gen = gsl_rng_alloc(par->gen_type);
        if (par->ran_gen == NULL) {
            if (lattice) {
                for (int i = 0; i < num_of_links; i++)
                    gsl_matrix_complex_free(lattice[i]);
                free(lattice);
            }
            printf("Error: Failed allocating memory for the random number generator! Exiting...\n");
            return 1;
        }
        if (!(par->seed)) 
            par->seed = (long)time(NULL);
        gsl_rng_set(par->ran_gen, par->seed);
        
        /* running the simulation with already given parameters */
        success = !simulate(par, lattice);
        
        gsl_rng_free(par->ran_gen);
        
        for (int i = 0; i < num_of_links; i++)
            gsl_matrix_complex_free(lattice[i]);
        free(lattice);

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
        int num_of_links, tempL;
   
        /* checking if the given lattice size is a power of 2 */
        tempL = strtod(s, NULL);
        if ((tempL > 0) && !(tempL & (tempL - 1))) par->L = tempL;
        else {
            fprintf(stderr, "L has to be a power of 2. Exiting...\n");
            return 1;
        }
        
        /* allocate memory for the lattice (link matrices) */
        num_of_links = par->L * par->L * par->L * par->L * 8;
        
        lattice = realloc(lattice, num_of_links * sizeof(gsl_matrix_complex *));
        if (lattice == NULL) {
            printf("Error: Failed allocating memory for the lattice! Exiting...\n");
            return 1;
        }
        for (int i = 0; i < num_of_links; i++) {
            lattice[i] = gsl_matrix_complex_alloc(2, 2);
            if (lattice[i] == NULL) {
                printf("Error: Failed allocating memory for the lattice-matrices! Exiting...\n");
                for (int j = 0; j < i; j++)
                    gsl_matrix_complex_free(lattice[i]);
                free(lattice);
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

int simulate(PAR *par, gsl_matrix_complex **lattice) {
    double acceptance = 0.;
    gsl_matrix_complex **su2;
    gsl_matrix **results;
    
    /* allocate memory for a specified number of random SU(2) matrices */
    su2 = malloc(par->n_su2 * sizeof(gsl_matrix_complex *));
    if (su2 == NULL) {
        printf("Error: Failed allocating memory for array of SU(2)-matrix addresses. Exiting...\n");
        return 1;
    }
    for (int i = 0; i < par->n_su2; i++) {
        su2[i] = gsl_matrix_complex_calloc(2, 2);
        if (su2[i] == NULL) {
            printf("Error: Failed allocating memory for SU(2)-matrix. Exiting...\n");
            for (int j = 0; j < i; j++) 
                gsl_matrix_complex_free(su2[j]);
            free(su2);
            return 1;
        }
    }
    /* allocate memory for a specified number of results */
    results = malloc(par->n_configs * sizeof(gsl_matrix *));
    if (results == NULL) {
        printf("Error: Failed allocating memory for array of results. Exiting...\n");
        for (int i = 0; i < par->n_su2; i++) 
            gsl_matrix_complex_free(su2[i]);
        free(su2);
        return 1;
    }
    for (int i = 0; i < par->n_configs; i++) {
        results[i] = gsl_matrix_alloc(par->L - 1, par->L - 1);
        if (results[i] == NULL) {
            printf("Error: Failed allocating memory for results matrix. Exiting...\n");
            for (int j = 0; j < par->n_su2; j++) 
                gsl_matrix_complex_free(su2[j]);
            free(su2);
            for (int j = 0; j < i; j++) 
                gsl_matrix_free(results[j]);
            free(results);
            return 1;
        }
    }

    
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
    if (init_su2(par, su2)) {
        for (int i = 0; i < par->n_su2; i++)
            gsl_matrix_complex_free(su2[i]);
        free(su2);
        for (int i = 0; i < par->n_configs; i++) 
            gsl_matrix_free(results[i]);
        free(results);
        return 1;
    }

    /* DEBUG: 
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
    printf("Initializing %d x %d x %d x %d lattice...\n", par->L, par->L, par->L, par->L);
    init_lattice(par, lattice, su2);
    
    /* thermalize the lattice */
    printf("Thermalizing (%d sweeps)...\n", par->n_therm);
    for (int i = 0; i < par->n_therm; i++) {
        if(update_lattice(par, lattice, su2, &acceptance)) {
            for (int j = 0; j < par->n_su2; j++) 
                gsl_matrix_complex_free(su2[j]);
            free(su2);
            for (int j = 0; j < par->n_configs; j++) 
                gsl_matrix_free(results[j]);
            free(results);
            return 1;
        }
        if (unitarize_lattice(par, lattice)) {
            for (int j = 0; j < par->n_su2; j++) 
                gsl_matrix_complex_free(su2[j]);
            free(su2);
            for (int j = 0; j < par->n_configs; j++) 
                gsl_matrix_free(results[j]);
            free(results);
            return 1;
        }
    }
    acceptance /= (double)(par->n_therm * par->L * par->L * par->L * par->L * 4 * par->n_hits);
    printf("Thermalization-acceptance: %3.2f\n", acceptance);
    acceptance = 0.;

    /* take measurements after every n_corr-th update-sweep */
    printf(
        "Taking measurements of %d configurations(n_corr=%d, n_hits=%d):\na x a\ta x 2a\n",
        par->n_configs, 
        par->n_corr, 
        par->n_hits
    );
    for (int i = 0; i < par->n_configs; i++) {
        for (int j = 0; j < par->n_corr; j++) {
            if(update_lattice(par, lattice, su2, &acceptance)) {
                for (int j = 0; j < par->n_su2; j++) 
                    gsl_matrix_complex_free(su2[j]);
                free(su2);
                for (int j = 0; j < par->n_configs; j++) 
                    gsl_matrix_free(results[j]);
                free(results);
                return 1;
            }
            if (unitarize_lattice(par, lattice)) {
                for (int j = 0; j < par->n_su2; j++) 
                    gsl_matrix_complex_free(su2[j]);
                free(su2);
                for (int j = 0; j < par->n_configs; j++) 
                    gsl_matrix_free(results[j]);
                free(results);
                return 1;
            }

            /* DEBUG: print out ReTr of some random link on the lattice / print out the matrix product of a
             * link and his conjugate transpose
            gsl_complex z_temp = gsl_complex_rect(0., 0.);
            gsl_matrix_complex *m_temp = NULL;
            m_temp = gsl_matrix_complex_alloc(3, 3);
            for (int k = 0; k < 3; k++) 
                z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(lattice[ind(4, 4, 4, 4, 0, par->L)], k, k));
            printf("DEBUG: ReTr = %7.2e\n", GSL_REAL(z_temp)); 
            gsl_blas_zgemm(
                CblasNoTrans, 
                CblasNoTrans, 
                gsl_complex_rect(1., 0.),
                lattice[ind(3, 3, 3, 3, 0, par->L)], 
                lattice[ind(4, 3, 3, 3, 4, par->L)], 
                gsl_complex_rect(0., 0.),
                m_temp
            );
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    z_temp = gsl_matrix_complex_get(m_temp, k, l);
                    printf("%3.2f+%3.2f\t", GSL_REAL(z_temp), GSL_IMAG(z_temp));
                }
                printf("\n");
            }
            gsl_matrix_complex_free(m_temp);
            */
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

        if (measure_aa_a2a(par, results, lattice)) {
            for (int j = 0; j < par->n_su2; j++) 
                gsl_matrix_complex_free(su2[j]);
            free(su2);
            for (int j = 0; j < par->n_configs; j++) 
                gsl_matrix_free(results[j]);
            free(results);
            return 1;
        }
        printf("%3.2f\t%3.2f\n", results[0], results[1]);
        /* measure tadpole and print out the results
        double tadpole_result;
        measure_tadpole(par, lattice, &tadpole_result);
        printf("u_0 = %3.2f\n", tadpole_result);
        */
        /* printf("%g\n", gauge_inv(par, lattice)); */

        /* DEBUG: gauge transform the lattice and look at the effect on measurements 
        double action;
        measure_action_r(par, lattice, &action);
        printf("action = %7.2e\n", action); 
        double loop;
        if (lattice_loop_rect(
            par, 
            lattice, 
            3, 
            3, 
            3, 
            3, 
            1, 
            2, 
            3, 
            4, 
            &loop
        )) {
            for (int j = 0; j < par->n_su2; j++) 
                gsl_matrix_complex_free(su2[j]);
            free(su2);
            return 1;
        }
        printf("DEBUG: loop = %7.2e\n", loop);
        if (gauge_transform_lattice(par, lattice)) {
            for (int j = 0; j < par->n_su2; j++) 
                gsl_matrix_complex_free(su2[j]);
            free(su2);
            return 1;
        }
        if (lattice_loop_rect(
            par, 
            lattice, 
            3, 
            3, 
            3, 
            3, 
            1, 
            2, 
            3, 
            4, 
            &loop
        )) {
            for (int j = 0; j < par->n_su2; j++) 
                gsl_matrix_complex_free(su2[j]);
            free(su2);
            return 1;
        }
        printf("DEBUG: loop = %7.2e\t<--after gauge transf\n", loop);
        

        if (measure_aa_a2a(par, results, lattice)) {
            for (int j = 0; j < par->n_su2; j++) 
                gsl_matrix_complex_free(su2[j]);
            free(su2);
            return 1;
        }
        printf("%3.2f\t%3.2f\t<--gauge transformed\n", results[0], results[1]);
        measure_action_r(par, lattice, &action);
        printf("action = %7.2e\n", action);
        DEBUG END */
    }
    acceptance /= (double)(par->n_configs * par->n_corr * par->L * par->L * par->L * par->L * 4 * par->n_hits);
    printf("Acceptance: %3.2f\n", acceptance);
    
    for (int i = 0; i < par->n_su2; i++)
        gsl_matrix_complex_free(su2[i]);
    free(su2);
    for (int i = 0; i < par->n_configs; i++) 
        gsl_matrix_free(results[i]);
    free(results);

    return 0;
}

/* calculate the position of a specific link in the array */
static inline int ind(int i, int j, int k, int l, int dir, int le) {
    return (((le * i + j) * le + k) * le + l) * 8 + dir;
}

/* calculate the position of a specific link in the array, while applying periodic boundary conditions to all
 * coordinates */
static inline int periodic_ind(int i, int j, int k, int l, int dir, int le) {
    return (((le * (i & (le - 1)) + (j & (le - 1))) * le + (k & (le - 1))) * le + (l & (le - 1))) * 8 + dir;
}

/* function to get the conjugate transpose matrix */
void psl_matrix_complex_dagger(gsl_matrix_complex *m) {
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++)
            gsl_matrix_complex_set(m, i, j, gsl_complex_conjugate(gsl_matrix_complex_get(m, j, i)));
    }
}

/* function to save the conjugate transpose matrix of src in dest */
void psl_matrix_complex_dagger_memcpy(gsl_matrix_complex *dest, gsl_matrix_complex *src) {
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++)
            gsl_matrix_complex_set(dest, i, j, gsl_complex_conjugate(gsl_matrix_complex_get(src, j, i)));
    }
}

/* void psl_vector_complex_cross(
        const gsl_vector_complex *vec_1, 
        const gsl_vector_complex *vec_2, 
        gsl_vector_complex *vec_result
) {
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3, 
            k = (i + 2) % 3;
        gsl_vector_complex_set(
            vec_result, 
            i, 
            gsl_complex_sub(
                gsl_complex_mul(
                    gsl_vector_complex_get(vec_1, j), 
                    gsl_vector_complex_get(vec_2, k)
                ), 
                gsl_complex_mul(
                    gsl_vector_complex_get(vec_1, k),
                    gsl_vector_complex_get(vec_2, j)
                )
            )
        );
    }
} */


/* this function calculates the matrix-product of three matrices and -->ADDs<-- the result onto the */
void psl_matrix_complex_product_3_add(
    PAR *par,
    const gsl_matrix_complex *matrix_1, 
    const gsl_matrix_complex *matrix_2, 
    const gsl_matrix_complex *matrix_3, 
    gsl_matrix_complex *m_sum
) {
    const gsl_complex z_zero = gsl_complex_rect(0., 0.), 
                    z_one = gsl_complex_rect(1., 0.);

    gsl_blas_zgemm(
        CblasNoTrans, 
        CblasNoTrans, 
        z_one, 
        matrix_1, 
        matrix_2, 
        z_zero, 
        par->m_workspace
    );
    gsl_blas_zgemm(
        CblasNoTrans, 
        CblasNoTrans, 
        z_one, 
        par->m_workspace, 
        matrix_3, 
        z_one, 
        m_sum
    );
}

/* calculates the product of 4 matrices */
void psl_matrix_complex_product_4(
    PAR *par,
    const gsl_matrix_complex *matrix_1, 
    const gsl_matrix_complex *matrix_2, 
    const gsl_matrix_complex *matrix_3, 
    const gsl_matrix_complex *matrix_4,
    gsl_matrix_complex *m_result
) {
    const gsl_complex z_zero = gsl_complex_rect(0., 0.), 
                    z_one = gsl_complex_rect(1., 0.);

    gsl_blas_zgemm(
        CblasNoTrans, 
        CblasNoTrans, 
        z_one, 
        matrix_1, 
        matrix_2, 
        z_zero, 
        m_result
    );
    gsl_blas_zgemm(
        CblasNoTrans, 
        CblasNoTrans, 
        z_one, 
        m_result, 
        matrix_3, 
        z_zero, 
        par->m_workspace
    );
    gsl_blas_zgemm(
        CblasNoTrans, 
        CblasNoTrans, 
        z_one, 
        par->m_workspace, 
        matrix_4, 
        z_zero, 
        m_result
    );
}

/* calculates the product of 6 matrices */
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
}

/* calculates a planar Wilson-loop */
int lattice_loop_rect(
    const PAR *par, 
    gsl_matrix_complex **lattice, 
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
    gsl_matrix_complex *m_temp = NULL;
    gsl_complex z_temp = gsl_complex_rect(0., 0.);

    m_temp = gsl_matrix_complex_alloc(2, 2);
    if (m_temp == NULL) {
        printf("Error: Failed allocating memory for temporary matrix used for calculating a rectangle loop. Exiting...\n");
        return 1;
    }

    gsl_matrix_complex_set_identity(m_temp);
    lattice_line_product(
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
    lattice_line_product(
        par, 
        lattice, 
        i_start + L_1 * (dir_1 == 0), 
        j_start + L_1 * (dir_1 == 1), 
        k_start + L_1 * (dir_1 == 2), 
        l_start + L_1 * (dir_1 == 3), 
        dir_2, 
        L_2, 
        m_temp
    );
    lattice_line_product(
        par, 
        lattice, 
        i_start + L_1 * (dir_1 == 0) + L_2 * (dir_2 == 0), 
        j_start + L_1 * (dir_1 == 1) + L_2 * (dir_2 == 1), 
        k_start + L_1 * (dir_1 == 2) + L_2 * (dir_2 == 2), 
        l_start + L_1 * (dir_1 == 3) + L_2 * (dir_2 == 3), 
        dir_1 + 4, 
        L_1, 
        m_temp
    );
    lattice_line_product(
        par, 
        lattice, 
        i_start + L_2 * (dir_2 == 0), 
        j_start + L_2 * (dir_2 == 1), 
        k_start + L_2 * (dir_2 == 2), 
        l_start + L_2 * (dir_2 == 3), 
        dir_2 + 4, 
        L_2, 
        m_temp
    );

    for (int m = 0; m < 2; m++) {
        z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp, m, m));
    }
    *result = GSL_REAL(z_temp) / 2.;

    gsl_matrix_complex_free(m_temp);

    return 0;
}

/* calculate the product of a given matrix m_product with n_matrices links along the direction dir from a
 * starting point and save it under m_product*/
void lattice_line_product(
    const PAR *par, 
    gsl_matrix_complex **lattice, 
    int i_start, 
    int j_start, 
    int k_start, 
    int l_start, 
    int dir, 
    int n_matrices, 
    gsl_matrix_complex *m_product
) {
    int sign, 
        dir_abs;
    const gsl_complex z_zero = gsl_complex_rect(0., 0.), 
                    z_one = gsl_complex_rect(1., 0.);
    
    if (dir < 4) { 
        sign = 1;
        dir_abs = dir;
    } else {
        sign = -1;
        dir_abs = dir - 4;
    }

    for (int x = 0; x < n_matrices; x++) {
        gsl_blas_zgemm(
            CblasNoTrans, 
            CblasNoTrans, 
            z_one, 
            m_product, 
            lattice[periodic_ind(
                i_start + sign * x * (dir_abs == 0), 
                j_start + sign * x * (dir_abs == 1), 
                k_start + sign * x * (dir_abs == 2), 
                l_start + sign * x * (dir_abs == 3), 
                dir, 
                par->L
            )], 
            z_zero, 
            par->m_workspace
        );
        gsl_matrix_complex_memcpy(m_product, par->m_workspace);
    }
}

int psl_matrix_complex_unitarize(gsl_matrix_complex *matrix) {
    gsl_complex z_temp;
    gsl_vector_complex_view vec[2];
    gsl_vector_complex *v_temp = NULL;

    v_temp = gsl_vector_complex_alloc(2);
    if (v_temp == NULL) {
        printf("Error: Failed allocating memory for temporary vector used for unitarization. Exiting...\n");
        return 1;
    }
    
    for (int i = 0; i < 2; i++) 
        vec[i] = gsl_matrix_complex_column(matrix, i);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < i; j++) {   /* i,j each go over the column vectors of the matrix */
            gsl_blas_zdotc(&(vec[j].vector), &(vec[i].vector), &z_temp);
            gsl_vector_complex_memcpy(v_temp, &(vec[j].vector));
            gsl_vector_complex_scale(v_temp, z_temp);
            gsl_vector_complex_sub(&(vec[i].vector), v_temp);
        }
        gsl_blas_zdscal(1. / gsl_blas_dznrm2(&(vec[i].vector)), &(vec[i].vector));  /* normalization */
    }
    
    /* DEBUG: 
    gsl_blas_zdotc(&(vec[1].vector), &(vec[2].vector), &z_temp);
    printf("DEBUG: dot(1,2) = %4.3f + %4.3fi\n", GSL_REAL(z_temp), GSL_IMAG(z_temp)); */

    /* cross product: 
    psl_vector_complex_cross(&(vec[0].vector), &(vec[1].vector), &(vec[2].vector));
    */
    gsl_vector_complex_free(v_temp);
    return 0;
}

int measure(PAR *par, gsl_matrix *results, gsl_matrix_complex **lattice) {
    const gsl_complex z_zero = gsl_complex_rect(0., 0.);
    double result;
    
    gsl_matrix_set_zero(results);

    for (int i; i < par->L; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    for (int dir_1 = 0; dir_1 < 4; dir_1++) {
                        for (int dir_2 = dir_1 + 1; dir_2 < 4; dir_2++) {
                            for (int L_1 = 1; L_1 < par->L; L_1++) {
                                for (int L_2 = L_1; L_2 < par->L; L_2++) {
                                    if(lattice_loop_rect(
                                        par, 
                                        lattice, 
                                        i, 
                                        j, 
                                        k, 
                                        l, 
                                        dir_1, 
                                        dir_2, 
                                        L_1, 
                                        L_2, 
                                        &result
                                    )) 
                                        return 1;
                                    gsl_matrix_set(
                                        results, 
                                        L_1 - 1, 
                                        L_2 - 1, 
                                        gsl_matrix_get(results, L_1 - 1, L_2 - 1) + result
                                    );
                                    /* save to file here */

                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    /* normalization missing! */

    return 0;
}

void measure_tadpole(PAR *par, gsl_matrix_complex **lattice, double *tadpole_result) {
    const gsl_complex z_zero = gsl_complex_rect(1., 0.);
    gsl_complex z_temp;

    *tadpole_result = 0.;

    for (int i = 0; i < par->L; i++) {
        for (int j = 0; j < par->L; j++) {
            for(int k = 0; k < par->L; k++) {
                for(int l = 0; l < par->L; l++) {
                    for (int dir = 0; dir < 4; dir++) {
                        z_temp = z_zero;
                        for (int m = 0; m < 2; m++) {
                            z_temp = gsl_complex_add(
                                z_temp, 
                                gsl_matrix_complex_get(
                                    lattice[ind(i, j, k, l, dir, par->L)], 
                                    m, 
                                    m
                                )
                            );
                        }
                        *tadpole_result += GSL_REAL(z_temp);
                    }
                }
            }
        }
    }
    *tadpole_result /= (double)(2 * par->L * par->L * par->L * par->L * 4);
}

int measure_aa_a2a(PAR *par, double *results, gsl_matrix_complex **lattice) {
    const gsl_complex z_zero = gsl_complex_rect(0., 0.); 
    gsl_complex z_temp;
    gsl_matrix_complex *m_temp_1 = NULL; 

    results[0] = 0.;
    results[1] = 0.;

    m_temp_1 = gsl_matrix_complex_alloc(2, 2);
    if (m_temp_1 == NULL) {
        printf("Error: Failed allocating memory for temporary matrices used in the measurement function. Exiting...\n");
        return 1;
    }
    
    for (int i = 0; i < par->L; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    for (int dir_1 = 0; dir_1 < 4; dir_1++) {
                        for (int dir_2 = dir_1 + 1; dir_2 < 4; dir_2++) {
                            /* calculate a x a Wilson-loop */
                            psl_matrix_complex_product_4(
                                par,
                                lattice[ind(i, j, k, l, dir_1, par->L)],
                                lattice[periodic_ind(
                                    i + (dir_1 == 0), 
                                    j + (dir_1 == 1), 
                                    k + (dir_1 == 2), 
                                    l + (dir_1 == 3), 
                                    dir_2,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_1 == 0) + (dir_2 == 0), 
                                    j + (dir_1 == 1) + (dir_2 == 1), 
                                    k + (dir_1 == 2) + (dir_2 == 2), 
                                    l + (dir_1 == 3) + (dir_2 == 3), 
                                    dir_1 + 4,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_2 == 0), 
                                    j + (dir_2 == 1), 
                                    k + (dir_2 == 2), 
                                    l + (dir_2 == 3), 
                                    dir_2 + 4,
                                    par->L
                                )],
                                m_temp_1
                            );
                            z_temp = z_zero;
                            for (int m = 0; m < 2; m++)     /* trace */
                                z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp_1, m, m));
                            results[0] += GSL_REAL(z_temp);

                            /* calculate first a x 2a Wilson-loop, while reusing the contentof m_workspace from the a x a Wilson-loop (it's kind of black magic... DON'T TRY AT HOME) */
                            psl_matrix_complex_product_4(
                                par,
                                par->m_workspace,   /* <-- dangerous stuff */
                                lattice[periodic_ind(
                                    i + (dir_2 == 0), 
                                    j + (dir_2 == 1), 
                                    k + (dir_2 == 2), 
                                    l + (dir_2 == 3), 
                                    dir_1 + 4,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_2 == 0) - (dir_1 == 0), 
                                    j + (dir_2 == 1) - (dir_1 == 1), 
                                    k + (dir_2 == 2) - (dir_1 == 2), 
                                    l + (dir_2 == 3) - (dir_1 == 3), 
                                    dir_2 + 4,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i - (dir_1 == 0), 
                                    j - (dir_1 == 1), 
                                    k - (dir_1 == 2), 
                                    l - (dir_1 == 3), 
                                    dir_1,
                                    par->L
                                )],
                                m_temp_1
                            );

                            z_temp = z_zero;
                            for (int m = 0; m < 2; m++)     /* trace */
                                z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp_1, m, m));
                            results[1] += GSL_REAL(z_temp);

                            /* calculate the second a x 2a Wilson-loop (no black magic here) */
                            psl_matrix_complex_product_6(
                                par,
                                lattice[ind(i, j, k, l, dir_2, par->L)],
                                lattice[periodic_ind(
                                    i + (dir_2 == 0), 
                                    j + (dir_2 == 1), 
                                    k + (dir_2 == 2), 
                                    l + (dir_2 == 3), 
                                    dir_1,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_1 == 0) + (dir_2 == 0), 
                                    j + (dir_1 == 1) + (dir_2 == 1), 
                                    k + (dir_1 == 2) + (dir_2 == 2), 
                                    l + (dir_1 == 3) + (dir_2 == 3), 
                                    dir_2 + 4,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_1 == 0), 
                                    j + (dir_1 == 1), 
                                    k + (dir_1 == 2), 
                                    l + (dir_1 == 3), 
                                    dir_2 + 4,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_1 == 0) - (dir_2 == 0), 
                                    j + (dir_1 == 1) - (dir_2 == 1), 
                                    k + (dir_1 == 2) - (dir_2 == 2), 
                                    l + (dir_1 == 3) - (dir_2 == 3), 
                                    dir_1 + 4,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i - (dir_2 == 0), 
                                    j - (dir_2 == 1), 
                                    k - (dir_2 == 2), 
                                    l - (dir_2 == 3), 
                                    dir_2,
                                    par->L
                                )],
                                m_temp_1
                            );
                            
                            z_temp = z_zero;
                            for (int m = 0; m < 2; m++)     /* trace */
                                z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp_1, m, m));
                            results[1] += GSL_REAL(z_temp);
                        }
                    }
                }
            }
        }
    }
    /* factor 1/2 from Wilson-loop and normalization */
    results[0] /= (double)(2 * par->L * par->L * par->L * par->L * 6); 
    results[1] /= (double)(2 * par->L * par->L * par->L * par->L * 12);
    return 0;
}

int measure_action_r(PAR *par, gsl_matrix_complex **lattice, double *action) {
    const gsl_complex z_zero = gsl_complex_rect(0., 0.); 
    gsl_complex z_temp;
    gsl_matrix_complex *m_temp_1 = NULL; 

    *action = 0.;

    m_temp_1 = gsl_matrix_complex_alloc(2, 2);
    if (m_temp_1 == NULL) {
        printf("Error: Failed allocating memory for temporary matrices used in the action-measurement function. Exiting...\n");
        return 1;
    }
    
    for (int i = 0; i < par->L; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    for (int dir_1 = 0; dir_1 < 4; dir_1++) {
                        for (int dir_2 = dir_1 + 1; dir_2 < 4; dir_2++) {
                            /* calculate a x a Wilson-loop */
                            psl_matrix_complex_product_4(
                                par,
                                lattice[ind(i, j, k, l, dir_1, par->L)],
                                lattice[periodic_ind(
                                    i + (dir_1 == 0), 
                                    j + (dir_1 == 1), 
                                    k + (dir_1 == 2), 
                                    l + (dir_1 == 3), 
                                    dir_2,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_1 == 0) + (dir_2 == 0), 
                                    j + (dir_1 == 1) + (dir_2 == 1), 
                                    k + (dir_1 == 2) + (dir_2 == 2), 
                                    l + (dir_1 == 3) + (dir_2 == 3), 
                                    dir_1 + 4,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_2 == 0), 
                                    j + (dir_2 == 1), 
                                    k + (dir_2 == 2), 
                                    l + (dir_2 == 3), 
                                    dir_2 + 4,
                                    par->L
                                )],
                                m_temp_1
                            );
                            z_temp = z_zero;
                            for (int m = 0; m < 2; m++)     /* trace */
                                z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp_1, m, m));
                            *action += GSL_REAL(z_temp);
                        }
                    }
                }
            }
        }
    }
    /* factor 1/2 from Plaquette operator and -beta from Wilson action */
    *action *= -par->beta / 2. / (double)(par->tadpole * par->tadpole * par->tadpole * par->tadpole); 
    return 0;
}

int measure_action_l(PAR *par, gsl_matrix_complex **lattice, double *action) {
    const gsl_complex z_zero = gsl_complex_rect(0., 0.); 
    gsl_complex z_temp;
    gsl_matrix_complex *m_temp_1 = NULL; 

    *action = 0.;

    m_temp_1 = gsl_matrix_complex_alloc(2, 2);
    if (m_temp_1 == NULL) {
        printf("Error: Failed allocating memory for temporary matrices used in the action-measurement function. Exiting...\n");
        return 1;
    }
    
    for (int i = 0; i < par->L; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    for (int dir_2 = 0; dir_2 < 4; dir_2++) {
                        for (int dir_1 = dir_2 + 1; dir_1 < 4; dir_1++) {
                            /* calculate a x a Wilson-loop */
                            psl_matrix_complex_product_4(
                                par,
                                lattice[ind(i, j, k, l, dir_1, par->L)],
                                lattice[periodic_ind(
                                    i + (dir_1 == 0), 
                                    j + (dir_1 == 1), 
                                    k + (dir_1 == 2), 
                                    l + (dir_1 == 3), 
                                    dir_2,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_1 == 0) + (dir_2 == 0), 
                                    j + (dir_1 == 1) + (dir_2 == 1), 
                                    k + (dir_1 == 2) + (dir_2 == 2), 
                                    l + (dir_1 == 3) + (dir_2 == 3), 
                                    dir_1 + 4,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_2 == 0), 
                                    j + (dir_2 == 1), 
                                    k + (dir_2 == 2), 
                                    l + (dir_2 == 3), 
                                    dir_2 + 4,
                                    par->L
                                )],
                                m_temp_1
                            );
                            z_temp = z_zero;
                            for (int m = 0; m < 2; m++)     /* trace */
                                z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp_1, m, m));
                            *action += GSL_REAL(z_temp);
                        }
                    }
                }
            }
        }
    }
    /* factor 1/2 from Plaquette operator and -beta from Wilson action */
    *action *= -par->beta / 2. / (double)(par->tadpole * par->tadpole * par->tadpole * par->tadpole); 
    return 0;
}

int unitarize_lattice(PAR *par, gsl_matrix_complex **lattice) {
    for (int i = 0; i < par->L; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    for (int dir = 0; dir < 4; dir++) {
                        if (psl_matrix_complex_unitarize(lattice[ind(i, j, k, l, dir, par->L)])) 
                            return 1;
                    }
                }
            }
        }
    }
    for (int i = 0; i < par->L; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    for (int dir = 0; dir < 4; dir++) {
                        psl_matrix_complex_dagger_memcpy(
                            lattice[ind(i, j, k, l, dir + 4, par->L)], 
                            lattice[periodic_ind(
                                i - (dir == 0), 
                                j - (dir == 1), 
                                k - (dir == 2), 
                                l - (dir == 3), 
                                dir, 
                                par->L
                            )]
                        );
                    }
                }
            }
        }
    }

    return 0;
}

int update_lattice(PAR *par, gsl_matrix_complex **lattice, gsl_matrix_complex **su2, double * acceptance) {
    double Delta_S;
    const gsl_complex z_zero = gsl_complex_rect(0., 0.), 
          z_one = gsl_complex_rect(1., 0.);
    gsl_complex z_temp;
    gsl_matrix_complex *m_temp_1 = NULL, 
                       *m_temp_2 = NULL, 
                       *m_temp_3 = NULL, 
                       *Gamma = NULL;

    m_temp_1 = gsl_matrix_complex_alloc(2, 2);
    m_temp_2 = gsl_matrix_complex_alloc(2, 2);
    m_temp_3 = gsl_matrix_complex_alloc(2, 2);
    Gamma = gsl_matrix_complex_alloc(2, 2);
    if ((m_temp_1 == NULL) || (m_temp_2 == NULL) || (m_temp_3 == NULL) || (Gamma == NULL)) {
        printf("Error: Failed allocating memory for temporary matrix/matrices and/or Gamma matrix (lattice update). Exiting...\n");
        return 1;
    } 

    /* Loop through the whole lattice (all independent links) */
    for (int i = 0; i < par->L; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    for (int dir_1 = 0; dir_1 < 4; dir_1++) {
                        /* calculate Gamma_mu(x) which is needed to calculate the differences in action later */
                        gsl_matrix_complex_set_zero(Gamma);
                        for (int dir_2 = 0; dir_2 < 4; dir_2++) {
                            if (dir_1 == dir_2) 
                                continue;
                            psl_matrix_complex_product_3_add(
                                par,
                                lattice[periodic_ind(
                                    i + (dir_1 == 0), 
                                    j + (dir_1 == 1),
                                    k + (dir_1 == 2),
                                    l + (dir_1 == 3),
                                    dir_2,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_1 == 0) + (dir_2 == 0),
                                    j + (dir_1 == 1) + (dir_2 == 1),
                                    k + (dir_1 == 2) + (dir_2 == 2),
                                    l + (dir_1 == 3) + (dir_2 == 3),
                                    dir_1 + 4,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_2 == 0), 
                                    j + (dir_2 == 1),
                                    k + (dir_2 == 2),
                                    l + (dir_2 == 3),
                                    dir_2 + 4,
                                    par->L
                                )],
                                Gamma
                            );
                            psl_matrix_complex_product_3_add(
                                par,
                                lattice[periodic_ind(
                                    i + (dir_1 == 0), 
                                    j + (dir_1 == 1),
                                    k + (dir_1 == 2),
                                    l + (dir_1 == 3),
                                    dir_2 + 4,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i + (dir_1 == 0) - (dir_2 == 0),
                                    j + (dir_1 == 1) - (dir_2 == 1),
                                    k + (dir_1 == 2) - (dir_2 == 2),
                                    l + (dir_1 == 3) - (dir_2 == 3),
                                    dir_1 + 4,
                                    par->L
                                )],
                                lattice[periodic_ind(
                                    i - (dir_2 == 0), 
                                    j - (dir_2 == 1),
                                    k - (dir_2 == 2),
                                    l - (dir_2 == 3),
                                    dir_2,
                                    par->L
                                )],
                                Gamma
                            );
                        }
                        /* do the specified number of MC-updates of this link */
                        for (int n = 0; n < par->n_hits; n++) {
                            /* DEBUG: 
                            double action_1, action_2;
                            if (measure_action_r(par, lattice, &action_1)) {
                                gsl_matrix_complex_free(m_temp_1);
                                gsl_matrix_complex_free(m_temp_2);
                                gsl_matrix_complex_free(m_temp_3);
                                gsl_matrix_complex_free(Gamma);
                                return 1;
                            }  */
                            
                            /* multiply with a random SU(2)-matrix and save in a temporary matrix, until we
                             * know if we will use it */
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                su2[(int)(par->n_su2 * gsl_rng_uniform(par->ran_gen))], 
                                lattice[ind(i, j, k, l, dir_1, par->L)],
                                z_zero, 
                                m_temp_1
                            );
                            /* measure the difference in action the new link operator does for the lattice */
                            gsl_matrix_complex_memcpy(m_temp_2, m_temp_1);
                            gsl_matrix_complex_sub(m_temp_2, lattice[ind(i, j, k, l, dir_1, par->L)]);
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans,
                                z_one,
                                m_temp_2,
                                Gamma,
                                z_zero,
                                m_temp_3
                            );
                            z_temp = z_zero;
                            for (int m = 0; m < 2; m++)     /* trace */
                                z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp_3, m, m));
                            Delta_S = -par->beta * GSL_REAL(z_temp) / 2. / (double)(par->tadpole * par->tadpole * par->tadpole * par->tadpole);
                            /* accept/reject step */
                            if (Delta_S <= 0.) {
                                /* save updated matrix and it's conjugate transpose */
                                gsl_matrix_complex_memcpy(lattice[ind(i, j, k, l, dir_1, par->L)], m_temp_1);
                                psl_matrix_complex_dagger_memcpy(
                                    lattice[periodic_ind(
                                        i + (dir_1 == 0), 
                                        j + (dir_1 == 1), 
                                        k + (dir_1 == 2), 
                                        l + (dir_1 == 3),
                                        dir_1 + 4, 
                                        par->L
                                    )],
                                    m_temp_1
                                );
                                *acceptance += 1.;

                                /* DEBUG: 
                                if (measure_action_r(par, lattice, &action_2)) {
                                    gsl_matrix_complex_free(m_temp_1);
                                    gsl_matrix_complex_free(m_temp_2);
                                    gsl_matrix_complex_free(m_temp_3);
                                    gsl_matrix_complex_free(Gamma);
                                    return 1;
                                }
                                printf("DEBUG: DS = %7.2e\nDEBUG: SS = %7.2e\n", Delta_S, action_2 - action_1);
                                if (fabs(Delta_S) > 1e-12) {
                                    for (int a = 0; a < 2; a++) {
                                        for (int b = 0; b < 2; b++) {
                                            z_temp = gsl_matrix_complex_get(lattice[ind(i, j, k, l, dir_1, par->L)], a, b);
                                            printf("%7.2e+i%7.2e\t", GSL_REAL(z_temp), GSL_IMAG(z_temp));
                                        }
                                        printf("\n");
                                    }
                                } */
                                
                            } else {
                                if (gsl_rng_uniform(par->ran_gen) < exp(-Delta_S)) {
                                    /* save updated matrix and it's conjugate transpose */
                                    gsl_matrix_complex_memcpy(lattice[ind(i, j, k, l, dir_1, par->L)], m_temp_1);
                                    psl_matrix_complex_dagger_memcpy(
                                        lattice[periodic_ind(
                                            i + (dir_1 == 0), 
                                            j + (dir_1 == 1), 
                                            k + (dir_1 == 2), 
                                            l + (dir_1 == 3),
                                            dir_1 + 4, 
                                            par->L
                                        )],
                                        m_temp_1
                                    );
                                    *acceptance += 1.;
                                    /* DEBUG: 
                                    if (measure_action_r(par, lattice, &action_2)) {
                                        gsl_matrix_complex_free(m_temp_1);
                                        gsl_matrix_complex_free(m_temp_2);
                                        gsl_matrix_complex_free(m_temp_3);
                                        gsl_matrix_complex_free(Gamma);
                                        return 1;
                                    }
                                    printf("DEBUG: DS = %7.2e\nDEBUG: SS = %7.2e\n", Delta_S, action_2 - action_1);
                                    if (fabs(Delta_S) > 1e-12) {
                                        for (int a = 0; a < 2; a++) {
                                            for (int b = 0; b < 2; b++) {
                                                z_temp = gsl_matrix_complex_get(lattice[ind(i, j, k, l, dir_1, par->L)], a, b);
                                                printf("%7.2e+i%7.2e\t", GSL_REAL(z_temp), GSL_IMAG(z_temp));
                                            }
                                            printf("\n");
                                        }
                                    } */
                                
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    gsl_matrix_complex_free(m_temp_1);
    gsl_matrix_complex_free(m_temp_2);
    gsl_matrix_complex_free(m_temp_3);
    gsl_matrix_complex_free(Gamma);

    return 0;
}


/* int init_su3(PAR *par, gsl_matrix_complex **su3) {
    * generate random matrices *
    for (int n = 0; n < par->n_su2 / 2; n++) {
        for (int i = 0; i < 2; i++) {
            for (int j = i; j < 2; j++) {   * su2_{i,j} are the matrix components *
                gsl_matrix_complex_set(
                    su2[n], 
                    i, 
                    j, 
                    gsl_complex_rect(
                        (double)(i == j), 
                        par->eps * (gsl_rng_uniform_pos(par->ran_gen) * 2. - 1.)
                    )
                );
            }
        }
        * make the matrix symmetric *
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < i; j++) {
                gsl_matrix_complex_set(
                    su2[n], 
                    i, 
                    j, 
                    gsl_matrix_complex_get(su2[n], j, i)
                );
            }
        }
        
        if (psl_matrix_complex_unitarize(su2[n])) {
            return 1;
        }

        * create inverse matrices by getting the conjugate transpose *
        psl_matrix_complex_dagger_memcpy(su2[n + par->n_su2 / 2], su2[n]);
    }

    * DEBUG:  
    gsl_complex z_temp;
    const gsl_complex z_zero = gsl_complex_rect(0., 0.), 
                    z_one = gsl_complex_rect(1., 0.);
    gsl_matrix_complex *m_temp = NULL, 
                       *m_temp_0 = NULL;
    gsl_permutation *p = NULL;
    
    m_temp = gsl_matrix_complex_alloc(2, 2);
    m_temp_0 = gsl_matrix_complex_alloc(2, 2);
    p = gsl_permutation_alloc(3);
    if ((m_temp == NULL) || (p == NULL) || (m_temp_0 == NULL)) {
        printf("Error: Failed allocating memory for a temporary matrix and/or permutation used in the debugging process. Exiting...\n");
        return 1;
    }
    for (int i = 0; i < par->n_su2; i++) {
        int signum;
        gsl_complex det;
        
        gsl_matrix_complex_memcpy(m_temp, su2[i]);
        gsl_linalg_complex_LU_decomp(m_temp, p, &signum);
        det = gsl_linalg_complex_LU_det(m_temp, signum);
        printf("DEBUG: det = %4.3f + %4.3fi\n", GSL_REAL(det), GSL_IMAG(det));
        printf("DEBUG: abs(det) = %4.3f\n", gsl_complex_abs(det));
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
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
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
    for (int i = 0; i < 3; i++) 
        z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp_0, i, i));
    printf("DEBUG: Trace = %e\n", GSL_REAL(z_temp));

    gsl_matrix_complex_free(m_temp);
    gsl_matrix_complex_free(m_temp_0);
    gsl_permutation_free(p); *

    return 0;
} */

void init_lattice(PAR *par, gsl_matrix_complex **lattice, gsl_matrix_complex **su2) {
    /* initialize all independent links with random SU(3) matrices */
    for (int i = 0; i < par->L; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    for (int dir = 0; dir < 4; dir++) {
                        gsl_matrix_complex_memcpy(
                            lattice[ind(i, j, k, l, dir, par->L)], 
                            su2[(int)(gsl_rng_uniform(par->ran_gen) * par->n_su2)]
                        );
                    }
                }
            }
        }
    }

    /* generate the daggered matrices */
    for (int i = 0; i < par->L; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    for (int dir_dagger = 4; dir_dagger < 8; dir_dagger++) {
                        psl_matrix_complex_dagger_memcpy(
                            lattice[ind(i, j, k, l, dir_dagger, par->L)], 
                            lattice[periodic_ind(
                                i - (dir_dagger == 4), 
                                j - (dir_dagger == 5), 
                                k - (dir_dagger == 6), 
                                l - (dir_dagger == 7), 
                                dir_dagger - 4, 
                                par->L
                            )]
                        );
                    }
                }
            }
        }
    }
}





