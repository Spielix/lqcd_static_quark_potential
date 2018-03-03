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

#include "static_quark_potential.h"

/* calculate the position of a specific link in the array */
static inline int ind(int i, int j, int k, int l, int dir, int le) {
    return (((le * i + j) * le + k) * le + l) * le + dir;
}

/* calculate the position of a specific link in the array, while applying periodic boundary conditions to all
 * coordinates */
static inline int periodic_ind(int i, int j, int k, int l, int dir, int le) {
    return (((le * (i & (le - 1)) + (j & (le - 1))) * le + (k & (le - 1))) * le + (l & (le - 1))) * le + dir;
}

/* function to get the conjugate transpose matrix */
void psl_matrix_complex_dagger(gsl_matrix_complex *m) {
    gsl_matrix_complex_transpose(m);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            gsl_matrix_complex_set(m, i, j, gsl_complex_conjugate(gsl_matrix_complex_get(m, i, j)));
    }
}

/* function to save the conjugate transpose matrix of src under dest */
void psl_matrix_complex_dagger_memcpy(gsl_matrix_complex *dest, gsl_matrix_complex *src) {
    gsl_matrix_complex_transpose_memcpy(dest, src); 
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            gsl_matrix_complex_set(dest, i, j, gsl_complex_conjugate(gsl_matrix_complex_get(dest, i, j)));
    }
}

int measure(Par *par, double *results, gsl_matrix_complex **lattice) {
    const gsl_complex z_zero = gsl_complex_rect(0., 0.), 
          z_one = gsl_complex_rect(1., 0.);
    gsl_complex z_temp;
    gsl_matrix_complex *m_temp_1 = NULL, 
                       *m_temp_2 = NULL;

    results[0] = 0.;
    results[1] = 0.;

    m_temp_1 = gsl_matrix_complex_alloc(3, 3);
    m_temp_2 = gsl_matrix_complex_alloc(3, 3);
    if ((m_temp_1 == NULL) || (m_temp_2 == NULL)) {
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
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                lattice[ind(i, j, k, l, dir_1, par->L)],
                                lattice[periodic_ind(
                                    i + (dir_1 == 0), 
                                    j + (dir_1 == 1), 
                                    k + (dir_1 == 2), 
                                    l + (dir_1 == 3), 
                                    dir_2,
                                    par->L
                                )],
                                z_zero,
                                m_temp_1
                            );
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                m_temp_1,
                                lattice[periodic_ind(
                                    i + (dir_1 == 0) + (dir_2 == 0), 
                                    j + (dir_1 == 1) + (dir_2 == 0), 
                                    k + (dir_1 == 2) + (dir_2 == 0), 
                                    l + (dir_1 == 3) + (dir_2 == 0), 
                                    dir_1 + 4,
                                    par->L
                                )],
                                z_zero,
                                m_temp_2
                            );
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                m_temp_2,
                                lattice[periodic_ind(
                                    i + (dir_2 == 0), 
                                    j + (dir_2 == 1), 
                                    k + (dir_2 == 2), 
                                    l + (dir_2 == 3), 
                                    dir_2 + 4,
                                    par->L
                                )],
                                z_zero,
                                m_temp_1
                            );
                            z_temp = z_zero;
                            for (int m = 0; m < 3; m++)     /* trace */
                                z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp_1, m, m));
                            results[0] += GSL_REAL(z_temp);

                            /* calculate first a x 2a Wilson-loop, while reusing the contentof m_temp_2 from the a x a Wilson-loop */
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                m_temp_2,
                                lattice[periodic_ind(
                                    i + (dir_2 == 0), 
                                    j + (dir_2 == 1), 
                                    k + (dir_2 == 2), 
                                    l + (dir_2 == 3), 
                                    dir_1 + 4,
                                    par->L
                                )],
                                z_zero,
                                m_temp_1
                            );
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                m_temp_1,
                                lattice[periodic_ind(
                                    i + (dir_2 == 0) - (dir_1 == 0), 
                                    j + (dir_2 == 1) - (dir_1 == 1), 
                                    k + (dir_2 == 2) - (dir_1 == 2), 
                                    l + (dir_2 == 3) - (dir_1 == 3), 
                                    dir_2 + 4,
                                    par->L
                                )],
                                z_zero,
                                m_temp_2
                            );
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                m_temp_2,
                                lattice[periodic_ind(
                                    i - (dir_1 == 0), 
                                    j - (dir_1 == 1), 
                                    k - (dir_1 == 2), 
                                    l - (dir_1 == 3), 
                                    dir_1,
                                    par->L
                                )],
                                z_zero,
                                m_temp_1
                            );

                            z_temp = z_zero;
                            for (int m = 0; m < 3; m++)     /* trace */
                                z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp_1, m, m));
                            results[1] += GSL_REAL(z_temp);

                            /* calculate the second a x 2a Wilson-loop */
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                lattice[ind(i, j, k, l, dir_2, par->L)],
                                lattice[periodic_ind(
                                    i + (dir_2 == 0), 
                                    j + (dir_2 == 1), 
                                    k + (dir_2 == 2), 
                                    l + (dir_2 == 3), 
                                    dir_1,
                                    par->L
                                )],
                                z_zero,
                                m_temp_1
                            );
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                m_temp_1,
                                lattice[periodic_ind(
                                    i + (dir_1 == 0) + (dir_2 == 0), 
                                    j + (dir_1 == 1) + (dir_2 == 1), 
                                    k + (dir_1 == 2) + (dir_2 == 2), 
                                    l + (dir_1 == 3) + (dir_2 == 3), 
                                    dir_2 + 4,
                                    par->L
                                )],
                                z_zero,
                                m_temp_2
                            );
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                m_temp_2,
                                lattice[periodic_ind(
                                    i + (dir_1 == 0), 
                                    j + (dir_1 == 1), 
                                    k + (dir_1 == 2), 
                                    l + (dir_1 == 3), 
                                    dir_2 + 4,
                                    par->L
                                )],
                                z_zero,
                                m_temp_1
                            );
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                m_temp_1,
                                lattice[periodic_ind(
                                    i + (dir_1 == 0) - (dir_2 == 0), 
                                    j + (dir_1 == 1) - (dir_2 == 1), 
                                    k + (dir_1 == 2) - (dir_2 == 2), 
                                    l + (dir_1 == 3) - (dir_2 == 3), 
                                    dir_1 + 4,
                                    par->L
                                )],
                                z_zero,
                                m_temp_2
                            );
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                m_temp_2,
                                lattice[periodic_ind(
                                    i - (dir_2 == 0), 
                                    j - (dir_2 == 1), 
                                    k - (dir_2 == 2), 
                                    l - (dir_2 == 3), 
                                    dir_2,
                                    par->L
                                )],
                                z_zero,
                                m_temp_1
                            );
                            
                            z_temp = z_zero;
                            for (int m = 0; m < 3; m++)     /* trace */
                                z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp_1, m, m));
                            results[1] += GSL_REAL(z_temp);
                        }
                    }
                }
            }
        }
    }
    /* factor 1/3 from Wilson-loop and normalization */
    results[0] /= (double)(3 * par->L * par->L * par->L * par->L * 6); 
    results[1] /= (double)(3 * par->L * par->L * par->L * par->L * 12);
    return 0;
}


void print_results(Par *par, double *results) {
}

int update_lattice(Par *par, gsl_matrix_complex **lattice, gsl_matrix_complex **su3, double * acceptance) {
    double Delta_S;
    const gsl_complex z_zero = gsl_complex_rect(0., 0.), 
          z_one = gsl_complex_rect(1., 0.);
    gsl_complex z_temp;
    gsl_matrix_complex *m_temp_1 = NULL, 
                       *m_temp_2 = NULL, 
                       *m_temp_3 = NULL, 
                       *Gamma = NULL;

    m_temp_1 = gsl_matrix_complex_alloc(3, 3);
    m_temp_2 = gsl_matrix_complex_alloc(3, 3);
    m_temp_3 = gsl_matrix_complex_alloc(3, 3);
    Gamma = gsl_matrix_complex_alloc(3, 3);
    if ((m_temp_1 == NULL) || (m_temp_2 == NULL) || (m_temp_3 == NULL) || (Gamma == NULL)) {
        printf("Error: Failed allocating memory for temporary matrix/matrices and/or Gamma matrix (lattice update). Exiting...\n");
        return 1;
    }

    *acceptance = 0.;

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
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans,
                                z_one,
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
                                z_zero,
                                m_temp_1
                            );
                            gsl_blas_zgemm(
                                CblasNoTrans,
                                CblasNoTrans,
                                z_one,
                                m_temp_1,
                                lattice[periodic_ind(
                                    i + (dir_2 == 0), 
                                    j + (dir_2 == 1),
                                    k + (dir_2 == 2),
                                    l + (dir_2 == 3),
                                    dir_2 + 4,
                                    par->L
                                )],
                                z_one,
                                Gamma
                            );
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans,
                                z_one,
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
                                z_zero,
                                m_temp_1
                            );
                            gsl_blas_zgemm(
                                CblasNoTrans,
                                CblasNoTrans,
                                z_one,
                                m_temp_1,
                                lattice[periodic_ind(
                                    i - (dir_2 == 0), 
                                    j - (dir_2 == 1),
                                    k - (dir_2 == 2),
                                    l - (dir_2 == 3),
                                    dir_2,
                                    par->L
                                )],
                                z_one,
                                Gamma
                            );
                        }
                        /* do the specified number of MC-updates of this link */
                        for (int n = 0; n < par->n_hits; n++) {
                            gsl_blas_zgemm(
                                CblasNoTrans, 
                                CblasNoTrans, 
                                z_one, 
                                su3[(int)(par->n_su3 * gsl_rng_uniform(par->ran_gen))], 
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
                            for (int m = 0; m < 3; m++)     /* trace */
                                z_temp = gsl_complex_add(z_temp, gsl_matrix_complex_get(m_temp_3, m, m));
                            Delta_S = -par->beta * GSL_REAL(z_temp);     /* ------ Faktoren??? ------ */
                            /* accept/reject step */
                            if (Delta_S < 0.) {
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
                                continue;
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
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    *acceptance /= (double)(par->L * par->L * par->L * par->L * 4 * par->n_hits);

    gsl_matrix_complex_free(m_temp_1);
    gsl_matrix_complex_free(m_temp_2);
    gsl_matrix_complex_free(m_temp_3);
    gsl_matrix_complex_free(Gamma);

    return 0;
}


int init_su3(Par *par, gsl_matrix_complex **su3) {
    gsl_vector_complex *vec[3], 
                       *tempv;
    gsl_complex tempz;

    /* allocating memory for three vectors which will temporarily store the column vectors of the matrices for
     * unitarization */
    for (int i = 0; i < 3; i++) {
        vec[i] = gsl_vector_complex_alloc(3);
    }
    tempv = gsl_vector_complex_alloc(3);
    if ((vec[0] == NULL) || (vec[1] == NULL) || (vec[2] == NULL) || (tempv == NULL)) {
        printf("Error: Allocating memory for vector in init_su3 failed. Exiting...\n");
        return 1;
    }

    for (int n = 0; n < par->n_su3 / 2; n++) {
        for (int i = 0; i < 3; i++) {
            for (int j = i; j < 3; j++) {   /* su3_{i,j} are the matrix components */
                gsl_matrix_complex_set(
                    su3[n], 
                    i, 
                    j, 
                    gsl_complex_rect(
                        (double)(i == j), 
                        par->eps * (gsl_rng_uniform_pos(par->ran_gen) * 2. - 1.)
                    )
                );
            }
        }
        /* make the matrix symmetric */
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < i; j++) {
                gsl_matrix_complex_set(
                    su3[n], 
                    i, 
                    j, 
                    gsl_matrix_complex_get(su3[n], j, i)
                );
            }
        }
        /* unitarize the matrix (Gram-Schmidt) */
        for (int i = 0; i < 3; i++) {
            gsl_matrix_complex_get_col(vec[i], su3[n], i);
            for (int j = 0; j < i; j++) {   /* i,j each go over the column vectors of the matrix */
                gsl_blas_zdotu(vec[j], vec[i], &tempz);
                gsl_vector_complex_memcpy(tempv, vec[j]);
                gsl_vector_complex_scale(tempv, tempz);
                gsl_vector_complex_sub(vec[i], tempv);
            }
            gsl_blas_zdscal(1. / gsl_blas_dznrm2(vec[i]), vec[i]);  /* normalization */
            gsl_matrix_complex_set_col(su3[n], i, vec[i]);
        }
        /* for (int j = 0; j < 3; j++) {    cross product 
            gsl_vector_complex_set(
                vec[2], 
                j, 
                gsl_complex_add(
                    gsl_complex_mul(
                        gsl_vector_complex_get(vec[1], (j + 1) % 3), 
                        gsl_vector_complex_get(vec[2], (j + 2) % 3)
                    ), 
                    gsl_complex_mul(
                        gsl_vector_complex_get(vec[1], (j + 2) % 3),
                        gsl_vector_complex_get(vec[2], (j + 1) % 3)
                    )
                )
            );
        } */
    }

    /* create inverse matrices by getting the conjugate transpose */
    for (int n = par->n_su3 / 2; n < par->n_su3; n++) {
        psl_matrix_complex_dagger_memcpy(su3[n], su3[n - par->n_su3 / 2]);
    }

    for (int i = 0; i < 3; i++)
        gsl_vector_complex_free(vec[i]);
    gsl_vector_complex_free(tempv);

    /* debugging: */
    const gsl_complex z_zero = gsl_complex_rect(0., 0.), z_one = gsl_complex_rect(1., 0.);
    gsl_matrix_complex *m_temp = NULL;
    gsl_permutation *p = NULL;
    
    m_temp = gsl_matrix_complex_alloc(3, 3);
    p = gsl_permutation_alloc(3);
    if ((m_temp == NULL) || (p == NULL)) {
        printf("Error: Failed allocating memory for a temporary matrix and/or permutation used in the debugging process. Exiting...\n");
        return 1;
    }
    for (int i = 0; i < par->n_su3; i++) {
        int signum;
        gsl_complex det;
        
        gsl_matrix_complex_memcpy(m_temp, su3[i]);
        gsl_linalg_complex_LU_decomp(m_temp, p, &signum);
        det = gsl_linalg_complex_LU_det(m_temp, signum);
        printf("DEBUG: det = %4.3f + %4.3fi\n", GSL_REAL(det), GSL_IMAG(det));
    }
    for (int i = 0; i < par->n_su3 / 2; i++) {
        gsl_complex z_temp;
        gsl_blas_zgemm(
            CblasNoTrans, 
            CblasNoTrans, 
            z_one, 
            su3[i], 
            su3[i + par->n_su3 / 2], 
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

    gsl_matrix_complex_free(m_temp);
    gsl_permutation_free(p);
    

    return 0;
}

void init_lattice(Par *par, gsl_matrix_complex **lattice, gsl_matrix_complex **su3) {
    /* initialize all independent links with random SU(3) matrices */
    for (int i = 0; i < par->L; i++) {
        for (int j = 0; j < par->L; j++) {
            for (int k = 0; k < par->L; k++) {
                for (int l = 0; l < par->L; l++) {
                    for (int dir = 0; dir < 4; dir++) {
                        gsl_matrix_complex_memcpy(
                            lattice[ind(i, j, k, l, dir, par->L)], 
                            su3[(int)(gsl_rng_uniform(par->ran_gen) * par->n_su3)]
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



int sim(Par *par, gsl_matrix_complex **lattice) {
    double results[2], 
           acceptance = 0.;
    gsl_matrix_complex **su3;
    
    /* allocate memory for a specified number of random SU(3) matrices */
    su3 = malloc(par->n_su3 * sizeof(gsl_matrix_complex *));
    if (su3 == NULL) {
        printf("Error: Failed allocating memory for array of SU(3)-matrix addresses. Exiting...\n");
        return 1;
    }
    for (int i = 0; i < par->n_su3; i++) {
        su3[i] = gsl_matrix_complex_alloc(3, 3);
        if (su3[i] == NULL) {
            printf("Error: failed allocating memory for su3 matrices. Exiting...\n");
            for (int j = 0; j < i; j++) 
                gsl_matrix_complex_free(su3[j]);
            free(su3);
            return 1;
        }
    }
    
    printf(
        "Starting Lattice-QCD simulation with parameters: L=%d, seed=%ld, eps=%4.3f, beta=%4.3f\n", 
        par->L, 
        par->seed,
        par->eps,
        par->beta
    );

    /* generate the SU(3) matrices */
    printf("Generating %d random SU(3)-matrices...", par->n_su3);
    if (init_su3(par, su3)) {
        for (int i = 0; i < par->n_su3; i++)
            gsl_matrix_complex_free(su3[i]);
         return 1;
    }
    printf("\n");

    /* initialize all links in the lattice */
    printf("Initializing %d x %d x %d x %d lattice...", par->L, par->L, par->L, par->L);
    init_lattice(par, lattice, su3);
    printf("\n");
    
    /* thermalize the lattice */
    printf("Thermalizing (%d sweeps)...", par->n_therm);
    for (int i = 0; i < par->n_therm; i++) {
        if(update_lattice(par, lattice, su3, &acceptance)) {
            for (int j = 0; j < par->n_su3; j++) 
                gsl_matrix_complex_free(su3[j]);
            free(su3);
            return 1;
        }
    }
    printf("\n");
    acceptance /= (double)(par->n_therm);
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
            if(update_lattice(par, lattice, su3, &acceptance)) {
                for (int j = 0; j < par->n_su3; j++) 
                    gsl_matrix_complex_free(su3[j]);
                free(su3);
                return 1;
            }
        }
        if (measure(par, results, lattice)) {
            for (int j = 0; j < par->n_su3; j++) 
                gsl_matrix_complex_free(su3[j]);
            free(su3);
            return 1;
        }
        printf("%3.2f\t%3.2f\n", results[0], results[1]);
    }
    acceptance /= (double)(par->n_configs * par->n_corr);
    printf("Acceptance: %3.2f\n", acceptance);
    
    for (int i = 0; i < par->n_su3; i++)
        gsl_matrix_complex_free(su3[i]);
    free(su3);

    return 0;
}


int read_args(Par *par, char *arg) {
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
        if (par->seed) gsl_rng_set(par->ran_gen, par->seed);
        else gsl_rng_set(par->ran_gen, (long)time(NULL));
        
        /* running the simulation with already given parameters */
        success = !sim(par, lattice);
        
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
            lattice[i] = gsl_matrix_complex_alloc(3, 3);
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
    
    if (!strcmp(arg, "n_su3")) {
        par->n_su3 = strtol(s, NULL, 0);
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

/*    if (!strcmp(arg, "gen_type")) {
        par->gen_type = strtol(s, NULL, 0);
        return 0;
    } */

    fprintf(stderr, "No such variable name: '%s'. Exiting...\n", arg);
    return 1;
}


int main(int argc, char *argv[])
{
    Par *par = NULL;
    /* allocate memory for simulation parameters */
     par = malloc(sizeof(Par));
     if (par == NULL) {
        printf("Error: Failed allocating memory for simulation parameters. Exiting...\n");
        exit(EXIT_FAILURE);
     }

    
    /* initialize simulation parameters with statndard values */
    par->L = 0;
    par->seed = 0;
    par->n_configs = 10;
    par->gen_type = gsl_rng_ranlxs0;
    par->n_su3 = 100;
    par->eps = 1.;
    par->n_hits = 10;
    par->n_therm = 500;
    par->n_corr = 50;
    par->beta = 0.;
  
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
