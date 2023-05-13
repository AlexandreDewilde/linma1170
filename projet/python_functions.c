#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmshc.h>
#include "matrix.h"
#include "elasticity.h"
#include "math.h"
#include "design.h"
#include "eigen.h"

static double E = 0.7e11;  // Young's modulus for Aluminum
static double nu = 0.3;    // Poisson coefficient
static double rho = 3000;  // Density of Aluminum

void init_gmsh() {
    int ierr;
    char* args[3] = {"", "-v", "0"};
    gmshInitialize(3, args, 0, 0, &ierr);
}

double* get_k_freq_tuning_fork(size_t const k, double const r1, double const r2, double const e, double const l) {
    designTuningForkSymmetric(r1, r2, e, l, 0.3, NULL);
    Matrix* A = compute_matrix_km(E, nu, rho);  

    // Power iteration + deflation to find k largest eigenvalues
    double *v = malloc(A->m * sizeof(double));
    double* freqs = malloc(k * sizeof(double));
    for(int ki = 0; ki < k; ki++) {
        double lambda = power_iteration(A, v);
        freqs[ki] = 1. / (2 * M_PI * sqrt(lambda));

        // Deflate matrix
        for(int i = 0; i < A->m; i++)
            for(int j = 0; j < A->n; j++)
                A->a[i][j] -= lambda * v[i] * v[j];
    }

    free_matrix(A);
    free(v);
    A = v = NULL;
    return freqs;
}