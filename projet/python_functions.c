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

double* get_k_freq_tuning_fork(size_t const k, double const r1, double const r2, double const e, double const l) {
    // Initialize Gmsh and create geometry
    int ierr;
    // char* arg1 = "-v";
    // char* arg2 = "0";
    // char* argv[2] = {arg1, arg2};
    gmshInitialize(0, NULL, 0, 0, &ierr);

    designTuningForkSymmetric(r1, r2, e, l, 0.3, NULL);
    Matrix *K, *M;
    size_t* boundary_nodes;
    size_t n_boundary_nodes;
    double* coord;
    assemble_system(&K, &M, &coord, &boundary_nodes, &n_boundary_nodes, E, nu, rho);

    // Remove lines from matrix that are boundary
    Matrix *K_new;
    Matrix *M_new;
    remove_bnd_lines(K, M, boundary_nodes, n_boundary_nodes, &K_new, &M_new, NULL);
    // free_matrix(K); free_matrix(M); K = M = NULL;
    
    inverse_matrix_permute(K_new, M_new);
    Matrix *A = M_new;
    free_matrix(K_new); K_new = NULL;
    

    // Power iteration + deflation to find k largest eigenvalues
    double * v = malloc(A->m * sizeof(double));
    double lambda;
    double* freqs = malloc(k * sizeof(double));
    for(int ki = 0; ki < k; ki++) {
        lambda = power_iteration(A, v);
        freqs[ki] = 1./(2 * M_PI * sqrt(lambda));

        // Deflate matrix
        for(int i = 0; i < A->m; i++)
            for(int j = 0; j < A->n; j++)
                A->a[i][j] -= lambda * v[i] * v[j];
    }

    free_matrix(M_new);
    free(boundary_nodes);
    return freqs;
}