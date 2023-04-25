#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmshc.h>
#include "matrix.h"
#include "elasticity.h"
#include "math.h"
#include "lu.h"
#include "design.h"
#include "eigen.h"
#include "rmck.h"
#include <lapacke.h> 
#include <cblas.h>

 #define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
 #define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define MESHSIZE 0.3
#define DELTA 0.00000001
#define LR .000001
#define FREQ 1567.98

double E = 0.7e11;  // Young's modulus for Aluminum
double nu = 0.3;    // Poisson coefficient
double rho = 3000;  //

void compute_gradient(double r1, double r2, double e, double l, double *gradient) {
    designTuningFork(r1, r2, e, l, MESHSIZE, NULL);
    Matrix *K, *M;
    size_t* boundary_nodes;
    size_t n_boundary_nodes;
    double * coord;
    assemble_system(&K, &M, &coord, &boundary_nodes, &n_boundary_nodes, E, nu, rho);

    // Remove lines from matrix that are boundary
    Matrix *K_new;
    Matrix *M_new;
    remove_bnd_lines(K, M, boundary_nodes, n_boundary_nodes, &K_new, &M_new, NULL);
    free_matrix(K); free_matrix(M); K = M = NULL;
    
    inverse_matrix_permute(K_new, M_new);
    Matrix *A = M_new;
    free_matrix(K_new); K_new = NULL;

    double v[A->m];
    double lambda = power_iteration(A, v);

    free_matrix(M_new); M_new = NULL;
    double factor[4] = {0.,0.,0.,0.};
    double freq = 1./(2*M_PI*sqrt(lambda));
    printf("%lf %lf\n", freq, FREQ);
    for (int i = 0; i < 4; i++) {
        factor[i] = DELTA;
        Matrix *K, *M;
        // int err;
        // gmshInitialize(0, NULL, 0, 0, &err);
        designTuningFork(r1+factor[0], r2+factor[1], e+factor[2], l+factor[3], MESHSIZE, NULL);
        assemble_system(&K, &M, &coord, &boundary_nodes, &n_boundary_nodes, E, nu, rho);

        // Remove lines from matrix that are boundary
        Matrix *K_new;
        Matrix *M_new;
        remove_bnd_lines(K, M, boundary_nodes, n_boundary_nodes, &K_new, &M_new, NULL);
        free_matrix(K); free_matrix(M); K = M = NULL;
        
        inverse_matrix_permute(K_new, M_new);
        Matrix *A = M_new;
        free_matrix(K_new); K_new = NULL;

        double v[A->m];
        double lambda_eps = power_iteration(A, v);
        gradient[i] = 2.*(freq -  FREQ) * (lambda - lambda_eps) / DELTA;
        free_matrix(M_new); M_new = NULL;
        factor[i] = 0.;
    }
}

int main (int argc, char *argv[]) {
    int ierr;
    gmshInitialize(argc, argv, 0, 0, &ierr);
    double params[4] = {0.001851, 0.022364, 0.038419, 0.076507};
    printf("%lf, %lf, %lf, %lf\n", params[0], params[1], params[2], params[3]);
    for (int it = 0; it < 1000; it++) {
        double gradient[4];
        compute_gradient(params[0], params[1], params[2], params[3], gradient);
        for (int i = 0; i < 4; i++) {
            params[i] = params[i] - LR * gradient[i];
        }
        // printf("%lf, %lf, %lf, %lf\n", gradient[0], gradient[1], gradient[2], gradient[3]);
        printf("%lf, %lf, %lf, %lf\n", params[0], params[1], params[2], params[3]);
    }
    printf("%lf, %lf, %lf, %lf\n", params[0], params[1], params[2], params[3]);
    return 0;
}
