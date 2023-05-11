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

#define MESHSIZE 0.2
#define DELTA 0.000001
#define LR .1
#define FREQ 1567.98

double E = 0.7e11;  // Young's modulus for Aluminum
double nu = 0.3;    // Poisson coefficient
double rho = 3000;  //

Matrix *compute_matrix(double r1, double r2, double e, double l, double factor[4])  {
    Matrix *K, *M;
    double *coord; size_t *boundary_nodes; size_t n_boundary_nodes;

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
    return A;
}

void compute_gradient(double r1, double r2, double e, double l, double *gradient) {
    double factor[4] = {0.,0.,0.,0.};
    Matrix *A = compute_matrix(r1, r2, e, l, factor);

    double v[A->m];
    double lambda[4];
    for (int k = 0; k < 4; k++) {
        lambda[k] = power_iteration(A, v);
        for(int i = 0; i < A->m; i++)
            for(int j = 0; j < A->n; j++)
                A->a[i][j] -= lambda[k] * v[i] * v[j];

    }
    free_matrix(A); A = NULL;
    
    double freq = 1./(2*M_PI*sqrt(lambda[1]));
    double freq2 = 1./(2*M_PI*sqrt(lambda[3]));
    printf("%lf %lf %lf %lf\n", freq, FREQ, freq2, 2*FREQ);
    for (int i = 0; i < 4; i++) {
        double lambda[2][4];
        for (int j = 0; j < 2; j++) {
            double fact = j;
            if (j == 2) fact = -1.;
            factor[i] = DELTA*fact;
            A = compute_matrix(r1, r2, e, l, factor);
            double v[A->m];
            for (int k = 0; k < 4; k++) {
                lambda[j][k] = power_iteration(A, v);
                for(int ii = 0; ii < A->m; ii++)
                    for(int jj = 0; jj < A->n; jj++)
                        A->a[ii][jj] -= lambda[j][k] * v[ii] * v[jj];
            }
            free_matrix(A); A = NULL;
            factor[i] = 0.;
        }
        // double f1 = 1./(2*M_PI*sqrt(lambda[1][1]));
        // double f2 = 1./(2*M_PI*sqrt(lambda[0][1]));
        // gradient[i] = (fabs(f1 - FREQ) - fabs(f2 - FREQ)) / 2 / DELTA;
        gradient[i] = (freq -  FREQ) * (lambda[1][1] - lambda[0][1]) / 2 / DELTA; // + (freq2-2*FREQ) * (lambda[1][3] - lambda[0][3]) / 2 / DELTA;
    }
}

int main (int argc, char *argv[]) {
    int ierr;
    gmshInitialize(argc, argv, 0, 0, &ierr);
    double params[4] = {0.004112, 0.022369, 0.038419, 0.076508};
    printf("%lf, %lf, %lf, %lf\n", params[0], params[1], params[2], params[3]);
    for (int it = 0; it < 1000; it++) {
        double gradient[4];
        compute_gradient(params[0], params[1], params[2], params[3], gradient);
        for (int i = 0; i < 4; i++) {
            
            params[i] = params[i] - LR * gradient[i];
        }
        printf("%lf, %lf, %lf, %lf\n", gradient[0], gradient[1], gradient[2], gradient[3]);
        printf("%lf, %lf, %lf, %lf\n", params[0], params[1], params[2], params[3]);
    }
    printf("%lf, %lf, %lf, %lf\n", params[0], params[1], params[2], params[3]);
    return 0;
}
