#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmshc.h>
#include "matrix.h"
#include "elasticity.h"
#include "math.h"
#include "design.h"
#include "eigen.h"

#define MESHSIZE 0.2
#define DELTA 0.000001
#define LR .1
#define FREQ 1567.98

double E = 0.7e11;  // Young's modulus for Aluminum
double nu = 0.3;    // Poisson coefficient
double rho = 3000;

int compute_gradient(double const r1, double const r2, double const e, double const l, double *gradient) {
    double factor[4] = {0.,0.,0.,0.};
    designTuningForkSymmetric(r1 + factor[0], r2 + factor[1], e + factor[2], l + factor[3], MESHSIZE, NULL);
    Matrix *A = compute_matrix_km(E, nu, rho);

    double v[A->m];
    double lambda = power_iteration(A, v);
    free_matrix(A); A = NULL;
    
    double freq = 1. / (2 * M_PI * sqrt(lambda));
    if (fabs(freq - FREQ) < 1.) {
        printf("%lf %lf\n", freq, FREQ);
        return 1;
    }
    for (int i = 0; i < 4; i++) {
        double lambda[2];
        for (int j = 0; j < 2; j++) {
            double fact = j;
            if (j == 2)
                fact = -1.;
            factor[i] = DELTA * fact;

            designTuningForkSymmetric(r1 + factor[0], r2 + factor[1], e + factor[2], l + factor[3], MESHSIZE, NULL);
            A = compute_matrix_km(E, nu, rho);

            double v[A->m];
            lambda[j] = power_iteration(A, v);
            free_matrix(A); A = NULL;
            factor[i] = 0.;
        }
        gradient[i] = (freq -  FREQ) * (lambda[0] - lambda[1]) / 2 / DELTA;
    }
    return 0;
}

int main (int argc, char *argv[]) {
    int ierr;
    gmshInitialize(argc, argv, 0, 0, &ierr);
    double params[4] = {6e-3, 11e-3, 38e-3, 82e-3};
    // printf("%lf, %lf, %lf, %lf\n", params[0], params[1], params[2], params[3]);
    for (int it = 0; it < 10000; it++) {
        double gradient[4];
        int res = compute_gradient(params[0], params[1], params[2], params[3], gradient);
        if (res) {
            printf("Converged: %d iterations\n", it);
            break;
        }
        for (int i = 0; i < 4; i++) {
            params[i] = fmin(fmax(params[i] - LR * gradient[i], 1e-3), 1.);
        }

        // printf("%lf, %lf, %lf, %lf\n", gradient[0], gradient[1], gradient[2], gradient[3]);
        // printf("%lf, %lf, %lf, %lf\n", params[0], params[1], params[2], params[3]);
    }
    printf("%.20f, %.20f, %.20lf, %.20lf\n", params[0], params[1], params[2], params[3]);
    return 0;
}
