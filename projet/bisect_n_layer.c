#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmshc.h>
#include "matrix.h"
#include "elasticity.h"
#include "math.h"
#include "design.h"
#include "eigen.h"
#include <cblas.h>

#define MESHSIZE 0.2
#define TARGET_FREQ 1567.98

double E = 0.7e11;  // Young's modulus for Aluminum
double nu = 0.3;    // Poisson coefficient
double rho = 3000;  // Mass density

// Fixed parameter
double const lh = 0.0659556622575007;
double const wh = 0.009164195816166532;
double const d[6] = {0.05929252977121609, 0.021409268854274743, 11e-3, 11e-3, 11e-3, 11e-3};
double const dec[6] = {.0, 0.0, 0, 0, 0, 0};
double const h[6] = {0.06581162483843764, 0.021227934873379272, 11e-3, 12e-3, 12e-3, 12e-3};
double length[6] = {50e-3, 50e-3, 50e-3, 82e-3, 82e-3, 82e-3};

int main (int argc, char *argv[]) {
  // Initialize Gmsh and create geometry
    int ierr;
    gmshInitialize(argc, argv, 0, 0, &ierr);
    for (int i = 1; i >= 0; i--) {
        double l = 1e-3;
        double r = 1e-1;
        size_t it = 0;
        while (r - l > 1e-16) {
            double const mid = l + (r - l) / 2;
            length[i] = mid;
            designTuningForkSymmetricNLayer(lh, wh, d, dec, h, length, 2, MESHSIZE, NULL);
            Matrix* A = compute_matrix_km(E, nu, rho);

            double* v = malloc(A->m * sizeof(double));
            double lambda = power_iteration(A, v);
            for (int j = 0; j < i; j++) {
                cblas_dger(CblasRowMajor, A->m, A->m, -lambda, v, 1, v, 1, A->data, A->m);
                lambda = power_iteration(A, v);
            }
            double const freq = 1. /(2 * M_PI * sqrt(lambda));
            if (it % 100 == 0)
                printf("Current freq %d is %lf\n", i + 1, freq);

            // printf("f1 = %.3lf\n", freq);
            if (fabs(freq - (i * 1) * TARGET_FREQ) < 1.) {
                printf("Convergence : %ld\n", it + 1);
                break;
            }

            free(v);
            free_matrix(A);

            if (freq > TARGET_FREQ * (i + 1))
                l = mid;
            else
                r = mid;
            it++;
        }
    }
    printf("%lf %lf\n", length[0], length[1]);
    // printf("To get f=%lf, r1=%lf r2=%lf e=%lf l=%lf\n", TARGET_FREQ, r1, r2, e, l);

    return 0;
}
