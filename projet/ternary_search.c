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
#define FREQ 1567.98

double E = 0.7e11;  // Young's modulus for Aluminum
double nu = 0.3;    // Poisson coefficient
double rho = 3000;

int main (int argc, char *argv[]) {
  // Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(argc, argv, 0, 0, &ierr);

  double l = 10e-3;
  double r = 1e-1;
  size_t it = 0;
  while (r - l > 1e-17) {
	double m1 = l + (r - l) / 3;
    double m2 = r - (r - l) / 3;
    double mids[2] = {m1, m2};
	double freqs[2];
    for (int it = 0; it < 2; it++) {
        designTuningForkSymmetric(6e-3, 11e-3, 38e-3, mids[it], MESHSIZE, NULL);
        Matrix* A = compute_matrix_km(E, nu, rho);

        // Power iteration + deflation to find k largest eigenvalues
        double* v = malloc(A->m * sizeof(double));
        double lambda = power_iteration(A, v);
        freqs[it] = 1. / (2 * M_PI * sqrt(lambda));
        
        free(v);
        free_matrix(A);
    }
    double f1 = (freqs[0] - FREQ) * (freqs[0] - FREQ);
    double f2 = (freqs[1] - FREQ) * (freqs[1] - FREQ);
    if (fabs(freqs[0] - FREQ) < 1.) {
        printf("Convergence : %ld\n", it + 1);
        break;
    }
    if (f1 > f2)
        l = mids[0];
    else
        r = mids[1];
    it++;
  }
  printf("Those parameters will produce the desired frequency r1=6e-3, r2=11e-3, e=28e-3, l=%lf\n", l);

  return 0;
}
