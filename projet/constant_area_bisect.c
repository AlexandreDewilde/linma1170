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
#define TARGET_FREQ 1567.98
#define TARGET_AREA 5e-4

double E = 0.7e11;  // Young's modulus for Aluminum
double nu = 0.3;    // Poisson coefficient
double rho = 3000;  // Mass density

// Fixed parameter
double r1 = 6e-3;
double r2 = 11e-3;
double e = 38e-3;

double compute_area(double const r1, double const r2, double const e, double const l) {
    return (r2 - r1) * e + 2 * (r2 - r1) * l + M_PI * (r2 * r2 - r1 * r1) / 4;
}

int main (int argc, char *argv[]) {
  // Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(argc, argv, 0, 0, &ierr);

  double l = 1e-3;
  double r = 1e-1;
  size_t it = 0;
  while (r - l > 1e-15) {
	double const mid = l + (r - l) / 2;

	designTuningForkSymmetric(r1, r2, e, mid, MESHSIZE, NULL);
	Matrix* A = compute_matrix_km(E, nu, rho);

	double* v = malloc(A->m * sizeof(double));
	double const lambda = power_iteration(A, v);
	double const freq = 1. /(2 * M_PI * sqrt(lambda));

	// printf("f1 = %.3lf\n", freq);
	if (fabs(freq - TARGET_FREQ) < 1.) {
		printf("Convergence for freq : %ld\n", it + 1);
		break;
	}

	free(v);
	free_matrix(A);

	if (freq > TARGET_FREQ)
		l = mid;
	else
		r = mid;
	it++;
  }
  double const new_length = l;
  l = 1e-3;
  r = 100.;
  it = 0;
  while (r - l > 1e-30) {
	double const mid = l + (r - l) / 2;
    double const area = compute_area(r1, r2, mid, new_length);
	// printf("f1 = %.3lf\n", freq);
	if (fabs(TARGET_AREA - area) < 1e-9) {
		printf("Convergence for area: %ld\n", it + 1);
		break;
	}

	if (area < TARGET_AREA)
		l = mid;
	else
		r = mid;
	it++;
  }


    printf("To get f=%lf with area=%lf, r1=%lf r2=%lf e=%lf l=%lf\n", TARGET_FREQ, TARGET_AREA, r1, r2, l, new_length);

  return 0;
}
