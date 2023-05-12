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

#define MESHSIZE 0.2
#define TARGET_FREQ 1567.98

double E = 0.7e11;  // Young's modulus for Aluminum
double nu = 0.3;    // Poisson coefficient
double rho = 3000;  // Mass density

// Fixed parameter
double r1 = 6e-3;
double r2 = 11e-3;
double e = 38e-3;

Matrix *compute_matrix(double r1, double r2, double e, double l)  {
	Matrix* K;
	Matrix* M;
	double* coord;
	size_t* boundary_nodes;
	size_t n_boundary_nodes;

	designTuningForkSymmetric(r1, r2, e, l, MESHSIZE, NULL);
	assemble_system(&K, &M, &coord, &boundary_nodes, &n_boundary_nodes, E, nu, rho);

	// Remove lines from matrix that are boundary
	Matrix* K_new;
	Matrix* M_new;
	remove_bnd_lines(K, M, boundary_nodes, n_boundary_nodes, &K_new, &M_new, NULL);
	free_matrix(K); free_matrix(M); K = M = NULL;
	
	inverse_matrix_permute(K_new, M_new);
	free_matrix(K_new); K_new = NULL;
	return M_new;
}

int main (int argc, char *argv[]) {
  // Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(argc, argv, 0, 0, &ierr);

  double l = 10e-3;
  double r = 1e-1;
  size_t it = 0;
  while (r - l > 1e-15) {
	double const mid = l + (r - l) / 2;

	Matrix* A = compute_matrix(r1, r2, e, mid);
	double* v = malloc(A->m * sizeof(double));
	double const lambda = power_iteration(A, v);
	double const freq = 1. /(2 * M_PI * sqrt(lambda));


	// printf("f1 = %.3lf\n", freq);
	if (fabs(freq - TARGET_FREQ) < 1.) {
		printf("Convergence : %ld\n", it);
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

  printf("To get f=%lf, r1=%lf r2=%lf e=%lf l=%lf\n", TARGET_FREQ, r1, r2, e, l);

  return 0;
}
