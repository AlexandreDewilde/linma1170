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
#define FREQ 1567.98

double E = 0.7e11;  // Young's modulus for Aluminum
double nu = 0.3;    // Poisson coefficient
double rho = 3000;  //

Matrix *compute_matrix(double r1, double r2, double e, double l)  {
	Matrix *K, *M;
	double *coord; size_t *boundary_nodes; size_t n_boundary_nodes;

	designTuningFork(r1, r2, e, l, MESHSIZE, NULL);
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

int main (int argc, char *argv[]) {
  // Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(argc, argv, 0, 0, &ierr);

  double l = 10e-3;
  double r = 1e-1;
  int it = 0;
  while (r - l > 1e-15) {
	double mid = l + (r - l) / 2;
	double freqs[4];
	Matrix* A = compute_matrix(6e-3, 11e-3, 38e-3, mid);

	// Power iteration + deflation to find k largest eigenvalues
	double *v = malloc(A->m * sizeof(double));
	double lambda, freq;
	for(int ki = 0; ki < 4; ki++) {
		lambda = power_iteration(A, v);
		freq = 1./(2*M_PI*sqrt(lambda));
		freqs[ki] = freq;
		
		// Deflate matrix
		for(int i = 0; i < A->m; i++)
			for(int j = 0; j < A->n; j++)
				A->a[i][j] -= lambda * v[i] * v[j];
	}
	printf("f1 = %.3lf, f2 = %.3lf\n", freqs[1], freqs[3]);
	if (fabs(freqs[1] - FREQ) < 0.5) {
		printf("Convergence : %d\n", it);
	}
	free(v);
	free_matrix(A);
	double f = freqs[1] - FREQ;
	if (f > 0)
		l = mid;
	else
		r = mid;
	it++;
  }
  printf("%lf\n", l);

  return 0;
}
