#include <stdio.h>
#include <stdlib.h>
#include <gmshc.h>
#include "matrix.h"
#include "elasticity.h"
#include "math.h"
#include "design.h"
#include "eigen.h"
#include <cblas.h>

#define max(a,b) \
   	({ __typeof__ (a) _a = (a); \
		__typeof__ (b) _b = (b); \
	_a > _b ? _a : _b; })

#define min(a,b) \
   	({ __typeof__ (a) _a = (a); \
	   __typeof__ (b) _b = (b); \
		_a < _b ? _a : _b; })


int main (int argc, char *argv[]) {

  	if (argc < 2) {
		printf("Usage: \n"
		"./project <k> <out>\n" 
		"---------------------------- \n\n"
		"- k is the number of frequencies to compute. \n "
		"- out is the output file to write the frequencies. \n "
		"\n");
		return -1;
  	}

	// Number of vibration modes to find
	size_t const k = atoi(argv[1]);

	// Define physical constantssss
	double const E = 0.7e11;  // Young's modulus for Aluminum
	double const nu = 0.3;    // Poisson coefficient
	double const rho = 3000;  // Density of Aluminum

	// Initialize Gmsh and create geometry
	int ierr;
	gmshInitialize(argc, argv, 0, 0, &ierr);
 
	// Create geometry.
	double const d[6] = {0.05929252977121609, 0.021409268854274743, 11e-3, 11e-3, 11e-3, 11e-3};
	double const dec[6] = {.0, 0.0, 0, 0, 0, 0};
	double const h[6] = {0.06581162483843764, 0.021227934873379272, 11e-3, 12e-3, 12e-3, 12e-3};
	double const l[6] = {0.040410, 0.056977, 50e-3, 82e-3, 82e-3, 82e-3};
	designTuningForkNLayer(0.06595566225750074, 0.009164195816166532, d, dec, h, l, 2, 0.3, NULL);
  
	// Assemble the 2 matrices of the linear elasticity problem: 
	// M is the mass matrix && K is the stiffness matrix
	Matrix* K;
	Matrix* M;
	size_t* boundary_nodes;
	size_t n_boundary_nodes;
	double* coord;
	assemble_system(&K, &M, &coord, &boundary_nodes, &n_boundary_nodes, E, nu, rho);
	free(coord);
	coord = NULL;
	size_t const m = K->m;

	// Remove lines from matrix that are boundary
	Matrix* K_new;
	Matrix* M_new;
	remove_bnd_lines(K, M, boundary_nodes, n_boundary_nodes, &K_new, &M_new, NULL);
  	free_matrix(K);
	free_matrix(M);
	// Avoid dangling pointers
	K = M = NULL;
  
	inverse_matrix_permute(K_new, M_new);
	Matrix* A = M_new;
	free_matrix(K_new);
	K_new = NULL;

	// Power iteration + deflation to find k largest eigenvalues
	double* const v = malloc(A->m * sizeof(double));
	double lambda;
	double freq;
	FILE * file = NULL;
	if (argc >= 3)
		file = fopen(argv[2], "w"); // open file to write frequencies
	
	double* const freq_vector = calloc(m, sizeof(double));
	for(size_t ki = 0; ki < k; ki++) {
		lambda = power_iteration(A, v);
		freq = 1. / (2 * M_PI *sqrt(lambda));
		
		if (file != NULL)
			fprintf(file, "%.9lf ", freq);

		printf("lambda = %.9e, f = %.3lf\n", lambda, freq);

		// Deflate matrix
		cblas_dger(CblasRowMajor, A->m, A->m, -lambda, v, 1, v, 1, A->data, A->m);
		// for(size_t i = 0; i < A->m; i++)
			// for(size_t j = 0; j < A->n; j++)
				// A->a[i][j] -= lambda * v[i] * v[j];

		// Put appropriate BC and plot
		size_t iv = 0;
		size_t i_bnd = 0; 
		for(size_t i = 0; i < m/2; i++) {
			if(i_bnd < n_boundary_nodes && i == boundary_nodes[i_bnd]) {
				i_bnd++;
				continue;
			}
			freq_vector[2 * i + 0] = v[2 * iv + 0];
			freq_vector[2 * i + 1] = v[2 * iv + 1];
			iv++;
		}
		visualize_in_gmsh(freq_vector, m / 2);
	}

	if (file != NULL)
  		fclose(file);

  	// gmshFltkRun(&ierr);

	// Don't need to free the os is faster...
	// free(freq_vector);
	// free_matrix(M_new);
	// free(boundary_nodes);
	// free(v);
  	return 0;
}