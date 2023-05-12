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

#define TARGET_FREQ 1567.98

int main (int argc, char *argv[]) {

  // Define physical constants
  double E = 0.7e11;  // Young's modulus for Aluminum
  double nu = 0.3;    // Poisson coefficient
  double rho = 3000;  // Density of Aluminum

  // Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(argc, argv, 0, 0, &ierr);

  double e = atof(argv[1]);
  double d1 = atof(argv[2]);
  double d2 = atof(argv[3]);
  double h1 = atof(argv[4]);
  double h2 = atof(argv[5]);
  double l = atof(argv[6]);
  double l2 = atof(argv[7]);
  designTuningForkSymmetric2Layer(e, d1, d2, h1, h2, l, l2, 0.2, NULL);
  
  // Number of vibration modes to find
  int k = 2;

  // Assemble the 2 matrices of the linear elasticity problem: 
  // M is the mass matrix
  // K is the stiffness matrix
  Matrix *K, *M;
  size_t* boundary_nodes;
  size_t n_boundary_nodes;
  double * coord;
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
  double lambda, freq[4];
  for(int ki = 0; ki < k; ki++) {
    lambda = power_iteration(A, v);
    freq[ki] = 1./(2*M_PI*sqrt(lambda));
    // fprintf(file, "%.9lf ", freq);

    // printf("%.3lf\n", freq);
    // Deflate matrix
    for(int i = 0; i < A->m; i++)
      for(int j = 0; j < A->n; j++)
        A->a[i][j] -= lambda * v[i] * v[j];
  }
  printf("%lf\n", fabs(freq[1] - 2. * TARGET_FREQ) + fabs(freq[0] - TARGET_FREQ));

  // free_matrix (K);
  // free_matrix (M);
  // free_matrix (K_new);
  // free_matrix (M_new);
  // free(boundary_nodes);

  return 0;
}