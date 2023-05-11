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

  if (argc < 2){
    printf("Usage: \n"
      "./project <k> <out>\n" 
      "---------------------------- \n\n"
      "- k is the number of frequencies to compute. \n "
      "- out is the output file to write the frequencies. \n "
      "\n");
    return -1;
  }

  // Define physical constants
  double E = 0.7e11;  // Young's modulus for Aluminum
  double nu = 0.3;    // Poisson coefficient
  double rho = 3000;  // Density of Aluminum

  // Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(argc, argv, 0, 0, &ierr);

  double r1 = atof(argv[2]);
  double r2 = atof(argv[3]);
  double e = atof(argv[4]);
  double l = atof(argv[5]);
  designTuningFork(r1, r2, e, l, 0.3, NULL);
  
  // Number of vibration modes to find
  int k = atoi(argv[1]);

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
  printf("%lf\n", fabs(freq[1] - TARGET_FREQ));

  // free_matrix (K);
  // free_matrix (M);
  // free_matrix (K_new);
  // free_matrix (M_new);
  // free(boundary_nodes);

  return 0;
}