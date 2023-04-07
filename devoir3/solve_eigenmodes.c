#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmshc.h>
#include "math.h"
#include "functions/headers/power_method.h"
#include "functions/headers/matrix.h"
#include "functions/headers/elasticity.h"
#include "functions/headers/lu.h"
#include "functions/headers/rayleigh_quotient.h"


const double PI = 3.14159265359;

void set_vector(double *vec, size_t size) {
    double sqrt_size = sqrt(size);
    for (int i = 0; i < size; i++) {
        vec[i] = 1. / sqrt_size;
    }
}

int main (int argc, char *argv[]){

    if (argc < 5){
        printf("Usage: \n"
                "./solve_deformation <geo_file.geo> <meshSizeFactor> <k> <out>\n" 
                "---------------------------- \n\n"
                "- Use this .geo file : \n"
                "      square.geo or tuningFork.geo\n"
                "- meshSizeFactor sets a size factor on the mesh size; for the square: \n "
                "      0. < meshSizeFactor < 1.\n" 
                "- k is the number of frequencies to compute. \n "
                "- out is the output file to write the frequencies. \n "
        "\n");
            return -1;
    } 
    int k = atoi(argv[3]);

    // Define physical constants
    double E = 0.7e11;  // Young's modulus for Aluminum
    double nu = 0.3;    // Poisson coefficient
    double rho = 3000;  // Density of Aluminum


    // Initialize Gmsh and load geometry
    int ierr;
    gmshInitialize(argc, argv, 0, 0, &ierr);
    gmshOpen(argv[1], &ierr);
        
    // Set mesh size factor and generate mesh
    double meshSizeFactor;
    sscanf(argv[2],"%lf",&meshSizeFactor);
    gmshOptionSetNumber("Mesh.MeshSizeFactor", meshSizeFactor, &ierr);
    gmshModelMeshGenerate(2, &ierr);

    // Assemble the 2 matrices of the linear elasticity problem: 
    // M is the mass matrix
    // K is the stiffness matrix
    // boundary_nodes contains the index of the boundary nodes inside the matrix
    // coordinates contains the coordinates of the nodes, useful for your permutation :)
    Matrix *K, *M;
    double *coord;
    size_t* boundary_nodes;
    size_t n_boundary_nodes;
    assemble_system(&K, &M, &boundary_nodes, &n_boundary_nodes, &coord, E, nu, rho);
    size_t matrix_size = K->m;

    // 1. Remove lines from matrix that correspond to boundary nodes
    reduce_matrix(&K, &M, boundary_nodes, n_boundary_nodes);

    // print_matrix(M);
    // K^-1
    Matrix *K_inv = inverse_matrix_permute(K);
    // Matrix *K_inv = inverse_matrix(K);
    free_matrix(K);
    K = NULL;

    // K^-1 M
    Matrix *KM = mult_matrix(K_inv, M);
    free_matrix(K_inv);
    free_matrix(M);
    K_inv = M = NULL;

    // 2. Power iteration to find the largest eigenvalues
    
    double *vec = malloc(sizeof(double)*KM->m);
    double *vec_zeros = malloc(sizeof(double)*matrix_size);
    double eigen_values[k];

    // printf("Start to search for eigenvectors\n");
    for (int it = 0; it < k; it++) {
        set_vector(vec, KM->m);
        double lambda = 0.;
        power_iteration(KM, &lambda, vec, 1e-20);
        eigen_values[it] = lambda;
        
        // Deflation matrix A = A - lambda * uu^t
        for (int i = 0; i < KM->m; i++) {
            for (int j = 0; j < KM->n; j++) {
                KM->a[i][j] -= lambda * vec[i] * vec[j];
            }
        }

        // Recreate base vector
        int current = 0;
        for (int i = 0; i < matrix_size; i++) {
            int skip = 0;
            for (int j = 0; j < n_boundary_nodes; j++) {
                if (boundary_nodes[j] == i/2) {
                    skip = 1; break;
                }
            }
            if (skip) {
                vec_zeros[i] = 0.;
            }
            else {
                vec_zeros[i] = vec[current++];
            }
        }
        // Add it to gmsh
        visualize_in_gmsh(vec_zeros, matrix_size/2);
    }
   
    // gmshFltkRun(&ierr);
    free(vec);
    free(vec_zeros);
    free_matrix(KM);
    free(boundary_nodes);
    free(coord);

    // 3. Write frequencies to file out
    FILE *f = fopen(argv[4], "w");
    if (f == NULL) {
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < k-1; i++) {
        fprintf(f, "%.20f ", 1/sqrt(eigen_values[i])/2/PI);
    }
    if (k) {
        fprintf(f, "%.20f", 1/sqrt(eigen_values[k-1])/2/PI);
    }
    fclose(f);
    
    return 0;
}
