#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmshc.h>
#include "math.h"
#include "functions/headers/matrix.h"
#include "functions/headers/elasticity.h"
#include "functions/headers/power_method.h"
#include "functions/headers/rayleigh_quotient.h"

int main(int argc, char *argv[]) {
    // Define physical constants
    double E = 0.7e11;  // Young's modulus for Aluminum
    double nu = 0.3;    // Poisson coefficient
    double rho = 3000;  // Density of Aluminum


    // Initialize Gmsh and load geometry
    int ierr;
    gmshInitialize(argc, argv, 0, 0, &ierr);
    gmshOpen("tuningFork.geo", &ierr);
    
    FILE *f = fopen("chart/convergence.csv", "w");
    // Set mesh size factor and generate mesh
    
    srand(time(NULL));
    fprintf(f, "meshsize;power_iteration;rayleigh_quotient_iteration\n");
    for (int it = 0; it < 50; it+=5) {
        double meshSizeFactor = ((double) (it*0.01 + 0.4));
        
        gmshOptionSetNumber("Mesh.MeshSizeFactor", meshSizeFactor, &ierr);
        gmshModelMeshGenerate(2, &ierr);

        Matrix *K, *M;
        double *coord;
        size_t* boundary_nodes;
        size_t n_boundary_nodes;
        assemble_system(&K, &M, &boundary_nodes, &n_boundary_nodes, &coord, E, nu, rho);
        size_t matrix_size = K->m;

        char *boundary_bool = calloc(K->m/2, 1);
        for (int i = 0; i < n_boundary_nodes; ++i) {
            boundary_bool[boundary_nodes[i]] = 1;
        }

        reduce_matrix(&K, &M, boundary_bool, n_boundary_nodes);
        free(boundary_bool);

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
        double *vec = malloc(sizeof(double)*KM->m);
        for (int i = 0; i < KM->m; i++) {
            vec[i] = 1. / KM->m;
        }
        double lambda = 0.;
        power_iteration(KM, &lambda, vec, 1e-10);
        for (int i = 0; i < KM->m; i++) {
                if (!(rand() % 4))
                    vec[i] += (((double) (rand() % 100))+1.) / (rand() % 10 + 1.);
            }
        double *vec2 = malloc(sizeof(double)*KM->m); 
        memcpy(vec2, vec, KM->m*sizeof(double));
        int it_pi = power_iteration(KM, &lambda, vec, 1e-10);
        int it_rqi = rayleigh_quotient_iteration(KM, &lambda, vec2, 1e-10);
        fprintf(f, "%lf;%d;%d\n", meshSizeFactor, it_pi, it_rqi);
        
        free_matrix(KM);
        free(vec);
        free(vec2);
    }
    fclose(f);
}
