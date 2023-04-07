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
    
    FILE *f = fopen("chart/hybrid.csv", "w");
    // Set mesh size factor and generate mesh
    
    srand(time(NULL));
    fprintf(f, "meshsize;power_iteration;hybrid_pi;hybrid_rqi\n");
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
        double *vec = malloc(sizeof(double)*KM->m);
        for (int i = 0; i < KM->m; i++) {
            vec[i] = 1. / KM->m;
        }
        double lambda = 0.;
        int it_pi = power_iteration(KM, &lambda, vec, 1e-10);
        for (int i = 0; i < KM->m; i++) {
            vec[i] = 1. / KM->m;
        }
        // printf("Hybrid method\n");
        int it_hybrid_pi = power_iteration(KM, &lambda, vec, 1.);
        int it_hybrid_rqi = rayleigh_quotient_iteration(KM, &lambda, vec, 1e-10);
        fprintf(f, "%lf;%d;%d;%d\n", meshSizeFactor, it_pi, it_hybrid_pi, it_hybrid_rqi);
        free_matrix(KM);
        free(vec);
    }
    fclose(f);
}
