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

int main (int argc, char *argv[]){

    if (argc < 4){
        printf("Usage: \n"
                "./solve_deformation <geo_file.geo> <meshSizeFactor> <k> <out>\n" 
                "---------------------------- \n\n"
                "- Use this .geo file : \n"
                "      square.geo or tuningFork.geo\n"
                "- meshSizeFactor sets a size factor on the mesh size; for the square: \n "
                "      0. < meshSizeFactor < 1.\n" 
                "- p 1 or 0, 1 mean permuted matrix for inversion 0 not \n "
        "\n");
            return -1;
    } 
    int p = atoi(argv[3]);

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
    // K^-1
    Matrix *K_inv;
    if (p) {
        K_inv = inverse_matrix_permute(K);
    }
    else {
        K_inv = inverse_matrix(K);
    }

    // OS will free by itself faster
    // free_matrix(K);
    // free_matrix(K_inv);
    // free_matrix(M);
    // gmshFltkRun(&ierr);
    // free(boundary_nodes);
    // free(coord);


}
