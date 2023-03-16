#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmshc.h>
#include <time.h>
#include "datastructures/headers/matrix.h"
#include "headers/elasticity.h"
#include "headers/lu.h"

// base on this SO post https://stackoverflow.com/questions/3557221/how-do-i-measure-time-in-c
double get_time()
{
    struct timespec now;
    timespec_get(&now, TIME_UTC);
    return now.tv_sec + now.tv_nsec*1e-9;
}

int main (int argc, char *argv[]) {

	if (argc != 3) {
		printf("Usage: \n"
			"./solve_deformation <geo_file.geo> <meshSizeFactor>\n" 
			"---------------------------- \n\n"
			"- Use this .geo file : \n"
			"      square.geo \n"
			"- meshSizeFactor sets a size factor on the mesh size; for the square: \n "
			"      0. < meshSizeFactor < 1.\n \n");
		return -1;
	}

	int ierr; // Gmsh error code

	// Initialize Gmsh and load geometry
	gmshInitialize(argc, argv, 0, 0, &ierr);
	gmshOpen(argv[1], &ierr);

	// Set mesh size factor and generate mesh
	double meshSizeFactor;
	sscanf(argv[2],"%lf",&meshSizeFactor);
	gmshOptionSetNumber("Mesh.MeshSizeFactor", meshSizeFactor, &ierr);
	gmshModelMeshGenerate(2, &ierr);

	// Setup linear system
	int n_nodes, n_triplets;
	size_t *gmsh_num;  // For final renumbering in Gmsh
	double *coord;     // Vector :  Coordinates of nodes (size 2*n_nodes)
	double *RHS;       // Vector : right-hand side of equation (size 2*n_nodes)
	Triplet* triplets; // Triplets of values to insert in matrix; each entry is a struct containing (i, j, val)
	assemble_system(&triplets, &n_triplets, &RHS, &coord, &gmsh_num, &n_nodes);

	// Build matrix
	int *perm = malloc(2*n_nodes*sizeof(int));
	for (int i = 0; i < 2*n_nodes; i++) perm[i] = i;
	compute_permutation(perm, coord, n_nodes, triplets, n_triplets);
	int *match_perm = malloc(2*n_nodes*sizeof(int));
	for (int i = 0; i < 2*n_nodes; i++) {
		match_perm[perm[i]] = i;
	}
	Matrix * K = allocate_matrix(2*n_nodes, 2*n_nodes);
	memset(K->data, 0, sizeof(double)*2*n_nodes*2*n_nodes);
	int bandwidth = 0;

	for (int t=0; t<n_triplets; t++){
		// K->a[match_perm[triplets[t].i]][match_perm[triplets[t].j]] += triplets[t].val;
		int current_band = match_perm[triplets[t].i] - match_perm[triplets[t].j];
		if (current_band < 0) current_band = -current_band;
		if (current_band > bandwidth) bandwidth = current_band;
	}

	BandMatrix *band_matrix =  allocate_band_matrix(2*n_nodes, bandwidth);
	memset(band_matrix->data, 0, 2*n_nodes*(2*bandwidth+1));
	for (int t=0; t<n_triplets; t++) {
		band_matrix->a[match_perm[triplets[t].i]][match_perm[triplets[t].j]] += triplets[t].val;
	}


	double *RHS_perm = malloc(sizeof(double)*n_nodes*2);

	for (int i = 0; i < 2*n_nodes; i++) {
		RHS_perm[match_perm[i]] = RHS[i];
	}



	// print_matrix(K);

	// Solve linear system
	// double start = get_time();
	// lu(K);
	// double end = get_time();
	// double diff = end - start;
	// printf("lu\t%lf\n",diff);
	// start = get_time();
	// solve(K, RHS_perm);
	// end = get_time();
	// diff = end - start;
	// printf("solve\t%lf\n",diff);

	// print_matrix(K);

	double start = get_time();
	lu_band(band_matrix);
	double end = get_time();
	double diff = end - start;
	printf("lu_band\t%lf\n",diff);
	start = get_time();
	solve_band(band_matrix, RHS_perm);
	end = get_time();
	diff = end - start;
	printf("solve_band\t%lf\n",diff);

	for (int i = 0; i < 2*n_nodes; i++) {
		RHS[i] = RHS_perm[match_perm[i]];
	}
	// Visualization in Gmsh
	visualize_in_gmsh(RHS, gmsh_num, n_nodes);

	// // Run the Gmsh GUI; Comment this line if you do not want the Gmsh window to launch
	gmshFltkRun(&ierr);
	
	// Free stuff
	free_band_matrix(band_matrix);
	free(perm);
	free(match_perm);
	free(RHS_perm);
	free(RHS);
	free(gmsh_num);
	free(coord);
	free(triplets);
	// free_matrix(K);
	gmshFinalize(&ierr);
	return 0;
}