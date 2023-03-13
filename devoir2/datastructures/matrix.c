#include "headers/matrix.h"
#include "../headers/rmck.h"

Matrix * allocate_matrix(int m, int n) {
	Matrix *mat = malloc(sizeof(Matrix));
	mat->data = calloc(n*m, sizeof(double));
	mat->a = malloc(sizeof(double*)*m);
	mat->n = n; mat->m = m;
	for (int i = 0; i < m; i++) {
		(mat->a)[i] = mat->data + i*n;
	}
	return mat;
}


int compute_permutation(int * perm, double * coord, int n_nodes, Triplet * triplets, int n_triplets) {
	rmck(perm, triplets, n_triplets, n_nodes);
	return 0;
}

void free_matrix(Matrix * mat) {
	free(mat->data);
	free(mat->a);
	free(mat);
}

BandMatrix * allocate_band_matrix(int m, int k) {
	BandMatrix *mat = malloc(sizeof(BandMatrix));
	mat->m = m; mat->k = k;
	int width = (2*k+1);
	mat->data = malloc(sizeof(double)*width*m);
	mat->a = malloc(sizeof(double*)*m);
	for (int i = 0; i < m; i++) {
		mat->a[i] = mat->data + width*i - i;
	}
	return mat;
}

void free_band_matrix(BandMatrix * mat) {
	free(mat->a);
	free(mat->data);
	free(mat);
}


void print_matrix(Matrix * A) {
	for (int i = 0; i < A->m; i++) {
		for (int j = 0; j < A->n; j++) printf("%lf ", A->a[i][j]);
		printf("\n");
	}
}

void print_vector(double *v, int n) {
	for (int i = 0; i < n; i++) {
		printf("%lf ", v[i]);
	}
	printf("\n");
}