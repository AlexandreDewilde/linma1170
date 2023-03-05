#include "matrix.h"

Matrix * allocate_matrix(int m, int n) {
	Matrix *mat = malloc(sizeof(Matrix));
	mat->data = malloc(sizeof(double)*n*m);
	mat->a = malloc(sizeof(double*)*m);
	mat->n = n; mat->m = m;
	for (int i = 0; i < m; i++) {
		(mat->a)[i] = mat->data + i*n;
	}
	return mat;
}

void free_matrix(Matrix * mat) {
	free(mat->data);
	free(mat->a);
	free(mat);
}
