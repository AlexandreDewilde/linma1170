#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"

Matrix * allocate_matrix(int m, int n) {
	Matrix *matrix = malloc(sizeof(matrix));
	matrix->m = m;
	matrix->n = n;
	matrix->data = malloc(m*n*sizeof(double));
	matrix->a = malloc(m*sizeof(double*));

	for(int i = 0; i < m; i++){
		matrix->a[i] = matrix->data + i*n;
	}
	return matrix;
}

void free_matrix(Matrix * mat) {
	free(mat->data);
	free(mat->a);
	free(mat);
}

void print_matrix(Matrix *mat, double cutoff){
	printf("[");
	for(int i = 0; i < mat->m; i++){
		if (i > 0) printf(" ");
		printf("[");
		for(int j = 0; j < mat->n; j++){
			double val = mat->a[i][j];
			val = fabs(val) > cutoff ? val : 0;
			printf("%10.5g", val);
			if (j != mat->n-1) printf(", ");
		}
		if(i == mat->m-1) printf("]");
		printf("]\n");
	}
}

int mult_matrix(Matrix * A, Matrix * B, Matrix * C) {
	if((A->m != C->m) || (B->m != C->m) || (A->n != B->m)){
		printf("ERROR : wrong matrix dimension\n");
		return -1;
	}
	for(int i = 0; i < A->m; i++){
		for(int j = 0; j < B->n; j++){
			double s = 0;
			for(int k = 0; k < A->n; k++){
				s += A->a[i][k]*B->a[k][j];
			}
			C->a[i][j] = s;
		}
	}
	return 0;
}
