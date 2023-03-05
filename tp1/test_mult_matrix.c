#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

int main() {
	Matrix * A = allocate_matrix(2,3);
	Matrix * B = allocate_matrix(3,2);
	Matrix * C = allocate_matrix(2,2);
	double ** a = A->a;
	a[0][0] = 1; a[0][1] = 2; a[0][2] = 1;
	a[1][0] = 1; a[1][1] = 1; a[1][2] = 2;

	double ** b = B->a;
	b[0][0] = 2; b[0][1] = 2;
	b[1][0] = 1; b[1][1] = 1;
	b[2][0] = 1; b[2][1] = 2;

	mult_matrix(A,B,C);
	double **c = C->a;
	for (int i=0; i<2; i++){
		for (int j=0; j<2; j++){
			printf("%.2lf ", c[i][j]);
		}
		printf("\n");
	}
	free_matrix(A);
	free_matrix(B);
	free_matrix(C);
}