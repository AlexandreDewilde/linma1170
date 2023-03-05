#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

int main() {
	int nrow = 100, ncol = 200;
	Matrix * mat = allocate_matrix(nrow, ncol);
	for(int i = 0; i < nrow; i++) {
		for(int j = 0; j < ncol; j++) {
			mat->a[i][j] = (double)rand()/RAND_MAX; // nombre aléatoire entre 0 et 1
		}
	}
	printf("A[12][116] = %.3lf\n", mat->a[12][116]);
	free_matrix(mat);
}