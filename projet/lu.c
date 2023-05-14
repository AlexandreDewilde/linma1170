#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lu.h"
#include <cblas.h>

int lu(Matrix* A) {
	size_t const m = A->m;
	for(size_t k = 0; k < m - 1; k++) {
		if(fabs(A->a[k][k]) < EPS)
			return -1;

		for(size_t j = k + 1; j < m; j++) {
			A->a[j][k] /= A->a[k][k];
			for(size_t l = k + 1; l < m; l++)
				A->a[j][l] -= A->a[j][k] * A->a[k][l];
		}
	}
	return 0;
}

int solve(Matrix* LU, double* y) {
	size_t const m = LU->m;
	double* x = y;

	// Résolution de L*x = y par substitution avant
	for(size_t i = 0; i < m; i++) {
		for(size_t j = 0; j < i; j++)
			x[i] -= LU->a[i][j] * x[j];
	}

	// Résolution de U*x = L^{-1} y par substitution arrière
	for(int i = m - 1; i >= 0; i--) {
		for(size_t j = i + 1; j < m; j++)
			x[i] -= LU->a[i][j] * x[j];
		x[i] /= LU->a[i][i];
	}

	return 0;
}

int lu_band(BandMatrix * A) {
	for (int k = 0; k < A->m - 1; k++) {
		double pivot = A->a[k][k];
		double pivot_abs = pivot;
		if (pivot < 0) pivot_abs = 0 - pivot;
		if (pivot_abs < EPS) return -1;
		int max_i = k + A->k + 1;
		if (max_i > A->m) max_i = A->m;
		for (int i = k + 1; i < max_i; i++) {
			A->a[i][k] /= pivot;
			int max_j = A->k + k + 1;
			if (max_j > A->m) max_j = A->m;
			for (int j = k + 1; j < max_j; j++) {
				A->a[i][j] -= A->a[i][k]*A->a[k][j];
			}
		}
	}
	return 0;
}

int solve_band(BandMatrix * LU, double * y) {
	for (int k = 0; k < LU->m; k++) {
		int min = k - LU->k;
		if(min < 0) min = 0;
		for (int i = min; i < k; i++) {
			y[k] -= LU->a[k][i] * y[i];
		}
	}
	// cblas_dtbsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasUnit, LU->m, LU->k, LU->data, 2*LU->k + 1, y, 1);
	cblas_dtbsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, LU->m, LU->k, LU->data, 2*LU->k+1, y, 1);

	// for (int k = LU->m - 1; k >= 0; k--) {
	//     int mx = k + LU->k + 1;
	//     if (mx > LU->m) mx = LU->m;
	//     for (int i = k + 1; i < mx; i++) {
	//         y[k] -= LU->a[k][i] * y[i];
	//     }
	//     y[k] /= LU->a[k][k];
	// }

	return 0;
}