#include "headers/lu.h"

int lu(Matrix *A) {
    for (int k = 0; k < fmin(A->n, A->m); k++) {
        if(fabs(A->data[k*A->n + k]) < EPS) {
            return -1;
        }
        double pivot = A->a[k][k];
        for (int i = k + 1; i < A->m; i++) {
            A->a[i][k] /= pivot;
            
            for (int j = k + 1; j < A->n; j++) {
                A->a[i][j] -= A->a[i][k]*A->a[k][j];
            }
        }
    }
    return 0;
}

int forward_substitution(Matrix *A, double *y) {
    for (int k = 0; k < A->m; k++) {
        for (int i = 0; i < k; i++) {
            y[k] -= A->a[k][i] * y[i];
        }
    }
    return 0;
}

int backward_substitution(Matrix *A, double *y) {
    for (int k = A->m - 1; k >= 0; k--) {
        for (int i = k + 1; i < A->m; i++) {
            y[k] -= A->a[k][i]*y[i];
        }
        y[k] /= A->a[k][k];
    }
    return 0;
}


int solve(Matrix * LU, double * y) {
    int res = forward_substitution(LU, y);
    if (res == -1) return -1;
    res = backward_substitution(LU, y);
    if (res == -1) return -1;
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

    for (int k = LU->m - 1; k >= 0; k--) {
        int mx = k + LU->k + 1;
        if (mx > LU->m) mx = LU->m;
        for (int i = k + 1; i < mx; i++) {
            y[k] -= LU->a[k][i] * y[i];
        }
        y[k] /= LU->a[k][k];
    }
	
	return 0;
}