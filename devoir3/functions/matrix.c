#include "headers/matrix.h"
#include "headers/rmck.h"
#include "headers/lu.h"

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
	for (int i = 0; i < n_nodes*2; i++) {
		perm[i] = i;
	}
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
	mat->data = calloc(width*m, sizeof(double));
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
		for (int j = 0; j < A->n; j++) printf("%0.20f ", A->a[i][j]);
		printf("\n");
	}
}

void print_vector(double *v, int n) {
	for (int i = 0; i < n; i++) {
		printf("%lf ", v[i]);
	}
	printf("\n");
}

int is_symmetric(Matrix * K) {
    for (int i = 0; i < K->m; i++) {
        for (int j = i+1; i < K->n; j++) {
            if (fabs(K->a[i][j]-K->a[j][i]) > 1e15) return 0;
        }
    }
    return 1;
}

Matrix *inverse_matrix(Matrix *A) {
	lu(A);
    Matrix *K_inv = allocate_matrix(A->m, A->n);
    double *vec = malloc(sizeof(double)*A->m);
    for (int i = 0; i < A->n; i++) {
        memset(vec, 0, sizeof(double)*A->m);
        vec[i] = 1.;
        solve(A, vec);
        for (int j = 0; j < A->m; j++) K_inv->a[j][i] = vec[j];
    }
	free(vec);
	return K_inv;
}



Matrix *inverse_matrix_permute(Matrix *A) {
    int *perm = malloc(sizeof(int)*A->m);
    int *inverse_perm = malloc(sizeof(int)*A->m);
    rmck_matrix(perm, A);
    Matrix *A_inv = allocate_matrix(A->m, A->n);
    for (int i = 0; i < A->m; i++) inverse_perm[perm[i]] = i;
    size_t bandwidth = 0;
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->m; j++) {
            if (fabs(A->a[i][j]) > 1e-20) {
                int current_band = (inverse_perm[i] - inverse_perm[j]);
                if (current_band < 0) current_band = -current_band;
                if (current_band > bandwidth) {
                    bandwidth = current_band;
                }
            }
        }
    }

    BandMatrix *band_matrix = allocate_band_matrix(A->m, bandwidth);
    for (int i = 0; i < A->m; i++) {
        int m = i - bandwidth;
        if (m < 0) m = 0;
        int M = i + bandwidth + 1;
        if (M > A->m) M = A->m;
        for (int j = m; j < M; j++) {
            band_matrix->a[i][j] = A->a[perm[i]][perm[j]];
        }
    }

    double *vec = malloc(sizeof(double)*A->m);
    lu_band(band_matrix);
    for (int i = 0; i < A->m; i++) {
        memset(vec, 0, A->m*sizeof(double));
        vec[i] = 1.;
        
        solve_band(band_matrix, vec);
        for (int j = 0; j < A->m; j++) {
            A_inv->a[perm[j]][perm[i]] = vec[j];
        }
    }
    free(vec);
    free_band_matrix(band_matrix);
    free(perm);
    free(inverse_perm);
	return A_inv;
}


Matrix *mult_matrix(Matrix *A, Matrix *B) {
	Matrix *mat = allocate_matrix(A->m, B->n);
	for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < B->n; j++) {
            double sum = 0.;
            for (int k = 0; k < A->n; k++) sum += A->a[i][k] * B->a[k][j];
            mat->a[i][j] = sum;
        }
    }
	return mat;
}


void reduce_matrix(Matrix **K, Matrix **M, char *boundary_bool, size_t n_boundary_nodes) {

    size_t new_size = (*K)->m - 2*n_boundary_nodes;
    Matrix *K_reduced = allocate_matrix(new_size, new_size);
    Matrix *M_reduced = allocate_matrix(new_size, new_size);
    int current_row = 0;
    for (int i = 0; i < (*K)->m; i++) {
        if (boundary_bool[i/2]) continue;
        int current_col = 0;
        for (int j = 0; j < (*K)->m; j++) {
            if (boundary_bool[j/2]) continue;
            K_reduced->a[current_row][current_col] = (*K)->a[i][j];
            M_reduced->a[current_row][current_col] = (*M)->a[i][j];
            ++current_col;
        }
        ++current_row;
    }
    free_matrix(*K);
    free_matrix(*M);
    *M = M_reduced;
    *K = K_reduced;
}