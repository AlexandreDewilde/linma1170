#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "matrix.h"
#include "rmck.h"
#include "lu.h"

Matrix * allocate_matrix(int m, int n) {
	Matrix * mat = (Matrix*) malloc(sizeof(Matrix));
	mat->m = m, mat->n = n;
	mat->data = (double*) calloc(m*n,sizeof(double));
	mat->a = (double**) calloc(m, sizeof(double*));
	for(int i = 0; i < m; i++)
		mat->a[i] = mat->data+i*n;
	return mat;
}

void free_matrix(Matrix * mat) {
	free(mat->a);
	free(mat->data);
	free(mat);
}

void print_vector(double * v, int n) {
	for(int i = 0; i < n; i++)
		printf("%.3e ", v[i]);
	printf("\n");
}

void print_matrix(Matrix * A) {
	for(int i = 0; i < A->m; i++)
		print_vector(A->a[i], A->n);
}

int is_symmetric(Matrix * K) {
	int symmetric = 1;
	for(int i = 0; i < K->m; i++)
		for(int j = i+1; j < K->n; j++)
			if(fabs((K->a[i][j] - K->a[j][i]) / K->a[i][j]) > 1e-12) {
				printf("%d %d\n", i, j);
				printf("%lf %lf\n", K->a[i][j], K->a[j][i]);
				symmetric = 0;
			}
	return symmetric;
}

typedef struct {
	int i;		// index
	double x,y;	// coordinates
} Node;

// Comparateur
static int cmp(const void * a, const void * b) {
	Node * na = (Node *) a;
	Node * nb = (Node *) b;
	if (na->x > nb->x) return 1;
	else return -1;
}

int compute_permutation(int * perm, double * coord, int n_nodes) {

	// // We assume perm is allocated but not initialized
	// for(int i = 0; i < n_nodes; i++) 
	// 	perm[i] = i;

	// qsort_r(perm, n_nodes, sizeof(int), coord, cmp);

	// Create Node structs
	Node * nodes = malloc(n_nodes * sizeof(Node));
	for(int i = 0; i < n_nodes; i++) {
		nodes[i].i = i;
		nodes[i].x = coord[2*i];
		nodes[i].y = coord[2*i+1];
	}

	// Sort nodes
	qsort(nodes, n_nodes, sizeof(Node), cmp);

	// Fetch permutation (we assume perm is allocated)
	for(int i = 0; i < n_nodes; i++)
		perm[i] = nodes[i].i;

	return 0;
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

int compute_bandwitdh(Matrix *A, int *inverse_perm) {
	size_t bandwidth = 0;
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->m; j++) {
            if (fabs(A->a[i][j]) > 1e-20) {
                int current_band = (inverse_perm[i] - inverse_perm[j]);
				// Absolute value
                if (current_band < 0)
					current_band = -current_band;
                if (current_band > bandwidth) {
                    bandwidth = current_band;
                }
            }
        }
    }
	return bandwidth;
}



void inverse_matrix_permute(Matrix *A, Matrix *M) {
	size_t m = A->m;

	// Compute permutation
    int perm[m], inverse_perm[m];
    rmck_matrix(perm, A);
    for (int i = 0; i < A->m; i++)
		inverse_perm[perm[i]] = i;

	// Compute bandwidth
    size_t bandwidth = compute_bandwitdh(A, inverse_perm);

    BandMatrix *band_matrix = allocate_band_matrix(m, bandwidth);
    for (int i = 0; i < A->m; i++) {
        int m = i - bandwidth;
        if (m < 0) m = 0;
        int M = i + bandwidth + 1;
        if (M > A->m) M = A->m;
        for (int j = m; j < M; j++) {
            band_matrix->a[i][j] = A->a[perm[i]][perm[j]];
        }
    }

	double vec[m];
    lu_band(band_matrix);
    for (int i = 0; i < A->m; i++) {
		for (int j = 0; j < A->m; j++)
        	vec[j] = M->a[perm[j]][perm[i]];
        
        solve_band(band_matrix, vec);
        for (int j = 0; j < A->m; j++) {
            M->a[perm[j]][perm[i]] = vec[j];
        }
    }
    free_band_matrix(band_matrix);
}