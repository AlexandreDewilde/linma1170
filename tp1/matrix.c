#include <stdlib.h>
#include "matrix.h"

// NB: voir tp1.pdf pour l'énoncé complet du TP.

Matrix * allocate_matrix(int m, int n) {
	// Exercice 2 : À vous de jouer ! Il faut :
	// 1. Allouer une structure Matrix
	// 2. Assigner m et n
	// 3. Allouer `data`
	// 4. Allouer `a` et y remplir les pointeurs vers les lignes
	Matrix *mat = malloc(sizeof(Matrix));
	mat->data = malloc(sizeof(double)*n*m);
	mat->n = n; mat->m = m;
	mat->a = malloc(sizeof(double*)*m);
	for (int i = 0; i < m; i++) {
		(mat->a)[i] = mat->data + i*n;
	}
	return mat;
}

void free_matrix(Matrix * mat) {
	// Exercice 2 : À vous de jouer ! Il faut :
	// 1. Libérer ce qui a été alloué dans `mat`
	// 2. Libérer `mat`
	free(mat->data);
	free(mat->a);
	free(mat);
}

int mult_matrix(Matrix * A, Matrix * B, Matrix * C) {
	// Exercice 3
	for (int i = 0; i < A->m; i++) {
		for (int j = 0; j< B->n; j++) {
			int sum = 0;
			for (int k = 0; k < A->n; k++) sum += A->data[i*A->n + k] * B->data[k*B->n + j];
			C->data[i*C->n + j] = sum;
		}
	}

}
