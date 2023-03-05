#include <stdlib.h>

#ifndef MATRIX_H // le header guard
#define MATRIX_H 

typedef struct Matrix {
	int m, n;       // dimensions de la matrice
	double * data;  // tableau 1D de taille m*n contenant les entrées de la matrice
	double ** a;    // tableau 1D de m pointeurs vers chaque ligne, pour pouvoir appeler a[i][j]
} Matrix;

Matrix * allocate_matrix(int m, int n); // alloue une matrice de dimensions données
void free_matrix(Matrix * mat);         // libère la mémoire de la structure Matrix donnée

#endif