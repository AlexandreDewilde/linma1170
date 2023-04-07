#include "vector.h"
#include <string.h>
#include <stdlib.h>
#ifndef RMCK_H
#define RMCK_H
#include "matrix.h"

void rmck(int *perm, Triplet *triplets, int n_triplets, int n_nodes);

void rmck_matrix(int *perm, Matrix *A);


#endif