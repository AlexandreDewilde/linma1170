#ifndef RAYLEIGH_QUOTIENT_H
#define RAYLEIGH_QUOTIENT_H
#include "matrix.h"
#include "lu.h"

int rayleigh_quotient_iteration(Matrix *A, double *lambda, double *vec, double eps);

#endif