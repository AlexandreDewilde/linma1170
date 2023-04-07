#ifndef VECTOR_H
#define VECTOR_H
#include <stdlib.h>

typedef struct {
    int *array;
    int n;
    int N;
} vector;

vector *create_vector();

void delete_vector(vector *vec);

void push_vector(vector *vec, int val);

int pop_vector(vector *vec);

int get_vector(vector *vec, int idx);

void swap_vector(vector *vec, int idx_1, int idx_2);

#endif