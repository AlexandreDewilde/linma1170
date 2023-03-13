#include "headers/vector.h"
#include <stdlib.h>

vector *create_vector() {
    vector *vec = malloc(sizeof(vector));    
    vec->n = 0;
    vec->N = 2;
    vec->array = malloc(sizeof(int) * 2);
    return vec;
}

void delete_vector(vector *vec) {
    free(vec->array);
    free(vec);
}

void push_vector(vector *vec, int val) {
    if (vec->N <= vec->n) {
        vec->N = vec->N * 2;
        vec->array = realloc(vec->array, vec->N * sizeof(int));
    }
    vec->array[vec->n++] = val; 
}

int pop_vector(vector *vec) {
    return vec->array[vec->n-- - 1];
}

int get_vector(vector *vec, int idx) {
    return vec->array[idx];
}

void swap_vector(vector *vec, int idx_1, int idx_2) {
    int temp = vec->array[idx_1];
    vec->array[idx_1] = vec->array[idx_2];
    vec->array[idx_2] = temp;
}