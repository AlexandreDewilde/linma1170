#include "headers/power_method.h"


int power_iteration(Matrix *A, double *lambda, double *vec, double eps) {
    double *vec2 = malloc(sizeof(double)*A->m);
    int iterations = 10000;
    *lambda = 0.;
    for (int it = 0; it < iterations; it++) {
        for (int i = 0; i < A->m; i++) {
            vec2[i] = 0.;
            for (int j = 0; j < A->m; j++ ) {
                vec2[i] += A->a[i][j]*vec[j];
            }
        }
        double norm = 0.;
        for (int k = 0; k < A->m; k++) {
            norm += vec2[k]*vec2[k];
        }
        norm = sqrt(norm);

        for (int k = 0; k < A->m; k++) {
            vec[k] = vec2[k] / norm;
        }
        double prev_lambda = *lambda;
        *lambda = 0.;
        for (int i = 0; i < A->m; i++) {
            vec2[i] = 0.;
            for (int j = 0; j < A->m; j++) {
                vec2[i] += vec[j] * A->a[j][i];
            }
        }
        for (int i = 0; i < A->m; i++) {
            *lambda += vec2[i] * vec[i];
        }
        double den = 0.;
        for (int i = 0; i < A->m; i++) {
            den += vec[i]*vec[i];
        }
        *lambda /= den;
        if (fabs((*lambda)-prev_lambda) < eps) {
            // printf("Power Method, number of iterations to converge : %d\n", it+1);
            return it + 1;
        }
    }
    
    free(vec2);
    return iterations;
}