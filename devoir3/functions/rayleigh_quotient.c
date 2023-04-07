#include "headers/rayleigh_quotient.h"

double mult(Matrix *A, double *vec) {
    double *vec2 = malloc(sizeof(double)*A->m);
    for (int i = 0; i < A->m; i++) {
        vec2[i] = 0.;
        for (int j = 0; j < A->m; j++) {
            vec2[i] += vec[j] * A->a[j][i];
        }
    }
    double lambda = 0.;
    for (int i = 0; i < A->m; i++) {
        lambda += vec2[i] * vec[i];
    }
    free(vec2);
    return lambda;
}

int rayleigh_quotient_iteration(Matrix *A, double *lambda, double *vec, double eps) {

    *lambda = mult(A, vec);
    Matrix *mat = allocate_matrix(A->m, A->n);
    int iterations = 10000;
    for (int it = 0; it < iterations; it++) {
        for (int i = 0; i < A->m; i++) {
            for (int j = 0; j < A->n; j++) {
                mat->a[i][j] = A->a[i][j];
            }
            mat->a[i][i] -= *lambda;
        }
        int res = lu(mat);
        if (res == -1) {
            printf("error with rayleigh\n");
            exit(EXIT_FAILURE);
        }
        solve(mat, vec);
        double sum = 0.;
        for (int i = 0; i < A->m; i++) sum += vec[i]*vec[i];
        sum = sqrt(sum);
        for (int i = 0; i < A->m; i++) vec[i] /= sum;
        double prev_lambda = *lambda;
        *lambda = mult(A, vec);
        if (fabs(*lambda - prev_lambda) < eps) {
            // printf("Rayleigh quotient iteration, number of iterations to converge : %d\n", it+1);
            return it+1;
        }
    }
    free_matrix(mat);
    return iterations;
}