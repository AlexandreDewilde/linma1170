#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "math.h"
#include "functions/headers/matrix.h"
#include "functions/headers/power_method.h"
#include "functions/headers/rayleigh_quotient.h"

int main() {
    srand(time(NULL));
    printf("size;power_iteration;hybrid_method\n");
    for (int it = 0; it < 100; it++) {
        size_t sz = (it + 1) * 10;
        Matrix *mat = allocate_matrix(sz, sz);
        for (int k = 0; k < 5; k++) {
            for (int i = 0; i < sz; i++) {
                mat->a[i][i] = 1.;
                for (int j = i + 1; j < sz; j++) {
                    mat->a[i][j] = mat->a[j][i] = ((double) (rand() % 100 + 1)) / (rand() % 5 + 1);
                }
            }
            double vec[sz];
            for (int i = 0; i < sz; i++) {
                vec[i] = 1. / sz;
            }
            double lambda = 0.;
            int it_pi = power_iteration(mat, &lambda, vec, 1e-20);
            for (int i = 0; i < sz; i++) {
                vec[i] = 1. / sz;
            }
            // printf("Hybrid method\n");
            int it_hybrid = power_iteration(mat, &lambda, vec, 1e-5);
            it_hybrid += rayleigh_quotient_iteration(mat, &lambda, vec, 1e-20);
            printf("%d;%d;%d\n", sz, it_pi, it_hybrid);
        }
        free_matrix(mat);
    }
}
