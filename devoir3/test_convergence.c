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
    printf("size;power_iteration;rayleigh_quotient\n");
    for (int it = 0; it < 10; it++) {
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
            power_iteration(mat, &lambda, vec, 1e-20);
            for (int i = 0; i < sz; i++) {
                if (!(rand() % 4))
                    vec[i] += 1. / (rand() % 100 + 100000.);
            }
            double vec2[sz]; memcpy(vec2, vec, sz*sizeof(double));
            int it_pi = power_iteration(mat, &lambda, vec, 1e-20);
            int it_rqi = rayleigh_quotient_iteration(mat, &lambda, vec2, 1e-20);
            printf("%d;%d;%d\n", sz, it_pi, it_rqi);
        }
        free_matrix(mat);
    }
}
