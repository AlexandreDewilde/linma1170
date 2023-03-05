#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "matrix.h"
#include "lu.h"

void mult_lu(Matrix *LU, Matrix *ret) {
    for (int i = 0; i < LU->m; i++) {
        for (int j = 0; j < LU->n; j++) {
            for (int k = 0; k < fmin(i,j)+1; k++) {
                double L = LU->a[i][k];
                if (i == k) L = 1.;
                ret->a[i][j] += L*LU->a[k][j];
            }
        }
    }
}

int check_same_matrix(double *m1, double *m2, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            int diff = (m1[i*n+j] - m2[i*n+j]);
            if (diff < 0) diff = - diff;
            if (diff > EPS) {
                return -1; 
            }
        }
    }
    return 0;
}

Matrix *random_matrix(int m, int n) {
    Matrix *mat = allocate_matrix(m, n);
    for (int i = 0; i < mat->m; i++) {
        for (int j = 0; j < mat->n; j++) {
            double nb = (double) (rand() % 1000);
            mat->a[i][j] = nb;
        }
        //Avoid pivot < EPS
        if (i < n && i < m && mat->a[i][i] < EPS)
            mat->a[i][i] += 1.; 
    }
    return mat;
}

int test_random_lu(Matrix *mat) {
    int m = mat->m; int n = mat->n;
    size_t size = sizeof(double) * n * m;
    double *copy = malloc(size);
    memcpy(copy, mat->data, size);

    int res = lu(mat);
    if (res == -1) {
        // means that pivot is null but the generated matrix has no null pivot so WA
        free(copy);
        return -1;
    }

    Matrix *mult = allocate_matrix(m, n);
    memset(mult->data, 0, n*m*sizeof(double));
    
    mult_lu(mat, mult);
    res = check_same_matrix(mult->data, copy, m, n);
    if (res != 0) {
        free(copy);
        return -1;
    }
    free_matrix(mult);
    free(copy);
    return 0;
}

int test_solve(Matrix *LU, double *A) {
    int m = LU->m, n = LU->n;
    double y[m], y_copy[m];
    for (int i = 0; i < m; i++) {
        double nb = (double) (rand() % 1000);
        y[i] = y_copy[i] = nb;
    }

    solve(LU, y);

    double y_calc[m];
    for (int i = 0; i < m; i++) {
        y_calc[i] = 0;
        for (int j = 0; j < n; j++) {
            y_calc[i] += y[j]*A[i*m+j]; 
        }
    }

    for (int i = 0; i < m; i++) {
        double diff = y_calc[i] - y_copy[i];
        if (diff < 0) diff = -diff;
        // Large espsilon because error can be big
        // Change by a calcul?
        if (diff > 0.1) {
            printf("%lf\n", diff);
            return -3;
        }
    }
    return 0;
}

// base on this SO post https://stackoverflow.com/questions/3557221/how-do-i-measure-time-in-c
double get_time()
{
    struct timespec now;
    timespec_get(&now, TIME_UTC);
    return now.tv_sec + now.tv_nsec*1e-9;
}

int main(int argc, char *argv[]) {
    // Test square matrices
    for (int i = 0; i < 1000; i++) {
        for (int j = 0; j < 5; j++) {
            Matrix *mat = random_matrix(i+1, i+1);
            size_t size = (i+1)*(i+1)*sizeof(double);
            double *copy = malloc(size);
            memcpy(copy, mat->data, size);
            int res = test_random_lu(mat);
            if (res != 0) {
                free(copy);
                free_matrix(mat);
                printf("%d\n",res);
                exit(EXIT_FAILURE);
            }
            res = test_solve(mat, copy);
            free(copy);
            free_matrix(mat);
            if (res != 0) {
                printf("%d\n", res);
                exit(EXIT_FAILURE);
            }
        }
    }

    

    // Test all type of matrices
    for (int i = 0; i < 1000; i++) {
        int n = rand()%100+1; int m = rand()%100+1;
        Matrix *mat = random_matrix(n, m);
        size_t size = n*m*sizeof(double);
        double *copy = malloc(size);
        memcpy(copy, mat->data, size);
        int res = test_random_lu(mat);
        free(copy);
        free_matrix(mat);
        if (res != 0) {
            printf("%d\n",res);
            exit(EXIT_FAILURE);
        }
    }

    // Test time
    printf("function\tsize\ttime\n");
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 1000; j++) {
            Matrix *mat = random_matrix(i+1, i+1);
            double y[i+1];
            for (int k = 0; k < i+1; k++) {
                double nb = (double) (rand() % 1000);
                y[k] = nb;
            }

            double start = get_time();
            lu(mat);
            double end = get_time();
            double diff = end - start;
            printf("lu\t%d\t%lf\n",i+1,diff);
            
            start = get_time();
            solve(mat, y);
            end = get_time();
            diff = end - start;
            printf("solve\t%d\t%lf\n",i+1,diff);
            free_matrix(mat);
        }
    }
    return 0;
}