#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "matrix.h"

double get_wtime(){
    struct timeval timecheck;
    gettimeofday(&timecheck, NULL);
    long usecs = (long)timecheck.tv_sec * 1000000 + (long)timecheck.tv_usec;
    return usecs * 1e-6;
}

int main() {
	int n = 10;

	// initialiser A de manière aleatoire
	Matrix * A = allocate_matrix(n, n);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			A->a[i][j] = (double)rand()/RAND_MAX; // nombre aléatoire entre 0 et 1
		}
	}

	// initialiser Q a l'identite
	Matrix * Q = allocate_matrix(n, n);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			Q->a[i][j] = (i==j);
		}
	}

	// initialiser R egal a A
	Matrix * R = allocate_matrix(n, n);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			R->a[i][j] = A->a[i][j];
		}
	}

    double tic = get_wtime();
    //-------------- A completer : calculer la QR avec Givens ----------------
	


















	
    //-------------------------------------------------------------------------
    double toc = get_wtime();
    printf("Time spent : %.3f ms\n", (toc - tic)*1e3);

	// print_matrix(A, 1e-15);
	// print_matrix(Q, 1e-15);
	// print_matrix(R, 1e-15);

	//check that A = QR
    double qr_error = 0;
	Matrix * QR = allocate_matrix(n, n);
	Matrix * QT = allocate_matrix(n, n);
    mult_matrix(Q,R,QR);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			qr_error += (QR->a[i][j] - A->a[i][j])*(QR->a[i][j] - A->a[i][j]);
			QT->a[j][i] = Q->a[i][j];
		}
	}
	free_matrix(QR);
    qr_error = sqrt(qr_error);
    printf("||QR - A||  = %.7e\n", qr_error);

	//check that Q is orthogonal
	double q_error = 0;
	Matrix * QQT = allocate_matrix(n, n);
    mult_matrix(Q,QT,QQT);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			q_error += (QQT->a[i][j] - (i == j))*(QQT->a[i][j] - (i == j));
		}
	}
	free_matrix(QQT);
	free_matrix(QT);
    q_error = sqrt(q_error);
    printf("||QQ* - I|| = %.7e\n", q_error);


	// Check that R is upper triangular
	double r_error = 0;
	for(int i = 1; i < n; i++) {
		for(int j = 0; j < i; j++) {
			r_error += R->a[i][j]*R->a[i][j];
		}
	}
    r_error = sqrt(r_error);
    printf("R error     = %.7e\n", r_error);


	free_matrix(R);
	free_matrix(Q);
	free_matrix(A);
}