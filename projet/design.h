#ifndef DESIGN_H // le header guard
#define DESIGN_H 

void designTuningFork(double r1, double r2, double e, double l, double meshSizeFactor, char * filename);


void designTuningForkSymmetric(double r1, double r2, double e, double l, double meshSizeFactor, char * filename);

void designTuningForkSymmetric2Layer(double e, double d1, double d2, double h1, double h2, double l1, double l2, double meshSizeFactor, char * filename);

void designTuningForkNLayer(double* d, double* h, double* l, size_t n, double meshSizeFactor, char * filename);

void designTuningForkSymmetricNLayer(double* d, double* h, double* l, size_t n, double meshSizeFactor, char * filename);

#endif