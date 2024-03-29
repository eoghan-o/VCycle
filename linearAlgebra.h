#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void matrixMultiply(double* a, double* b, int m, int n, int k, double* matrix);
void matrixAdd(double* a, double* b, int m, int n, double* result);
void matVec(double* matrix, double* vec, int m, int n, double* result);
void lower(double* matrix, int m, double* result);
void upper(double* matrix, int m, double* result);
void diag(double* matrix, int m, double* result);
double dotProduct(double u[], double v[], int size);
void printMatrix(double* matrix, int m, int n);