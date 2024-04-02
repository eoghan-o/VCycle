#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double norm(double* v, int n) {
    double sum = 0;
    int i = 0; 
    for(i = 0; i < n; i++) {
        sum = sum + pow(v[i], 2);
    }
    return pow(sum, 0.5);
}
// a is size m x k
// b is size k x n
void matrixMultiply(double* a, double* b, int m, int n, int k, double* result) {
    int i = 0; // row
    int j = 0; // col
    int l = 0; // vec position
    for(i = 0; i < m; i ++) {
        for(j = 0; j < n; j++) {
            double sum = 0;
            for(l = 0; l < k; l++) {
                sum = sum + a[k*i + l] * b[n*l + j];
            }
            result[n * i + j] = sum;
        }
    }
}

void matrixAdd(double* a, double* b, int m, int n, double* result) {
    int i = 0;
    int j = 0;
    for(i = 0; i < m; i++) {
        for(j = 0; j < n; j++) {
            result[n*i + j] = a[n*i + j] + b[n*i + j];
        }
    }
}

void matVec(double* matrix, double* vec, int m, int n, double* result) {
    int i = 0;
    int j = 0;
    for(i = 0; i < m; i++) {
        double sum = 0;
        for(j = 0; j < n; j++) {
            sum = sum + matrix[n * i + j] * vec[j];
        }
        result[i] = sum;
    }
}

void lower(double* matrix, int m, double* result) {
    int i = 0;
    int j = 0;
    for(i = 0; i < m; i++) {
        for(j = 0; j < i; j++) {
            result[m * i + j] = matrix[m * i + j];
        }
    }
}

void upper(double* matrix, int m, double* result) {
    int i = 0;
    int j = 0;
    for(i = 0; i < m; i++) {
        for(j = i + 1; j < m; j++) {
            result[m * i + j] = matrix[m * i + j];
        }
    }
}

void diag(double* matrix, int m, double* result) {
    int i = 0;
    for(i = 0; i < m; i++) {
        result[m *i + i] = matrix[m * i + i];
    }
} 

//
void transpose(double* A, int m, int n, double* result) {
    int i = 0;
    int j = 0;
    for(i = 0; i < m; i++) {
        for(j = 0; j < n; j++) {
            result[m*j + i] = A[n * i + j];
        }
    }
}

double dotProduct(double u[], double v[], int size) {
    double sum = 0;
    int i = 0;
    for(i = 0; i < size; i++) {
        sum = sum + u[i] * v[i];
    }
    return sum;
}

void printMatrix(double* matrix, int m, int n) {
    int i = 0;
    int j = 0;
    for(i = 0; i < m; i++) {
        for(j = 0; j < n; j++) {
            printf("%f   ", matrix[n*i+j]);
        }
        printf("\n");
    }
}
/*
int main() {
    double A[4] = {2, 7, 1, 2};
    double B[4] = {3, 1, 1, 0};

    double u[2] = {1, 2};
    double v[2] = {-2, -1};
    double dot = dotProduct(u, v, 2);

    double result[4] = {0, 0, 0, 0};
    matrixMultiply(A, B, 2, 2, 2, result);

    double matvec[2] = {0, 0};
    matVec(A, u, 2, 2, matvec);
    
    printf("\nMatrix:\n%f %f\n%f %f\nDot: %f\nMatvec:\n%f\n%f", result[0], result[1], result[2], result[3], dot, matvec[0], matvec[1]);
}
*/