#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearAlgebra.h"

double func(double x) {
    return x;
}

// Takes R^{n/2 - 1} to R^{n - 1}
void prolongation(double* v2h, int n, double* vh) {
    int i = 0;
    int a = 0;
    int b = 0;
    vh[0] = v2h[0]/2.0; // Update first
    vh[n-2] = v2h[n/2 - 2]; // Update last
    for(i = 1; i < n-2; i++) {
        vh[i] = (v2h[a] + v2h[b])/2.0; // Update middle
        if(i%2) {
            a++;
        } else {
            b++;
        }
    }
}

// Depracated
/*
void prolongation(double* v2h, int n, double* vh) {
    double I2h[(n-1)*(n/2-1)];
    int i = 0;
    int j = 0;
    for(j = 0; j < n/2-1; j++) {
        for(i = 0; i < n-1; i++) {
            if(2*j == i) {
                I2h[(n/2-1)*i + j] = 1;
            } else if(i == 2*j + 1) {
                I2h[(n/2-1)*i + j] = 2;
            } else if(i == 2*j + 2) {
                I2h[(n/2-1)*i + j] = 1;
            } else {
                I2h[(n/2-1)*i + j] = 0;
            }
        }
    }
    matVec(I2h, v2h, n-1, n/2 - 1, vh);
}
*/

// Takes R^{n-1} to R^{n/2 - 1}
void restriction(double* vh, int n, double* v2h) {
    int i = 0;
    int j = 0;
    for(i = 0; i < n/2 - 1; i++) {
        v2h[i] = (vh[j] + 2*vh[j+1] + vh[j+2])/4.0;
        j = j + 2;
    }
}

// Depracated
/*
void restriction(double* vh, int n, double* v2h) {
    double Ih[(n/2-1)*(n-1)];
    int i = 0;
    int j = 0;
    for(i = 0; i < n/2-1; i++) {
        for(j = 0; j < n-1; j++) {
            if(j == 2*i) {
                Ih[(n-1)*i + j] = 1;
            } else if(j == 2*i+1) {
                Ih[(n-1)*i + j] = 2;
            } else if(j == 2*i + 2) {
                Ih[(n-1)*i + j] = 1;
            } else {
                Ih[(n-1)*i + j] = 0;
            }
        }
    }
    matVec(Ih, vh, (n/2 - 1), n-1, v2h);
}
*/

void gaussSeidel(double* A, double* x, double* b, int n, int nu) {
    int i = 0;
    int j = 0;
    int k = 0;
    for(k = 0; k < nu; k++) {
        for(i = 0; i < n; i++) {
            double sum = 0;
            for(j = 0; j < n; j++) {
                if(j != i) {
                    sum = sum + A[n*i + j] * x[j];
                }
            }
            x[i] = (b[i] - sum)/A[n*i + i];
        }
    }
}

void generatePoissonMatrix(double* result, double sigma, double h, int n) {
    int i = 0;
    int j = 0;
    for(i = 0; i < n-1; i++) {
        for(j = 0; j < n-1; j++) {
            if(i == j) {
                result[(n-1)*i + j] = 2 + sigma*pow(h,2);
            } else if(j == i + 1 || i == j + 1) {
                result[(n-1)*i + j] = -1;
            } else {
                result[(n-1)*i + j] = 0;
            }
        }
    }
}

// Solve Au = f
void vCycle(double* A, double* u, double* f, int n, int nu1, int nu2, double h, double sigma) {
    if(n != 2) {
        gaussSeidel(A, u, f, n-1, nu1);
        printMatrix(u, n-1, 1);

        double residual[n-1];
        matVec(A, u, n-1, n-1, residual); // r = Au

        int i = 0;
        for(i = 0; i < n-1; i++) {
            residual[i] = f[i] - residual[i]; // r <- f - r = f - Au
        }

        double residual2h[n/2-1];
        restriction(residual, n, residual2h); // form r2h

        double A2h[(n/2-1)*(n/2-1)];
        generatePoissonMatrix(A2h, sigma, 2*h, n/2); // form A2h

        double error2h[n/2-1]; // Generate the initial guess e2h
        for(i = 0; i < n/2-1; i++) {
            error2h[i] = 0;
        }

        vCycle(A2h, error2h, residual2h, n/2, nu1, nu2, 2*h, sigma); // Do it again! A2h e2h = r2h

        double errorh[n-1]; 
        prolongation(error2h, n, errorh);  // Interpolate e2h to eh

        for(i = 0; i < n-1; i++) {
            u[i] = u[i] + errorh[i];
        }
        
        gaussSeidel(A, u, f, n-1, nu2); // Ah uh = fh w/ initial guess uh, nu2 times
    } else {
        // A is a 1 by 1, f is a 1 x 1
        // Data goes into initialGuess
        u[0] = f[0]/A[0];
    }
}


int main() { 
    int k = 3;
    int n = (int)(pow(2, k));

    int nu1 = 2;
    int nu2 = 1;

    double h = 1.0/n;    

    // Define f vector
    double sigma = 1.0;

    // Poisson problem Matrix
    double A[(n-1)*(n-1)];
    generatePoissonMatrix(A, sigma, h, n);

    double u[n-1];
    double f[n-1];
    int i = 0;
    int j = 0;

    for(i = 0; i < n-1; i++) {
        u[i] = 0;
        f[i] = func(i*h) * pow(h,2);
    }

    double residualNorm = 0.0;
    double tolerance = pow(10, -6);

    i = 0;
    while(i == 0 || residualNorm > tolerance) {
        double temp[n-1];
        for(j = 0; j < n-1; j++) {
            temp[j] = u[j];
        }

        vCycle(A, u, f, n, nu1, nu2, h, sigma);
        
        for(j = 0; j < n; j++) {
            temp[j] = u[j] - temp[j];
        }
        residualNorm = norm(temp, n-1);
        printf("\nCycle %d\nNorm of u(i+1)-u(i) %f\n", i+1, residualNorm);
        i = i+1;      
    }
    printf("Solution vector u:\n");
    for(i = 0; i < n-1; i++) {
        printf("%f\n", u[i]);
    }
    return 0;
}