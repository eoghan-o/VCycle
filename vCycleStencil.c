#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearAlgebra.h"

double stencilValue(int i, int j, int n, double sigma, double h) {
    if(i == j) { 
        return 2 + sigma * pow(h,2.0);
    } else if(i == -1 || j == -1 || i == n + 1 || j == n + 1) {
        return 0;
    } else if(i == j + 1 || j == i + 1) {
        return -1.0;
    } else {
        return 0;
    }
}

void gaussSeidel(double* x, double* b, int n, int nu, double sigma, double h) {
    int i = 0;
    int j = 0;
    int k = 0;
    for(k = 0; k < nu; k++) {
        for(i = 0; i < n - 1; i++) {
            double sum = 0;
            for(j = 0; j < n - 1; j++) {
                if(j != i) {
                    sum = sum + stencilValue(i, j, n - 1, sigma, h) * x[j];
                }
            }
            x[i] = (b[i] - sum)/stencilValue(i, i, n - 1, sigma, h);
        }
    }
}

void prolongation(double* v2h, int n, double* vh) {
    int i = 0;
    int a = 0;
    int b = 0;
    vh[0] = v2h[0]/2.0;
    vh[n-2] = v2h[n/2 - 2];
    for(i = 1; i < n-2; i++) {
        vh[i] = (v2h[a] + v2h[b])/2.0;
        if(i%2) {
            a++;
        } else {
            b++;
        }
    }
}

void restriction(double* vh, int n, double* v2h) {
    int i = 0;
    int j = 0;
    for(i = 0; i < n/2 - 1; i++) {
        v2h[i] = (vh[j] + 2*vh[j+1] + vh[j+2])/4.0;
        j = j + 2;
    }
}

void vCycle(double* u, double* f, int n, int nu1, int nu2, double sigma, double h) {
    if(n != 2) {
        gaussSeidel(u, f, n, nu1, sigma, h);

        double residual[n-1];
        int i = 0;
        int j = 0;
        for(i = 0; i < n-1; i++) {
            double sum = 0;
            for(j = i - 1; j < i + 2; j++) {
                sum = sum + stencilValue(i, j, n - 1, sigma, h) * u[j];
            }
            residual[i] = f[i] - sum;
        }

        double residual2h[n/2-1];
        restriction(residual, n, residual2h);

        double error2h[n/2-1];
        for(i = 0; i < n/2-1; i++) {
            error2h[i] = 0;
        }

        vCycle(error2h, residual2h, n/2, nu1, nu2, sigma, 2*h);

        double errorh[n-1];
        prolongation(error2h, n, errorh);

        for(i = 0; i < n-1; i++) {
            u[i] = u[i] + errorh[i];
        }

        gaussSeidel(u, f, n, nu2, sigma, h);
    } else {
        u[0] = f[0]/(2+sigma*pow(h,2));
    }
}     

int simulate(int k, FILE* fptr) {

    int n = (int)pow(2,k);

    int nu1 = 2;
    int nu2 = 1;

    double h = 1.0/n;

    double sigma = 1.0;

    double u[n-1];
    double f[n-1];

    int i = 0;
    int j = 0;

    for(i = 0; i < n; i++) {
        u[i] = 0;
        f[i] = i*h * pow(h,2);
    }

    double residualNorm = 0.0;
    double tolerance = pow(10,-6);

    i = 0;
    int a = 0;
    int b = 0;
    double temp[n-1];
    double residual[n-1];
    double sum = 0;
    while(i == 0 || residualNorm > tolerance) {
        vCycle(u, f, n, nu1, nu2, sigma, h);

        for(a = 0; a < n-1; a++) {
            sum = 0.0;
            for(b = a - 1; b < a + 2; b++) {
                sum = sum + stencilValue(a, b, n - 1, sigma, h) * u[b];
            }
            residual[a] = f[a] - sum;
        }

        residualNorm = norm(residual, n-1);
        //printf("\nCycle %d\nNorm of residual: %f\n", i+1, residualNorm);
        i = i+1;

        fprintf(fptr, "%d,", i);
        for(j = 0; j < n-2; j++) {
            fprintf(fptr, "%.10f,", u[j]);
        }
        fprintf(fptr, "%.10f\n", u[n-2]);
    }

    /*
    printf("Solution vector u:\n");
    for(i = 0; i < n-1; i++) {
        printf("%f\n", u[i]);
    }
    */
    return i;
}

int main() {
    FILE* fptr;
    fptr = fopen("solution.csv", "w");

    int k = 3;
    int num = simulate(k, fptr);

    //for(k = 3; k <= 10; k++) {
    //    printf("k = %d -> Num iterations: %d\n", k, simulate(k, fptr));
    //}
    return 0;
}