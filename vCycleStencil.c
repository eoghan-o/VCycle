#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearAlgebra.h"
#include <string.h>

double stencilValue(int i, int j, int n, double sigma, double h) {
    if(i == j) { 
        return 2.0/pow(h,2) + sigma;
    } else if(i == -1 || j == -1 || i == n-1 || j == n-1) {
        return 0;
    } else if(i == j + 1 || j == i + 1) {
        return -1.0/pow(h,2);
    } else {
        return 0;
    }
}

void gaussSeidel(double* x, double* b, int n, int nu, double sigma, double h) {
    int i = 0;
    int j = 0;
    int k = 0;
    double sum = 0;
    for(k = 0; k < nu; k++) {
        for(i = 0; i < n - 1; i++) {
            sum = 0;
            for(j = 0; j < n - 1; j++) {
                if(j != i) {
                    sum = sum + stencilValue(i, j, n, sigma, h) * x[j]; 
                }
            }
            x[i] = (b[i] - sum)/stencilValue(i, i, n, sigma, h); // x = D^(-1)(b - (U+L)x)
        }
    }
}

void prolongation(double* v2h, int n, double* vh) {
    int i = 0;
    int a = 0;
    int b = 0;
    vh[0] = v2h[0]/2.0;
    vh[n-2] = v2h[n/2 - 2]/2.0;
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
    int i = 0;
    int j = 0;
    double sum = 0.0;
    double residual[n-1];
    double residual2h[n/2-1];
    double error2h[n/2-1];
    double errorh[n-1];
    if(n != 2) {
        gaussSeidel(u, f, n, nu1, sigma, h);

        for(i = 0; i < n-1; i++) {
            sum = 0;
            for(j = i - 1; j < i + 2; j++) {
                sum = sum + stencilValue(i, j, n, sigma, h) * u[j];
            }
            residual[i] = f[i] - sum;
        }

        restriction(residual, n, residual2h);

        for(i = 0; i < n/2-1; i++) {
            error2h[i] = 0;
        }

        vCycle(error2h, residual2h, n/2, nu1, nu2, sigma, 2*h);
        
        prolongation(error2h, n, errorh);

        for(i = 0; i < n-1; i++) {
            u[i] = u[i] + errorh[i];
        }

        gaussSeidel(u, f, n, nu2, sigma, h);
    } else {
        u[0] = f[0]/(2+sigma*pow(h,2))*pow(h,2);
        //printf("Coarsest grid level, norm of u: %.16e\n", norm(u, n-1)*pow(h,0.5));
    }
}     

void muCycle(int mu, double* u, double* f, int n, int nu1, int nu2, double sigma, double h) {
    int i = 0;
    int j = 0;
    double sum = 0.0;
    double residual[n-1];
    double residual2h[n/2-1];
    double errorh[n-1];
    double error2h[n/2-1];
    if(n != 2) {
        gaussSeidel(u, f, n, nu1, sigma, h);

        for(i = 0; i < n-1; i++) {
            sum = 0;
            for(j = i - 1; j < i + 2; j++) {
                sum = sum + stencilValue(i, j, n, sigma, h) * u[j];
            }
            residual[i] = f[i] - sum;
        }

        restriction(residual, n, residual2h);

        for(i = 0; i < n/2-1; i++) {
            error2h[i] = 0;
        }

        if(n != 4) {
            for(i = 0; i < mu; i++) {
                muCycle(mu, error2h, residual2h, n/2, nu1, nu2, sigma, 2*h);
            }
        } else {
            muCycle(mu, error2h, residual2h, n/2, nu1, nu2, sigma, 2*h);
        }

        prolongation(error2h, n, errorh);

        for(i = 0; i < n-1; i++) {
            u[i] = u[i] + errorh[i];
        }

        gaussSeidel(u, f, n, nu2, sigma, h);
    } else {
        u[0] = f[0]/(2+sigma*pow(h,2))*pow(h,2);
    }
}

void FMG(double* u, double* f, int n, int nu1, int nu2, double sigma, double h) {
    double f2h[n/2-1];
    double u2h[n/2-1];
    int i = 0;
    if(n != 2) {
        restriction(f, n, f2h);
        for(i = 0; i < n/2-1; i++) {
            u2h[i] = 0;
        }
        FMG(u2h, f2h, n, nu1, nu2, sigma, h);
        prolongation(u2h, n, u);
    }
    vCycle(u, f, n, nu1, nu2, sigma, h);
}

int simulate(int k, int nu1, int nu2, double tolerance, FILE* fptr) {

    fprintf(fptr, "Simulating with k = %d\n", k);

    int n = (int)pow(2,k);
    double h = 1.0/n;
    double sigma = 0.0;
    double u[n-1];
    double f[n-1];
    int a = 0;
    int b = 0;
    double temp[n-1];
    double residual[n-1];
    double residualNorm = 0.0;
    double sum = 0;

    int i = 0;
    int j = 0;
    for(i = 0; i < n; i++) {
        u[i] = (double)(rand()%(10001))/10000.0; //sin(M_PI * 1 * i * h);
        f[i] = 0;
    }
    i = 0;
    
    while(i == 0 || residualNorm > tolerance) {
        //vCycle(u, f, n, nu1, nu2, sigma, h);
        //muCycle(2, u, f, n, nu1, nu2, sigma, h);
        FMG(u, f, n, nu1, nu2, sigma, );

        // ======= Residual ========
        for(a = 0; a < n-1; a++) {
            sum = 0.0;
            for(b = a - 1; b < a + 2; b++) {
                sum = sum + stencilValue(a, b, n, sigma, h) * u[b]; // Au
            }
            residual[a] = f[a] - sum; // r = f - Au
        }
        double fNorm = norm(f, n-1);
        if(fNorm != 0) {
            residualNorm = norm(residual, n-1)/norm(f, n-1); // ||r||/||f||
        } else {
            residualNorm = norm(residual, n-1);
        }
        //printf("After iteration %d, norm of the residual is %.16e\n", i+1, residualNorm);
        // ======= Residual ========

        i = i+1;
        // ==== Write output to file ====
        fprintf(fptr, "%d,", i);
        for(j = 0; j < n-2; j++) {
            fprintf(fptr, "%.10f,", u[j]);
        }
        fprintf(fptr, "%.10f\n", u[n-2]);

        printf("\nLOOP\n");
    }

    return i; // Return number of iterations
}

int main(int argc, char **argv) {
    FILE* fptr;
    fptr = fopen("solution.csv", "w");

    double tolerance = pow(10.0, -6);
    double nu1 = 1; // Number of Gauss-Seidel iterations before coarsening
    double nu2 = 1; // Number of Gauss Seidel iterations after coarsening

    if(argc == 4) {
        nu1 = atoi(argv[1]);
        nu2 = atoi(argv[2]);
        tolerance = pow(10, -1*atof(argv[3]));
    }

    int k = 4;
    for(k = 4; k <= 4; k++) {
        printf("k = %d -> Num iterations: %d\n", k, simulate(k, nu1, nu2, tolerance, fptr));
    }
    return 0;
}