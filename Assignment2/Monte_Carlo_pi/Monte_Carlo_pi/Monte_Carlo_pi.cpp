#include <iostream>
#include <omp.h>
#include <cstdio>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <cstdio>
#define _USE_MATH_DEFINES
#include <math.h>
#define L 2
#define N 10
#define T 2



void AllocateVector(double** Vector, int size)
{
    (*Vector) = new double[size];
}

void FreeVector(double** Vector)
{
    delete[](*Vector);
}

void PrintVector(double* Vector, int size)
{
    for (int i = 0; i < size; i++) {
        printf("%f", Vector[i]);
        printf("\n");
    }
}


double mu(double t) { // at x = 0 or x = N
    return 0;
}

double u0(double x) {  // at t = 0 
    return 0.1 * sin(M_PI * x);
}

double alpha(double x) {
    return exp(cos(x)) / 10;
}

void Solve(double* solution, double tau) {
    double x;
    double t;
    double h = (double)1 / N;

    for (int j = 0; j < L; j++) {
        solution[j * N] = mu(j);
        solution[((j + 1) * N)-1] = mu(j);
    }
    for (int i = 1; i < N-1; i++) {
        solution[i] = u0(i*h);
        solution[N + i] = u0(i*h);
    }
    PrintVector(solution, N * L);
    printf("=============================\n\n");
    for (int j = 1; j < L-1; j++) {

        for (int i = 1; i < N - 1; i++) {
            x = i * h;
            t = j * tau;
            solution[(j + 1) * N + i] = 2 * solution[N * j + i] - solution[N * (j - 1) + i] + (tau * tau) * alpha(x) * (1 / (h*h)) * 
                (solution[N * j + i - 1] - 2 * solution[N * j + i] + solution[N * j + i + 1]);
        }
    }
}

int main()
{
    double tau = (double)T/L;
    double* solution;
    AllocateVector(&solution, N*L);
    Solve(solution, tau);
    PrintVector(solution, N * L);
    return 0;
}