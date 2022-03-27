#include "math.h"
#include <omp.h>
#include <cstdio>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <cstdio>



#define RIGHT_BOUND 1.0
#define LEFT_BOUND 1.0
#define N_MAX 100000
#define _USE_MATH_DEFINES
#define EPSILON 0.0000000001
double f(double x, double y) {
	return -4;
}
double mu1(double y) {
	return (1 + pow(y+2, 2));
}
double mu2(double y) {
	return (4 + pow(y+2, 2));
}
double mu3(double x) {
	return (pow(x+1, 2) + 4);
}
double mu4(double x) {
	return (pow(x+1, 2) + 9);
}
double FuncU(double x, double y) {
	return (pow(x+1, 2) + pow(y+2, 2));
}

void plotExactSolution(int N) {
	double h = LEFT_BOUND / (double)N;
	double k = RIGHT_BOUND / (double)N;
	FILE* gnuplot = _popen("gnuplot -persist", "w");
	fprintf(gnuplot, "set view map\n");
	fprintf(gnuplot, "set dgrid3d\n");
	fprintf(gnuplot, "set pm3d interpolate 0,0\n");
	fprintf(gnuplot, "set title 'Heatmap'\n");
	fprintf(gnuplot, "set xrange [-0.001:1.001]\n");
	fprintf(gnuplot, "splot '-' using 1:2:3 with pm3d\n");
	for (int i = 0; i <= N; ++i) {
		double x = h * i;
		for (int j = 0; j <= N; ++j) {
			double y = k * j;
			//printf("(%f, %f, %f)\n", x, y + 2, FuncU(x, y));
			fprintf(gnuplot, "%f %f %f\n", x, y+2, FuncU(x, y));
		}
	}
	fprintf(gnuplot, "e\n");
	_pclose(gnuplot);
}

// Vector allocation
//void AllocateVector(double** Vector, int size)
//{
//	(*Vector) = new double[size];
//}
// Vector release
//void FreeVector(double** Vector)
//{
//	delete[](*Vector);
//}
// Vector printing
//void PrintVector(double* Vector, double* Vector2,  int size)
//{
//	//printf("%s\n", VectorName);
//	for (int i = 0; i < size; i++){
//		printf("%f", Vector[i]);
//		printf("      ,%f", Vector2[i]);
//		printf("\n");
//	}
//}
//--------------------------------------------------------------------------
// Initialization of the data by the partial differential equation
//--------------------------------------------------------------------------
// matrix allocation (n, m)
// Diagonals of the matrix are stored by rows
// the shift of the diagonal elements is stored in the index array
void CreateDUMatrix(int n, int m, double** Matrix, double** Index)
{
	double hsqr = (double)n * n / LEFT_BOUND / LEFT_BOUND; // 1/h
	double ksqr = (double)m * m / RIGHT_BOUND / RIGHT_BOUND; // 1/k
	double A = 2 * (hsqr + ksqr);
	int size = (n - 1) * (m - 1), bandWidth = 5;
	AllocateVector(Matrix, size * bandWidth);
	AllocateVector(Index, bandWidth);
	(*Index)[0] = -n + 1; (*Index)[1] = -1; (*Index)[2] = 0; (*Index)[3] = 1; (*Index)[4] = n - 1;
	for (int i = 0; i < size; i++)
	{
		if (i >= n - 1) (*Matrix)[i * bandWidth] = -ksqr;
		else (*Matrix)[i * bandWidth] = 0.0;
		if (i % (n - 1) != 0) (*Matrix)[i * bandWidth + 1] = -hsqr;
		else (*Matrix)[i * bandWidth + 1] = 0.0;
		(*Matrix)[i * bandWidth + 2] = A;
		if ((i + 1) % (n - 1) != 0) (*Matrix)[i * bandWidth + 3] = -hsqr;
		else (*Matrix)[i * bandWidth + 3] = 0.0;
		if (i < (n - 1) * (m - 2)) (*Matrix)[i * bandWidth + 4] = -ksqr;
		else (*Matrix)[i * bandWidth + 4] = 0.0;
	}
}


// The right-hand side of the equation
void CreateDUVector(int n, int m, double** Vector)
{
	double h = LEFT_BOUND / (double)n;
	double k = RIGHT_BOUND / (double)m;
	double hsqr = (double)n * n / LEFT_BOUND / LEFT_BOUND;
	double ksqr = (double)m * m / RIGHT_BOUND / RIGHT_BOUND;
	AllocateVector(Vector, (n - 1) * (m - 1));
	for (int j = 0; j < m - 1; j++)
	{
		for (int i = 0; i < n - 1; i++)
			(*Vector)[j * (n - 1) + i] = f((double)(i) * h, (double)(j) * k);
		(*Vector)[j * (n - 1)] += hsqr * mu1((double)(j) * k);
		(*Vector)[j * (n - 1) + n - 2] += hsqr * mu2((double)(j) * k);
	}
	for (int i = 0; i < n - 1; i++)
	{
		(*Vector)[i] += ksqr * mu3((double)(i) * h);
		(*Vector)[(m - 2) * (n - 1) + i] += ksqr * mu4((double)(i) * h);
	}
}


void GetFirstApproximation(double* Result, int size)
{
	for (int i = 0; i < size; i++)
		Result[i] = 0.0;
}


//double GetWParam(double Step);

double BandOverRelaxation(double* Matrix, double* Vector, double* Result, double* Index, int size, int bandWidth, double WParam, double Accuracy, int& StepCount)
{
	double CurrError;//achieved accuracy
	double sum, TempError;
	int ii, index = Index[bandWidth - 1], bandHalf = (bandWidth - 1) / 2;
	StepCount = 0;
	do
	{
		CurrError = -1.0;
		for (int i = index; i < size + index; i++)
		{
			ii = i - index;
			TempError = Result[i];
			sum = 0.0;
			for (int j = 0; j < bandWidth; j++) {
				int myidx = i + Index[j];
				sum += Matrix[ii * bandWidth + j] * Result[myidx];
			}
			Result[i] = (Vector[ii] - sum) * WParam / Matrix[ii * bandWidth + bandHalf] + Result[i];
			TempError = fabs(Result[i] - TempError);
			if (TempError > CurrError) CurrError = TempError;
		}
		StepCount++;
	} while ((CurrError > Accuracy) && (StepCount < N_MAX));
	return CurrError;
}


double SolvePoisson(int n, int m, double* Solution, double Accuracy, double& ORAccuracy, int& StepCount)
{
	double* Matrix, * Vector, * Result;
	double* Index;
	int size = (n - 1) * (m - 1), ResSize = size + 2 * (n - 1), bandWidth = 5;
	double start, finish; double time;
	double WParam, step = (n / LEFT_BOUND > m / RIGHT_BOUND) ? (double)LEFT_BOUND / n : (double)RIGHT_BOUND / m;
	CreateDUMatrix(n, m, &Matrix, &Index);
	CreateDUVector(n, m, &Vector);
	AllocateVector(&Result, ResSize);
	GetFirstApproximation(Result, ResSize);
	WParam = 1;//GetWParam(step);
	start = omp_get_wtime();
	ORAccuracy = BandOverRelaxation(Matrix, Vector, Result, Index, size, bandWidth, WParam, Accuracy, StepCount);
	finish = omp_get_wtime();
	time = (finish - start);
	memcpy(Solution, Result + n - 1, sizeof(double) * size);
	FreeVector(&Matrix);
	FreeVector(&Index);
	FreeVector(&Vector);
	FreeVector(&Result);
	return time;
}


double SolutionCheck(double* solution, int n, int m)
{
	double h = LEFT_BOUND / (double)n, k = RIGHT_BOUND / (double)m;
	double err = 0, temp;
	double* temp2;
	AllocateVector(&temp2, (n - 1) * (m - 1));
	for (int j = 0; j < m - 1; j++)
		for (int i = 0; i < n - 1; i++)
		{
			temp2[j * (n - 1) + i] = FuncU((double)(i) * h, (double)(j) * k);
			temp = fabs(solution[j * (n - 1) + i] - FuncU((double)(i) * h, (double)(j) * k));
			if (temp > err)
				err = temp;
		}
	//PrintVector(solution, temp2, (n - 1) * (m - 1));
	return err;
}



//int main(int argc, char* argv[])
//{
	//int n[5] = { 10, 50, 100, 500, 1000 };
	//int m[5] = { 10, 50, 100, 500, 1000 };
	//double Epsilon[3] = { 0.001, 0.0001, 0.00001 };
	////printf("|%-15s|%-20s|%-15s|%-20s|\n", "N & M ", "EPSILON", "Iterations", "Achieved Epsilon");
	////printf("%-25s%-20s%-10.2f%-10.2f%-10.2f\n", name.c_str(), title.c_str(), gross, tax, net);
	////for (int j = 0; j < 3; j++) {
	//printf("|  N & M  |  EPSILON  |  Iterations  |   Achieved Epsilon   |\n");
	//	//for (int i = 0; i < 5; i++) {
	//		int StepCount;
	//		int size;
	//		double time;
	//		double Accuracy = Epsilon[1];
	//		double AcAccuracy;
	//		double Correctness;
	//		double* Solution;
	//		size = (n[1] - 1) * (m[1] - 1);
	//		AllocateVector(&Solution, size);
	//		time = SolvePoisson(n[1], m[1], Solution, Accuracy, AcAccuracy, StepCount);

	//		//printf("UpperRelaxation:\ntime = %.15f\n", time);
	//		//printf("Accuracy = %.15f, stepCount = %d\n", AcAccuracy, StepCount);
	//		Correctness = SolutionCheck(Solution, n[1], m[1]);
	//		//printf("Exact and Approximate solution comparison = %.15f\n", Correctness);

	//		printf("| %5i   | %5f  | %7i      | %15f      |\n", n[1], Epsilon[1], StepCount, AcAccuracy);
	//		//printf(file, "%d;%d;%.15f;%.15f;%.15f;%d;%.15f;\n", n, m, Accuracy, AcAccuracy,StepCount, time);
	//		FreeVector(&Solution);

	//		plotExactSolution(n[1]);

	////	}
	////}

	//return 0;
//}