#pragma once
#ifndef HEADER_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define HEADER_H

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <omp.h>

struct crsMatrix
{
    int N;
    int NZ;
    double* Value;
    int* Col;
    int* row_index;
};

struct regMatrix
{
    int N;
    double** Value;
};

struct Vec
{
    int N;
    double** Value;
};

void AllocateCRSMatrix(int N, int NZ, crsMatrix& mtx);
void AllocateRegMatrix(int N, regMatrix& mtx);
void AllocateVector(int N, Vec& vec);

void FreeRegMatrix(regMatrix& mtx);
void FreeCRSMatrix(crsMatrix& mtx);
void FreeVector(Vec& vec);

double next();

void GenerateRegularCRS(int seed, int N, int cntInRow, crsMatrix& mtx);
void GenerateSpecialCRS(int seed, int N, int cntInRow, crsMatrix& mtx);
void GenerateVector(int seed, int N, Vec& vec);

double Transpose2(crsMatrix imtx, crsMatrix& omtx);

int CompareVectors(Vec& vec1, Vec& vec2, double THRESHOLD, double& diff);
int CompareMatrix(regMatrix& reg1, regMatrix& reg2, double THRESHOLD, double& diff);

int SeqMult(crsMatrix A, crsMatrix B, crsMatrix& C, double& time);
int ParMult(crsMatrix A, crsMatrix B, crsMatrix& C, int ThN, double& time);
void SeqRegMult(regMatrix A, regMatrix B, regMatrix& C, double& time);
void ParRegMult(regMatrix A, regMatrix B, regMatrix& C, int ThN, double& time);


void RegVecMult(regMatrix A, Vec& iVec, Vec& oVec, double& time);
void CRSVecMult(crsMatrix A, Vec& iVec, Vec& oVec, double& time);
void ParRegVecMult(regMatrix A, Vec& iVec, Vec& oVec, int ThN, double& time);
void ParCRSVecMult(crsMatrix A, Vec& iVec, Vec& oVec, int ThN, double& time);


void CRStoReg(crsMatrix iM, regMatrix& oM);

void printReg(regMatrix& mtx);
void printCRS(crsMatrix& mtx);
void printVec(Vec& vec);
void printMult(crsMatrix& mtx, Vec& iVec, Vec& oVec);

#endif