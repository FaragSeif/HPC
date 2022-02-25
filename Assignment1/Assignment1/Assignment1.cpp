// Assignment1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <iomanip>
#include "Header.h" 

using std::vector;
using std::cout;
using std::endl;
using namespace std::chrono;


//========================================================
//==================== Allocation ========================
//========================================================
void AllocateCRSMatrix(int N, int NZ, crsMatrix& mtx) {
    mtx.N = N;
    mtx.NZ = NZ;
    mtx.Value = new double[NZ]();
    mtx.Col = new int[NZ]();
    mtx.row_index = new int[N + 1]();
}
void AllocateRegMatrix(int N, regMatrix& mtx) {
    mtx.N = N;
    mtx.Value = new double* [N];
    for (int i = 0; i < N; i++) {
        mtx.Value[i] = new double[N]();
    }
}
void AllocateVector(int N, Vec& vec) {
    vec.N = N;
    vec.Value = new double* [N];
    for (int i = 0; i < N; i++) {
        vec.Value[i] = new double[1]();
    }
}

//========================================================
//===================== Freeing ==========================
//========================================================
void FreeRegMatrix(regMatrix& mtx) {
    for (int i = 0; i < mtx.N; ++i)
        delete[] mtx.Value[i];
    delete[] mtx.Value;
}
void FreeCRSMatrix(crsMatrix& mtx) {
    delete[] mtx.Value;
    delete[] mtx.Col;
    delete[] mtx.row_index;
}
void FreeVector(Vec& vec) {
    for (int i = 0; i < vec.N; ++i)
        delete[] vec.Value[i];
    delete[] vec.Value;
}

//========================================================
//==================== Generation ========================
//========================================================
void GenerateRegularCRS(int seed, int N, int cntInRow, crsMatrix& mtx)
{
    int i, j, k, f, tmp, notNull, c;

    srand(seed);

    notNull = cntInRow * N;
    AllocateCRSMatrix(N, notNull, mtx);

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < cntInRow; j++)
        {
            do
            {
                mtx.Col[i * cntInRow + j] = rand() % N;
                f = 0;
                for (k = 0; k < j; k++)
                    if (mtx.Col[i * cntInRow + j] == mtx.Col[i * cntInRow + k])
                        f = 1;
            } while (f == 1);
        }
        for (j = 0; j < cntInRow - 1; j++)
            for (k = 0; k < cntInRow - 1; k++)
                if (mtx.Col[i * cntInRow + k] > mtx.Col[i * cntInRow + k + 1])
                {
                    tmp = mtx.Col[i * cntInRow + k];
                    mtx.Col[i * cntInRow + k] = mtx.Col[i * cntInRow + k + 1];
                    mtx.Col[i * cntInRow + k + 1] = tmp;
                }
    }

    for (i = 0; i < cntInRow * N; i++)
        mtx.Value[i] = next(); //* MAX_VAL;

    c = 0;
    for (i = 0; i <= N; i++)
    {
        mtx.row_index[i] = c;
        c += cntInRow;
    }
}
void GenerateSpecialCRS(int seed, int N, int cntInRow, crsMatrix& mtx)
{
    srand(seed);
    double end = pow((double)cntInRow, 1.0 / 3.0);
    double step = end / N;

    vector<int>* columns = new vector<int>[N];
    int NZ = 0;

    for (int i = 0; i < N; i++)
    {
        int rowNZ = int(pow((double(i + 1) * step), 3) + 1);
        NZ += rowNZ;
        int num1 = (rowNZ - 1) / 2;
        int num2 = rowNZ - 1 - num1;

        if (rowNZ != 0)
        {
            if (i < num1)
            {
                num2 += num1 - i;
                num1 = i;
                for (int j = 0; j < i; j++)
                    columns[i].push_back(j);
                columns[i].push_back(i);
                for (int j = 0; j < num2; j++)
                    columns[i].push_back(i + 1 + j);
            }
            else
            {
                if (N - i - 1 < num2)
                {
                    num1 += num2 - (N - 1 - i);
                    num2 = N - i - 1;
                }
                for (int j = 0; j < num1; j++)
                    columns[i].push_back(i - num1 + j);
                columns[i].push_back(i);
                for (int j = 0; j < num2; j++)
                    columns[i].push_back(i + j + 1);
            }
        }
    }

    AllocateCRSMatrix(N, NZ, mtx);

    int count = 0;
    int sum = 0;
    for (int i = 0; i < N; i++)
    {
        mtx.row_index[i] = sum;
        sum += columns[i].size();
        for (unsigned int j = 0; j < columns[i].size(); j++)
        {
            mtx.Col[count] = columns[i][j];
            mtx.Value[count] = next();
            count++;
        }
    }
    mtx.row_index[N] = sum;

    delete[] columns;
}
void GenerateVector(int seed, int N, Vec& vec) {
    srand(seed);
    for (int i = 0; i < vec.N; i++)
        vec.Value[i][0] = next();
}


//========================================================
//================== Multiplication ======================
//========================================================
int SeqMult(crsMatrix A, crsMatrix B, crsMatrix& C, double& time)
{
    if (A.N != B.N)
        return 1;

    int N = A.N;
    vector<int> columns;
    vector<double> values;
    vector<int> row_index;
    int ZERO_IN_CRS = 0;

    auto start = high_resolution_clock::now();
    int rowNZ;
    row_index.push_back(0);
    for (int i = 0; i < N; i++)
    {
        rowNZ = 0;
        for (int j = 0; j < N; j++)
        {
            double sum = 0;
            for (int k = A.row_index[i]; k < A.row_index[i + 1]; k++)
            {
                for (int l = B.row_index[j]; l < B.row_index[j + 1]; l++)
                {
                    if (A.Col[k] == B.Col[l])
                    {
                        sum += A.Value[k] * B.Value[l];
                        break;
                    }
                }
            }
            if (fabs(sum) > ZERO_IN_CRS)
            {
                columns.push_back(j);
                values.push_back(sum);
                rowNZ++;
            }
        }
        row_index.push_back(rowNZ + row_index[i]);
    }

    AllocateCRSMatrix(N, columns.size(), C);

    for (unsigned int j = 0; j < columns.size(); j++)
    {
        C.Col[j] = columns[j];
        C.Value[j] = values[j];
    }
    for (int i = 0; i <= N; i++)
        C.row_index[i] = row_index[i];

    auto finish = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(finish - start);
    time = duration.count();

    return 0;
}
void SeqRegMult(regMatrix A, regMatrix B, regMatrix& C, double& time) {
    AllocateRegMatrix(A.N, C);
    auto start = high_resolution_clock::now();
    for (int i = 0; i < A.N; i++)
        for (int j = 0; j < A.N; j++)
            for (int k = 0; k < A.N; k++)
                C.Value[j][i] += A.Value[j][k] * B.Value[k][i];
    auto finish = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(finish - start);
    time = duration.count();
}

void RegVecMult(regMatrix A, Vec& iVec, Vec& oVec, double& time) {
    auto start = high_resolution_clock::now();
    for (int i = 0; i < A.N; i++)
    {
        // Computing the i-th component of the vector y
        double sum = 0;
        for (int j = 0; j < A.N; j++)
            sum = sum + A.Value[i][j] * iVec.Value[j][0];
        oVec.Value[i][0] = sum;
    }
    auto finish = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(finish - start);
    time = duration.count();
}
void CRSVecMult(crsMatrix A, Vec& iVec, Vec& oVec, double& time) {
    auto start = high_resolution_clock::now();
    for (int i = 0; i < A.N; i++)
    {
        // Computing the i-th component of the vector y
        double sum = 0;
        int j1 = A.row_index[i];
        int j2 = A.row_index[i + 1];
        for (int j = j1; j < j2; j++)
            sum = sum + A.Value[j] * iVec.Value[A.Col[j]][0];
        oVec.Value[i][0] = sum;
    }
    auto finish = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(finish - start);
    time = duration.count();
}

//==================== Parallel ==========================
int ParMult(crsMatrix A, crsMatrix B, crsMatrix& C, int& ThN, double& time)
{
    if (A.N != B.N)
        return 1;

    int N = A.N;
    vector<int> columns;
    vector<double> values;
    vector<int> row_index;
    int ZERO_IN_CRS = 0;

    auto start = high_resolution_clock::now();
    int rowNZ;
    row_index.push_back(0);
    omp_set_num_threads(ThN);
    for (int i = 0; i < N; i++)
    {
        rowNZ = 0;
        #pragma omp parallel for
        for (int j = 0; j < N; j++)
        {
            double sum = 0;
            for (int k = A.row_index[i]; k < A.row_index[i + 1]; k++)
            {
                for (int l = B.row_index[j]; l < B.row_index[j + 1]; l++)
                {
                    if (A.Col[k] == B.Col[l])
                    {
                        sum += A.Value[k] * B.Value[l];
                        break;
                    }
                }
            }
            if (fabs(sum) > ZERO_IN_CRS)
            {
                columns.push_back(j);
                values.push_back(sum);
                rowNZ++;
            }
        }
        row_index.push_back(rowNZ + row_index[i]);
    }

    AllocateCRSMatrix(N, columns.size(), C);

    for (unsigned int j = 0; j < columns.size(); j++)
    {
        C.Col[j] = columns[j];
        C.Value[j] = values[j];
    }
    for (int i = 0; i <= N; i++)
        C.row_index[i] = row_index[i];

    auto finish = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(finish - start);
    time = duration.count();

    return 0;
}
void ParRegMult(regMatrix A, regMatrix B, regMatrix& C, int ThN, double& time) {
    AllocateRegMatrix(A.N, C);
    auto start = high_resolution_clock::now();
    omp_set_num_threads(ThN);
    #pragma omp parallel for
    for (int i = 0; i < A.N; i++)
        for (int j = 0; j < A.N; j++)
            for (int k = 0; k < A.N; k++)
                C.Value[j][i] += A.Value[j][k] * B.Value[k][i];
    auto finish = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(finish - start);
    time = duration.count();
}

void ParRegVecMult(regMatrix A, Vec& iVec, Vec& oVec, int ThN, double& time) {
    auto start = high_resolution_clock::now();
    //int n = ThN;
    omp_set_num_threads(ThN);
    #pragma omp parallel for
    for (int i = 0; i < A.N; i++)
    {
        // Computing the i-th component of the vector y
        double sum = 0;
        for (int j = 0; j < A.N; j++)
            sum = sum + A.Value[i][j] * iVec.Value[j][0];
        oVec.Value[i][0] = sum;
    }
    auto finish = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(finish - start);
    time = duration.count();
}
void ParCRSVecMult(crsMatrix A, Vec& iVec, Vec& oVec, int ThN, double& time) 
{
    int n = ThN;
    omp_set_num_threads(ThN);
    auto start = high_resolution_clock::now();
    #pragma omp parallel for
    for (int i = 0; i < A.N; i++)
    {
        // Computing the i-th component of the vector y
        double sum = 0;
        int j1 = A.row_index[i];
        int j2 = A.row_index[i + 1];
        for (int j = j1; j < j2; j++)
            sum = sum + A.Value[j] * iVec.Value[A.Col[j]][0];
        oVec.Value[i][0] = sum;
    }
    auto finish = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(finish - start);
    time = duration.count();
}


//========================================================
//==================== UTILITIES =========================
//========================================================
double next()
{
    return ((double)rand() / (double)RAND_MAX);
}
void CRStoReg(crsMatrix iM, regMatrix& oM) {
    AllocateRegMatrix(iM.N, oM);
    for (int i = 0; i < iM.N; i++) {
        for (int j = iM.row_index[i]; j < iM.row_index[i + 1]; j++) {
            oM.Value[i][iM.Col[j]] = iM.Value[j];
        }
    }
}
int CompareVectors(Vec& vec1, Vec& vec2, double THRESHOLD, double& diff)
{
    diff = 0.0;
    if (vec1.N != vec2.N)
        return -1;
    for (int i = 0; i < vec1.N; i++)
    {
        if (diff < fabs(vec1.Value[i][0] - vec2.Value[i][0]))
        {
            diff = fabs(vec1.Value[i][0] - vec2.Value[i][0]);
        }
    }
    if (diff > THRESHOLD) return 1;
    else return 0;
}
int CompareMatrix(regMatrix& reg1, regMatrix& reg2, double THRESHOLD, double& diff)
{
    diff = 0.0;
    if (reg1.N != reg2.N)
        return -1;
    for (int i = 0; i < reg1.N; i++)
    {
        for (int j = 0; j < reg1.N; j++) {
            double test = fabs(reg1.Value[i][j] - reg2.Value[i][j]);
            //cout << test << endl;
            if (diff < fabs(reg1.Value[i][j] - reg2.Value[i][j]))
            {
                diff = fabs(reg1.Value[i][j] - reg2.Value[i][j]);
            }
        }
    }
    if (diff <= THRESHOLD) return 1;
    else return 0;
}
double Transpose2(crsMatrix imtx, crsMatrix& omtx)
{
    int i, j;

    auto start = high_resolution_clock::now();

    AllocateCRSMatrix(imtx.N, imtx.NZ, omtx);

    memset(omtx.row_index, 0, (imtx.N + 1) * sizeof(int));
    for (i = 0; i < imtx.NZ; i++)
        omtx.row_index[imtx.Col[i] + 1]++;

    int S = 0;
    for (i = 1; i <= imtx.N; i++)
    {
        int tmp = omtx.row_index[i];
        omtx.row_index[i] = S;
        S = S + tmp;
    }

    for (i = 0; i < imtx.N; i++)
    {
        int j1 = imtx.row_index[i];
        int j2 = imtx.row_index[i + 1];
        int Col = i;
        for (j = j1; j < j2; j++)
        {
            double V = imtx.Value[j];
            int RIndex = imtx.Col[j];
            int IIndex = omtx.row_index[RIndex + 1];
            omtx.Value[IIndex] = V;
            omtx.Col[IIndex] = Col;
            omtx.row_index[RIndex + 1]++;
        }
    }

    auto finish = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(finish - start);
    return duration.count();
}


//========================================================
//===================== Printing =========================
//========================================================
void printReg(regMatrix& mtx) {
    cout << "\n--------------------------------- In Regualr Form ------------------------------------" << endl;
    for (int i = 0; i < mtx.N; i++)
    {
        for (int j = 0; j < mtx.N; j++)
            cout << std::setw(12) << mtx.Value[i][j] << "  ";
        cout << endl;
    }
    cout << "--------------------------------------------------------------------------------------" << endl;
}
void printCRS(crsMatrix& mtx) {
    for (int i = 0; i < mtx.NZ - 1; i++)
    {
        if (i == 0) cout << "Value Matrix: [ ";
        cout << mtx.Value[i] << " , ";
    }
    cout << mtx.Value[mtx.NZ - 1] << " ]" << endl;


    for (int i = 0; i < mtx.NZ - 1; i++)
    {
        if (i == 0) cout << "Column Matrix: [ ";
        cout << mtx.Col[i] << " , ";
    }
    cout << mtx.Col[mtx.NZ - 1] << " ]" << endl;


    for (int i = 0; i < mtx.N; i++)
    {
        if (i == 0) cout << "Row Idx Matrix: [ ";
        cout << mtx.row_index[i] << " , ";
    }
    cout << mtx.row_index[mtx.N] << " ]" << endl;
    regMatrix regMtx;
    CRStoReg(mtx, regMtx);
    printReg(regMtx);
}
void printVec(Vec& vec) {
    for (int i = 0; i < vec.N; i++)
        cout << vec.Value[i][0] << "  ";
    cout << endl;
}
void printMult(crsMatrix& mtx, Vec& iVec, Vec& oVec) {
    regMatrix regMtx;
    CRStoReg(mtx, regMtx);
    for (int i = 0; i < regMtx.N; i++)
    {
        cout << "|";
        for (int j = 0; j < regMtx.N; j++)
            cout << std::setw(10) << regMtx.Value[i][j] << "  ";
        cout << "|" << std::setw(10) << "|" << std::setw(10) << iVec.Value[i][0] << std::setw(1) << "|";
            
        cout << std::setw(10) << "|" << std::setw(10) << oVec.Value[i][0] << std::setw(1) << "|" << endl;
    }
}
//=======================================================