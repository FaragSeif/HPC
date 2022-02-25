#include "Header.h"
#include "gnuplot-iostream.h"


using std::vector;
using std::cout;
using std::endl;
using namespace std::chrono;

int main() {
    system("echo ==================================== Seif Farag ==========================================");
    system("echo ================================== Assignment #1 ========================================= & echo.");
    system("echo CPU INFO:");
    system("echo ----------");
    system("wmic cpu get deviceid, name, numberofcores, maxclockspeed");
    system("echo RAM INFO:");
    system("echo ----------");
    system("wmic MemoryChip get BankLabel, Capacity, Speed");
    system("echo ===============================================================");

    Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");

    int N[] = { 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000 };

    int ThN[] = { 2,4,8 };
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 5000 / 500 / 2; i++) {
            crsMatrix cA, cB, cBT, cC, PcC;
            regMatrix rA, rB, rC, rcC, PrC;
            Vec iVec, oVec;
            int cntInRow;

            do {
                cntInRow = rand();
            } while (cntInRow >= N[i] / 3);

            AllocateVector(N[i], iVec);
            AllocateVector(N[i], oVec);
            GenerateVector(1, N[i], iVec);
            GenerateRegularCRS(1, N[i], cntInRow, cA);
            GenerateSpecialCRS(2, N[i], cntInRow, cB);
            //Transpose2(cB, cBT);
            double timeM1, timeM2, timeM3, timeM4, VtimeM1, VtimeM2, VtimeM3, VtimeM4;
            vector<vector<double>> v;

            //SeqMult(cA, cBT, cC, timeM1);
            //ParMult(cA, cBT, PcC, ThN[i], timeM2);

            CRStoReg(cA, rA);
            //CRStoReg(cB, rB);
            //CRStoReg(cC, rcC);

            //SeqRegMult(rA, rB, rC, timeM3);
            //ParRegMult(rA, rB, PrC, ThN[i],timeM4);

            //CRSVecMult(cA, iVec, oVec, VtimeM1);
            ParCRSVecMult(cA, iVec, oVec, ThN[j], VtimeM1);
            //RegVecMult(rA, iVec, oVec, VtimeM3);
            ParRegVecMult(rA, iVec, oVec, ThN[j], VtimeM4);
            cout << "( " << i << " , " << j << " )";
            //cout << "Initializing a CRS matrix with " << N << " rows and colums..." << endl;
            //cout << "===============================================================\n";
            //cout << "CRS Matrix Initialized.!!" << endl;
            //cout << "===============================================================\n";
            //cout << "Initializing a CRS matrix with " << N*2 << " rows and colums..." << endl;
            //cout << "===============================================================\n";
            //cout << "Cubic CRS Matrix Initialized.!!" << endl;
            //cout << "===============================================================";
            //cout << endl << "Converting CRS matrix to Regular..." << endl;
            //cout << "===============================================================";

            //double diff;
            //int code = CompareMatrix(rC, rcC, 0.00000002, diff);
            //if (code == -1) cout << "wrong dimension";
            //else if (code == 0) cout << "For Error threshold of 0.00000002 and N = " << N * i << "  ..INCORRECT.!!\n";
            //else cout << "For Error threshold of 0.00000002 and N = " << N * i << "  ..CORRECT.!!\n";

            /*double diff;
            int code = CompareVectors(iVec, oVec, 0.00000002, diff);
            if (code == -1) cout << "wrong dimension";
            else if (code == 0) cout << "For Error threshold of 0.00000002 and N = " << N*i << "  ..INCORRECT.!!\n";
            else cout << "For Error threshold of 0.00000002 and N = " << N*i << "  ..CORRECT.!!\n";*/
            //cout << endl << "Sequential CRS time: " << timeM1 << endl;
            //cout << endl << "Parallel CRS time: " << timeM3 << endl;
            //cout << endl << "Sequential Reg time: " << timeM2 << endl;
            //cout << endl << "Parallel Reg time: " << timeM2 << endl;


            /*FreeCRSMatrix(cA);
            FreeCRSMatrix(cB);
            FreeCRSMatrix(cC);
            FreeCRSMatrix(PcC);
            FreeRegMatrix(rA);
            FreeRegMatrix(rB);
            FreeRegMatrix(rC);
            FreeRegMatrix(PrC);
            FreeRegMatrix(rcC);
            FreeVector(iVec);
            FreeVector(oVec);*/
        }
    }
    //vector<vector<double>> v{
    //            { 1, 2},
    //            { 3, 3},
    //            { 5, 4},
    //            { 7, 5},
    //            { 9, 6},
    //            { 11, 7},
    //            { 13, 8},
    //            { 15, 10},
    //            { 17, 11},
    //            { 19, 15}
    //};
    //double v2[10][2] = { {1,2},{3,4},{3,6}, {3,8}, {3,10}, {3,12}, {3,14}, {3,16}, {3,18}, {3,20} };
    //gp << "plot '-' with lines\n";
    //gp.send(v2);
    //std::cin.get();
	return 0;
}