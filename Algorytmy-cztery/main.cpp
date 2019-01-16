#include <iostream>
#include <algorithm>
#include <vector>
#include "MyMatrix.h"
#include "Generator.h"
#include "MonteCarlo.h"
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include "Eigen/Dense"
#include "Eigen/SparseLU"
#include <Eigen/Sparse>
using namespace std;
using namespace Eigen;
std::vector<double> odejmowanie_wektorow(std::vector<double> B, std::vector<double> X, int const MatrixSize){
    std::vector<double> wynik(MatrixSize);
    for (int i =0 ; i< MatrixSize ; i++){
        wynik[i] = abs(B[i] - X[i]);
    }
    return wynik;
}

double liczenie_normy(std::vector<double> X, int const MatrixSize){
    double suma=0;
    double norma=0;
    for (int i =0 ; i< MatrixSize ; i++){
        suma += (double)(X[i] * X[i]);
    }
    norma = (double)(sqrt(suma) / MatrixSize);
    //cout << norma << endl;
    return norma;
}

int main()
{
	double d=0;
    string wyniki = "bladJ;bladS\n";

	int const MatrixSize = ((NOTN + 1)*(NOTN + 2)) / 2;
	vector<double> wektorZerowy(MatrixSize);
	vector<double> wektorZerowyJ(MatrixSize);
	vector<double> wektorZerowyS(MatrixSize);
	vector<double> wektorWynikowy(MatrixSize);
	vector<double> wektorWynikowyJ(MatrixSize);
	vector<double> wektorWynikowyS(MatrixSize);
	/*vector<double> wektorMonte(MatrixSize);
	vector<double> wektorMonte2(MatrixSize);
	vector<double> wektorMonte3(MatrixSize);
	vector<double> wektorMonte4(MatrixSize);*/
	cout<<"b\n";
	Generator gen;

	double bladJ, bladS, bladwybory;
    std::chrono::duration<double> start,stop,roznica_czasow, srednia;


	//---------------------------testowanie-----------------------------------//


    int k=0;
	for (int s = 1 ; s < 2 ; s ++){

        //generujemy puste wektory i macierze
        for (int i = 0; i < MatrixSize - 1; i++) {
		wektorZerowy[i] = 0;
        }
        wektorZerowy[MatrixSize - 1] = -1;
        for (int i = 0; i < MatrixSize - 1; i++) {
            wektorZerowyJ[i] = 0;
        }
        wektorZerowyJ[MatrixSize - 1] = -1;
        for (int i = 0; i < MatrixSize - 1; i++) {
            wektorZerowyS[i] = 0;
        }
        wektorZerowyS[MatrixSize - 1] = -1;


        //generujemy macierze
        MyMatrix macierz(gen.generujMacierz(NOTN),MatrixSize);
        MyMatrix macierzJ(gen.generujMacierz(NOTN),MatrixSize);
        MyMatrix macierzS(gen.generujMacierz(NOTN),MatrixSize);

        //wyliczamy wektory i macierze
        //wektorWynikowy = macierz.czesc_wybor(wektorZerowy,MatrixSize);
        wektorWynikowy = macierz.czesc_wybor_ulepszony(wektorZerowy,MatrixSize);
        //wektorWynikowyJ = macierzJ.jacobiIteracje(wektorZerowyJ,MatrixSize,500);
        //wektorWynikowyS = macierzS.GSeidelIteracje(wektorZerowyS,MatrixSize,500);


        //liczenie roznicy czasow pomiedzy gausami
        /*auto start = std::chrono::system_clock::now();
        wektorWynikowy = macierz.czesc_wybor(wektorZerowy,MatrixSize);
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        cout << "czas " << elapsed_seconds.count() << endl;

        auto start_u = std::chrono::system_clock::now();
        wektorWynikowyJ = macierzJ.czesc_wybor_ulepszony(wektorZerowyJ,MatrixSize);
        auto end_u = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds_u = end_u-start_u;
        cout << "czas_u " << elapsed_seconds_u.count() << endl;

        roznica_czasow += elapsed_seconds - elapsed_seconds_u;
*/
        //liczymy normy
        //bladJ = liczenie_normy(odejmowanie_wektorow(wektorWynikowy,wektorWynikowyJ,MatrixSize),MatrixSize);
        //bladS = liczenie_normy(odejmowanie_wektorow(wektorWynikowy,wektorWynikowyS,MatrixSize),MatrixSize);
        //std::cout << std::setprecision(12) << bladJ << " bj" << endl;
        //std::cout << std::setprecision(12) << bladS << " bs" << endl;
        /*bladwybory = liczenie_normy(odejmowanie_wektorow(wektorWynikowy,wektorMonte,MatrixSize),MatrixSize);
        cout << "1000;" << bladwybory << endl;
        bladwybory += liczenie_normy(odejmowanie_wektorow(wektorWynikowy,wektorMonte2,MatrixSize),MatrixSize);
        cout << "50000;" << bladwybory << endl;
        bladwybory += liczenie_normy(odejmowanie_wektorow(wektorWynikowy,wektorMonte3,MatrixSize),MatrixSize);
        cout << "200000;" << bladwybory << endl;
        bladwybory += liczenie_normy(odejmowanie_wektorow(wektorWynikowy,wektorMonte4,MatrixSize),MatrixSize);
        cout << "1000000;" << bladwybory << endl;*/

	}
SparseMatrix<double> sm1(MatrixSize,MatrixSize);
SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<Index> >   solver;
MyMatrix EigenMatrix(gen.generujMacierz(NOTN),MatrixSize);
//EigenMatrix.WyswietlMacierz(MatrixSize);
for(int i=0;i<MatrixSize;i++){
    for(int j=0;j<MatrixSize;j++){
            sm1.insert(i,j)=EigenMatrix.tab[i][j];
    }
}
//std::cout << sm1 << endl;
        VectorXd x(MatrixSize), b(MatrixSize);
                for (int i = 0; i < MatrixSize - 1; i++) {
		b(i)=0;
        }

        b(MatrixSize - 1) = -1;
// fill A and b;
// Compute the ordering permutation vector from the structural pattern of A
solver.compute(sm1);
if(solver.info()!=Success){
    return -1;
}
//Use the factors to solve the linear system
  x = solver.solve(b);

	return 0;
}
