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
using namespace std;

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

	int const MatrixSize = ((N + 1)*(N + 2)) / 2;
    MonteCarlo monte;
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

        //liczymy wektor dla Monte Carlo
        /*for(int i=0;i<=N;i++){
            for (int j = 0; j <= N; j++) {
                if (i + j <= N) {
                    wektorMonte[k]=monte.monteCarlo(10000,i, j, N);
                    k++;
                    if(k==45) cout<< "wieksze pol" <<endl;
                    //cout<< k <<endl;
                }
            }
        }

        cout << "drugi" << endl;
        k=0;
        for(int i=0;i<=N;i++){
            for (int j = 0; j <= N; j++) {
                if (i + j <= N) {
                    wektorMonte2[k]=monte.monteCarlo(50000,i, j, N);

                    //cout << wektorMonte[k] << endl;
                    k++;
                    if(k==45) cout<< "wieksze pol" <<endl;
                    //if(k==15) cout<< "ponadpol" <<endl;
                }
            }
        }
        cout << "trzeci" << endl;
        k=0;
        for(int i=0;i<=N;i++){
            for (int j = 0; j <= N; j++) {
                if (i + j <= N) {
                    wektorMonte3[k]=monte.monteCarlo(200000,i, j, N);

                    //cout << wektorMonte[k] << endl;
                    k++;
                    if(k==45) cout<< "wieksze pol" <<endl;
                    //if(k==15) cout<< "ponadpol" <<endl;
                }
            }
        }
        cout << "czwarty" << endl;
        k=0;
        for(int i=0;i<=N;i++){
            for (int j = 0; j <= N; j++) {
                if (i + j <= N) {
                    wektorMonte4[k]=monte.monteCarlo(1000000,i, j, N);

                    //cout << wektorMonte[k] << endl;
                    k++;
                    if(k==45) cout<< "wieksze pol" <<endl;
                }
            }
        }*/

        //generujemy macierze
        MyMatrix macierz(gen.generujMacierz(N),MatrixSize);
        MyMatrix macierzJ(gen.generujMacierz(N),MatrixSize);
        MyMatrix macierzS(gen.generujMacierz(N),MatrixSize);

        //wyliczamy wektory i macierze
        //wektorWynikowy = macierz.czesc_wybor(wektorZerowy,MatrixSize);
        wektorWynikowy = macierz.czesc_wybor_ulepszony(wektorZerowy,MatrixSize);
        wektorWynikowyJ = macierzJ.jacobiIteracje(wektorZerowyJ,MatrixSize,500);
        wektorWynikowyS = macierzS.GSeidelIteracje(wektorZerowyS,MatrixSize,500);


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
        bladJ = liczenie_normy(odejmowanie_wektorow(wektorWynikowy,wektorWynikowyJ,MatrixSize),MatrixSize);
        bladS = liczenie_normy(odejmowanie_wektorow(wektorWynikowy,wektorWynikowyS,MatrixSize),MatrixSize);
        std::cout << std::setprecision(12) << bladJ << " bj" << endl;
        std::cout << std::setprecision(12) << bladS << " bs" << endl;
        /*bladwybory = liczenie_normy(odejmowanie_wektorow(wektorWynikowy,wektorMonte,MatrixSize),MatrixSize);
        cout << "1000;" << bladwybory << endl;
        bladwybory += liczenie_normy(odejmowanie_wektorow(wektorWynikowy,wektorMonte2,MatrixSize),MatrixSize);
        cout << "50000;" << bladwybory << endl;
        bladwybory += liczenie_normy(odejmowanie_wektorow(wektorWynikowy,wektorMonte3,MatrixSize),MatrixSize);
        cout << "200000;" << bladwybory << endl;
        bladwybory += liczenie_normy(odejmowanie_wektorow(wektorWynikowy,wektorMonte4,MatrixSize),MatrixSize);
        cout << "1000000;" << bladwybory << endl;*/




	}



	return 0;
}
