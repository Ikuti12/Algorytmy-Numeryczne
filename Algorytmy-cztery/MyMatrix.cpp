#include <iostream>
#include <cmath>
#include <vector>
#include "MyMatrix.h"

#define eps 0.00000000000001

MyMatrix::MyMatrix(double** arr, int const MatrixSize)
{
    for (int i = 0; i < MatrixSize; i++)
    {
        tab.push_back(std::vector<double>(MatrixSize));
    }
    for (int i = 0; i < MatrixSize; i++)
    {
        for (int j = 0; j < MatrixSize; j++)
            tab[i][j] = arr[i][j];
    }
}

MyMatrix::~MyMatrix() {}

void MyMatrix::WyswietlMacierz(int const MatrixSize)
{
    for (int i = 0; i < MatrixSize; i++)
    {
        for (int j = 0; j < MatrixSize; j++)
            std::cout << tab[i][j] << " ";

		std::cout << std::endl;
    }
}

std::vector<double>& MyMatrix::operator[](int i) { return tab[i]; }


void MyMatrix::zamien_wiersz(int i, int k, int j)
{
    double c;
    for (int z = j; z < N; z++)
    {
        c = tab[i][z];
        tab[i][z] = tab[k][z];
        tab[k][z] = c;
    }
}

std::vector<double> MyMatrix::wylicz(std::vector<double> B, int const MatrixSize)
{
    std::vector<double> wynik(MatrixSize);
    double suma;
    for (int i = MatrixSize - 1; i >= 0; i--)
    {
        suma = 0;
        for (int j = MatrixSize - 1; j > i; j--)
        {
            suma += tab[i][j] * wynik[j];
        }
        wynik[i] = (B[i] - suma) / tab[i][i];
        if (wynik[i] == -0.0) {
            wynik[i] = 0.0;
        }
    }
    return wynik;
}

std::vector<double> MyMatrix::czesc_wybor(std::vector<double> B, int const MatrixSize)
{
    double pom, maximum, modul, wspolcz;
    std::vector<double> wynik(MatrixSize);
    int wiersz;
    for (int n = 0; n < MatrixSize - 1; n++)
    {
        // sortujemy wiersze wed³ug kolejnoœci w n-tej kolumnie
        wiersz = n;
        maximum = std::abs(tab[n][n]);
        for (int i = n + 1; i < MatrixSize; i++)
        {
            modul = std::abs(tab[i][n]);
            if (maximum < modul && modul !=0)
            {
                maximum = modul;
                wiersz = i;
            }
        }

        if (wiersz != n)
        {
            zamien_wiersz(n, wiersz, n);
            pom = B[n];
            B[n] = B[wiersz];
            B[wiersz] = pom;
        }

        // wykonujemy zerowanie wszystkich wierszy poni¿ej n-tego
        for (int i = n + 1; i < MatrixSize; i++)
        {
            wspolcz = tab[i][n] / tab[n][n];
            for (int j = n; j < MatrixSize; j++)
            {
                tab[i][j] = tab[i][j] - wspolcz * tab[n][j];
            }
            B[i] = B[i] - wspolcz * B[n];
        }
    }
    wynik = wylicz(B, MatrixSize);
    return wynik;
}

std::vector<double> MyMatrix::czesc_wybor_ulepszony(std::vector<double> B, int const MatrixSize){
        double pom, maximum, modul, wspolcz;
        std::vector<double> wynik(MatrixSize);
        int wiersz;
        for (int n = 0; n < MatrixSize - 1; n++)
        {
            // sortujemy wiersze według kolejności w n-tej kolumnie
            wiersz = n;
            maximum = std::abs(tab[n][n]);
            for (int i = n + 1; i < MatrixSize; i++)
            {
                modul = std::abs(tab[i][n]);
                if (maximum < modul && modul !=0)
                {
                    maximum = modul;
                    wiersz = i;
                }
            }

            if (wiersz != n)
            {
                zamien_wiersz(n, wiersz, n);
                pom = B[n];
                B[n] = B[wiersz];
                B[wiersz] = pom;
            }

            // wykonujemy zerowanie wszystkich wierszy poniżej n-tego
            for (int i = n + 1; i < MatrixSize; i++)
            {
                if(tab[i][n] == 0)
                {
                    continue;
                }
                wspolcz = tab[i][n] / tab[n][n];
                for (int j = n; j < MatrixSize; j++)
                {
                    tab[i][j] = tab[i][j] - wspolcz * tab[n][j];
                }
                B[i] = B[i] - wspolcz * B[n];
            }
        }
        wynik = wylicz(B, MatrixSize);
        return wynik;
}

std::vector<double> MyMatrix::jacobi(std::vector<double> B, int const MatrixSize) {
		std::vector<double> X(MatrixSize);
		std::vector<double> X1(MatrixSize);
		double pom,result;
		double sum;
		int ile=0;

		for(int i=0; i<MatrixSize; i++)
        {
            X[i] = 0;
        }
        do
        {
            ile++;
            X1 = X;
            for(int i=1; i<MatrixSize; i++)
                    {
                        sum = 0;
                        for(int j=1;j<MatrixSize;j++)
                        {
                            if(i != j)
                                sum += tab[i][j]*X[j];
                        }
                        X[i] = (-sum + B[i]) / tab[i][i];
                        if(X[i] == -0.0)
                            X[i] = 0.0;
                    }

                    pom = 0;
                    result=0;
                    for(int i=0;i<MatrixSize;i++)
                    {
                        pom = fabs(X[i]-X1[i]);
                        result += pom*pom;
                    }
                    result = sqrt(result);
        }while(result > eps);

		std::cout<<ile<<std::endl;
		return X1;
}

std::vector<double> MyMatrix::jacobiIteracje(std::vector<double> B, int const MatrixSize, int iteracje) {
		std::vector<double> X(MatrixSize);
		std::vector<double> X1(MatrixSize);
		double pom,result;
		double sum;
		int ile=0;

		for(int i=0; i<MatrixSize; i++)
        {
            X[i] = 0;
        }
        do
        {
            ile++;
            X1 = X;
            for(int i=1; i<MatrixSize; i++)
                    {
                        sum = 0;
                        for(int j=1;j<MatrixSize;j++)
                        {
                            if(i != j)
                                sum += tab[i][j]*X[j];
                        }
                        X[i] = (-sum + B[i]) / tab[i][i];
                        if(X[i] == -0.0)
                            X[i] = 0.0;
                    }

                    pom = 0;
                    result=0;
                    for(int i=0;i<MatrixSize;i++)
                    {
                        pom = fabs(X[i]-X1[i]);
                        result += pom*pom;
                    }
                    result = sqrt(result);
        }while(ile < iteracje);
		return X1;
}

std::vector<double> MyMatrix::GSeidel(std::vector<double> B, int const MatrixSize) {
        std::vector<double> X(MatrixSize);
        std::vector<double> X1(MatrixSize);
        double pom,result;
        double sum;
        double sumdwa;
        int ile=0;

        for(int i=0; i<MatrixSize; i++)
        {
            X[i] = 0;
        }
        do
        {
            ile++;
            X1 = X;
            for(int i=1; i<MatrixSize; i++)
                    {
                        sum = 0;
                        sumdwa=0;
                        for(int j=1;j<i-1;j++){
                            sum+=tab[i][j]*X[j];
                        }
                        for(int j=i+1;j<MatrixSize;j++){
                            sumdwa+= tab[i][j] * X[j];
                        }
                        X[i] = ((-sum-sumdwa) + B[i]) / tab[i][i];
                        if(X[i] == -0.0)
                            X[i] = 0.0;
                    }

                    pom = 0;
                    result=0;
                    for(int i=0;i<MatrixSize;i++)
                    {
                        pom = fabs(X[i]-X1[i]);
                        result += pom*pom;
                    }
                    result = sqrt(result);
        }while(result > eps);

        std::cout<<ile<<std::endl;
        return X1;
}

std::vector<double> MyMatrix::GSeidelIteracje(std::vector<double> B, int const MatrixSize, int iteracje) {
        std::vector<double> X(MatrixSize);
        std::vector<double> X1(MatrixSize);
        double pom,result;
        double sum;
        double sumdwa;
        int ile=0;

        for(int i=0; i<MatrixSize; i++)
        {
            X[i] = 0;
        }
        do
        {
            ile++;
            X1 = X;
            for(int i=1; i<MatrixSize; i++)
                    {
                        sum = 0;
                        sumdwa=0;
                        for(int j=1;j<i-1;j++){
                            sum+=tab[i][j]*X[j];
                        }
                        for(int j=i+1;j<MatrixSize;j++){
                            sumdwa+= tab[i][j] * X[j];
                        }
                        X[i] = ((-sum-sumdwa) + B[i]) / tab[i][i];
                        if(X[i] == -0.0)
                            X[i] = 0.0;
                    }

                    pom = 0;
                    result=0;
                    for(int i=0;i<MatrixSize;i++)
                    {
                        pom = fabs(X[i]-X1[i]);
                        result += pom*pom;
                    }
                    result = sqrt(result);
        }while(ile < iteracje);

        return X1;
}
