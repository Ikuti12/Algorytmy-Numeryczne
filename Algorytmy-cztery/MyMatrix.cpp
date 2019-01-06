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
		std::vector<double> x1(MatrixSize);
		std::vector<double> x2(MatrixSize);
		bool koniec;
		std::vector<std::vector<double> > M(MatrixSize);
		for (int i = 0; i < MatrixSize; i++)
			M[i].resize(MatrixSize);
		for (int i = 0; i < MatrixSize; i++) {
			X[i] = (double) (1 / tab[i][i]);

		}
		for (int i = 0; i < MatrixSize; i++) {
			for (int j = 0; j < MatrixSize; j++) {
				if (i == j) {
					M[i][j] = 0;
				}
				else {
					M[i][j] = -(tab[i][j] * X[i]);
				}
			}
		}
		for (int i = 0; i < MatrixSize; i++) {
			x1[i] = 0;
		}
		for (int i = 0; i < MatrixSize; i++) {
                x2[i] = (double)(X[i] * B[i]);
                for (int j = 0; j < MatrixSize; j++) {
                    x2[i] += (double)(M[i][j] * x1[j]);
                }
            }
        koniec=true;

		do {
		    for (int i = 0; i < MatrixSize; i++) {
				x1[i] = x2[i];
			}
			for (int i = 0; i < MatrixSize; i++) {
				x2[i] = (double)(X[i] * B[i]);
				for (int j = 0; j < MatrixSize; j++) {
					x2[i] += (double)(M[i][j] * x1[j]);
				}
			}
            for (int i = 0; i < MatrixSize; i++)
                {
                    if (fabs(x2[i] - x1[i]) > eps) { koniec = true; break; }
                    else koniec = false;
                }

		}while(koniec);
		//std::cout << "Wynik Jacobi" << std::endl;
		//for (int i = 0; i < MatrixSize; i++) {
		//	std::cout << x1[i] << std::endl;
		//}
		return x1;
	}

std::vector<double> MyMatrix::jacobiPoprawione(std::vector<double> B, int const MatrixSize) {
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

std::vector<double> MyMatrix::GSeidelPoprawione(std::vector<double> B, int const MatrixSize) {
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

std::vector<double> MyMatrix::GSeidel(std::vector<double> B, int const MatrixSize) {
		double dzielnik;
		double norma1, norma2;
		bool koniec;
		std::vector<double> B2(MatrixSize);
		std::vector<double> x(MatrixSize);
		std::vector<double> x2(MatrixSize);
		std::vector<std::vector<double> > U(MatrixSize);
		for (int i = 0; i < MatrixSize; i++)
			U[i].resize(MatrixSize);
		std::vector<std::vector<double> > L(MatrixSize);
		for (int i = 0; i < MatrixSize; i++)
			L[i].resize(MatrixSize);
		std::vector<std::vector<double> > D(MatrixSize);
		for (int i = 0; i < MatrixSize; i++)
			D[i].resize(MatrixSize);
		std::vector<std::vector<double> > I(MatrixSize);
		for (int i = 0; i < MatrixSize; i++)
			I[i].resize(MatrixSize);
		for (int i = 0; i < MatrixSize; i++) {
			for (int j = 0; j < MatrixSize; j++)
			{
				L[i][j] = 0;
				U[i][j] = 0;
				D[i][j] = 0;
				I[i][j] = 0;
			}
		}
		for (int i = 0; i < MatrixSize; i++)
		{
			for (int j = 0; j < MatrixSize; j++)
			{
				// utworzenie macierzy dolnotrojkatnej
				if (i > j)
					L[i][j] = tab[i][j];
				// utworzenie macierzy gornotrojkatnej
				else if (i < j)
					U[i][j] = tab[i][j];
				else
				{
					// utworzenie macierzy z przekatnej macierzy glownej
					D[i][i] = tab[i][i];
					// jedynki na przekatnej
					I[i][i] = 1;
				}
			}
		}
		// wyznaczenie p = min{||-(L+D)^-1 * N||:||-(L+D)^-1 * N||}

        // obliczenie D + L (wynik w L)
		for (int i = 0; i < MatrixSize; i++)
			L[i][i] = D[i][i];
		// obliczenie (D + L)^-1 (wynik w I)
		for (int i = 0; i < MatrixSize; i++)
		{
			dzielnik = L[i][i];
			for (int j = 0; j < MatrixSize; j++)
			{
				L[i][j] /= (double)dzielnik;
				I[i][j] /= (double)dzielnik;
			}
			for (int k = 0; k < MatrixSize; k++)
			{
				if (k == i) continue;
				dzielnik = L[k][i];
				for (int j = 0; j < MatrixSize; j++)
				{
					L[k][j] -= (double)(L[i][j] * dzielnik);
					I[k][j] -= (double)(I[i][j] * dzielnik);
				}
			}
		}
		// obliczenie (L + D)^-1 * N (wynik w L)
		for (int i = 0; i < MatrixSize; i++)
		{
			for (int j = 0; j < MatrixSize; j++)
			{
				L[i][j] = 0;
				for (int k = 0; k < MatrixSize; k++)
					L[i][j] += (double)(I[i][k] * U[k][j]);
			}
		}

		// obliczenie -(L+D)^-1 * N
		for (int i = 0; i < MatrixSize; i++)
			for (int j = 0; j < MatrixSize; j++)
				L[i][j] *= (double)(-1);
		// norma "jeden"
		norma1 = 0;
		for (int i = 0; i < MatrixSize; i++)
		{
			dzielnik = 0;
			for (int j = 0; j < MatrixSize; j++)
				dzielnik += fabs(L[j][i]);
			if (dzielnik > norma1) norma1 = dzielnik;
		}
		// norma "nieskonczonosc"
		norma2 = 0;
		for (int i = 0; i < MatrixSize; i++)
		{
			dzielnik = 0;
			for (int j = 0; j < MatrixSize; j++)
				dzielnik += fabs(L[i][j]);
			if (dzielnik > norma2) norma2 = dzielnik;
		}

		// p = min(norma1, norma2)
		if (norma1 > norma2) norma1 = norma2;

		// ciag nie jest zbiezny do rozwiazania ukladu rownan
		if (norma1 >= 1)
			norma1 = 0.5;

		// inicjalizacja wektora wynikow
		for (int i = 0; i < MatrixSize; i++)
			x[i] = 0;

		koniec = true;
		do
		{
			// przepisanie x - aktualnych wynikow do x2 - wynikow z poprzedniej iteracji
			for (int i = 0; i < MatrixSize; i++)
				x2[i] = x[i];

			// wykonanie kolejnej iteracji
			for (int i = 0; i < MatrixSize; i++)
			{
				B2[i] = B[i];

				for (int j = 0; j < i; j++)
					B2[i] -= (tab[i][j] * x[j]);

				for (int j = i + 1; j < MatrixSize; j++)
					B2[i] -= tab[i][j] * x2[j];

				x[i] = (double)(B2[i] / tab[i][i]);
			}

			// sprawdzenie warunku zakonczenia: ||x(k)-x(k-1)|| <= epsilon
			for (int i = 0; i < MatrixSize; i++)
			{
				if (fabs(x[i] - x2[i]) > eps) { koniec = true; break; }
				else koniec = false;
			}
		} while (koniec);

		// wyswietlenie wyniku
		//std::cout << "Wynik GSeidel:\n";
		for (int i = 0; i < MatrixSize; i++) {
			if (x[i] == -0.0) {
				x[i] = 0.0;
			}
		}
		return x;

}
