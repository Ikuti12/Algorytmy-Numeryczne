#include <iostream>
#include <algorithm>
#include <vector>
#include "Generator.h"
#include "Probability.h"

double** Generator::generujMacierz(int allCount)
{
    int matrixSize = (allCount + 1)*(allCount + 2) / 2;
    double** wynik = new double*[matrixSize];

    for (int i = 0; i < matrixSize; i++)
        wynik[i] = new double[matrixSize];

    int nCount = 0;
    int yCount = 0;
    int uCount = 0;
    int Y = 0;
    int N = 0;
    int x = 0;
    int y = 0;

    for (int temp = allCount + 1; temp >= 0; temp--)
    {
        for (int i = temp; i > 0; i--)
        {
            nCount = N;
            uCount = i - 1;
            yCount = allCount - nCount - uCount;
            N++;
            std::vector<int> v = { yCount,nCount,uCount };
            lista.push_back(v);
        }
        N = 0;
        Y++;
    }

    //for (int i = 0; i < lista.size(); i++)
     //   std::cout << lista[i][0] << lista[i][1] << lista[i][2] << std::endl;

    //std::cout <<"To jest miejsce" << znajdzMiejsce(1,2,0) << std::endl;
    int yWithOne = 0;

    for (int i = 0; i < matrixSize; i++)
        for (int j = 0; j < matrixSize; j++)
            wynik[i][j] = 0;

    for (int temp = allCount + 1; temp >= 0; temp--)
    {
        for (int i = temp; i > 0; i--)
        {
            nCount = N;
            uCount = i - 1;
            yCount = allCount - nCount - uCount;
            Probability probability(yCount, nCount, allCount);
            probability.countProbability();
            //std::cout <<"[Y: " << Y << ", N: " << N << "] , U: " << uCount << ", N: " << nCount << ", Y: " << yCount << " | " <<"Probabilities: plusN: " << probability.plusN << ", plusY: " << probability.plusY << ", minusYandN: " << probability.minusYandN << ", doNothing: " << probability.doNothing << std::endl;
            N++;
            if (probability.plusN != 0) {
                int indeks = znajdzMiejsce(yCount, nCount+1, uCount-1);
                wynik[x][indeks] = probability.plusN;
            }
            if (probability.plusY != 0) {
                int indeks = znajdzMiejsce(yCount+1, nCount, uCount - 1);
                wynik[x][indeks] = probability.plusY;
            }
            if (probability.minusYandN != 0) {
                int indeks = znajdzMiejsce(yCount-1, nCount - 1, uCount +2);
                wynik[x][indeks] = probability.minusYandN;
            }
            if ((x==0 && y==0)|| (x == matrixSize - 1 && y == matrixSize - 1) || (x==allCount && y==allCount)) {
                wynik[x][y] = probability.doNothing - 2;
            }
            else {
            wynik[x][y] = probability.doNothing - 1;
            }
            x++;
            y++;
        }
        N = 0;
        Y++;
    }
    return wynik;
}

int Generator::znajdzMiejsce(int Y, int N, int U) {
		std::vector<int> szukanie = { Y,N,U };
		auto miejsce = find(lista.begin(), lista.end(), szukanie);
		if (miejsce == lista.end()) {
			std::cout << "Nie znalazlem!";
		}
		return miejsce - lista.begin();
	}
