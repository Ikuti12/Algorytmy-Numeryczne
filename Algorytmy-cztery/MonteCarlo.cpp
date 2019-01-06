#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <ctime>
#include "MonteCarlo.h"

double MonteCarlo::monteCarlo(int iteracje, int yCount, int nCount, int allCount) {

	srand(time(0));
	double allYesCount = 0;
	double wynik;
	int uCount = allCount - yCount - nCount;
	for (int j = 0; j < iteracje; j++) {
		int uCountTemp = uCount;
		int yCountTemp = yCount;
		int nCountTemp = nCount;

		std::vector<Agent> agents(allCount);
		for (int i = 0; i < yCountTemp; i++) {
			agents[i] = Agent("Y");
		}
		for (int i = yCountTemp; i < nCountTemp + yCountTemp; i++) {
			agents[i] = Agent("N");
		}
		for (int i = nCountTemp + yCountTemp; i < allCount; i++) {
			agents[i] = Agent("U");
		}
		std::random_shuffle(agents.begin(), agents.end());
		int tescik = 0;
		while (!(yCountTemp == allCount || nCountTemp == allCount || uCountTemp == allCount)) {
			int rand1 = rand() % allCount;
			int rand2 = rand() % allCount;
			while (rand1 == rand2) {
				rand2 = rand() % allCount;
			}
			//cout << "Working" << endl;
			/*if (tescik < 20) {
				tescik++;
				cout << tescik << " "<< yCountTemp << nCountTemp << uCountTemp << endl;
				cout << rand1 << rand2 << endl;
				cout << tescik << " " << agents[rand1].getState() << agents[rand2].getState() << endl;
				for (int i = 0; i < allCount; i++) {
					cout << agents[i].getState();
				}
				cout << endl;
			}
			if (tescik < 20) {
				cout << tescik << " " << agents[rand1].getState() << agents[rand2].getState() << endl;
			}*/
			int result = agents[rand1].compare((&agents[rand2]));

			if (result == 1) {
				yCountTemp++;
				uCountTemp--;

			}
			else if (result == 2) {
				yCountTemp--;
				nCountTemp--;
				uCountTemp += 2;

			}
			else if (result == 3) {
				nCountTemp++;
				uCountTemp--;
			}
		}
		//cout << j << endl;

		if (yCountTemp == allCount) {
			allYesCount++;
		}
	}
	//cout << allYesCount << endl;
	//cout << allCount << endl;
	wynik = (double)allYesCount/(double)iteracje;
	//std::cout << "wynik z monte " << wynik << std::endl;
	return wynik;
}
