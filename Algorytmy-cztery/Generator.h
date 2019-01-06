#ifndef GENERATOR_H_INCLUDED
#define GENERATOR_H_INCLUDED

class Generator {
public:
	double** generujMacierz(int);
private:
	int znajdzMiejsce(int, int, int);
	std::vector<std::vector<int> > lista;
};

/* X =3
i Y N U
0 0 0 3
1 0 1 2
2 0 2 1
3 0 3 0
4 1 0 2
5 1 1 1
6 1 2 0
7 2 0 1
8 2 1 0
9 3 0 0

*/

#endif // GENERATOR_H_INCLUDED
