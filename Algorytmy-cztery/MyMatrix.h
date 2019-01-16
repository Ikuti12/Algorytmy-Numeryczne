#ifndef MYMATRIX_H_INCLUDED
#define MYMATRIX_H_INCLUDED

#define NOTN 50 // Total amount of agents here!

class MyMatrix
{
    public:
	// konstruktor i dekonstruktor
	MyMatrix(double**, int const);
	~MyMatrix();

	void WyswietlMacierz(int const);

	// przeciazenia operatorów
	std::vector<double>& operator[](int);
	//MyMatrix& operator+(MyMatrix a);

	void zamien_wiersz(int, int, int);
	std::vector<double> wylicz(std::vector<double>, int const);
	std::vector<double> czesc_wybor(std::vector<double>, int const);
	std::vector<double> czesc_wybor_ulepszony(std::vector<double>, int const);
	std::vector<double> jacobi(std::vector<double>, int const);
	std::vector<double> jacobiIteracje(std::vector<double>, int const, int);
	std::vector<double> GSeidel(std::vector<double>, int const);
	std::vector<double> GSeidelIteracje(std::vector<double>, int const, int);
	std::vector<std::vector<double>> tab;
};



#endif // MYMATRIX_H_INCLUDED
