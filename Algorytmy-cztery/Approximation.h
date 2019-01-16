#ifndef APPROXIMATION_H_INCLUDED
#define APPROXIMATION_H_INCLUDED

class Approximation
{
    public:
	// konstruktor i dekonstruktor
	Approximation();
	~Approximation();

	std::vector<double> run(std::vector<double>, std::vector<double>, int);
};

#endif // APPROXIMATION_H_INCLUDED
