#ifndef PROBABILITY_H_INCLUDED
#define PROBABILITY_H_INCLUDED

class Probability {
public:
	int yCount;
	int nCount;
	int uCount;
	int allCount;

	double doNothing = 0;
	double plusN = 0;
	double plusY = 0;
	double minusYandN = 0;
	long omega;

	Probability(int, int, int);

	void countProbability();
private:
	long countPairs(long);
	static long getSumFromPartOfFactorial(long, long);
};

#endif // PROBABILITY_H_INCLUDED
