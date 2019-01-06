#include "Probability.h"

Probability::Probability(int yCount, int nCount, int allCount)
{
    this->yCount = yCount;
    this->nCount = nCount;
    this->allCount = allCount;
    this->uCount = allCount - yCount - nCount;
    if (yCount + nCount > allCount) {
        throw " za duze!!";
    }
}

void Probability::countProbability()
{
    omega = countPairs(allCount);

    // do nothing
    long YY = countPairs(yCount);
    long NN = countPairs(nCount);
    long UU = countPairs(uCount);

    // do nothing
    doNothing = (YY + NN + UU) / (double)omega;

    // plus N
    long NU = countPairs(nCount + uCount) - UU - NN;
    plusN = NU / (double)omega;

    // plus Y
    long YU = countPairs(yCount + uCount) - UU - YY;
    plusY = YU / (double)omega;

    // minus N and Y
    long YN = countPairs(yCount + nCount) - NN - YY;
    minusYandN = YN / (double)omega;

    //Y = (YU)/(double) omega;
    //N = (NU)/(double) omega;
    //U = (YY + NN + UU + YN)/(double) omega;

}

long Probability::countPairs(long n)
{
    if (n > 1) {
        return getSumFromPartOfFactorial(n - 1, n) / 2;
    }
    else return 0;
}

long Probability::getSumFromPartOfFactorial(long from, long to)
{
    long result = 1;
    for (long i = from; i <= to; i++) {
        result *= i;
    }
    return result;
}
