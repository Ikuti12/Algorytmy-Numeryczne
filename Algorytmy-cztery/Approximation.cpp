#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "Approximation.h"
#include "MyMatrix.h"

Approximation::Approximation() {}
Approximation::~Approximation() {}

std::vector<double> Approximation::run(std::vector<double> x, std::vector<double> y, int m)
{
    int n = x.size();
    if(m < n)
    {
        std::vector<double> b(m);
        MyMatrix F(n, m+1);

        for(int i=0;i<m;i++)
            for(int j=0;j<n;j++)
                F[j][i] = pow(x[j],i);

        MyMatrix trans_F = F.transposed(m,n);
        MyMatrix A = trans_F * F;
        b = trans_F * y;
        std::vector<double> c = A.inversion() * b;
        return c;
    }
    else
    {
        throw std::invalid_argument( "Bad value for m!\n" );
    }
}
