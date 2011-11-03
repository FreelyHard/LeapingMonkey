#include "Legendre.h"
#include "LegendreAbscissas.h"
#include "lapacke.h"
#include "cblas.h"
#include <stdio.h>
#include <iostream>

double Legendre::alpha(int n, double x) const {
  return -(2.0*n - 1.0)/(double)n*x;
}

double Legendre::beta(int n, double x) const {
  return (n - 1.0)/(double)n;
}

double Legendre::function(int N, double x) const {
  if (N == 0) {
    return 1.0;
  } else if (N == 1) {
    return x;
  } else {
    return  ((2.0*N - 1.0)*x*function(N-1,x) 
        - (N-1.0)*function(N-2,x))/(double)N;
  }
}

bool Legendre::setRank(int N) {
  delete[] differentiationMatrix;
  differentiationMatrix = NULL;
  delete[] valuesToCoefficients;
  valuesToCoefficients = NULL;
  abscissas = LegendreAbscissas::getAbscissas(N);
  if (abscissas == 0) {
    nBasis = 0;
    return false;
  }
  nBasis = N;
  return true;
}

double* Legendre::coefficientsOfDerivativeMatrix() const {
  double* differentiate = new double[nBasis*nBasis];
  for (int i = 0; i < nBasis*nBasis; i++) differentiate[i] = 0.0;
  for (int iRow = 0; iRow < nBasis; iRow++) {
    for (int iCol = iRow + 1; iCol < nBasis; iCol += 2) {
      int i = index(iRow, iCol);
      differentiate[i] = 2.0*iRow + 1.0;
    }
  }
  return differentiate;
}

double Legendre::evaluate(double x, const double* coeffs) const {
  double pMinus1 = 1.0;
  double p = x;
  double result = coeffs[0] + x*coeffs[1];
  for (int n = 1; n < nBasis-1; n++) {
    double pPlus1 = (double)(2*n+1)/(double)(n+1)*x*p
      - (double)n/(double)(n+1)*pMinus1;
    result += coeffs[n+1]*pPlus1;
    pMinus1 = p;
    p = pPlus1;
  }
  return result;
}
