#include "Chebyshev.h"
#include "lapacke.h"
#include "cblas.h"
#include <stdio.h>
#include <iostream>
#include <math.h>

Chebyshev::~Chebyshev() {
  delete[] abscissas;
}

double Chebyshev::alpha(int n, double x) {
  return -2.0*x;
}

double Chebyshev::beta(int n, double x) {
  return 1.;
}

double Chebyshev::function(int N, double x) const {
  // TODO: remove recursion for speedup.
  if (N == 0) {
    return 1.0;
  } else if (N == 1) {
    return x;
  } else {
    return  2.0*x*function(N-1,x) - function(N-2, x);
  }
}

double* Chebyshev::coefficientsOfDerivativeMatrix() {
  // From formula A.15 of Boyd.
  double* differentiate = new double[nBasis*nBasis];
  for (int i = 0; i < nBasis*nBasis; i++) differentiate[i] = 0.0;
  differentiate[index(nBasis-2, nBasis-1)] = 2.0*(nBasis-1);
  for (int iRow = nBasis-2; iRow >= 0; iRow--) {
    // The part that depends on the next order coefficient of f
    differentiate[index(iRow, iRow + 1)] = 2.0*(iRow + 1);
    // The part that depends on the next order coefficient of f'
    for (int iCol = iRow + 3; iCol < nBasis; iCol++) {
      differentiate[index(iRow, iCol)] += differentiate[index(iRow+2,iCol)];
    }
    // Take care of c_k
    for (int iCol = 0; iCol < nBasis; iCol++)
      differentiate[index(0,iCol)] *= 0.5;
  }
  return differentiate;
}

double Chebyshev::evaluate(double x, double* coeffs) {
  double tMinus1 = 1.0;
  double t = x;
  double result = coeffs[0] + x*coeffs[1];
  for (int n = 1; n < nBasis-1; n++) {
    double tPlus1 = 2.0*x*t - tMinus1;
    result += coeffs[n+1]*tPlus1;
    tMinus1 = t;
    t = tPlus1;
  }
  return result;
}

double* Chebyshev::chebyshevToTaylorCoefficientsMatrix() {
  double* matrix = new double[nBasis*nBasis];
  for (int i = 0; i < nBasis*nBasis; i++) matrix[i] = 0;
  matrix[0] = 1.; // T_0(x) = 1 = x^0.
  matrix[index(1,1)] = 1.; // T_1(x) = x = x^1
  for (int iCheb = 2; iCheb < nBasis; iCheb++) {
    matrix[index(0, iCheb)] = -matrix[index(0,iCheb-2)];
    for (int iTaylor = 1; iTaylor < nBasis; iTaylor++) {
      int iMinus = index(iTaylor-1, iCheb-1); // x*T_{n-1}(x)
      int iMinus2 = index(iTaylor, iCheb-2);
      int i = index(iTaylor, iCheb);
      matrix[i] = 2.*matrix[iMinus] - matrix[iMinus2];
    }
  }
  return matrix;
}
