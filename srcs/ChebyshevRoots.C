#include "ChebyshevRoots.h"
#include "lapacke.h"
#include "cblas.h"
#include <stdio.h>
#include <iostream>
#include <math.h>

bool ChebyshevRoots::setRank(int N) {
  if (N > 1) {
    delete[] differentiationMatrix;
    differentiationMatrix = NULL;
    delete[] valuesToCoefficients;
    valuesToCoefficients = NULL;
    nBasis = N;
    double* myAbscissas = new double[2*nBasis];

    // First the actual abscissas
    for (int i = 0; i < nBasis; i++) {
      myAbscissas[nBasis-1-i] = cos(M_PI*(double)(2*i+1)/(double)(2*nBasis));
    }
    // Now the Gauss-Chebyshev quadrature weights.
    double invN = 1.0/(double)nBasis;
    for (int i = 0; i < nBasis; i++) {
      myAbscissas[nBasis + i] = M_PI*invN;
    }

    // Set the pointer properly.
    delete[] abscissas;
    abscissas = myAbscissas;
    return true;
  }
  return false;
}

double *ChebyshevRoots::getAbscissas(int N) const {
  double *x = new double[N];
    for (int i = 0; i < N; i++) {
      x[N - 1 - i] = cos(M_PI*(double)(2*i+1)/(double)(2*N));
    }
  return x;
}
