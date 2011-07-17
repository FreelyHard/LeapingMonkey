#include "ChebyshevExtrema.h"
#include "mkl_lapack.h"
#include "mkl_cblas.h"
#include <stdio.h>
#include <iostream>
#include <math.h>

bool ChebyshevExtrema::setRank(int N) {
  if (N > 1) {
    delete[] differentiationMatrix;
    differentiationMatrix = NULL;
    delete[] valuesToCoefficients;
    valuesToCoefficients = NULL;
    nBasis = N;
    double* myAbscissas = new double[2*nBasis];

    // First the actual abscissas
    for (int i = 0; i < nBasis; i++) {
      myAbscissas[nBasis - 1 - i] = cos(M_PI*(double)i/(double)(nBasis - 1));
    }
    // Now the Chebyshev-Lobatto quadrature weights.
    double pi_N1 = M_PI/(double)(nBasis - 1);
    myAbscissas[nBasis] = 0.5*pi_N1;
    for (int i = 1; i < nBasis-1; i++) {
      myAbscissas[nBasis + i] = pi_N1;
    }
    myAbscissas[2*nBasis - 1] = 0.5*pi_N1;
    delete[] abscissas;
    abscissas = myAbscissas;
    return true;
  }
  return false;
}
