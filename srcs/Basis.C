#include "Basis.h"
#include "lapacke.h"
#include "cblas.h"
#include <stdio.h>
#include <iostream>

Basis::~Basis() {
  delete[] differentiationMatrix;
  delete[] valuesToCoefficients;
}

int Basis::getRank() const {
  return nBasis;
}

const double* Basis::getAbscissas() const {
  if (nBasis != 0) {
    return abscissas;
  }
  return 0;
}

double* Basis::coefficientsToValuesMatrix() const {
  double *coefficientsToValuesMatrix = new double[nBasis*nBasis];
  /*
  for (int iRow = 0; iRow < nBasis; iRow++) {
    for (int iCol = 0; iCol < nBasis; iCol++) {
      int i = index(iRow,iCol);
      coefficientsToValuesMatrix[i] = function(iCol,abscissas[iRow]);
    }
  }*/
  // First two columns need to be done manually.
  for (int iRow = 0; iRow < nBasis; iRow++) {
    int i = index(iRow, 0);
    coefficientsToValuesMatrix[i] = function(0,abscissas[iRow]);
    i = index(iRow, 1);
    coefficientsToValuesMatrix[i] = function(1,abscissas[iRow]);
  }
  // Recurse on everything else...
  for (int iRow = 0; iRow < nBasis; iRow++) {
    for (int iCol = 2; iCol < nBasis; iCol++) {
      double a = alpha(iCol, abscissas[iRow]);
      double b = beta(iCol, abscissas[iRow]);
      int i = index(iRow, iCol);
      int iMinus1 = index(iRow, iCol - 1);
      int iMinus2 = index(iRow, iCol - 2);
      coefficientsToValuesMatrix[i] = -a*coefficientsToValuesMatrix[iMinus1]
        -b*coefficientsToValuesMatrix[iMinus2];
    }
  }
  return coefficientsToValuesMatrix;
}

const double* Basis::getValuesToCoefficientsMatrix() const {
#warning "Should maybe do some error checking."
  if (valuesToCoefficients == NULL) {
    valuesToCoefficients = coefficientsToValuesMatrix();
    int pivots[nBasis];
    // Lapack LU decomposition.
    int error = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, nBasis, nBasis,
        valuesToCoefficients, nBasis, pivots);
    // Compute the inverse from the LU decomposition.
    error = LAPACKE_dgetri(LAPACK_ROW_MAJOR, nBasis, valuesToCoefficients,
        nBasis, pivots);
  }
  return valuesToCoefficients;
}

int Basis::index(int iRow, int iCol) const {
  return nBasis*iRow + iCol;
}

const double* Basis::getDifferentiationMatrix() const {
  if (differentiationMatrix == NULL) {
    // First we transform values to coefficients, then we compute the
    // derivative using the coefficients.
    const double* constv2c = getValuesToCoefficientsMatrix();
    double v2c[nBasis*nBasis];
    for (int i = 0; i < nBasis*nBasis; i++) v2c[i] = constv2c[i];
    double* c2dc = coefficientsOfDerivativeMatrix();
    double v2dc[nBasis*nBasis];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        nBasis, nBasis, nBasis, 1.0, c2dc, nBasis, v2c, nBasis,
        0.0, v2dc, nBasis);
    delete[] c2dc;

    // Now the matrix will return the coefficients of the derivative and
    // we need to get the values of the derivative.
    double* c2v = coefficientsToValuesMatrix();
    differentiationMatrix = new double[nBasis*nBasis];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        nBasis, nBasis, nBasis, 1.0, c2v, nBasis, v2dc, nBasis,
        0.0, differentiationMatrix, nBasis);
    delete[] c2v;
  }
  return differentiationMatrix;
}

void Basis::fillCoefficients(const double* values, double* coeffs) const {
  const double* v2c = getValuesToCoefficientsMatrix();
  for (int iRow = 0; iRow < nBasis; iRow++) {
    coeffs[iRow] = 0.0;
    for (int iCol = 0; iCol < nBasis; iCol++) {
      coeffs[iRow] += v2c[index(iRow,iCol)]*values[iCol];
    }
  }
}

double* Basis::interpolate(const double* values, 
    const double* x, int nPoints) const {
  double* fAtX = new double[nPoints];
  double coeffs[nBasis];
  fillCoefficients(values, coeffs);
  for (int i = 0; i < nPoints; i++) {
    fAtX[i] = evaluate(x[i], coeffs);
  }
  return fAtX;
}

double Basis::integrate(const double* values) const {
  double integral = 0.0;
  const double* weights = getAbscissas();
  for (int i = 0; i < nBasis; i++) integral += weights[i + nBasis]*values[i];
  return integral;
}
