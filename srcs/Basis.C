#include "Basis.h"
#include "mkl_lapack.h"
#include "mkl_cblas.h"
#include <stdio.h>
#include <iostream>

Basis::~Basis() {
  delete[] differentiationMatrix;
  delete[] valuesToCoefficients;
}

int Basis::getRank() {
  return nBasis;
}

const double* Basis::getAbscissas() {
  if (nBasis != 0) {
    return abscissas;
  }
  return 0;
}

double* Basis::coefficientsToValuesMatrix() {
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

const double* Basis::getValuesToCoefficientsMatrix() {
  if (valuesToCoefficients == NULL) {
    valuesToCoefficients = coefficientsToValuesMatrix();
    int pivots[nBasis];
    int error;
    // Lapack LU decomposition.
    dgetrf(&nBasis,&nBasis,valuesToCoefficients,&nBasis,pivots,&error);
    // Compute the inverse from the LU decomposition.
    int lwork = nBasis*nBasis;
    double scratchSpace[lwork]; // Should compute this optimally...
    dgetri(&nBasis, valuesToCoefficients, &nBasis, pivots,
        scratchSpace, &lwork, &error);
  }
  return valuesToCoefficients;
}

int Basis::index(int iRow, int iCol) {
  return nBasis*iRow + iCol;
}

const double* Basis::getDifferentiationMatrix() {
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

void Basis::fillCoefficients(const double* values, double* coeffs) {
  const double* v2c = getValuesToCoefficientsMatrix();
  for (int iRow = 0; iRow < nBasis; iRow++) {
    coeffs[iRow] = 0.0;
    for (int iCol = 0; iCol < nBasis; iCol++) {
      coeffs[iRow] += v2c[index(iRow,iCol)]*values[iCol];
    }
  }
}

double* Basis::interpolate(const double* values, 
    const double* x, int nPoints) {
  double* fAtX = new double[nPoints];
  double coeffs[nBasis];
  fillCoefficients(values, coeffs);
  for (int i = 0; i < nPoints; i++) {
    fAtX[i] = evaluate(x[i], coeffs);
  }
  return fAtX;
}

double Basis::integrate(const double* values) {
  double integral = 0.0;
  const double* weights = getAbscissas();
  for (int i = 0; i < nBasis; i++) integral += weights[i + nBasis]*values[i];
  return integral;
}
