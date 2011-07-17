#include <cstdlib>
#include <iostream>
#include "Basis2D.h"

void Basis2D::setRanks(int n1, int n2) {
  basis1->setRank(n1);
  basis2->setRank(n2);
  delete[] differentiationMatrices;
  differentiationMatrices = NULL;
}

void Basis2D::fillCoefficients(const double* u, double* coefficients) {
  int N1 = basis1->getRank();
  int N2 = basis2->getRank();
  // The easy coordinates...
  for (int iCoord2 = 0; iCoord2 < N2; iCoord2++) {
    int index = functionIndex(0, iCoord2);
    // Pointer arithmetic.
    basis1->fillCoefficients(u + index, coefficients + index);
  }
  // The complicated coordinates...
  double coeffs1d[N2];
  double coeffs2d[N2];
  for (int iCoord1 = 0; iCoord1 < N1; iCoord1++) {
    // Copy 1d coefficients
    for (int iCoord2 = 0; iCoord2 < N2; iCoord2++) {
      int index = functionIndex(iCoord1, iCoord2);
      coeffs1d[iCoord2] = coefficients[index];
    }
    basis2->fillCoefficients(coeffs1d, coeffs2d);
    // Copy the 2d coefficients into the coefficients array.
    for (int iCoord2 = 0; iCoord2 < N2; iCoord2++) {
      int index = functionIndex(iCoord1, iCoord2);
      coefficients[index] = coeffs2d[iCoord2];
    }
  }
}

double* Basis2D::tensorInterpolate(const double* u, const double* coord1,
    int n1, const double* coord2, int n2) {
  int nCoord1 = basis1->getRank();
  int nCoord2 = basis2->getRank();
  double* uInterpolatedAlongCoord2[nCoord1];
  for (int iCoord1 = 0; iCoord1 < nCoord1; iCoord1++) {
    double uConstCoord1[nCoord2];
    for (int iCoord2 = 0; iCoord2 < nCoord2; iCoord2++) {
      int index = functionIndex(iCoord1, iCoord2);
      uConstCoord1[iCoord2] = u[index];
    }
    uInterpolatedAlongCoord2[iCoord1] = 
      basis2->interpolate(uConstCoord1, coord2, n2);
  }
  double* uInterpolated = new double[n1*n2];
  for (int i2 = 0; i2 < n2; i2++) {
    double uConstCoord2[nCoord1];
    for (int iCoord1 = 0; iCoord1 < nCoord1; iCoord1++) {
      uConstCoord2[iCoord1] = uInterpolatedAlongCoord2[iCoord1][i2];
    }
    double* uInterpolatedAlongCoord1 = 
      basis1->interpolate(uConstCoord2, coord1, n1);
    for (int i1 = 0; i1 < n1; i1++) {
      int index = i2*n1 + i1; // Copied from functionIndex...
      uInterpolated[index] = uInterpolatedAlongCoord1[i1];
    }
    //delete[] uInterpolatedAlongCoord1;
  }
  for (int i = 0; i < nCoord1; i++) delete[] uInterpolatedAlongCoord2[i];
  //delete[] uInterpolatedAlongCoord2;
  return uInterpolated;
}

const double* Basis2D::getDifferentiationMatrices() {
  if (differentiationMatrices == NULL) {
    int n1 = basis1->getRank();
    int n2 = basis2->getRank();
    int n2D = n1*n2;
    //Instantiate and zero the matrices.
    delete[] differentiationMatrices;
    differentiationMatrices = new double[2*n2D*n2D];
    for (int i = 0; i < 2*n2D*n2D; i++) differentiationMatrices[i] = 0.;

    const double* diff1 = basis1->getDifferentiationMatrix();
    const double* diff2 = basis2->getDifferentiationMatrix();
    for (int iCoord1 = 0; iCoord1 < n1; iCoord1++) {
      for (int iCoord2 = 0; iCoord2 < n2; iCoord2++) {
        int iRow = functionIndex(iCoord1,iCoord2);
        // The derivative along COORD1.
        // In this case the derivative only depends on function values that
        // share the same value of coord2.
        for (int jCoord1 = 0; jCoord1 < n1; jCoord1++) {
          int iCol = functionIndex(jCoord1, iCoord2);
          int iMatrix = matrixIndex(iRow, iCol);
          int i1D = basis1->index(iCoord1, jCoord1);
          differentiationMatrices[iMatrix] = diff1[i1D];
        }
        // The derivative along COORD2.
        // In this case the derivative only depends on function values that
        // share the same value of coord1.
        for (int jCoord2 = 0; jCoord2 < n2; jCoord2++) {
          int iCol = functionIndex(iCoord1, jCoord2);
          int iMatrix = matrixIndex(iRow, iCol);
          int i1D = basis2->index(iCoord2, jCoord2);
          differentiationMatrices[iMatrix + n2D*n2D] = diff2[i1D];
        }
      }
    }
  }
  return differentiationMatrices;
}

Basis2D::~Basis2D() {
  delete[] differentiationMatrices;
}

int Basis2D::functionIndex(int iCoord1, int iCoord2) {
  return iCoord2*basis1->getRank() + iCoord1;
}

int Basis2D::matrixIndex(int outputFunctionIndex, int inputFunctionIndex) {
  return outputFunctionIndex*(basis1->getRank()*basis2->getRank()) 
    + inputFunctionIndex;
}
