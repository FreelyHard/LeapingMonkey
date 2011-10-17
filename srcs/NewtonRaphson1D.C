#include <iostream>
#include <cmath>
#include <string>

#include "NewtonRaphson1D.h"

#include "lapacke.h"
#include "cblas.h"

bool NewtonRaphson1D::computeSolution(double* u, void* constants,
    double tolerance, int maxIterations) {
  setPointerToSolution(u);
  return computeSolution(tolerance, maxIterations);
}

void NewtonRaphson1D::setPointerToSolution(double* u) {
  sigma = u;
}

const double* NewtonRaphson1D::getResidue() {
  return residue;
}

const double* NewtonRaphson1D::getDerivative() {
  return sigmap;
}

bool NewtonRaphson1D::computeSolution(double tolerance, int maxIterations) {
  // Set up the guts. By assumption, the guess has been set in sigma.
  delete[] sigmap;
  delete[] sigmapp;
  delete[] residue;
  int nBasis = basis->getRank();
  sigmap = new double[basis->getRank()];
  sigmapp = new double[basis->getRank()];
  residue = new double[basis->getRank()];
  computeDifferentiationMatrices();
  computeDerivatives();
  // Start grinding.
  double totalResidue = getTotalResidue();
  int nIterations = 0;
  for (int i = 0; i < basis->getRank(); i++) {
  }
  while (totalResidue > tolerance && nIterations < maxIterations) {
    totalResidue = singleNewtonRaphson();
    nIterations++;
  }
  if (totalResidue < tolerance && totalResidue > 0.0) return true;
  return false;
}

double NewtonRaphson1D::map(double xi) {
  return xi;
}

double NewtonRaphson1D::dInverseMap(double xi) {
  return 1.0;
}

double NewtonRaphson1D::getTotalResidue() {
  double totalResidue = 0.;
  const double* xi = basis->getAbscissas();
  for (int i = 0; i < basis->getRank(); i++) {
    residue[i] = -f(map(xi[i]), sigma[i], sigmap[i], sigmapp[i]);
    if (i == 0 && gType != NONE) {
      residue[i] = -g(sigma[i], sigmap[i]);
    } else if (i == basis->getRank()-1 && hType != NONE) {
      if (hType == INITIAL) {
        residue[i] = -h(sigma[0], sigmap[0]);
      } else {
        residue[i] = -h(sigma[i], sigmap[i]);
      }
    }
    totalResidue += pow(residue[i], 2);
  }
  return sqrt(totalResidue);
}

void NewtonRaphson1D::computeDerivatives() {
  int nBasis = basis->getRank();
  for (int iRow = 0; iRow < nBasis; iRow++) {
    sigmap[iRow] = 0;
    sigmapp[iRow] = 0;
    for (int iCol = nBasis-1; iCol >= 0; iCol--) {
      int i = iRow*nBasis + iCol;
      sigmap[iRow] += diff[i]*sigma[iCol];
      sigmapp[iRow] += doubleDiff[i]*sigma[iCol];
    }
  }
}

void NewtonRaphson1D::registerBoundaries(
    BoundaryType gtype, BoundaryType htype) {
  gType = gtype;
  hType = htype;
}

void NewtonRaphson1D::computeDifferentiationMatrices() {
  int nBasis = basis->getRank();
  // The differentiation matrix follows by chain rule.
  const double* dbydxi = basis->getDifferentiationMatrix();
  diff = new double[nBasis*nBasis];
  const double* xi = basis->getAbscissas();
  for (int iRow = 0; iRow < nBasis; iRow++) {
    double dxidx = dInverseMap(xi[iRow]);
    for (int iCol = 0; iCol < nBasis; iCol++) {
      int index = iRow*nBasis + iCol;
      diff[index] = dxidx*dbydxi[index];
    }
  }
  // Double differentiation is the differentiation twice!
  doubleDiff = new double[nBasis*nBasis];
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
      nBasis, nBasis, nBasis, 1.0, diff, nBasis, diff, nBasis,
      0.0, doubleDiff, nBasis);
}

double NewtonRaphson1D::singleNewtonRaphson() {
  int nBasis = basis->getRank();
  double* jac = getJacobian();

  // The following is exactly lapack's dgesv written for the transpose.
  //  Because lapacke wants one to invert ldb, I had to figure this out...
  int pivots[nBasis];
  int info = 0, one = 1;

  // Compute LU decomposition.
  LAPACK_dgetrf(&nBasis, &nBasis, jac, &nBasis, pivots, &info);

  // Only fortran thinks that it is the transpose
  std::string transpose = "Transpose"; 
  char trans = transpose.c_str()[0];
  // Solve matrix equation given the LU decomposition of the transpose.
  LAPACK_dgetrs(&trans, &nBasis, &one, jac, &nBasis, pivots, 
      residue, &nBasis, &info);
  delete[] jac;
  if (info != 0) {
    std::cerr << "Error occured in inversion.\n";
    return -1.0;
  }

  // Update the solution.
  for (int i = 0; i < nBasis; i++) {
    sigma[i] += residue[i];
  }
  computeDerivatives();
  return getTotalResidue();
}

double kroneckerDelta(int i, int j) {
  if (i == j) return 1.0;
  return 0.0;
}

double* NewtonRaphson1D::getJacobian() {
  int nBasis = basis->getRank();
  const double* xi = basis->getAbscissas();
  double* jac = new double[nBasis*nBasis];
  for (int iRow = 0; iRow < nBasis; iRow++) {
    double x = map(xi[iRow]);
    // Compute the constant derivatives.
    double dfdsigma, dfdsigmap, dfdsigmapp;
    if ( (iRow != 0 || gType == NONE) 
         && (iRow != nBasis-1 || hType == NONE)) {
      dfdsigma =
        dfdu(x, sigma[iRow], sigmap[iRow], sigmapp[iRow]);
      dfdsigmap =
        dfdup(x, sigma[iRow], sigmap[iRow], sigmapp[iRow]);
      dfdsigmapp =
        dfdupp(x, sigma[iRow], sigmap[iRow], sigmapp[iRow]);
    } else if (iRow == 0 && gType != NONE) {
      dfdsigma = dgdu(sigma[iRow], sigmap[iRow]);
      dfdsigmap = dgdup(sigma[iRow], sigmap[iRow]);
      dfdsigmapp = 0.0;
    } else if (iRow == nBasis-1 && hType != NONE) {
      int i = iRow;
      if (hType == INITIAL) i = 0;
      dfdsigma = dhdu(sigma[i], sigmap[i]);
      dfdsigmap = dhdup(sigma[i], sigmap[i]);
      dfdsigmapp = 0.0;
    }
    for (int iCol = 0; iCol < nBasis; iCol++) {
      int iRowMod = iRow;
      if (iRow == nBasis-1 && hType == INITIAL) iRowMod == 0;
      int iMatrix = iRowMod*nBasis + iCol;
      // The aesthetically pleasing spectral Jacobian matrix:
      //  a logical sum of vector-Matrix row products.
      jac[iMatrix] =
        dfdsigma*kroneckerDelta(iRowMod,iCol) +
        dfdsigmap*diff[iMatrix] +
        dfdsigmapp*doubleDiff[iMatrix];
    }
  }
  return jac;
}

NewtonRaphson1D::~NewtonRaphson1D() {
  delete[] residue;
  delete[] sigmap;
  delete[] sigmapp;
  delete[] diff;
  delete[] doubleDiff;
}

double NewtonRaphson1D::g(double u, double up) {
  return 0.0;
}

double NewtonRaphson1D::dgdu(double u, double up) {
  return 0.0;
}

double NewtonRaphson1D::dgdup(double u, double up) {
  return 0.0;
}

double NewtonRaphson1D::h(double u, double up) {
  return 0.0;
}

double NewtonRaphson1D::dhdu(double u, double up) {
  return 0.0;
}

double NewtonRaphson1D::dhdup(double u, double up) {
  return 0.0;
}
