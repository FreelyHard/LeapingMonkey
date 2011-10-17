#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cassert>
#include <ctime>

#define protected public
#define private public
#include "Hamiltonian.h"
#include "lapacke.h"
#define N1 13
#define N2 7
#define REGPOWER 2.

double abs(double x) {
  if (x < 0.0) return -x;
  return x;
}

class TestHamiltonian: public Hamiltonian {
  public:
    TestHamiltonian(ExtrinsicCurvature &kij, SphericalBasis* aBasis):
      Hamiltonian(kij,aBasis) {};
    void runTests();
    void showSpectra();
  private:
    void runOperatorKernelTest();
    void runDoubleFreeTest();
    void runJacobianTest();
    void runResidueTest();
};

void TestHamiltonian::runResidueTest() {
  // Set up the fake test data.
  int nRadial = basis->getRank(SphericalBasis::COORD1);
  int nTheta = basis->getRank(SphericalBasis::COORD2);
  const double* r = basis->getCoord(SphericalBasis::COORD1);
  const double* theta = basis->getCoord(SphericalBasis::COORD2);
  delete[] sigma; delete[] u; delete[] kappa; delete[] betaOverSigma7;
  sigma = new double[nTheta];
  betaOverSigma7 = new double[nTheta];
  u = new double[nTheta*nRadial];
  kappa = new double[nTheta*nRadial];
  regPower = 0.;
  mass = 0.0;
  for (int iTheta = 0; iTheta < nTheta; iTheta++) {
    sigma[iTheta] = -5. + 1.0*(double)rand()/(double)RAND_MAX;
    betaOverSigma7[iTheta] = -5. + 1.0*(double)rand()/(double)RAND_MAX;
    for (int iR = 0; iR < nRadial; iR++) {
      int iFunc = basis->functionIndex(iR, iTheta);
      u[iFunc] = -1.;
      kappa[iFunc] = betaOverSigma7[iTheta];
    }
    betaOverSigma7[iTheta] *= pow(sigma[iTheta], -7.);
  }
  int nPoints = nRadial*nTheta;
  double residue[nPoints];
  double totalResidue = fillResidue(residue);
  for (int iR = 0; iR < nRadial - 1; iR++) {
    for (int iTheta = 0; iTheta < nTheta; iTheta++) {
      int iResidue = basis->functionIndex(iR, iTheta);
      assert(abs(residue[iResidue]) < 1.0e-11);
    }
  }
}

void TestHamiltonian::runJacobianTest() {
  // Set up the fake test data.
  int nRadial = basis->getRank(SphericalBasis::COORD1);
  int nTheta = basis->getRank(SphericalBasis::COORD2);
  const double* r = basis->getCoord(SphericalBasis::COORD1);
  const double* theta = basis->getCoord(SphericalBasis::COORD2);
  delete[] sigma; delete[] u; delete[] kappa; delete[] betaOverSigma7;
  sigma = new double[nTheta];
  betaOverSigma7 = new double[nTheta];
  u = new double[nTheta*nRadial];
  kappa = new double[nTheta*nRadial];
  for (int iTheta = 0; iTheta < nTheta; iTheta++) {
    sigma[iTheta] = 0.;
    betaOverSigma7[iTheta] = 0.;
    for (int iR = 0; iR < nRadial; iR++) {
      int iFunc = basis->functionIndex(iR, iTheta);
      u[iFunc] = r[iR];
      kappa[iFunc] = 0.0;
    }
  }
  regPower = 0.;
  int nPoints = nRadial*nTheta;
  double jacobian[nPoints*nPoints];
  fillJacobian(jacobian);
  // This part tests the boundary condition.
  for (int iTheta = 0; iTheta < nTheta; iTheta++) {
    int iR = nRadial-1;
    int iGrad = basis->functionIndex(iR, iTheta);
    double gradient = 0;
    for (int iu = 0; iu < nPoints; iu++) {
      int iTranspose = basis->matrixIndex(iu, iGrad);
      gradient += jacobian[iTranspose]*u[iu];
    }
    assert(abs(2. - gradient) < 5.0e-14);
  }
  // This part tests the laplacian part.
  const double* lap = basis->getLaplacian();
  for (int igR = 0; igR < nRadial-1; igR++) {
    for (int iT = 0; iT < nTheta; iT++) {
      int iGrad = basis->functionIndex(igR, iT);
      double gradient = 0.;
      for (int iTheta = 0; iTheta < nTheta; iTheta++) {
        for (int iR = 0; iR < nRadial; iR++) {
          int iu = basis->functionIndex(iR, iTheta);
          int iTranspose = basis->matrixIndex(iu, iGrad);
          gradient += jacobian[iTranspose]*r[iR]*r[iR];
        }
      }
      assert(abs(gradient/pow(r[igR], regPower+2.)-6) < 5.e-12);
    }
  }
}

void TestHamiltonian::runTests() {
  std::cout << "Testing hamiltonian equation.\n";
  int nTests = 0;

  runDoubleFreeTest();
  std::cout << ".\n"; nTests++;

  runJacobianTest();
  std::cout << ".\n"; nTests++;

  runResidueTest();
  std::cout << ".\n"; nTests++;

  std::cout << "Ok. Ran " << nTests << " tests successfully.\n";
}

void TestHamiltonian::runDoubleFreeTest() {
  freeFields();
}

int main() {
  srand(time(NULL));
  ExtrinsicCurvature kij;
  SphericalBasis* basis = new SphericalBasis();
  TestHamiltonian test(kij, basis);
  test.setRanks(N1,N2);
  test.runTests();
}
