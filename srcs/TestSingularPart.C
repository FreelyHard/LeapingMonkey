#include <iostream>
#include <cassert>
#include <cmath>

#define private public
#include "SingularPart.h"
#include "ChebyshevRoots.h"

class TestSingularPart: public SingularPart {
  public:
    TestSingularPart(Basis* base, ExtrinsicCurvature &kij): 
      SingularPart(base, kij) {};
    void runTests();
  private:
    void testSphericallySymmetric();
    void testResidueAndOperator();
    void testConvergence();
};

double abs(double x) {
  if (x < 0.) return -x;
  return x;
}

void TestSingularPart::testConvergence() {
  curvature.setSpin(0.5);
  double* sigma = new double[basis->getRank()];
  bool converged = fillSigma(sigma);
  delete[] sigma;
  assert(converged);
}

void TestSingularPart::testResidueAndOperator() {
  int nTheta = basis->getRank();
  double* beta = new double[nTheta];
  double* sigma = new double[nTheta];
  double* theta = new double[nTheta];
  const double* u = basis->getAbscissas();
  for (int i = 0; i < nTheta; i++) {
    theta[i] = acos(-1.)*0.5*(u[i]+1.);
    beta[i] = -8.*(-25./4.*pow(cos(theta[i]), 2.) + 2.)*pow(cos(theta[i]), 14.);
    sigma[i] = pow(cos(theta[i]), 2.);
  }
  double* lap = getOperator(theta);
  double* residue = new double[nTheta];
  fillResidue(sigma, beta, lap, residue);
  for (int i = 0; i < nTheta; i++) assert(abs(residue[i]) < 1.e-10);
  delete[] residue;
  delete[] lap;
  delete[] sigma;
  delete[] beta;
  delete[] theta;
}

void TestSingularPart::testSphericallySymmetric() {
  double actualSigma = sqrt(1.5);
  double* sigma = new double[basis->getRank()];
  fillSigma(sigma);
  for (int i = 0; i < basis->getRank(); i++) 
    assert(abs(actualSigma - sigma[i]) < 1.0e-14);
  delete[] sigma;
}

void TestSingularPart::runTests() {
  std::cout << "Testing SingularPart class.\n";
  int nTests = 0;

  testSphericallySymmetric();
  std::cout << ".\n"; nTests++;

  testResidueAndOperator();
  std::cout << ".\n"; nTests++;

  testConvergence();
  std::cout << ".\n"; nTests++;

  std::cout << "Ok. Ran " << nTests << " tests.\n";
};

int main() {
  ExtrinsicCurvature kij;
  kij.setMass(1.0);
  ChebyshevRoots zenithalBasis;
  zenithalBasis.setRank(20);
  TestSingularPart test(&zenithalBasis, kij);
  test.runTests();
}
