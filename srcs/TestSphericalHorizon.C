#include <iostream>
#include <cassert>

#include "ChebyshevRoots.h"
#include "SphericalHorizon.h"

double abs(double x) {
  if (x < 0.) return -x;
  return x;
}

class TestSphericalHorizon: public SphericalHorizon {
  public:
    TestSphericalHorizon(Basis* aBasis, Fields* initialData,
        int n, double tol, int maxIter):
      SphericalHorizon(aBasis, initialData), nPoints(n),
      tolerance(tol), maxIterations(maxIter) {};
    void runTests();
  private:
    int nPoints;
    int maxIterations;
    double tolerance;
    double* runFindHorizonTest();
    void testHorizonMass(double* u);
    void testHorizonSpin(double* u);
};

void TestSphericalHorizon::runTests() {
  std::cout << "Testing SphericalHorizon object.\n";
  int nTests = 0;

  double* u = runFindHorizonTest();
  nTests++; std::cout << ".\n";

  testHorizonMass(u);
  nTests++; std::cout << ".\n";

  testHorizonSpin(u);
  delete[] u;
  nTests++; std::cout << ".\n";

  std::cout << "OK. Ran " << nTests << " successfully.\n";
}

void TestSphericalHorizon::testHorizonMass(double* u) {
  double mass = computeHorizonMass(u);
  assert(abs(mass-1.) < 9.0e-6);
}

void TestSphericalHorizon::testHorizonSpin(double* u) {
  double spin = computeHorizonSpin(u);
  // Measuring 0 is a horrible test... the coefficient may be wrong.
  assert(abs(spin) < 1.0e-16);
}

double* TestSphericalHorizon::runFindHorizonTest() {
  double* u = new double[nPoints];
  // Test the outer horizon finder.
  for (int i = 0; i < nPoints; i++) u[i] = 0.6;
  double chi = 1.0;
  computeSolution(u, &chi, tolerance, maxIterations);
  for (int i = 0; i < nPoints; i++) {
    assert( abs(u[i] - 0.5) < 2.0e-15);
  }
  // Test the inner horizon finder.
  chi = -1.0;
  for (int i = 0; i < nPoints; i++) u[i] = 0.6;
  computeSolution(u, &chi, tolerance, maxIterations);
  for (int i = 0; i < nPoints; i++) {
    assert( abs(u[i] - 0.5) < 2.0e-15);
  }
  return u;
}

int main() {
  std::cout.precision(16);
  int nR = 11;
  int nTheta = 7;
  double rMax = 10;
  double tolerance = 1.5e-15;
  int maxIterations = 100;
  double bareMass = 1.;
  int nHorizon = 15;

  // Create a Schwarzschild solution.
  ExtrinsicCurvature kij;
  SphericalBasis* sBasis = new SphericalBasis();
  Hamiltonian ham(kij, sBasis);
  ham.setBareMass(bareMass);
  ham.setMaximumRadius(rMax);
  ham.setRanks(nR, nTheta);
  ham.solve(tolerance, maxIterations);
  Fields* fields = new Fields(ham, kij);

  ChebyshevRoots* basis = new ChebyshevRoots();
  basis->setRank(nHorizon);
  
  TestSphericalHorizon test(basis, fields, nHorizon, tolerance, maxIterations);
  test.runTests();
}
