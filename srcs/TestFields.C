#include <iostream>
#include <cassert>
#include <ctime>

#include "Fields.h"

class TestFields: public Fields {
  public:
    TestFields(Hamiltonian &ham, ExtrinsicCurvature &kij) :
      Fields(ham, kij) {};
};

int main() {
  int nR = 11;
  int nTheta = 7;
  double rMax = 10;
  double tolerance = 1.5e-10;
  int maxIterations = 100;
  double bareMass = 1.;

  ExtrinsicCurvature kij;
  SphericalBasis* basis = new SphericalBasis();
  Hamiltonian ham(kij,basis);
  ham.setBareMass(bareMass);
  ham.setMaximumRadius(rMax);
  ham.setRanks(nR, nTheta);
  ham.solve(tolerance, maxIterations);
  TestFields test(ham, kij);

  std::cout << test.psi(0.5, 1.0) << "\n";
  std::cout << test.psi(1.0, 1.0) << "\n";
}
