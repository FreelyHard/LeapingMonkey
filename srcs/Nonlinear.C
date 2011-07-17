#include <iostream>
#include <cmath>

#include "Hamiltonian.h"

double differenceBetweenSolutions(double* high, double* low, int n) {
  double euclideanNorm = 0.0;
  for (int i = 0; i < n; i++) {
    euclideanNorm += (high[i] - low[i])*(high[i] - low[i]);
  }
  return sqrt(euclideanNorm);
}


int main() {
  // Parameters...
  // Numerical...
  int nSimple = 65;
  int nRadial = nSimple;
  int nTheta = nSimple;
  double rMax = 10;
  double tolerance = 1.5e-10;
  int maxIterations = 10000;
  bool verbose = true;
  // Physical...
  double bareMass = 1.;
  double trumpetMass = 0.;
  double momentum = 0.;
  double spin = 0.25;

  std::cout << "Creating basic extrinsic curvature and hamiltonian...\n";
  ExtrinsicCurvature kij;
  kij.setMass(trumpetMass);
  kij.setMomentum(momentum);
  kij.setSpin(spin);
  std::cout << "Creating basic Hamiltonian.\n";
  SphericalBasis* basis = new SphericalBasis();
  Hamiltonian ham(kij, basis);
  ham.setMaximumRadius(rMax);
  ham.setBareMass(bareMass);
  ham.setVerbose(verbose);

  int n1 = nRadial;
  int n2 = nTheta;
 
  ham.setRanks(nRadial, nTheta);
  std::cout << "Solving Hamiltonian.\n";
  ham.solve(tolerance, maxIterations);
  ham.printToFile();

  delete basis;
}
