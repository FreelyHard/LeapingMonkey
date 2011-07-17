#include "SpectralHorizon.h"
#include <iostream>
#include <cstdlib>
#include <assert.h>
#include <fstream>
#include <math.h>

double abs( double x) {
  if (x < 0) return -x;
  return x;
}

bool iterate(SpectralHorizon* horizon, double momentMin,
    double momentMax, double momentStep, double mass, double otherMoment,
    void (*momentSelect)(SpectralHorizon*, double, double, double),
    char* filename) {
  double tolerance = 1.0e-12;
  double maxIterations = 100;

  std::ofstream massFile;
  massFile.open(filename, std::ios::app);
  massFile << "# moment horizonMass\n";
  massFile.precision(14);
  massFile.flush();

  printf("---> Iterating over moments...\n");
  int N = 1000;
  for (double moment = momentMin; 
      abs(moment) < abs(momentMax); 
      moment += momentStep) {
    printf("Setting moment to %g.\n",moment);
    // Set moments and compute.
    momentSelect(horizon, mass, moment, otherMoment);
    bool converged = horizon->findHorizon(tolerance, maxIterations);
    if (converged) {
      printf("Solution found.\n");
      horizon->interpolateToFile(1000);
      massFile << moment << " " << horizon->getHorizonMass() << "\n";
      massFile.flush();
    } else {
      printf("Solution failed.\n");
      return false;
    }
  }
  massFile.close();
  return true;
}

void changeQuadrupole(SpectralHorizon* horizon, double mass, 
    double quadrupole, double dipole) {
  horizon->setMoments(mass, dipole, quadrupole);
}

void changeDipole(SpectralHorizon* horizon, double mass,
    double dipole, double quadrupole) {
  horizon->setMoments(mass, dipole, quadrupole);
}

int main() {
  printf("Building Schwarzschild solution.\n");
  double m = 1.0;
  double dipole = 0.0;
  double quadrupole = 0.0;
  Legendre* basis = new Legendre();
  SpectralHorizon* horizon = new SpectralHorizon(m,dipole,quadrupole,basis);
  int nBasis = 20;
  horizon->initializeRepresentation(nBasis);
  double guess[nBasis];
  for (int i = 0; i < nBasis; i++) guess[i] = 0.6;
  horizon->setGuess(guess);
  printf("Mass = %f\n", horizon->getHorizonMass());

  iterate(horizon, 0.0, 0.2, 0.0025, m, quadrupole, &changeDipole,
     "/dev/null");
  horizon->setGuess(guess);
  /*iterate(horizon, 0.0, 0.13, 0.01, m, dipole, &changeQuadrupole,
    "quadrupole.mass");
  horizon->setGuess(guess);
  iterate(horizon, 0.0, -0.19, -0.01, m, dipole, &changeQuadrupole,
    "quadrupole.mass");
*/
  delete horizon;
}
