#include <iostream>
#include <cassert>

#define protected public
#include "SpectralHorizon.h"

double abs( double x) {
  if (x < 0) return -x;
  return x;
}

int main() {
  printf("Testing SpectralHorizon object.\n");
  Legendre* basis = new Legendre();
  SpectralHorizon sphere(1.0,0.0,0.0,basis);
  int nBasis = 5;
  sphere.initializeRepresentation(nBasis);
  double guess[nBasis];
  for (int i = 0; i < nBasis; i++) guess[i] = 0.5;
  sphere.setGuess(guess);
  sphere.findHorizon(1.0e-12, 10);
  const double* residue = sphere.getResidue();
  for (int i = 0; i < nBasis; i++) assert(abs(residue[i]) < 7.0e-15);
  printf(".\n");
  printf("OK.\n");
}
