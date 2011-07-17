#include <iostream>
#include <cstdlib>
#include <ctime>
#include <assert.h>

#include "ConformalFactor.h"

double abs(double x) {
  if (x < 0.) return -x;
  return x;
}

int main() {
  printf("Testing conformal factor and Apparent horizon equations.\n");
  ConformalFactor psi(1.0,0.0, 0.0);
  double r0 = 0.5;
  double z0 = 0.0;
  double curvature = 4.0*psi.dpsidr(r0, z0)/psi.psi(r0, z0);
  curvature = -curvature - 1./r0;
  assert(curvature == 2.0);
  printf(".\n");

  srand(time(NULL));
  psi.setMoments(1.0, 0.5, 0.5);
  for (int i = 0; i < 100; i++) {
    double r = 0.001 + (double)rand()/(double)RAND_MAX*10.0;
    double z = -5.0 + (double)rand()/(double)RAND_MAX*10.0;
    double psirr = psi.d2psidr2(r,z);
    double psirOverr = psi.dpsidr(r,z)/r;
    double psizz = psi.d2psidz2(r,z);
    double lap = psirr + psizz + psirOverr;
    assert(abs(lap) < 1.0e-12); // This will sometimes randomly fail if one gets too close to the singularity.
  }
  printf(".\n");
  printf("OK.\n");
}
