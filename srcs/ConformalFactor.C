#include "ConformalFactor.h"
#include <math.h>

double pow(double x, int n) {
  double y = 1;
  for (int i = 0; i<n; i++) y *= x;
  return y;
}

double ConformalFactor::getMass() {
  return m;
}

double ConformalFactor::getQuadrupoleMoment() {
  return q;
}

double ConformalFactor::getDipoleMoment() {
  return d;
}

void ConformalFactor::setMoments(double mass, double dipole, 
    double quadrupole) {
  m = mass;
  d = dipole;
  q = quadrupole;
}

double ConformalFactor::psi(double r, double z) const {
  return 0.1e1 + m * pow(r * r + z * z, -0.1e1 / 0.2e1) / 0.2e1 + d * z * pow(r * r + z * z, -0.3e1 / 0.2e1) / 0.2e1 + q * (r * r - 0.2e1 * z * z) * pow(r * r + z * z, -0.5e1 / 0.2e1) / 0.8e1;
}

double ConformalFactor::dpsidr(double r, double z) const {
  return -m * pow(r * r + z * z, -0.3e1 / 0.2e1) * r / 0.2e1 - 0.3e1 / 0.2e1 * d * z * pow(r * r + z * z, -0.5e1 / 0.2e1) * r + q * r * pow(r * r + z * z, -0.5e1 / 0.2e1) / 0.4e1 - 0.5e1 / 0.8e1 * q * (r * r - 0.2e1 * z * z) * pow(r * r + z * z, -0.7e1 / 0.2e1) * r;
}

double ConformalFactor::dpsidz(double r, double z) const {
  return -m * pow(r * r + z * z, -0.3e1 / 0.2e1) * z / 0.2e1 + d * pow(r * r + z * z, -0.3e1 / 0.2e1) / 0.2e1 - 0.3e1 / 0.2e1 * d * z * z * pow(r * r + z * z, -0.5e1 / 0.2e1) - q * z * pow(r * r + z * z, -0.5e1 / 0.2e1) / 0.2e1 - 0.5e1 / 0.8e1 * q * (r * r - 0.2e1 * z * z) * pow(r * r + z * z, -0.7e1 / 0.2e1) * z;
}

double ConformalFactor::d2psidr2(double r, double z) const {
  return 0.3e1 / 0.2e1 * m * pow(r * r + z * z, -0.5e1 / 0.2e1) * r * r - m * pow(r * r + z * z, -0.3e1 / 0.2e1) / 0.2e1 + 0.15e2 / 0.2e1 * d * z * pow(r * r + z * z, -0.7e1 / 0.2e1) * r * r - 0.3e1 / 0.2e1 * d * z * pow(r * r + z * z, -0.5e1 / 0.2e1) + q * pow(r * r + z * z, -0.5e1 / 0.2e1) / 0.4e1 - 0.5e1 / 0.2e1 * q * r * r * pow(r * r + z * z, -0.7e1 / 0.2e1) + 0.35e2 / 0.8e1 * q * (r * r - 0.2e1 * z * z) * pow(r * r + z * z, -0.9e1 / 0.2e1) * r * r - 0.5e1 / 0.8e1 * q * (r * r - 0.2e1 * z * z) * pow(r * r + z * z, -0.7e1 / 0.2e1);
}

double ConformalFactor::d2psidrdz(double r, double z) const {
  return 0.3e1 / 0.2e1 * m * pow(r * r + z * z, -0.5e1 / 0.2e1) * r * z - 0.3e1 / 0.2e1 * d * pow(r * r + z * z, -0.5e1 / 0.2e1) * r + 0.15e2 / 0.2e1 * d * z * z * pow(r * r + z * z, -0.7e1 / 0.2e1) * r + 0.5e1 / 0.4e1 * q * r * pow(r * r + z * z, -0.7e1 / 0.2e1) * z + 0.35e2 / 0.8e1 * q * (r * r - 0.2e1 * z * z) * pow(r * r + z * z, -0.9e1 / 0.2e1) * r * z;
}

double ConformalFactor::d2psidz2(double r, double z) const {
  return 0.3e1 / 0.2e1 * m * pow(r * r + z * z, -0.5e1 / 0.2e1) * z * z - m * pow(r * r + z * z, -0.3e1 / 0.2e1) / 0.2e1 - 0.9e1 / 0.2e1 * d * z * pow(r * r + z * z, -0.5e1 / 0.2e1) + 0.15e2 / 0.2e1 * d * pow(z, 0.3e1) * pow(r * r + z * z, -0.7e1 / 0.2e1) - q * pow(r * r + z * z, -0.5e1 / 0.2e1) / 0.2e1 + 0.5e1 * q * z * z * pow(r * r + z * z, -0.7e1 / 0.2e1) + 0.35e2 / 0.8e1 * q * (r * r - 0.2e1 * z * z) * pow(r * r + z * z, -0.9e1 / 0.2e1) * z * z - 0.5e1 / 0.8e1 * q * (r * r - 0.2e1 * z * z) * pow(r * r + z * z, -0.7e1 / 0.2e1);
}
