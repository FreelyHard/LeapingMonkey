#include <cmath>

#include "SphericalHorizon.h"

// Macros, so that Maple output works directly.
#define Axx fields->Axx
#define Axz fields->Axz
#define Azz fields->Azz
#define Axxr fields->Axxr
#define Axzr fields->Axzr
#define Azzr fields->Azzr
#define psi fields->psi
#define psir fields->psir
#define psit fields->psit
#define psirr fields->psirr
#define psirt fields->psirt

SphericalHorizon::SphericalHorizon(Basis* aBasis, Fields* initialData):
  NewtonRaphson1D(aBasis), fields(initialData)
{
}

bool SphericalHorizon::computeSolution(double* u, void* selector,
    double tolerance, int maxIterations) {
  setPointerToSolution(u);
  chi = *((double*)selector);
  return NewtonRaphson1D::computeSolution(tolerance, maxIterations);
}

double SphericalHorizon::f(double theta, double u, double up, double upp) {
  return psi(u, theta) * (0.2e1 - up / u * 0.1e1 / tan(theta) - 0.1e1 / (up * up + u * u) * upp * u) / 0.4e1 + chi * (Axx(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * pow(-u * sin(theta) + cos(theta) * up, 0.2e1) - 0.2e1 * Axz(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * (-sin(theta) * u * u * cos(theta) - u * up + 0.2e1 * pow(cos(theta), 0.2e1) * up * u + cos(theta) * up * up * sin(theta)) + Azz(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * pow(u * cos(theta) + sin(theta) * up, 0.2e1)) * pow(psi(u, theta), -0.3e1) / 0.4e1 + u * psir(u, theta) - up / u * psit(u, theta);
}

double SphericalHorizon::dfdu(double theta, double u, double up, double upp) {
  return psir(u, theta) * (0.2e1 - up / u * 0.1e1 / tan(theta) - 0.1e1 / (up * up + u * u) * upp * u) / 0.4e1 + psi(u, theta) * (up * pow(u, -0.2e1) * 0.1e1 / tan(theta) + 0.2e1 * pow(up * up + u * u, -0.2e1) * upp * u * u - 0.1e1 / (up * up + u * u) * upp) / 0.4e1 + chi * (Axxr(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * pow(-u * sin(theta) + cos(theta) * up, 0.2e1) - Axx(u, theta) * pow(up * up + u * u, -0.3e1 / 0.2e1) * pow(-u * sin(theta) + cos(theta) * up, 0.2e1) * u - 0.2e1 * Axx(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * (-u * sin(theta) + cos(theta) * up) * sin(theta) - 0.2e1 * Axzr(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * (-sin(theta) * u * u * cos(theta) - u * up + 0.2e1 * pow(cos(theta), 0.2e1) * up * u + cos(theta) * up * up * sin(theta)) + 0.2e1 * Axz(u, theta) * pow(up * up + u * u, -0.3e1 / 0.2e1) * (-sin(theta) * u * u * cos(theta) - u * up + 0.2e1 * pow(cos(theta), 0.2e1) * up * u + cos(theta) * up * up * sin(theta)) * u - 0.2e1 * Axz(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * (-0.2e1 * sin(theta) * u * cos(theta) - up + 0.2e1 * pow(cos(theta), 0.2e1) * up) + Azzr(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * pow(u * cos(theta) + sin(theta) * up, 0.2e1) - Azz(u, theta) * pow(up * up + u * u, -0.3e1 / 0.2e1) * pow(u * cos(theta) + sin(theta) * up, 0.2e1) * u + 0.2e1 * Azz(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * (u * cos(theta) + sin(theta) * up) * cos(theta)) * pow(psi(u, theta), -0.3e1) / 0.4e1 - 0.3e1 / 0.4e1 * chi * (Axx(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * pow(-u * sin(theta) + cos(theta) * up, 0.2e1) - 0.2e1 * Axz(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * (-sin(theta) * u * u * cos(theta) - u * up + 0.2e1 * pow(cos(theta), 0.2e1) * up * u + cos(theta) * up * up * sin(theta)) + Azz(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * pow(u * cos(theta) + sin(theta) * up, 0.2e1)) * pow(psi(u, theta), -0.4e1) * psir(u, theta) + psir(u, theta) + u * psirr(u, theta) + up * pow(u, -0.2e1) * psit(u, theta) - up / u * psirt(u, theta);
}

double SphericalHorizon::dfdup(double theta, double u, double up, double upp) {
  return psi(u, theta) * (-0.1e1 / u * 0.1e1 / tan(theta) + 0.2e1 * pow(up * up + u * u, -0.2e1) * upp * u * up) / 0.4e1 + chi * (-Axx(u, theta) * pow(up * up + u * u, -0.3e1 / 0.2e1) * pow(-u * sin(theta) + cos(theta) * up, 0.2e1) * up + 0.2e1 * Axx(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * (-u * sin(theta) + cos(theta) * up) * cos(theta) + 0.2e1 * Axz(u, theta) * pow(up * up + u * u, -0.3e1 / 0.2e1) * (-sin(theta) * u * u * cos(theta) - u * up + 0.2e1 * pow(cos(theta), 0.2e1) * up * u + cos(theta) * up * up * sin(theta)) * up - 0.2e1 * Axz(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * (-u + 0.2e1 * pow(cos(theta), 0.2e1) * u + 0.2e1 * cos(theta) * up * sin(theta)) - Azz(u, theta) * pow(up * up + u * u, -0.3e1 / 0.2e1) * pow(u * cos(theta) + sin(theta) * up, 0.2e1) * up + 0.2e1 * Azz(u, theta) * pow(up * up + u * u, -0.1e1 / 0.2e1) * (u * cos(theta) + sin(theta) * up) * sin(theta)) * pow(psi(u, theta), -0.3e1) / 0.4e1 - 0.1e1 / u * psit(u, theta);
}

double SphericalHorizon::dfdupp(double theta, double u, double up, double upp) {
  return -psi(u, theta) / (up * up + u * u) * u / 0.4e1;
}

double SphericalHorizon::map(double xi) {
  return (xi + 1.0)*0.5*M_PI;
}

double SphericalHorizon::dInverseMap(double xi) {
  return M_2_PI;
}

double SphericalHorizon::computeHorizonMass(double* r) {
  // Set up some variables.
  int nBasis = basis->getRank();
  double integrand[nBasis];
  const double *xi = basis->getAbscissas();
  const double *diff = basis->getDifferentiationMatrix();
  double rPrime[nBasis];

  // compute the derivative.
  for (int iDerivative = 0; iDerivative < nBasis; iDerivative++) {
    rPrime[iDerivative] = 0.;
    for (int iFunction = 0; iFunction < nBasis; iFunction++) {
      int index = basis->index(iDerivative, iFunction);
      rPrime[iDerivative] += diff[index]*r[iFunction];
    }
  }

  // Compute the integrand. See ApparentHorizon.pdf
  for (int i = 0; i < nBasis; i++) {
    double theta = map(xi[i]);
    double myPsi = psi(r[i], theta);
    integrand[i] = 2.0*M_PI*pow(psi(r[i], theta), 4.)*r[i]*sin(theta)*
      sqrt(rPrime[i]*rPrime[i] + r[i]*r[i]);
    // Need to convert to Chebyshev coordinates and multiply by weight
    integrand[i] *= 0.5*M_PI*sqrt(1. - xi[i]*xi[i]);
  }

  return sqrt(basis->integrate(integrand)*M_1_PI*0.0625);
}

double SphericalHorizon::computeHorizonSpin(double* r) {
  // Set up some variables.
  int nBasis = basis->getRank();
  double integrand[nBasis];
  const double *xi = basis->getAbscissas();
  const double *diff = basis->getDifferentiationMatrix();
  double rPrime[nBasis];

  // compute the derivative.
  for (int iDerivative = 0; iDerivative < nBasis; iDerivative++) {
    rPrime[iDerivative] = 0.;
    for (int iFunction = 0; iFunction < nBasis; iFunction++) {
      int index = basis->index(iDerivative, iFunction);
      rPrime[iDerivative] += diff[index]*r[iFunction];
    }
  }

  // Compute the integrand. See ApparentHorizon.pdf
  for (int i = 0; i < nBasis; i++) {
    double theta = map(xi[i]);
    double axy = fields->Axy(r[i], theta);
    double ayz = fields->Ayz(r[i], theta);
    double sTheta = sin(theta);
    double cTheta = cos(theta);
    integrand[i] = -0.25*r[i]*r[i]*(
             r[i]*(sTheta*axy + cTheta*ayz) -
        rPrime[i]*(cTheta*axy - sTheta*ayz))*sTheta*sTheta;
    // Need to convert to Chebyshev coordinates and multiply by weight
    integrand[i] *= 0.5*M_PI*sqrt(1. - xi[i]*xi[i]);
  }

  return basis->integrate(integrand);
}
