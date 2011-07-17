#include <cmath>
#include <iostream>

#include "Fields.h"
#include "ChebyshevRoots.h"

Fields::Fields(Hamiltonian &ham, ExtrinsicCurvature kij) : k(kij) {
  // Get data from the Hamiltonian.
  mass = ham.getBareMass();
  sigma = ham.getSingularAngularPart();
  regPower = ham.getSingularPower();
  u = ham.getRemainder();
  // Set up the basis.
  int* ranks = ham.getRanks();
  basis.setRanks(ranks[0], ranks[1]);
  basis.setMaximumRadius(ham.getMaximumRadius());
  computeDerivatives();
  delete[] ranks;
}

double Fields::computeADMMass() {
  Basis* thetaBasis = basis.getBasis(SphericalBasis::COORD2);
  const double* xi = thetaBasis->getAbscissas();
  int nR = basis.getRank(SphericalBasis::COORD1);
  const double* r = basis.getCoord(SphericalBasis::COORD1);
  int nTheta = basis.getRank(SphericalBasis::COORD2);
  const double* theta = basis.getCoord(SphericalBasis::COORD2);
  double integrand[nTheta];
  for (int iTheta = 0; iTheta < nTheta; iTheta++) {
    // This formula holds because of the boundary condition.
    double myPsi = psi(r[nR-1], theta[iTheta]);
    integrand[iTheta] = (myPsi - 1.)*
      r[nR-1]*sin(theta[iTheta]);
    // Change to the spectral coordinate to integrate...
    integrand[iTheta] *= 0.5*M_PI*sqrt(1.0 - xi[iTheta]*xi[iTheta]);
  }
  return thetaBasis->integrate(integrand);
}

Fields::~Fields() {
  delete[] u;
  delete[] ur;
  delete[] ut;
  delete[] urr;
  delete[] urt;
  delete[] sigma;
  delete[] sigmap;
}

void Fields::computeDerivatives() {
  int nR = basis.getRank(SphericalBasis::COORD1);
  int nTheta = basis.getRank(SphericalBasis::COORD2);
  ur = new double[nR*nTheta];
  ut = new double[nR*nTheta];
  urr = new double[nR*nTheta];
  urt = new double[nR*nTheta];
  sigmap = new double[nTheta];

  const double* ddr = basis.getDiff(SphericalBasis::COORD1);
  const double* ddtheta = basis.getDiff(SphericalBasis::COORD2);

  for (int iDerivative = 0; iDerivative < nR*nTheta; iDerivative++) {
    ur[iDerivative] = 0.;
    ut[iDerivative] = 0.;
    for (int iFunction = 0; iFunction < nR*nTheta; iFunction++) {
      int iMatrix = basis.matrixIndex(iDerivative, iFunction);
      ur[iDerivative] += ddr[iMatrix]*u[iFunction];
      ut[iDerivative] += ddtheta[iMatrix]*u[iFunction];
    }
  }
  for (int iDerivative = 0; iDerivative < nR*nTheta; iDerivative++) {
    urr[iDerivative] = 0.;
    urt[iDerivative] = 0.;
    for (int iFunction = 0; iFunction < nR*nTheta; iFunction++) {
      int iMatrix = basis.matrixIndex(iDerivative, iFunction);
      urr[iDerivative] += ddr[iMatrix]*ur[iFunction];
      urt[iDerivative] += ddtheta[iMatrix]*ur[iFunction];
    }
  }

  // This is white boxish. If we go this far, we might as well not
  // compute it if it is zero.
  ChebyshevRoots thetaBasis;
  thetaBasis.setRank(nTheta);
  const double* ddxi = thetaBasis.getDifferentiationMatrix();
  for (int iDerivative = 0; iDerivative < nTheta; iDerivative++) {
    sigmap[iDerivative] = 0.;
    for (int iFunction = 0; iFunction < nTheta; iFunction++) {
      int iMatrix = thetaBasis.index(iDerivative, iFunction);
      // dxi/dtheta = 2/pi and this is chain rule
      sigmap[iDerivative] += M_2_PI*ddxi[iMatrix]*sigma[iFunction];
    }
  }
}

double Fields::psi(double r, double theta) {
  // These are kind of hackish.
  double* mySigma = 
    basis.interpolate1d(SphericalBasis::COORD2, sigma, &theta, 1);
  double* remainder = basis.tensorInterpolate(u, &r, 1, &theta, 1);
  double psi = 1. + 0.5*mass/r + mySigma[0]*pow(r, -regPower) + remainder[0];
  delete[] remainder;
  delete[] mySigma;
  return psi;
}

double Fields::psir(double r, double theta) {
  // These are kind of hackish.
  double* mySigma = 
    basis.interpolate1d(SphericalBasis::COORD2, sigma, &theta,1);
  double* remainder = basis.tensorInterpolate(ur, &r, 1, &theta, 1);
  double psir = -0.5*mass/(r*r) - regPower*mySigma[0]*pow(r, -regPower-1.)
    + remainder[0];
  delete[] remainder;
  delete[] mySigma;
  return psir;
}

double Fields::psirr(double r, double theta) {
  // These are kind of hackish.
  double* mySigma = 
    basis.interpolate1d(SphericalBasis::COORD2, sigma, &theta,1);
  double* remainder = basis.tensorInterpolate(urr, &r, 1, &theta, 1);
  double psirr = mass/(r*r*r)
    + regPower*(regPower + 1.)*mySigma[0]*pow(r, -regPower-2.)
    + remainder[0];
  delete[] remainder;
  delete[] mySigma;
  return psirr;
}

double Fields::psit(double r, double theta) {
  // These are kind of hackish.
  double* mySigma = 
    basis.interpolate1d(SphericalBasis::COORD2, sigmap, &theta,1);
  double* remainder = basis.tensorInterpolate(ut, &r, 1, &theta, 1);
  double psit = mySigma[0]*pow(r, -regPower) + remainder[0];
  delete[] remainder;
  delete[] mySigma;
  return psit;
}

double Fields::psirt(double r, double theta) {
  // These are kind of hackish.
  double* mySigma = 
    basis.interpolate1d(SphericalBasis::COORD2, sigmap, &theta,1);
  double* remainder = basis.tensorInterpolate(urt, &r, 1, &theta, 1);
  double psirt = - regPower*mySigma[0]*pow(r, -regPower-1)
    + remainder[0];
  delete[] remainder;
  delete[] mySigma;
  return psirt;
}

double Fields::Axx(double r, double theta) {
  return k.Axx(r,theta);
}

double Fields::Axy(double r, double theta) {
  return k.Axy(r,theta);
}

double Fields::Ayz(double r, double theta) {
  return k.Ayz(r,theta);
}

double Fields::Axz(double r, double theta) {
  return k.Axz(r,theta);
}

double Fields::Azz(double r, double theta) {
  return k.Azz(r,theta);
}

double Fields::Axxr(double r, double theta) {
  return k.Axxr(r,theta);
}

double Fields::Axzr(double r, double theta) {
  return k.Axzr(r,theta);
}

double Fields::Azzr(double r, double theta) {
  return k.Azzr(r,theta);
}
