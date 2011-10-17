#include "SpectralHorizon.h"
#include "SpectralHorizonMacros.h"
#include "NewtonRaphson1D.h"

#include "lapacke.h"
#include "cblas.h"

#include <math.h>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

SpectralHorizon::SpectralHorizon(double mass, double dipole, double quadrupole,
    Legendre* aBasis):
        psi(mass, dipole, quadrupole), nBasis(0), 
        basis(aBasis), NewtonRaphson1D(aBasis)
{
};

void SpectralHorizon::printHorizon() {
  const double* xi = basis->getAbscissas();
  double totalResidue = getTotalResidue();
  const double* residue = getResidue();
  for (int i = 0; i < nBasis; i++) {
    printf("%10.8e %10.8e %10.8e\n",
        acos(xi[i]), horizonRadius[i], residue[i]);
  }
  printf("Total residual error = %16.14f\n", totalResidue);
}

double SpectralHorizon::dInverseMap(double xi) {
  return -sqrt(1.0 - xi*xi);
}

void SpectralHorizon::printToFile() {
  std::ostringstream o;
  o << "m=" << psi.getMass() << "_d=" << psi.getDipoleMoment() << "_q=" <<
   psi.getQuadrupoleMoment() << ".dat";
  std::string filename = o.str();
  std::cout << "Writing to " + filename + ".\n";
  std::ofstream myFile;
  myFile.open(filename.c_str());
  myFile << "# theta\tradius\tresidual\n";
  myFile.precision(14);
  const double* xi = basis->getAbscissas();
  const double* residue = getResidue();
  for (int i = 0; i < nBasis; i++) {
    myFile << acos(xi[i]) << " " << horizonRadius[i] << " "
      << residue[i] << "\n";
  }
  myFile.close();
}

void SpectralHorizon::interpolateToFile(int N) {
  std::ostringstream o;
  o << "m=" << psi.getMass() << "_d=" << psi.getDipoleMoment() << "_q=" <<
   psi.getQuadrupoleMoment() << "N=" << nBasis << ".dat";
  std::string filename = o.str();
  std::cout << "Writing to " + filename + ".\n";
  std::ofstream myFile;
  myFile.open(filename.c_str());
  myFile << "# theta\tradius\tresidual\n";
  myFile.precision(14);

  double h = acos(-1.0)/(double)(N-1);
  double x[N];
  for (int i = 0; i < N; i++) x[i] = cos(h*(double)i);

  const double* residue = getResidue();
  double* interpSigma = basis->interpolate(horizonRadius,x,N);
  double* interpResidue = basis->interpolate(residue,x,N);
  for (int i = 0; i < N; i++) {
    myFile << acos(x[i]) << " " << interpSigma[i] << " "
      << interpResidue[i] << "\n";
  }
  myFile.close();
  delete[] interpSigma;
  delete[] interpResidue;
}

bool SpectralHorizon::initializeRepresentation(int n) {
  // Free up memory if we have already set the basis->
  if (nBasis != 0) {
    freeMemory();
  }
  bool basisMade = basis->setRank(n);
  if (basisMade) {
    nBasis = n;
    horizonRadius = new double[nBasis];
  } else {
    nBasis = 0;
  }
  return basisMade;
}

// Garbage collection.
void SpectralHorizon::freeMemory() {
  if (nBasis != 0) {
    delete[] horizonRadius;
  }
}

SpectralHorizon::~SpectralHorizon() {
  freeMemory();
  delete basis;
}

double SpectralHorizon::f(double xi, double u, double up, double upp) {
  double v = sqrt(1. - xi*xi);
  return phi(u * v, u * xi) * (-(-0.2e1 * up * up + u * upp - u * u) / (up * up + u * u) + (-up * xi + u * v) / u / v) / 0.4e1 + (-up * xi + u * v) * phir(u * v, u * xi) + (up * v + u * xi) * phiz(u * v, u * xi);
}

double SpectralHorizon::dfdu(double xi, double u, double up, double upp) {
  double v = sqrt(1. - xi*xi);
  return (v * phir(u * v, u * xi) + xi * phiz(u * v, u * xi)) * (-(-0.2e1 * up * up + u * upp - u * u) / (up * up + u * u) + (-up * xi + u * v) / u / v) / 0.4e1 + phi(u * v, u * xi) * (-(upp - 0.2e1 * u) / (up * up + u * u) + 0.2e1 * (-0.2e1 * up * up + u * upp - u * u) * pow(up * up + u * u, -0.2e1) * u + 0.1e1 / u - (-up * xi + u * v) * pow(u, -0.2e1) / v) / 0.4e1 + v * phir(u * v, u * xi) + (-up * xi + u * v) * (phirr(u * v, u * xi) * v + phirz(u * v, u * xi) * xi) + xi * phiz(u * v, u * xi) + (up * v + u * xi) * (phirz(u * v, u * xi) * v + phizz(u * v, u * xi) * xi);
}

double SpectralHorizon::dfdup(double xi, double u, double up, double upp) {
  double v = sqrt(1. - xi*xi);
  return phi(u * v, u * xi) * (0.4e1 * up / (up * up + u * u) + 0.2e1 * (-0.2e1 * up * up + u * upp - u * u) * pow(up * up + u * u, -0.2e1) * up - xi / u / v) / 0.4e1 - xi * phir(u * v, u * xi) + v * phiz(u * v, u * xi);
}

double SpectralHorizon::dfdupp(double xi, double u, double up, double upp) {
  double v = sqrt(1. - xi*xi);
  return -phi(u * v, u * xi) / (up * up + u * u) * u / 0.4e1;
}

void SpectralHorizon::setGuess(double* guess) {
  for (int i = 0; i < nBasis; i++) horizonRadius[i] = guess[i];
}

bool SpectralHorizon::findHorizon(double tolerance, int maxIterations) {
  return computeSolution(horizonRadius, NULL, tolerance, maxIterations);
}


void SpectralHorizon::setMoments(double m, double d, double q) {
  psi.setMoments(m, d, q);
}

double SpectralHorizon::getHorizonMass() {
  double* sigma = horizonRadius;
  const double* sigmap = getDerivative();
  double integrand[nBasis];
  const double* xi = basis->getAbscissas();
  // A = 2\pi\int{psi^4*r*ds/d\theta d\theta}^\pi_0
  // ds/d\theta = \sqrt{\sigma^2 + \sigma'^2}
  for (int i = 0; i < nBasis; i++) {
    double theta = acos(xi[i]);
    double r = sigma[i]*sin(theta);
    double z = sigma[i]*cos(theta);
    double psi4 = pow(psi.psi(r,z), 4.0);
    integrand[i] = psi4*sigma[i]*
      sqrt(sigma[i]*sigma[i] + sigmap[i]*sigmap[i]);
  }
  return sqrt(basis->integrate(integrand)/8.0);
}
