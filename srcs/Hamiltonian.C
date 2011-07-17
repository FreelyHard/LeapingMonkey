#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "mkl_cblas.h"
#include "mkl_lapack.h"

#include "Hamiltonian.h"
#include "SingularPart.h"

Hamiltonian::Hamiltonian(ExtrinsicCurvature &kij, SphericalBasis* aBasis) : 
  k(kij), diffOperator(NULL), mass(1.0), sigma(NULL),
  kappa(NULL), u(NULL), residue(NULL), betaOverSigma7(NULL),
  verbose(false), basis(aBasis), severe(true) {
}

int* Hamiltonian::getRanks() {
  int* ranks = new int[2];
  ranks[0] = basis->getRank(SphericalBasis::COORD1);
  ranks[1] = basis->getRank(SphericalBasis::COORD2);
  return ranks;
}

void Hamiltonian::printToFile() {
  std::ostringstream o;
  int nR = basis->getRank(SphericalBasis::COORD1);
  int nTheta = basis->getRank(SphericalBasis::COORD2);
  // Filename...
  o << "hamiltonian_NR=" << nR
    << "_NTheta=" << nTheta
    << "_m=" << mass
    << "_M=" << k.getMass()
    << "_S=" << k.getSpin()
    << "_P=" << k.getMomentum()
    << ".dat";
  std::string filename = o.str();
  std::cout << "Writing to " + filename + ".\n";
  std::ofstream myFile;
  myFile.open(filename.c_str()); 
  // The header contains number of collocations
  myFile << "# NR = " << nR
    << "\n# NTheta = " << nTheta << "\n";
  // The coordinates of the collocations
  const double* r = basis->getCoord(SphericalBasis::COORD1);
  myFile << "# Radius values:\n";
  for (int iR = 0; iR < nR; iR++) myFile << r[iR] << " ";
  myFile << "\n";
  const double* theta = basis->getCoord(SphericalBasis::COORD2);
  myFile << "# Theta values:\n";
  for (int iTheta = 0; iTheta < nTheta; iTheta++) {
    myFile << theta[iTheta] << " ";
  }
  myFile << "\n";
  // The function values at the collocation points.
  myFile << "# Remainder (u) values, radius changes along row:\n";
  for (int iTheta = 0; iTheta < nTheta; iTheta++) {
    for (int iR = 0; iR < nR; iR++) {
      int iu = basis->functionIndex(iR, iTheta);
      myFile << u[iu] << " ";
    }
    myFile << "\n";
  }
  // The coefficients of the function
  double coeffs[nR*nTheta];
  basis->fillCoefficients(u, coeffs);
  myFile << "# Coefficients values, order of radius basis changes along row:\n";
  for (int iTheta = 0; iTheta < nTheta; iTheta++) {
    for (int iR = 0; iR < nR; iR++) {
      int iCoeff = basis->functionIndex(iR, iTheta);
      myFile << coeffs[iCoeff] << " ";
    }
    myFile << "\n";
  }
  myFile.close();
}

double* Hamiltonian::interpolateRemainder(int nR, int nTheta) {
  // The coordinates of the collocations
  double r[nR];
  double dr = rMax/(double)(nR-1);
  for (int i = 0; i < nR; i++) r[i] = dr*i;
  double theta[nTheta];
  double dTheta = acos(-1.0)/(double)(nTheta-1);
  for (int i = 0; i < nTheta; i++) theta[i] = dTheta*i;
  return basis->tensorInterpolate(u, r, nR, theta, nTheta);
}

void Hamiltonian::interpolateToFile(int nR, int nTheta) {
  std::ostringstream o;
  int realnR = basis->getRank(SphericalBasis::COORD1);
  int realnTheta = basis->getRank(SphericalBasis::COORD2);
  // Filename...
  o << "hamiltonian_interpolated_NR=" << nR << "(" << realnR << ")"
    << "_NTheta=" << nTheta << "(" << realnTheta << ")"
    << "_m=" << mass
    << "_M=" << k.getMass()
    << "_S=" << k.getSpin()
    << "_P=" << k.getMomentum()
    << ".dat";
  std::string filename = o.str();
  std::cout << "Writing to " + filename + ".\n";
  std::ofstream myFile;
  myFile.open(filename.c_str()); 
  myFile.precision(16);
  // The header contains number of collocations
  myFile << "# NR = " << nR
    << "\n# NTheta = " << nTheta << "\n";
  // The coordinates of the collocations
  double r[nR];
  double dr = rMax/(double)(nR-1);
  for (int i = 0; i < nR; i++) r[i] = dr*i;
  myFile << "# Radius values:\n";
  for (int iR = 0; iR < nR; iR++) myFile << r[iR] << " ";
  myFile << "\n";
  double theta[nTheta];
  double dTheta = acos(-1.0)/(double)(nTheta-1);
  for (int i = 0; i < nTheta; i++) theta[i] = dTheta*i;
  myFile << "# Theta values:\n";
  for (int iTheta = 0; iTheta < nTheta; iTheta++) {
    myFile << theta[iTheta] << " ";
  }
  myFile << "\n";
  double* uInterp = basis->tensorInterpolate(u, r, nR, theta, nTheta);
  // The function values at the collocation points.
  myFile << "# Remainder (u) values, radius changes along row:\n";
  for (int iTheta = 0; iTheta < nTheta; iTheta++) {
    for (int iR = 0; iR < nR; iR++) {
      int iu = nR*iTheta + iR;
      myFile << uInterp[iu] << " ";
    }
    myFile << "\n";
  }
  delete[] uInterp;
}

double Hamiltonian::getMaximumRadius() {
  return rMax;
}

void Hamiltonian::setMaximumRadius(double radius) {
  rMax = radius;
  basis->setMaximumRadius(radius);
}

void Hamiltonian::setRanks(int N1, int N2) {
  basis->setRanks(N1,N2);
  freeFields();
}

double Hamiltonian::getSingularPower() {
  return regPower;
}

void Hamiltonian::freeFields() {
  delete[] diffOperator; diffOperator = NULL;
  delete[] kappa; kappa = NULL;
  delete[] u; u = NULL;
  delete[] residue; residue = NULL;
  delete[] sigma; sigma = NULL;
  delete[] betaOverSigma7; betaOverSigma7 = NULL;
}

Hamiltonian::~Hamiltonian() {
  freeFields();
}

const double* Hamiltonian::getOperator() {
  if (diffOperator == NULL) {
    // Get and copy in the Laplacian.
    int nRadial = basis->getRank(SphericalBasis::COORD1);
    int nTheta = basis->getRank(SphericalBasis::COORD2);
    diffOperator = new double[nRadial*nRadial*nTheta*nTheta];
    const double* laplacian = basis->getLaplacian();
    const double* radii = basis->getCoord(SphericalBasis::COORD1);
    const double* thetas = basis->getCoord(SphericalBasis::COORD2);
    for (int iR = 0; iR < nRadial; iR++) {
      double r = radii[iR];
      for (int iTheta = 0; iTheta < nTheta; iTheta++) {
        double theta = thetas[iTheta];
        int iDerivative = basis->functionIndex(iR,iTheta);
        for (int iFunction = 0; iFunction < nRadial*nTheta; iFunction++) {
          int index = basis->matrixIndex(iDerivative,iFunction);
          int indexTrans = basis->matrixIndex(iFunction, iDerivative);
          diffOperator[indexTrans] = pow(r, regPower+2)*laplacian[index];
        }
      }
    }

    // Now we set up some boundary conditions.
    // Boundary condition. u' + u/r = ...
    const double* ddr = basis->getDiff(SphericalBasis::COORD1);
    int iR = nRadial-1;
    for (int iTheta = 0; iTheta < nTheta; iTheta++) {
      int iBC = basis->functionIndex(iR, iTheta);
      double sum = 0;
      for (int iFunction = 0; iFunction < nRadial*nTheta; iFunction++) {
        int iTranspose = basis->matrixIndex(iFunction, iBC);
        int iMatrix = basis->matrixIndex(iBC, iFunction);
        diffOperator[iTranspose] = ddr[iMatrix];
      }
      int iDiag = basis->matrixIndex(iBC, iBC);
      diffOperator[iDiag] += 1./radii[iR];
    }
  }
  return diffOperator;
}

double Hamiltonian::fillResidue(double* residue) {
  const double* laplacianBC = getOperator();
  int nRadial = basis->getRank(SphericalBasis::COORD1);
  int nTheta = basis->getRank(SphericalBasis::COORD2);
  int nPoints = nRadial*nTheta;
  const double* r = basis->getCoord(SphericalBasis::COORD1);
  double totalResidue = 0.;
  for (int iR = 0; iR < nRadial; iR++) {
    for (int iTheta = 0; iTheta < nTheta; iTheta++) {
      int iResidue = basis->functionIndex(iR, iTheta);
      residue[iResidue] = 0.0;
      for (int iu = 0; iu < nPoints; iu++) {
        int iLap = basis->matrixIndex(iu, iResidue);
        // Note minus sign to prepare for Jacobian inversion.
        residue[iResidue] -= laplacianBC[iLap]*u[iu];
      }
      // The bulk.
      if (iR < nRadial - 1) {
        // NOTE: this only has the delta function term in it.
        //Note minus sign.
        residue[iResidue] -= -0.125*betaOverSigma7[iTheta]
          + 0.125*kappa[iResidue]
          *pow(pow(r[iR], regPower)*(1. + 0.5*mass/r[iR] + u[iResidue])
              + sigma[iTheta], -7.);
      } else {
        residue[iResidue] -= -(regPower - 1.)*sigma[iTheta]
          *pow(r[iR], -(regPower + 1.));
      }
      totalResidue += pow(residue[iResidue], 2.);
    }
  }
  return sqrt(totalResidue);
}

void Hamiltonian::fillJacobian(double* jacobian) {
  const double* laplacianBC = getOperator();
  int nRadial = basis->getRank(SphericalBasis::COORD1);
  int nTheta = basis->getRank(SphericalBasis::COORD2);
  int nPoints = nRadial*nTheta;
  const double* r = basis->getCoord(SphericalBasis::COORD1);
  for (int iR = 0; iR < nRadial; iR++) {
    for (int iTheta = 0; iTheta < nTheta; iTheta++) {
      int iLap = basis->functionIndex(iR, iTheta);
      for (int iFunction = 0; iFunction < nPoints; iFunction++) {
        int iTranspose = basis->matrixIndex(iFunction, iLap);
        jacobian[iTranspose] = laplacianBC[iTranspose];
      }
      // NOTE: this only contains the rho = delta term.
      if (iR != (nRadial -1)) {
        int iDiag = basis->matrixIndex(iLap, iLap);
        jacobian[iDiag] -= 7./8.*kappa[iLap]*pow(
            pow(r[iR], regPower)*(1. + 0.5*mass/r[iR] + u[iLap])
               + sigma[iTheta], -8.)*pow(r[iR], regPower);
      }
    }
  }
}

void Hamiltonian::solve(double tolerance, int maxIterations) {
  // The fields need to be set up.
  setUpProblemData();
  // Then we perform Newton-Raphson iterations.
  bool converged = newtonRaphson(tolerance, maxIterations);
  if (!converged) {
    std::cerr << "Failed to converge to solution.\n";
  }
}

bool Hamiltonian::newtonRaphson(double tolerance, int maxIterations) {
  // Acquire spectral information.
  int nRadial = basis->getRank(SphericalBasis::COORD1);
  int nTheta = basis->getRank(SphericalBasis::COORD2);
  int nPoints = nRadial*nTheta;
  const double* r = basis->getCoord(SphericalBasis::COORD1);
  const double* theta = basis->getCoord(SphericalBasis::COORD2);
  // Set up the data.
  double jacobian[nPoints*nPoints];
  double residue[nPoints];
  double totalResidue = fillResidue(residue);
  int nIterations = 0;
  int one = 1, pivots[nPoints], info;
  while (totalResidue > tolerance && nIterations++ < maxIterations) {
    fillJacobian(jacobian);
    dgesv(&nPoints, &one, jacobian, &nPoints, pivots, residue, &nPoints, &info);
    if (info == 0) {
      for (int i = 0; i < nPoints; i++) u[i] += residue[i];
      totalResidue = fillResidue(residue);
      if (verbose) {
        std::cout << "Iteration " << nIterations <<
          " residue = " << totalResidue << "\n";
      }
    } else {
      // Inversion failed...
      std::cerr << "Inversion failed.\n";
      nIterations = maxIterations;
    }
  }

  if (totalResidue < tolerance) {
    return true;
  } else {
    if (severe) {
      delete[] u;
      u = NULL;
    }
    return false;
  }
}

double* Hamiltonian::getSingularAngularPart() {
  double* sigmaCopy = NULL;
  if (sigma != NULL) {
    int nTheta = basis->getRank(SphericalBasis::COORD2);
    sigmaCopy = new double[nTheta];
    for (int i = 0; i < nTheta; i++) sigmaCopy[i] = sigma[i];
  }
  return sigmaCopy;
}

double* Hamiltonian::getRemainder() {
  int nRadial = basis->getRank(SphericalBasis::COORD1);
  int nTheta = basis->getRank(SphericalBasis::COORD2);
  double* uCopy = new double[nRadial*nTheta];
  for (int i = 0; i < nRadial*nTheta; i++) {
    uCopy[i] = u[i];
  }
  return uCopy;
}

void Hamiltonian::setUpProblemData() {
  int nRadial = basis->getRank(SphericalBasis::COORD1);
  int nTheta = basis->getRank(SphericalBasis::COORD2);
  const double* r = basis->getCoord(SphericalBasis::COORD1);
  const double* theta = basis->getCoord(SphericalBasis::COORD2);
  delete[] u;
  u = new double[nRadial*nTheta];
  for (int i = 0; i < nRadial*nTheta; i++) u[i] = 0.0;
  // If the necessary singularity is less strong than the 1/r term from
  // the delta function, ignore it.
  SingularPart singularity(basis->getBasis(SphericalBasis::COORD2), k);
  delete[] sigma;
  sigma = new double[nTheta];
  regPower = singularity.getRegularityPower();
  delete[] betaOverSigma7;
  betaOverSigma7 = k.getBeta(theta, nTheta);
  if (mass > 0.0 && regPower < 1) {
    for (int i = 0; i < nTheta; i++) {
      sigma[i] = 0.0;
      betaOverSigma7[i] = 0.0;
    }
  } else {
    singularity.fillSigma(sigma);
    for (int i = 0; i < nTheta; i++) betaOverSigma7[i] *= pow(sigma[i], -7.);
  }
  // The extrinsic curvature terms.
  delete[] kappa;
  kappa = new double[nTheta*nRadial];
  for (int iTheta = 0; iTheta < nTheta; iTheta++) {
    for (int iR = 0; iR < nRadial; iR++) {
      int iKappa = basis->functionIndex(iR, iTheta);
      kappa[iKappa] = k.kappa(r[iR], theta[iTheta]);
    }
  }
}

void Hamiltonian::setBareMass(double m) {
  mass = m;
}

double Hamiltonian::getBareMass() {
  return mass;
}

void Hamiltonian::setSevere(bool severity) {
  severe = severity;
}

void Hamiltonian::setVerbose(bool verbosity) {
  verbose = verbosity;
}
