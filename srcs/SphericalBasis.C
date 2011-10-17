#include <cmath>
#include <iostream>
#include "SphericalBasis.h"
#include "ChebyshevRoots.h"
#include "cblas.h"

double* SphericalBasis::tensorInterpolate(const double* u, const double* r,
    int nR, const double* theta, int nTheta) {
  double coord1[nR];
  for (int iR = 0; iR < nR; iR++) {
    coord1[iR] = -1. + 2.0/rMax*r[iR];
  }
  double coord2[nTheta];
  for (int iTheta = 0; iTheta < nTheta; iTheta++) {
    coord2[iTheta] = -1. + 2.0/acos(-1.)*theta[iTheta];
  }
  return basis->tensorInterpolate(u, coord1, nR, coord2, nTheta);
}

double* SphericalBasis::interpolate1d(Direction dir, const double* u, 
    const double* x, int nX) {
  if (dir == COORD1) {
    double xi[nX];
    for (int i = 0; i < nX; i++) {
      xi[i] = 2./rMax*x[i] - 1.;
    }
    return radialBasis->interpolate(u, xi, nX);
  } else if (dir == COORD2) {
    double xi[nX];
    for (int i = 0; i < nX; i++) {
      xi[i] = M_2_PI*x[i] - 1.;
    }
    return thetaBasis->interpolate(u, xi, nX);
  }
  return NULL;
}

void SphericalBasis::fillCoefficients(const double* u, double* coeffs) {
  basis->fillCoefficients(u,coeffs);
}

SphericalBasis::SphericalBasis(): 
  laplacian(NULL), rMax(1.0), r(NULL), theta(NULL), dBydtheta(NULL), 
  dBydr(NULL)
{
  radialBasis = new ChebyshevRoots();
  thetaBasis = new ChebyshevRoots();
  basis = new Basis2D(radialBasis, thetaBasis);
}

SphericalBasis::~SphericalBasis() {
  delete basis;
  delete thetaBasis;
  delete radialBasis;
  delete[] laplacian;
  delete[] r;
  delete[] theta;
  delete[] dBydtheta;
  delete[] dBydr;
}

void SphericalBasis::setRanks(int radialRank, int zenithalRank) {
  int N1 = radialBasis->getRank();
  int N2 = thetaBasis->getRank();
  bool newRanks = !(N1 == radialRank && N2 == zenithalRank);
  if (radialRank > 1 && zenithalRank > 1 && newRanks) {
    delete[] laplacian; laplacian = NULL;
    basis->setRanks(radialRank, zenithalRank);
    delete[] dBydr; dBydr = NULL;
    delete[] dBydtheta; dBydtheta = NULL;
  }
}

void SphericalBasis::setMaximumRadius(double maximumRadius) {
  if (maximumRadius > 0.0 && maximumRadius != rMax) {
    delete[] r; r = NULL;
    delete[] laplacian; laplacian = NULL;
    delete[] dBydr; dBydr = NULL;
    delete[] dBydtheta; dBydtheta = NULL;
    rMax = maximumRadius;
  }
}

const double* SphericalBasis::getLaplacian() {
  if (laplacian == NULL) {
    int N1 = radialBasis->getRank();
    int N2 = thetaBasis->getRank();

    const double* dBydr = getdbydr();
    const double* dBydtheta = getdbydtheta();

    // For the second derivatives, we multiply the matrices...
    double d2r[N1*N1*N2*N2], d2theta[N1*N1*N2*N2];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        N1*N2, N1*N2, N1*N2, 1.0, dBydr, N1*N2, dBydr, N1*N2,
        0.0, d2r, N1*N2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        N1*N2, N1*N2, N1*N2, 1.0, dBydtheta, N1*N2,
        dBydtheta, N1*N2, 0.0, d2theta, N1*N2);

    // Now we combine the pieces into the laplacian.
    laplacian = new double[N1*N1*N2*N2];
    const double* x = radialBasis->getAbscissas();
    const double* u = thetaBasis->getAbscissas();
    for (int idR = 0; idR < N1; idR++) {
      double r = 0.5*rMax*(x[idR] + 1.);
      double invR = 1./r;
      for (int idTheta = 0; idTheta < N2; idTheta++) {
        int iDerivative = basis->functionIndex(idR, idTheta);
        double theta = (u[idTheta]+1.)*acos(-1.)*0.5;
        for (int iFunction = 0; iFunction < N1*N2; iFunction++) {
          int index = basis->matrixIndex(iDerivative, iFunction);
          laplacian[index] = d2r[index] + invR*invR*d2theta[index]
            + 2.0*invR*dBydr[index]
            + cos(theta)/sin(theta)*invR*invR*dBydtheta[index];
        }
      }
    }
  }
  return laplacian;
}

double* SphericalBasis::getdbydr() {
  if (dBydr == NULL) {
    const double* coordDiff = basis->getDifferentiationMatrices();
    int N1 = radialBasis->getRank();
    int N2 = thetaBasis->getRank();
    dBydr = new double[N1*N1*N2*N2];
    // Differentiation is just chain rule.
    for (int idTheta = 0; idTheta < N2; idTheta++) {
      for (int idR = 0; idR < N1; idR++) {
        int iDerivative = basis->functionIndex(idR, idTheta);
        // We multiply each row by the corresponding derivative
        // appearing in the chain rule.
        double dxdr = 2./rMax;
        for (int iFunction = 0; iFunction < N1*N2; iFunction++) {
          int index = basis->matrixIndex(iDerivative, iFunction);
          dBydr[index] = coordDiff[index]*dxdr;
        }
      }
    }
  }
  return dBydr;
}

double* SphericalBasis::getdbydtheta() {
  if (dBydtheta == NULL) {
    const double* coordDiff = basis->getDifferentiationMatrices();
    int N1 = radialBasis->getRank();
    int N2 = thetaBasis->getRank();
    dBydtheta = new double[N1*N1*N2*N2];
    const double* u = thetaBasis->getAbscissas();
    // Differentiation is just chain rule.
    for (int idR = 0; idR < N1; idR++) {
      for (int idTheta = 0; idTheta < N2; idTheta++) {
        int iDerivative = basis->functionIndex(idR, idTheta);
        // We multiply each row by the corresponding derivative
        // appearing in the chain rule.
        double dudtheta = 2./acos(-1.);
        for (int iFunction = 0; iFunction < N1*N2; iFunction++) {
          int index = basis->matrixIndex(iDerivative, iFunction);
          dBydtheta[index] = coordDiff[index + N1*N1*N2*N2]*dudtheta;
        }
      }
    }
  }
  return dBydtheta;
}

const double* SphericalBasis::getR() {
  if (r == NULL) {
    int N1 = radialBasis->getRank();
    r = new double[N1];
    const double* x = radialBasis->getAbscissas();
    for (int i = 0; i < N1; i++) {
      r[i] = 0.5*rMax*(x[i] + 1.);
    }
  }
  return r;
}

const double* SphericalBasis::getTheta() {
  if (theta == NULL) {
    int N2 = thetaBasis->getRank();
    theta = new double[N2];
    const double* u = thetaBasis->getAbscissas();
    for (int i = 0; i < N2; i++) {
      theta[i] = 0.5*acos(-1.)*(u[i] + 1.);
    }
  }
  return theta;
}

Basis* SphericalBasis::getBasis(Direction dir) {
  if (dir == COORD1) {
    return radialBasis;
  } else if (dir == COORD2) {
    return thetaBasis;
  }
  return NULL;
}

const double* SphericalBasis::getCoord(Direction dir) {
  if (dir == COORD1) {
    return getR();
  } else if (dir == COORD2) {
    return getTheta();
  }
  return NULL;
}

const double* SphericalBasis::getDiff(Direction dir) {
  if (dir == COORD1) {
    return getdbydr();
  } else if (dir == COORD2) {
    return getdbydtheta();
  }
  return NULL;
}

int SphericalBasis::getRank(Direction dir) {
  if (dir == COORD1) {
    return radialBasis->getRank();
  } else if (dir == COORD2) {
    return thetaBasis->getRank();
  }
  return 0;
}

int SphericalBasis::functionIndex(int i, int j) {
  return basis->functionIndex(i,j);
}

int SphericalBasis::matrixIndex(int iOutput, int iInput) {
  return basis->matrixIndex(iOutput, iInput);
}
