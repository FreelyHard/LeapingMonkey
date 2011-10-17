#include <cmath>
#include <iostream>

#include "cblas.h"
#include "lapacke.h"
#include "SingularPart.h"

bool SingularPart::newtonRaphson(double* sigma, double* beta,
    double* diffOp) const {
  int nTheta = basis->getRank();
  double* residue = new double[nTheta];
  double* jacobian = new double[nTheta*nTheta];
  double totalResidue = fillResidue(sigma, beta, diffOp, residue);
  int maxIterations = 10, nIterations = 0;
  double tolerance = 1.0e-12*nTheta;
  double p = getRegularityPower();
  while (totalResidue > tolerance && nIterations++ < maxIterations) {
    // Make the Jacobian
    for (int i = 0; i < nTheta; i++) {
      for (int j = 0; j < nTheta; j++) {
        int iTrans = basis->index(j,i);
        int iDirect = basis->index(i,j);
        // We transpose because of usage of Lapack (fortran).
        jacobian[iTrans] = diffOp[iDirect];
      }
      int iDelta = basis->index(i,i);
      jacobian[iDelta] += p*(p-1.) - 0.875*beta[i]*pow(sigma[i], -8);
    }
    // Solve the linear system.
    int one = 1, pivots[nTheta];
#warning "TODO: Don't need to compute transpose."
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, nTheta, one, jacobian, 
        nTheta, pivots, residue, nTheta);
    for (int i = 0; i < nTheta; i++) {
      sigma[i] += residue[i];
    }
    totalResidue = fillResidue(sigma, beta, diffOp, residue);
  }
  return (maxIterations >= nIterations) && (totalResidue < tolerance);
}

bool SingularPart::fillSigma(double* sigma) const {
  // Set up the variables.
  int nTheta = basis->getRank();
  const double* u = basis->getAbscissas();
  double* theta = new double[nTheta];
  for (int i = 0; i < nTheta; i++) theta[i] = acos(-1.)*0.5*(u[i] + 1.);
  double* beta =  curvature.getBeta(theta, nTheta);
  double p = getRegularityPower();

  // Set sigma to the solution if beta were constant.
  for (int i = 0; i < nTheta; i++) sigma[i] = p < 1. ?
    pow(beta[i]*0.125/(p*(1.-p)), 0.125) :
    pow(beta[i]*0.125/(p*(p - 1.)), 0.125) ;

  double* lap = getOperator(theta);
  bool converged = newtonRaphson(sigma, beta, lap);
  if (!converged) {
    std::clog << "Newton-Raphson did not converge.\n";
  }

  delete[] lap;
  delete[] theta;
  delete[] beta;
  return converged;
}

double* SingularPart::getOperator(double* theta) const  {
  int nTheta = basis->getRank();
  const double* ddu = basis->getDifferentiationMatrix(); // actually d/du
  double dudtheta = 2./acos(-1.);
  double* laplacian = new double[nTheta*nTheta];
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
      nTheta, nTheta, nTheta, 1.0, ddu, nTheta, ddu, nTheta,
      0.0, laplacian, nTheta);
  for (int iLaplacian = 0; iLaplacian < nTheta; iLaplacian++) {
    double cotTheta = cos(theta[iLaplacian])/sin(theta[iLaplacian]);
    for (int iFunction = 0; iFunction < nTheta; iFunction++) {
      int i = basis->index(iLaplacian, iFunction);
      laplacian[i] = (laplacian[i]*dudtheta + cotTheta*ddu[i])*dudtheta;
    }
  }
  return laplacian;
}

double SingularPart::fillResidue(double* sigma, double* beta, 
    double* differentialOperator, double* residue) const {
  int nTheta = basis->getRank();
  double p = getRegularityPower();
  double totalResidue = 0.0;
  for (int iLap = 0; iLap < nTheta; iLap++) {
    double lap = 0;
    for (int iSigma = 0; iSigma < nTheta; iSigma++) {
      int i = basis->index(iLap, iSigma);
      lap += differentialOperator[i]*sigma[iSigma];
    }
    residue[iLap] = -(lap + p*(p-1.)*sigma[iLap] + 
      beta[iLap]*0.125*pow(sigma[iLap], -7.));
    totalResidue += pow(residue[iLap], 2.);
  }
  return sqrt(totalResidue);
}

double SingularPart::getRegularityPower() const {
  int n = curvature.getSingularityPower();
  return (double)(n-2)*0.125;
}
