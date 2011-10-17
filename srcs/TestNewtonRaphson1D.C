#include <iostream>
#include <cassert>
#include <ctime>
#include <cmath>

#define private public
#include "NewtonRaphson1D.h"
#include "ChebyshevExtrema.h"

#define RANK 31
#define TOLERANCE 1.0e-16
#define MAXITERATIONS 1

double abs(double x) {
  if (x < 0.0) {
    return -x;
  }
  return x;
}

class TestNewtonRaphson1D: public NewtonRaphson1D {
  public:
    TestNewtonRaphson1D(Basis* aBasis);
    void runTests();
    void runJacobianTest();
    void runResidueTest();
  protected:
    virtual double f(double x, double u, double up, double upp);
    virtual double dfdu(double x, double u, double up, double upp);
    virtual double dfdup(double x, double u, double up, double upp);
    virtual double dfdupp(double x, double u, double up, double upp);
    virtual double g(double u, double up);
    virtual double dgdu(double u, double up);
    virtual double h(double u, double up);
    virtual double dhdu(double u, double up);

  private:
    void runBCTest();
    double a;
    double b;
    double c;
    double g0;
    double h0;
    double alphabc;
    double betabc;
    double alphaiv;
    double betaiv;
    double rplus;
    double rminus;
};

TestNewtonRaphson1D::TestNewtonRaphson1D(Basis* aBasis) : 
    NewtonRaphson1D(aBasis) {
  srand(time(NULL));
  a = -0.5 + 0.5*(double)rand()/(double)RAND_MAX;
  b = 2.0 + (double)rand()/(double)RAND_MAX;
  c = -0.5 + 0.5*(double)rand()/(double)RAND_MAX;
  g0 = -0.5 + 0.5*(double)rand()/(double)RAND_MAX;
  h0 = -0.5 + 0.5*(double)rand()/(double)RAND_MAX;
  double discriminant = sqrt(b*b - 4.*a*c);
  rplus  = 0.5*(-b + discriminant)/a;
  rminus = 0.5*(-b - discriminant)/a;

  alphabc = - (-h0*exp(-rplus) + g0*exp(rplus))/
    (-exp(-rminus+rplus) + exp(rminus - rplus));
  betabc  =   (-h0*exp(-rminus) + g0*exp(rminus))/
    (-exp(-rminus+rplus) + exp(rminus - rplus));

  alphaiv = -(g0*rplus  - h0)*exp(rminus)/(rminus - rplus);
  betaiv  =  (g0*rminus - h0)*exp(rplus )/(rminus - rplus);
}

void TestNewtonRaphson1D::runResidueTest() {
  registerBoundaries(BOUNDARY, BOUNDARY);
  double u[RANK];
  for (int i = 0; i < RANK; i++) u[i] = 0;
  setPointerToSolution(u);
  double up[RANK], upp[RANK];
  sigmap = up;
  sigmapp = upp;
  computeDifferentiationMatrices();
  computeDerivatives();
  residue = new double[RANK];
  double totalResidue = getTotalResidue();

  assert(abs(residue[0] - g0) < 1.0e-15);
  for (int i = 1; i < RANK-1; i++) {
    assert(abs(residue[i]) < 1.0e-15);
  }
  assert(abs(residue[RANK-1] - h0) < 1.0e-15);

  delete[] residue;
  sigma = NULL;
  sigmap = NULL;
  sigmapp = NULL;
  delete[] diff;
  delete[] doubleDiff;
}

void TestNewtonRaphson1D::runJacobianTest() {
  registerBoundaries(BOUNDARY, BOUNDARY);
  double u[RANK];
  for (int i = 0; i < RANK; i++) u[i] = 0;
  setPointerToSolution(u);
  double up[RANK], upp[RANK];
  sigmap = up;
  sigmapp = upp;
  computeDifferentiationMatrices();
  computeDerivatives();
  residue = new double[RANK];
  double* jac = getJacobian();

  double fdjac[RANK*RANK];
  
  double h = 1.0e-5;
  double residueplus[RANK];
  for (int iCollocation = 0; iCollocation < RANK; iCollocation++) {
    u[iCollocation] = h;
    computeDerivatives();
    getTotalResidue();
    for (int i = 0; i < RANK; i++) residueplus[i] = residue[i];

    u[iCollocation] = -h;
    computeDerivatives();
    getTotalResidue();
    for (int i = 0; i < RANK; i++) {
      int ijac = i*RANK + iCollocation;
      fdjac[ijac] = (residueplus[i] - residue[i])/(2*h);
    }
    
    u[iCollocation] = 0;
  }

  for (int i = 0; i < RANK; i++) {
    for (int j = 0; j < RANK; j++) {
      int ijac = i*RANK + j;
      assert(fdjac[ijac]+jac[ijac] < 2.0e-12);
    }
  }

  delete[] residue;
  sigma = NULL;
  sigmap = NULL;
  sigmapp = NULL;
  delete[] diff;
  delete[] doubleDiff;
  delete[] jac;
}

void TestNewtonRaphson1D::runTests() {
  std::cout << "Testing Newton-Raphson 1D class.\n";
  int nTests = 0;

  runResidueTest();
  nTests++; std::cout << ".\n";

  runJacobianTest();
  nTests++; std::cout << ".\n";

  runBCTest();
  nTests++; std::cout << ".\n";

  std::cout << "Done. Ran " << nTests << " tests successfully.\n";
}

void TestNewtonRaphson1D::runBCTest() {
  registerBoundaries(BOUNDARY, BOUNDARY);
  double u[RANK];
  const double* x = basis->getAbscissas();
  for (int i = 0; i < RANK; i++) {
    u[i] = 0;
  }
  computeSolution(u, NULL, TOLERANCE, MAXITERATIONS);
  for (int i = 0; i < RANK; i++) {
    double solution = alphabc*exp(rminus*x[i]) + betabc*exp(rplus*x[i]);
    if (!(abs( u[i] - solution) < 1.0e-10)) {
      std::cerr << "Bad luck. Sometimes the random number make the";
      std::cerr << " exponential solution converge very poorly. Try again.\n";
      bool solutionWithinTolerance = (abs( u[i] - solution) < 1.0e-10);
      std::cout << u[i] - solution << "\n";
      assert(solutionWithinTolerance);
    }
  }
}

double TestNewtonRaphson1D::f(double x, double u, 
    double up, double upp) {
  return a*upp + b*up + c*u;
}


double TestNewtonRaphson1D::dfdu(double x, double u, 
    double up, double upp) {
  return c;
}

double TestNewtonRaphson1D::dfdup(double x, double u, 
    double up, double upp) {
  return b;
}

double TestNewtonRaphson1D::dfdupp(double x, double u, 
    double up, double upp) {
  return a;
}

double TestNewtonRaphson1D::g(double u, double up) {
  return u - g0;
}

double TestNewtonRaphson1D::dgdu(double u, double up) {
  return 1.;
}

double TestNewtonRaphson1D::h(double u, double up) {
  return u - h0;
}

double TestNewtonRaphson1D::dhdu(double u, double up) {
  return 1.;
}

int main() {
  int rank = RANK;
  ChebyshevExtrema* basis = new ChebyshevExtrema();
  basis->setRank(RANK);
  TestNewtonRaphson1D test(basis);
  test.runTests();
  delete basis;
}
