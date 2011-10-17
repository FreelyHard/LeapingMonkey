#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cassert>
#include "cblas.h"
// No, I'm not proud of myself.
#define private protected
#include "SphericalBasis.h"

#define N1 15
#define N2 17

double abs(double x) {
  if (x < 0.) return -x;
  return x;
}

class TestSphericalBasis: public SphericalBasis {
  public:
    TestSphericalBasis(): SphericalBasis() {};
    void runTests();
  private:
    void runDifferentiateTest();
    void runLaplacianTest();
    void runRGradTest();
};

void TestSphericalBasis::runDifferentiateTest() {
  double f[N1*N2];        // f = 3x^2 + 4*z^2 or r^2(3 + cos^2\theta)
  double dfdr[N1*N2];     // dfdr = 2r(3 + cos^2\theta)
  double dfdtheta[N1*N2]; // dfdtheta = -2r^2*cos*sin
  const double* r = getCoord(COORD1);
  const double* theta = getCoord(COORD2);
  for (int iR = 0; iR < N1; iR++) {
    for (int iTheta = 0; iTheta < N2; iTheta++) {
      int index = basis->functionIndex(iR, iTheta);
      f[index] = r[iR]*r[iR]*(3. + cos(theta[iTheta])*cos(theta[iTheta]));
      dfdr[index] = 2.*f[index]/r[iR];
      dfdtheta[index] = -2*r[iR]*r[iR]*cos(theta[iTheta])*sin(theta[iTheta]);
    }
  }
  // Compute the derivatives spectrally.
  double* ddr = getdbydr();
  double* ddtheta = getdbydtheta();
  double sdfdr[N1*N2], sdfdtheta[N1*N2];
  for (int iDerivative = 0; iDerivative < N1*N2; iDerivative++) {
    sdfdr[iDerivative] = 0.;
    sdfdtheta[iDerivative] = 0.;
    for (int iFunction = 0; iFunction < N1*N2; iFunction++) {
      int index = basis->matrixIndex(iDerivative, iFunction);
      sdfdr[iDerivative] += ddr[index]*f[iFunction];
      sdfdtheta[iDerivative] += ddtheta[index]*f[iFunction];
    }
    assert(abs(sdfdr[iDerivative] - dfdr[iDerivative]) < 9.0e-9);
    assert(abs(sdfdtheta[iDerivative] - dfdtheta[iDerivative]) < 9.0e-9);
  }

  double d2r[N1*N1*N2*N2];
  double d2theta[N1*N1*N2*N2];
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
       N1*N2, N1*N2, N1*N2, 1.0, ddtheta, N1*N2,
       ddtheta, N1*N2, 0.0, d2theta, N1*N2);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
       N1*N2, N1*N2, N1*N2, 1.0, ddr, N1*N2,
       ddr, N1*N2, 0.0, d2r, N1*N2);
  for (int iTheta = 0; iTheta < N2; iTheta++) {
    for (int iR = 0; iR < N1; iR++) {
      int iDerivative = basis->functionIndex(iR, iTheta);
      double sd2ft = 0.0;
      double sd2fr = 0.0;
      for (int iFunction = 0; iFunction < N1*N2; iFunction++) {
        int index = basis->matrixIndex(iDerivative, iFunction);
        sd2ft += d2theta[index]*f[iFunction];
        sd2fr += d2r[index]*f[iFunction];
      }
      double sfakelap = 2.*sd2fr + sd2ft/r[iR]/r[iR];
      assert(abs(sfakelap - 14.) < 1.0e-8);
    }
  }

}

void TestSphericalBasis::runLaplacianTest() {
  double f[N1*N2];        // f = 3x^2 + 4*z^2 or r^2(3 + cos^2\theta)
  const double* x = radialBasis->getAbscissas();
  const double* u = thetaBasis->getAbscissas();
  for (int iR = 0; iR < N1; iR++) {
    for (int iTheta = 0; iTheta < N2; iTheta++) {
      double r = (x[iR] + 1.)*0.5*rMax;
      int index = basis->functionIndex(iR, iTheta);
      double theta = 0.5*acos(-1.)*(u[iTheta] + 1.);
      f[index] = r*r*(3. + cos(theta)*cos(theta));
    }
  }
  const double* lap = getLaplacian();
  double lapF[N1*N2];
  for (int iTheta = 0; iTheta < N2; iTheta++) {
    for (int iR = 0; iR < N1-1; iR++) {
      int iDerivative = basis->functionIndex(iR, iTheta);
      lapF[iDerivative] = 0.;
      for (int iFunction = 0; iFunction < N1*N2; iFunction++) {
        int index = basis->matrixIndex(iDerivative, iFunction);
        lapF[iDerivative] += lap[index]*f[iFunction];
      }
      assert(abs(lapF[iDerivative] - 20.) < 1.0e-8);
    }
  }
}

void TestSphericalBasis::runTests() {
  int ntests = 0;
  std::cout << "Testing SphericalBasis class.\n";
  
  runDifferentiateTest();
  ntests++; std::cout << ".\n";

  runLaplacianTest();
  ntests++; std::cout << ".\n";

  std::cout << "OK. Ran " << ntests << " tests.\n";
}

int main() {
  TestSphericalBasis test;
  test.setRanks(N1,N2);
  test.runTests();
}
