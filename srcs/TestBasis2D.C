#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <assert.h>
#define private public
#define protected public
#include "Basis2D.h"
#include "ChebyshevRoots.h"
#include "Legendre.h"

#define N1 5
#define N2 6

double abs(double x) {
  if (x < 0.) return -x;
  return x;
}

class TestBasis2D : public Basis2D {
  public:
    /**
     * Constructs the test object for the Basis2D class.
     */
    TestBasis2D(Basis* basisA, Basis* basisB): Basis2D(basisA,basisB) {};

    /**
     * Runs the tests on the Basis2D object.
     */
    void runTests();

  private:
    /**
     * Tests the construction of differentiation by computing the
     * derivative of a polynomial function and checking the result
     * against the exact result.
     */
    void testDifferentiate();

    /**
     * Tests the ability to determine 2D coefficients.
     */
    void testCoefficients();
};

void TestBasis2D::testDifferentiate() {
  // Create a random polynomial.
  double a[N1*N2];
  for (int i = 0; i < N1*N2; i++) 
    a[i] = (-1.0 + 2.0*(double)rand()/(double)RAND_MAX)/
      (double)(i*i+1);
  // Compute the values at the collocation points.
  const double* x = basis1->getAbscissas();
  const double* y = basis2->getAbscissas();
  double f[N1*N2];
  for (int iX = 0; iX < N1; iX++) {
    for (int iY = 0; iY < N2; iY++) {
      f[functionIndex(iX, iY)] = 0.0;
      for (int iPowX = 0; iPowX < N1; iPowX++) {
        for (int iPowY = 0; iPowY < N2; iPowY++) {
          f[functionIndex(iX,iY)] +=
            a[functionIndex(iPowX,iPowY)]*pow(x[iX], iPowX)*pow(y[iY], iPowY);
        }
      }
    }
  }
  // Compute the exact derivatives at the collocation points.
  double dfdx[N1*N2];
  double dfdy[N1*N2];
  for (int iX = 0; iX < N1; iX++) {
    for (int iY = 0; iY < N2; iY++) {
      dfdx[functionIndex(iX, iY)] = 0.0;
      for (int iPowX = 1; iPowX < N1; iPowX++) {
        for (int iPowY = 0; iPowY < N2; iPowY++) {
          dfdx[functionIndex(iX,iY)] +=
            iPowX*a[functionIndex(iPowX,iPowY)]
              *pow(x[iX], iPowX-1)*pow(y[iY], iPowY);
        }
      }
      dfdy[functionIndex(iX, iY)] = 0.0;
      for (int iPowX = 0; iPowX < N1; iPowX++) {
        for (int iPowY = 1; iPowY < N2; iPowY++) {
          dfdy[functionIndex(iX,iY)] +=
            iPowY*a[functionIndex(iPowX,iPowY)]
              *pow(x[iX], iPowX)*pow(y[iY], iPowY-1);
        }
      }
    }
  }
  // Now compute the approximate derivatives.
  double sdfdx[N1*N2], sdfdy[N1*N2];
  for (int i = 0; i < N1*N2; i++) {
    sdfdx[i] = 0.0;
    sdfdy[i] = 0.0;
  }
  const double* diff = getDifferentiationMatrices();
  for (int iDerivative = 0; iDerivative < N1*N2; iDerivative++) {
    for (int iFunction = 0; iFunction < N1*N2; iFunction++) {
      sdfdx[iDerivative] +=
        diff[matrixIndex(iDerivative,iFunction)]*f[iFunction];
      sdfdy[iDerivative] +=
        diff[matrixIndex(iDerivative,iFunction) + N1*N2*N1*N2]*f[iFunction];
    }
  }
  // Assert equality.
  for (int iX = 0; iX < N1; iX++) {
    for (int iY = 0; iY < N2; iY++) {
      int iFunction = functionIndex(iX,iY);
      assert( abs(dfdx[iFunction] - sdfdx[iFunction]) < 1.0e-14);
      assert( abs(dfdy[iFunction] - sdfdy[iFunction]) < 1.0e-14);
    }
  }
}

void TestBasis2D::testCoefficients() {
  double actualCoeffs[N1*N2];
  for (int iN1 = 0; iN1 < N1; iN1++) {
    for (int iN2 = 0; iN2 < N2; iN2++) {
      int i = functionIndex(iN1, iN2);
      actualCoeffs[i] = 1.0/(double)(iN1+1)/(double)(iN2+1);
    }
  }
  double u[N1*N2];
  const double* coord1 = basis1->getAbscissas();
  const double* coord2 = basis2->getAbscissas();
  for (int iCoord1 = 0; iCoord1 < N1; iCoord1++) {
    double x = coord1[iCoord1];
    for (int iCoord2 = 0; iCoord2 < N2; iCoord2++) {
      double y = coord2[iCoord2];
      int index = functionIndex(iCoord1, iCoord2);
      u[index] = 0.;
      for (int iN1 = 0; iN1 < N1; iN1++) {
        for (int iN2 = 0; iN2 < N2; iN2++) {
          int coeffIndex = functionIndex(iN1, iN2);
          u[index] += actualCoeffs[coeffIndex]*basis1->function(iN1, x)*
            basis2->function(iN2, y);
        }
      }
    }
  }
  // Spectrally compute the coefficients;
  double newCoeffs[N1*N2];
  fillCoefficients(u, newCoeffs);
  for (int iN1 = 0; iN1 < N1; iN1++) {
    for (int iN2 = 0; iN2 < N2; iN2++) {
      int index = functionIndex(iN1, iN2);
      assert(abs(actualCoeffs[index] - newCoeffs[index]) < 1.0e-15);
    }
  }
}

void TestBasis2D::runTests() {
  int nTests = 0;
  std::cout << "Running tests on Basis2D class.\n";

  testDifferentiate();
  nTests++; std::cout << ".\n";

  testCoefficients();
  nTests++; std::cout << ".\n";

  std::cout << "Ok. Ran " << nTests << " tests.\n";
}

int main() {
  // Create the 2d basis.
  ChebyshevRoots* basis1 = new ChebyshevRoots();
  basis1->setRank(N1);
  Legendre* basis2 = new Legendre();
  basis2->setRank(N2);
  TestBasis2D test(basis1, basis2);

  // Initialize the random number generator.
  srand(time(NULL));
  test.runTests();
  return 0;
}
