#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <time.h>
#include "ChebyshevRoots.h"

/**
 * A class to test the ChebyshevRoots class.
 */
class TestChebyshevRoots: public ChebyshevRoots {
  public:
    /**
     * The constructor does nothing.
     */
    TestChebyshevRoots() {};

    /**
     * Runs the various tests on the ChebyshevRoots class.
     */
    void runTests();

  private:

    /**
     * Tests the integration ability by integrating some polynomials.
     */
    void testIntegrate();

    /**
     * Tests the basis function. Calls the first 5 of them and tests the
     * output versus the known output. Asserts that all the basis
     * functions are computed properly.
     */
    void testBasisFunction();

    /**
     * Tests the interpolation functionality.
     */
    void testInterpolate();

    /**
     * Tests the ability to evaluate a function from it's coefficients
     * using Clenshaw. Picks random coefficients and 10 random points
     * and compares to the basisFunction evaluation method (which should
     * pass the tests independently).
     */
    void testEvaluate();

    /**
     * Tests the coefficientsToValues matrix. Uses a set of coefficients
     * and computes the corresponding values. Tests against known
     * values.
     */
    void testCoefficientsToValuesMatrices();

    /**
     * Gets the coefficientsOfDerivativeMatrix from the super class and
     * then tests it. Asserts that the tests are passed to a given
     * internal tolerance
     */
    void testDifferentiateCoefficients();

    /**
     * Gets the differentiation matrix and then differentiates a
     * function. It asserts that the values of the derivatives are equal
     * to the known values to an internal tolerance.
     */
    void testDifferentiate();

    /**
     * Tests the ability to regain the coefficients from the values.
     */
    void testValuesToCoefficients();

    /**
     * Tests the ability to transform coefficients;
     */
    void testChebyshevToTaylor();
};

double abs(double x) {
  if (x < 0.0) return -x;
  return x;
}

void TestChebyshevRoots::runTests() {
  printf("Running tests on ChebyshevRoots polynomial class.\n");
  int tests = 0;
  assert(setRank(5) == true);

  testBasisFunction();
  printf(".\n"); tests++; 

  testEvaluate();
  printf(".\n"); tests++; 

  testInterpolate();
  printf(".\n"); tests++; 

  testCoefficientsToValuesMatrices();
  printf(".\n"); tests++; 

  testDifferentiateCoefficients();
  printf(".\n"); tests++; 

  testDifferentiate();
  printf(".\n"); tests++; 

  testValuesToCoefficients();
  printf(".\n"); tests++; 

  testChebyshevToTaylor();
  printf(".\n"); tests++; 

  testIntegrate();
  printf(".\n"); tests++; 

  printf("OK. Ran %d tests successfully.\n",tests);
}

void TestChebyshevRoots::testIntegrate() {
  int rankHere = 55;
  setRank(rankHere);
  double f[rankHere];
  double a[4];
  a[0] = 1.0;
  a[1] = 2.0;
  a[2] = 3.0;
  a[3] = 4.0; // f = 1 + 2x + 3x^2 + 4x^3
  
  const double* x = getAbscissas();
  for (int i = 0; i < rankHere; i++) {
    f[i] = a[0];
    for (int j = 1; j < 5; j++) f[i] += a[j]*pow(x[i], j);
  }
  double integral = integrate(f);
  assert(abs(integral - 2.5*M_PI) < 1.0e-15);
  setRank(5);
}

void TestChebyshevRoots::testChebyshevToTaylor() {
  double f[5];
  const double* x = getAbscissas();
  double* c2t = chebyshevToTaylorCoefficientsMatrix();
  for (int ipow = 0; ipow < 5; ipow++) {
    for (int i = 0; i < 5; i++) f[i] = pow(x[i], ipow);
    double coeffs[5];
    fillCoefficients(f, coeffs);
    for(int iTaylor = 0; iTaylor < 5; iTaylor++) {
      double taylorCoeff = 0.0;
      for (int iCheb = 0; iCheb < 5; iCheb++) {
        taylorCoeff += c2t[index(iTaylor, iCheb)]*coeffs[iCheb];
      }
      if (iTaylor == ipow) {
        assert(abs(taylorCoeff-1.0) < 1.0e-15);
      } else {
        assert(abs(taylorCoeff) < 1.0e-15);
      }
    }
  }
  delete[] c2t;
}

void TestChebyshevRoots::testBasisFunction() {
  double x = 0.375;
  assert(function(0,x) == 1);
  assert(function(1,x) == x);
  assert(function(2,x) == -.71875);
  assert(function(3,x) == -.9140625);
  assert(function(4,x) == 0.033203125);
  assert(abs(function(5,x) - 0.93896484375) < 2.0E-16);
}

void TestChebyshevRoots::testCoefficientsToValuesMatrices() {
  double coefficients[5] = {1.0, 0.5, 0.25, 0.125, 0.0625};
  double values[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  const double* matrix = coefficientsToValuesMatrix();
  for (int iRow = 0; iRow < 5; iRow++) {
    for (int iCol = 0; iCol < 5; iCol++) {
      int i = index(iRow, iCol);
      values[iRow] += coefficients[iCol]*matrix[i];
    }
  }
  double solution[5] = {1.7705692254263069972, 1.0471927508671712981,
    .81250000000000000000, .69717162764848656190, .67256639605803514283};
  for (int i = 0; i < 5; i++) {
    assert( abs(values[i] - solution[4-i]) < 2.5e-16);
  }
  delete[] matrix;
  matrix = getValuesToCoefficientsMatrix();
  double newCoefficients[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  for (int iRow = 0; iRow < 5; iRow++) {
    for (int iCol = 0; iCol < 5; iCol++) {
      int i = index(iRow, iCol);
      newCoefficients[iRow] += matrix[i]*values[iCol];
    }
  }
  for (int i = 0; i < 5; i++) 
    assert( abs(newCoefficients[i] - coefficients[i]) < 5.0e-16);
}

void TestChebyshevRoots::testDifferentiateCoefficients() {
  double* matrix = coefficientsOfDerivativeMatrix();
  double coefficients[5] = {1., 2., 3., 4., 5.};
  double derivativeCoeffs[5] = {0., 0., 0., 0., 0.};
  for (int iRow = 0; iRow < 5; iRow++) {
    for (int iCol = 0; iCol < 5; iCol++) {
      int i = index(iRow, iCol);
      derivativeCoeffs[iRow] += matrix[i]*coefficients[iCol];
    }
  }
  double actualCoeffs[5] = {14.0, 52.0, 24.0, 40.0, 0.};
  for (int i = 0; i < 5; i++)
    assert(derivativeCoeffs[i] == actualCoeffs[i]);
}

void TestChebyshevRoots::testDifferentiate() {
  const double* differentiate = getDifferentiationMatrix();
  // These are the values of the function used in testCoefficientsToValues()
  double values[5] = {1.7705692254263069972, 1.0471927508671712981,
    .81250000000000000000, .69717162764848656190, .67256639605803514283};
  double derivative[5] = {0., 0., 0., 0., 0.};
  for (int iRow = 0; iRow < 5; iRow++) {
    for (int iCol = 0; iCol < 5; iCol++) {
      derivative[iRow] += differentiate[index(iRow, iCol)]*values[4-iCol];
    }
  }
  double solution[5] = {3.2022401463701774909, 1.0493868745099223395,
    .12500000000000000001, .23708763392765652423, -.23871465480775635474};
  for (int i = 0; i < 5; i++) 
    assert(abs(derivative[i]-solution[4-i]) < 4.E-15);
}

void TestChebyshevRoots::testValuesToCoefficients() {
  // These were used in the test of the inverse functionality.
  double values[5] = {1.7705692254263069972, 1.0471927508671712981,
    .81250000000000000000, .69717162764848656190, .67256639605803514283};
  double revValues[5] = {values[4], values[3], values[2], values[1], values[0]};
  double solution[5] = {1.0, 0.5, 0.25, 0.125, 0.0625};
  double coeffs[5];
  fillCoefficients(revValues, coeffs);
  for (int i = 0; i < 5; i++)
    assert(abs(coeffs[i]-solution[i]) < 4.5e-16);
}

void TestChebyshevRoots::testEvaluate() {
  srand(time(NULL));
  double coeffs[5];
  for (int i = 0; i < 5; i++) 
    coeffs[i] = -5.0 + 10.0*(double)rand()/(double)RAND_MAX;
  // Now we actually test some points.
  for (int i = 0; i < 10; i++) {
    double x = -1.0 + 2.0*(double)rand()/(double)RAND_MAX;
    double f = 0.0;
    for (int n = 0; n < 5; n++) f += coeffs[n]*function(n,x);
    double f2 = evaluate(x,coeffs);
    assert(abs(f-f2) < 1.0e-14); // I want to be conservative on the assert.
  }
}

void TestChebyshevRoots::testInterpolate() {
  double x[25];
  for (int i = 0; i < 25; i++) 
    x[i] = -1.0 + 2.0*(double)rand()/(double)RAND_MAX;
  double values[5];
  for (int i = 0; i < 5; i++)
    values[i] = -5.0 + 10.0*(double)rand()/(double)RAND_MAX;
  double coeffs[5];
  fillCoefficients(values, coeffs);
  double* fAtX = interpolate(values,x,25);
  for (int i = 0; i < 25; i++) {
    double f2 = evaluate(x[i],coeffs);
    assert(abs(fAtX[i]-f2) < 1.0e-14);
  }
  delete[] fAtX;
}

int main() {
  TestChebyshevRoots test;
  test.runTests();
  return 0;
}
