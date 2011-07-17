#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cassert>

#include "ChebyshevExtrema.h"

/**
 * A class to test the ChebyshevExtrema class.
 */
class TestChebyshevExtrema: public ChebyshevExtrema {
  public:
    /**
     * The constructor does nothing.
     */
    TestChebyshevExtrema() {};

    /**
     * Runs the various tests on the ChebyshevExtrema class.
     */
    void runTests();

  private:

    /**
     * Tests the quadrature formulae, mostly the weight formula.
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
};

double abs(double x) {
  if (x < 0.0) return -x;
  return x;
}

void TestChebyshevExtrema::runTests() {
  printf("Running tests on ChebyshevExtrema polynomial class.\n");
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

  testIntegrate();
  printf(".\n"); tests++; 

  printf("OK. Ran %d tests successfully.\n",tests);
}

void TestChebyshevExtrema::testIntegrate() {
  int rankHere = 5;
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
}

void TestChebyshevExtrema::testBasisFunction() {
  double x = 0.375;
  assert(function(0,x) == 1);
  assert(function(1,x) == x);
  assert(function(2,x) == -.71875);
  assert(function(3,x) == -.9140625);
  assert(function(4,x) == 0.033203125);
  assert(abs(function(5,x) - 0.93896484375) < 2.0E-16);
}

void TestChebyshevExtrema::testCoefficientsToValuesMatrices() {
  double coefficients[5] = {1.0, 0.5, 0.25, 0.125, 0.0625};
  double values[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  const double* matrix = coefficientsToValuesMatrix();
  for (int iRow = 0; iRow < 5; iRow++) {
    for (int iCol = 0; iCol < 5; iCol++) {
      int i = index(iRow, iCol);
      values[iRow] += coefficients[iCol]*matrix[i];
    }
  }
  double solution[5] = {1.9375, 1.2026650429449553216,
    .8125, .67233495705504467835, .6875};
  for (int i = 0; i < 5; i++) 
    assert( abs(values[i] - solution[4-i]) < 2.5e-16);
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

void TestChebyshevExtrema::testDifferentiateCoefficients() {
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

void TestChebyshevExtrema::testDifferentiate() {
  const double* differentiate = getDifferentiationMatrix();
  // These are the values of the function used in testCoefficientsToValues()
  // The ordering is backwards though, so the derivative has a minus
  // sign.
  double values[5] = {1.9375000000000000000, 1.2026650429449553216,
    .81250000000000000000, .67233495705504467835, .68750000000000000000};
  double derivative[5] = {0., 0., 0., 0., 0.};
  for (int iRow = 0; iRow < 5; iRow++) {
    for (int iCol = 0; iCol < 5; iCol++) {
      derivative[iRow] += differentiate[index(iRow, iCol)]*values[iCol];
    }
  }
  double solution[5] = {3.6250000000000000000, 1.5821067811865475244,
    .12500000000000000001, .16789321881345247559, -.37500000000000000003};
  for (int i = 0; i < 5; i++) 
    assert(abs(derivative[i]+solution[i]) < 4.E-15);
}

void TestChebyshevExtrema::testValuesToCoefficients() {
  // These were used in the test of the inverse functionality.
  // The change in signs to the abscissas picks up a minus sign on the
  // odd coefficients.
  double values[5] = {1.9375000000000000000, 1.2026650429449553216,
    .81250000000000000000, .67233495705504467835, .68750000000000000000};
  double solution[5] = {1.0, -0.5, 0.25, -0.125, 0.0625};
  double coeffs[5];
  fillCoefficients(values, coeffs);
  for (int i = 0; i < 5; i++)
    assert(abs(coeffs[i]-solution[i]) < 4.0e-16);
}

void TestChebyshevExtrema::testEvaluate() {
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

void TestChebyshevExtrema::testInterpolate() {
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
  TestChebyshevExtrema test;
  test.runTests();
  return 0;
}
