#include "TestLegendre.HH"
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <assert.h>

double abs(double x) {
  if (x < 0.0) return -x;
  return x;
}

void TestLegendre::runTests() {
  std::cout << "Running tests on Legendre polynomial class.\n";
  int tests = 0;
  assert(setRank(5) == true);

  testBasisFunction();
  std::cout << ".\n"; tests++; 

  testEvaluate();
  std::cout << ".\n"; tests++; 

  testInterpolate();
  std::cout << ".\n"; tests++; 

  testCoefficientsToValuesMatrices();
  std::cout << ".\n"; tests++; 

  testDifferentiateCoefficients();
  std::cout << ".\n"; tests++; 

  testDifferentiate();
  std::cout << ".\n"; tests++; 

  testValuesToCoefficients();
  std::cout << ".\n"; tests++; 

  std::cout << "OK. Ran " << tests << " tests successfully.\n";
}

void TestLegendre::testBasisFunction() {
  double x = 0.5;
  assert(function(0,x) == 1);
  assert(function(1,x) == x);
  assert(function(2,x) == -0.125);
  assert(function(3,x) == -0.4375);
  assert(function(4,x) == -0.2890625);
  assert(function(5,x) == 0.08984375);
}

void TestLegendre::testCoefficientsToValuesMatrices() {
  double coefficients[5] = {1.0, 0.5, 0.25, 0.125, 0.0625};
  double values[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  const double* matrix = coefficientsToValuesMatrix();
  for (int iRow = 0; iRow < 5; iRow++) {
    for (int iCol = 0; iCol < 5; iCol++) {
      int i = index(iRow, iCol);
      values[iRow] += coefficients[iCol]*matrix[i];
    }
  }
  double solution[5] = {.68257536428805005463, .74513775143552098982,
    .89843750000000000000, 1.1792615356095370567, 1.7140130029878795533};
  for (int i = 0; i < 5; i++) 
    assert( abs(values[i] - solution[i]) < 2.5e-16);
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

void TestLegendre::testDifferentiateCoefficients() {
  double* matrix = coefficientsOfDerivativeMatrix();
  double coefficients[5] = {1., 2., 3., 4., 5.};
  double derivativeCoeffs[5] = {0., 0., 0., 0., 0.};
  for (int iRow = 0; iRow < 5; iRow++) {
    for (int iCol = 0; iCol < 5; iCol++) {
      int i = index(iRow, iCol);
      derivativeCoeffs[iRow] += matrix[i]*coefficients[iCol];
    }
  }
  double actualCoeffs[5] = {6., 24., 20., 35., 0.};
  for (int i = 0; i < 5; i++) 
    assert(derivativeCoeffs[i] == actualCoeffs[i]);
}

void TestLegendre::testDifferentiate() {
  const double* differentiate = getDifferentiationMatrix();
  // These are the values of the function used in testCoefficientsToValues()
  double values[5] = {.68257536428805005463, .74513775143552098982,
    .89843750000000000000, 1.1792615356095370567, 1.7140130029878795533};
  double derivative[5] = {0., 0., 0., 0., 0.};
  for (int iRow = 0; iRow < 5; iRow++) {
    for (int iCol = 0; iCol < 5; iCol++) {
      derivative[iRow] += differentiate[index(iRow, iCol)]*values[iCol];
    }
  }
  double solution[5] = {0.1359455071179487986e-1, .26211706521167405295,
    .31250000000000000000, .90653768089899526373, 2.1510840365108691368};
  for (int i = 0; i < 5; i++) 
    assert(abs(derivative[i]-solution[i]) < 2.75E-15);
}

void TestLegendre::testValuesToCoefficients() {
  // These were used in the test of the inverse functionality.
  double values[5] = {.68257536428805005463, .74513775143552098982,
    .89843750000000000000, 1.1792615356095370567, 1.7140130029878795533};
  double solution[5] = {1.0, 0.5, 0.25, 0.125, 0.0625};
  double coeffs[5];
  fillCoefficients(values, coeffs);
  for (int i = 0; i < 5; i++)
    assert(abs(coeffs[i]-solution[i]) < 4.0e-16);
}

void TestLegendre::testEvaluate() {
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

void TestLegendre::testInterpolate() {
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
