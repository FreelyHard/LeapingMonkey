#define private protected
#include "NewtonRaphson1DGeneral.h"
#include "ChebyshevExtrema.h"
#include <iostream>
#include <cassert>
#include <cmath>

// The fields I am using are:
//    u1[i] = x^4
//    u2[i] = x
//    u3[i] = x^2/2
//    u4[i] = x^5
//    u5[i] = x^5/5

double abs(double x) {
  return x < 0 ? -x : x;
}

class TestClass : public Monkey::NewtonRaphson1DGeneral {
  public:
    TestClass();
    virtual ~TestClass();
    void runTests();
  protected:
    virtual double equation(int, int);
    virtual double jacobian(int, int, int, int);
    void runGetFieldsTest();
    void runComputeDerivativesTest();
    void runComputeResidueTest();
    void runComputeJacobianTest();
    double bulk1(int);
    double bulk2(int);
    double bulk3(int);
    double bulk4(int);
    double bulk5(int);
};

//    u1[i] = x^4
double TestClass::bulk1(int iCollocation) {
  double* u = getField(2,0);   //    x^4
  double* up = getField(2,1);  // 4 *x^3
  double* upp = getField(2,2); // 12*x^2
  return upp[iCollocation]*pow(coordinates[0][iCollocation], 2)
    - pow(coordinates[0][iCollocation], 1)*up[iCollocation]
    - 8*u[iCollocation];
}

//    u2[i] = x
double TestClass::bulk2(int iCollocation) {
  double* upp = getField(3,2); // 0
  return upp[iCollocation];;
}

//    u3[i] = x^2/2
double TestClass::bulk3(int iCollocation) {
  double* upp = getField(4,2); // 1
  return upp[iCollocation] - 1;
}

//    u4[i] = x^5
double TestClass::bulk4(int iCollocation) {
  double* u = getField(5,0); // x^5
  double* up = getField(5,1); // x^4*5
  double* upp = getField(5,2); // x^3*20
  return upp[iCollocation] - 20*pow(coordinates[1][iCollocation], 3);
}

//    u5[i] = x^5/5
double TestClass::bulk5(int iCollocation) {
  double* u = getField(6,0); // x^5/5
  double* up = getField(6,1); // x^4
  double* upp = getField(6,2); // x^3*4
  return upp[iCollocation] - up[iCollocation]*4/coordinates[1][iCollocation];
}

void TestClass::runTests() {
  int nTests = 0;

  runGetFieldsTest();
  nTests++; std::cout << ".\n";

  runComputeDerivativesTest();
  nTests++; std::cout << ".\n";

  runComputeResidueTest();
  nTests++; std::cout << ".\n";

  runComputeJacobianTest();
  nTests++; std::cout << ".\n";

  std::cout << "Performed " << nTests << " tests successfully.\n";
}

void TestClass::runComputeDerivativesTest() {
  // Some simple polynomial tests will be exact.
  int nBasis = bases[0]->getRank();
  double *u1 = getField(2, 0);
  double *u2 = getField(3, 0);
  double *u3 = getField(4, 0);
  for (int i = 0; i < nBasis; i++) {
    u1[i] = pow(coordinates[0][i], 4);
    u2[i] = coordinates[0][i];
    u3[i] = pow(coordinates[0][i], 2)*0.5;
  }
  nBasis = bases[1]->getRank();
  double *u4 = getField(5, 0);
  double *u5 = getField(6, 0);
  for (int i = 0; i < nBasis; i++) {
    u4[i] = pow(coordinates[1][i], 5);
    u5[i] = pow(coordinates[1][i], 5)/5.;
  }
  computeDerivatives();

  double *u4p = getField(5, 1);
  double *u4pp = getField(5, 2);
  double *u5p = getField(6, 1);
  double *u5pp = getField(6, 2);
  for (int i = 0; i < nBasis; i++) {
    assert(abs(u4p[i] - pow(coordinates[1][i], 4)*5) < 1.0e-13);
    assert(abs(u4pp[i] - pow(coordinates[1][i], 3)*20) < 1.0e-13);
    assert(abs(u5p[i] - pow(coordinates[1][i], 4)) < 1.0e-13);
    assert(abs(u5pp[i] - pow(coordinates[1][i], 3)*4) < 1.0e-13);
  }
  nBasis = bases[0]->getRank();
  double *u1p = getField(2, 1);
  double *u2p = getField(3, 1);
  double *u3p = getField(4, 1);
  double *u1pp = getField(2, 2);
  double *u2pp = getField(3, 2);
  double *u3pp = getField(4, 2);
  for (int i = 0; i < nBasis; i++) {
    assert(abs(u1p[i] - pow(coordinates[0][i], 3)*4) < 1.0e-13);
    assert(abs(u1pp[i] - pow(coordinates[0][i], 2)*12) < 1.0e-13);
    assert(abs(u2p[i] - 1.) < 1.0e-13);
    assert(abs(u2pp[i]) < 1.0e-13);
    assert(abs(u3p[i] - pow(coordinates[0][i], 1)) < 1.0e-13);
    assert(abs(u3pp[i] - 1.0) < 1.0e-13);
  }
}

void TestClass::runComputeResidueTest() {
  computeResidue();
  for (int i = 0; i < problemDimension; i++) {
    assert(abs(residue[i]) < 1.5e-14);
  }
}

void TestClass::runGetFieldsTest() {
  double *myFieldsPointer = getField(2,0);
  int i = 0;
  int iField = numberOfAuxiliaries;
  for (int iDomain = 0; iDomain < numberOfDomains; iDomain++) {
    int nBasis = bases[iDomain]->getRank();
    for (int iField = 0; iField < fieldsPerDomain[iDomain]; iField++) {
      for (int id = 0; id < 3; id++) {
        for (int j = 0; j < nBasis; j++) {
          myFieldsPointer[i] = i++;
        }
      }
    }
  }
  // Assert that i is the number of elements in fields.
  assert(i == 3*(5*3 + 6*2));
  auxiliaries[0] = 1000;
  auxiliaries[1] = 1111;

  int old = -5;
  // the difference in initial elements is the size of the vector.
  for (iField = 2; iField < 7; iField++) {
    for (int id = 0; id < 3; id++) {
      int inew = (int)getField(iField, id)[0];
      assert(inew - old == 5 || inew - old == 6);
      old = inew;
    }
  }
  // The last element is 80.
  assert((int)getField(6, 2)[5] == 80);
}

double TestClass::equation(int iEquation, int iCollocation) {
  switch(iEquation) {
    case 2:
      return bulk1(iCollocation);
    case 3:
      return bulk2(iCollocation);
    case 4:
      return bulk3(iCollocation);
    case 5:
      return bulk4(iCollocation);
    case 6:
      return bulk5(iCollocation);
  }
  return 0;
}

double TestClass::jacobian(int iDomain, int iVariable,
    int iDerivative, int iCollocation)  {
  return 0;
}

TestClass::~TestClass() {
  delete[] fieldsPerDomain;
  for (int i = 0; i < numberOfDomains; i++) {
    delete bases[i];
  }
  delete[] bases;
}

void TestClass::runComputeJacobianTest() {
  computeJacobian();
  for (int iOut = 0; iOut < problemDimension; iOut++) {
    double notZero = 0.0;
    for (int iIn = 0; iIn < problemDimension; iIn++) {
      int iMatrix = iOut*problemDimension + iIn;
      notZero += totalJacobian[iMatrix];
    }
    assert(abs(notZero) > 1.0e-16);
  }
}

TestClass::TestClass(): NewtonRaphson1DGeneral()
{
  // The number of domains.
  numberOfDomains = 2;

  // The number of fields per domain.
  fieldsPerDomain = new int[2];
  fieldsPerDomain[0] = 3;
  fieldsPerDomain[1] = 2;

  // The bases for the individual domains.
  bases = new Basis*[2];
  bases[0] = new ChebyshevExtrema();
  bases[0]->setRank(5);
  bases[1] = new ChebyshevExtrema();
  bases[1]->setRank(6);

  // The number of auxiliary conditions.
  numberOfAuxiliaries = 2;

  // Initialise template method internals.
  initialise();
}

int main() {
  std::cout << "Testing NewtonRaphson1DGeneral class.\n";

  TestClass test;
  test.runTests();
}
