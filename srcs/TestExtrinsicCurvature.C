#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cassert>
#include "ExtrinsicCurvature.h"

#define MYL 0.1
#define MYP 0.2

double abs(double x) {
  if (x < 0.) return -x;
  return x;
}

class TestExtrinsicCurvature: public ExtrinsicCurvature {
  public:
    TestExtrinsicCurvature(): ExtrinsicCurvature() {};
    void runTests();
  private:
    void testASquared();
};

void TestExtrinsicCurvature::runTests() {
  std::cout << "Running tests on extrinsic curvature.\n";
  int nTests = 0;

  std::cout << "OK. Ran " << nTests << " tests successfully.\n";
}

int main() {
  TestExtrinsicCurvature k;
  k.setSpin(MYL);
  k.setMomentum(MYP);
  srand(time(NULL));
  k.runTests();
}
