#ifndef TESTLEGENDRE_H_
#define TESTLEGENDRE_H_

#include "Legendre.h"
/**
 * A class to test the Legendre class.
 */
class TestLegendre: public Legendre {
  public:
    /**
     * The constructor does nothing.
     */
    TestLegendre() {};

    /**
     * Runs the various tests on the Legendre class.
     */
    void runTests();

  private:

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
#endif
