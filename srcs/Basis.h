#ifndef BASIS_H_
#define BASIS_H_
#include <cstdlib>

/**
 * \brief A 1d polynomial basis for the representation of truncated
 * approximations of functions.
 *
 * A class to represent a truncated set of Basis polynomial
 * basis functions. The abscissas and weights are pulled out of the
 * BasisAbscissas namespace where they are saved from maple output.
 */
class Basis {
  public:
    /**
     * The constructor does nothing. This is necessary as I don't want
     * to throw an exception from the constructor in the case that the
     * chosen number of basis functions does not have collocation points
     * implemented in the BasisAbscissas namespace.
     */
    Basis(): nBasis(0), valuesToCoefficients(NULL),
      differentiationMatrix(NULL), abscissas(NULL) {};

    /**
     * Deletes the differentiation matrix, if it has been computed.
     */
    virtual ~Basis();

    /**
     * Sets the pointers to the abscissas and weights. Returns a
     * non-zero return code if the abscissas are not implemented in
     * BasisAbscissas.
     * \param N The number of basis functions to use.
     * \retval Returns 0 if the abscissas exist, -1 otherwise.
     */
    virtual bool setRank(int N) = 0;

    /**
     * Returns the abscissas, provided they have been set. Otherwise
     * returns NULL.
     * \retval The array of abscissas.
     */
    const double* getAbscissas();

    /**
     * Returns a pointer to the differentiation matrix. The
     * differentiation matrix takes a vector of function values at the
     * collocation points and returns a vector containing the function's
     * derivative values at the collocation points.
     * \retval The matrix which differentiates functions (represented by
     * vectors).
     */
    const double* getDifferentiationMatrix();
    
    /**
     * The matrix operator which transforms between values at the
     * collocation points and the coefficients of the corresponding
     * Basis polynomial interpolant. Remember to delete[] the matrix
     * when you are finished with it!
     */
    const double* getValuesToCoefficientsMatrix();

    /**
     * Interpolates the function represented by values at the nPoints
     * points in x. The returned array is a new double[] so remember to
     * delete[] it.
     * \param values The function values at the collocation points.
     * \param x The points to interpolate the functioin at.
     * \param nPoints The number of points to interpolate.
     * \retval The function evaluated at the points x.
     */
    double* interpolate(const double* values, const double* x, int nPoints);

    /**
     * Returns the rank of the basis. The number of basis functions, the
     * number of collocation points, etc.
     * \retval The dimension of the function space spanned by the basis.
     */
    int getRank() const;

    /**
     * Returns the index into the matrices.
     * \param iRow The row index.
     * \param iCol The column index.
     * \retval The index of coefficientsToValues etc.
     */
    int index(int iRow, int iCol);

    /**
     * From a vector of values, fills the vector of coefficients.
     * \param values The values of the funciton at the collocation
     * points.
     * \param coefficients The vector of coefficients to be filled.
     */
    void fillCoefficients(const double* values, double* coefficients);

    /**
     * Integrates the function whose value at the collocation points are
     * given by values.
     * \param values The value of the function at the abscissas.
     * \retval The integral of the function, using Gaussian quadrature.
     */
    double integrate(const double* values);

  protected:

    /**
     * The parameter alpha in the recursion relation for the family of
     * polynomials:
     * \f[
     *   \phi_n(x) + \alpha_n(x)\phi_{n-1}(x) + \beta_n(x)\phi_{n-2}(x)
     *   = 0
     * \f]
     * \param n The order of the highest order basis function.
     * \param x The coordinate to evaluate the coefficient at.
     * \retval The parameter alpha.
     */
    virtual double alpha(int n, double x) = 0;

    /**
     * The parameter beta in the recursion relation for the family of
     * polynomials:
     * \f[
     *   \phi_n(x) + \alpha_n(x)\phi_{n-1}(x) + \beta_n(x)\phi_{n-2}(x)
     *   = 0
     * \f]
     * \param n The order of the highest order basis function.
     * \param x The coordinate to evaluate the coefficient at.
     * \retval The parameter beta.
     */
    virtual double beta(int n, double x) = 0;

    /**
     * Evaluates the function represented by the coefficients at the
     * point x using the recursion relation between the polynomials to
     * minimize the number of computations at each step.
     * \param x The coordinate to evaluate the point at.
     * \param coeffs The coefficients defining the function.
     * \retval The function evaluated at the point x.
     */
    virtual double evaluate(double x, double* coeffs) = 0;

    /**
     * The individual basis functions. Not meant to be used
     * directly as they are low-level for doing spectral methods.
     * \param n The order of the basis function.
     * \param x The coordinate to evaluate the basis function at.
     * \retval \phi_n(x) the Basis polynomial of n'th order evaluated at
     * x.
     */
    virtual double function(int n, double x) const = 0;

    /**
     * Computes and returns the matrix operator which transforms
     * coefficients of the polynomials to values at the collocation
     * points (which are the coefficients of the corresponding cardinal
     * functions). Remember to delete[] the matrix when you are finished
     * with it.
     * \retval The coefficient-to-values matrix operator.
     */
    virtual double* coefficientsToValuesMatrix();

    /**
     * Builds and returns the matrix that takes a vector of coefficients
     * of a function and returns the coefficients of the derivative of
     * that function.
     * \retval The coefficients-to-coefficients-of-derivative operator.
     */
    virtual double* coefficientsOfDerivativeMatrix() = 0;

    /**
     * The number of basis functions.
     */
    int nBasis;

    /**
     * The abscissas (first half) and weights (second half).
     */
    const double* abscissas;

    /**
     * The values to coefficients matrix.
     */
    double* valuesToCoefficients;

    /**
     * The matrix which differentiates (vector represented) functions.
     */
    double* differentiationMatrix;
};
#endif
