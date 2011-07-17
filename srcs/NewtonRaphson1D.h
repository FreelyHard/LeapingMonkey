#ifndef NEWTONRAPHSON1D_H_
#define NEWTONRAPHSON1D_H_

#include "Basis.h"

/**
 * \brief An abstract object that encapsulates the logic necessary to
 * solve one dimensional second order boundary value problems.
 *
 * One may note that this is an example of the "Template Method" design
 * pattern.
 *
 * The constructor and Jacobian methods must be supplied by the derived
 * class. The actual solver function can be overridden if there are
 * problem dependent constants that the user may want to at solution
 * time. The rest is handled internally. The matrix is set up, the
 * residue is computed, the iterations are performed, etc.
 *
 * It is possible to set up the problem so that the Jacobians involve
 * the spectral coordinates directly but the derivatives are computed in
 * terms of the physical coordinates. In this case, one defines the
 * functions f, dfdu, dfdup, dfdupp in terms of \f$\xi\f$, leaves the
 * mapping untouched (so that calls to f and so on are not mapped) but
 * overrides dInverseMap (so that the differentiation matrices are
 * computed in terms of the physical coordinates). This is perhaps
 * ungainly and should perhaps only be used if the use of the spectral
 * coordinates makes computations much easier.
 *
 * There is a fairly strong ordering of the private routines. Before the
 * residue is called, the derivatives must have been computed. Before
 * the derivatives can be computed the differentiation matrices must
 * have been computed. Before the a single Newton-Raphson step, the
 * residue must have been computed.
 *
 * The abstract set of equations that this object solves is
 * \f[
 *  \begin{array}{rcl}
 *    f\left(x, u, \frac{du}{dx}, \frac{d^2u}{dx^2}\right) = 0
 *    & \mathrm{ for } & x \in [a,b]\\
 *    g\left(u, \frac{du}{dx}\right) = 0 &\mathrm{ for }& x = a\\
 *    h\left(u, \frac{du}{dx}\right) = 0 &\mathrm{ for }& x = b
 *    \ \mathrm{or}\ a
 *  \end{array}
 * \f]
 * here the functions \f$g\f$ and \f$h\f$ should be considered boundary
 * conditions, (or initial conditions if \f$h\f$ is defined at
 * \f$x=a\f$).
 */
class NewtonRaphson1D {
  public:

    /**
     * Used to specify whether the function g and h are boundary
     * conditions, initial conditions or unused.
     */
    enum BoundaryType {
      NONE = 0,     /**< Uses collocation normally. */
      INITIAL = 1,  /**< Equation is an initial condition. */
      BOUNDARY = 3, /**< Equation is a boundary condition. */
    };

    /**
     * \brief The only constructor. Sets the pointer to a basis.
     *
     * \param aBasis The polynomial basis to use for the solution.
     */
    NewtonRaphson1D(Basis* aBasis): basis(aBasis), diff(NULL),
      doubleDiff(NULL), sigmap(NULL), sigmapp(NULL), residue(NULL),
      gType(NONE), hType(NONE) {};

    /**
     * Delete[]'s the fields the differentiation
     * matrices) which were created with new double[].
     */
    virtual ~NewtonRaphson1D();

    /**
     * Performs multiple Newton Raphson iterations until the solution is
     * acquired. The default assumes that there are no problem dependent
     * constants and just calls computeSolution(double, int).
     * \param u On input, the initial guess, on output, the solution.
     * \param constants An object storing various problem dependent
     * constants.
     * \param tolerance The tolerance to converge to in the residual.
     * \param maxIterations The maximum number of iterations to perform.
     * \retval converged The boolean of whether the solution converged or not.
     */
    virtual bool computeSolution(double* u, void* constants, 
        double tolerance, int maxIterations);

    /**
     * Computes the total residue, the Euclidean norm of the residue
     * vector.
     * \retval residue The Euclidean norm of the residue vector.
     */
    double getTotalResidue();

    /**
     * Returns a pointer to the residue. Returns NULL if the solution
     * has not been computed.
     * \retval residue The residue of the solution with respect to the
     * differential equation /f$f/f$.
     */
    const double* getResidue();

    /**
     * Returns a pointer to the derivative or NULL if the solution has
     * not been computed.
     * \retval sigmap The derivative of the solution.
     */
    const double* getDerivative();

    /**
     * Specifies whether to use boundary values, initial conditions or
     * nothing.
     * \param gType Whether the function g is a boundary condition...
     * \param hType Whether the function h is a boundary condition...
     */
    void registerBoundaries(BoundaryType gType, BoundaryType hType);

  protected:

    /**
     * Sets the pointer to the solution vector. Should be used by
     * computeSolution.
     * \param u The vector to store the solution.
     */
    void setPointerToSolution(double* u);

    /**
     * Specifies the state of function g.
     */
    BoundaryType gType;

    /**
     * Specifies the state of function h.
     */
    BoundaryType hType;

    /**
     * Performs multiple Newton Raphson iterations until the solution is
     * acquired.
     * \param tolerance The tolerance to converge to in the residual.
     * \param maxIterations The maximum number of iterations to perform.
     * \retval converged The boolean of whether the solution converged or not.
     */
    bool computeSolution(double tolerance, int maxIterations);

    /**
     * The mapping between the spectral coordinates and the problem
     * dependent coordinates. The default assumes that the problem
     * coordinates are the spectral coordinates.
     * \param xi The spectral coordinates, for instance \f$\xi \in
     * [-1,1]\f$ for a Chebyshev basis.
     * \retval x The physical coordinates, whatever those might be.
     */
    virtual double map(double xi);

    /**
     * The derivative of the inverse mapping between the spectral
     * coordinates and the problem dependent coordinates. The default
     * assumes that the problem coordinates are the spectral coordinates.
     * This is necessary so that the differentiation matrices can be set
     * up by chain rule, since \f$\frac{d}{dx} =
     * \frac{d\xi}{dx}\frac{d}{d\xi}\f$ and the spectral basis can
     * compute \f$\frac{d}{d\xi}\f$.
     * \param xi The spectral coordinates, for instance \f$\xi \in
     * [-1,1]\f$ for a Chebyshev basis.
     * \retval dxidx The derirvative of the physical coordinates
     * \f$\frac{d\xi}{dx}\f$, whatever those might be.
     */
    virtual double dInverseMap(double xi);

    /**
     * The differential equation.
     * \param x The coordinate.
     * \param u The function to solve for.
     * \param up The derivative of the function u.
     * \param upp The second derivative of u.
     * \retval The functional derivative df/dup.
     */
    virtual double f(double x, double u, double up, double upp) = 0;

    /**
     * The functional derivative of the differential equation with
     * respect to the function u.
     * \param x The coordinate.
     * \param u The function to solve for.
     * \param up The derivative of the function u.
     * \param upp The second derivative of u.
     * \retval The functional derivative df/dup.
     */
    virtual double dfdu(double x, double u, double up, double upp) = 0;

    /**
     * The functional derivative of the differential equation with
     * respect to the function up.
     * \param x The coordinate.
     * \param u The function to solve for.
     * \param up The derivative of the function u.
     * \param upp The second derivative of u.
     * \retval The functional derivative df/dup.
     */
    virtual double dfdup(double x, double u, double up, double upp) = 0;

    /**
     * The functional derivative of the differential equation with
     * respect to the function upp.
     * \param x The coordinate.
     * \param u The function to solve for.
     * \param up The derivative of the function u.
     * \param upp The second derivative of u.
     * \retval The functional derivative df/dupp.
     */
    virtual double dfdupp(double x, double u, double up, double upp) = 0;

    /**
     * The left boundary condition.
     * \param u The function value at the left boundary.
     * \param up The derivative of the function at the left boundary.
     * \retval g The boundary condition on the left.
     */
    virtual double g(double u, double up);

    /**
     * The left boundary condition.
     * \param u The function value at the left boundary.
     * \param up The derivative of the function at the left boundary.
     * \retval g The boundary condition on the left.
     */
    virtual double dgdu(double u, double up);

    /**
     * The left boundary condition.
     * \param u The function value at the left boundary.
     * \param up The derivative of the function at the left boundary.
     * \retval g The boundary condition on the left.
     */
    virtual double dgdup(double u, double up);

    /**
     * The right boundary condition.
     * \param u The function value at the left boundary.
     * \param up The derivative of the function at the left boundary.
     * \retval h The boundary condition on the left.
     */
    virtual double h(double u, double up);

    /**
     * The right boundary condition.
     * \param u The function value at the left boundary.
     * \param up The derivative of the function at the left boundary.
     * \retval h The boundary condition on the left.
     */
    virtual double dhdu(double u, double up);

    /**
     * The left boundary condition.
     * \param u The function value at the left boundary.
     * \param up The derivative of the function at the left boundary.
     * \retval g The boundary condition on the left.
     */
    virtual double dhdup(double u, double up);

    /**
     * The representation space for the numerical solution.
     */
    Basis* basis;

  private:

    /**
     * The residue.
     */
    double* residue;

    /**
     * The pointer to the solution, owned by caller.
     */
    double* sigma;

    /**
     * The derivative of the solution.
     */
    double* sigmap;

    /**
     * The second derivative of the solution.
     */
    double* sigmapp;

    /**
     * The matrix to differentiate functions. Specifically up = diff*u.
     */
    double* diff;

    /**
     * The matrix to twice differentiate functions. Specifically, upp =
     * doubleDiff*u. 
     */
    double* doubleDiff;

    /**
     * Computes the differentiation matrices.
     */
    void computeDifferentiationMatrices();

    /**
     * From the function, computes derivatives.
     */
    void computeDerivatives();

    /**
     * Takes a single Newton-Raphson step.
     * \retval The total residue. Negative if an error ocurred.
     */
    double singleNewtonRaphson();

    /**
     * Creates a new double[nBasis*nBasis] containing the current
     * Jacobian, \f$df_i/du_j\f$, for (i,j) = 0..nBasis-1. Remember to
     * delete[] it.
     * \retval jac The Jacobian \f$df_i/du_j\f$.
     */
    double* getJacobian();

    /**
     * Delete[]s the arrays that were once new double[]s.
     */
    void freeMemory();

};

#endif
