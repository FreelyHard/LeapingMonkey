#ifndef SPECTRALHORIZON_H_
#define SPECTRALHORIZON_H_
#include "Legendre.h"
#include "ConformalFactor.h"
#include "NewtonRaphson1D.h"

/**
 * \brief An apparent horizon finding class, uses Legendre polynomials
 * to spectrally approximate the horizon.
 *
 * A class that computes apparent horizon using axisymmetric spherical
 * harmonics (the Legendre polynomials). Also computes horizon masses,
 * and can save the horizon to file for plotting. Uses conformal factor
 * to represent the metric, so only conformally flat metrics are
 * supported.
 */
class SpectralHorizon : public NewtonRaphson1D {
  public:
    /**
     * \brief The only constructor.
     *
     * Simply initializes the reference to the conformal factor.
     * \param mass The mass of the spacetime.
     * \param dipole The dipole moment of the spacetime.
     * \param quadrupole The quadrupole moment of the spacetime.
     * \param aBasis The basis to use.
     */
    SpectralHorizon(double mass, double dipole, double quadrupole,
        Legendre* aBasis);

    /**
     * Delete[]'s the fields sigma(p,pp) (and the differentiation
     * matrices) which were created with new double[].
     */
    ~SpectralHorizon();

    /**
     * Changes the moments of the conformal factor.
     * \param m The new mass.
     * \param d The new dipole moment.
     * \param q The new quadrupole moment.
     */
    void setMoments(double m, double d, double q);

    /**
     * Initializes the Legendre basis representation. First it tries to
     * set the basis up, then assuming that passed, initialises space
     * for the function, its derivatives and the differentiation
     * operators.
     * \param n The dimension of the linear representation to use.
     * \retval Returns the boolean of whether the representation was
     * initialized.
     */
    bool initializeRepresentation(int n);

    /**
     * Sets the apparent horizon to some guess.
     * \param guess The vector containing the radii of the horizon.
     */
    void setGuess(double* guess);

    /**
     * Finds the horizon.
     * \param tolerance The tolerance to converge to in the residual.
     * \param maxIterations The maximum number of iterations to perform.
     * \retval The boolean of whether the solution converged or not.
     */
    bool findHorizon(double tolerance, int maxIterations);

    /**
     * Prints the horizon to screen.
     */
    void printHorizon();

    /**
     * Prints the horizon data to a file.
     */
    void printToFile();

    /**
     * Interpolates the horizon and then prints the results to file.
     * \param N The number of points to have on the interpolated curve.
     */
    void interpolateToFile(int N);

    /**
     * Returns the horizon mass: \f$M = \sqrt{A/16\pi}\f$. In this case,
     * \f$A = 2\pi\int_0^\pi{\psi^4r\frac{ds}{d\theta} d\theta}\f$ with
     * \f$\frac{ds}{d\theta} = \sqrt{\sigma^2 +{\sigma'}^2}\f$. The
     * trick is converting to the Legendre polynomial's u coordinate.
     * \retval The horizon mass.
     */
    double getHorizonMass();

  protected:

    /**
     * Here \f$\xi = \cos\theta\f$ so \f$\frac{d\xi}{d\theta} =
     * -\sin\theta = -\sqrt{1 - \xi^2}\f$.
     * \param xi The spectral coordinate \f$\xi\f$.
     * \retval The derivative of the inverse mapping,
     *  \f$\frac{d\xi}{d\theta}\f$. 
     */
    virtual double dInverseMap(double xi);

  private:

    /**
     * The dimension of the representation space.
     */
    int nBasis;

    /**
     * The radius field.
     */
    double* horizonRadius;

    /**
     * The representation space for the horizon.
     */
    Legendre* basis;

    /**
     * Delete[]s the arrays that were once new double[]s.
     */
    void freeMemory();

    /**
     * This computes the residual of the apparent horizon equation.
     * \param xi The cosine of the angle, cos(theta).
     * \param u The radius as a function of angle, theta.
     * \param up The derivative of u w.r.t. theta: \f$du/d\theta\f$
     * \param upp The second derivative of u.
     * \retval The residual of the apparent horizon equation.
     */
    virtual double f(double xi, double u, double up, double upp);

    /**
     * The functional derivative of the apparent horizon equation with
     * respect to the function u.
     * \param xi The cosine of the angle, cos(theta).
     * \param u The radius as a function of angle, theta.
     * \param up The derivative of u w.r.t. theta: \f$du/d\theta\f$
     * \param upp The second derivative of u.
     * \retval The functional derivative df/du.
     */
    virtual double dfdu(double xi, double u, double up, double upp);

    /**
     * The functional derivative of the apparent horizon equation with
     * respect to the function up.
     * \param xi The cosine of the angle, cos(theta).
     * \param u The radius as a function of angle, theta.
     * \param up The derivative of u w.r.t. theta: \f$du/d\theta\f$
     * \param upp The second derivative of u.
     * \retval The functional derivative df/dup.
     */
    virtual double dfdup(double xi, double u, double up, double upp);

    /**
     * The functional derivative of the apparent horizon equation with
     * respect to the function upp.
     * \param xi The cosine of the angle, cos(theta).
     * \param u The radius as a function of angle, theta.
     * \param up The derivative of u w.r.t. theta: \f$du/d\theta\f$
     * \param upp The second derivative of u.
     * \retval The functional derivative df/dupp.
     */
    virtual double dfdupp(double xi, double u, double up, double upp);

    /**
     * The conformal factor. I am too lazy to do some autoreplace on the
     * various phi, phir, phiz, etc. that occur in the Maple output, so
     * this is called indirectly through the above functions.
     */
    ConformalFactor psi;

    /**
     * Creates a new double[nBasis*nBasis] containing the current
     * Jacobian, \f$df_i/du_j\f$, for (i,j) = 0..nBasis-1. Remember to
     * delete[] it.
     * \retval The Jacobian \f$df_i/du_j\f$.
     */
    double* getJacobian();

    /**
     * Takes a single Newton-Raphson step.
     * \retval The total residue. Negative if an error ocurred.
     */
    double singleNewtonRaphson();

};

#endif
