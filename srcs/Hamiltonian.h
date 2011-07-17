#ifndef HAMILTONIAN_H_
#define HAMILTONIAN_H_
#include "Basis2D.h"
#include "SphericalBasis.h"
#include "ExtrinsicCurvature.h"

/**
 * \brief Conformal Transverse-Traceless Hamiltonian.
 *
 * A class which is responsible for solving and handling the Hamiltonian
 * after the Conformally flat, maximally sliced, transverse traceless
 * decomposition. Uses spherical coordinates \f$(r,\theta)\f$ and maps
 * them onto spectral coordinates \f$(x,u)\f$ using
 * \f[
 *  x = \frac{2r}{r_\mathrm{max}} - 1
 * \f]
 * and
 * \f[
 *  u = \frac{2\theta}{\pi}-1
 * \f]
 * The decomposition used to solve the problem regularly is
 * \f[
 *  \psi = 1 + \psi_\rho + \frac{\sigma(\theta)}{r^p} + u(r,\theta)
 * \f]
 * The differential equation in this case is
 * \f[
 *  r^{p+2}\triangle u - \frac{\beta}{8\sigma^7} + 
 *    \frac{\kappa}{8\left(r^p(1+\psi_\rho + u) + \sigma\right)^7} = 0
 * \f]
 * provided that the function \f$\sigma(\theta)\f$ is coming from the
 * SingularPart object. See documentation of that class for information
 * about the differential equation governing it's behavior.
 */
class Hamiltonian {
  public:
    /**
     * Constructs the basic Hamiltonian object. Takes a SphericalBasis
     * so that derivatives don't need to be recomputed.
     * \param kij The extrinsic curvature associated with the initial
     * data.
     * \param aBasis A pointer to the spherical basis.
     */
    Hamiltonian(ExtrinsicCurvature &kij, SphericalBasis* aBasis);

    /**
     * Constructs the basic Hamiltonian object. Initializes the various
     * basis objects.
     * \param kij The extrinsic curvature associated with the initial
     * data.
     */
    Hamiltonian(ExtrinsicCurvature &kij);

    /**
     * Destructor, needed to delete some new objects.
     */
    ~Hamiltonian();

    /**
     * Sets the bare mass of the spacetime.
     * \param m The bare mass to set.
     */
    void setBareMass(double m);

    /**
     * Gets the bare mass of the spacetime.
     */
    double getBareMass();

    /**
     * Solves the Hamiltonian constraint using Newton Raphson.
     */
    void solve(double tolerance, int maxIterations);

    /**
     * Set the ranks of the basis.
     */
    void setRanks(int N1, int N2);

    /**
     * Prints the remainder term to file.
     */
    void printToFile();

    /**
     * Interpolates the function, then prints the result to file.
     * \param nR The number of radial points to interpolate onto.
     * \param nTheta The number of azimuthal points to interpolate onto.
     */
    void interpolateToFile(int nR, int nTheta);

    /**
     * Sets the maximum radius.
     * \param radius The maximum radius to set.
     */
    void setMaximumRadius(double radius);

    /**
     * Gets the maximum radius.
     * \retval radius The maximum radius.
     */
    double getMaximumRadius();

    /**
     * Interpolates the function, then returns the result.
     * \param nR The number of radial points to interpolate onto.
     * \param nTheta The number of azimuthal points to interpolate onto.
     * \retval A new double[] containing the interpolated result.
     */
    double* interpolateRemainder(int nR, int nTheta);

    /**
     * Toggles whether to be verbose or not.
     * \param verbosity The boolean state of verbosity to set.
     */
    void setVerbose(bool verbosity);

    /**
     * Toggles whether to delete the solution when the NR fails to
     * converge.
     * \param severity The boolean state of severity to set.
     */
    void setSevere(bool severity);

    /**
     * Returns the angular part of the singular term \f$\sigma\f$. The
     * result is a new double[] so remember to delete[] it.
     * \retval A new double[] containing \f$\sigma\f$.
     */
    double* getSingularAngularPart();

    /**
     * Returns the ranks of the basis. Remember to delete[] the returned
     * array.
     * \retval The ranks of the basis.
     */
    int* getRanks();

    /**
     * Returns the power of the singular part of the equation \f$p\f$.
     * \retval The power of the singular part.
     */
    double getSingularPower();

    /**
     * Returns a copy of the remainder, u.
     * \retval u The remainder in the Hamiltonian decomposition.
     */
    double* getRemainder();

  private:
    /**
     * Sets up the data for solving the problem.
     */
    void setUpProblemData();

    /**
     * Performs Newton-Raphson iterations on the Nonlinear constraint
     * equation for the regular function u. Returns true when and if it
     * has converged to a solution, other returns false.
     * \retval Whether the solution converged to the tolerance or not.
     */
    bool newtonRaphson(double tolerance, int maxIterations);

    /**
     * Fills the supplied vector with the residue evaluated at the
     * collocation point. Returns the Euclidean norm of the residue
     * vector.
     * \param residue The residue vector to be filled.
     * \retval The norm of the residue vector.
     */
    double fillResidue(double* residue);

    /**
     * Fills the Jacobian matrix for the nonlinear constraint problem
     * for the regular part of the conformal factor. The indexing is
     * transposed, to prepare for inversion with dgesv.
     * \param jacobian The jacobian matrix to fill.
     */
    void fillJacobian(double* jacobian);

    // TODO: Remember the boundary condition
    /**
     * Computes (if needed) and returns the differential operator + boundary
     * condition operator. In this case, the differential operator is
     * given by \f$r^{p+2}\triangle\f$. The transpose is computed, since
     * this will be used by dgesv for Newton-Raphson iterations.
     * \retval The Laplacian/boundary condition operator.
     */
    const double* getOperator();

    /**
     * The 2D spherical basis.
     */
    SphericalBasis* basis;

    /**
     * The mass.
     */
    double mass;

    /**
     * The verbosity flag.
     */
    bool verbose;

    /**
     * Deletes the solution if it is non-convergent in the
     * Newton-Raphson.
     */
    bool severe;

    /**
     * The maximum radius.
     */
    double rMax;

    /**
     * The power of \f$r\f$ in the singular part of the decomposition.
     */
    double regPower;

    /**
     * The angular part of the singularity: \f$\sigma(\theta)\f$!
     */
    double* sigma;

    /**
     * The numerator of the conformal extrinsic curvature squared:
     * \f$\kappa\f$! See ExtrinsicCurvature for explanation.
     */
    double* kappa;

    /**
     * The O(1) part of kappa, depends only on angle, divided by sigma
     * to the seventh power. See ExtrinsicCurvature for explanation of
     * the beta term.
     */
    double* betaOverSigma7;

    /**
     * \brief The Laplacian operator plus the boundary conditions.
     *
     * The only differential operator here. It
     * is a bit of a misnomer, as it also contains the boundary
     * conditions.
     */
    double* diffOperator;

    /**
     * \brief The non-singular part of the conformal factor.
     */
    double* u;

    /**
     * The residue of the Hamiltonian equation.
     */
    double* residue;

    /**
     * The extrinsic curvature associated with the initial data.
     */
    ExtrinsicCurvature k;

    /**
     * Frees up the new[] doubles. Should be called whenever the rMax,
     * mass, or ranks of the bases change (or a similar change).
     */
    void freeFields();
};

#endif
