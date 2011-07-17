#ifndef SINGULARPART_H_
#define SINGULARPART_H_

#include "Basis.h"
#include "ExtrinsicCurvature.h"

/**
 * \brief Solves the differential equations governing the singular part of the
 * conformal factor, and returns it.
 *
 * One should recall that the conformal factor is decomposed as
 * \f[
 *  \psi(r,\theta) = 1 + \psi_\rho(r,\theta) + \frac{\sigma(\theta)}{r^p} + u(r,\theta)
 * \f]
 * choosing \f$p=\frac{n-2}{8}\f$ leads to the equation
 * \f[
 *  \frac{\partial^2\sigma}{\partial\theta^2} +
 *      \cot\theta\frac{\partial\sigma}{\partial\theta} 
 *      + p(p-1)\sigma + \frac{\beta}{8\sigma} = 0
 * \f]
 * with \f$\beta\f$ coming from the extrinsic curvature.
 */
class SingularPart {
  public:
    /**
     * Constructs the singular part solver. Needs the zenithal basis and
     * the extrinsic curvature.
     * \param base The zenithal basis.
     * \param kij The extrinsic curvature.
     */
    SingularPart(Basis* base, const ExtrinsicCurvature &kij): 
      basis(base), curvature(kij) {};

    /**
     * Computes the sigma function in the decomposition of the conformal
     * factor. Returns a boolean determining whether or not it converged
     * to the solution.
     * \param sigma On output, the function sigma evaluated at 
     * the zenithal collocation points.
     * \retval Whether or not the solution converged (true).
     */
    bool fillSigma(double* sigma) const;

    /**
     * Returns the power p, such that 1/r^p is the strength of the
     * singularity in the conformal factor.
     * \retval The regularity power.
     */
    double getRegularityPower() const;

  private:

    /**
     * Computes the residue and stores it in the supplied array.
     * \param sigma The current best guess at the solution.
     * \param beta The extrinsic curvature's source function.
     * \param differentialOperator The differential operator of the
     * equation.
     * \param residue A pointer to an array to store the residue.
     * \retval The total residue, the Euclidean norm of the residue.
     */
    double fillResidue(double* sigma, double* beta, 
        double* differentialOperator, double* residue) const;

    /**
     * Returns what is essentially the angular part of the Laplacian, in
     * axisymmetry. The result is a new double[], so remember to
     * delete[] it.
     * \param theta The values of theta at the collocation points.
     * \retval The differential operator.
     */
    double* getOperator(double* theta) const;

    /**
     * Performs Newton-Raphson iterations until sigma converges to the
     * solution.
     * \param sigma On input, the initial guess, on output the actual
     * function sigma;
     * \param beta The beta function defined by the extrinsic curvature.
     * \param diffOp The differential operator part of the PDE.
     * \retval The boolean value of whether it converged (true) or not.
     */
    bool newtonRaphson(double* sigma, double* beta, double* diffOp) const;

    /**
     * The zenithal basis.
     */
    Basis* basis;

    /**
     * The extrinsic curvature.
     */
    ExtrinsicCurvature curvature;
};

#endif
