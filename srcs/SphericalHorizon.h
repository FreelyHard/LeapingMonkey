#ifndef SPHERICALHORIZON_H_
#define SPHERICALHORIZON_H_

#include "NewtonRaphson1D.h"
#include "Fields.h"

/**
 * \brief Encapsulates the code necessary to find an apparent horizon on a
 * maximally sliced, conformally flat, transverse tracelless 3-manifold.
 */
class SphericalHorizon: public NewtonRaphson1D {
  public:

    /**
     * Constructor, initializes the NewtonRaphson solver and sets the
     * fields.
     * \param aBasis The basis to use for finding the horizon.
     * \param initialData The conformal factor and extrinsic curvature.
     */
    SphericalHorizon(Basis* aBasis, Fields* initialData);

    /**
     * Computes the selected (inner or outer) apparent horizon.
     * \param u On input, the initial guess, on output the horizon.
     * \param selector The selector for the inner (-1.0) or outer (+1.)
     * horizon. Should be a single double.
     * \param tolerance The residual on the Newton-Raphson to tolerate.
     * \param maxIterations The maximum number of iterations to perform
     * in the Newton-Raphson iteration.
     * \retval converged Whether the Newton-Raphson converged or not.
     */
    virtual bool computeSolution(double* u, void* selector, double tolerance,
        int maxIterations);

    /**
     * \brief Computes the irreducible horizon mass.
     *
     * In this context the horizon mass \f$M = \sqrt{A/16\pi}\f$ 
     * is given in terms of the area:
     * \f[
     *  A =  2\pi \int_0^\pi {r\psi^4\sqrt{\left(\sigma'\right)^2
     *  + r^2}\;\sin{\theta} \; d\theta}.
     * \f]
     * \param r The horizon found using computeSolution().
     * \retval M_H The horizon mass \f$M = \sqrt{A/16\pi}\f$.
     */
    double computeHorizonMass(double* r);

    /**
     * \brief Computes the horizon spin
     *
     * In this context the horizon spin is given as
     * \f[
     *  J =  2\pi \int_0^\pi {r\psi^4\sqrt{\left(\sigma'\right)^2
     *  + r^2}\;\sin{\theta} \; d\theta}.
     * \f]
     * \param r The horizon found using computeSolution().
     * \retval J The angular momentum of the horizon.
     */
    double computeHorizonSpin(double* r);

  protected:

    /**
     * Computes the residual of the apparent horizon equation. Produced
     * from maple, see the latex documentation also.
     * \param u The radial function, \f$\sigma\f$.
     * \param up The derivative of the radial function
     *  \f$\frac{d\sigma}{d\theta}\f$
     * \param upp The second derivative of the radial function
     *  \f$\frac{d^2\sigma}{d\theta^2}\f$
     * \param theta The zenithal angle.
     * \retval The residual of the apparent horizon equation.
     */
    virtual double f(double theta, double u, double up, double upp);

    /**
     * Computes the derivative of the apparent horizon equation with
     * respect to \f$\sigma\f$. Produced from maple.
     * \param u The radial function, \f$\sigma\f$.
     * \param up The derivative of the radial function
     *  \f$\frac{d\sigma}{d\theta}\f$
     * \param upp The second derivative of the radial function
     *  \f$\frac{d^2\sigma}{d\theta^2}\f$
     * \param theta The zenithal angle.
     * \retval The derivative of the apparent horizon equation with
     * respect to \f$\sigma\f$.
     */
    virtual double dfdu(double theta, double u, double up, double upp);

    /**
     * Computes the derivative of the apparent horizon equation with
     * respect to \f$\sigma'\f$. Produced from maple.
     * \param u The radial function, \f$\sigma\f$.
     * \param up The derivative of the radial function
     *  \f$\frac{d\sigma}{d\theta}\f$
     * \param upp The second derivative of the radial function
     *  \f$\frac{d^2\sigma}{d\theta^2}\f$
     * \param theta The zenithal angle.
     * \retval The derivative of the apparent horizon equation with
     * respect to \f$\sigma' = \frac{d\sigma}{d\theta}\f$.
     */
    virtual double dfdup(double theta, double u, double up, double upp);

    /**
     * Computes the derivative of the apparent horizon equation with
     * respect to \f$\sigma''\f$. Produced from maple.
     * \param u The radial function, \f$\sigma\f$.
     * \param up The derivative of the radial function
     *  \f$\frac{d\sigma}{d\theta}\f$
     * \param upp The second derivative of the radial function
     *  \f$\frac{d^2\sigma}{d\theta^2}\f$
     * \param theta The zenithal angle.
     * \retval The derivative of the apparent horizon equation with
     * respect to \f$\sigma'' = \frac{d^2\sigma}{d\theta^2}\f$.
     */
    virtual double dfdupp(double theta, double u, double up, double upp);

    /**
     * Maps the spectral coordinates onto the physical ones. In this
     * case, the interval \f$[-1,1]\f$ is mapped onto \f$[0,\pi]\f$.
     * \param xi The spectral coordinate \f$\xi\in [-1,1]\f$.
     * \retval theta The zenithal angle \f$\theta \in [0,\pi]\f$.
     */
    virtual double map(double xi);

    /**
     * The derivative of the inverse map. A constant.
     * \param xi The spectral coordinate, unused.
     * \retval dxidtheta The derivative of the inverse mapping 
     * \f$\frac{d\xi}{d\theta}\f$
     */
    virtual double dInverseMap(double xi);

  private:

    /**
     * Stores which horizon is being found.
     */
    double chi;

    /**
     * The initial data fields: the conformal factor \f$\psi\f$ and the
     * extrinsic curvature \f$k_{ij}\f$ and \f$A_{ij}\f$.
     */
    Fields* fields;
};

#endif
