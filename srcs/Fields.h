#ifndef FIELDS_H_
#define FIELDS_H_

#include "ExtrinsicCurvature.h"
#include "Hamiltonian.h"
#include "SphericalBasis.h"

/**
 * \brief
 * Encapsulates all the methods necessary for the initial data in the
 * maximally sliced, transverse traceless axisymmetric formulation. 

 * In this case the metric is given by
 * \f[
 *  g_{ij} = \psi^4 \delta_{ij}
 * \f]
 * with \f$f_{ij}\f$ the metric of flat space. The extrinsic curvature
 * is given by
 * \f[
 *  k_{ij} = \psi^2A_{ij}
 * \f]
 * and since the background is flat, one does not need to distinguish
 * between the covariant and contravariant conformal extrinsic
 * curvature. That is
 * \f[
 *  A_{ij} = A^{ij}
 * \f]
 * in Cartesian coordinates. In this case, the coordinates are spherical
 * (although vector and tensor bases are Cartesian). Recall that the
 * conformal factor is decomposed as
 * \f[
 *  \psi(r,\theta) = 1 + \frac{m}{2r} + \frac{\sigma(\theta)}{r^p} + u(r,\theta)
 * \f]
 *
 * This object knows that the Hamiltonian solution uses the
 * ChebyshevRoots in the angular direction...
 */
class Fields {
  public:
    /**
     * Sets up the internal data structures for the particular fields.
     * \param ham The object encapsulating the Hamiltonian constraint.
     * \param kij The extrinsic curvature object.
     */
    Fields(Hamiltonian &ham, ExtrinsicCurvature kij);

    /**
     * Object destructor delete[]'s some new double[]'s.
     */
    ~Fields();

    /**
     * The conformal factor, psi.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The conformal factor.
     */
    double psi(double r, double theta);

    /**
     * The derivative of the conformal factor psi with respect to r.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The specified derivative of the conformal factor.
     */
    double psir(double r, double theta);

    /**
     * The derivative of the conformal factor psi with respect to theta.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The specified derivative of the conformal factor.
     */
    double psit(double r, double theta);

    /**
     * The second derivative of the conformal factor psi with respect to r.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The specified derivative of the conformal factor.
     */
    double psirr(double r, double theta);

    /**
     * The second derivative of the conformal factor psi with respect
     * to theta and r.
     * The conformal factor, psi.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The specified derivative of the conformal factor.
     */
    double psirt(double r, double theta);

    /**
     * The conformal extrinsic curvature A_{xx}.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The conformal extrinsic curvature.
     */
    double Axx(double r, double theta);

    /**
     * The conformal extrinsic curvature A_{xy}.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The conformal extrinsic curvature.
     */
    double Axy(double r, double theta);

    /**
     * The conformal extrinsic curvature A_{yz}.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The conformal extrinsic curvature.
     */
    double Ayz(double r, double theta);

    /**
     * The conformal extrinsic curvature A_{xz}.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The conformal extrinsic curvature.
     */
    double Axz(double r, double theta);

    /**
     * The conformal extrinsic curvature A_{zz}.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The conformal extrinsic curvature.
     */
    double Azz(double r, double theta);

    /**
     * The radial derivative of the conformal extrinsic curvature A_{xx}.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The conformal extrinsic curvature.
     */
    double Axxr(double r, double theta);

    /**
     * The radial derivative of the conformal extrinsic curvature A_{xz}.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The conformal extrinsic curvature.
     */
    double Axzr(double r, double theta);

    /**
     * The radial derivative of the conformal extrinsic curvature A_{zz}.
     * \param r The radial coordinate.
     * \param theta The zenithal coordinate.
     * \retval The conformal extrinsic curvature.
     */
    double Azzr(double r, double theta);

    /**
     * \brief Computes the ADM mass.
     *
     * Used to compute the mass of trumpet spacetimes and also to
     * compute a sort of consistency. In general the ADM
     * mass is given by 
     * \f[
     *  M = \int_{\mathcal{S}_\infty \left(
     *  \frac{\partial g_{ij}}{\partial x^i} - 
     *  \frac{\partial g_{ii}}{\partial x^j}\right) s^j dA}
     * \f]
     * but in spherical symmetry we have.
     * \f[
     *  M = -\int{\psi^3 \frac{\partial\psi}{\partial r} 
     *  r^2\,\sin\theta\,d\theta}
     * \f]
     * \retval The ADM mass evaluated at the maximum radius.
     */
    double computeADMMass();

  private:

    /**
     * The remainder part of the conformal factor.
     */
    double* u;

    /**
     * The radial derivative of the remainder part of the conformal factor.
     */
    double* ur;

    /**
     * The angular derivative of the remainder part of the conformal factor.
     */
    double* ut;

    /**
     * The second radial derivative of the remainder part of the conformal factor.
     */
    double* urr;

    /**
     * The mixed second derivative of the remainder part of the conformal factor.
     */
    double* urt;

    /**
     * The extrinsic curvature.
     */
    ExtrinsicCurvature k;

    /**
     * The mass part of the conformal factor decomposition.
     */
    double mass;

    /**
     * The singular part of the conformal factor decomposition.
     */
    double* sigma;

    /**
     * The derivative of the singular part of the conformal factor decomposition.
     */
    double* sigmap;

    /**
     * The spherical basis.
     */
    SphericalBasis basis;

    /**
     * The power of the singular part \f$p\f$.
     */
    double regPower;

    /**
     * Computes the derivatives of the remainder and sigma.
     */
    void computeDerivatives();
};

#endif
