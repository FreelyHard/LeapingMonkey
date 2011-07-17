#ifndef CONFORMALFACTOR_H_
#define CONFORMALFACTOR_H_ 

/**
 * \brief A class that encapsulates linear conformal factors.
 *
 * The exact form of the conformal factor defined through this class is
 * \f[
 *  \psi_\mathrm{singular} = 1 + \frac{1}{2r} + \frac{dz}{2r^3} + \frac{q}{8r^5}(x^2 + y^2 - 2z^2)
 * \f]
 * it contains psi and its derivatives. All formulae are analytic. So this class only stores information
 * about the moments.
 */
class ConformalFactor {
  public:
    /**
     * The constructor for the conformal factor.
     * \param mass The ADM mass of the spacetime.
     * \param dipole The dipole moment of the spacetime.
     * \param quadrupole The quadrupole moment of the spacetime.
     */
    ConformalFactor(double mass, double dipole, double quadrupole) :
      m(mass), d(dipole), q(quadrupole) {};

    /**
     * Changes the moments of the conformal factor.
     * \param mass The ADM mass of the spacetime.
     * \param dipole The dipole moment of the spacetime.
     * \param quadrupole The quadrupole moment of the spacetime.
     */
    void setMoments(double mass, double dipole, double quadrupole);

    /**
     * The conformal factor itself.
     * \param r The radial coordinate r.
     * \param z The cylindrical height coordinate z.
     * \retval The conformal factor.
     */
    double psi(double r, double z) const;

    /**
     * The derivative of the conformal factor with respect to r.
     * \param r The radial coordinate r.
     * \param z The cylindrical height coordinate z.
     * \retval The derivative of the conformal factor: dpsi/dr.
     */
    double dpsidr(double r, double z) const;

    /**
     * The derivative of the conformal factor with respect to z.
     * \param r The radial coordinate r.
     * \param z The cylindrical height coordinate z.
     * \retval The derivative of the conformal factor: dpsi/dz.
     */
    double dpsidz(double r, double z) const;

    /**
     * The second derivative of the conformal factor with respect to r.
     * \param r The radial coordinate r.
     * \param z The cylindrical height coordinate z.
     * \retval The derivative of the conformal factor: d^2psi/dr^2.
     */
    double d2psidr2(double r, double z) const;

    /**
     * The second derivative of the conformal factor with respect to z.
     * \param r The radial coordinate r.
     * \param z The cylindrical height coordinate z.
     * \retval The derivative of the conformal factor: d^2psi/dz^2.
     */
    double d2psidz2(double r, double z) const;

    /**
     * The second derivative of the conformal factor with respect to rz.
     * \param r The radial coordinate r.
     * \param z The cylindrical height coordinate z.
     * \retval The derivative of the conformal factor: d^2psi/drdz.
     */
    double d2psidrdz(double r, double z) const;

    /**
     * Returns the mass.
     * \retval The mass.
     */
    double getMass();

    /**
     * Returns the dipole moment of the spacetime.
     * \retval The dipole moment.
     */
    double getDipoleMoment();

    /**
     * Returns the quadrupole moment of the spacetime.
     * \retval The quadrupole moment.
     */
    double getQuadrupoleMoment();

  private:
    /**
     * The ADM mass parameter.
     */
    double m;

    /**
     * The dipole moment.
     */
    double d;

    /**
     * The quadrupole moment.
     */
    double q;
};
#endif
