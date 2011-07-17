#ifndef TRUMPETCONFORMAL_H_
#define TRUMPETCONFORMAL_H_

/**
 * \brief Codes up the "trumpet" solution.
 */
class TrumpetConformal {
  public:
    /**
     * Constructor, sets the mass.
     */
    TrumpetConformal(double mass): m(mass) {};

    /**
     * Returns the conformal factor.
     */
    double psi(double r);
  private:
    /**
     * The mass.
     */
    double m;

    /**
     * The trumpet areal like coordinate in terms of the isotropic
     * coordinate. Uses Newton's method.
     * \param r The isotropic radius
     * \retval The areal coordinate.
     */
    double trumpetR(double r);

    /**
     * The derivative of the isotropic coordinate with respect to the
     * trumpet coordinate.
     * \param R The trumpet areal like coordinate.
     * \retval The derivative of the istropic radius.
     */
    double drdR(double R);

    /**
     * The isotropic coordinate in terms of the trumpet coordinate.
     * \param R The trumpet areal like coordinate.
     * \retval The istropic radius.
     */
    double isotropicr(double R);
};

#endif
