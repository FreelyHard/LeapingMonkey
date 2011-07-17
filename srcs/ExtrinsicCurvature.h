#ifndef EXTRINSICCURVATURE_H_
#define EXTRINSICCURVATURE_H_

/**
 * \brief Encodes the methods for computing the extrinsic curvature in
 * the conformal transverse tracelss decomposition.
 *
 * The class assumes axisymmetry, so things are simplified a bit.
 */
class ExtrinsicCurvature {
  public:
    /**
     * Constructor, sets the linear and angular momentum to zero.
     */
    ExtrinsicCurvature(): sz(0.0), pz(0.0), M(0.0), Q(0.0), cQ(0.0) {};

    /**
     * Sets the spin of the spacetime.
     * \param S The spin to set.
     */
    void setSpin(double S);

    /**
     * Returns the spin of the extrinsic curvature.
     */
    double getSpin();

    /**
     * Sets the spin curl moment.
     * \param cq The curl moment to set.
     */
    void setSpinCurl(double cq);

    /**
     * Gets the spin curl moment.
     * \retval cQ The spin curl moment.
     */
    double getSpinCurl();

    /**
     * Sets the spin quadrupole moment.
     * \param q The quadrupole moment to set.
     */
    void setQuadrupole(double q);

    /**
     * Gets the spin quadrupole moment.
     * \retval Q The spin quadrupole moment.
     */
    double getQuadrupole();

    /**
     * Sets the momentum of the spacetime.
     * \param P The momentum to set.
     */
    void setMomentum(double P);

    /**
     * Returns the momentum of the extrinsic curvature.
     */
    double getMomentum();

    /**
     * Sets the mass of the spacetime, through the extrinsic curvature.
     * \param m The mass of the spacetime to set.
     */
    void setMass(double m);

    /**
     * Returns the mass of the extrinsic curvature.
     */
    double getMass();

    /**
     * Get the conformal extrinsic curvature.
     * \param x The coordinate x, the first coordinate.
     * \param y The coordinate y, the second coordinate.
     * \param z The coordinate z, the third coordinate.
     * \retval A new double[6] containing Aij use indexAij to get
     * quantities from this. Remember to delete[] the result.
     */
    double* getAij(double x, double y, double z);

    /**
     * A helper function to abstract the indexing. Will return garbage
     * and you will overwrite/access faulty memory if the indices are
     * out of bounds!
     * \param i The first tensor index, should be in [1,2,3].
     * \param j The second tensor index, should be in [1,2,3].
     * \retval The index into the returned Aij from getAij.
     */
    int indexAij(int i, int j);

    /**
     * Returns the value of \f$A^{ij}A_{ij}r^n\f$ at the requested point.
     * \param r The radius of the point to evaluate at.
     * \param theta The zenithal angle to evaluate at.
     * \retval The numerator of the conformal extrinsic curvature squared.
     */
    double kappa(double r, double theta);

    /**
     * Returns the power of the 1/r^n type singularity in A^2.
     * \retval The power n of the 1/r^n singularity.
     */
    int getSingularityPower() const;

    /**
     * Returns a vector containing the values of Beta. The function beta
     * is implicitly defined as
     * \f[
     *  A^{ij}A_{ij} \equiv \frac{\kappa(r,\theta)}{r^n} \equiv
     *    \frac{\sum_{i=1}^{n-1}r^i\alpha(\theta) + \beta(\theta)}{r^n}
     * \f]
     * with n chosen such that \f$\beta\neq 0\f$. Remember to delete[]
     * the returned array.
     * \param theta The values of theta at which to compute beta.
     * \param nTheta The number of collocations: the rank of the basis.
     * \retval A new double[] containing the values of beta.
     */
    double* getBeta(const double* theta, int nTheta) const;

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

  private:
    /**
     * The angular momentum.
     */
    double sz;

    /**
     * The linear momentum.
     */
    double pz;

    /**
     * The mass for trumpet type initial data.
     */
    double M;

    /**
     * The sign of the mass term.
     */
    double xi;

    /**
     * The spin quadrupole moment.
     */
    double Q;

    /**
     * The spin curl term.
     */
    double cQ;

};

#endif
