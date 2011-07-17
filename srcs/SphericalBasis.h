#ifndef SPHERICALBASIS_H_
#define SPHERICALBASIS_H_
#include "Basis2D.h"
#include "Chebyshev.h"

/**
 * \brief Provides a Chebyshev basis for (2D) spherical coordinates.
 *
 * A class which is responsible for handling a spectral basis on
 * spherical coordinates and building differential operators for those
 * coordinates.
 */
class SphericalBasis {
  public:

    /**
     * The coordinate directions.
     */
    enum Direction {
      COORD1 = 0, /**< Specifies the first coordinate direction. */
      COORD2 = 1  /**< Specifies the second coordinate direction. */ 
    };

    /**
     * Constructs the basic SphericalBasis object. Initializes the various
     * basis objects.
     */
    SphericalBasis();

    /**
     * Destructor, needed to delete some new objects.
     */
    ~SphericalBasis();

    /**
     * Returns the basis responsible for the direction dir.
     * \param dir The desired coordinate direction.
     * \retval The basis used in that coordinate direction.
     */
    Basis* getBasis(Direction dir);

    /**
     * Sets the rank of the two bases.
     * \param radialRank The rank to set in the radial basis.
     * \param zenithalRank The rank to set in the zenithal basis.
     */
    void setRanks(int radialRank, int zenithalRank);

    /**
     * Sets the maximum radius of the domain \f$r_\mathrm{max}\f$.
     */
    void setMaximumRadius(double maximumRadius);

    /**
     * Computes (if needed) and returns the Laplacian + boundary
     * condition operator.
     * In spherical coordinates, the Laplacian is given by
     * \f[
     *  \triangle \equiv \frac{\partial^2}{\partial r^2} +
     *  \frac{1}{r^2}\frac{\partial^2}{\partial\theta^2} +
     *  \frac{2}{r}\frac{\partial}{\partial r} +
     *  \frac{\tan(\theta)}{r^2}\frac{\partial}{\partial\theta}\mathrm{.}
     * \f]
     * The boundary condition in this case is 
     * \f[
     *  \frac{\partial \psi}{\partial r} + \frac{\psi}{r} = \frac{1}{r}
     * \f]
     * although we only set the operator part in this case.
     * \retval The Laplacian/boundary condition operator.
     */
    const double* getLaplacian();

    /**
     * Abstract way of getting the i'th coordinate's collocation points.
     */
    const double* getCoord(Direction dir);

    /**
     * Abstract way of getting the i'th coordinate's differential
     * operator.
     */
    const double* getDiff(Direction dir);

    /**
     * Abstract way of getting the ranks of the individual bases.
     */
    int getRank(Direction dir);

    /**
     * Gives the index into the differential operator matrix given two
     * function indexes.
     * \param iOutput The function index of the output vector.
     * \param iInput The function index of the input vector.
     * \retval The index of the corresponding element in the matrix.
     */
    int matrixIndex(int iOutput, int iInput);

    /**
     * Gives the index of the function given two collocation indexes.
     * \param i The collocation index of the first coordinate.
     * \param j The collocation index of the second coordiante.
     * \retval The corresponding index of the function at the specified
     * collocation point.
     */
    int functionIndex(int i, int j);

    /**
     * For a function u, fills the supplied vector with the coefficients
     * of the function, indexed by function index, where
     * functionIndex(i,j) corresponds to the coefficient of
     * phi_i(x)*psi_j(y) and phi and psi are the basis functions for
     * coordinate 1 and 2 respectively.
     * \param u The function to get the coefficients of.
     * \param coeffs The vector to fill with coefficients.
     */
    void fillCoefficients(const double* u, double* coeffs);

    /**
     * Interpolates the supplied function u onto a tensor product of
     * coordinates, i.e. onto the product set coord1*coord2. Returns a
     * new double[] containing those values.
     * \param u The function to interpolate.
     * \param r The values of the radius to interpolate onto.
     * \param nR The number of points in the radius.
     * \param theta The values of theta to interpolate onto.
     * \param nTheta The number of points in theta.
     * \retval The function u interpolated onto the tensor product
     * (meshgrid) of points.
     */
    double* tensorInterpolate(const double* u, const double* r,
        int nR, const double* theta, int nTheta);

    /**
     * Interpolates a 1d function (a function of direction 1 or 2) onto
     * the supplied desired coordinates. Remember to delete the returned
     * array.
     * \param dir The direction to interpolate in (selects coordinate
     * direction).
     * \param u The function to interpolate.
     * \param x The coordinates to interpolate the function onto.
     * \param nX The number of points to interpolate onto.
     * \retval A new double[] containing the interpolated point.
     */
    double* interpolate1d(Direction dir, const double* u, 
        const double* x, int nX);

  private:

    /**
     * Returns the values of r at the collocation points.
     * \retval An array containing the values of r.
     */
    const double* getR();

    /**
     * Returns the collocation values of theta.
     * \retval An array containing the values of theta.
     */
    const double* getTheta();

    /**
     * Returns the matrix operator \f$\partial/\partial r\f$ given by
     * \f[
     *  \frac{\partial}{\partial r} =
     *  \frac{dx}{dr}\frac{\partial}{\partial x} =
     *  \frac{2}{r_\mathrm{max}}\frac{\partial}{\partial x}
     * \f]
     * where \f$\partial/\partial x\f$ comes from \ref Basis
     * "radialBasis". Since a new double[] is returned, remember to
     * delete it.
     * \retval The matrix operator \f$\partial/\partial r\f$.
     */
    double* getdbydr();

    /**
     * Returns the matrix operator \f$\partial/\partial\theta\f$ given by
     * \f[
     *  \frac{\partial}{\partial\theta} =
     *  \frac{du}{d\theta}\frac{\partial}{\partial u} =
     *  -\sin(\theta)\frac{\partial}{\partial u}
     * \f]
     * where \f$\partial/\partial \theta\f$ comes from \ref Basis
     * "thetaBasis". Since a new double[] is returned, remember to
     * delete it.
     * \retval The matrix operator \f$\partial/\partial \theta\f$.
     */
    double* getdbydtheta();

    /**
     * The 2D basis, tensor product of \f$T_n(x)\f$ with 
     * \f$x = 2r/r_\mathrm{max}-1\f$ and 
     * \f$T_n(u)\f$ with \f$u = 2\theta/\pi - 1\f$.
     */
    Basis2D* basis;

    /**
     * The zenithal basis, with functions \f$P_n(u)\f$ with 
     * \f$u = 2\theta/\pi - 1\f$.
     */
    Chebyshev* thetaBasis;

    /**
     * The radial basis, with functions \f$T_n(x)\f$ with 
     * \f$x = 2r/r_\mathrm{max}-1\f$.
     */
    Chebyshev* radialBasis;

    /**
     * The maximum radius of the domain.
     */
    double rMax;

    /**
     * The radial derivative operator.
     */
    double* dBydr;

    /**
     * The angular derivative operator.
     */
    double* dBydtheta;

    /**
     * The Laplacian operator.
     */
    double* laplacian;

    /**
     * The collocation values of r.
     */
    double* r;

    /**
     * The collocation values of theta.
     */
    double* theta;
};

#endif
