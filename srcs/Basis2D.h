#ifndef BASIS2D_H_
#define BASIS2D_H_
#include "Basis.h"
//TODO: find out if a forward declaration suffices here...

/**
 * \brief A tensor product of 2 Basis objects, forms a 2d basis.
 *
 * A class to represent a basis for 2d functions. Uses tensor products
 * of 1d basis functions.
 */
class Basis2D {
  public:
    /**
     * Constructs the 2d basis from 2 1d bases. This class doesn't think
     * that it owns the 1d bases, so the caller is responsible for their
     * deletion.
     * \param basisA A pointer to the basis for the first coordinate.
     * \param basisB A pointer to the basis for the second coordinate.
     */
    Basis2D(Basis* basisA, Basis* basisB): basis1(basisA), basis2(basisB),
        differentiationMatrices(NULL) {};

    /**
     * The destructor, to clean up some new double[]'s.
     */
    ~Basis2D();

    /**
     * Returns an array containing the 2 differentiation matrices. The
     * first is d/du1, the second is d/du2. Where u1 and u2 are the
     * first and second coordinates.
     */
    const double* getDifferentiationMatrices();

    /**
     * The coordinate directions.
     */
    enum Direction {
      COORD1 = 0, /**< Specifies the first coordinate direction. */
      COORD2 = 1  /**< Specifies the second coordinate direction. */ 
    };

    /**
     * Differentiates the function in the chosen direction. The
     * direction should be one of the Basis2D::Direction enums.
     * \param direction The direction to differentiate in.
     * \param function The function values at the collocation points.
     * \retval The derivative values at the collocation points.
     */
    double* differentiateCoord(Direction direction, const double* function);

    /**
     * Returns the index into a matrix operator. Where the function is a
     * vector representing the function values at collocation points and
     * its dependence on indexing those collocation points has been
     * delegated through the functionIndex() function.
     * \param outputFunctionIndex The row index of the matrix operator.
     * \param inputFunctionIndex The column index of the matrix
     * operator.
     * \retval The index into the matrix such that if b = Ax,
     * A[matrixIndex] =
     * db_{outputFunctionIndex}/dx_{inputFunctionIndex}.
     */
    int matrixIndex(int outputFunctionIndex, int inputFunctionIndex);

    /**
     * Returns an index into a function vector so that
     * f[functionIndex(iCoord1, iCoord2)] corresponds to the function
     * evaluated at the collocation point (x1_iCoord1, x2_iCoord2). The
     * first coordinate index is the fastest growing, so the indexing is
     * fortran style.
     * \param iCoord1 The index of the first coordinate's collocation
     * point.
     * \param iCoord2 The index of the second coordinate's collocation
     * point.
     * \retval The index into a function vector corresponding to the
     * correct collocations.
     */
    int functionIndex(int iCoord1, int iCoord2);

    /**
     * Interpolates the supplied function u onto a tensor product of
     * coordinates, i.e. onto the product set coord1*coord2. Returns a
     * new double[] containing those values.
     * \param u The function to interpolate.
     * \param coord1 The values of coord1 to interpolate onto.
     * \param n1 The number of points in coord1.
     * \param coord2 The values of coord2 to interpolate onto.
     * \param n2 The number of points in coord2.
     * \retval The function u interpolated onto the tensor product
     * (meshgrid) of points.
     */
    double* tensorInterpolate(const double* u, const double* coord1, int n1,
        const double* coord2, int n2);

    /**
     * Returns the coefficients of the function u, indexed using
     * functionIndex() where the arguments are the power of the
     * coordinate. This is a two stage process, first the values are
     * converted to a collection of 1d coefficients:
     * \f[
     *  \{ u_n \} = \{ u(x_i, y_j) \} \rightarrow \{ \sum_n U_n(y_j) T_n(x) \}
     * \f]
     * the 1d coefficients are treated as 1d function values and then
     * converted to coefficients.
     * \param u The function to return the coefficients of.
     * \param coeffs The vector of coefficients to be filled.
     */
    void fillCoefficients(const double* u, double* coeffs);

    /**
     * Sets the rank of the individual basis objects.
     * \param n1 The value of the rank of basis 1 to set.
     * \param n2 The value of the rank of basis 2 to set.
     */
    void setRanks(int n1, int n2);

  protected:
    
    /**
     * The basis for coordinate direction 1.
     */
    Basis* basis1;

    /**
     * The basis for coordinate direction 2.
     */
    Basis* basis2;

    /**
     * The derivative operators.
     */
    double* differentiationMatrices;
};

#endif
