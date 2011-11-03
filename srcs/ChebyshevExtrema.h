#ifndef CHEBYSHEVEXTREMA_H_
#define CHEBYSHEVEXTREMA_H_
#include <cstdlib>
#include "Chebyshev.h"

/**
 * \brief The Chebyshev polynomials using the extrema for collocations.
 *
 * A class to represent a truncated set of ChebyshevExtrema polynomial
 * basis functions.
 */
class ChebyshevExtrema: public Chebyshev {
  public:
    /**
     * The constructor does nothing.
     */
    ChebyshevExtrema(): Chebyshev() {};

    /**
     * Sets the pointers to the abscissas and weights.
     * \param N The number of basis functions to use.
     * \retval Returns 0 if the abscissas exist, -1 otherwise.
     */
    virtual bool setRank(int N);

    virtual double *getAbscissas(int N) const;

    using Basis::getAbscissas;
};
#endif
