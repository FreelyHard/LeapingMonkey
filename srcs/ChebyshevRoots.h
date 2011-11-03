#ifndef CHEBYSHEVROOTS_H_
#define CHEBYSHEVROOTS_H_
#include <cstdlib>
#include "Chebyshev.h"

/**
 * \brief The Chebyshev polynomials using the roots for collocations.
 *
 * A class to represent a truncated set of ChebyshevRoots polynomial
 * basis functions.
 */
class ChebyshevRoots: public Chebyshev {
  public:
    /**
     * The constructor does nothing.
     */
    ChebyshevRoots(): Chebyshev() {};

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
