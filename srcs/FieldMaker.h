#ifndef FIELDMAKER_H_
#define FIELDMAKER_H_

#include "Fields.h"
#include "SphericalBasis.h"

/**
 * \brief Produces Fields objects.
 *
 * Useful so that users do not need to call individual
 * routines responsible for solving differential equations.
 */
class FieldMaker {
  public:
    /**
     * \brief Constructs the FieldMaker object. 
     *
     * Needs to know the numerical parameters associated with the 
     * Newton-Raphson iterations.
     * \param tolerance_ The magnitude of the residual to tolerate in
     * the Newton-Raphson iteration.
     * \param maxIterations_ The maximum number of iterations to perform
     * in the Newton-Raphson iteration.
     * \param nRadial_ The number of basis functions to use in the radial
     * direction;
     * \param nTheta_ The number of basis functions to use in the
     * zenithal direction.
     * \param rMax_ The maximum radius.
     */
    FieldMaker(double tolerance_, double maxIterations_, 
        int nRadial_, int nTheta_, double rMax_);

    /**
     * Destructor, kills the basis.
     */
    ~FieldMaker();

    /**
     * Sets the verbosity level of computations.
     * \param verbosity Whether to print statements (true) or not.
     */
    void setVerbosity(bool verbosity);

    /**
     * Sets the severity of the Hamiltonian.
     * \param severity The severity to set.
     */
    void setSevere(bool severity);

    /**
     * \brief Sets the tolerance of the Newton-Raphson iteration.
     * \param tolerance_ The magnitude of the residual to tolerate in
     * the Newton-Raphson iteration.
     */
    void setTolerance(double tolerance_);

    /**
     * \brief Sets the maximum number of iterations in the
     * Newton-Raphson iterations.
     * \param maxIterations_ The maximum number of iterations to perform
     * in the Newton-Raphson iteration.
     */
    void setMaxIterations(double maxIterations_);

    /**
     * \brief Sets the number of radial basis functions.
     * \param nRadial_ The number of radial basis functions to set.
     */
    void setNRadial(int nRadial_);

    /**
     * \brief Sets the number of angular basis functions.
     * \param nTheta_ The number of angular basis functions to set.
     */
    void setNTheta(int nTheta_);

    /**
     * \brief Sets the maximum radius.
     * \param rMax_ The maximum radius to set.
     */
    void setMaximumRadius(double rMax_);

    /**
     * \brief Creates a Fields object with the prescribed parameters.
     * \todo Pass a pointer to a data object.
     * \param bareMass The bareMass of the Hamiltonian.
     * \param momentum The momentum of the ExtrinsicCurvature.
     * \param spin The spin of the ExtrinsicCurvature.
     * \param trumpetMass The mass term of the ExtrinsicCurvature.
     */
    Fields* createFields(double bareMass, double momentum, 
        double spin, double trumpetMass);

    /**
     * \brief Creates a Fields object with the prescribed parameters.
     * \todo Pass a pointer to a data object.
     * \param bareMass The bareMass of the Hamiltonian.
     * \param momentum The momentum of the ExtrinsicCurvature.
     * \param spin The spin of the ExtrinsicCurvature.
     * \param trumpetMass The mass of the ExtrinsicCurvature.
     * \param quadZ The z derivative of a curl.
     * \param quadPhi The curl of a curl quadrupole moment.
     */
    Fields* createFields(double bareMass, double momentum, 
        double spin, double trumpetMass, double quadZ, double quadPhi);
  private:

    /**
     * \brief The spherical basis used to solve the Hamiltonians.
     */
    SphericalBasis* basis;

    /**
     * \brief The magnitude of the residual to tolerate in
     * the Newton-Raphson iteration.
     */
    double tolerance;

    /**
     * \brief The maximum number of iterations to perform
     * in the Newton-Raphson iteration.
     */
    double maxIterations;

    /**
     * \brief The number of basis functions to use in the radial
     * direction.
     */
    int nRadial;

    /**
     * \brief The number of basis functions to use in the
     * zenithal direction.
     */
    int nTheta;

    /**
     * \brief The maximum radius.
     */
    double rMax;

    /**
     * Controls verbosity of some subroutines.
     */
    bool verbose;

    /**
     * Controls the severity of the Hamiltonian.
     */
    bool severe;

};

#endif
