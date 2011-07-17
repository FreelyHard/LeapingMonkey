#include <iostream>

#include "FieldMaker.h"
#include "Fields.h"

FieldMaker::FieldMaker(double tolerance_, double maxIterations_, 
    int nRadial_, int nTheta_, double rMax_):
      tolerance(tolerance_), maxIterations(maxIterations_),
      nRadial(nRadial_), nTheta(nTheta_), rMax(rMax_), verbose(false),
      severe(true) 
{
  basis = new SphericalBasis();
  basis->setMaximumRadius(rMax);
};

FieldMaker::~FieldMaker() {
  delete basis;
}

void FieldMaker::setSevere(bool severity) {
  severe = severity;
}

void FieldMaker::setVerbosity(bool verbosity) {
  verbose = verbosity;
}

void FieldMaker::setTolerance(double tolerance_) {
  tolerance = tolerance_;
}

void FieldMaker::setMaxIterations(double maxIterations_) {
  maxIterations = maxIterations_;
}

void FieldMaker::setNRadial(int nRadial_) {
  nRadial = nRadial_;
}

void FieldMaker::setNTheta(int nTheta_) {
  nTheta = nTheta_;
}

void FieldMaker::setMaximumRadius(double rMax_) {
  rMax = rMax_;
}

Fields* FieldMaker::createFields(double bareMass, double momentum,
    double spin, double trumpetMass) {
  return createFields(bareMass, momentum, spin, trumpetMass, 0.0, 0.0);
}

Fields* FieldMaker::createFields(double bareMass, double momentum,
    double spin, double trumpetMass, double quadZ, double quadPhi) {
  //TODO: store kij and ham so that we aren't always creating and
  //destroying
  ExtrinsicCurvature kij;
  kij.setSpin(spin);
  kij.setMomentum(momentum);
  kij.setMass(trumpetMass);
  kij.setSpinCurl(quadPhi);
  kij.setQuadrupole(quadZ);

  Hamiltonian ham(kij, basis);
  ham.setSevere(severe);
  ham.setBareMass(bareMass);
  ham.setMaximumRadius(rMax);
  ham.setVerbose(verbose);
  ham.setRanks(nRadial, nTheta);

  ham.solve(tolerance, maxIterations);

  return new Fields(ham, kij);
}
