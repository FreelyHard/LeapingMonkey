#include <iostream>
#include <sstream>
#include <cassert>
#include <string>
#include <fstream>
#include <cmath>

#include "FieldMaker.h"
#include "ChebyshevRoots.h"
#include "SphericalHorizon.h"
#include "TrumpetConformal.h"

// Sadly this code is highly imperative... but I am unsure how to
// abstract it. Plus it just needs to produce data for a paper.
int main() {

  //Select which data to compute...
  bool computePunctureSpins = false;
  bool testTrumpetSolution = false;
  bool computeTrumpetSpins = false;
  bool singularSpin = false;
  bool computeQuadrupoles = true;

  // Numerical parameters.
  int nTheta = 40;
  int nRadial = 130;
  int rMax = 10;
  double tolerance = 1.0e-8;
  int maxIterations = 15;
  int hIterations = 20;
  bool verbose = true;
  bool severe = false;

  // Base physical parameters.
  double bareMass = 1.;
  double trumpetMass = 0.;
  double momentum = 0.;
  double spin = 0.0015;

  // Family parameters.
  double dlogJ = 0.025;
  double logJMin = -2.0;
  double logJMax = 2. + dlogJ*0.5;

  std::cout.precision(16);

  ChebyshevRoots* aBasis = new ChebyshevRoots();
  aBasis->setRank(nTheta);

  FieldMaker solver(tolerance, maxIterations, nRadial, nTheta, rMax);
  solver.setVerbosity(verbose);
  solver.setSevere(severe);

  // Start building solutions...
  double chiInner = -1.0;
  double chiOuter = 1.0;
  double* outer = new double[nTheta];
  for (int i = 0; i < nTheta; i++) {
    outer[i] = 0.5;
  }

  if (computePunctureSpins) {
    std::cout << "Computing spinning puncture solutions.\n";
    std::stringstream filename;
    filename << "puncture." << nTheta << "." << nRadial << ".spins";
    std::ofstream myFile;
    myFile.open(filename.str().c_str());
    myFile << "# M_p M_ADM M_t J epsilon M_irr M zeta chi\n";
    myFile.precision(16);
    for (double logJ = logJMin; logJ < logJMax; logJ += dlogJ) {
      double J = pow(10.0, logJ);
      // Scale the tolerance with J.
      solver.setTolerance(tolerance*J);
      if (verbose) std::cout << "Puncture: J = " << J << "\n";
      Fields* fields = 
        solver.createFields(bareMass, momentum, J, trumpetMass);
      double mADM = fields->computeADMMass();
      SphericalHorizon* horizonFinder = new SphericalHorizon(aBasis, fields);
      horizonFinder->computeSolution(outer, &chiOuter,
          tolerance, hIterations);

      double mIrr = horizonFinder->computeHorizonMass(outer);
      double spinHO = horizonFinder->computeHorizonSpin(outer);
      double mass = sqrt(mIrr*mIrr + 0.25*spinHO*spinHO*pow(mIrr, -2));
      myFile << bareMass << " " << mADM << " "
        << trumpetMass << " " << J << " " << J/(mADM*mADM) << " " <<
        mIrr << " " << mass << " "<< 0.5*spinHO/(mIrr*mIrr) << " " << 
        spinHO/(mass*mass) << "\n";
      myFile.flush();

      delete fields;
      delete horizonFinder;
    }
    myFile.close();
  } // End loop where puncture spins are created.

  // Test that we can compute a trumpet and that the error is small.
  if (testTrumpetSolution) {
    TrumpetConformal trumpet(1.0);
    std::cout << "Computing normal trumpet solution.\n";

    Fields* fields = solver.createFields(0.0, 0.0, 0.0, 1.0);

    double numericalValue = fields->psi(1.0, 1.0);
    double actualValue = trumpet.psi(1.0);
    double r = rMax*0.5;
    double bound = fields->psir(r, 1.0) + (fields->psi(r, 1.0)-1.)/r;

    assert( -(numericalValue - actualValue) < 1.6e-4);
    assert(bound < 8.0e-5);

    delete fields;

  }

  solver.setTolerance(tolerance);
  if (computeTrumpetSpins) {
    for (int i = 0; i < nTheta; i++) {
      outer[i] = 0.5;
    }
    std::cout << "Computing spinning trumpet solutions.\n";
    bareMass = 0.0;
    trumpetMass = 1.0;
    std::stringstream filename;
    filename << "trumpet." << nTheta << "." << nRadial << ".spins";
    std::ofstream myFile;
    myFile.open(filename.str().c_str());
    myFile.precision(16);
    myFile << "# M_p M_ADM M_t J epsilon M_irr M zeta chi\n";
    for (double logJ = logJMin; logJ < logJMax; logJ += dlogJ) {
      double J = pow(10.0, logJ);
      if (verbose) std::cout << "Trumpet: J = " << J << "\n";
      // Scale mass so mADM ~ 1.
      Fields* fields =
        solver.createFields(bareMass, momentum, J, trumpetMass);
      double mADM = fields->computeADMMass();
      trumpetMass = trumpetMass/mADM;
      delete fields;
      fields = solver.createFields(bareMass, momentum, J, trumpetMass);
      mADM = fields->computeADMMass();

      SphericalHorizon* horizonFinder = new SphericalHorizon(aBasis, fields);
      horizonFinder->computeSolution(outer, &chiOuter,
            tolerance, hIterations);

      double mIrr = horizonFinder->computeHorizonMass(outer);
      double spinHO = horizonFinder->computeHorizonSpin(outer);
      double mass = sqrt(mIrr*mIrr + 0.25*spinHO*spinHO*pow(mIrr, -2));
      myFile << bareMass << " " << mADM << " "
        << trumpetMass << " " << J << " " << J/(mADM*mADM) << " " <<
        mIrr << " " << mass << " "<< 0.5*spinHO/(mIrr*mIrr) << " " << 
        spinHO/(mass*mass) << "\n";
      myFile.flush();

      delete fields;
      delete horizonFinder; 
    }
    bareMass = 1.0;
    trumpetMass = 0.0;
    myFile.close();
  }

  if (singularSpin) {
    std::cout << "Creating singular spin source.\n";
    Fields* fields = solver.createFields(0., 0., 1., 0.);
    for (int i = 0; i < nTheta; i++) {
      outer[i] = 0.5;
    }
    SphericalHorizon* horizonFinder = new SphericalHorizon(aBasis, fields);
    double spinHO = horizonFinder->computeHorizonSpin(outer);
    double mIrr = horizonFinder->computeHorizonMass(outer);
    double mass = sqrt(mIrr*mIrr + 0.25*spinHO*spinHO*pow(mIrr, -2));
    double mADM = fields->computeADMMass();
    std::cout << mADM << " " << 1./(mADM*mADM) << " " << spinHO << " " 
      << mIrr << " " << mass << " " << spinHO/(mass*mass) << "\n";
  }

  if (computeQuadrupoles) {
    std::cout << "Creating quadrupole sources.\n";
    bareMass = 1.0;
    trumpetMass = 0.0;
    double J = 0.5;
    for (int i = 0; i < nTheta; i++) {
      outer[i] = 0.5;
    }
    // Set an initial guess for the apparent horizon...
    Fields* fields = 
        solver.createFields(bareMass, momentum, J, trumpetMass);
    SphericalHorizon* horizonFinder = new SphericalHorizon(aBasis, fields);
    horizonFinder->computeSolution(outer, &chiOuter, tolerance, hIterations);
    delete fields;
    delete horizonFinder;

    std::stringstream filename;
    filename << "quadrupole.1." << nTheta << "." << nRadial << ".spins";
    std::ofstream myFile;
    myFile.open(filename.str().c_str());
    myFile.precision(16);
    myFile << "# q M_p M_ADM M_t J epsilon M_irr M zeta chi\n";
    double dq = 0.01;
    double qMax = 0.5;
    double qMin = dq;
    for (double q = qMin; q < qMax; q += dq) {
      fields = solver.createFields(bareMass, momentum, J, trumpetMass, q, 0.0);
      horizonFinder = new SphericalHorizon(aBasis, fields);
      horizonFinder->computeSolution(outer, &chiOuter, tolerance, hIterations);

      double mADM = fields->computeADMMass();
      double mIrr = horizonFinder->computeHorizonMass(outer);
      double spinHO = horizonFinder->computeHorizonSpin(outer);
      double mass = sqrt(mIrr*mIrr + 0.25*spinHO*spinHO*pow(mIrr, -2));
      myFile << q << " " << bareMass << " " << mADM << " "
        << trumpetMass << " " << J << " " << J/(mADM*mADM) << " " <<
        mIrr << " " << mass << " "<< 0.5*spinHO/(mIrr*mIrr) << " " << 
        spinHO/(mass*mass) << "\n";
      myFile.flush();

      delete fields;
      delete horizonFinder; 
    }
    myFile.close();

    // The computations for the second type of quadrupole moment.
    std::stringstream filename2;
    filename2 << "quadrupole.2." << nTheta << "." << nRadial << ".spins";
    std::ofstream myFile2;
    myFile2.open(filename2.str().c_str());
    myFile2.precision(16);
    myFile2 << "# q M_p M_ADM M_t J epsilon M_irr M zeta chi\n";
    for (int i = 0; i < nTheta; i++) {
      outer[i] = 0.5;
    }
    for (double q = qMin; q < qMax; q += dq) {
      fields = solver.createFields(bareMass, momentum, J, trumpetMass, 0.0, q);
      horizonFinder = new SphericalHorizon(aBasis, fields);
      horizonFinder->computeSolution(outer, &chiOuter, tolerance, hIterations);

      double mADM = fields->computeADMMass();
      double mIrr = horizonFinder->computeHorizonMass(outer);
      double spinHO = horizonFinder->computeHorizonSpin(outer);
      double mass = sqrt(mIrr*mIrr + 0.25*spinHO*spinHO*pow(mIrr, -2));
      myFile2 << q << " " << bareMass << " " << mADM << " "
        << trumpetMass << " " << J << " " << J/(mADM*mADM) << " " <<
        mIrr << " " << mass << " "<< 0.5*spinHO/(mIrr*mIrr) << " " << 
        spinHO/(mass*mass) << "\n";
      myFile2.flush();

      delete fields;
      delete horizonFinder; 
    }
    myFile2.close();
  }

}
