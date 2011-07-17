#include <cmath>
#include "ExtrinsicCurvature.h"

void ExtrinsicCurvature::setSpinCurl(double cq) {
  cQ = cq;
}

double ExtrinsicCurvature::getSpinCurl() {
  return cQ;
}

void ExtrinsicCurvature::setQuadrupole(double q) {
  Q = q;
}

double ExtrinsicCurvature::getQuadrupole() {
  return Q;
}

void ExtrinsicCurvature::setSpin(double S) {
  sz = S;
}

double ExtrinsicCurvature::getSpin() {
  return sz;
}

void ExtrinsicCurvature::setMomentum(double P) {
  pz = P;
}

double ExtrinsicCurvature::getMomentum() {
  return pz;
}

void ExtrinsicCurvature::setMass(double m) {
  if (m < 0) {
    xi = -1.0;
  } else {
    xi = 1.0;
  }
  M = m;
}

double ExtrinsicCurvature::getMass() {
  return M;
}

double ExtrinsicCurvature::kappa(double R, double theta) {
  if (Q != 0.0 || cQ != 0.0) {
    return 0.288e3 * cQ * cQ - 0.288e3 * sz * R * pow(cos(theta), 0.3e1) * Q + 0.576e3 * cQ * cQ * pow(cos(theta), 0.2e1) + 0.1008e4 * Q * Q * pow(cos(theta), 0.2e1) - 0.1080e4 * Q * Q * pow(cos(theta), 0.4e1) + 0.18e2 * sz * sz * R * R + 0.9e1 / 0.2e1 * pz * pz * pow(R, 0.4e1) + 0.81e2 / 0.8e1 * xi * xi * pow(M, 0.4e1) * R * R + 0.72e2 * pz * R * R * cQ - 0.18e2 * R * R * sz * sz * pow(cos(theta), 0.2e1) + 0.9e1 * pow(R, 0.4e1) * pz * pz * pow(cos(theta), 0.2e1) + 0.72e2 * Q * Q - 0.27e2 / 0.2e1 * pow(R, 0.3e1) * pz * xi * sqrt(0.3e1) * M * M * cos(theta) + 0.288e3 * R * sz * cos(theta) * Q - 0.288e3 * pz * R * R * pow(cos(theta), 0.2e1) * cQ + 0.108e3 * xi * sqrt(0.3e1) * M * M * R * cQ * cos(theta);
  } else if (sz != 0.0 || M != 0.0) {
    return 0.9e1 * R * R * pz * pz * pow(cos(theta), 0.2e1) - 0.18e2 * sz * sz * pow(cos(theta), 0.2e1) + 0.18e2 * sz * sz + 0.9e1 / 0.2e1 * pz * pz * R * R + 0.81e2 / 0.8e1 * pow(M, 0.4e1) * xi * xi - 0.27e2 / 0.2e1 * R * pz * xi * sqrt(0.3e1) * M * M * cos(theta);
  } else if (pz != 0.0) {
    return 0.9e1 / 0.2e1 * pz * pz * (0.2e1 * pow(cos(theta), 0.2e1) + 0.1e1);
  }
  return 0.0;
}

int ExtrinsicCurvature::getSingularityPower() const {
  if (Q != 0.0 || cQ != 0.0) {
    return 8;
  } else if (sz != 0.0 || M != 0.0) {
    return 6;
  } else if (pz != 0.0) {
    return 4;
  }
  return 0;
}

double* ExtrinsicCurvature::getBeta(const double* theta_, int nTheta) const {
  double* beta = new double[nTheta];
  for (int i = 0; i < nTheta; i++) {
    double theta = theta_[i];
    if (Q != 0.0 || cQ != 0.0) {
      beta[i] = 0.288e3 * cQ * cQ + 0.576e3 * cQ * cQ * pow(cos(theta), 0.2e1) + 0.1008e4 * Q * Q * pow(cos(theta), 0.2e1) - 0.1080e4 * Q * Q * pow(cos(theta), 0.4e1) + 0.72e2 * Q * Q;
    } else if (sz != 0.0 || M != 0.0) {
      beta[i] = -0.18e2 * sz * sz * pow(cos(theta), 0.2e1) + 0.18e2 * sz * sz + 0.81e2 / 0.8e1 * pow(M, 0.4e1) * xi * xi;
    } else if (pz != 0.0) {
      beta[i] = 0.9e1 / 0.2e1 * pz * pz * (0.2e1 * pow(cos(theta), 0.2e1) + 0.1e1);
    } else {
      beta[i] = 0.0;
    }
  }
  return beta;
}

int ExtrinsicCurvature::indexAij(int i, int j) {
  //TODO: throw exception if j is faulty...
  if (i > j) {
    return indexAij(j,i);
  } else if (i == 1) {
    return j-1;
  } else if (i == 2) {
    return 3 + j-2;
  } else {
    return 5;
  }
}

double ExtrinsicCurvature::Axx(double R, double theta) {
  return 0.3e1 / 0.4e1 * (-0.2e1 * R * R * pow(cos(theta), 0.3e1) * pz - 0.2e1 * R * xi * sqrt(0.3e1) * M * M + 0.3e1 * xi * sqrt(0.3e1) * M * M * R * pow(cos(theta), 0.2e1) - 0.64e2 * cQ * cos(theta) + 0.80e2 * cQ * pow(cos(theta), 0.3e1)) * pow(R, -0.4e1);
}

double ExtrinsicCurvature::Axy(double R, double theta) {
  return -0.3e1 * pow(R, -0.4e1) * pow(sin(theta), 0.2e1) * (R * sz + 0.10e2 * Q * cos(theta));
}

double ExtrinsicCurvature::Axz(double R, double theta) {
  return -0.3e1 / 0.4e1 * sin(theta) * (-0.2e1 * R * R * pz - 0.2e1 * R * R * pow(cos(theta), 0.2e1) * pz + 0.3e1 * R * xi * sqrt(0.3e1) * M * M * cos(theta) - 0.16e2 * cQ + 0.80e2 * cQ * pow(cos(theta), 0.2e1)) * pow(R, -0.4e1);
}

double ExtrinsicCurvature::Ayz(double R, double theta) {
  return -0.3e1 * pow(R, -0.4e1) * sin(theta) * (R * sz * cos(theta) + 0.10e2 * Q * pow(cos(theta), 0.2e1) - 0.2e1 * Q);
}

double ExtrinsicCurvature::Azz(double R, double theta) {
  return -0.3e1 / 0.4e1 * (-0.2e1 * pz * R * R * cos(theta) - 0.2e1 * R * R * pow(cos(theta), 0.3e1) * pz + 0.3e1 * xi * sqrt(0.3e1) * M * M * R * pow(cos(theta), 0.2e1) - R * xi * sqrt(0.3e1) * M * M + 0.80e2 * cQ * pow(cos(theta), 0.3e1) - 0.48e2 * cQ * cos(theta)) * pow(R, -0.4e1);
}

double ExtrinsicCurvature::Axxr(double R, double theta) {
  return 0.3e1 / 0.4e1 * (-0.4e1 * R * pow(cos(theta), 0.3e1) * pz - 0.2e1 * xi * sqrt(0.3e1) * M * M + 0.3e1 * xi * sqrt(0.3e1) * M * M * pow(cos(theta), 0.2e1)) * pow(R, -0.4e1) - 0.3e1 * (-0.2e1 * R * R * pow(cos(theta), 0.3e1) * pz - 0.2e1 * R * xi * sqrt(0.3e1) * M * M + 0.3e1 * xi * sqrt(0.3e1) * M * M * R * pow(cos(theta), 0.2e1) - 0.64e2 * cQ * cos(theta) + 0.80e2 * cQ * pow(cos(theta), 0.3e1)) * pow(R, -0.5e1);
}

double ExtrinsicCurvature::Axzr(double R, double theta) {
  return -0.3e1 / 0.4e1 * sin(theta) * (-0.4e1 * pz * R - 0.4e1 * R * pow(cos(theta), 0.2e1) * pz + 0.3e1 * xi * sqrt(0.3e1) * M * M * cos(theta)) * pow(R, -0.4e1) + 0.3e1 * sin(theta) * (-0.2e1 * R * R * pz - 0.2e1 * R * R * pow(cos(theta), 0.2e1) * pz + 0.3e1 * R * xi * sqrt(0.3e1) * M * M * cos(theta) - 0.16e2 * cQ + 0.80e2 * cQ * pow(cos(theta), 0.2e1)) * pow(R, -0.5e1);
}

double ExtrinsicCurvature::Azzr(double R, double theta) {
  return -0.3e1 / 0.4e1 * (-0.4e1 * pz * R * cos(theta) - 0.4e1 * R * pow(cos(theta), 0.3e1) * pz + 0.3e1 * xi * sqrt(0.3e1) * M * M * pow(cos(theta), 0.2e1) - xi * sqrt(0.3e1) * M * M) * pow(R, -0.4e1) + 0.3e1 * (-0.2e1 * pz * R * R * cos(theta) - 0.2e1 * R * R * pow(cos(theta), 0.3e1) * pz + 0.3e1 * xi * sqrt(0.3e1) * M * M * R * pow(cos(theta), 0.2e1) - R * xi * sqrt(0.3e1) * M * M + 0.80e2 * cQ * pow(cos(theta), 0.3e1) - 0.48e2 * cQ * cos(theta)) * pow(R, -0.5e1);
}
