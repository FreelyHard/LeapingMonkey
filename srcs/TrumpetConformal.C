#include <cmath>
#include <cassert>
#include "TrumpetConformal.h"

double abs(double x) {
  if (x < 0.0) return -x;
  return x;
}

double TrumpetConformal::psi(double r) {
  return sqrt(trumpetR(r)/r);
}

double TrumpetConformal::trumpetR(double r) {
  double guessR = r;
  if (r < 1.5) guessR = 1.51;
  double error = r - isotropicr(guessR);
  int count = 0;
  int countMax = 1000;
  while (abs(error) > 1.0e-12*r && count < countMax) {
    count++;
    double delta = error/drdR(guessR);
    guessR += delta;
    error = r - isotropicr(guessR);
  }
  assert(guessR == guessR); // Test for nan.
  assert(count < countMax);
  return guessR;
}

double TrumpetConformal::isotropicr(double R) {
  return 0.2500000000e0 * (0.2e1 * R + m + sqrt(0.4e1 * R * R + 0.4e1 * m * R + 0.3e1 * m * m)) * pow(0.8242640686e1 * (0.2e1 * R - 0.3e1 * m) / (0.8e1 * R + 0.6e1 * m + 0.3e1 * sqrt(0.8e1 * R * R + 0.8e1 * m * R + 0.6e1 * m * m)), sqrt(0.2e1) / 0.2e1);
}

double TrumpetConformal::drdR(double R) {
  return 0.2203665914e-13 * (0.672792206e9 * sqrt(0.8e1 * R * R + 0.8e1 * m * R + 0.6e1 * m * m) * sqrt(0.4e1 * R * R + 0.4e1 * m * R + 0.3e1 * m * m) * R * m + 0.1600000000e10 * sqrt(0.8e1 * R * R + 0.8e1 * m * R + 0.6e1 * m * m) * pow(R, 0.3e1) + 0.1009188309e10 * sqrt(0.8e1 * R * R + 0.8e1 * m * R + 0.6e1 * m * m) * pow(m, 0.3e1) + 0.6788225099e10 * pow(R, 0.3e1) * m + 0.7079393920e10 * R * R * m * m + 0.2982337650e10 * pow(m, 0.3e1) * R + 0.2400000000e10 * sqrt(0.4e1 * R * R + 0.4e1 * m * R + 0.3e1 * m * m) * pow(R, 0.3e1) - 0.1427207794e10 * sqrt(0.4e1 * R * R + 0.4e1 * m * R + 0.3e1 * m * m) * pow(m, 0.3e1) + 0.800000000e9 * sqrt(0.8e1 * R * R + 0.8e1 * m * R + 0.6e1 * m * m) * sqrt(0.4e1 * R * R + 0.4e1 * m * R + 0.3e1 * m * m) * R * R - 0.263603897e9 * sqrt(0.8e1 * R * R + 0.8e1 * m * R + 0.6e1 * m * m) * sqrt(0.4e1 * R * R + 0.4e1 * m * R + 0.3e1 * m * m) * m * m + 0.2145584412e10 * sqrt(0.8e1 * R * R + 0.8e1 * m * R + 0.6e1 * m * m) * m * R * R + 0.145584412e9 * sqrt(0.8e1 * R * R + 0.8e1 * m * R + 0.6e1 * m * m) * R * m * m + 0.2194112549e10 * sqrt(0.4e1 * R * R + 0.4e1 * m * R + 0.3e1 * m * m) * m * R * R + 0.2442640686e10 * sqrt(0.4e1 * R * R + 0.4e1 * m * R + 0.3e1 * m * m) * R * m * m + 0.4800000000e10 * pow(R, 0.4e1) + 0.1118376617e10 * pow(m, 0.4e1)) * pow((0.1648528137e10 * R - 0.2472792206e10 * m) / (0.8e1 * R + 0.6e1 * m + 0.3e1 * sqrt(0.8e1 * R * R + 0.8e1 * m * R + 0.6e1 * m * m)), 0.707106781e9 / 0.1000000000e10) / (0.2e1 * R - 0.3e1 * m) * pow(0.8e1 * R * R + 0.8e1 * m * R + 0.6e1 * m * m, -0.1e1 / 0.2e1) / (0.8e1 * R + 0.6e1 * m + 0.3e1 * sqrt(0.8e1 * R * R + 0.8e1 * m * R + 0.6e1 * m * m)) * pow(0.4e1 * R * R + 0.4e1 * m * R + 0.3e1 * m * m, -0.1e1 / 0.2e1);
}
