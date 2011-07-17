#include <cstdlib>
#include <iostream>
#include <cmath>
#include "TrumpetConformal.h"

int main() {
  double mass = 1.;
  double r = 0.3;
  double h = 0.00001;
  TrumpetConformal trumpet(mass);
  double psi0 = trumpet.psi(r);
  double psip = trumpet.psi(r+h);
  double psim = trumpet.psi(r-h);
  double lapPsi = (psip + psim - 2.*psi0)/(h*h) + 1./r*(psip - psim)/h;
  double source = 81./64.*pow(mass, 4.)*pow(r, -6.)*pow(psi0, -7);
  double zero = lapPsi + source;
  std::cout << " We have D^2 psi + A^2/8 psi^7 = " << zero << 
    " using finite differences.\n";
}
