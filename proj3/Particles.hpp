#include <armadillo>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <math.h>

class Particle
{
  public:
    double charge;
    double mass;
    arma::vec position;
    arma::vec velocity;

    Particle(double c, double m, arma::vec pos, arma::vec vel){
      charge = c;
      mass = m;
      position = pos;
      velocity = vel;
    }
};
