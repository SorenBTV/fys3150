#include <armadillo>
#include <iostream>


class Particle
{
    public:

    double q;
    double m;
    arma::vec pos;
    arma::vec vel;

    Particle(double charge, double mass, arma::vec position, arma::vec velocity)
    {
        q = charge;
        m = mass;
        pos = position;
        vel = velocity;
    }
};