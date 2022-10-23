#include "Particle.hpp"

Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity)
{
    q = charge;
    m = mass;
    pos = position;
    vel = velocity;
}

double Particle::charge()
{
    return q;
}

double Particle::mass()
    {
        return m;
    }

arma::vec Particle::position()
{
    return pos;
}

arma::vec Particle::velocity()
{
    return vel;
}