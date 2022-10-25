#include "Particle.hpp"

Particle::Particle(double charge_in, double mass_in, arma::vec position_in, arma::vec velocity_in)
{
    q = charge_in;
    m = mass_in;
    pos = position_in;
    vel = velocity_in;
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

void Particle::new_position(arma::vec new_r)
{
    r = new_r;
}

void Particle::new_velocity(arma::vec new_v)
{
    v = new_v;
}