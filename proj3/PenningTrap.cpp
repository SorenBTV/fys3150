#include "Particle.hpp"
#include "PenningTrap.hpp"

#include <armadillo>
#include <iostream>

//Constructor
PenningTrap::PenningTrap(double B_0, double V_0, double d_in)
{
    B0 = B_0;
    V0 = V_0;
    d = d_in;
    ke = 1.38935333e5;
};

//Add a particle to the trap
void PenningTrap::add_Particle(Particle &p_in)
{
    particles.push_back(p_in);
};

// External electric field at point r=(x,y,z)
arma::vec external_E_field(arma::vec r)
{
    arma::vec V(3);
    v(0) = -2*r(0);
    V(1) = -2*r(1);
    V(2) = 4*r(2);
    V = V * (V0/(2*pow(d,2)));
    return V;
}

// External magnetic field at point r=(x,y,z)
arma::vec external_B_field(arma::vec r)
{
    arma::vec B(3);
    B(0) = 0;
    B(1) = 0;
    B(2) = B_0;
    return B;
}


// Force on particle_i from particle_j
arma::vec force_particle(int i, int j)
{
    arma::vec force;
    Particle p_i = p_in[i];
    Particle p_j = p_in[j];
    arma::vec pos_i = p_i.position;
    arma::vec pos_j = p_j.position;
    arma::vec dist = pos_i - pos_j;

    double abs_dist = abs(dist(0) + dist(1) + dist(2));
    force = ke * p_i.q * p_j.q * pow(abs_dist,-3) * dist;
    return force;
}

// The total force on particle_i from the external fields
arma::vec total_force_external(int i);

// The total force on particle_i from the other particles
arma::vec total_force_particles(int i);

// The total force on particle_i from both external fields and other particles
arma::vec total_force(int i);

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void evolve_RK4(double dt);

// Evolve the system one time step (dt) using Forward Euler
void evolve_forward_Euler(double dt);