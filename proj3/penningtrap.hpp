#ifndef __PenningTrap_hpp__  
#define __PenningTrap_hpp__

#include <string>
#include <vector>
#include <armadillo>
#include <cmath>
#include "Particle.hpp"

class PenningTrap
{

    public:

    double B_0;
    double V_0;
    double d;
    double ke;
    double freq;
    double amp;
    bool particle_interactions = true;

    std::vector<Particle> particles;


    // Constructor
    PenningTrap(double B0, double V0, double d_in);

    // Add a particle to the trap
    void add_particle(Particle p_in);

    //Particle ineractions true or false
    void particle_interaction(bool particle_interactions);

    //Time dependant V
    double V_time(double time);

    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r);

    // External electric field at point r=(x,y,z) and time (s)
    arma::vec external_E_field(arma::vec r, double time);

    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r);

    // Force on particle_i from particle_j
    arma::vec force_particle(int i, int j);

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i);

    // The total force on particle_i from the external fields with time dependant electric field
    arma::vec total_force_external(int i, double time);

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i);

    // The total force on particle_i from both external fields and other particles and time dependant electric field
    arma::vec total_force(int i, double time);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order with time dependant electric field
    void evolve_RK4(double dt, double time);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt);
};


#endif