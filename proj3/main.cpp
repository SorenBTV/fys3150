#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "Particle.hpp"
#include "PenningTrap.hpp"


void one_particle_trap(PenningTrap trap, Particle p, int N_steps, double time_end);
void two_particle_trap(PenningTrap trap, Particle p1, Particle p2, bool particle_interactions, double dt, double time_end);

int main()
{   
    //Particle 1
    double q = 1;
    double m = 39.96;
    arma::vec r1 = {20, 0, 20};
    arma::vec v1 = {0, 25, 0};
    Particle p1(q, m, r1, v1);

    //Particle 2
    arma::vec r2 = {25, 25, 0};
    arma::vec v2 = {0, 40, 5};
    Particle p2(q, m, r1, v1);

    //Constructing trap
    double B_0 = 9.65 * 10;
    double V_0 = 2.41e6;
    double d = 500;
    PenningTrap trap(B_0, V_0, d);

    two_particle_trap(trap, p1, p2, true, 0.01, 50);

    return 0;
}


void two_particle_trap(PenningTrap trap, Particle p1, Particle p2, bool particle_interactions, double dt, double time_end)
{
    trap.add_particle(p1);
    trap.add_particle(p2);
    trap.particle_interaction(particle_interactions);

    int N_steps = time_end/dt;

    //position data (p1)
    arma::vec x_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec y_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec z_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    //velocity data (p1)
    arma::vec v_x_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec v_z_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    
    //position data (p2)
    arma::vec x_vals_p2 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec y_vals_p2 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec z_vals_p2 = arma::vec(N_steps, arma::fill::zeros);
    //velocity data (p2)
    arma::vec v_x_vals_p2 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec v_z_vals_p2 = arma::vec(N_steps, arma::fill::zeros);

    // initial values
    x_vals_p1(0) = p1.position()(0);
    y_vals_p1(0) = p1.position()(1);
    z_vals_p1(0) = p1.position()(2);

    x_vals_p2(0) = p2.position()(0);
    y_vals_p2(0) = p2.position()(1);
    z_vals_p2(0) = p2.position()(2);

    v_x_vals_p1(0) = p1.velocity()(0);
    v_z_vals_p1(0) = p1.velocity()(2);
    v_x_vals_p2(0) = p2.velocity()(0);
    v_z_vals_p2(0) = p2.velocity()(2);
    
    // compute position with Runge-Kutta 4th order
    for(int i = 1; i < N_steps; i++)
    {
        trap.evolve_RK4(dt);

        // p1 values
        arma::vec r = trap.particles[0].position();
        arma::vec v = trap.particles[0].velocity();
        
        x_vals_p1(i) = r(0);
        y_vals_p1(i) = r(1);
        z_vals_p1(i) = r(2);

        v_x_vals_p1(i) = v(0);
        v_z_vals_p1(i) = v(2);

        // p2 values
        r = trap.particles[1].position();
        v = trap.particles[1].velocity();

        x_vals_p2(i) = r(0);
        y_vals_p2(i) = r(1);
        z_vals_p2(i) = r(2);
        
        v_x_vals_p2(i) = v(0);
        v_z_vals_p2(i) = v(2);

        std::cout << r << "\n";
    }

}