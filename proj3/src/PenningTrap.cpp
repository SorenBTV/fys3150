#include <vector>      // For vector
#include <string>      // For string
#include <stdlib.h>    // For rand (from C). For more powerful random number generation in C++ we should use <random>
#include <stdexcept>   // For runtime_error
#include <cmath>

#include "PenningTrap.hpp"
#include "Particle.hpp"



    // Constructor
    PenningTrap::PenningTrap(double B0, double V0, double d_in)
    {
        B_0 = B0;
        V_0 = V0;
        d = d_in;
    }

    // Add a particle to the trap
    void PenningTrap::add_particle(Particle p_in)
    {
        PenningTrap::particles.push_back(p_in);
    }

    //Particle interactions
    void PenningTrap::particle_interaction(bool particle_interactions)
    {
        PenningTrap::particle_interactions = particle_interactions;
    }


    double PenningTrap::V_time(double time)
    {
        return V_0*(1*amp*cos(freq*time));
    }


    // External electric field at point r=(x,y,z)
    arma::vec PenningTrap::external_E_field(arma::vec r)
    {
        arma::vec E(3);

        if(arma::norm(r)>d)
        {
            E = {0,0,0};
        }
        else
        {
            E(0) = -2*r(0);
            E(1) = -2*r(1);
            E(2) = 4*r(2);
        }
        return -V_0/(2*pow(d,2)) * E;
    }

    // External electric field at point r=(x,y,z) and time(s)
    arma::vec PenningTrap::external_E_field(arma::vec r, double time)
    {
        arma::vec E(3);
        double V_t = V_time(time);

        if(arma::norm(r)>d)
        {
            E = {0,0,0};
        }
        else
        {
            E(0) = -2*r(0);
            E(1) = -2*r(1);
            E(2) = 4*r(2);
        }
        return -V_t/(2*pow(d,2)) * E;
    }

    // External magnetic field at point r=(x,y,z)
    arma::vec PenningTrap::external_B_field(arma::vec r)
    {
        arma::vec B(3);
        if(arma::norm(r)<d)
        {
            B = {0,0,0};
        }
        else
        {
            B(0) = 0;
            B(1) = 0;
            B(2) = B_0;
        }
        return B;
    }

    // Force on particle_i from particle_j
    arma::vec PenningTrap::force_particle(int i, int j)
    {
        double ke = 1.38935333e5;
        double q_i = particles[i].charge();
        arma::vec r_i = particles[i].position();
        double q_j = particles[j].charge();
        arma::vec r_j = particles[j].position();

        double konst = ke * q_i * q_j;
        arma::vec F = konst * (r_i - r_j)/pow(arma::norm(r_i-r_j),3);

        return F;
    }

    // The total force on particle_i from the external fields
    arma::vec PenningTrap::total_force_external(int i)
    {
        arma::vec r = particles[i].position();
        arma::vec v = particles[i].velocity();
        double q = particles[i].charge();
        
        arma::vec B = external_B_field(r);
        arma::vec E = external_E_field(r);

        arma::vec F(3);
        F(0) = q * E(0) + q * (v(1)*B(2)-v(2)*B(1));
        F(1) = q * E(1) + q * (v(0)*B(2)-v(2)*B(0));
        F(2) = q * E(2) + q * (v(0)*B(1)-v(1)*B(0));
    
        return F;
    }

    // The total force on particle_i from the external fields and time dependant electric field
    arma::vec PenningTrap::total_force_external(int i, double time)
    {
        arma::vec r = particles[i].position();
        arma::vec v = particles[i].velocity();
        double q = particles[i].charge();
        
        arma::vec B = external_B_field(r);
        arma::vec E = external_E_field(r, time);

        arma::vec F(3);
        F(0) = q * E(0) + q * (v(1)*B(2)-v(2)*B(1));
        F(1) = q * E(1) + q * (v(0)*B(2)-v(2)*B(0));
        F(2) = q * E(2) + q * (v(0)*B(1)-v(1)*B(0));
    
        return F;
    }

    // The total force on particle_i from the other particles
    arma::vec PenningTrap::total_force_particles(int i)
    {
        arma::vec F(3);
        F(0) = 0;
        F(1) = 0;
        F(2) = 0;

        for (int j=0; j<particles.size(); j++)
        {
            if(i != j)
            {
                F(0) += force_particle(i,j)(0);
                F(1) += force_particle(i,j)(1);
                F(2) += force_particle(i,j)(2);
            }
        }
        return F;
    }

    // The total force on particle_i from both external fields and other particles
    arma::vec PenningTrap::total_force(int i)
    {
        arma::vec F(3);
        if(particle_interactions)
        {
            F = total_force_external(i) + total_force_particles(i);
        }
        else
        {
            F = total_force_external(i);
        }
        

        return F;
    }

    // The total force on particle_i from both external fields and other particles and time dependant electric field
    arma::vec PenningTrap::total_force(int i, double time)
    {
        arma::vec F(3);
        if(particle_interactions)
        {
            F = total_force_external(i, time) + total_force_particles(i);
        }
        else
        {
            F = total_force_external(i, time);
        }
        

        return F;
    }

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void PenningTrap::evolve_RK4(double dt)
    {
        
        //Temporary copy of position and velocity of particles in the trap
        std::vector<Particle> p_tmp = particles;

        int N = p_tmp.size();

        std::vector<arma::vec> k1_r(N);
        std::vector<arma::vec> k1_v(N);

        for (int i=0; i<N; i++)
        {
            arma::vec v = particles[i].velocity();
            double m = particles[i].mass();
            arma::vec F = total_force(i);
            k1_r[i] = dt * v;
            k1_v[i] = dt * F/m;
        }

        for (int i=0; i<N; i++)
        {
            arma::vec r = p_tmp[i].position();
            arma::vec v = p_tmp[i].velocity();
            particles[i].new_position(r + k1_r[i]/2);
            particles[i].new_velocity(v + k1_v[i]/2);
        }

        std::vector<arma::vec> k2_r(N);
        std::vector<arma::vec> k2_v(N);

        for (int i=0; i<N; i++)
        {
            arma::vec v = particles[i].velocity();
            double m = particles[i].mass();
            arma::vec F = total_force(i);
            k2_r[i] = dt * v;
            k2_v[i] = dt * F/m;
        }
        for (int i=0; i<N; i++)
        {
            arma::vec r = p_tmp[i].position();
            arma::vec v = p_tmp[i].velocity();
            particles[i].new_position(r + k2_r[i]/2);
            particles[i].new_velocity(v + k2_v[i]/2);
        }


        std::vector<arma::vec> k3_r(N);
        std::vector<arma::vec> k3_v(N);

        for (int i=0; i<N; i++)
        {
            arma::vec v = particles[i].velocity();
            double m = particles[i].mass();
            arma::vec F = total_force(i);
            k3_r[i] = dt * v;
            k3_v[i] = dt * F/m;
        }
        for (int i=0; i<N; i++)
        {
            arma::vec r = p_tmp[i].position();
            arma::vec v = p_tmp[i].velocity();
            particles[i].new_position(r + k3_r[i]/2);
            particles[i].new_velocity(v + k3_v[i]/2);
        }


        std::vector<arma::vec> k4_r(N);
        std::vector<arma::vec> k4_v(N);
        for (int i=0; i<N; i++)
        {
            arma::vec v = particles[i].velocity();
            double m = particles[i].mass();
            arma::vec F = total_force(i);
            k4_r[i] = dt * v;
            k4_v[i] = dt * F/m;
        }

        for (int i=0; i<N; i++)
        {
            arma::vec r = p_tmp[i].position();
            arma::vec v = p_tmp[i].velocity();

            arma::vec new_r = r +  1/6 * (k1_r[i] + 2*k2_r[i] + 2*k3_r[i] + k4_r[i]);
            arma::vec new_v = v +  1/6 * (k1_v[i] + 2*k2_v[i] + 2*k3_v[i] + k4_v[i]);

            particles[i].new_position(new_r);
            particles[i].new_velocity(new_v);
        }
    }


    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void PenningTrap::evolve_RK4(double dt, double time)
    {
        
        //Temporary copy of position and velocity of particles in the trap
        std::vector<Particle> p_tmp = particles;

        int N = p_tmp.size();

        std::vector<arma::vec> k1_r(N);
        std::vector<arma::vec> k1_v(N);

        for (int i=0; i<N; i++)
        {
            arma::vec v = particles[i].velocity();
            double m = particles[i].mass();
            arma::vec F = total_force(i, time);
            k1_r[i] = dt * v;
            k1_v[i] = dt * F/m;
        }

        for (int i=0; i<N; i++)
        {
            arma::vec r = p_tmp[i].position();
            arma::vec v = p_tmp[i].velocity();
            particles[i].new_position(r + k1_r[i]/2);
            particles[i].new_velocity(v + k1_v[i]/2);
        }

        std::vector<arma::vec> k2_r(N);
        std::vector<arma::vec> k2_v(N);

        for (int i=0; i<N; i++)
        {
            arma::vec v = particles[i].velocity();
            double m = particles[i].mass();
            arma::vec F = total_force(i, time + dt/2);
            k2_r[i] = dt * v;
            k2_v[i] = dt * F/m;
        }
        for (int i=0; i<N; i++)
        {
            arma::vec r = p_tmp[i].position();
            arma::vec v = p_tmp[i].velocity();
            particles[i].new_position(r + k2_r[i]/2);
            particles[i].new_velocity(v + k2_v[i]/2);
        }


        std::vector<arma::vec> k3_r(N);
        std::vector<arma::vec> k3_v(N);

        for (int i=0; i<N; i++)
        {
            arma::vec v = particles[i].velocity();
            double m = particles[i].mass();
            arma::vec F = total_force(i, time + dt/2);
            k3_r[i] = dt * v;
            k3_v[i] = dt * F/m;
        }
        for (int i=0; i<N; i++)
        {
            arma::vec r = p_tmp[i].position();
            arma::vec v = p_tmp[i].velocity();
            particles[i].new_position(r + k3_r[i]/2);
            particles[i].new_velocity(v + k3_v[i]/2);
        }


        std::vector<arma::vec> k4_r(N);
        std::vector<arma::vec> k4_v(N);
        for (int i=0; i<N; i++)
        {
            arma::vec v = particles[i].velocity();
            double m = particles[i].mass();
            arma::vec F = total_force(i, time + dt);
            k4_r[i] = dt * v;
            k4_v[i] = dt * F/m;
        }

        for (int i=0; i<N; i++)
        {
            arma::vec r = p_tmp[i].position();
            arma::vec v = p_tmp[i].velocity();

            arma::vec new_r = r +  1/6 * (k1_r[i] + 2*k2_r[i] + 2*k3_r[i] + k4_r[i]);
            arma::vec new_v = v +  1/6 * (k1_v[i] + 2*k2_v[i] + 2*k3_v[i] + k4_v[i]);

            particles[i].new_position(new_r);
            particles[i].new_velocity(new_v);
        }
    }

    // Evolve the system one time step (dt) using Forward Euler
    //void PenningTrap::evolve_forward_Euler(double dt);

