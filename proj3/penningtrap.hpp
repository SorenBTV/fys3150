#include "Particle.hpp"

class PenningTrap
{

    public:

    double B_0;
    double V_0;
    double d_in;
    double ke;

    std::vector<Particle> particles;


    // Constructor
    PenningTrap(double B0, double V0, double d)
    {
        B_0 = B0;
        V_0 = V0;
        d_in = d;
    }

    // Add a particle to the trap
    void add_particle(Particle p_in)
    {
        particles.push_back(p_in);
    }


    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r)
    {
        arma::vec E(3);
        E(0) = 9.65*r(0);
        E(1) = 9.65*r(1);
        E(2) = -9.65*2*r(2);
        return E;
    }

    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field()
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
        double ke = 1.38935333e5;
        Particle p_i = particles[i];
        Particle p_j = particles[j];

        double pos_ix = p_i.pos(0);
        double pos_iy = p_i.pos(1);
        double pos_iz = p_i.pos(2);

        double pos_jx = p_j.pos(0);
        double pos_jy = p_j.pos(1);
        double pos_jz = p_j.pos(2);

        double konst = ke * p_i.q * p_j.q;
        arma::vec pos = {pos_ix-pos_jx, pos_iy-pos_jy, pos_iz-pos_jz};

        double param = sqrt(pow(pos(0),2) + pow(pos(2),2) + pow(pos(2),2));
        double parampow = pow(param, 3);

        arma::vec F(3);
        F(0) = konst * pos(0)/parampow;
        F(1) = konst * pos(1)/parampow;
        F(2) = konst * pos(2)/parampow;

        return F;
    }

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i)
    {

        Particle p_i = particles[i];
        arma::vec v = p_i.vel;
        arma::vec B = external_B_field();
        arma::vec E = external_E_field(p_i.pos);

        arma::vec F(3);
        F(0) = p_i.q * E(0) + p_i.q * (v(1)*B(2)-v(2)*B(1));
        F(1) = p_i.q * E(1) + p_i.q * (v(0)*B(2)-v(2)*B(0));
        F(2) = p_i.q * E(2) + p_i.q * (v(0)*B(1)-v(1)*B(0));
    
        return F;
    }

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i)
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
    arma::vec total_force(int i)
    {
        arma::vec F(3);
        F = total_force_external(i) + total_force_particles(i);

        return F;
    }

    arma::mat krv(int i, arma::mat kr, arma::mat kv, int step, double dt, arma::vec r, arma::vec v)
    {
        step = 0;
        for(int j=0; j<particles.size(); j++)
        {  
            double charge = particles[j].q;
            double mass = particles[j].m;
            arma::vec temp_r = particles[j].pos;
            arma::vec temp_v = particles[j].vel;
            arma::vec F = total_force(j);

            kr(0,j) = temp_v(0);
            kr(1,j) = temp_v(1);
            kr(2,j) = temp_v(2);

            kv(0,j) = F(0)/mass;
            kv(1,j) = F(1)/mass;
            kv(2,j) = F(2)/mass;

            temp_r = {r(i,step) + kr(0,j)*dt/2, r(i,step+1) + kr(1,j)*dt/2, r(i,step+2) + kr(2,j)*dt/2};
            temp_v = {v(i,step) + kv(0,j)*dt/2, v(i,step+1) + kv(1,j)*dt/2, v(i,step+2) + kv(2,j)*dt/2};
            particles[j] = {charge, mass, temp_r, temp_v};
            step += 3;
        }
    return kr, kv;

    }


    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt, int N)
    {
        int t_end = 50;
        arma::mat r(N, particles.size()*3);
        arma::mat v(N, particles.size()*3);

        int step = 0;
        for(int i=0; i<particles.size(); i++)
        {
            Particle p = particles[i];
            v(0, step) = p.vel(0);
            v(0, step+1) = p.vel(1);
            v(0, step+2) = p.vel(2);

            r(0, step) = p.pos(0);
            r(0, step+1) = p.pos(1);
            r(0, step+2) = p.pos(2);

        step += 3;
        }

        for(int i=0; i<N-1; i++)
        {
            arma::mat kr1(3, particles.size());
            arma::mat kv1(3, particles.size());

            arma::mat kr2(3, particles.size());
            arma::mat kv2(3, particles.size());

            arma::mat kr3(3, particles.size());
            arma::mat kv3(3, particles.size());

            arma::mat kr4(3, particles.size());
            arma::mat kv4(3, particles.size());


            kr1, kv1 = krv(i, kr1, kv1, 0, dt, r, v);
            kr2, kv2 = krv(i, kr2, kv2, 0, dt, r, v);
            kr3, kv3 = krv(i, kr3, kv3, 0, dt, r, v);
            kr4, kv4 = krv(i, kr4, kv4, 0, dt, r, v);

                    
            step = 0;
            for(int j=0; j<particles.size(); j++)
            {
                r(i+1,step) = r(i,step) + (dt / 6) * (kr1(0,j) + 2*kr2(0,j) + 2*kr3(0,j) + kr4(0,j));
                r(i+1,step+1) = r(i,step+1) + (dt / 6) * (kr1(1,j) + 2*kr2(1,j) + 2*kr3(1,j) + kr4(1,j));
                r(i+1,step+2) = r(i,step+2) + (dt / 6) * (kr1(2,j) + 2*kr2(2,j) + 2*kr3(2,j) + kr4(2,j));

                v(i+1,step) = v(i,step) + (dt / 6) * (kv1(0,j) + 2*kv2(0,j) + 2*kv3(0,j) + kv4(0,j));
                v(i+1,step+1) = v(i,step+1) + (dt / 6) * (kv1(1,j) + 2*kv2(1,j) + 2*kv3(1,j) + kv4(1,j));
                v(i+1,step+2) = v(i,step+2) + (dt / 6) * (kv1(2,j) + 2*kv2(2,j) + 2*kv3(2,j) + kv4(2,j));
                
                step += 3;
            }
        }


        
    }

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt);
};