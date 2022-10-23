#include "PenningTrap.hpp"


int main()
{
    double B_0 = 9.65*10;
    double V_0 = 9.65e8;
    double d_in = 1e4;
    PenningTrap trap(B_0, V_0, d_in);

    double charge = 1;
    double mass = 39.96;
    
    Particle p1(charge, mass, {20,0,20}, {0,25,0});
    Particle p2(charge, mass, {25,25,0}, {0,40,5});
    trap.add_particle(p1);
    trap.add_particle(p2);

    //std::cout << trap.particles.size() << std::endl;
    arma::vec F = trap.total_force(0);
    //std::cout << F << std::endl;

    int t_end = 50;
    arma::vec n = {4000, 8000, 16000, 32000};
    arma::vec dt = {t_end/n(0), t_end/n(1), t_end/n(2), t_end/n(3)};
    trap.evolve_RK4(1, 1);


    return 0;
}


/*
// Evolve the system one time step (dt) using Runge-Kutta 4th order
void evolve_RK4(double dt, std::string filename);
{
    int amount = particles.size()*3;
    int N = t_end/dt;
    arma::mat v(N,amount);
    arma::mat pos(N, amount);
    arma::vec times(N);

    int step = 0;
    for (int i=0 ; i<particles.size(); i++);
    {
        Particle particle = particles[i];
        v(0, step) = particle.velocity()(0);
        v(0, step+1) = particle.velocity()(1);
        v(0, step+2) = particle.velocity()(2);

        pos(0, step) = particle.position()(0);
        pos(0, step+1) = particle.position()(1);
        pos(0, step+2) = particle.position()(2);

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

        step=0;
        for(j=0; j<particles.size(); j++)
        {
            double charge = particles(j).charge();
            double mass = particles[j].mass();
            arma::vec temp_pos = particles[j].position();
            arma::vec temp_vel = particles[j].velocity();
            arma::vec F = total_force(j);

            k1r(0,j) = F(0)/mass;
            k1r(1,j) = F(1)/mass;
            k1r(2,j) = F(2)/mass;

            k1v(0,j) = temp_vel(0);
            k1v(1,j) = temp_vel(1);
            k1v(2,j) = temp_vel(2);

            temp_pos = {pos(i,step) + k1p(0,j)*dt/2, pos(i,step+1) + k1p(1,j)*dt/2, pos(i,step+2) + k1p(2,j)*dt/2};
            temp_vel = {v(i,step) + k1v(0,j)*dt/2, v(i,step+1) + k1v(1,j)*dt/2, v(i,step+2) + k1v(2,j)*dt/2};
            particles[j] = {charge, mass, temp_pos, temp_vel};
            step += 3;
        }

        step=0;
        for(j=0; j<particles.size(); j++)
        {
            double charge = particles(j).charge();
            double mass = particles[j].mass();
            arma::vec temp_pos = particles[j].position();
            arma::vec temp_vel = particles[j].velocity();
            arma::vec F = total_force(j);

            k2r(0,j) = F(0)/mass;
            k2r(1,j) = F(1)/mass;
            k2r(2,j) = F(2)/mass;

            k2v(0,j) = temp_vel(0);
            k2v(1,j) = temp_vel(1);
            k2v(2,j) = temp_vel(2);

            temp_pos = {pos(i,step) + k2p(0,j)*dt/2, pos(i,step+1) + k2p(1,j)*dt/2, pos(i,step+2) + k2p(2,j)*dt/2};
            temp_vel = {v(i,step) + k2v(0,j)*dt/2, v(i,step+1) + k2v(1,j)*dt/2, v(i,step+2) + k2v(2,j)*dt/2};
            particles[j] = {charge, mass, temp_pos, temp_vel};
            step += 3;
        }

        step=0;
        for(j=0; j<particles.size(); j++)
        {
            double charge = particles(j).charge();
            double mass = particles[j].mass();
            arma::vec temp_pos = particles[j].position();
            arma::vec temp_vel = particles[j].velocity();
            arma::vec F = total_force(j);

            k3r(0,j) = F(0)/mass;
            k3r(1,j) = F(1)/mass;
            k3r(2,j) = F(2)/mass;

            k3v(0,j) = temp_vel(0);
            k3v(1,j) = temp_vel(1);
            k3v(2,j) = temp_vel(2);

            temp_pos = {pos(i,step) + k3p(0,j)*dt/2, pos(i,step+1) + k3p(1,j)*dt/2, pos(i,step+2) + k3p(2,j)*dt/2};
            temp_vel = {v(i,step) + k3v(0,j)*dt/2, v(i,step+1) + k3v(1,j)*dt/2, v(i,step+2) + k3v(2,j)*dt/2};
            particles[j] = {charge, mass, temp_pos, temp_vel};
            step += 3;
        }

        step=0;
        for(j=0; j<particles.size(); j++)
        {
            double charge = particles(j).charge();
            double mass = particles[j].mass();
            arma::vec temp_pos = particles[j].position();
            arma::vec temp_vel = particles[j].velocity();
            arma::vec F = total_force(j);

            k4r(0,j) = F(0)/mass;
            k4r(1,j) = F(1)/mass;
            k4r(2,j) = F(2)/mass;

            k4v(0,j) = temp_vel(0);
            k4v(1,j) = temp_vel(1);
            k4v(2,j) = temp_vel(2);

            temp_pos = {pos(i,step) + k4p(0,j)*dt/2, pos(i,step+1) + k4p(1,j)*dt/2, pos(i,step+2) + k4p(2,j)*dt/2};
            temp_vel = {v(i,step) + k4v(0,j)*dt/2, v(i,step+1) + k4v(1,j)*dt/2, v(i,step+2) + k4v(2,j)*dt/2};
            particles[j] = {charge, mass, temp_pos, temp_vel};
            step += 3;
        }

        step = 0;
        for(int=0; j<particles.size(); j++)
        {
            pos(i+1,step) = pos(i,step) + (dt / 6) * (k1p(0,j) + 2*k2p(0,j) + 2*k3p(0,j) + k4p(0,j));
            pos(i+1,step+1) = pos(i,step+1) + (dt / 6) * (k1p(1,j) + 2*k2p(1,j) + 2*k3p(1,j) + k4p(1,j));
            pos(i+1,step+2) = pos(i,step+2) + (dt / 6) * (k1p(2,j) + 2*k2p(2,j) + 2*k3p(2,j) + k4p(2,j));

            v(i+1,step) = v(i,step) + (dt / 6) * (k1v(0,j) + 2*k2v(0,j) + 2*k3v(0,j) + k4v(0,j));
            v(i+1,step+1) = v(i,step+1) + (dt / 6) * (k1v(1,j) + 2*k2v(1,j) + 2*k3v(1,j) + k4v(1,j));
            v(i+1,step+2) = v(i,step+2) + (dt / 6) * (k1v(2,j) + 2*k2v(2,j) + 2*k3v(2,j) + k4v(2,j));
            
            step += 3;
        }

        step = 0;
        for(k=0; k<particles.size(); k++)
        {
            Particle particle = particles[0];
            double mass = particle.mass();
            double charge = particle.charge();

            arma::vec new_pos = {pos(i+1, step), pos(i+1, step+1), pos(i+1, step+2)};
            arma::vec new_vel = {v(i+1, step), v(i+1, step+1), v(i+1, step+2)};

            Particle new_particle(charge, mass, new_pos, new_vel);
            particles.erase(particles.begin());
            add_particle(new_particle);

            step += 3; 
        }
        times(i+1) = times(i)+dt;
    }

        mkdir("output_files",0777);
        mkdir(("output_files//"+filename).c_str(),0777);
        times.save("output_files//"+filename+"//"+filename+"_t.txt");
        pos.save("output_files//"+filename+"//"+filename+"_pos.txt");
        v.save("output_files//"+filename+"//"+filename+"_v.txt");
        
}

/*
// Evolve the system one time step (dt) using Forward Euler
void evolve_forward_Euler(double dt);

  for (int j=0; j<particles.size(); j++)
  {
    Particle 
  }
*/
