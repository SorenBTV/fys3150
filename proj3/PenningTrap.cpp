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

    double t_end = 50;
    double n = 4000;
    double dt = t_end/n;
    trap.evolve_RK4(dt, n);


    return 0;
}

/*
// Evolve the system one time step (dt) using Forward Euler
void evolve_forward_Euler(double dt);

  for (int j=0; j<particles.size(); j++)
  {
    Particle 
  }
*/
