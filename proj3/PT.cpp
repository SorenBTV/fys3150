#include "PenningTrap.hpp"


int main(){
  double B_0 = 9.65 * 10;
  double V_0 = 2.41 * 1000000;
  double d = 500;

  //Particle

  double q = 1;
  double m = 39.96;
  arma::vec posvec1 = {20,0,20};
  arma::vec posvec2 = {25,25,0};
  arma::vec velvec1 = {0,25,0};
  arma::vec velvec2 = {0,40,5};



  PenningTrap trap(B_0, V_0, d);
  Particle p1(q,m,posvec1,velvec1);
  Particle p2(q,m,posvec2,velvec2);


  trap.add_particle(p1);
  trap.add_particle(p2);


  double n = 4500;
  double t_end = 50;
  double dt = t_end/n;

  //cout << trap.external_E_field(particle1.position);


  //trap.particle_show();

    arma::vec x_p1 = arma::vec(n, arma::fill::zeros);
    arma::vec y_p1 = arma::vec(n, arma::fill::zeros);
    arma::vec z_p1 = arma::vec(n, arma::fill::zeros);
    //velocity data (p1)
    arma::vec vx_p1 = arma::vec(n, arma::fill::zeros);
    arma::vec vy_p1 = arma::vec(n, arma::fill::zeros);
    arma::vec vz_p1 = arma::vec(n, arma::fill::zeros);
    
    //position data (p2)
    arma::vec x_p2 = arma::vec(n, arma::fill::zeros);
    arma::vec y_p2 = arma::vec(n, arma::fill::zeros);
    arma::vec z_p2 = arma::vec(n, arma::fill::zeros);
    //velocity data (p2)
    arma::vec vx_p2 = arma::vec(n, arma::fill::zeros);
    arma::vec vy_p2 = arma::vec(n, arma::fill::zeros);
    arma::vec vz_p2 = arma::vec(n, arma::fill::zeros);

    x_p1(0) = p1.position(0);
    y_p1(0) = p1.position(1);
    z_p1(0) = p1.position(2);

    x_p2(0) = p2.position(0);
    y_p2(0) = p2.position(1);
    z_p2(0) = p2.position(2);

    vx_p1(0) = p1.velocity(0);
    vy_p1(0) = p1.velocity(1);
    vz_p1(0) = p1.velocity(2);
    vx_p2(0) = p2.velocity(0);
    vy_p1(0) = p1.velocity(1);
    vz_p2(0) = p2.velocity(2);
    

for (int j=0; j<n; j++)
{
  trap.evolve_RK4(dt,n);

  arma::vec r = trap.particles[0].position;
  arma::vec v = trap.particles[0].velocity;

  x_p1(j) = r(0);
  y_p1(j) = r(1);
  z_p1(j) = r(2);

  vx_p1(j) = v(0);
  vy_p1(j) = v(1);
  vz_p1(j) = v(2);

  arma::vec r1 = trap.particles[1].position;
  arma::vec v1 = trap.particles[1].velocity;

  x_p2(j) = r1(0);
  y_p2(j) = r1(1);
  z_p2(j) = r1(2);

  vx_p2(j) = v1(0);
  vy_p2(j) = v1(1);
  vz_p2(j) = v1(2);
}
    x_p1.save("data/x_values_RK4_p1.txt", arma::raw_ascii);
    y_p1.save("data/y_values_RK4_p1.txt", arma::raw_ascii);
    z_p1.save("data/z_values_RK4_p1.txt", arma::raw_ascii);

    x_p2.save("data/x_values_RK4_p2.txt", arma::raw_ascii);
    y_p2.save("data/y_values_RK4_p2.txt", arma::raw_ascii);
    z_p2.save("data/z_values_RK4_p2.txt", arma::raw_ascii);

    vx_p1.save("data/v_x_values_RK4_p1.txt", arma::raw_ascii);
    vy_p1.save("data/v_x_values_RK4_p1.txt", arma::raw_ascii);
    vz_p1.save("data/v_z_values_RK4_p1.txt", arma::raw_ascii);

    vx_p2.save("data/v_x_values_RK4_p2.txt", arma::raw_ascii);
    vy_p2.save("data/v_x_values_RK4_p2.txt", arma::raw_ascii);
    vz_p2.save("data/v_z_values_RK4_p2.txt", arma::raw_ascii);

  //trap.particle_show();

  //trap.write_file_position_p1("Position_x4500.txt", n);


  return 0;
}

