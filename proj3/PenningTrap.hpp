#include "Particles.hpp"

class PenningTrap
{
public:
  double B_0;
  double V_0;
  double d;
  std::vector <Particle> particles;

  PenningTrap(double B0, double V0, double de){
    B_0 = B0;
    V_0 = V0;
    d = de;
  }

  // Add a particle to the trap
  void add_particle(Particle p_in){
    particles.push_back(p_in);
  }

  // External electric field at point r=(x,y,z)
  arma::vec external_E_field(arma::vec r){
    arma::vec E(3);
    if (arma::norm(r) < d){
      E(0) = -2*r(0);
      E(1) = -2*r(1);
      E(2) = 4*r(2);
      return -V_0/(2*d*d) * E;
  }

  //std::cout << arma::norm(r) << d;
  E = {0,0,0};
  return E;
  }

  // External magnetic field at point r=(x,y,z)
  arma::vec external_B_field(arma::vec r){
    if (arma::norm(r) < d){
      arma::vec B ={0,0,B_0};
    return B;}

    arma::vec BB = {0,0,0};
    return BB;
  }

  // Force on particle_i from particle_j
  arma::vec force_particle(int i, int j)
  {
      double ke = 1.38935333e5;
      Particle p_i = particles[i];
      Particle p_j = particles[j];

      double konst = ke * (p_i.charge * p_j.charge)/p_i.mass;
      double norm = pow(arma::norm(p_i.position - p_j.position),3);

      arma::vec F = konst*(p_i.position-p_j.position)/norm;
      return F;
  }

  // The total force on particle_i from the external fields
  arma::vec total_force_external(int i)
  {

      Particle p_i = particles[i];
      arma::vec B = external_B_field(p_i.position);
      arma::vec E = external_E_field(p_i.position);

      arma::vec F = p_i.charge * E + p_i.charge*arma::cross(p_i.velocity, B);
      return F;
    }

  // The total force on particle_i from the other particles
  arma::vec total_force_particles(int i){

    Particle p_i = particles[i];
    arma::vec force = {0,0,0};

    for (int j=0; j<particles.size();j++){

      if (i!=j){
        Particle p_j = particles[j];
        force += force_particle(i,j);
      }
    }
    return force;
  }

  // The total force on particle_i from both external fields and other particles
  arma::vec total_force(int i){
    arma::vec totforce = {0,0,0};
    //arma::vec g = {0,0,-9.81};
    arma::vec ext_fields = total_force_external(i);
    arma::vec otherparticles = total_force_particles(i);

    totforce += ext_fields + otherparticles;
    return totforce;
  }

  arma::vec f(double wz,double w0, arma::vec r, arma::vec v){
    arma::vec rr = {w0*v[1] + wz*r[0], -w0*v[0] + wz*r[1], -wz*r[2]};
    return rr;
  }

  void particle_show(){
    std::cout << "Position" << std::endl << particles[0].position << particles[1].position << std::endl <<"Velocity" << particles[0].velocity << particles[1].velocity <<std::endl;
  }

  // Evolve the system one time step (dt) using Runge-Kutta 4th order
  void evolve_RK4(double dt, int N)
  {
      for (int i = 0; i<particles.size();i++)
      {
        arma::vec r = particles[i].position;
        arma::vec v = particles[i].velocity;
        arma::vec F = total_force_particles(i);
        double m = particles[i].mass;

        double wz = 2*particles[i].charge*V_0/(particles[i].mass*d*d);
        double w0 = particles[i].charge*B_0/particles[i].mass;

        arma::vec k1v = dt*(f(wz, w0, r,v)/m + F);
        arma::vec k1r = dt*v;

        arma::vec k2v = dt*(f(wz, w0, r + k1r/2, v + k1v/2)/m + F);
        arma::vec k2r = dt*(v + k1v/2);

        arma::vec k3v = dt*(f(wz, w0, r + k2r/2, v + k2v/2)/m + F);
        arma::vec k3r = dt*(v + k2v/2);

        arma::vec k4v = dt*(f(wz, w0, r + k3r, v + k3v)/m + F);
        arma::vec k4r = dt*(v + k3v);

        double ty = 1.0/6.0;
        particles[i].position = r + ty*(k1r + 2*k2r + 2*k3r + k4r);
        particles[i].velocity = v + ty*(k1v + 2*k2v + 2*k3v + k4v);
        
        /*
        if (j % 10 == 0 & i == 0){
        std::cout << "Velocity " << i << std::endl << particles[i].velocity << "position " << i << std::endl << particles[i].position << std::endl;
      }*/
      }   
  }

  // Evolve the system one time step (dt) using Forward Euler
  void evolve_forward_Euler(double dt, int i);

void write_file_position_p1(std::string filename)
{
  std::ofstream ofile;
    ofile.open(filename);
    int width = 12;
    int prec = 4;
    
    
    
      ofile << std::setw(width) << std::setprecision(prec) << std::scientific << particles[0].position << std::endl;
    
    
    ofile.close();
}
/*
void write_file_velocity(std::string filename)
{
  
  std::ofstream ofile;
    ofile.open(filename);
    int width = 12;
    int prec = 4;

  
    
    ofile.close();
}*/
};
