#ifndef __Particle_hpp__  
#define __Particle_hpp__

#include <string>
#include <armadillo>

class Particle
{
    public:

    double q;
    double m;
    arma::vec pos;
    arma::vec vel;

    Particle(double charge, double mass, arma::vec position, arma::vec velocity);

    double charge();
    double mass();
    arma::vec position();
    arma::vec velocity();
    void new_position(arma::vec new_r);
    void new_velocity(arma::vec new_v);



};

#endif