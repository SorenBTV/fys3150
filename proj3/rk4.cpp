
arma::mat krv(int j, arma::mat kr, arma::mat kv, int step)
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