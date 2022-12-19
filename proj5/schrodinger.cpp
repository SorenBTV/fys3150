#include "schrodinger.hpp"



void schrodinger::initialize(double h, double dt, double T, double xc, double yc, double sigmax, double sigmay,
                             double px, double py, double v0, int nslits, bool save_prob)
{
    
    save_prob_ = save_prob;
    h_ = h;
    M_ = 1./h_ + 1;
    dt_ = dt;
    T_ = T;
    x_c = xc;
    y_c = yc;
    sigma_x = sigmax;
    sigma_y = sigmay;
    p_x = px;
    p_y = py;
    v0_ = v0;
    nslits_ = nslits;


    n = T_/dt_;


    a_ = arma::cx_vec((M_-2)*(M_-2));
    b_ = arma::cx_vec((M_-2)*(M_-2));

    u_ = arma::cx_vec((M_-2)*(M_-2), arma::fill::zeros);
    initilalize_u();


    V_ = arma::zeros(M_, M_);
    if(nslits_ > 0)
    {
        initialize_V();
    }

    A_ = arma::sp_cx_mat((M_-2)*(M_-2), (M_-2)*(M_-2));
    B_ = arma::sp_cx_mat((M_-2)*(M_-2), (M_-2)*(M_-2));


}

void schrodinger::initilalize_u()
{
    arma::cx_double val;
    arma::cx_double sum = 0;
    double x;
    double y;

    for(int i=1; i<=M_-2; i++)
    {
        x = i*h_;

        for(int j=1; j<=M_-2; j++)
        {   
            y = j*h_;

            val = exp(-(pow(x - x_c, 2)/(2 * pow(sigma_x, 2))) - (pow(y - y_c, 2) / (2 * pow(sigma_y, 2)))
                + 1.0i * p_x * (x - x_c) + 1.0i * p_y * (y - y_c));

            sum += val* conj(val);
            u_(ij_k(i, j)) = val;
        }
    }
    
    for (int k=0; k<(M_-2)*(M_-2); k++)
    {
        u_(k)/=sqrt(sum);
    }
    
}


int schrodinger::ij_k(int i, int j)
{
    return (i - 1) + (M_ - 2) * (j - 1);  
    
}




void schrodinger::A_B()
{
    r_ = arma::cx_double(0,dt_/2./h_/h_);
    arma::cx_double constant;

    for (int i=1; i<=M_-2; i++)
    {
        for (int j=1; j<=M_-2; j++)
        {
            constant = 4.0*r_ + arma::cx_double(0,dt_/2.*V_(i,j));
            int k = ij_k(i, j);
            a_(k) = 1.0 + constant;
            b_(k) = 1.0 - constant;
        }
    }
    A_.diag() = a_;
    A_.diag(1) += -r_;
    A_.diag(-1) += -r_;
    A_.diag(M_-2) += -r_;
    A_.diag(-M_+2) += -r_;

    for (int i=M_-2; i<(M_-2)*(M_-2); i+=M_-2)
    {
        A_(i, i-1) = 0;
        A_(i-1, i) = 0;
    }

    B_ = A_ * -1;
    B_.diag() = b_;

}

void schrodinger::update_u()
{
    b_ = B_ * u_;
    u_ = arma::spsolve(A_, b_);
}


void schrodinger::initialize_V()
{

    double th = 0.02;
    double sep = 0.05;
    double ap = 0.05;
    double mid = 0.5;
    int odd = nslits_ % 2;

    // runs over all indices
    for (int j = 0; j < M_; j++)
    {
        for (int i = 0; i < M_; i++)
        {
            // Fill for middle x-values
            if (i * h_ >= mid - th / 2 && i * h_ <= mid + th / 2)
            {
                // opening for odd numbered amount of slits
                if (odd == 1 && j * h_ < mid + ap / 2 && j * h_ > mid - ap / 2)
                {
                    V_(i, j) = 0;
                }
                else
                {
                    V_(i, j) = v0_;
                }
                // Fill out rest of slits
                for (int k = 0; k < (nslits_ - odd) / 2; k++)
                {
                    if (
                        // Setting up boundaries and slits
                        j * h_ >= (k + 0.5 * odd) * ap + (k + 1 + 0.5 * (odd - 1)) * sep + mid           
                            && j * h_ < (k + 1 + 0.5 * odd) * ap + (k + 1 + 0.5 * (odd - 1)) * sep + mid 
                        
                        || (j * h_ < -(k + 0.5 * odd) * ap - (k + 1 + 0.5 * (odd - 1)) * sep + mid           
                            && j * h_ >= -(k + 1 + 0.5 * odd) * ap - (k + 1 + 0.5 * (odd - 1)) * sep + mid)) 
                    {
                        V_(i, j) = 0; 
                    }
                }
            }
        }
    }
   if(nslits_ == 1)
   {
        V_.save(string("data/V_1.dat"));
   }
   else if(nslits_ == 2)
   {
        V_.save(string("data/V_2.dat"));
   }
   else if(nslits_ == 3)
   {
        V_.save(string("data/V_3.dat"));
   }
}


void schrodinger::simulation()
{
    res = arma::cx_cube(n+1, M_, M_);
    res.zeros();

    for (int i=1; i<M_-1; i++)
    {
        for (int j=1; j<M_-1; j++)
            {
                res(0,i,j) = u_(ij_k(i,j));
            }
    }


    for (int k=1; k<=n; k++)
    {
        cout << k << endl;
        update_u();
        
        for (int i=1; i<M_-1; i++)
        {
            for (int j=1; j<M_-1; j++)
            {
                res(k,i,j) = u_(ij_k(i,j));
            }
        }
    }


}   

