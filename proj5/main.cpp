#include "schrodinger.hpp"



void write_prob(schrodinger test, double sigmay, double v0,int slits, string filename)
{
    test.initialize(0.005, 2.5e-5, 0.008, 0.25, 0.5, 0.05, sigmay, 200, 0,
                    v0, slits, true);
    test.A_B();
    test.simulation();

    test.res.save(filename);
}



void write_simulation(schrodinger test, int slits)
{
    test.initialize(0.005, 2.5e-5, 0.002, 0.25, 0.5, 0.05, 0.2, 200, 0,
                    1e10, slits, true);

    test.A_B();
    test.simulation();

    if(slits==1)
    {
        test.res.save(string("data/single_slit.dat"));
    }
    else if(slits==2)
    {
        test.res.save(string("data/double_slit.dat"));
    }
    else if(slits==3)
    {
        test.res.save(string("data/triple_slit.dat"));
    } 
}


int main()
{

    // Problem 7

    schrodinger test1;
    write_prob(test1, 0.05, 0, 0, "data/test7.1.dat");

    schrodinger test2;
    write_prob(test2, 0.1, 1.0e10, 2, "data/test7.2.dat");


    //Problem 8/9
    schrodinger single_slit;
    write_simulation(single_slit, 1);

    schrodinger double_slit;
    write_simulation(double_slit, 2);

    schrodinger triple_slit;
    write_simulation(triple_slit, 3);



    return 0;
}
