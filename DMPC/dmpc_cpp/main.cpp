#include <iostream>
#include <Eigen/Dense>
#include <eigen-quadprog/QuadProg.h>
#include "dmpc.h"


using namespace Eigen;
using namespace std;
int main()

{
    std::vector<MatrixXd> hola;
    MatrixXd m = MatrixXd::Random(3,3);
    m = (m + MatrixXd::Constant(3,3,1.2)) * 100;
    cout << "m =" << endl << m << endl;
    hola.push_back(m);
    hola.push_back(MatrixXd::Constant(3,3,1.2));
    cout << "hola[1] = " << hola.at(0) << endl;
    cout << "hola[2] = " << hola.size() << endl;
    Vector3d pmin;
    pmin << -2.5, -2.5, 0.2;
    Vector3d pmax;
    pmax << 2.5, 2.5, 2.2;
    Params p = {0.4,20,15,2,1.5,0.5,0.5};
    SoftDMPC test(p);
    int N = 4;
    float rmin = 0.91;
    test.gen_rand_pts(N,pmin,pmax,rmin);

}