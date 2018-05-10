#include <iostream>
#include <Eigen/Dense>
#include <eigen-quadprog/src/QuadProg.h>
#include "dmpc.h"
#include <iomanip>
#include <chrono>

using namespace Eigen;
using namespace std;
using namespace std::chrono;

int main()

{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    Vector3d pmin;
    pmin << -2.5, -2.5, 0.2;
    Vector3d pmax;
    pmax << 2.5, 2.5, 2.2;
//    Params p = {0.4,20,15,2,1.5,0.5,0.5};
    DMPC test;
    int N = 70;
    float rmin = 0.91;
    MatrixXd po = test.gen_rand_pts(N,pmin,pmax,rmin);
    MatrixXd pf = test.gen_rand_pts(N,pmin,pmax,rmin);
//    cout << "po = " << endl << po << endl;
//    cout << "pf = " << endl << pf << endl;

//    Vector3d po1(0,0,1.5);
//    Vector3d po2(0,2,1.5);
//    Vector3d pf1(0,2,1.5);
//    Vector3d pf2(0,0,1.5);
//    MatrixXd po(3,2);
//    po << po1,po2;
//    MatrixXd pf(3,2);
//    pf << pf1,pf2;

    std::vector<Trajectory> solution = test.solveDMPC(po,pf);
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();

    cout << "Execution time = " << duration/1000000.0 << "s" << endl;

}