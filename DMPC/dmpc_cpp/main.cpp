#include <iostream>
#include "dmpc.h"

using namespace Eigen;
using namespace std;
using namespace std::chrono;

int main()

{
    Vector3d pmin;
    pmin << -1.0, -1.0, 0.2;
//    pmin << -4, -4, 0.2;
    Vector3d pmax;
    pmax << 1.0, 1.0, 2.2;
//    pmax << 4, 4, 3.2;
//    Params p = {0.4,20,15,2,1.5,0.5,0.5};
    DMPC test;
    DMPC test2;
    int N = 25;
    float rmin_init = 0.75;
    MatrixXd po = test.gen_rand_pts(N, pmin, pmax, rmin_init);
    MatrixXd pf = test.gen_rand_perm(po);

//    cout << "po = " << endl << po << endl;
//    cout << "pf = " << endl << pf << endl;

//    Vector3d po1(-1.0, 1.0, 1.0);
//    Vector3d po2( 0.0, 1.0, 0.8);
//    Vector3d po3(1.0, 1.0, 1.5);
//    Vector3d po4( -1.0, 0.0, 0.4);
//    Vector3d po5 (0.0, 0.0 , 1.3);
//    Vector3d po6 (1.0, 0.0 , 0.7);
//    Vector3d po7( -1.0, -1.0 , 0.9);
//    Vector3d po8( 0.0, -1.0 , 1.4);
//    Vector3d po9 (1.0, -1.0 , 0.6);

//    Vector3d po1(1.0,1.0,1.5);
//    Vector3d po2(-1,-1,1.5);
//    Vector3d po3(-1,1,1.5);
//    Vector3d po4(1,-1,1.5);
//
//    Vector3d pf1(-1.0,-1.0,1.5);
//    Vector3d pf2(1,1,1.5);
//    Vector3d pf3(1,-1,1.5);
//    Vector3d pf4(-1,1,1.5);
//
//    MatrixXd po(3,4);
//    po << po1,po2,po3,po4;
//    MatrixXd pf(3,4);
//    pf << pf1,pf2,pf3,pf4;
//    pf = test.gen_rand_perm(po);

    test.set_final_pts(pf);
    test.set_initial_pts(po);
    test2.set_final_pts(pf);
    test2.set_initial_pts(po);

    cout << "OOQP" << endl;
    cout << "----------------------" << endl;
    test.use_OOQP();
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    std::vector<Trajectory> sol_parallel = test.solveParallelDMPCv2(po, pf);
    std::vector<Trajectory> sol_para_short = test.solution_short;
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(t2 - t1).count();
    cout << "Total Parallel Execution Computation time = "
         << duration / 1000000.0 << "s" << endl << endl;

    cout << "GUROBI" << endl;
    cout << "----------------------" << endl;
//    test2.use_OOQP();
    t1 = high_resolution_clock::now();
    std::vector<Trajectory> sol_parallel2 = test2.solveParallelDMPCv2(po,pf);
    std::vector<Trajectory> sol_para2_short = test2.solution_short;
    t2 = high_resolution_clock::now();
    duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Total Parallel Execution Computation time = "
         << duration/1000000.0 << "s" << endl;

    // Write result to txt file (to be read by MATLAB)

    char const *file = "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/CPP_results/trajectories.txt";
    test.trajectories2file(sol_para_short, file);
}