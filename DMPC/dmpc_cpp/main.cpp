#include <iostream>
#include "dmpc.h"

using namespace Eigen;
using namespace std;
using namespace std::chrono;

int main()

{
    Vector3d pmin;
    pmin << -4, -4, 0.2;
    Vector3d pmax;
    pmax << 4, 4, 3.2;
//    Params p = {0.4,20,15,2,1.5,0.5,0.5};
    DMPC test;
    int N = 100;
    float rmin_init = 0.91;
    MatrixXd po = test.gen_rand_pts(N,pmin,pmax,rmin_init);
    MatrixXd pf = test.gen_rand_perm(po);

//    cout << "po = " << endl << po << endl;
//    cout << "pf = " << endl << pf << endl;

//    Vector3d po1(0.01,0,1.5);
//    Vector3d po2(0,2,1.5);
//    Vector3d pf1(0.01,0,1.5);
//    Vector3d pf2(0,-0.5,1.5);
//    MatrixXd po(3,2);
//    po << po1,po2;
//    MatrixXd pf(3,2);
//    pf << pf1,pf2;

    test.set_final_pts(pf);
    test.set_initial_pts(po);

    cout << "PARALLEL EXECUTION" << endl;
    cout << "----------------------" << endl;

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    std::vector<Trajectory> sol_parallel = test.solveParallelDMPC(po,pf);
    std::vector<Trajectory> sol_para_short = test.solution_short;
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Total Parallel Execution Computation time = "
         << duration/1000000.0 << "s" << endl << endl;

    cout << "SEQUENTIAL EXECUTION" << endl;
    cout << "----------------------" << endl;

    t1 = high_resolution_clock::now();
    std::vector<Trajectory> sol_seq = test.solveDMPC(po,pf);
    std::vector<Trajectory> sol_seq_short = test.solution_short;
    t2 = high_resolution_clock::now();
    duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Total Sequential Execution Computation time = "
         << duration/1000000.0 << "s" << endl;

    // Write result to txt file (to be read by MATLAB)

    char const *file = "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/CPP_results/trajectories.txt";
    test.trajectories2file(sol_para_short,file);
}





