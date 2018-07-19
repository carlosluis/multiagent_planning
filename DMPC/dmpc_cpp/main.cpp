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
    DMPC test("cplex");
    DMPC test2("ooqp");
    int N = 1;
    float rmin_init = 0.75;
//    MatrixXd po = test.gen_rand_pts(N,pmin,pmax,rmin_init);
//    MatrixXd pf = test.gen_rand_perm(po);

    Vector3d po1(-1.0, 1.0, 1.0);
    Vector3d po2(-0.5, 1.0, 1.0);
    Vector3d po3(0.0, 1.0, 1.0);
    Vector3d po4(0.5, 1.0, 1.0);
    Vector3d po5( 1.0, 1.0, 1.0);
    Vector3d po6(-1.0, 0.5, 1.0);
    Vector3d po7( -0.5, 0.5, 1.0);
    Vector3d po8( 0.0, 0.5, 1.0);
    Vector3d po9( 0.5, 0.5, 1.0);
    Vector3d po10( 1.0, 0.5, 1.0);
    Vector3d po11( -1.0, 0.0, 1.0);
    Vector3d po12( -0.5, 0.0, 1.0);
    Vector3d po13( 0.0, 0.0, 1.0);
    Vector3d po14( 0.5, 0.0, 1.0);
    Vector3d po15( 1.0, 0.0, 1.0);
    Vector3d po16(-1.0, -0.5, 1.0);
    Vector3d po17( -0.5, -0.5, 1.0);
    Vector3d po18( 0.0, -0.5, 1.0);
    Vector3d po19( 0.5, -0.5, 1.0);
    Vector3d po20( 1.0, -0.5, 1.0);
    Vector3d po21(-1.0, -1.0, 1.0);
    Vector3d po22(-0.5, -1.0, 1.0);
    Vector3d po23(0.0, -1.0, 1.0);
    Vector3d po24(0.5, -1.0, 1.0);
    Vector3d po25(1.0, -1.0, 1.0);

    MatrixXd po(3,25);
//    po << po1,po2;
    po << po1,po2,po3,po4,po5,po6,po7,po8,po9,po10,
            po11,po12,po13,po14,po15,po16,po17,po18,po19,po20,
            po21,po22,po23,po24,po25;

    MatrixXd pf(3,25);
//    pf << po2,po1;
    pf << po25,po24,po23,po22,po21,po20,po19,po18,po17,po16,po15,
            po14,po13,po12,po11,po10,po9,po8,po7,po6,po5,
            po4,po3,po2,po1;
//    pf = test.gen_rand_perm(po);

    test.set_final_pts(pf);
    test.set_initial_pts(po);
    test2.set_final_pts(pf);
    test2.set_initial_pts(po);

    cout << "CPLEX" << endl;
    cout << "----------------------" << endl;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    std::vector<Trajectory> sol_parallel = test.solveParallelDMPCv2(po,pf);
    std::vector<Trajectory> sol_para_short = test.solution_short;
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Total Parallel Execution Computation time = "
         << duration/1000000.0 << "s" << endl << endl;

    cout << "EIGEN-QUADPROG" << endl;
    cout << "----------------------" << endl;
    t1 = high_resolution_clock::now();
    std::vector<Trajectory> sol_parallel2 = test2.solveParallelDMPCv2(po,pf);
    std::vector<Trajectory> sol_para2_short = test2.solution_short;
    t2 = high_resolution_clock::now();
    duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Total Parallel Execution Computation time = "
         << duration/1000000.0 << "s" << endl;

    // Write result to txt file (to be read by MATLAB)

    char const *file = "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/CPP_results/trajectories.txt";
    test.trajectories2file(sol_para_short,file);
}