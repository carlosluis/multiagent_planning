#include <iostream>
#include "dmpc.h"


using namespace Eigen;
using namespace std;
using namespace std::chrono;

int main()

{
    Vector3d pmin;
    pmin << -2.5, -2.5, 0.2;
    Vector3d pmax;
    pmax << 2.5, 2.5, 5.2;
    Params p = {0.2,50,15,2,2.0,0.35,1.0,2.0,100,0.01,0.05,1};
    DMPC test("ooqp",p);
    DMPC test2("ooqp",p);
//    DMPC test3("cplex",p);

    test.set_boundaries(pmin,pmax);
    test2.set_boundaries(pmin,pmax);
//    test3.set_boundaries(pmin,pmax);

    int N = 20;
    float rmin_init = 0.75;
    MatrixXd po = test.gen_rand_pts(N,pmin,pmax,rmin_init);
    MatrixXd pf = test.gen_rand_pts(N,pmin,pmax,rmin_init);

    test.set_final_pts(pf);
    test.set_initial_pts(po);
    test2.set_final_pts(pf);
    test2.set_initial_pts(po);
//    test3.set_final_pts(pf);
//    test3.set_initial_pts(po);

    test.set_k_factor(0);
    test2.set_k_factor(-1);

    cout << "FACTOR K_CTR = K" << endl;
    cout << "----------------------" << endl;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    std::vector<Trajectory> sol_parallel = test.solveParallelDMPCv2();
    std::vector<Trajectory> sol_para_short = test.solution_short;
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Total Parallel Execution Computation time = "
         << duration/1000000.0 << "s" << endl << endl;

    cout << "FACTOR K_CTR = K-1" << endl;
    cout << "----------------------" << endl;
    t1 = high_resolution_clock::now();
    std::vector<Trajectory> sol_parallel2 = test2.solveParallelDMPCv2();
    std::vector<Trajectory> sol_para2_short = test2.solution_short;
    t2 = high_resolution_clock::now();
    duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Total Parallel Execution Computation time = "
         << duration/1000000.0 << "s" << endl << endl;
////
//    cout << "CPLEX" << endl;
//    cout << "----------------------" << endl;
//    t1 = high_resolution_clock::now();
//    std::vector<Trajectory> sol_parallel3 = test3.solveParallelDMPCv2();
//    std::vector<Trajectory> sol_para3_short = test3.solution_short;
//    t2 = high_resolution_clock::now();
//    duration = duration_cast<microseconds>( t2 - t1 ).count();
//    cout << "Total Parallel Execution Computation time = "
//         << duration/1000000.0 << "s" << endl;

    // Write result to txt file (to be read by MATLAB)

    char const *file = "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/dmpc/cpp_results/trajectories.txt";
    test.trajectories2file(sol_para_short,file);
}