#include <iostream>
#include "dmpc.h"


using namespace Eigen;
using namespace std;
using namespace std::chrono;

void test2file(const std::vector<MatrixXd> &time_vec,
               const VectorXd &cluster_size,
               const VectorXd &num_vehicles, char const* pathAndName){
    ofstream file(pathAndName, ios::out | ios::trunc);
    if(file)  // succeeded at opening the file
    {
        // instructions
        cout << "Writing solution to text file..." << endl;

        // Write the average time
        file << cluster_size.size() << " " <<  num_vehicles.size()
             << " " << time_vec.at(0).cols() <<  endl;
        file << cluster_size.transpose() << " "
             << num_vehicles.transpose() << endl;
        for (int i = 0; i < cluster_size.size(); ++i){
            file << time_vec.at(i) << endl;
        }
        file.close();  // close the file after finished
    }

    else
    {
        cerr << "Error while trying to open file" << endl;
    }
}

int main()

{
    // Test definitions
    VectorXd cluster_size(6);
    cluster_size << 1, 2, 4, 6, 8 , 10;
    VectorXd num_vehicles(2);
    num_vehicles << 10,20;
    int num_trials = 5;
    float rmin_init = 0.75;
    std::string solver = "ooqp";

    // Trial set up
    Vector3d pmin;
    pmin << -3.0, -3.0, 0.2;
    Vector3d pmax;
    pmax << 3.0, 3.0, 3.2;
    Params p = {0.2,30,15,2,2.0,0.35,1.0,2.0,100,0.01,0.05,1};

    int size_clust_arr = cluster_size.size();
    int size_vehic_arr = num_vehicles.size();

    // Create the instances of DMPC with different cluster size
    std::vector<DMPC> dmpc_vec;

    // Store all the times, to be processed by MATLAB
    std::vector<MatrixXd> time_vec;
    for(int i=0; i<size_clust_arr; ++i){
        DMPC temp(solver,p);
        temp.set_boundaries(pmin,pmax);
        temp.set_cluster_num(cluster_size[i]);
        dmpc_vec.push_back(temp);
        time_vec.push_back(MatrixXd::Zero(size_vehic_arr,num_trials));
    }

    for(int i=0; i<size_vehic_arr; ++i){
        int success = 0;
        for(int j=0; success<num_trials; ++j){
            cout << "Doing trial #" << success
                 << " with " << num_vehicles[i] << " agents" << endl;
            // Generate the test
            MatrixXd po = dmpc_vec.at(0).gen_rand_pts(num_vehicles[i],
                                                      pmin,pmax,rmin_init);
            MatrixXd pf = dmpc_vec.at(0).gen_rand_pts(num_vehicles[i],
                                                      pmin,pmax,rmin_init);

            for (int k=0; k<size_clust_arr; ++k){
                // Set initial and final pts
                dmpc_vec.at(k).set_initial_pts(po);
                dmpc_vec.at(k).set_final_pts(pf);
                int  l = 0;

                // Solve the problem and get the computation time
                high_resolution_clock::time_point t1 = high_resolution_clock::now();
                dmpc_vec.at(k).solveParallelDMPCv2();
                high_resolution_clock::time_point t2 = high_resolution_clock::now();
                auto duration = duration_cast<microseconds>( t2 - t1 ).count();

                // Check if it was successful before recording it
                if (dmpc_vec.at(k).successful){
                    time_vec.at(k)(i,success) = duration/1000000.0;
                }
                else continue;
            }
            if (dmpc_vec.at(0).successful)
                success++;
        }
    }

    // Write result to txt file (to be read by MATLAB)

    char const *file = "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/CPP_results/cluster_test.txt";
    test2file(time_vec,cluster_size,num_vehicles,file);
}