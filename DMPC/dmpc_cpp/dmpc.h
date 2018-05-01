//
// Created by carlos on 16/03/18.
//
#ifndef DMPC_CPP_DMPC_H
#define DMPC_CPP_DMPC_H

#include <Eigen/Dense>
#include <eigen-quadprog/QuadProg.h>

using namespace Eigen;
using namespace std;

struct Constraint {
    MatrixXd A;
    VectorXd b;
};

struct Trajectory {
    MatrixXd pos;
    MatrixXd vel;
    MatrixXd acc;
};

class SoftDMPC {
public:
    // Constructor
    SoftDMPC();
    // Destructor
    ~SoftDMPC(){};

    // Public methods

private:
    // Private Variables

    float _h; // time step, in seconds
    int _T; // Max time to complete trajectory
    float _K; // number of time steps for the trajectory
    int _k_hor; // length of the prediction horizon
    int _order; // order of the ellipsoid for collision constraint
    float _c; // multiplier for constraint in the Z direction
    float _rmin;
    float _alim;
    Matrix3d _E;
    Matrix3d _E1;
    Matrix3d _E2;
    Matrix<double, 6, 6> _A;
    Matrix<double, 6, 3> _b;
    MatrixXd _Lambda;
    MatrixXd _Delta;


    // Private Methods
    bool check_collisions(MatrixXd prev_p, std::vector<MatrixXd> obs, int n);
    MatrixXd get_lambda_mat(int h, int K);
    MatrixXd get_delta_mat (int k_hor);
    Trajectory init_dmpc (MatrixXd po, MatrixXd pf);
    Constraint build_constraint (MatrixXd prev_p,
                                 MatrixXd po,
                                 MatrixXd vo,
                                 std::vector<MatrixXd> obs);



};

#endif //DMPC_CPP_DMPC_H
