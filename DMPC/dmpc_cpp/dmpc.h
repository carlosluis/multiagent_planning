//
// Created by carlos on 16/03/18.
//
#ifndef DMPC_CPP_DMPC_H
#define DMPC_CPP_DMPC_H

#include <Eigen/Dense>
#include <eigen-quadprog/QuadProg.h>

using namespace Eigen;
using namespace std;

// Structure definitions
struct Constraint {
    MatrixXd A;
    VectorXd b;
};

struct Trajectory {
    MatrixXd pos;
    MatrixXd vel;
    MatrixXd acc;
};

struct Params {
    float h; // time step, in seconds
    int T; // Max time to complete trajectory
    int k_hor; // length of the prediction horizon
    int order; // order of the ellipsoid for collision constraint
    float c; // multiplier for constraint in the Z direction
    float rmin;
    float alim;
};

static const Params default_params = {0.2,20,15,2,1.5,0.5,0.5};

// Class definition
class SoftDMPC {
public:
    // Constructor
    SoftDMPC(Params params = default_params);
    // Destructor
    ~SoftDMPC(){};

    // Public variables

    // Public methods
    MatrixXd gen_rand_pts(int N, Vector3d pmin, Vector3d pmax, float rmin);
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
    MatrixXd _A0;

    // Private Methods
    bool check_collisions(MatrixXd prev_p, std::vector<MatrixXd> obs, int n);
    MatrixXd get_lambda_mat(int h, int K);
    MatrixXd get_delta_mat (int K);
    MatrixXd get_A0_mat (int K);
    Trajectory init_dmpc (MatrixXd po, MatrixXd pf);
    Constraint build_constraint (MatrixXd prev_p,
                                 MatrixXd po,
                                 MatrixXd vo,
                                 std::vector<MatrixXd> obs);
};

#endif //DMPC_CPP_DMPC_H
