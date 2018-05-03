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

static const Params default_params = {0.2,10,15,2,1.5,0.5,0.5};

// Class definition
class DMPC {
public:
    // Constructor
    DMPC(Params params = default_params);
    // Destructor
    ~DMPC(){};

    // Public variables

    // Public methods
    MatrixXd gen_rand_pts(int N, Vector3d pmin, Vector3d pmax, float rmin);
    void set_initial_pts(MatrixXd po);
    void set_final_pts(MatrixXd pf);\

private:
    // Private Variables

    // Algorithm parameters
    float _h; // time step, in seconds
    int _T; // Max time to complete trajectory
    float _K; // number of time steps for the trajectory
    int _k_hor; // length of the prediction horizon
    int _order; // order of the ellipsoid for collision constraint
    float _c; // multiplier for constraint in the Z direction
    float _rmin;
    float _alim;

    // Workspace boundaries
    Vector3d _pmin;
    Vector3d _pmax;

    // Goals
    MatrixXd _po; // set of initial positions
    MatrixXd _pf; // set of final positions

    // Ellipsoid variables
    Matrix3d _E;
    Matrix3d _E1;
    Matrix3d _E2;

    // Model-related matrices
    Matrix<double, 6, 6> _A;
    Matrix<double, 6, 3> _b;
    MatrixXd _Lambda; // Matrix to recover position from acceleration
    MatrixXd _Delta; // Used for input variation computation
    MatrixXd _A0; // Propagation of initial states in position

    // Private Methods
    bool check_collisions(Vector3d prev_p, std::vector<MatrixXd> obs, int n, int k);
    MatrixXd get_lambda_mat(int h, int K);
    MatrixXd get_delta_mat (int K);
    MatrixXd get_A0_mat (int K);
    Trajectory init_dmpc (Vector3d po, Vector3d pf);
    Constraint build_collconstraint (Vector3d prev_p,
                                 Vector3d po,
                                 Vector3d vo,
                                 std::vector<MatrixXd> obs, int n, int k);
    Trajectory solveDMPC(Vector3d po,Vector3d pf,
                                      Vector3d vo,Vector3d ao,
                                      int n, std::vector<MatrixXd> obs);
};

#endif //DMPC_CPP_DMPC_H
