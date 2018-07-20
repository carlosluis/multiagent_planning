//
// Created by carlos on 16/03/18.
//
#ifndef DMPC_CPP_DMPC_H
#define DMPC_CPP_DMPC_H

#include <Eigen/Dense>
#include <eigen-quadprog/src/QuadProg.h>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <chrono>
#include <thread>
#include <fstream>
#include <algorithm>
#include <ooqp_eigen_interface/OoqpEigenInterface.hpp>
#define IL_STD 1
#include <ilcplex/ilocplex.h>

using namespace Eigen;
using namespace std;

// Structure definitions
struct Constraint {
    MatrixXd A;
    VectorXd b;
    VectorXd prev_dist;
};

struct Trajectory {
    MatrixXd pos;
    MatrixXd vel;
    MatrixXd acc;
};

struct Params {
    float h; // time step, in seconds
    int maxT; // Max time to complete trajectory
    int k_hor; // length of the prediction horizon
    int order; // order of the ellipsoid for collision constraint
    float c; // multiplier for constraint in the Z direction
    float rmin;
    float alim;
};

static const Params default_params = {0.2,30,15,2,2.0,0.35,1.0};

// Class definition
class DMPC {
public:
    // Constructor
    DMPC(std::string solver_name = "quadprog", Params params = default_params);
    // Destructor
    ~DMPC(){};

    // Public variables

    // Public methods

    // Goal generation
    MatrixXd gen_rand_pts(const int &N,
                          const Vector3d &pmin,
                          const Vector3d &pmax,
                          const float &rmin);
    MatrixXd gen_rand_perm (const MatrixXd &po);

    // Setters
    void set_initial_pts(const MatrixXd &po);
    void set_final_pts(const MatrixXd &pf);
    void use_OOQP();

    // Top level method to be called by program
    std::vector<Trajectory> solveDMPC(const MatrixXd &po,
                                      const MatrixXd &pf);
    std::vector<Trajectory> solveParallelDMPC(const MatrixXd &po,
                                              const MatrixXd &pf);
    std::vector<Trajectory> solveParallelDMPCv2(const MatrixXd &po,
                                                const MatrixXd &pf);

    void trajectories2file(const std::vector<Trajectory> &src,
                           char const* pathAndName);

    // public variables

    std::vector<Trajectory> solution_short;


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
    float _h_scaled;

    std::string _solver_name;

    // Workspace boundaries
    Vector3d _pmin;
    Vector3d _pmax;

    // QP solver
    QuadProgDense _qp;

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
    MatrixXd _A_v; // Matrix to recover velocity from acceleration
    MatrixXd _Delta; // Used for input variation computation
    MatrixXd _A0; // Propagation of initial states in position

    int _fail; //keeps track if QP failed or not
    int _bla;
    bool _execution_ended;
    int _failed_i;

    //CPLEX variables
    std::vector<CPXENVptr> _env;
    std::vector<CPXLPptr> _lp;

    // Private Methods

    // Matrix building methods
    void get_lambda_A_v_mat(const int &K);
    void get_delta_mat (const int &K);
    void get_A0_mat (const int &K);

    // Initialization method of algorithm
    Trajectory init_dmpc (const Vector3d &po,
                          const Vector3d &pf);

    // Collision check and constraint construction
    bool check_collisions(const Vector3d &prev_p,
                          const std::vector<MatrixXd> &obs,
                          const int &n, const int &k);
    std::vector<bool> check_collisionsv2(const Vector3d &prev_p,
                                         const std::vector<MatrixXd> &obs,
                                         const int &n, const int &k);

    Constraint build_collconstraint (const Vector3d &prev_p,
                                     const Vector3d &po,
                                     const Vector3d &vo,
                                     const std::vector<MatrixXd> &obs,
                                     const int &n, const int &k);

    Constraint build_collconstraintv2 (const Vector3d &prev_p,
                                       const Vector3d &po,
                                       const Vector3d &vo,
                                       const std::vector<MatrixXd> &obs,
                                       const std::vector<bool> &violation_vec,
                                       const int &n, const int &k);

    // Optimization problem solving
    Trajectory solveQP(const Vector3d &po, const Vector3d &pf,
                       const Vector3d &vo, const Vector3d &ao,
                       const int &n, const std::vector<MatrixXd> &obs);

    Trajectory solveQPv2(const Vector3d &po, const Vector3d &pf,
                       const Vector3d &vo, const Vector3d &ao,
                       const int &n, const std::vector<MatrixXd> &obs,
                        const int &id_cluster);

    // Post checks
    bool reached_goal(const std::vector<Trajectory> &all_trajectories,
                      const MatrixXd &pf, const float &error_tol, const int &N);

    bool reached_goalv2(const std::vector<Trajectory> &all_trajectories,
                      const MatrixXd &pf, const float &error_tol, const int &N, const int &k);

    bool collision_violation(const std::vector<Trajectory> &solution);

    std::vector<Trajectory> interp_trajectory(const std::vector<Trajectory> &sol,
                                              const double &step_size);

    double get_trajectory_time(const std::vector<Trajectory> &solution);

    void scale_solution(std::vector<Trajectory> &sol,
                        const int &vmax, const int &amax);

    void cluster_solve(const int &k,
                       std::vector<Trajectory> &all_trajectories,
                       std::vector<MatrixXd> &obs,
                       const std::vector<int> &agents,
                       const std::vector<MatrixXd> &prev_obs);

    void cluster_solvev2(const int &k,
                         std::vector<Trajectory> &all_trajectories,
                         std::vector<MatrixXd> &obs,
                         const std::vector<int> &agents,
                         const std::vector<MatrixXd> &prev_obs,
                         const int id_cluster);
    void init_cplex(int id);
    void terminate_cplex(int id);
    void eigen_to_cplex(const Eigen::MatrixXd &H, int *&matbeg,
                        int *&matcnt, int *&matind, double *&matval,
                        int &numrows, int &numcols, int &numnz);
};

#endif //DMPC_CPP_DMPC_H
