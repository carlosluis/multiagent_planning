/*
 * dmpc.h
 *
 *  Created On : Mar 16, 2018
 *      Author : Carlos Luis
 *      Email  : carlos.luis@robotics.utias.utoronto.ca
 */

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

/**
 * Structure to specify a constraint of type Ax < b
 */

struct Constraint {
    MatrixXd A;           /**< Matrix A of dimensions (nxn) */
    VectorXd b;           /**< Vector b of dimensions (nx1) */
    VectorXd prev_dist;   /**< Distance to predicted colliding neighbours */
};

/**
 * Structure to specify a time-parameterized trajectory of position,
 * velocity and acceleration for a double-integrator agent in 3D
 */
struct Trajectory {
    MatrixXd pos;  /**< Position in XYZ  */
    MatrixXd vel;  /**< Velocity in XYZ */
    MatrixXd acc;  /**< Acceleration in XYZ */
};

/**
 * Structure to hold several parameters used in the DMPC algorithm
 * All of these are configured in the controller_param_example YAML file
 */
struct Params {
    float h;             /**< Time step for the computation */
    int T;               /**< Max time to complete trajectory */
    int k_hor;           /**< Length of the prediction horizon */
    int order;           /**< Order of the ellipsoid for collision constraint */
    float c;             /**< Multiplier for constraint in the Z direction */
    float rmin;          /**< Min inter-agent distance for collision in XY */
    float alim;          /**< Acceleration limit for every component XYZ */
    float vlim;          /**< Velocity limit to scale the final solution */
    int freq;            /**< Frequency for trajectory interpolation */
    float goal_tol;      /**< Allowed distance error to goal */
    float collision_tol; /**< Allowed violation of collision constraint */
    int speed;           /**< Number of time steps for error calculation */
};

// Default parameters
static const Params default_params = {0.2,10,12,2,1.5,0.5,2.0,2.0,
                                      100,0.05,0.05,1};

// Class definition
class DMPC {
public:
    /**
     * Load parameters into private variables of the class, initialize
     * variables, construct MPC matrices
     * @param params DMPC parameters, read in YAML file
     * @param solver_name choose which solver to use
     */
    DMPC(std::string solver_name = "quadprog", Params params = default_params);
    ~DMPC(){};

    // Public variables

    std::vector<Trajectory> solution_short; /**< Solution pre-interpolation */
    bool successful;                        /**< True if transition was OK */

    // Public methods

    /**
     * Generate random N points, within boundaries, separated by a minimum
     * distance
     * @param N number of vehicles
     * @param pmin lower limit of workspace
     * @param pmax upper limit of workspace
     * @param rmin minimum distance between generated points
     * @return Matrix with points following the convention
     *  [x1 x2 x3 ... xn]
     *  |y1 y2 y3 ... yn|
     *  [z1 z2 z3 ... zn]
     */
    MatrixXd gen_rand_pts(const int &N,
                          const Vector3d &pmin,
                          const Vector3d &pmax,
                          const float &rmin);

    /**
     * Generate a random permutation from a set of initial points:
     * the permutation is such that every vehicle needs to move, i.e.
     * po(i) != pf(i)
     * @param po set of initial locations, following the convention
     * @return set of final locations, following the convention
     */
    MatrixXd gen_rand_perm (const MatrixXd &po);

    /*
     * Public setters
     */

    /**
     * Verifies that the set of initial points is within the workspace and
     * non-colliding w.r.t the variable rmin.
     * If everything is consistent, write the values to private variables. If
     * some of the initial locations are out of bounds, it will relocate it
     * within the workspace.
     * @param po set of initial positions
     */
    void set_initial_pts(const MatrixXd &po);

    void set_k_factor(const int &k_factor);

    /**
     * Verifies consistency of the set of final locations (inbounds). Will
     * relocate if final position if out of bounds
     * @param pf set of final locations
     */
    void set_final_pts(const MatrixXd &pf);

    /**
     * Set workspace boundaries, copy to private variables
     * @param pmin lower limit of workspace
     * @param pmax upper limit of workspace
     */
    void set_boundaries (const Vector3d &pmin, const Vector3d &pmax);

    void set_cluster_num (const int &num);

    /**
     * Solve the DMPC algorithm sequentially for each agent of the transition.
     * It will create a transition trajectory for a fixed amount of time, and
     * after finishing it will do a sanity check to verify if the agents
     * arrived at their goals
     * @return the position, velocity and acceleration profiles for each agent
     */
    std::vector<Trajectory> solveDMPC();

    /**
     * First version of the DMPC algorithm using parallelization. The agents
     * are separated into clusters, that are solved in parallel (max of 8
     * clusters running in separate threads). This version will simulate a
     * trajectory for a fixed amount of time and do the sanity check after
     * solving
     * @return the position, velocity and acceleration profiles for each agent
     */
    std::vector<Trajectory> solveParallelDMPC();

    /**
     * Second version of the DMPC algorithm using parallelization. After each
     * time step of the trajectory, we check if the agents arrived at their
     * goal to end the execution. After solving there's still sanity checks
     * on collision violation.
     * @return the position, velocity and acceleration profiles for each agent
     */
    std::vector<Trajectory> solveParallelDMPCv2();

    /**
     * Write a .txt file to be read by MATLAB, specifying certain parameters
     * of the algorithm and the solution. Visualization and debugging purposes
     * @param src
     * @param pathAndName
     */
    void trajectories2file(const std::vector<Trajectory> &src,
                           char const* pathAndName);

private:
    // Private Variables

    /*
     * Algorithm parameters
     */

    std::string _solver_name; /**< Indicate what solver to use */
    float _h;             /**< Time step for the computation */
    int _T;               /**< Max time to complete trajectory */
    int _k_hor;           /**< Length of the prediction horizon */
    int _order;           /**< Order the ellipsoid for collision constraint */
    float _c;             /**< Multiplier for constraint in the Z direction */
    float _rmin;          /**< Min inter-agent distance for collision in XY */
    float _alim;          /**< Acceleration limit for every component XYZ */
    float _vlim;          /**< Velocity limit to scale the final solution */
    int _freq;            /**< Frequency for trajectory interpolation */
    float _goal_tol;      /**< Allowed distance error to goal */
    float _collision_tol; /**< Allowed violation of collision constraint */
    int _speed;           /**< Number of time steps for error calculation */

    float _K;             /**< Max number of time steps to complete transit*/
    float _h_scaled;      /**< Resulting time step after scaling solution */
    int _num_clusters;    /**< # of clusters to divide the problem in */

    // Workspace boundaries
    Vector3d _pmin;  /**< Lower limit of workspace [x,y,z] */
    Vector3d _pmax;  /**< Upper limit of workspace [x,y,z] */

    // Goals
    MatrixXd _po;  /**< Set of initial positions */
    MatrixXd _pf;  /**< Set of final positions */

    // Ellipsoid variables
    /*
     * Modifying a sphere to become a sphere is done by multiplying with a
     * diagonal matrix of the form:
     *     [a 0 0]
     * E = |0 b 0| ; if a = b = 1 & c > 1, you get an ellipsoid elongated in z
     *     [0 0 c]
     */
    Matrix3d _E;   /**< Scaling matrix to get an ellipsoid */
    Matrix3d _E1;  /**< E^(-1) */
    Matrix3d _E2;  /**< E^(-2) */

    // Model-related matrices
    Matrix<double, 6, 6> _A;  /**< 3D double integrator dynamics */
    Matrix<double, 6, 3> _b;  /**< 3D */
    MatrixXd _Lambda;         /**< Get position from acceleration input */
    MatrixXd _A_v;            /**< Get velocity from acceleration input */
    MatrixXd _Delta;          /**< Used for input variation computation */
    MatrixXd _A0;             /**< Propagation of initial states in position */

    int _fail; //keeps track if QP failed or not
    int _k_factor;


    //CPLEX variables
    std::vector<CPXENVptr> _env;  /**< CPLEX environment vector */
    std::vector<CPXLPptr> _lp;    /**< CPLEX problem vector */

    /*
     * Private methods
     */

    /**
     * Construct matrices Lambda and A_v
     * @param K horizon length
     */
    void get_lambda_A_v_mat(const int &K);

    /**
     * Construct Delta matrix
     * @param K horizon length
     */
    void get_delta_mat (const int &K);

    /**
     * Construct A0 matrix for initial condition propagation
     * @param K horizon length
     */
    void get_A0_mat (const int &K);

    /**
     * Initialize DMPC algorithm by specifying a straight line from po to pf
     * @param po initial location of i-th agent
     * @param pf desired final location of i-th agent
     * @return initial trajectory in position (vel and acc set to zero)
     */
    Trajectory init_dmpc (const Vector3d &po,
                          const Vector3d &pf);

    /**
     * Used in solveDMPC and solveParallelDMPC. Check if there're predicted
     * collisions with neighbours.
     * @param prev_p previous position horizon of i-th agent
     * @param obs previous position horizon of the N agents
     * @param n ID of the the i-th agent
     * @param k check collision at the k-th time step of horizon
     * @return true if any constraint is violateed
     */
    bool check_collisions(const Vector3d &prev_p,
                          const std::vector<MatrixXd> &obs,
                          const int &n, const int &k);

    /**
     * Build the collision constraint. Used with Parallel DMPC v2
     * @param prev_p previous position horizon of i-th agent
     * @param po initial location of i-th agent at the current time step
     * @param vo initial velocity of i-th agent at the current time step
     * @param obs previous position horizon of the N agents
     * @param violation_vec output og check_collisions method
     * @param n ID of the i-th agent
     * @param k k-th time step of the horizon
     * @return
     */
    std::vector<bool> check_collisionsv2(const Vector3d &prev_p,
                                         const std::vector<MatrixXd> &obs,
                                         const int &n, const int &k);

    /**
     * Build the collision constraint. Used in solveDMPC and Parallel DMPC.
     * It will include constraints with EVERY neighbour.
     * @param prev_p previous position horizon of i-th agent
     * @param po initial location of i-th agent at the current time step
     * @param vo initial velocity of i-th agent at the current time step
     * @param obs previous position horizon of the N agents
     * @param n ID of the i-th agent
     * @param k k-th time step of the horizon
     * @return pair A,b of a linear constraint Ax < b
     */
    Constraint build_collconstraint (const Vector3d &prev_p,
                                     const Vector3d &po,
                                     const Vector3d &vo,
                                     const std::vector<MatrixXd> &obs,
                                     const int &n, const int &k);

    /**
     * Build the collision constraint. Used in Parallel DMPC v2. It will
     * include constraints only with the neighbours in violation_vec
     * @param prev_p previous position horizon of i-th agent
     * @param po initial location of i-th agent at the current time step
     * @param vo initial velocity of i-th agent at the current time step
     * @param obs previous position horizon of the N agents
     * @param violation_vec output of check_collisions method
     * @param n ID of the i-th agent
     * @param k k-th time step of the horizon
     * @return pair A,b of a linear constraint Ax < b
     */
    Constraint build_collconstraintv2 (const Vector3d &prev_p,
                                       const Vector3d &po,
                                       const Vector3d &vo,
                                       const std::vector<MatrixXd> &obs,
                                       const std::vector<bool> &violation_vec,
                                       const int &n, const int &k);

    /**
     * Build and solve the QP formulation within the DMPC algorithm.
     * @param po initial location of the i-th agent at the current time step
     * @param pf final desired location of i-th agent
     * @param vo initial velocity of the i-th agent at the current time step
     * @param ao initial acceleration of the i-th agent at the current time step
     * @param n ID of the i-th agent
     * @param obs previous position horizon of the N agents
     * @return position, velocity and acceleration over horizon for i-th agent
     */
    Trajectory solveQP(const Vector3d &po, const Vector3d &pf,
                       const Vector3d &vo, const Vector3d &ao,
                       const int &n, const std::vector<MatrixXd> &obs);

    /**
     * Build and solve the QP formulation within the DMPC algorithm. Include
     * the ID of the cluster (for parallel implementation), to use the
     * appropriate CPLEX environment and problem.
     * @param po initial location of the i-th agent at the current time step
     * @param pf final desired location of i-th agent
     * @param vo initial velocity of the i-th agent at the current time step
     * @param ao initial acceleration of the i-th agent at the current time step
     * @param n ID of the i-th agent
     * @param obs previous position horizon of the N agents
     * @param id_cluster ID of the cluster agent 'n' is in
     * @return position, velocity and acceleration over horizon for i-th agent
     */
    Trajectory solveQPv2(const Vector3d &po, const Vector3d &pf,
                       const Vector3d &vo, const Vector3d &ao,
                       const int &n, const std::vector<MatrixXd> &obs,
                        const int &id_cluster);

    /**
     * This function is called within solveDMPC and solveParallelDMPC to
     * verify, at the end of the execution, if all agents arrived at their goals
     * @param all_trajectories solution vector, populated as the algorithm runs
     * @param pf desired final locations for the N agetns
     * @param goal_tol tolerance to desired location, in meters
     * @param N number of agents
     * @return true if all vehicles reached their goals
     */
    bool reached_goal(const std::vector<Trajectory> &all_trajectories,
                      const MatrixXd &pf, const float &goal_tol, const int &N);

    /**
     * Called within Parallel v2 to check, after solving each time step of
     * the transition, if the agents arrived at their desired location
     * @param all_trajectories solution vector, populated as the algorithm runs
     * @param pf desired final locations for the N agetns
     * @param goal_tol tolerance to desired location, in meters
     * @param N number of agents
     * @param k time step of the trajectory
     * @return
     */
    bool reached_goalv2(const std::vector<Trajectory> &all_trajectories,
                        const MatrixXd &pf, const float &goal_tol,
                        const int &N, const int &k);

    /**
     * After solving and interpolating the solution, verify that no
     * collisions were violated in the transition
     * @param solution interpolated solution
     * @return true if collisions will occur
     */
    bool collision_violation(const std::vector<Trajectory> &solution);

    /**
     * Interpolate the solution vector, to within the specified step size
     * @param sol solution vector for all the agents
     * @param step_size desired step size of the interpolation
     * @return interpolated solution
     */
    std::vector<Trajectory> interp_trajectory(const std::vector<Trajectory> &sol,
                                              const double &step_size);

    /**
     * Calculate the predicted time to complete the transition
     * @param solution interpolated solution
     * @return time to complete transition, in seconds
     */
    double get_trajectory_time(const std::vector<Trajectory> &solution);


    /**
     * Solve each agent cluster in a separate thread. Calls all the v1 functions
     * @param k current time step of the trajectory to solve
     * @param all_trajectories vector of solution of the transition
     * @param obs updated position horizon of the N agents
     * @param agents ID of agents within the cluster
     * @param prev_obs previous position horizon of the N agents
     */
    void cluster_solve(const int &k,
                       std::vector<Trajectory> &all_trajectories,
                       std::vector<MatrixXd> &obs,
                       const std::vector<int> &agents,
                       const std::vector<MatrixXd> &prev_obs);

    /**
     * Solve each agent cluster in a separate thread
     * @param k current time step of the trajectory to solve
     * @param all_trajectories vector of solution of the transition
     * @param obs updated position horizon of the N agents
     * @param agents ID of agents within the cluster
     * @param prev_obs previous position horizon of the N agents
     * @param id_cluster ID of the cluster being solved
     */
    void cluster_solvev2(const int &k,
                         std::vector<Trajectory> &all_trajectories,
                         std::vector<MatrixXd> &obs,
                         const std::vector<int> &agents,
                         const std::vector<MatrixXd> &prev_obs,
                         const int id_cluster);

    /**
     * Scale the resulting solution to push the velocity or acceleration limits.
     * The output will only modify the velocity and acceleration profiles,
     * while compressing in time the position profile.
     * @param sol solution vector, before interpolation
     * @param vmax maximum allowed velocity in x, y or z
     * @param amax maximum allowed acceleration in x, y, or z
     */
    void scale_solution(std::vector<Trajectory> &sol,
                        const float &vmax, const float &amax);

    /**
     * Initialize the CPLEX environment for a specified cluster ID
     * @param id ID corresponding to the cluster
     */
    void init_cplex(int id);

    /**
     * Terminate the CPLEX environment of a specified cluster
     * @param id ID corresponding to the cluster
     */
    void terminate_cplex(int id);

    /**
     * Convert eigen matrices to the CPLEX sparse format
     * @param H square matrix (either cost or constraint)
     * @param matbeg beginning index of columns
     * @param matcnt number of non-zero elements of columns
     * @param matind row number of a specified matval
     * @param matval non-zero value of matrix H
     * @param numrows rows of H
     * @param numcols columns of H
     * @param numnz non-zero elements
     */
    void eigen_to_cplex(const Eigen::MatrixXd &H, int *&matbeg,
                        int *&matcnt, int *&matind, double *&matval,
                        int &numrows, int &numcols, int &numnz);
};

#endif //DMPC_CPP_DMPC_H
