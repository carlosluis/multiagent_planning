/*
 * dmpc.cpp
 *
 *  Created On : Mar 16, 2018
 *      Author : Carlos Luis
 *      Email  : carlos.luis@robotics.utias.utoronto.ca
 */

#include <iostream>
#include "dmpc.h"

using namespace Eigen;
using namespace std;
using namespace std::chrono;

bool execution_ended;
int failed_i_global;

DMPC::DMPC(std::string solver_name,Params params)
{
    // Load parameters into private variables
    _solver_name = solver_name;
    _h = params.h;
    _T = params.T;
    _k_hor = params.k_hor;
    _order = params.order;
    _c = params.c;
    _rmin = params.rmin;
    _alim = params.alim;
    _vlim = params.vlim;
    _freq = params.freq;
    _goal_tol = params.goal_tol;
    _collision_tol = params.collision_tol;
    _speed = params.speed;
    _K = _T / _h + 1;

    successful = false;

    // Initialize srand DO THIS ONLY ONCE
    srand((unsigned int) time(0));

    // Ellipsoid definitions
    Vector3d v(3);
    v << 1, 1, _c;
    _E = v.asDiagonal();
    _E1 = _E.inverse();
    _E2 = _E1.array().pow(_order);

    // Double integrator model
    _A << 1,0,0,_h,0,0,
          0,1,0,0,_h,0,
          0,0,1,0,0,_h,
          0,0,0,1,0,0,
          0,0,0,0,1,0,
          0,0,0,0,0,1;
    double h2 = pow(_h,2);
    _b << h2/2,0,0,
          0,h2/2,0,
          0,0,h2/2,
          _h,0,0,
          0,_h,0,
          0,0,_h;

    // Get matrices for MPC computation
    get_lambda_A_v_mat(_k_hor);
    get_delta_mat(_k_hor);
    get_A0_mat(_k_hor);

    // Vicon Room boundaries
    _pmin << -1.0, -1.0, 0.2;
    _pmax << 1.0, 1.0, 2.2;

    // Default number of clusters, in case the setter is never called
    _num_clusters = 8;
}

/***************************************************************************
    *
    *  MPC Matrices
    *
    **************************************************************************/

void DMPC::get_lambda_A_v_mat(const int &K) {
    /*
     * Matrices of the form:
     *         [b           0        ...   0]
     * Lambda =|Ab          b        ...   0|
     *         |..                   ...   0|
     *         [A^(K-1)b   A^(K-2)b  ...   b]
     */
    MatrixXd Apos = MatrixXd::Zero(3 * K, 3 * K);
    MatrixXd Avel = MatrixXd::Zero(3 * K, 3 * K);
    MatrixXd prev_row = MatrixXd::Zero(6, 3 * K);
    MatrixXd new_row = MatrixXd::Zero(6, 3 * K);
    MatrixXd add_b = MatrixXd::Zero(6, 3 * K);
    int idx = 0;
    for (int k = 0; k < K; ++k) {
        add_b << MatrixXd::Zero(_b.rows(), _b.cols() * (k)), _b,
                MatrixXd::Zero(_b.rows(), _b.cols() * (K - k - 1));
        new_row = _A * prev_row + add_b;
        Apos.middleRows(idx, 3) = new_row.middleRows(0, 3);
        Avel.middleRows(idx, 3) = new_row.middleRows(3, 3);
        prev_row = new_row;
        idx += 3;
    }
    _Lambda = Apos;
    _A_v = Avel;
}

void DMPC::get_delta_mat(const int &K) {
    /*
     * Matrix of the form:
     *        [ I   0   0 ...   0   0]
     * Delta =|-I   I   0 ...   0   0|
     *        | 0  -I   I ...   0   0|
     *        | ...  ...  ...   0   0|
     *        [ 0   0   0 ...  -I   I]
     */
    MatrixXd Delta = MatrixXd::Zero(3 * K, 3 * K);
    MatrixXd new_row = MatrixXd::Zero(3, 3 * K);
    Delta.topRows(3) << MatrixXd::Identity(3, 3),
            MatrixXd::Zero(3, 3 * (K - 1));
    MatrixXd b = MatrixXd::Zero(3, 6);
    b << (-1) * MatrixXd::Identity(3, 3), MatrixXd::Identity(3, 3);
    int idx = 3;
    for (int k = 0; k < K - 1; ++k) {
        new_row << MatrixXd::Zero(3, 3 * k), b, MatrixXd::Zero(3,
                                                               3 * (K - k - 2));
        Delta.middleRows(idx, 3) = new_row;
        idx += 3;
    }

    _Delta = Delta;
}

void DMPC::get_A0_mat(const int &K) {
    /*
     * Matrix of the form:
     *
     * A0 = [A  A^2  ...  A^K]'
     *
     */

    MatrixXd A0 = MatrixXd::Zero(3 * K, 6);
    MatrixXd new_row = MatrixXd::Zero(6, 6);
    MatrixXd prev_row = MatrixXd::Identity(6, 6);
    int idx = 0;
    for (int k = 0; k < K; ++k) {
        new_row = _A * prev_row;
        A0.middleRows(idx, 3) = new_row.middleRows(0, 3);
        prev_row = new_row;
        idx += 3;
    }
    _A0 = A0;
}

/***************************************************************************
    *
    *  Initialization
    *
    **************************************************************************/

Trajectory DMPC::init_dmpc(const Vector3d &po, const Vector3d &pf) {

    Trajectory init;

    // The diference tells us the direction of movement to go from po to pf
    Vector3d diff = pf - po;

    // Matlab equivalent to t = 0:h:(K-1)*h
    VectorXd t = VectorXd::LinSpaced(_K, 0, (_K - 1) * _h);

    // Restrict solution to the first k_hor steps
    VectorXd t_hor = t.head(_k_hor);

    // Initialize all matrices at zero
    init.pos = MatrixXd::Zero(3, t_hor.size());
    init.vel = MatrixXd::Zero(3, t_hor.size());
    init.acc = MatrixXd::Zero(3, t_hor.size());

    // Define a straight line from po to pf
    for (int i = 0; i < t_hor.size(); ++i) {
        init.pos.col(i) = po + t_hor[i] * diff / 10;   //((_K-1)*_h);
    }
    return init;
}

MatrixXd DMPC::gen_rand_pts(const int &N,
                            const Vector3d &pmin,
                            const Vector3d &pmax,
                            const float &rmin) {
    MatrixXd pts = MatrixXd::Zero(3, N);
    Vector3d candidate = MatrixXd::Zero(3, 1);
    VectorXd dist;
    bool pass = false;

    // Generate first point
    pts.col(0) = pmin.array()
                 + (pmax - pmin).array() *
                   ((MatrixXd::Random(3, 1).array() + 1) / 2);

    for (int n = 1; n < N; ++n) {
        while (!pass) {
            // Candidate picked randomly within workspace boundaries
            candidate = pmin.array()
                        + (pmax - pmin).array() *
                          ((MatrixXd::Random(3, 1).array() + 1) / 2);

            // Calculate distance to every previous pts calculated
            dist = ((((pts.leftCols(n)).colwise()
                      -
                      candidate).array().square()).colwise().sum()).array().sqrt();

            // If the candidate is sufficiently separated from previous pts,
            // then we add it to the Matrix of valid pts
            for (int k = 0; k < n; ++k) {
                pass = dist[k] > rmin;
                if (!pass)
                    break;
            }
            if (pass)
                pts.col(n) = candidate.array();
        }
        pass = false;
    }
    return pts;
}

MatrixXd DMPC::gen_rand_perm(const MatrixXd &po) {
    int N = po.cols();
    std::vector<int> array;
    std::vector<int> array_aux;
    int perm[N];
    MatrixXd pf = MatrixXd::Zero(3, N);
    for (int i = 0; i < N; i++) array.push_back(i);
    array_aux = array;

    // Random permute the order, making sure that po(i) != pf(i)
    // Such that every agent has to move during the transition
    for (int i = 0; i < N; i++) {
        int j;
        array_aux.clear();
        array_aux = array;
        array_aux.erase(std::remove(array_aux.begin(), array_aux.end(), i),
                        array_aux.end());
        if (i == N - 1) {
            perm[i] = array.at(0);
        } else if (i == N - 2 && array_aux.back() == N - 1) {
            perm[i] = array_aux.back();
            array.erase(std::remove(array.begin(), array.end(), perm[i]),
                        array.end());
        } else {
            j = rand() % (N - i - 1);
            perm[i] = array_aux.at(j);
            array.erase(
                    std::remove(array.begin(), array.end(), array_aux.at(j)),
                    array.end());
        }
    }

    for (int i = 0; i < N; i++) {
        pf.col(i) = po.col(perm[i]);
    }
    return pf;
}

void DMPC::set_initial_pts(const MatrixXd &po) {
    // Verify if points are within workspace
    int N = po.cols();
    Vector3d differ;
    double dist;
    int r, c;
    _po = po;
    for (int i = 0; i < 3; ++i) {
        if (_po.row(i).maxCoeff(&r, &c) > _pmax(i)) {
            _po(i, c) = _pmax(i);
            cout << "WARNING: initial position of vehicle #" << c
                 << " is out of bounds" << endl;
            cout << "Rewriting initial position to be inside of workspace"
                 << endl;
        }

        if (_po.row(i).minCoeff(&r, &c) < _pmin(i)) {
            _po(i, c) = _pmin(i);
            cout << "WARNING: initial position of vehicle #" << c
                 << " is out of bounds" << endl;
            cout << "Rewriting initial position to be inside of workspace"
                 << endl;
        }
    }

    // After clipping (if any) check if there's no initial collisions
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            if (i != j) {
                differ = _E1 * (po.col(i) - po.col(j));
                dist = pow((differ.array().pow(_order)).sum(),1.0/_order);
                if (dist < _rmin-0.05) // we add 5cm tolerance factor to it
                {
                    cout << "Collision constraint violation: ";
                    cout << "Vehicles " << i << " and " << j;
                    cout << " will be " << dist << "m";
                    cout << " apart at their starting locations" << endl;
                }
            }
        }
    }
}

void DMPC::set_final_pts(const MatrixXd &pf) {
    // Verify if points are within workspace
    int N = pf.cols();
    Vector3d differ;
    double dist;
    int r, c;
    _pf = pf;
    for (int i = 0; i < 3; ++i) {
        if (_pf.row(i).maxCoeff(&r, &c) > _pmax(i)) {
            _pf(i, c) = _pmax(i);
            cout << "WARNING: final position of vehicle #" << c
                 << " is out of bounds" << endl;
            cout << "Rewriting final position to be inside of workspace"
                 << endl;
        }

        if (_pf.row(i).minCoeff(&r, &c) < _pmin(i)) {
            _pf(i, c) = _pmin(i);
            cout << "WARNING: initial position of vehicle #" << c
                 << " is out of bounds" << endl;
            cout << "Rewriting initial position to be inside of workspace"
                 << endl;
        }
    }

    // After clipping (if any) check if there's no initial collisions
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            if (i != j) {
                differ = _E1 * (pf.col(i) - pf.col(j));
                dist = pow((differ.array().pow(_order)).sum(),1.0/_order);
                if (dist < _rmin-0.05) // we add 5cm tolerance factor to it
                {
                    cout << "Collision constraint violation: ";
                    cout << "Vehicles " << i << " and " << j;
                    cout << " will be " << dist << "m";
                    cout << " apart at their final locations" << endl;
                    cout << "Solver won't be able to draw agents to their final locations!"<<endl;
                }
            }
        }
    }
}

void DMPC::set_boundaries(const Vector3d &pmin, const Vector3d &pmax) {
    _pmin = pmin;
    _pmax = pmax;
}

void DMPC::set_cluster_num(const int &num) {
    _num_clusters = num;
}

/***************************************************************************
    *
    *  Collision checking and constraint building
    *
    **************************************************************************/

bool DMPC::check_collisions(const Vector3d &prev_p,
                            const std::vector<MatrixXd> &obs,
                            const int &n, const int &k)
{
    bool violation = false;
    Vector3d pj;
    Vector3d diff;
    double dist;
    for (int i = 0; i < obs.size(); ++i)
    {
        if (i != n)
        {
            pj = obs[i].col(k);
            diff = _E1*(prev_p - pj);
            dist = pow(((diff.array().pow(_order)).sum()),1.0/_order);
            if (dist < _rmin)
                return violation = true;
        }
    }
    return violation;
}

std::vector<bool> DMPC::check_collisionsv2(const Vector3d &prev_p,
                                           const std::vector<MatrixXd> &obs,
                                           const int &n, const int &k)
{
    Vector3d pj;
    Vector3d diff;
    int N = obs.size();
    std::vector<bool> violation(N);
    std::vector<bool> viol_constr(N);
    std::vector<double> dist(N-1);
    int idx = 0;

    // Check the distance from i-th agent to each of the N-1 neighbours
    // Populate the viol_constr array with the conflicting neighbours
    for (int i = 0; i < N; ++i)
    {
        if (i != n)
        {
            pj = obs[i].col(k);
            diff = _E1*(prev_p - pj);
            dist[idx] = pow(((diff.array().pow(_order)).sum()),1.0/_order);
            violation[i] = (dist[idx] < _rmin);
            viol_constr[i] = (dist[idx] < _rmin*(1+(float)k/_k_hor));
            idx++;
        }
        else
        {
            violation[i] = false;
            viol_constr[i] = false;
        }
    }

    // In case no violations occured, we clear the array to check outside
    // this function its size and make a decision of what to do next
    if (std::find(violation.begin(),violation.end(),true) == violation.end())
        viol_constr.clear();

    return viol_constr;
}

Constraint DMPC::build_collconstraint(const Vector3d &prev_p,
                                      const Vector3d &po,
                                      const Vector3d &vo,
                                      const std::vector<MatrixXd> &obs,
                                      const int &n, const int &k)
{
    int N_obs = obs.size();
    Vector3d pj;
    MatrixXd Ain_total = MatrixXd::Zero(N_obs-1,3*_k_hor);
    VectorXd bin_total = MatrixXd::Zero(N_obs-1,1);
    Vector3d diff;
    double dist;
    double r;
    int idx = 0;
    Matrix <double, 6,1> initial_states;
    MatrixXd diff_row = MatrixXd::Zero(1,3*_k_hor);
    initial_states << po,vo;
    Constraint collision;

    // Build the collision constraint based on its linearization around the
    // previous solution (prev_p & pj are both past information)
    for (int i=0; i < N_obs; ++i)
    {
        if(i!=n)
        {
            pj = obs[i].col(k);
            diff = _E1*(prev_p - pj);
            dist = pow(((diff.array().pow(_order)).sum()),1.0/_order);
            diff = (_E2*(prev_p - pj)).array().pow(_order - 1);

            r = pow(dist,_order-1)*(_rmin - dist) + diff.transpose()*prev_p
                - diff.transpose()*_A0.middleRows(3*(k-1),3)*initial_states;

            diff_row << MatrixXd::Zero(1,3*(k-1)),
                        diff.transpose(),
                        MatrixXd::Zero(1,3*(_k_hor-(k-1)-1));
            Ain_total.row(idx) = -diff_row*_Lambda;
            bin_total[idx] = -r;
            idx++;
        }
    }
    collision.A = Ain_total;
    collision.b = bin_total;
    return collision;
}

Constraint DMPC::build_collconstraintv2(const Vector3d &prev_p,
                                        const Vector3d &po,
                                        const Vector3d &vo,
                                        const std::vector<MatrixXd> &obs,
                                        const std::vector<bool> &violation_vec,
                                        const int &n, const int &k)
{
    int N_obs = obs.size();
    int N_violation = std::count (violation_vec.begin(), violation_vec.end(), true);
    Vector3d pj;
    MatrixXd Ain_total = MatrixXd::Zero(N_violation,3*_k_hor);
    VectorXd bin_total = MatrixXd::Zero(N_violation,1);
    Vector3d diff;
    double dist;
    double r;
    int idx = 0;
    Matrix <double, 6,1> initial_states;
    MatrixXd diff_row = MatrixXd::Zero(1,3*_k_hor);
    initial_states << po,vo;
    Constraint collision;
    int k_ctr = k;
    collision.prev_dist = VectorXd::Zero(N_violation);

    // Build the collision constraint based on its linearization around the
    // previous solution (prev_p & pj are both past information)
    for (int i=0; i < N_obs; ++i)
    {
        if(i!=n && violation_vec[i])
        {
            pj = obs[i].col(k);
            diff = _E1*(prev_p - pj);
            dist = pow(((diff.array().pow(_order)).sum()),1.0/_order);
            diff = (_E2*(prev_p - pj)).array().pow(_order - 1);

            collision.prev_dist[idx] = pow(dist,_order-1);

            r = pow(dist,_order-1)*(_rmin - dist) + diff.transpose()*prev_p
                - diff.transpose()*_A0.middleRows(3*k_ctr,3)*initial_states;

            diff_row << MatrixXd::Zero(1,3*k_ctr),
                    diff.transpose(),
                    MatrixXd::Zero(1,3*(_k_hor-k_ctr-1));
            Ain_total.row(idx) = -diff_row*_Lambda;
            bin_total[idx] = -r;
            idx++;
        }
    }
    collision.A = Ain_total;
    collision.b = bin_total;
    return collision;
}

/***************************************************************************
    *
    *  QP building and solving
    *
    **************************************************************************/

Trajectory DMPC::solveQP(const Vector3d &po, const Vector3d &pf,
                         const Vector3d &vo, const Vector3d &ao,
                         const int &n, const std::vector<MatrixXd> &obs)
{
    int N = obs.size();          // Number of vehicles in transition
    bool violation;              // True if collision constraint is violated
    Trajectory solution;         // Pos, vel, acc profiles
    Constraint coll_constraint;  // Linear inequality constraint Ax <= b
    MatrixXd prev_p = obs.at(n); // Previous position horizon of n-th vehicle
    VectorXd a0_1;               // Vector with initial acceleration

    // Initial state joint vector (pos + vel)
    Matrix <double, 6,1> initial_states;
    initial_states << po,vo;

    // Amount of variables/ ineq constraints: no collision case
    int n_var = 3*_k_hor;                   // 3D acc vector over the horizon
    int n_ineq = 2*(3*_k_hor)+2*(3*_k_hor); // Workspace boundaries + acc limits

    // Amount of variables/ ineq constraints: collision case
    int n_var_aug = n_var + N -1;      // Add slack variable for relaxation
    int n_ineq_aug = 2*(N-1) + n_ineq; // Add collisions + slack variable < 0

    // Number of variables / ineq constraints to be used by the qp
    int qp_nvar;
    int qp_nineq;
    int qp_neq = 0;  // No equality constraints in this problem

    // Collision constraints augmented with relaxation variable
    MatrixXd collconstrA_aug = MatrixXd::Zero(2*(N-1), n_var_aug);
    VectorXd collconstrb_aug = VectorXd::Zero(2*(N-1));

    // Penalty matrices for DMPC calculation
    MatrixXd Q = MatrixXd::Zero(n_var,n_var);
    MatrixXd Q_aug = MatrixXd::Zero(n_var_aug,n_var_aug);
    MatrixXd R = MatrixXd::Zero(n_var,n_var);
    MatrixXd R_aug = MatrixXd::Zero(n_var_aug,n_var_aug);
    MatrixXd S = MatrixXd::Zero(n_var,n_var);
    MatrixXd S_aug = MatrixXd::Zero(n_var_aug,n_var_aug);

    // Penalty matrices added in case of collision violation
    VectorXd f_w = VectorXd::Zero(n_var_aug);
    MatrixXd W = MatrixXd::Zero(n_var_aug,n_var_aug);

    // Final Quadratic and linear matrices for QP solver
    MatrixXd H;
    VectorXd f;

    // Tuning factor of speed
    int spd = 1;

    // Augmented model matrices in case of collisions
    MatrixXd Lambda_aug = MatrixXd::Zero(n_var_aug,n_var_aug);
    MatrixXd Lambda_aug_in = MatrixXd::Zero(n_var,n_var_aug);
    MatrixXd Delta_aug = MatrixXd::Zero(n_var_aug,n_var_aug);
    MatrixXd A0_aug = MatrixXd::Zero(n_var_aug,6);

    // Auxiliary variables
    VectorXd init_propagation = _A0*initial_states;
    VectorXd init_propagation_aug;
    VectorXd pf_rep;
    VectorXd alim_rep = _alim*VectorXd::Ones(3*_k_hor);

    // Constraint matrices and vectors to pass to the QP solver
    MatrixXd Ain;
    VectorXd bin;

    // QP results
    VectorXd x; // complete result
    VectorXd p; // position vector solution
    VectorXd v; // velocity vector solution
    VectorXd a; // acceleration vector solution

    for (int k=0; k < _k_hor; ++k)
    {
        violation = check_collisions(prev_p.col(k),obs,n,k);
        if (violation)
        {
            coll_constraint = build_collconstraint(prev_p.col(k),po,vo,obs,n,k);
            collconstrA_aug << coll_constraint.A, MatrixXd::Identity(N-1,N-1),
                      MatrixXd::Zero(N-1,n_var), MatrixXd::Identity(N-1,N-1);
            collconstrb_aug << coll_constraint.b, MatrixXd::Zero(N-1,1);
            break;
        }
    }

    // Choose appropriate Q,S,R matrices depending on the scenario

    // Case of no collisions and farther than 1m from goal
    if (!violation && (po-pf).norm() >= 1)
    {
        Q.block(3*(_k_hor-spd),3*(_k_hor-spd),3*spd,3*spd) = 1000*MatrixXd::Identity(3*spd,3*spd);
        R = 1*MatrixXd::Identity(n_var,n_var);
        S = 10*MatrixXd::Identity(n_var,n_var);
    }

    // Case of no collisions and close to goal
    else if (!violation && (po-pf).norm() < 1)
    {
        Q.block(3*(_k_hor-spd),3*(_k_hor-spd),3*spd,3*spd) = 10000*MatrixXd::Identity(3*spd,3*spd);
        R = 1*MatrixXd::Identity(n_var,n_var);
        S = 10*MatrixXd::Identity(n_var,n_var);
    }

    // Case of collisions in the horizon
    else
    {
        Q.block(3*(_k_hor-spd),3*(_k_hor-spd),3*spd,3*spd) = 1000*MatrixXd::Identity(3*spd,3*spd);
        R = 1*MatrixXd::Identity(n_var,n_var);
        S = 100*MatrixXd::Identity(n_var,n_var);
    }

    // If there are collisions, use augmented matrices
    // required for collision constraints relaxation

    if (violation)
    {
        qp_nvar = n_var_aug;
        qp_nineq = n_ineq_aug;

        Q_aug.block(0,0,n_var,n_var) = Q;
        R_aug.block(0,0,n_var,n_var) = R;
        S_aug.block(0,0,n_var,n_var) = S;

        Lambda_aug.block(0,0,n_var,n_var) = _Lambda;
        Lambda_aug_in.block(0,0,n_var,n_var) = _Lambda;
        Delta_aug.block(0,0,n_var,n_var) = _Delta;
        A0_aug.block(0,0,n_var,6) = _A0;

        init_propagation_aug = A0_aug*initial_states;

        // Build complete inequality constraints:
        // 1) collision + relaxation variable
        // 2) workspace boundaries
        // 3) acceleration limits
        bin = VectorXd::Zero(n_ineq_aug);
        bin << collconstrb_aug,
               _pmax.replicate(_k_hor,1) - init_propagation,
               -_pmin.replicate(_k_hor,1) + init_propagation,
               alim_rep, alim_rep;

        Ain = MatrixXd::Zero(n_ineq_aug,n_var_aug);

        Ain << collconstrA_aug,
                Lambda_aug_in, -Lambda_aug_in,
                MatrixXd::Identity(n_var,n_var_aug),
                -MatrixXd::Identity(n_var,n_var_aug);

        // Build linear and quadratic cost matrices

        f_w << VectorXd::Zero(n_var),
               -pow(10,6)*VectorXd::Ones(N-1);

        //NOTE: W *NEEDS* to be different than zero, if not H has determinant 0

        W.block(n_var,n_var,N-1,N-1) = pow(10,0)*(MatrixXd::Identity(N-1,N-1));

        a0_1 = VectorXd::Zero(n_var_aug);
        a0_1 << ao, VectorXd::Zero(3*(_k_hor-1) + N - 1);

        pf_rep = VectorXd::Zero(n_var_aug);
        pf_rep << pf.replicate(_k_hor,1),MatrixXd::Zero(N-1,1);

        f = VectorXd::Zero(n_var_aug);
        H = MatrixXd::Zero(n_var_aug,n_var_aug);

        f = -2*(pf_rep.transpose()*Q_aug*Lambda_aug -
                init_propagation_aug.transpose()*Q_aug*Lambda_aug +
                a0_1.transpose()*S_aug*Delta_aug);

        f += f_w;

        H = 2*(Lambda_aug.transpose()*Q_aug*Lambda_aug
               + Delta_aug.transpose()*S_aug*Delta_aug
               + R_aug + W);
    }

    else // matrices when there're no collisions
    {
        qp_nvar = n_var;
        qp_nineq = n_ineq;

        bin = VectorXd::Zero(n_ineq);
        bin <<  _pmax.replicate(_k_hor,1) - init_propagation,
                -_pmin.replicate(_k_hor,1) + init_propagation,
                alim_rep, alim_rep;

        Ain = MatrixXd::Zero(n_ineq,n_var);

        Ain <<  _Lambda, -_Lambda,
                MatrixXd::Identity(n_var,n_var),
                -MatrixXd::Identity(n_var,n_var);

        a0_1 = VectorXd::Zero(n_var);
        a0_1 << ao, VectorXd::Zero(3*(_k_hor-1));
        pf_rep = VectorXd::Zero(n_var);
        pf_rep << pf.replicate(_k_hor,1);

        f = VectorXd::Zero(n_var);
        H = MatrixXd::Zero(n_var,n_var);

        f = -2*(pf_rep.transpose()*Q*_Lambda -
                init_propagation.transpose()*Q*_Lambda +
                a0_1.transpose()*S*_Delta);

        H = 2*(_Lambda.transpose()*Q*_Lambda
               + _Delta.transpose()*S*_Delta
               + R);
    }

    // Declare quadprog object and solve the QP
    QuadProgDense _qp(qp_nvar,qp_neq,qp_nineq);

    _qp.solve(H,f,MatrixXd::Zero(0, qp_nvar),VectorXd::Zero(0),Ain,bin);
    x = _qp.result();
    _fail = _qp.fail();
    if(_fail){
        cout << "Fail = " << _fail << endl;
        execution_ended = true;
        failed_i_global = n;
    }

    // Extract acceleration from the result
    a = x.head(n_var);

    // propagate states
    p = _Lambda*a + init_propagation;
    v = _A_v*a + vo.replicate(_k_hor,1);

    // Convert a 3*K vector into a 2D matrix of size (3,K)
    MatrixXd pos_aux(3,_k_hor);
    MatrixXd vel_aux(3,_k_hor);
    MatrixXd acc_aux(3,_k_hor);
    int idx = 0;
    for (int i = 0; i < _k_hor; ++i)
    {
        pos_aux.col(i) = p.segment(idx,3);
        vel_aux.col(i) = v.segment(idx,3);
        acc_aux.col(i) = a.segment(idx,3);
        idx += 3;
    }

    // Assign values to the trajectory solution
    solution.pos = pos_aux;
    solution.vel = vel_aux;
    solution.acc = acc_aux;
    return solution;
}

Trajectory DMPC::solveQPv2(const Vector3d &po, const Vector3d &pf,
                         const Vector3d &vo, const Vector3d &ao,
                         const int &n, const std::vector<MatrixXd> &obs,
                           const int &id_cluster)
{
    CPXENVptr env = _env.at(id_cluster); // CPLEX environment for this cluster
    CPXLPptr lp = _lp.at(id_cluster);    // CPLEX problem object

    int N = obs.size();              // Number of vehicles in transition
    int N_violation;                 // Number of conflicting neighbours
    bool violation;                  // True if collision constraint is violated
    Trajectory solution;             // Pos, vel, acc profiles
    Constraint coll_constraint;      // Linear inequality constraint Ax <= b
    MatrixXd prev_p = obs.at(n);     // Previous pos horizon of n-th vehicle
    VectorXd a0_1;                   // Vector with initial acceleration
    std::vector<bool> violation_vec; // IDs of conflicting neighbours
    int status = 0;                  // exit flag of any QP solver

    // Initial state joint vector (pos + vel)
    Matrix <double, 6,1> initial_states;
    initial_states << po,vo;

    // Amount of variables/ ineq constraints: no collision case
    int n_var = 3*_k_hor;                   // 3D acc over the horizon
    int n_ineq = 2*(3*_k_hor)+2*(3*_k_hor); // Workspace boundaries + acc limits

    // Number of variables / ineq constraints to be used by the QP
    int qp_nvar;
    int qp_nineq;
    int qp_neq = 0; // No equality constraints in this problem

    // Penalty matrices for DMPC calculation
    MatrixXd Q = MatrixXd::Zero(n_var,n_var);  // Trajectory error
    MatrixXd R = MatrixXd::Zero(n_var,n_var);  // Input magnitude
    MatrixXd S = MatrixXd::Zero(n_var,n_var);  // Input variation

    // Final Quadratic and linear matrices for QP solver
    MatrixXd H;
    VectorXd f;

    // Tuning factors
    int spd = _speed;
    float collision_tol = _collision_tol;
    int term = -1*pow(10,6);

    // Auxiliary variables
    VectorXd init_propagation = _A0*initial_states;
    VectorXd init_propagation_aug;
    VectorXd pf_rep;
    VectorXd alim_rep = _alim*VectorXd::Ones(3*_k_hor);

    // Constraint matrices and vectors to pass to the QP solver
    MatrixXd Ain;
    VectorXd bin;

    // Augmented variables
    int n_var_aug;
    int n_ineq_aug;
    MatrixXd collconstrA_aug;
    VectorXd collconstrb_aug;
    MatrixXd Q_aug;
    MatrixXd R_aug;
    MatrixXd S_aug;

    // Penalty matrices added in case of collision violation
    VectorXd f_w;
    MatrixXd W;

    // Augmented model matrices in case of collisions
    MatrixXd Lambda_aug;
    MatrixXd Lambda_aug_in;
    MatrixXd Delta_aug;
    MatrixXd A0_aug;

    // QP results
    VectorXd x;  // Complete result
    VectorXd p;  // Position vector solution
    VectorXd v;  // Velocity vector solution
    VectorXd a;  // Acceleration vector solution

    // Check for collisions in the prediction horizon
    for (int k=0; k < _k_hor; ++k)
    {
        violation_vec = check_collisionsv2(prev_p.col(k),obs,n,k);
        violation = violation_vec.size() > 0;
        if (violation) // violations occured
        {
//            if (k==0)
//                continue;
            N_violation = std::count (violation_vec.begin(),
                                      violation_vec.end(), true);

            // Amount of variables/ ineq constraints: collision case
            n_var_aug = n_var + N_violation;     // Add relaxation variables
            n_ineq_aug = 3*N_violation + n_ineq; // Add collisions + bounds

            // Collision constraints augmented with relaxation variable
            collconstrA_aug = MatrixXd::Zero(3*N_violation, n_var_aug);
            collconstrb_aug = VectorXd::Zero(3*N_violation);

            coll_constraint = build_collconstraintv2(prev_p.col(k),po,vo,obs,
                                                     violation_vec,n,k);
            MatrixXd diag_prevdist = coll_constraint.prev_dist.asDiagonal();

            collconstrA_aug << coll_constraint.A, diag_prevdist,
                               MatrixXd::Zero(N_violation,n_var),
                               MatrixXd::Identity(N_violation,N_violation),
                               MatrixXd::Zero(N_violation,n_var),
                              -MatrixXd::Identity(N_violation,N_violation);

            collconstrb_aug << coll_constraint.b, VectorXd::Zero(N_violation),
                               collision_tol*VectorXd::Ones(N_violation);
            break;
        }
    }

    // Choose appropriate Q,S,R matrices depending on the scenario

    // Case of no collisions and farther than 1m from goal
    if (!violation && (po-pf).norm() >= 1)
    {
        Q.block(3*(_k_hor-spd),3*(_k_hor-spd),3*spd,3*spd) =
                1000*MatrixXd::Identity(3*spd,3*spd);
        R = 1*MatrixXd::Identity(n_var,n_var);
        S = 10*MatrixXd::Identity(n_var,n_var);
    }

        // Case of no collisions and close to goal
    else if (!violation && (po-pf).norm() < 1)
    {
        Q.block(3*(_k_hor-spd),3*(_k_hor-spd),3*spd,3*spd) =
                10000*MatrixXd::Identity(3*spd,3*spd);
        R = 1*MatrixXd::Identity(n_var,n_var);
        S = 10*MatrixXd::Identity(n_var,n_var);
    }

        // Case of collisions in the horizon
    else
    {
        Q.block(3*(_k_hor-spd),3*(_k_hor-spd),3*spd,3*spd) =
                1000*MatrixXd::Identity(3*spd,3*spd);
        R = 1*MatrixXd::Identity(n_var,n_var);
        S = 100*MatrixXd::Identity(n_var,n_var);
    }

    // If there are collisions, use augmented matrices
    // required for collision constraints relaxation

    if (violation)
    {
        qp_nvar = n_var_aug;
        qp_nineq = n_ineq_aug;

        Q_aug = MatrixXd::Zero(n_var_aug,n_var_aug);
        R_aug = MatrixXd::Zero(n_var_aug,n_var_aug);
        S_aug = MatrixXd::Zero(n_var_aug,n_var_aug);

        // Penalty matrices added in case of collision violation
        f_w = VectorXd::Zero(n_var_aug);
        W = MatrixXd::Zero(n_var_aug,n_var_aug);

        // Augmented model matrices in case of collisions
        Lambda_aug = MatrixXd::Zero(n_var_aug,n_var_aug);
        Lambda_aug_in = MatrixXd::Zero(n_var,n_var_aug);
        Delta_aug = MatrixXd::Zero(n_var_aug,n_var_aug);
        A0_aug = MatrixXd::Zero(n_var_aug,6);

        Q_aug.block(0,0,n_var,n_var) = Q;
        R_aug.block(0,0,n_var,n_var) = R;
        S_aug.block(0,0,n_var,n_var) = S;

        Lambda_aug.block(0,0,n_var,n_var) = _Lambda;
        Lambda_aug_in.block(0,0,n_var,n_var) = _Lambda;
        Delta_aug.block(0,0,n_var,n_var) = _Delta;
        A0_aug.block(0,0,n_var,6) = _A0;

        init_propagation_aug = A0_aug*initial_states;

        // Build complete inequality constraints:
        // 1) collision + relaxation variable
        // 2) workspace boundaries
        // 3) acceleration limits
        bin = VectorXd::Zero(n_ineq_aug);
        bin <<  collconstrb_aug,
                _pmax.replicate(_k_hor,1) - init_propagation,
                -_pmin.replicate(_k_hor,1) + init_propagation,
                alim_rep, alim_rep;

        Ain = MatrixXd::Zero(n_ineq_aug,n_var_aug);

        Ain << collconstrA_aug,
                Lambda_aug_in, -Lambda_aug_in,
                MatrixXd::Identity(n_var,n_var_aug),
                -MatrixXd::Identity(n_var,n_var_aug);

        // Build linear and quadratic cost matrices
        f_w << VectorXd::Zero(n_var),
                term * VectorXd::Ones(N_violation);
//                term*coll_constraint.prev_dist.array().inverse() ;

        // W *NEEDS* to be different than zero, if not H has determinant 0

        W.block(n_var,n_var,N_violation,N_violation) =
                pow(10,0)*(MatrixXd::Identity(N_violation,N_violation));

        a0_1 = VectorXd::Zero(n_var_aug);
        a0_1 << ao, VectorXd::Zero(3*(_k_hor-1) + N_violation);

        pf_rep = VectorXd::Zero(n_var_aug);
        pf_rep << pf.replicate(_k_hor,1),MatrixXd::Zero(N_violation,1);

        f = VectorXd::Zero(n_var_aug);
        H = MatrixXd::Zero(n_var_aug,n_var_aug);

        f = -2*(pf_rep.transpose()*Q_aug*Lambda_aug -
                init_propagation_aug.transpose()*Q_aug*Lambda_aug +
                a0_1.transpose()*S_aug*Delta_aug);

        f += f_w;

        H = 2*(Lambda_aug.transpose()*Q_aug*Lambda_aug
               + Delta_aug.transpose()*S_aug*Delta_aug
               + R_aug + W);
    }

    else // matrices when there're no collisions
    {
        qp_nvar = n_var;
        qp_nineq = n_ineq;

        bin = VectorXd::Zero(n_ineq);
        bin <<  _pmax.replicate(_k_hor,1) - init_propagation,
                -_pmin.replicate(_k_hor,1) + init_propagation,
                alim_rep, alim_rep;

        Ain = MatrixXd::Zero(n_ineq,n_var);

        Ain <<  _Lambda, -_Lambda,
                MatrixXd::Identity(n_var,n_var),
                -MatrixXd::Identity(n_var,n_var);

        a0_1 = VectorXd::Zero(n_var);
        a0_1 << ao, VectorXd::Zero(3*(_k_hor-1));
        pf_rep = VectorXd::Zero(n_var);
        pf_rep << pf.replicate(_k_hor,1);

        f = VectorXd::Zero(n_var);
        H = MatrixXd::Zero(n_var,n_var);

        f = -2*(pf_rep.transpose()*Q*_Lambda -
                init_propagation.transpose()*Q*_Lambda +
                a0_1.transpose()*S*_Delta);

        H = 2*(_Lambda.transpose()*Q*_Lambda
               + _Delta.transpose()*S*_Delta
               + R);
    }

    /*
     * Differentiate between the different possible solvers to use
     */

    if (!strcmp(_solver_name.c_str(),"quadprog")) // Eigen-quadprog solver
    {
        // Declare quadprog object and solve the QP
        QuadProgDense qp(qp_nvar,qp_neq,qp_nineq);
        MatrixXd H_sym = MatrixXd::Zero(qp_nvar,qp_nvar);
        H_sym = 0.5*(H+H.transpose());
        qp.solve(H_sym,f,MatrixXd::Zero(0, qp_nvar),VectorXd::Zero(0),Ain,bin);
        x = qp.result();

        status = qp.fail();
        int tries = 0;
        float lim = 0.05;

        while (status && tries < 20){
            lim = 2*lim;
            term = 2*term;
            // Debug print
//            cout << "Infeasible - Retrying..." << endl;
            collconstrb_aug << coll_constraint.b, VectorXd::Zero(N_violation),
                               lim*VectorXd::Ones(N_violation);

            bin << collconstrb_aug,
                    _pmax.replicate(_k_hor,1) - init_propagation,
                    -_pmin.replicate(_k_hor,1) + init_propagation,
                    alim_rep, alim_rep;

            f_w << VectorXd::Zero(n_var),
//                    term * VectorXd::Ones(N_violation);
                    term*coll_constraint.prev_dist.array().inverse() ;

            f = -2*(pf_rep.transpose()*Q_aug*Lambda_aug -
                    init_propagation_aug.transpose()*Q_aug*Lambda_aug +
                    a0_1.transpose()*S_aug*Delta_aug);

            f += f_w;
            qp.solve(H_sym,f,MatrixXd::Zero(0, qp_nvar),VectorXd::Zero(0),Ain,
                     bin);
            x = qp.result();
            status = qp.fail();
            tries++;
        }
    }

    else if (!strcmp(_solver_name.c_str(),"ooqp"))
    {
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        MatrixXd H_sym = MatrixXd::Zero(qp_nvar,qp_nvar);
        H_sym = 0.5*(H+H.transpose());
        status = !ooqpei::OoqpEigenInterface::solve(H_sym.sparseView(),
                                                    f, Ain.sparseView(),
                                                    bin, x,true);
        int tries = 0;
        float lim = 0.05;
        while (status && tries < 20) {
            lim = 2 * lim;
            term = 2 * term;
            // Debug print
//            cout << "Infeasible - Retrying..." << endl;
            collconstrb_aug << coll_constraint.b, VectorXd::Zero(N_violation),
                               lim * VectorXd::Ones(N_violation);

            bin << collconstrb_aug,
                    _pmax.replicate(_k_hor, 1) - init_propagation,
                    -_pmin.replicate(_k_hor, 1) + init_propagation,
                    alim_rep, alim_rep;

            f_w << VectorXd::Zero(n_var),
                    term * VectorXd::Ones(N_violation);
//                    term*coll_constraint.prev_dist.array().inverse() ;

            f = -2 * (pf_rep.transpose() * Q_aug * Lambda_aug -
                      init_propagation_aug.transpose() * Q_aug * Lambda_aug +
                      a0_1.transpose() * S_aug * Delta_aug);

            f += f_w;
            status = !ooqpei::OoqpEigenInterface::solve(H_sym.sparseView(), f,
                                                        Ain.sparseView(),
                                                        bin, x);
            tries++;
        }

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>( t2 - t1 ).count();
//        cout << "QP time OOQP = "
//             << duration/1000000.0 << "s" << endl;
    }

    else if (!strcmp(_solver_name.c_str(),"cplex"))
    {

        CPXsetdblparam(env,CPX_PARAM_BAREPCOMP, 1e-8);
//        CPXsetintparam(env,CPX_PARAM_BARITLIM,1000);
//        CPXsetintparam(env,CPX_PARAM_BARORDER, 2);
//        CPXsetintparam(env,CPXPARAM_QPMethod,CPX_ALG_DUAL);
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        int lpstat;
        double objval;
        MatrixXd H_sym = MatrixXd::Zero(qp_nvar,qp_nvar);
        H_sym = 0.5*(H+H.transpose());

        // Variables for the sparse rep. of the quadratic cost term, H
        int hnumrows,hnumcols,hnumnz;
        int *hmatbeg,*hmatcnt,*hmatind;
        double *hmatval;

        // Variables for the sparse rep. of the constraint matrix Ain
        int tnumrows,tnumcols,tnumnz;
        int *tmatbeg,*tmatcnt,*tmatind;
        double *tmatval;

        // Convert matrices to CPLEX representation
        eigen_to_cplex(Ain, tmatbeg, tmatcnt, tmatind, tmatval,
                       tnumrows, tnumcols, tnumnz);
        eigen_to_cplex(H_sym, hmatbeg, hmatcnt, hmatind, hmatval,
                       hnumrows, hnumcols, hnumnz);

        // IBM information on how to use the functions here
        // http://www-01.ibm.com/support/knowledgecenter/SSSA5P_12.4.0/
        // ilog.odms.cplex.help/refcallablelibrary/html/functions/CPXcopylp.html
        double * objective = new double[f.size()];
        double * lb = new double[tnumcols];
        double * ub = new double[tnumcols];
        double * rhs = new double[tnumrows];
        double * sol = new double[tnumcols];
        char * sense = new char[tnumrows];

        // Array of length at least numcols with objective function coefficients
        for (uint i = 0; i < f.rows(); i++){
            objective[i] = f(i);
        }
        // Lower and upper bounds on each of the variables
        for (uint i = 0; i < f.rows(); i++) {
            lb[i] = -CPX_INFBOUND;
            ub[i] = CPX_INFBOUND;
        }
        // Inequality constraints
        for (uint i = 0; i < bin.size(); i++) {
            rhs[i] = bin(i);
            sense[i] = 'L';
        }

        // Copy linear problem, constraint matrix. specify minimization problem
        CPXcopylp(env, lp, tnumcols, tnumrows, CPX_MIN, objective, rhs, sense,
                  tmatbeg, tmatcnt, tmatind, tmatval, lb, ub, NULL);
        if (status) {
            printf("CPXcopylp failed.\n");
            terminate_cplex(id_cluster);
        }

        // Copy quadratic problem, H quadratic cost function matrix
        CPXcopyquad(env, lp, hmatbeg, hmatcnt, hmatind, hmatval);
        if (status) {
            printf("Failed to copy quadratic matrix. \n");
            terminate_cplex(id_cluster);
        }

        // Solve optimization
        status = CPXqpopt(env,lp);

        if (status) {
            printf("Failed to optimize QP: status= %i \n", status);
            terminate_cplex(id_cluster);
        }

        status = CPXsolution(env, lp, &lpstat, &objval, sol, NULL, NULL, NULL);
        // Check if solution was successfully returned
        if (status) {
            printf("Failed to retreive CPXsolution.\n");
            terminate_cplex(id_cluster);
        }

        x = VectorXd::Zero(qp_nvar);
        // Copy optimized solution to learned input vector xl
        for (uint i = 0; i < qp_nvar; i++) {
            x(i) = sol[i];
        }

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>( t2 - t1 ).count();
//        cout << "QP time CPLEX = "
//             << duration/1000000.0 << "s" << endl;
    }

    if(status){
        cout << "Fail completely" << endl;
        execution_ended = true;
        failed_i_global = n;
    }

    // Extract acceleration from the result
    a = x.head(n_var);

    // Propagate states
    p = _Lambda*a + init_propagation;
    v = _A_v*a + vo.replicate(_k_hor,1);

    // Convert a 3*K vector into a 2D matrix of size (3,K)
    MatrixXd pos_aux(3,_k_hor);
    MatrixXd vel_aux(3,_k_hor);
    MatrixXd acc_aux(3,_k_hor);
    int idx = 0;
    for (int i = 0; i < _k_hor; ++i)
    {
        pos_aux.col(i) = p.segment(idx,3);
        vel_aux.col(i) = v.segment(idx,3);
        acc_aux.col(i) = a.segment(idx,3);
        idx += 3;
    }

    // Assign values to the trajectory solution
    solution.pos = pos_aux;
    solution.vel = vel_aux;
    solution.acc = acc_aux;
//    cout << solution.pos << endl << endl;
//    cout << solution.vel << endl << endl;
//    cout << solution.acc << endl << endl;
    return solution;
}

/***************************************************************************
    *
    *  DMPC algorithm (all versions)
    *
    **************************************************************************/

std::vector<Trajectory> DMPC::solveDMPC()
{
    int N = _po.cols(); // Number of agents = number of rows of po
    int N_cmd = _pf.cols(); // Number of agents of transition = number of
    // rows of pf

    // Variables
    std::vector<Trajectory> all_trajectories(N_cmd);
    std::vector<MatrixXd> aux_trajectories(N-N_cmd);
    std::vector<Trajectory> solution(N_cmd);
    Vector3d poi;
    Vector3d pfi;
    Trajectory trajectory_i;
    Vector3d pok;
    Vector3d vok;
    Vector3d aok;
    int failed_i;
    bool arrived = false;

    high_resolution_clock::time_point t1;
    high_resolution_clock::time_point t2;

    std::vector<MatrixXd> prev_obs;
    std::vector<MatrixXd> obs;

    if (N_cmd < N)
    {
        VectorXd rep = VectorXd::Zero(3*_k_hor);
        MatrixXd pos_aux(3,_k_hor);

        for (int i=0; i < N-N_cmd; ++i)
        {
            rep = (_po.col(i+N_cmd)).replicate(_k_hor,1);
            int idx = 0;
            for (int l = 0; l < _k_hor; ++l)
            {
                pos_aux.col(l) = rep.segment(idx,3);
                idx += 3;
            }
            aux_trajectories.at(i) = pos_aux;
        }
    }

    t1 = high_resolution_clock::now();

    // Generate trajectory for each time step of trajectory, for each agent
    for (int k=0; k < _K; ++k)
    {
        obs.clear();
        for (int i=0; i < N; ++i)
        {
            if (i < N_cmd)
            {
                if (k == 0)
                {   // Initialize
                    all_trajectories.at(i).pos = MatrixXd::Zero(3,_K);
                    all_trajectories.at(i).vel = MatrixXd::Zero(3,_K);
                    all_trajectories.at(i).acc = MatrixXd::Zero(3,_K);
                    trajectory_i.pos = MatrixXd::Zero(3,_k_hor);
                    trajectory_i.vel = MatrixXd::Zero(3,_k_hor);
                    trajectory_i.acc = MatrixXd::Zero(3,_k_hor);
                    poi = _po.col(i);
                    pfi = _pf.col(i);
                    trajectory_i = init_dmpc(poi,pfi);
                    _fail = 0;
                }
                else
                {   // Update previous solution, solve current QP
//                t1 = high_resolution_clock::now();
                    pok = all_trajectories.at(i).pos.col(k-1);
                    vok = all_trajectories.at(i).vel.col(k-1);
                    aok = all_trajectories.at(i).acc.col(k-1);
                    trajectory_i = solveQP(pok,_pf.col(i),vok,aok,i,prev_obs);
//                t2 = high_resolution_clock::now();
//                auto duration = duration_cast<microseconds>( t2 - t1 ).count();
//                cout << "Computation time of solveQP = "
//                     << duration/1000000.0 << "s" << endl;
                }

                if (_fail)
                {
                    failed_i = i;
                    break;
                }

                // If QP didn't fail, then update the solution vector and obstacle list
                obs.push_back(trajectory_i.pos);
                all_trajectories.at(i).pos.col(k) = trajectory_i.pos.col(0);
                all_trajectories.at(i).vel.col(k) = trajectory_i.vel.col(0);
                all_trajectories.at(i).acc.col(k) = trajectory_i.acc.col(0);
            }

            else // uncommanded vehicles act as static obstacles
                obs.push_back(aux_trajectories.at(i-N_cmd));
        }

        if (_fail)
        {
            cout << "Failed - problem unfeasible @ k_T = " << k
                 << ", vehicle #" << failed_i << endl;
            break;
        }
        prev_obs = obs;
    }

    solution_short = all_trajectories;

    t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Solve all QPs computation time = " << duration/1000000.0 << "s"
         << endl;

    if(!_fail)
    {
        cout << "Optimization problem feasible: solution found" << endl;
        // Check if every agent reached its goal
        arrived = reached_goal(all_trajectories,_pf,0.05,N_cmd);
    }

    t1 = high_resolution_clock::now();

    if (arrived && !_fail)
    {
        cout << "All vehicles reached their goals" << endl;
        // Interpolate for better resolution (e.g. 100 Hz)
        solution = interp_trajectory(all_trajectories,1.0 / _freq);

        // Check if collision constraints were not violated
        bool violation = collision_violation(solution);

        // Calculate minimum time to complete trajectory, within 5cm of goals
        double time = get_trajectory_time(solution);
    }

    t2 = high_resolution_clock::now();
    duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Post-checks computation time = " << duration/1000000.0 << "s" << endl;
    return solution;
}

std::vector<Trajectory> DMPC::solveParallelDMPC()
{
    int N = _po.cols(); // Number of agents = number of rows of po
    int N_cmd = _pf.cols(); // Number of agents of transition = number of
    // rows of pf
    execution_ended = false; // reset variable before solving

    // Variables
    std::vector<Trajectory> all_trajectories(N_cmd);
    std::vector<MatrixXd> aux_trajectories(N-N_cmd);
    std::vector<Trajectory> solution(N_cmd);
    Vector3d poi;
    Vector3d pfi;
    Trajectory trajectory_i;
    Vector3d pok;
    Vector3d vok;
    Vector3d aok;
    bool arrived = false;

    // Separate N agents into 4 clusters to be solved in parallel
    int n_cluster = 8;
    if (n_cluster > N_cmd)
        n_cluster = N_cmd;
    std::vector<int> all_idx(n_cluster);
    std::vector<thread> all_threads(n_cluster);
    std::vector<std::vector<int>> all_clusters(n_cluster);

    int agentsXcluster = N_cmd/n_cluster;
    int residue = N_cmd%n_cluster;
    int cluster_capacity;
    int N_index = 0;

    for(int i=0; i<n_cluster; ++i)
    {
        if (residue!=0){
            cluster_capacity = agentsXcluster + 1;
            residue--;
        }
        else cluster_capacity = agentsXcluster;

        int curr_index = N_index;
        for(int j=curr_index; j<curr_index+cluster_capacity; ++j)
        {
            all_clusters.at(i).push_back(j);
            N_index = j+1;
        }
    }

    high_resolution_clock::time_point t1;
    high_resolution_clock::time_point t2;

    std::vector<MatrixXd> prev_obs(N);
    std::vector<MatrixXd> obs(N);

    if (N_cmd < N)
    {
        VectorXd rep = VectorXd::Zero(3*_k_hor);
        MatrixXd pos_aux(3,_k_hor);

        for (int i=N_cmd-1; i < N; ++i)
        {
            rep = (_po.col(i)).replicate(_k_hor,1);
            int idx = 0;
            for (int l = 0; l < _k_hor; ++l)
            {
                pos_aux.col(l) = rep.segment(idx,3);
                idx += 3;
            }
            obs.at(i) = pos_aux;
        }
    }

    t1 = high_resolution_clock::now();
    // Generate trajectory for each time step of trajectory, for each agent
    for (int k=0; k < _K; ++k)
    {
        // this is the loop that we want to parallelize
        for (int i=0; i<all_clusters.size(); ++i)
        {
            all_threads.at(i) = std::thread{&DMPC::cluster_solve, *this,
                                            ref(k),
                                            ref(all_trajectories),
                                            ref(obs),
                                            ref(all_clusters.at(i)),
                                            ref(prev_obs)};
        }
        for (int i=0; i<all_clusters.size(); ++i)
            all_threads.at(i).join();

        if (execution_ended)
        {
            cout << "Failed - problem unfeasible @ k_T = " << k
                 << ", vehicle #" << failed_i_global << endl;
            break;
        }

        prev_obs = obs;
    }

    solution_short = all_trajectories;
    t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Solve all QPs computation time = " << duration/1000000.0 << "s"
         << endl;

    if(!execution_ended)
    {
        cout << "Optimization problem feasible: solution found" << endl;
        // Check if every agent reached its goal
        arrived = reached_goal(all_trajectories,_pf,0.05,N_cmd);
    }

    t1 = high_resolution_clock::now();

    if (arrived && !execution_ended)
    {
        cout << "All vehicles reached their goals" << endl;

        scale_solution(solution_short,2.0,3.0);
        // Interpolate for better resolution (e.g. 100 Hz)
        solution = interp_trajectory(all_trajectories,0.0067);

        // Check if collision constraints were not violated
        bool violation = collision_violation(solution);

        // Calculate minimum time to complete trajectory, within 5cm of goals
        double time = get_trajectory_time(solution);
    }

    t2 = high_resolution_clock::now();
    duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Post-checks computation time = " << duration/1000000.0 << "s" << endl;
    return solution;
}

std::vector<Trajectory> DMPC::solveParallelDMPCv2()
{
    int N = _po.cols();     // Number of agents = number of rows of po
    int N_cmd = _pf.cols(); // Number of agents of transition = number of
                            // rows of pf
    execution_ended = false; // Reset variable before solving
    successful = false;      // Reset variable before solving

//    cout << "I'm solving the problem using " << _solver_name << endl;

    // Variables
    std::vector<Trajectory> all_trajectories(N_cmd);
    std::vector<MatrixXd> aux_trajectories(N-N_cmd);
    std::vector<Trajectory> solution(N_cmd);
    Vector3d poi;
    Vector3d pfi;
    Trajectory trajectory_i;
    Vector3d pok;
    Vector3d vok;
    Vector3d aok;
    bool arrived = false;

    // Separate N agents into clusters to be solved in parallel
    int n_cluster = _num_clusters;
    if (n_cluster > N_cmd)
        n_cluster = N_cmd;
    std::vector<int> all_idx(n_cluster);
    std::vector<thread> all_threads(n_cluster);
    std::vector<std::vector<int>> all_clusters(n_cluster);

    int agentsXcluster = N_cmd/n_cluster;
    int residue = N_cmd%n_cluster;
    int cluster_capacity;
    int N_index = 0;

    // initialize n_cluster CPLEX environments
    std::vector<CPXENVptr> env(n_cluster);
    std::vector<CPXLPptr> lp(n_cluster);
    _env = env;
    _lp = lp;
    for(int i=0; i<n_cluster; ++i)
    {
        init_cplex(i);
        if (residue!=0){
            cluster_capacity = agentsXcluster + 1;
            residue--;
        }
        else cluster_capacity = agentsXcluster;

        int curr_index = N_index;
        for(int j=curr_index; j<curr_index+cluster_capacity; ++j)
        {
            all_clusters.at(i).push_back(j);
            N_index = j+1;
        }
    }

    high_resolution_clock::time_point t1;
    high_resolution_clock::time_point t2;

    std::vector<MatrixXd> prev_obs(N);
    std::vector<MatrixXd> obs(N);

    if (N_cmd < N)
    {
        VectorXd rep = VectorXd::Zero(3*_k_hor);
        MatrixXd pos_aux(3,_k_hor);

        for (int i=N_cmd-1; i < N; ++i)
        {
            rep = (_po.col(i)).replicate(_k_hor,1);
            int idx = 0;
            for (int l = 0; l < _k_hor; ++l)
            {
                pos_aux.col(l) = rep.segment(idx,3);
                idx += 3;
            }
            obs.at(i) = pos_aux;
        }
    }

    t1 = high_resolution_clock::now(); // tic

    // Generate trajectory for each time step of trajectory, for each agent
    // Iterate until transition is completed or max steps is exceeded
    int k = 0;
    while (!arrived && k < _K )
    {
        // This is the loop that we want to parallelize
        for (int i=0; i<all_clusters.size(); ++i)
        {
            int id = i;
            all_threads.at(i) = std::thread{&DMPC::cluster_solvev2, *this,
                                            ref(k),
                                            ref(all_trajectories),
                                            ref(obs),
                                            ref(all_clusters.at(i)),
                                            ref(prev_obs),
                                            id};

        }
        for (int i=0; i<all_clusters.size(); ++i)
            all_threads.at(i).join();

        if (execution_ended)
        {
//            cout << "Failed - problem unfeasible @ k_T = " << k
//                 << ", vehicle #" << failed_i_global << endl;
            break;
        }

        prev_obs = obs;

        // check if the drones arrived at their goal
        arrived = reached_goalv2(all_trajectories,_pf,_goal_tol,N_cmd,k);
        k = k + 1;
    }

    // Routine to trim the resulting trajectory, which is a vector of size of
    // the maximum allowed time. So we trim it to k time steps only
    solution_short.clear(); // Just in case the variable was written before
    for (int i =0; i < N_cmd; ++i)
    {
        Trajectory aux;
        aux.pos = all_trajectories.at(i).pos.leftCols(k-1);
        aux.vel = all_trajectories.at(i).vel.leftCols(k-1);
        aux.acc = all_trajectories.at(i).acc.leftCols(k-1);
        solution_short.push_back(aux);
    }

    t2 = high_resolution_clock::now(); // toc
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
//    cout << "Solve all QPs computation time = " << duration/1000000.0 << "s"
//         << endl;

    // Sanity post-checks
    if(!execution_ended)
    {
//        cout << "Optimization problem feasible: solution found" << endl;
    }

    t1 = high_resolution_clock::now(); // tic

//    if (!execution_ended && !arrived)
//        cout << "The vehicles cannot finish the transition within the maximum allowed time"
//             << endl;

    if (arrived && !execution_ended)
    {
//        cout << "All vehicles reached their goals" << endl;
        // Scale solution to reach velocity and acceleration limits
        scale_solution(solution_short,_vlim,_alim);

        // Interpolate for better resolution (e.g. 100 Hz)
        solution = interp_trajectory(solution_short,1.0 / _freq);

        // Check if collision constraints were not violated
        bool violation = collision_violation(solution);
        successful = !violation;
//        cout << "Successful = " << successful << endl;

        // Calculate minimum time to complete trajectory, within 5cm of goals
        double time = get_trajectory_time(solution);
    }

    t2 = high_resolution_clock::now(); // toc
    duration = duration_cast<microseconds>( t2 - t1 ).count();
//    cout << "Post-checks computation time = "
//         << duration/1000000.0 << "s" << endl;
    return solution;
}

void DMPC::cluster_solve(const int &k,
                   std::vector<Trajectory> &all_trajectories,
                   std::vector<MatrixXd> &obs,
                   const std::vector<int> &agents,
                   const std::vector<MatrixXd> &prev_obs)
{
    Vector3d poi;
    Vector3d pfi;
    Trajectory trajectory_i;
    Vector3d pok;
    Vector3d vok;
    Vector3d aok;

    for (int i=agents.front(); i <= agents.back(); ++i)
    {
        if (k == 0)
        {   // Initialize
            all_trajectories.at(i).pos = MatrixXd::Zero(3,_K);
            all_trajectories.at(i).vel = MatrixXd::Zero(3,_K);
            all_trajectories.at(i).acc = MatrixXd::Zero(3,_K);
            obs.at(i) = MatrixXd::Zero(3,_k_hor);
            trajectory_i.pos = MatrixXd::Zero(3,_k_hor);
            trajectory_i.vel = MatrixXd::Zero(3,_k_hor);
            trajectory_i.acc = MatrixXd::Zero(3,_k_hor);
            poi = _po.col(i);
            pfi = _pf.col(i);
            trajectory_i = init_dmpc(poi,pfi);
            _fail = 0;
        }
        else
        {   // Update previous solution, solve current QP
            pok = all_trajectories.at(i).pos.col(k-1);
            vok = all_trajectories.at(i).vel.col(k-1);
            aok = all_trajectories.at(i).acc.col(k-1);
            trajectory_i = solveQP(pok,_pf.col(i),vok,aok,i,prev_obs);
        }

        if (_fail)
        {
            break;
        }

        // If QP didn't fail, then update the solution vector and obstacle list
        obs.at(i) = trajectory_i.pos;
        all_trajectories.at(i).pos.col(k) = trajectory_i.pos.col(0);
        all_trajectories.at(i).vel.col(k) = trajectory_i.vel.col(0);
        all_trajectories.at(i).acc.col(k) = trajectory_i.acc.col(0);
    }
}

void DMPC::cluster_solvev2(const int &k,
                           std::vector<Trajectory> &all_trajectories,
                           std::vector<MatrixXd> &obs,
                           const std::vector<int> &agents,
                           const std::vector<MatrixXd> &prev_obs,
                           const int id_cluster)
{
    Vector3d poi;
    Vector3d pfi;
    Trajectory trajectory_i;
    Vector3d pok;
    Vector3d vok;
    Vector3d aok;

    for (int i=agents.front(); i <= agents.back(); ++i)
    {
        if (k == 0)
        {   // Initialize
            all_trajectories.at(i).pos = MatrixXd::Zero(3,_K);
            all_trajectories.at(i).vel = MatrixXd::Zero(3,_K);
            all_trajectories.at(i).acc = MatrixXd::Zero(3,_K);
            obs.at(i) = MatrixXd::Zero(3,_k_hor);
            trajectory_i.pos = MatrixXd::Zero(3,_k_hor);
            trajectory_i.vel = MatrixXd::Zero(3,_k_hor);
            trajectory_i.acc = MatrixXd::Zero(3,_k_hor);
            poi = _po.col(i);
            pfi = _pf.col(i);
            trajectory_i = init_dmpc(poi,pfi);
            _fail = 0;
        }
        else
        {   // Update previous solution, solve current QP
            pok = all_trajectories.at(i).pos.col(k-1);
            vok = all_trajectories.at(i).vel.col(k-1);
            aok = all_trajectories.at(i).acc.col(k-1);
            trajectory_i = solveQPv2(pok,_pf.col(i),vok,aok,i,prev_obs, id_cluster);
        }

        if (execution_ended)
        {
            break;
        }

        // If QP didn't fail, then update the solution vector and obstacle list
        obs.at(i) = trajectory_i.pos;
        all_trajectories.at(i).pos.col(k) = trajectory_i.pos.col(0);
        all_trajectories.at(i).vel.col(k) = trajectory_i.vel.col(0);
        all_trajectories.at(i).acc.col(k) = trajectory_i.acc.col(0);
    }
}

/***************************************************************************
    *
    *  Post-solve checks, scaling and recording of solution
    *
    **************************************************************************/

bool DMPC::reached_goal(const std::vector<Trajectory> &all_trajectories,
                        const MatrixXd &pf, const float &goal_tol,
                        const int &N)
{   Vector3d diff;
    double dist;
    bool reached = true;
    for (int i=0; i < N; ++i)
    {
        diff = all_trajectories.at(i).pos.col(_K-1) - pf.col(i);
        dist = pow(((diff.array().pow(2)).sum()),1.0/2);
        if (dist > goal_tol){
            cout << "Vehicle " << i << " did not reach its goal by "
                 << dist << "m" << endl;
            reached = false;
        }
    }
    return reached;
}

bool DMPC::reached_goalv2(const std::vector<Trajectory> &all_trajectories,
                          const MatrixXd &pf, const float &goal_tol,
                          const int &N, const int &k)
{   Vector3d diff;
    double dist;
    bool reached = true;
    for (int i=0; i < N; ++i)
    {
        diff = all_trajectories.at(i).pos.col(k) - pf.col(i);
        dist = pow(((diff.array().pow(2)).sum()),1.0/2);
        if (dist > goal_tol)
            reached = false;
    }
    return reached;
}

double DMPC::get_trajectory_time(const std::vector<Trajectory> &solution)
{
    MatrixXd diff;
    int N = solution.size();
    int length_t = solution.at(0).pos.cols();
    VectorXd dist;
    VectorXd time(N);
    double min_traj_time = 0.0;

    for (int i=0; i < N; ++i)
    {
        diff = solution.at(i).pos - (_pf.col(i).replicate(1,length_t));
        dist = pow(((diff.array().pow(_order)).colwise().sum()),1.0/2);
        // Go backwards and find the last element greater than 5 cm
        for (int k = length_t-1; k >= 0 ; --k)
        {
            if(dist[k] >= 0.05)
            {
                time[i] = (float)k/100;
                break;
            }
        }
    }
    min_traj_time = time.maxCoeff();
//    cout << "The transition will be completed in " << min_traj_time << "s" << endl;
    return min_traj_time;
}

void DMPC::scale_solution(std::vector<Trajectory> &sol,
                          const float &vmax, const float &amax)
{
    int N = sol.size();
    int time_steps = sol.at(0).pos.cols();
    MatrixXd ak_mod = MatrixXd::Zero(time_steps,N);
    MatrixXd vk_mod = MatrixXd::Zero(time_steps,N);
    for (int i = 0; i<N; ++i)
    {
        ak_mod.col(i) = pow(((sol.at(i).acc.array().pow(2)).colwise().sum()),1.0/2.0);
        vk_mod.col(i) = pow(((sol.at(i).vel.array().pow(2)).colwise().sum()),1.0/2.0);
    }

    float r_factor = std::min(amax/ak_mod.maxCoeff(),vmax/vk_mod.maxCoeff());

    _h_scaled = _h/sqrt(r_factor);

    for(int i = 0; i < N; i++)
    {
        sol.at(i).acc = r_factor*sol.at(i).acc;
        for(int k=0; k< (time_steps-1); ++k)
        {
            sol.at(i).vel.col(k+1) = sol.at(i).vel.col(k) + _h_scaled*sol.at(i).acc.col(k);
        }
    }
}
std::vector<Trajectory> DMPC::interp_trajectory(const std::vector<Trajectory> &sol,
                                                const double &step_size)
{
    float T = (sol.at(0).pos.cols() - 1)*_h_scaled;
    int K = T/step_size + 1;
    int N = sol.size();
    double t0 = 0;
    std::vector<Trajectory> all_trajectories(N);
    VectorXd t = VectorXd::LinSpaced(K,0,(K-1)*step_size);
    std::vector<double> px;
    std::vector<double> py;
    std::vector<double> pz;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;
    std::vector<double> ax;
    std::vector<double> ay;
    std::vector<double> az;

    std::vector<double> px_inter;
    std::vector<double> py_inter;
    std::vector<double> pz_inter;
    std::vector<double> vx_inter;
    std::vector<double> vy_inter;
    std::vector<double> vz_inter;
    std::vector<double> ax_inter;
    std::vector<double> ay_inter;
    std::vector<double> az_inter;

    for (int n=0; n<N; ++n)
    {
        all_trajectories.at(n).pos = MatrixXd::Zero(3,K);
        all_trajectories.at(n).vel = MatrixXd::Zero(3,K);
        all_trajectories.at(n).acc = MatrixXd::Zero(3,K);
        px.clear();
        py.clear();
        pz.clear();
        vx.clear();
        vy.clear();
        vz.clear();
        ax.clear();
        ay.clear();
        az.clear();

        px_inter.clear();
        py_inter.clear();
        pz_inter.clear();
        vx_inter.clear();
        vy_inter.clear();
        vz_inter.clear();
        ax_inter.clear();
        ay_inter.clear();
        az_inter.clear();

        for(int i=0; i<sol.at(0).pos.row(0).cols(); ++i)
        {
            px.push_back(sol.at(n).pos.row(0).col(i).sum());
            py.push_back(sol.at(n).pos.row(1).col(i).sum());
            pz.push_back(sol.at(n).pos.row(2).col(i).sum());

            vx.push_back(sol.at(n).vel.row(0).col(i).sum());
            vy.push_back(sol.at(n).vel.row(1).col(i).sum());
            vz.push_back(sol.at(n).vel.row(2).col(i).sum());

            ax.push_back(sol.at(n).acc.row(0).col(i).sum());
            ay.push_back(sol.at(n).acc.row(1).col(i).sum());
            az.push_back(sol.at(n).acc.row(2).col(i).sum());
        }


        boost::math::cubic_b_spline<double> spline_px(px.begin(), px.end(),t0, _h_scaled);
        boost::math::cubic_b_spline<double> spline_py(py.begin(), py.end(),t0, _h_scaled);
        boost::math::cubic_b_spline<double> spline_pz(pz.begin(), pz.end(),t0, _h_scaled);

        boost::math::cubic_b_spline<double> spline_vx(vx.begin(), vx.end(),t0, _h_scaled);
        boost::math::cubic_b_spline<double> spline_vy(vy.begin(), vy.end(),t0, _h_scaled);
        boost::math::cubic_b_spline<double> spline_vz(vz.begin(), vz.end(),t0, _h_scaled);

        boost::math::cubic_b_spline<double> spline_ax(ax.begin(), ax.end(),t0, _h_scaled);
        boost::math::cubic_b_spline<double> spline_ay(ay.begin(), ay.end(),t0, _h_scaled);
        boost::math::cubic_b_spline<double> spline_az(az.begin(), az.end(),t0, _h_scaled);

        for(int i=0; i<t.size(); ++i)
        {
            px_inter.push_back(spline_px(t[i]));
            py_inter.push_back(spline_py(t[i]));
            pz_inter.push_back(spline_pz(t[i]));

            vx_inter.push_back(spline_vx(t[i]));
            vy_inter.push_back(spline_vy(t[i]));
            vz_inter.push_back(spline_vz(t[i]));

            ax_inter.push_back(spline_ax(t[i]));
            ay_inter.push_back(spline_ay(t[i]));
            az_inter.push_back(spline_az(t[i]));
        }

        all_trajectories.at(n).pos.row(0) = Map<VectorXd>(px_inter.data(),px_inter.size());
        all_trajectories.at(n).pos.row(1) = Map<VectorXd>(py_inter.data(),py_inter.size());
        all_trajectories.at(n).pos.row(2) = Map<VectorXd>(pz_inter.data(),pz_inter.size());

        all_trajectories.at(n).vel.row(0) = Map<VectorXd>(vx_inter.data(),vx_inter.size());
        all_trajectories.at(n).vel.row(1) = Map<VectorXd>(vy_inter.data(),vy_inter.size());
        all_trajectories.at(n).vel.row(2) = Map<VectorXd>(vz_inter.data(),vz_inter.size());

        all_trajectories.at(n).acc.row(0) = Map<VectorXd>(ax_inter.data(),ax_inter.size());
        all_trajectories.at(n).acc.row(1) = Map<VectorXd>(ay_inter.data(),ay_inter.size());
        all_trajectories.at(n).acc.row(2) = Map<VectorXd>(az_inter.data(),az_inter.size());
    }

    return all_trajectories;
}

bool DMPC::collision_violation(const std::vector<Trajectory> &solution)
{
    int N = solution.size();
    MatrixXd differ;
    VectorXd dist;
    bool violation = false;
    double min_dist;
    int pos;

    for (int i=0; i<N; ++i)
    {
        for(int j=i+1; j<N; ++j)
        {
            if (i!=j)
            {
                differ = _E1*(solution.at(i).pos-solution.at(j).pos);
                dist = pow(((differ.array().pow(_order)).colwise().sum()),1.0/_order);
                min_dist = dist.minCoeff(&pos);

                if (min_dist < _rmin-0.05) // we add 5cm tolerance factor to it
                {
                    violation = true;
//                    cout << "Collision constraint violation: ";
//                    cout << "Vehicles " << i << " and " << j;
//                    cout << " will be " << min_dist << "m";
//                    cout << " apart @ t = " << pos/100.0 << "s" << endl;
                }
            }
        }
    }

//    if (!violation)
//        cout << "No collisions found!" << endl;
    return violation;
}

void DMPC::trajectories2file(const std::vector<Trajectory> &solution,
                       char const* pathAndName)
{
    int N_cmd = solution.size();
    int N = _po.cols();
    ofstream file(pathAndName, ios::out | ios::trunc);
    if(file)  // succeeded at opening the file
    {
        // instructions
        cout << "Writing solution to text file..." << endl;

        // write a few simulation parameters used in the reading end
        file << N << " " <<  N_cmd << " "<< _h_scaled << " " << _pmin.transpose() << " " << _pmax.transpose() << endl;
        file << _po << endl;
        file << _pf << endl;

        for(int i=0; i < N_cmd; ++i)
        {
            file << solution.at(i).pos << endl;
        }

        for(int i=0; i < N_cmd; ++i)
        {
            file << solution.at(i).vel << endl;
        }

        for(int i=0; i < N_cmd; ++i)
        {
            file << solution.at(i).acc << endl;
        }

        file.close();  // close the file after finished
    }

    else
    {
        cerr << "Error while trying to open file" << endl;
    }
}

/***************************************************************************
    *
    *  CPLEX functions
    *
    **************************************************************************/

void DMPC::init_cplex(int id) {

    _env.at(id) = NULL;
    _lp.at(id) = NULL;
    int status;

    // Initialize CPLEX
    _env.at(id) = CPXopenCPLEX(&status);
    if (_env.at(id) == NULL) {
        char errmsg[1024];
        fprintf(stderr, "Could not open CPLEX .\n");
        CPXgeterrorstring(_env.at(id), status, errmsg);
        fprintf(stderr, "%s", errmsg);
        terminate_cplex(id);
    }
    status = CPXsetintparam(_env.at(id), CPXPARAM_ScreenOutput, CPX_OFF);
    if ( status ) {
        fprintf(stderr, "Failure to turn on screen indicator, error %d.\n", status);
        terminate_cplex(id);
    }
    status = CPXsetintparam(_env.at(id), CPXPARAM_Read_DataCheck, CPX_OFF);
    if ( status ) {
        fprintf(stderr,
                "Failure to turn on data checking, error %d.\n", status);
        terminate_cplex(id);
    }
    // Create problem
    _lp.at(id) = CPXcreateprob(_env.at(id), &status, "This_problem");
    if (_lp.at(id) == NULL) {
        printf("Failed to create problem. \n");
        status = 1;
        terminate_cplex(id);
    }
    return;
}

void DMPC::terminate_cplex(int id) {

    int status = 0;
    std::cout << "Terminating CPLEX" << std::endl;
    // Free lplem allocated by CPXcreatelp
    if (_lp.at(id) != NULL) {
        status = CPXfreeprob(_env.at(id), &_lp.at(id));
        if (status)
            fprintf(stderr, "CPXfreelp failed, error code %d. \n", status);
    }
    // Free cplex enviromnent
    if (_env.at(id) != NULL) {
        status = CPXcloseCPLEX(&_env.at(id));
        if (status) {
            char errmsg[1024];
            fprintf(stderr, "Could not close CPLEX enviroment.\n");
            CPXgeterrorstring(_env.at(id), status, errmsg);
            fprintf(stderr, "%s", errmsg);
        }
    }
}

void DMPC::eigen_to_cplex(const Eigen::MatrixXd &H, int *&matbeg,
                          int *&matcnt, int *&matind, double *&matval,
                          int &numrows, int &numcols, int &numnz) {
    // Process matrices to fit cplex matrix requirements
    int i, j, cnt, k;
    numrows = H.rows();
    numcols = H.cols();
    numnz = 0;

    // Find nonzero values
    for (i = 0; i < numrows; i++) {
        for (j = 0; j < numcols; j++) {
            if (H(i, j) != 0.0)
                numnz++;
        }
    }

    matbeg = new int[numcols];  // Contains index of beginning of column j
    matcnt = new int[numcols];  // Contains how many nz entries in columns j
    matind = new int[numnz];  // Contains row number of coefficient value val[k]
    matval = new double[numnz];
    numnz = 0;
    k = 0;
    // Arrange nonzero values into corresponding vectors
    for (j = 0; j < numcols; j++) {
        cnt = 0;
        for (i = 0; i < numrows; i++) {
            if (H(i, j) != 0.0) {
                cnt++;  // Number of nz values in this column
                numnz++;
                matind[k] = i;
                matval[k] = H(i, j);
                k++;
            }
        }
        matbeg[j] = numnz-cnt;
        matcnt[j] = cnt;
    }
}