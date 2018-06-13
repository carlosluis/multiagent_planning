//
// Created by carlos on 16/03/18.
//

#include <iostream>
#include "dmpc.h"

using namespace Eigen;
using namespace std;
using namespace std::chrono;

bool execution_ended;
int failed_i_global;

DMPC::DMPC(Params params)
{
    // Load parameters into private variables
    _h = params.h;
    _T = params.T;
    _k_hor = params.k_hor;
    _order = params.order;
    _c = params.c;
    _rmin = params.rmin;
    _alim = params.alim;
    _K =  _T/_h + 1; // number of time steps
    srand((unsigned int) time(0)); // initialize srand DO THIS ONLY ONCE

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
    _pmin << -4, -4, 0.2;
    _pmax << 4, 4, 3.2;
}

void DMPC::get_lambda_A_v_mat(const int &K)
{
    MatrixXd Apos = MatrixXd::Zero(3*K,3*K);
    MatrixXd Avel = MatrixXd::Zero(3*K,3*K);
    MatrixXd prev_row =  MatrixXd::Zero(6,3*K);
    MatrixXd new_row =  MatrixXd::Zero(6,3*K);
    MatrixXd add_b = MatrixXd::Zero(6,3*K);
    int idx = 0;
    for(int k=0; k < K; ++k)
    {
        add_b << MatrixXd::Zero(_b.rows(),_b.cols()*(k)),_b,
                MatrixXd::Zero(_b.rows(),_b.cols()*(K-k-1));
        new_row = _A*prev_row + add_b;
        Apos.middleRows(idx,3) = new_row.middleRows(0,3);
        Avel.middleRows(idx,3) = new_row.middleRows(3,3);
        prev_row = new_row;
        idx += 3;
    }
    _Lambda = Apos;
    _A_v = Avel;
}

void DMPC::get_delta_mat(const int &K)
{
    MatrixXd Delta = MatrixXd::Zero(3*K,3*K);
    MatrixXd new_row = MatrixXd::Zero(3,3*K);
    Delta.topRows(3) << MatrixXd::Identity(3,3),
                        MatrixXd::Zero(3,3*(K-1));
    MatrixXd b = MatrixXd::Zero(3,6);
    b << (-1)*MatrixXd::Identity(3,3), MatrixXd::Identity(3,3);
    int idx = 3;
    for(int k=0; k<K-1; ++k)
    {
       new_row <<  MatrixXd::Zero(3,3*k),b,MatrixXd::Zero(3,3*(K-k-2));
       Delta.middleRows(idx,3) = new_row;
       idx += 3;
    }

    _Delta = Delta;
}

void DMPC::get_A0_mat(const int &K)
{
    MatrixXd A0 = MatrixXd::Zero(3*K,6);
    MatrixXd new_row = MatrixXd::Zero(6,6);
    MatrixXd prev_row = MatrixXd::Identity(6,6);
    int idx = 0;
    for(int k=0; k<K; ++k)
    {
        new_row = _A*prev_row;
        A0.middleRows(idx,3) = new_row.middleRows(0,3);
        prev_row = new_row;
        idx += 3;
    }
    _A0 = A0;
}

Trajectory DMPC::init_dmpc(const Vector3d &po, const Vector3d &pf)
{
    Trajectory init;
    Vector3d diff = pf - po;
    VectorXd t = VectorXd::LinSpaced(_K,0,(_K-1)*_h);
    VectorXd t_hor = t.head(_k_hor);
    init.pos = MatrixXd::Zero(3,t_hor.size());
    init.vel = MatrixXd::Zero(3,t_hor.size());
    init.acc = MatrixXd::Zero(3,t_hor.size());
    for (int i=0; i<t_hor.size(); ++i)
    {
        init.pos.col(i) = po + t_hor[i]*diff/((_K-1)*_h);
    }
    return init;
}

MatrixXd DMPC::gen_rand_pts(const int &N,
                            const Vector3d &pmin,
                            const Vector3d &pmax,
                            const float &rmin)
{
    MatrixXd pts = MatrixXd::Zero(3,N);
    Vector3d candidate = MatrixXd::Zero(3,1);
    VectorXd dist;

    // Generate first point
    pts.col(0) = pmin.array()
                 + (pmax-pmin).array()*((MatrixXd::Random(3,1).array() + 1)/2);
    bool pass = false;
    for (int n=1; n < N; ++n)
    {
        while(!pass)
        {
            candidate = pmin.array()
                        + (pmax-pmin).array()*((MatrixXd::Random(3,1).array() + 1)/2);
            dist = ((((pts.leftCols(n)).colwise()
                   - candidate).array().square()).colwise().sum()).array().sqrt();

            for (int k=0; k < n; ++k){
                pass = dist[k] > rmin;
                if (!pass)
                    break;
            }
            if (pass)
                pts.col(n) = candidate.array();
        }
        pass = false;
    }
//    cout << "Random Points are:" << endl << pts << endl;
    return pts;
}

MatrixXd DMPC::gen_rand_perm(const MatrixXd &po)
{
    int N = po.cols();
    std::vector<int> array;
    std::vector<int> array_aux;
    int perm[N];
    MatrixXd pf = MatrixXd::Zero(3,N);
    for (int i = 0; i < N; i++) array.push_back(i);
    array_aux = array;
    cout << "perm = " << perm << endl;

    // Random permutation the order
    for (int i = 0; i < N; i++) {
        int j;
        array_aux.clear();
        array_aux = array;
        array_aux.erase(std::remove(array_aux.begin(), array_aux.end(), i), array_aux.end());
        if (i == N-1){
            perm[i] = array.at(0);
        }

        else if (i == N-2 && array_aux.back() == N-1){
            perm[i] = array_aux.back();
            array.erase(std::remove(array.begin(), array.end(), perm[i]), array.end());
        }
        else{
            j = rand() % (N-i-1);
            perm[i] = array_aux.at(j);
            array.erase(std::remove(array.begin(), array.end(), array_aux.at(j)), array.end());
        }
    }

    for (int i = 0; i < N; i++) {
        pf.col(i) = po.col(perm[i]);
    }
    return pf;
}

void DMPC::set_initial_pts(const MatrixXd &po)
{
    _po = po;
}

void DMPC::set_final_pts(const MatrixXd &pf)
{
    _pf = pf;
}

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

    for (int i=0; i < N_obs; ++i)
    {
        if(i!=n)
        {
            pj = obs[i].col(k);
            diff = _E1*(prev_p - pj);
            dist = pow(((diff.array().pow(_order)).sum()),1.0/_order);
            diff = (_E2*(prev_p - pj)).array().pow(_order - 1);

            r = pow(dist,_order-1)*(_rmin - dist) + diff.transpose()*prev_p
                - diff.transpose()*_A0.middleRows(3*k,3)*initial_states;

            diff_row << MatrixXd::Zero(1,3*k),
                        diff.transpose(),
                        MatrixXd::Zero(1,3*(_k_hor-k-1));
            Ain_total.row(idx) = -diff_row*_Lambda;
            bin_total[idx] = -r;
            idx++;
        }
    }
    collision.A = Ain_total;
    collision.b = bin_total;
    return collision;
}

Trajectory DMPC::solveQP(const Vector3d &po, const Vector3d &pf,
                   const Vector3d &vo, const Vector3d &ao,
                   const int &n, const std::vector<MatrixXd> &obs)
{
    int N = obs.size(); // number of vehicles for transition
    bool violation; // check if collision constraint is violated in horizon
    Trajectory solution;
    Constraint coll_constraint;
    MatrixXd prev_p = obs.at(n); // previous solution of n-th vehicle
    VectorXd a0_1;

    // Initial state joint vector (pos + vel)
    Matrix <double, 6,1> initial_states;
    initial_states << po,vo;

    // Define optimization problem constants

    // Amount of variables/ ineq constraints: no collision case
    int n_var = 3*_k_hor; // 3D acceleration vector for the horizon length
    int n_ineq = 2*(3*_k_hor)+2*(3*_k_hor); // workspace boundaries + acc limits

    // Amount of variables/ ineq constraints: collision case
    int n_var_aug = n_var + N -1; // add slack variable for constraint relaxation
    int n_ineq_aug = 2*(N-1) + n_ineq; // add collisions + slack variable < 0

    // Number of variables / ineq constraints to be used by the qp
    int qp_nvar;
    int qp_nineq;
    int qp_neq = 0; // No equality constraints in this problem

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
        Q.block(3*(_k_hor-1),3*(_k_hor-1),3,3) = 1000*MatrixXd::Identity(3,3);
        R = 1*MatrixXd::Identity(n_var,n_var);
        S = 10*MatrixXd::Identity(n_var,n_var);
    }

    // Case of no collisions and close to goal
    else if (!violation && (po-pf).norm() < 1)
    {
        Q.block(3*(_k_hor-1),3*(_k_hor-1),3,3) = 10000*MatrixXd::Identity(3,3);
        R = 1*MatrixXd::Identity(n_var,n_var);
        S = 10*MatrixXd::Identity(n_var,n_var);
    }

    // Case of collisions in the horizon
    else
    {
        Q.block(3*(_k_hor-1),3*(_k_hor-1),3,3) = 1000*MatrixXd::Identity(3,3);
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
               -pow(10,5)*VectorXd::Ones(N-1);

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

//    t2 = high_resolution_clock::now();
//    duration = duration_cast<microseconds>( t2 - t1 ).count();
//    cout << "Building the QP computation time = "
//         << duration/1000000.0 << "s" << endl;

//    t1 = high_resolution_clock::now();
    // Declare quadprog object and solve the QP
    QuadProgDense _qp(qp_nvar,qp_neq,qp_nineq);

    _qp.solve(H,f,MatrixXd::Zero(0, qp_nvar),VectorXd::Zero(0),Ain,bin);

//    t2 = high_resolution_clock::now();
//    duration = duration_cast<microseconds>( t2 - t1 ).count();
//    cout << "Solving the optimization computation time = "
//         << duration/1000000.0 << "s" << endl;

//    t1 = high_resolution_clock::now();
    x = _qp.result();
    _fail = _qp.fail();
    if(_fail){
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

//    t2 = high_resolution_clock::now();
//    duration = duration_cast<microseconds>( t2 - t1 ).count();
//    cout << "Formatting the solution computation time = "
//         << duration/1000000.0 << "s" << endl
//         << "-----------------------------" << endl;

    return solution;
}

std::vector<Trajectory> DMPC::solveDMPC(const MatrixXd &po,
                                        const MatrixXd &pf)
{
    int N = po.cols(); // Number of agents = number of rows of po
    int N_cmd = pf.cols(); // Number of agents of transition = number of rows of pf

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
            rep = (po.col(i+N_cmd)).replicate(_k_hor,1);
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
                    poi = po.col(i);
                    pfi = pf.col(i);
                    trajectory_i = init_dmpc(poi,pfi);
                    _fail = 0;
                }
                else
                {   // Update previous solution, solve current QP
//                t1 = high_resolution_clock::now();
                    pok = all_trajectories.at(i).pos.col(k-1);
                    vok = all_trajectories.at(i).vel.col(k-1);
                    aok = all_trajectories.at(i).acc.col(k-1);
                    trajectory_i = solveQP(pok,pf.col(i),vok,aok,i,prev_obs);
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
    cout << "Solve al QPs computation time = " << duration/1000000.0 << "s" << endl;

    if(!_fail)
    {
        cout << "Optimization problem feasible: solution found" << endl;
        // Check if every agent reached its goal
        arrived = reached_goal(all_trajectories,pf,0.05,N_cmd);
    }

    t1 = high_resolution_clock::now();

    if (arrived && !_fail)
    {
        cout << "All vehicles reached their goals" << endl;
        // Interpolate for better resolution (e.g. 100 Hz)
        solution = interp_trajectory(all_trajectories,0.01);

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

std::vector<Trajectory> DMPC::solveParallelDMPC(const MatrixXd &po,
                                          const MatrixXd &pf)
{
    int N = po.cols(); // Number of agents = number of rows of po
    int N_cmd = pf.cols(); // Number of agents of transition = number of rows of pf

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
    int residue = N%n_cluster;
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
            rep = (po.col(i)).replicate(_k_hor,1);
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
    cout << "Solve al QPs computation time = " << duration/1000000.0 << "s" << endl;

    if(!execution_ended)
    {
        cout << "Optimization problem feasible: solution found" << endl;
        // Check if every agent reached its goal
        arrived = reached_goal(all_trajectories,pf,0.05,N_cmd);
    }

    t1 = high_resolution_clock::now();

    if (arrived && !execution_ended)
    {
        cout << "All vehicles reached their goals" << endl;
        // Interpolate for better resolution (e.g. 100 Hz)
        solution = interp_trajectory(all_trajectories,0.01);

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
//            failed_i = i;
            break;
        }

        // If QP didn't fail, then update the solution vector and obstacle list
        obs.at(i) = trajectory_i.pos;
        all_trajectories.at(i).pos.col(k) = trajectory_i.pos.col(0);
        all_trajectories.at(i).vel.col(k) = trajectory_i.vel.col(0);
        all_trajectories.at(i).acc.col(k) = trajectory_i.acc.col(0);
    }
}

bool DMPC::reached_goal(const std::vector<Trajectory> &all_trajectories,
                        const MatrixXd &pf, const float &error_tol, const int &N)
{   Vector3d diff;
    double dist;
    bool reached = true;
    for (int i=0; i < N; ++i)
    {
        diff = all_trajectories.at(i).pos.col(_K-1) - pf.col(i);
        dist = pow(((diff.array().pow(2)).sum()),1.0/2);
        if (dist > error_tol){
            cout << "Vehicle " << i << " did not reach its goal by " << dist << "m" << endl;
            reached = false;
        }
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
                time[i] = k/100.0;
                break;
            }
        }
    }
    min_traj_time = time.maxCoeff();
    cout << "The transition will be completed in " << min_traj_time << "s" << endl;
    return min_traj_time;
}

std::vector<Trajectory> DMPC::interp_trajectory(const std::vector<Trajectory> &sol,
                                                const double &step_size)
{
    int K = _T/step_size + 1;
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


        boost::math::cubic_b_spline<double> spline_px(px.begin(), px.end(),t0, _h);
        boost::math::cubic_b_spline<double> spline_py(py.begin(), py.end(),t0, _h);
        boost::math::cubic_b_spline<double> spline_pz(pz.begin(), pz.end(),t0, _h);

        boost::math::cubic_b_spline<double> spline_vx(vx.begin(), vx.end(),t0, _h);
        boost::math::cubic_b_spline<double> spline_vy(vy.begin(), vy.end(),t0, _h);
        boost::math::cubic_b_spline<double> spline_vz(vz.begin(), vz.end(),t0, _h);

        boost::math::cubic_b_spline<double> spline_ax(ax.begin(), ax.end(),t0, _h);
        boost::math::cubic_b_spline<double> spline_ay(ay.begin(), ay.end(),t0, _h);
        boost::math::cubic_b_spline<double> spline_az(az.begin(), az.end(),t0, _h);

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
                    cout << "Collision constraint violation: ";
                    cout << "Vehicles " << i << " and " << j;
                    cout << " will be " << min_dist << "m";
                    cout << " apart @ t = " << pos/100.0 << "s" << endl;
                }
            }
        }
    }
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
        file << N << " " <<  N_cmd << " " << _pmin.transpose() << " " << _pmax.transpose() << endl;
        file << _po << endl;
        file << _pf << endl;

        for(int i=0; i < N_cmd; ++i)
        {
            file << solution.at(i).pos << endl;

        }

        file.close();  // close the file after finished
    }

    else
    {
        cerr << "Error while trying to open file" << endl;
    }
}