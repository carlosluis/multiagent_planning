//
// Created by carlos on 16/03/18.
//

#include <iostream>
#include "dmpc.h"

using namespace Eigen;
using namespace std;

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
    _Lambda = this->get_lambda_mat(_h,_k_hor);
    _Delta = this->get_delta_mat(_k_hor);
    _A0 = this->get_A0_mat(_k_hor);


    // Just testing some functions
    _pmin << -2.5, -2.5, 0.2;
    _pmax << 2.5, 2.5, 2.2;
    Vector3d po1(0,0,1.5);
    Vector3d po2(0,1,1.5);
    Vector3d pf1(0,1,1.5);
    Vector3d pf2(0,0,1.5);
    Trajectory agent1 = this->init_dmpc(po1,pf1);
    Trajectory agent2 = this->init_dmpc(po2,pf2);
    cout << "Init Trajectory agent 1:" << endl << agent1.pos << endl;
    cout << "Init Trajectory agent 2:" << endl << agent2.pos << endl;

    // build obstacle list to test different functions
    std::vector<MatrixXd> obs;
    obs.push_back(agent1.pos);
    obs.push_back(agent2.pos);
    bool violation = check_collisions(agent1.pos.col(13),obs,0,13);
    cout << "violation = " << violation << endl;
    Vector3d vo(0,0,0);
    Vector3d ao(0,0,0);
    Constraint collisions = build_collconstraint(agent1.pos.col(13),po1,vo,obs,0,13);
    cout << "Collision matrix A for agent 1" << endl << collisions.A << endl;
    cout << "Collision matrix b for agent 1" << endl << collisions.b << endl;
    Trajectory sol_agent1 = solveDMPC(po1,pf1,vo,ao,0,obs);

}

MatrixXd DMPC::get_lambda_mat(int h, int K)
{
    MatrixXd Apos = MatrixXd::Zero(3*K,3*K);
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
        prev_row = new_row;
        idx += 3;
    }
    return Apos;
}

MatrixXd DMPC::get_delta_mat(int K)
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
    return Delta;
}

MatrixXd DMPC::get_A0_mat(int K)
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
    return A0;
}

Trajectory DMPC::init_dmpc(Vector3d po, Vector3d pf)
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

MatrixXd DMPC::gen_rand_pts(int N, Vector3d pmin, Vector3d pmax, float rmin)
{
    MatrixXd pts = MatrixXd::Zero(3,N);
    Vector3d candidate = MatrixXd::Zero(3,1);
    VectorXd dist;

    // Generate first point
    std::srand((unsigned int) time(0));
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
    cout << "Random Points are:" << endl << pts << endl;
    return pts;
}

void DMPC::set_initial_pts(MatrixXd po)
{
    _po = po;
}

void DMPC::set_final_pts(MatrixXd pf)
{
    _pf = pf;
}

bool DMPC::check_collisions(Vector3d prev_p, std::vector<MatrixXd> obs,
                            int n, int k)
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

Constraint DMPC::build_collconstraint(Vector3d prev_p, Vector3d po,
                                      Vector3d vo, std::vector<MatrixXd> obs,
                                      int n, int k)
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

Trajectory DMPC::solveDMPC(Vector3d po,Vector3d pf,
                                  Vector3d vo,Vector3d ao,
                                  int n, std::vector<MatrixXd> obs)
{
    int N = obs.size(); // number of vehicles for transition
    bool violation;
    Trajectory solution;
    Constraint coll_constraint;
    Constraint total_collconstr;
    MatrixXd prev_p = obs.at(n); // previous solution of n-th vehicle
    MatrixXd collconstrA_aug(2*(N-1), 3*_k_hor + N -1);
    VectorXd collconstrb_aug(2*(N-1));
    VectorXd a0_1;

    // Initial state joint vector
    Matrix <double, 6,1> initial_states;
    initial_states << po,vo;

    // Define optimization problem constants

    // Amount of variables/ ineq constraints: no collision case
    const int n_var = 3*_k_hor; // 3D acceleration vector for the horizon length
    const int n_ineq = 2*(3*_k_hor)+2*(3*_k_hor); // workspace boundaries + acc limits

    // Amount of variables/ ineq constraints: collision case
    const int n_var_aug = n_var + N -1; // add slack variable for constraint relaxation
    const int n_ineq_aug = 2*(N-1) + n_ineq; // add collisions + slack variable < 0

    // Number of variables / ineq constraints to be used by the qp
    int qp_nvar;
    int qp_nineq;
    const int qp_neq = 0; // No equality constraints in this problem

    // Penalty matrices for DMPC calculation
    MatrixXd Q(n_var,n_var);
    MatrixXd Q_aug(n_var_aug,n_var_aug);
    MatrixXd R(n_var,n_var);
    MatrixXd R_aug(n_var_aug,n_var_aug);
    MatrixXd S(n_var,n_var);
    MatrixXd S_aug(n_var_aug,n_var_aug);

    // Penalty matrices added in case of collision violation
    VectorXd f_w(n_var_aug);
    MatrixXd W(n_var_aug,n_var_aug);

    // Final Quadratic and linear matrices for QP solver
    MatrixXd H;
    VectorXd f;

    // Augmented model matrices in case of collisions
    MatrixXd Lambda_aug(n_var_aug,n_var_aug);
    MatrixXd Lambda_aug_in(n_var,n_var_aug);
    MatrixXd Delta_aug(n_var_aug,n_var_aug);
    MatrixXd A0_aug(n_var_aug,6);

    // Auxiliary variables
    VectorXd init_propagation = _A0*initial_states;
    VectorXd init_propagation_aug;
    VectorXd pf_rep(n_var_aug);
    pf_rep << pf.replicate(_k_hor,1),MatrixXd::Zero(N-1,1);
    VectorXd alim_rep = _alim*MatrixXd::Ones(3*_k_hor,1);

    // Constraint matrices and vectors to pass to the QP solver
    MatrixXd Ain;
    MatrixXd bin;

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

        bin = MatrixXd::Zero(n_ineq_aug,1);
        bin << collconstrb_aug,
               _pmax.replicate(_k_hor,1) - init_propagation,
               -_pmin.replicate(_k_hor,1) + init_propagation,
               alim_rep, alim_rep; // input limit constraint

        f_w << MatrixXd::Zero(n_var,1),
               -pow(10,5)*MatrixXd::Ones(N-1,1);

        W.block(n_var,n_var,N-1,N-1) = pow(10,-9)*(MatrixXd::Identity(N-1,N-1));

        a0_1 = MatrixXd::Zero(n_var_aug,1);
        a0_1 << ao, MatrixXd::Zero(3*(_k_hor-1) + N - 1,1);
        f = MatrixXd::Zero(n_var_aug,1);

        f = -2*(pf_rep.transpose()*Q_aug*Lambda_aug -
                init_propagation_aug.transpose()*Q_aug*Lambda_aug +
                a0_1.transpose()*S_aug*Delta_aug);

        f += f_w;

        H = 2*(Lambda_aug.transpose()*Q_aug*Lambda_aug
               + Delta_aug.transpose()*S_aug*Delta_aug
               + R_aug + W);

        Ain = MatrixXd::Zero(n_ineq_aug,n_var_aug);

        Ain << collconstrA_aug,
               Lambda_aug_in, -Lambda_aug_in,
               MatrixXd::Identity(n_var,n_var_aug),
               -MatrixXd::Identity(n_var,n_var_aug);
    }

    else // matrices when there're no collisions
    {

    }

//    cout << H << endl;

    _qp = new QuadProgDense(qp_nvar,qp_neq,qp_nineq);
    _qp->solve(H,f,MatrixXd::Zero(0, qp_nvar),VectorXd::Zero(0),Ain,bin);
    cout << "Exit flag = " << _qp->fail() << endl;
    cout << "Result is = " << _qp->result() << endl;

    return solution;
}