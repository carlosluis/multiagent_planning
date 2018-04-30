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

struct State {
    MatrixXd pos;
    MatrixXd vel;
    MatrixXd acc;
};

bool CheckforColl(MatrixXd p, std::vector<MatrixXd> L, int k, int r_min);

MatrixXd getPosMat(int h, int K);

MatrixXd getPosVelMat(int h, int K);

MatrixXd initSolution(MatrixXd po, MatrixXd pf, int h, int K);

Constraint collConstr(MatrixXd p, MatrixXd po, int k, std::vector<MatrixXd> L, MatrixXd Ain, int r_min);

State propState(MatrixXd po, MatrixXd a, int h);

int maxDeviation(MatrixXd p, MatrixXd);







#endif //DMPC_CPP_DMPC_H
