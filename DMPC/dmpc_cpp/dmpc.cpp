//
// Created by carlos on 16/03/18.
//

#include <iostream>
#include "dmpc.h"

using namespace Eigen;
using namespace std;

SoftDMPC::SoftDMPC():
        _h(0.2),
        _T(20),
        _K(0),
        _k_hor(15),
        _order(2),
        _c(1.5),
        _rmin(0.5),
        _alim(0.5)
{
    _K =  _T/_h + 1; // number of time steps

    // Ellipsoid constants
    Vector3f v(3);
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
    float h2 = pow(_h,2);
    _b << h2/2,0,0,
          0,h2/2,0,
          0,0,h2/2,
          0,0,_h,
          0,_h,0,
          0,0,_h;
}

MatrixXd get_lambda_mat(int h, int K)
{
    // Kinematic model A,B matrices


}


