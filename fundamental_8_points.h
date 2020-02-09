#ifndef _FUNDAMENTAL_8_POINTS_
#define _FUNDAMENTAL_8_POINTS_

#include "math\matrix_svd.h"
#include "math\matrix.h"
#include "math\vector.h"

math::Matrix<double,3,3> fundamental_8_points(const math::Matrix<double,3,8>& points1
                                            ,const math::Matrix<double,3,8>& points2)
    {
        //direct linear transform
        math::Matrix<double,8,9> A;
        for(int i = 0; i < 8; i++)
        {
            math::Vec3d p1 = points1.col(i);
            math::Vec3d p2 = points2.col(i);

            A(i,0) = p1[0]*p2[0];
            A(i,1) = p1[0]*p2[1];
            A(i,2) = p1[0];
            A(i,3) = p1[1]*p2[0];
            A(i,4) = p1[1]*p2[1];
            A(i,5) = p1[1];
            A(i,6) = p2[0];
            A(i,7) = p2[1];
            A(i,8) = 1.0;

        }

        math::Matrix<double,9,9> vv;
        math::matrix_svd<double,8,9>(A,nullptr,nullptr,&vv);
        math::Vector<double,9> f = vv.col(8);

        //奇异值约束
        math::Matrix<double,3,3> F;
        F(0,0) = f[0];  F(0,1) = f[1];  F(0,2) = f[2];
        F(1,0) = f[3];  F(1,1) = f[4];  F(1,2) = f[5];
        F(2,0) = f[6];  F(2,1) = f[7];  F(2,2) = f[8];

        math::Matrix<double,3,3> U,S,V;
        math::matrix_svd(F,&U,&S,&V);
        S(2,2) = 0;
        F = (U*S*V.transposed()).transposed();

        return F;


    }

#endif