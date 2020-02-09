#ifndef _FUNDAMENTAL_RANSAC_
#define _FUNDAMENTAL_RANSAC_

#include <vector>
#include "math\functions.h"
#include "math\matrix_svd.h"
#include "fundamental_8_points.h"
#include <set>

typedef math::Matrix<double,3,3> FundamentalMatrix;

/*
存储匹配点
*/
struct correspondance2D2D
{
    double p1[2];
    double p2[2];
};
typedef std::vector<correspondance2D2D> correspondances2D2D;

/**
 * description 用于RANSAC采样成功所需要的采样次数
 * @param p -- 内点的概率
 * @param K --拟合模型需要的样本个数，对应基础矩阵num_samples=8
 * @param z  -- 预期的采样成功的概率
 *                          log(1-z)
 *       需要的采样次数 M = -----------
 *                          log(1-p^K)
 * Example: For p = 50%, z = 99%, n = 8: M = log(0.001) / log(0.99609) = 1176.
 * 需要采样1176次从而保证RANSAC的成功率不低于0.99.
 * @return

*/
int calc_ransac_iterations(double p,
                            int k,
                            double z = 0.99)
{
    double prob_all_good = math::fastpow(p,k);
    double num_iterations = std::log(1-z)
                            / std::log(1-prob_all_good);
    return static_cast<int>(math::round(num_iterations));
}

/**
 * \description 给定基础矩阵和一对匹配点，计算匹配点的sampson 距离，用于判断匹配点是否是内点,
 * 计算公式如下：
 *              SD = (x'Fx)^2 / ( (Fx)_1^2 + (Fx)_2^2 + (x'F)_1^2 + (x'F)_2^2 )
 * @param F-- 基础矩阵
 * @param m-- 匹配对
 * @return
 */
double calc_sampon_distance(const FundamentalMatrix& F,const correspondance2D2D& m)
{
    // double p2_F_p1 = 0;
    // p2_F_p1 += m.p2[0]*(m.p1[0]*F(0,0)+m.p1[1]*F(1,0)+m.p1[2]*F(2,0));
    // p2_F_p1 += m.p2[1]*(m.p1[0]*F(0,1)+m.p1[1]*F(1,1)+m.p1[2]*F(2,1));
    // p2_F_p1 += m.p2[2]*(m.p1[0]*F(0,2)+m.p1[1]*F(1,2)+m.p1[2]*F(2,2));
    // p2_F_p1 *= p2_F_p1;

    // double sum = 0.0;
    // sum += math::fastpow(F(0,0)*m.p1[0]+F(0,1)*m.p1[1]+F(0,2)*m.p1[2],2);
    // sum += math::fastpow(F(1,0)*m.p1[0]+F(1,1)*m.p1[1]+F(1,2)*m.p1[2],2);
    // sum += math::fastpow(m.p2[0]*F(0,0)+m.p2[1]*F(1,0)+m.p2[2]*F(2,0),2);
    // sum += math::fastpow(m.p2[0]*F(0,1)+m.p2[1]*F(1,1)+m.p2[2]*F(2,1),2);

    double p2_F_p1 = 0.0;
    p2_F_p1 += m.p2[0] * (m.p1[0] * F[0] + m.p1[1] * F[1] + F[2]);
    p2_F_p1 += m.p2[1] * (m.p1[0] * F[3] + m.p1[1] * F[4] + F[5]);
    p2_F_p1 +=     1.0 * (m.p1[0] * F[6] + m.p1[1] * F[7] + F[8]);
    p2_F_p1 *= p2_F_p1;

    double sum = 0.0;
    sum += math::fastpow(m.p1[0] * F[0] + m.p1[1] * F[1] + F[2], 2);
    sum += math::fastpow(m.p1[0] * F[3] + m.p1[1] * F[4] + F[5], 2);
    sum += math::fastpow(m.p2[0] * F[0] + m.p2[1] * F[3] + F[6], 2);
    sum += math::fastpow(m.p2[0] * F[1] + m.p2[1] * F[4] + F[7], 2);

    return p2_F_p1 / sum;
}

/*
* description 利用最小二乘法计算基础矩阵
* @param matches --输入的匹配对，必须大于8对
* @return F --基础矩阵
*/
FundamentalMatrix calc_fundamental_least_squares(const correspondances2D2D& matches)
{
    if(matches.size() < 8)
    {
        throw std::invalid_argument("At least 8 points required");
    }

    //创建 Nx9 维矩阵A，每一对匹配点生成矩阵A的一行。
    std::vector<double> A(matches.size()*9);
    for (int i = 0; i < matches.size(); i++)
    {
        const correspondance2D2D& p = matches[i];
            A[i*9+0] = p.p1[0]*p.p2[0];
            A[i*9+1] = p.p1[0]*p.p2[1];
            A[i*9+2] = p.p1[0];
            A[i*9+3] = p.p1[1]*p.p2[0];
            A[i*9+4] = p.p1[1]*p.p2[1];
            A[i*9+5] = p.p1[1];
            A[i*9+6] = p.p2[0];
            A[i*9+7] = p.p2[1];
            A[i*9+8] = 1.0;
    }

    //利用SVD计算基础矩阵
    std::vector<double> vv(9*9);
    math::matrix_svd<double>(&A[0],matches.size(),9,nullptr,nullptr,&vv[0]);

    FundamentalMatrix F;

    //将vv的最后一列作为解
    for(int i = 0; i<9; i++)
    {
        F[i] = vv[i*9+8];
    }

    //极限约束
    math::Matrix<double,3,3>U,S,V;
    math::matrix_svd(F,&U,&S,&V);
    S(2,2) = 0;
    F = -(U*S*V.transposed()).transposed();

    return F;
    
}

/*
* descripation 给定匹配对和基础矩阵，计算内点个数
* param matches
* param F
* return
*/
std::vector<int> find_inliers(const correspondances2D2D& matches,const FundamentalMatrix& F
                                ,const double& thresh)
{
    const double squared_thresh = thresh*thresh;
    std::vector<int> inliers;
    for(int i = 0; i < matches.size(); i++)
    {
        double error = calc_sampon_distance(F,matches[i]);
        if(error < squared_thresh)
        {
            inliers.push_back(i);
        }
    }
    return inliers;
}

/*
* descripation 利用ransac方法计算基础矩阵
* @matches -匹配点的集合
* @inlier_ratio -内点的概率
* @n_samples -求解模型需要的最少的点数
* @prob_acc -最终想要达到的准确率
* @inlier_thresh 属于内点的阈值
*/
FundamentalMatrix calc_fundamental_ransac(const correspondances2D2D& matches
                                        ,const float inlier_ratio = 0.5,
                                        const int n_samples = 8,
                                        const double prob_acc = 0.99,
                                        const double inlier_thresh = 0.0015)
{
    //计算采用次数
    int n_iterations = calc_ransac_iterations(inlier_ratio,n_samples);

    //ransac最终估计的内点
    std::vector<int> best_inliers;

    std::cout<<"RANSAC-F: Running for "<<n_iterations
             <<" iterations,threshold: "<<inlier_thresh
             <<"..."<<std::endl;
    for(int i = 0; i<n_iterations; i++)
    {
        // 1.随机找到8对不重复的点
        std::vector<int> indices;
        while(indices.size()<8)
        {
            
            int it = std::rand() % matches.size();
            
            indices.push_back(it);
        }

        math::Matrix<double,3,8> pset1,pset2;
        std::vector<int>::const_iterator iter = indices.cbegin();

        for(int j = 0; j <8; j++, iter++)
        {
            const correspondance2D2D& match = matches[*iter];
            pset1(0,j) = match.p1[0];
            pset1(1,j) = match.p1[1];
            pset1(2,j) = 1;

            pset2(0,j) = match.p2[0];
            pset2(1,j) = match.p2[1];
            pset2(2,j) = 1;
        } 

        //2. 8点法估计相机基础矩阵
        FundamentalMatrix F = fundamental_8_points(pset1,pset2);

        //3. 统计所有的内点个数
        std::vector<int> inlier_indices = find_inliers(matches,F,inlier_thresh);

        if(inlier_indices.size()> best_inliers.size())
        {
            best_inliers.swap(inlier_indices);
        }
    }

    correspondances2D2D corr_f;
    for(int i = 0; i< best_inliers.size(); i++)
    {
        corr_f.push_back(matches[best_inliers[i]]);
    }

    //利用所有内点进行最小二乘估计
    FundamentalMatrix F = calc_fundamental_least_squares(corr_f);

    std::cout<<"inlier number: "<< best_inliers.size()<<std::endl;
    std::cout<<"F\n: "<< F<<std::endl;

    std::cout<<"result should be: \n"
             <<"inliner number: 272\n"
             <<"F: \n"
             <<"-0.00961384 -0.0309071 0.703297\n"
             <<"0.0448265 -0.00158655 -0.0555796\n"
             <<"-0.703477 0.0648517 -0.0117791\n";

}

#endif