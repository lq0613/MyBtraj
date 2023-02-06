#ifndef _TRAJECTORY_GENERATOR_H_
#define _TRAJECTORY_GENERATOR_H_


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
// #include "mosek.h"
#include "/usr/local/include/osqp/osqp.h"//自己电脑中的osqp
#include "bezier_base.h"
#include "data_type.h"
#include"timeallocation.h"

using namespace std;
using namespace Eigen;

class TrajectoryGenerator {
private:

public:
        TrajectoryGenerator(){}
        ~TrajectoryGenerator(){}

        /* Use Bezier curve for the trajectory */
       int BezierPloyCoeffGeneration(
            const vector<Square> &corridor,
            const MatrixXd &MQM_j,
            const MatrixXd &MQM_a,
            const MatrixXd &MQM_v,
            const double kj,
            const double ka,
            const double kv,
            const MatrixXd &pos,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,
            const double margin,
            const bool & isLimitVel,
            const bool & isLimitAcc,
            MatrixXd & PolyCoeff);
        Vector2d getPosFromBezier(const MatrixXd & polyCoeff,Square square, double t_now, int seg_now,int _traj_order,VectorXd _C );
        void GetandOutPath(const vector<Square> &corridor,TrajectoryGenerator traj,const MatrixXd &_bezier_coeff,int _traj_order,VectorXd _C);
        
protected:

        virtual OSQPSettings* SolverDefaultSettings();
        
        void FreeData(OSQPData* data);
        

        template <typename T> T* CopyData(const std::vector<T>& vec) {
        T* data = new T[vec.size()];
        memcpy(data, vec.data(), sizeof(T) * vec.size());
        return data;
    }
};
vector<vector<double>> CSCReduction(vector<c_float> A_data,vector<c_int> A_indices,vector<c_int> A_indptr);

void OutputToCSV(vector<vector<double>>A);

#endif
