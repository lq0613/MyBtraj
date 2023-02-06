//#include"map_generator.hpp"
#include"data_type.h"
#include"A_star.h"
#include"corridor.h"
#include"timeallocation.h"
#include"trajectory_generator.h"
#include<iostream>
#include"bezier_base.h"
#include "signed_distance_field_2d.h"
#include <time.h> 
using namespace Eigen;
//#include"map_generator.hpp"
using namespace std;
using namespace planning;
using OccupancyMap = GridMap2D<uint8_t>;
void LinkGridMap(vector<vector<GridNode>> &GridMap,vector<vector<GridNode>> const &MapWithObs);
int main(int argc,char**argv)
{
    double start=clock();
    A_star astar;
    // GridMap2D<uint8_t> grid_map;
    // grid_map.set_cell_number({8,160});
    // grid_map.set_origin(std::array<double, 2>{0.0, 0.0});
    // grid_map.set_resolution(std::array<double, 2>{1, 1});
    // SignedDistanceField2D sign(std::move(grid_map));
    
    // sign.TestBasic();
    // sign.TestSignedDistance();
    double kj=1.000;//优化函数权重
    double ka=1.0;
    double kv=1.0;
    int _traj_order=5;//阶数  
    double _cube_margin=0.0;//安全距离

    MatrixXd pos = MatrixXd::Zero(2,2);//初始点和终止点的p
    MatrixXd vel = MatrixXd::Zero(2,2);//v
    MatrixXd acc = MatrixXd::Zero(2,2);//a

    pos.row(0)<<0,0;
    pos.row(1)<<10,198;

    vector<vector<GridNode>> GridMap;
    astar.InitMap();
    astar.AddObs();
     start=clock();
    vector<GridNode*>AstarPath=astar.AstarSearch(pos);
    double start2=clock();
    astar.OutputAstarPath(AstarPath); 
    LinkGridMap(GridMap,astar.GetMap());
    vector<Square>corridor=CorridorGenerate(AstarPath,GridMap);
    timeAllocation(corridor,pos,false);
    //PrintCorridor(corridor);

    TrajectoryGenerator traj;
    MatrixXd _bezier_coeff;
    Bernstein _bernstein_j;

    _bernstein_j.setParam(3,12,3.0);
    MatrixXd _MQM_j=_bernstein_j.getMQM()[_traj_order];
    VectorXd _C_j   = _bernstein_j.getC()[_traj_order];
    MatrixXd _Q_j=_bernstein_j.getQ()[_traj_order];
    // cout<<"MQMJ="<< _MQM_j<<endl;
    // cout<<"QJ="<< _Q_j<<endl;

    Bernstein _bernstein_a;
    _bernstein_a.setParam(3,12,2.0);
    MatrixXd _MQM_a=_bernstein_a.getMQM()[_traj_order];
    VectorXd _C_a   = _bernstein_a.getC()[_traj_order];
    MatrixXd _Q_a=_bernstein_a.getQ()[_traj_order];
    // cout<<"MQMa="<< _MQM_a<<endl;
    // cout<<"Qa="<< _Q_a<<endl;

    Bernstein _bernstein_v;
    _bernstein_v.setParam(3,12,1.0);
    MatrixXd _MQM_v=_bernstein_v.getMQM()[_traj_order];
    VectorXd _C_v   = _bernstein_v.getC()[_traj_order];
    MatrixXd _Q_v=_bernstein_v.getQ()[_traj_order];
    // cout<<"MQMv="<< _MQM_v<<endl;
    // cout<<"Qv="<< _Q_v<<endl;

    int i=traj.BezierPloyCoeffGeneration
        ( corridor, _MQM_j,_MQM_a,_MQM_v,kj,ka,kv, pos, vel, acc, 20, 6, _traj_order, 
         _cube_margin, 1,1, _bezier_coeff );

    traj.GetandOutPath(corridor,traj,_bezier_coeff,_traj_order,_C_a);
    cout<<"end"<<endl;
    double end=clock();
    cout<<"time total cost"<<(end-start)/1000000<<endl;
    cout<<"time except A* cost"<<(end-start2)/1000000<<endl;
    return 0;
}
void LinkGridMap(vector<vector<GridNode>> &GridMap,vector<vector<GridNode>> const &MapWithObs){
    for(int i=0;i<2*X_size;i+=1)
    {
        //先插入，否则vector不能按下标访问
        vector<GridNode>tmp;
        GridMap.push_back(tmp);
        for(int j=0;j<2*Y_size;j+=1)
        {
           GridNode tmp;
           GridMap[i].push_back(tmp);
           if(i%2!=0||j%2!=0) GridMap[i][j].index=0;
           else  GridMap[i][j].index=MapWithObs[step_length*i][step_length*j].index;
           GridMap[i][j].coord=Eigen::Vector2d(step_length*i,step_length*j);

        }
    }

}