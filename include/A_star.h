#ifndef A_STAR_H
#define A_STAR_H
#include"data_type.h"
#include"bezier_base.h"
using namespace std;
class A_star
{
private:
    std::ofstream outfile;
    const int kCost1 = 10; //直移一格消耗
    const int kCost2 = 14; //斜移一格消耗
    std::multimap<double, GridNodePtr> openList;  //开启列表
    std::vector<GridNodePtr> gridPath;
    vector<vector<GridNode>> MapWithObs;
public:
    vector<vector<GridNode>> GetMap();
    int calcG(GridNode *temp_start, GridNode *gridnode);
    int calcH(GridNode *gridnode, GridNode *end);
    int calcF(GridNode *gridnode);
    vector<GridNodePtr> retrievePath(GridNodePtr current); 
    std::vector<GridNode *> getSurroundNodes(const GridNode *gridnode, bool isIgnoreCorner) ;
    void findPath(GridNode *startPoint, GridNode *endPoint, bool isIgnoreCorner);
    std::vector<GridNode*> GetPath(GridNode &startPoint, GridNode &endPoint, bool isIgnoreCorner);
    bool isCanreach(const GridNode *girdnode, const GridNode *target, bool isIgnoreCorner) ;
   
    void InitMap();
    void AddObs();
    void AddCarToMap(int xmin,int xmax,int ymin,int ymax);
    void AddSpecialCarToMap(int xinit,int yinit);
    void PrintMap();
    vector<GridNode*> AstarSearch(MatrixXd pos);
    void OutputAstarPath(vector<GridNode*>AstarPath);
    A_star(/* args */);
    ~A_star();
};



#endif