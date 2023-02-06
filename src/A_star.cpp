#include"A_star.h"
A_star::A_star(){};
A_star::~A_star(){};

void A_star::InitMap(){
    for(double i=0;i<X_size;i+=1)
    {
        //先插入，否则vector不能按下标访问
        vector<GridNode>tmp;
        MapWithObs.push_back(tmp);
        for(double j=0;j<Y_size;j+=1)
        {
           GridNode tmp;
           MapWithObs[i].push_back(tmp);
           MapWithObs[i][j].index=0;
           MapWithObs[i][j].coord=Eigen::Vector2d(i,j);

        }
    }
    
}
void A_star:: AddCarToMap(int xmin,int xmax,int ymin,int ymax){
    if(xmin<0||xmax>=X_size||ymin<0||ymax>=Y_size){
        cout<<"the obstacle is out of the map, please try another parameter";
        return ;
    }
        for(int i=xmin;i<xmax;i++)
            {
                for(int j=ymin;j<ymax;j++)
                {
                    MapWithObs[i][j].index=1;
             }
         }
    }
void A_star::AddSpecialCarToMap(int xinit,int yinit){
    if(xinit<0||xinit+4>=X_size||yinit-4<0||yinit>=Y_size){
        cout<<"the obstacle is out of the map, please try another parameter";
        return ;    
    }
    int count=0;
    while(count<4){
        for(int i=xinit,j=yinit;i<xinit+8;i++,j++)
        {
            MapWithObs[i][j].index=1;
        }
        xinit++;
        yinit--;
        count++;
    }
}
void A_star::AddObs(){
     AddCarToMap(10,15,20,30);
     AddCarToMap(5,9,5,15);
      AddCarToMap(7,11,60,70);
      //AddCarToMap(0,10,86,90);
     AddCarToMap(0,5,5,15);
      AddCarToMap(8,13,95,105);
     AddCarToMap(0,5,140,150);
     AddCarToMap(0,5,185,195);
     AddSpecialCarToMap(4,120);  

    outfile.open("/home/liuqiao/project/Btraj/MyBtraj 1-12/csv/obstacle.csv");
    for(double i=0;i<X_size;i+=1)
    {
        for(double j=0;j<Y_size;j+=1)
        {
            if(MapWithObs[i][j].index!=0)
            {
                outfile<<i<<","<<j<<endl;
            }
        }
    }

}
void A_star::PrintMap(){
  for(double i=0;i<X_size;i+=1)
    {
        for(double j=0;j<Y_size;j+=1)
        {
           cout<<MapWithObs[i][j].index;
        }
        cout<<endl;
    }
    return;

}

void A_star::OutputAstarPath(vector<GridNode*>AstarPath){
    ofstream outfile;
    vector<double> OutputAstarPath;
    for(auto pathnode:AstarPath)
    {
        Eigen::Vector2d tmpcoord=pathnode->coord;
        MapWithObs[tmpcoord[0]][tmpcoord[1]].index=2;
        OutputAstarPath.push_back(tmpcoord[0]);
        OutputAstarPath.push_back(tmpcoord[1]);
    }
    outfile.open("/home/liuqiao/project/Btraj/MyBtraj 1-12/csv/astar.csv");
    for(int i=0;i<OutputAstarPath.size()-1;i+=2) outfile<<OutputAstarPath[i]<<","<<OutputAstarPath[i+1]<<endl;
    outfile.close();
    //PrintMap(MapWithObs);
    return;
}
int  A_star::calcG(GridNode *temp_start, GridNode *gridnode){
    int extraG = (abs(gridnode->coord[0] - temp_start->coord[0]) + 
                        abs(gridnode->coord[1] - temp_start->coord[1])) ==1 ? kCost1 : kCost2;
	int parentG = gridnode->parent == NULL ? 0 : gridnode->parent->G; //如果是初始节点，则其父节点是空
	return parentG + extraG;
}
int  A_star::calcH(GridNode *gridnode, GridNode *end){
    // return sqrt((double)(end->coord[0] - gridnode->coord[0])*(double)(end->coord[0] - gridnode->coord[0]) + 
    //             (double)(end->coord[1] - gridnode->coord[1])*(double)(end->coord[1] - gridnode->coord[1]))*kCost1;
    return (abs((double)(end->coord[0] - gridnode->coord[0]))+abs((double)(end->coord[1] - gridnode->coord[1])))*kCost1;
}
int  A_star::calcF(GridNode *gridnode){
	return gridnode->G + gridnode->H;
}
bool A_star::isCanreach(const GridNode *girdnode, const GridNode *target, bool isIgnoreCorner) {
    double x=target->coord[0],y=target->coord[1];
    if (x<0 || x>X_size - 1
		|| y<0 || y>Y_size - 1
		|| target->index == 1
		|| x == girdnode->coord[0]&&y == girdnode->coord[1]
		|| target->id==-1) //如果点与当前节点重合、超出地图、是障碍物、或者在关闭列表中，返回false
        return false;
	else
	{
		if (abs(girdnode->coord[0] - target->coord[0]) + abs(girdnode->coord[1] - target->coord[1]) == 1) //非斜角可以
			return true;
		else
		{
          //斜对角要判断是否绊住
         // int x=
			if (MapWithObs[girdnode->coord[0]][target->coord[1]].index == 0 
                        && MapWithObs[target->coord[0]][girdnode->coord[1]].index == 0)
				return true;
			else
				return isIgnoreCorner;
		}
        return true;
	}
}
void A_star::findPath(GridNode *startPoint, GridNode *endPoint, bool isIgnoreCorner){
    GridNode*tmp=startPoint;
    GridNode*endPointptr=endPoint;
    GridNodePtr curPoint = NULL;
    tmp->parent=NULL;
    tmp->H=calcH(tmp,endPointptr);
    tmp->F=calcF(tmp);
    tmp->id=1;
    openList.insert(make_pair(tmp->F,tmp)); //置入起点,拷贝开辟一个节点，内外隔离
    
	while (!openList.empty())
	{

		curPoint = openList.begin() -> second; //找到F值最小的点
        //cout<<openList.size()<<endl;
        if(curPoint->coord[0]==endPoint->coord[0]&&curPoint->coord[1]==endPoint->coord[1])
        {
            gridPath=retrievePath(curPoint);
            cout<<"overoverover__________________"<<endl;
            return;
        }
        
        openList.erase(openList.begin());
        curPoint->id=-1;//放到关闭列表
		//1,找到当前周围八个格中可以通过的格子
		auto surroundPoints = getSurroundNodes(curPoint, isIgnoreCorner);
       // cout<<surroundPoints.size()<<endl;
		for (auto &target : surroundPoints)
		{
			//2,对某一个格子，如果它不在开启列表中，加入到开启列表，设置当前格为其父节点，计算F G H
            if(target->id!=1)
			{
                target->id=1;
				target->parent = curPoint;
				target->G = calcG(curPoint, target);//这里这个G的计算也可以优化
				target->H = calcH(target, endPoint);
				target->F = calcF(target);
                openList.insert(make_pair(target->F,target)); 
                continue;
			}
			//3，对某一个格子，它在开启列表中，计算G值, 如果比原来的大, 就什么都不做, 否则设置它的父节点为当前点,并更新G和F
			else
			{
				int tempG = calcG(curPoint, target);
				if (tempG<=target->G)
				{
					target->parent = curPoint;
                    target->H=calcH(curPoint,target);
					target->G = tempG;
					target->F = calcF(target);
				}
			}				
		}
	}
     
	return ;
}
std::vector<GridNode *>A_star:: getSurroundNodes(const GridNode *gridnode, bool isIgnoreCorner) {

    std::vector<GridNode *> surroundPoints; 
	for (int x = gridnode->coord[0] - 1; x <= gridnode->coord[0] + 1; x+=1)
    {
        if(x<0||x>X_size-1) continue;
        for (int y = gridnode->coord[1] - 1; y <= gridnode->coord[1] + 1; y+=1)
            { 
                if(y<0||y>Y_size-1) continue;
                //Eigen::Vector2d coord={x,y};
                GridNode*target=&MapWithObs[x][y];
	            if (isCanreach(gridnode, target, isIgnoreCorner))
                {
                    surroundPoints.push_back(target);
                }
		        
            }

    }
	
 
	return surroundPoints;
}
vector<GridNode*> A_star::AstarSearch(MatrixXd pos){

    GridNode* start=&MapWithObs[pos(0,0)][pos(0,1)];
    GridNode* end=&MapWithObs[pos(1,0)][pos(1,1)];
    if(MapWithObs[start->coord[0]][start->coord[1]].index==1||
        MapWithObs[end->coord[0]][end->coord[1]].index==1)
        {
            cout<<"primal path generate failed because start or end point is the obstacle point";
            cout<<endl;
            cout<<"please change start or end point"<<endl;
        }
    findPath(start, end, false);
    
 
    return gridPath;

}
vector<GridNodePtr> A_star::retrievePath(GridNodePtr current)
{   
    vector<GridNodePtr> path;
    path.push_back(current);

    while(current->parent != NULL)
    {
        current = current -> parent;
        path.push_back(current);
    }
    reverse(path.begin(),path.end());
    openList.clear();

    return path;
}
vector<vector<GridNode>>A_star:: GetMap(){
    return MapWithObs;
}