#ifndef _DATA_TYPE_H
#define _DATA_TYPE_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <map>
#include <utility>

#define inf 1>>30
#define X_size 16
#define Y_size 200
#define _max_inflate_iter_x 20//x方向膨胀最大迭代次数
#define _max_inflate_iter 200//迭代最大次数
#define step_length 0.5
#define corrsidorstep 1//有了判断node是不是在corridor内部的逻辑后，就不需要这个step了
#define vecTimeAllocation  20
#define accTimeAllocation 4.0
#define min_center_distance 4.0
# define _GLIBCXX_USE_CXX11_ABI 1
// #define _minimize_j 3.0//公式中的k，到jerk
// #define _minimize_a 2.0//公式中的k，到a
// #define _minimize_v 1.0//公式中的k，到v

struct Square;
struct GridNode;
typedef GridNode* GridNodePtr;
struct GridNode
{     
   int id=0;        // 1--> open set, -1 --> closed set
   Eigen::Vector2d coord;//坐标
   int  index;//是否被占用，1=有障碍物，0=无障碍物
   std::multimap<double, GridNodePtr>::iterator nodeMapIt;
   double G,F,H;
   GridNode*parent;
   double occupancy;
   GridNode( Eigen::Vector2d _coord)
   {  
      id = 0;
      index = 0;
      coord = _coord;

      G = 0;
      F = 0;
      parent = NULL;
   }
   bool operator==(GridNode gridnode){
      return coord==gridnode.coord;
   }
   GridNode(){};
   
   ~GridNode(){};
};

struct Square
{   
   /*
   p1 _ _ _ _ _ _ _ _ _ _p2
      |                  |
      |                  |
      |                  |
      |                  | 
      |                  |
   p4 |__ _ _ _ _ _ _ _ _| p3  -

y  |
   |
   |
   |
   |
   |-------------->x*/
      std::pair<double,double>p1;  
      std::pair<double,double>p2;  
      std::pair<double,double>p3;  
      std::pair<double,double>p4;  
      std::pair<double,double>center;
      std::vector< std::pair<double, double> > box;
      bool valid;    // indicates whether this cube should be deleted

      double t; // time allocated to this square                                                                        
      void PrintCenter()
      {
         //std::cout<<"center of the suqare: \n"<<center<<std::endl;            
      }
      bool operator==(Square sq){
         bool xequal=(p1.first==sq.p1.first&&p2.first==sq.p2.first&&p3.first==sq.p3.first&&p4.first==sq.p4.first);
         bool yequal=(p1.second==sq.p1.second&&p2.second==sq.p2.second&&p3.second==sq.p3.second&&p4.second==sq.p4.second);
         return yequal&&xequal;

      }
      void SetBox(){//用与设置square中心点
        center.first=(double)(p1.first+p2.first+p3.first+p4.first)/4;
        center.second=(double)(p1.second+p2.second+p3.second+p4.second)/4;
        box.clear();
        box.resize(2);
        box[0]=std::make_pair(p4.first,p2.first);
        box[1]=std::make_pair(p4.second,p2.second);
      }

      Square()
      {  
         center ={0,0};
         p1={0,0};p2={0,0};p3={0,0};p4={0,0};
         valid = true;
         t = 0.0;
      }
      Square(GridNode* gridnode)
      {
         p1={gridnode->coord[0],gridnode->coord[1]};
         p2={gridnode->coord[0],gridnode->coord[1]};
         p3={gridnode->coord[0],gridnode->coord[1]};
         p4={gridnode->coord[0],gridnode->coord[1]};
         center={gridnode->coord[0],gridnode->coord[1]};
         valid = true;
         t = 0.0;
      }

      ~Square(){}
};



#endif