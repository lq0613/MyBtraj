#include"timeallocation.h"

using namespace std;
//匀速分配时间，和Btraj有点不同，目前没看懂他的逻辑
//起点时间不能为0,会影响后续求解！每个点中存储的时间都是上一个节点到这里需要的时间（不叠加）
//第一个corridor的t是从起点到第二个square的时间
void timeAllocation(vector<Square>&corridor,const MatrixXd &pos,bool IsAverage){
    int corridorsize=corridor.size();
    double _start_x=pos(0,0);
    double _start_y=pos(0,1);
    double _end_x  =pos(1,0);
    double _end_y  =pos(1,1);
    double distance=0.0;
if(IsAverage){
    if(corridorsize==0){
        cout<<"corridor num=0"<<endl;
        return ;
    }
    else if(corridorsize==1){
        distance=(_start_x-_end_x)*(_start_x-_end_x)+(_start_y-_end_y)*(_start_y-_end_y);
        distance=sqrt(distance);
        corridor[0].t=distance/vecTimeAllocation;
        
    }
    else{
        distance=(_start_x-corridor[1].center.first)*(_start_x-corridor[1].center.first)
                +(_start_y-corridor[1].center.second)*(_start_y-corridor[1].center.second);
        distance=sqrt(distance);
        corridor[0].t=distance/vecTimeAllocation;
        int _end_idx=corridorsize-1;
        distance=(_end_x-corridor[_end_idx].center.first)*(_end_x-corridor[_end_idx].center.first)
                +(_end_y-corridor[_end_idx].center.second)*(_end_y-corridor[_end_idx].center.second);
        distance=sqrt(distance);
        corridor[_end_idx].t=distance/vecTimeAllocation;
        for(int i=1;i<_end_idx;i++)
        {
            distance=(corridor[i].center.first-corridor[i+1].center.first)*
                        (corridor[i].center.first-corridor[i+1].center.first)+
                        (corridor[i].center.second-corridor[i+1].center.second)*
                        (corridor[i].center.second-corridor[i+1].center.second);
            distance=sqrt(distance);
            corridor[i].t=distance/vecTimeAllocation;
        }
    }

    }
else{//梯形分配，加速度恒定
    double d0=vecTimeAllocation*vecTimeAllocation/accTimeAllocation;
    if(corridorsize==0){
    cout<<"corridor num=0"<<endl;
    return ;
    }
    else if(corridorsize==1){
        distance=(_start_x-_end_x)*(_start_x-_end_x)+(_start_y-_end_y)*(_start_y-_end_y);
        distance=sqrt(distance);
        //下面这两句用vt图推导的，加速度恒定
        if(distance<=d0) corridor[0].t=2*sqrt(distance/accTimeAllocation);
        else             corridor[0].t=vecTimeAllocation/accTimeAllocation+distance/vecTimeAllocation; 
    }
    else{
        distance=(_start_x-corridor[1].center.first)*(_start_x-corridor[1].center.first)
                +(_start_y-corridor[1].center.second)*(_start_y-corridor[1].center.second);
        distance=sqrt(distance);
        if(distance<=d0) corridor[0].t=2*sqrt(distance/accTimeAllocation);
        else             corridor[0].t=vecTimeAllocation/accTimeAllocation+distance/vecTimeAllocation; 

        int _end_idx=corridorsize-1;
        distance=(_end_x-corridor[_end_idx].center.first)*(_end_x-corridor[_end_idx].center.first)
                +(_end_y-corridor[_end_idx].center.second)*(_end_y-corridor[_end_idx].center.second);
        distance=sqrt(distance);
        if(distance<=d0) corridor[_end_idx].t=2*sqrt(distance/accTimeAllocation);
        else             corridor[_end_idx].t=vecTimeAllocation/accTimeAllocation+distance/vecTimeAllocation; 

        for(int i=1;i<_end_idx;i++)
        {
            distance=(corridor[i].center.first-corridor[i+1].center.first)*
                        (corridor[i].center.first-corridor[i+1].center.first)+
                        (corridor[i].center.second-corridor[i+1].center.second)*
                        (corridor[i].center.second-corridor[i+1].center.second);
            distance=sqrt(distance);
        if(distance<=d0) corridor[i].t=2*sqrt(distance/accTimeAllocation);
        else             corridor[i].t=vecTimeAllocation/accTimeAllocation+distance/vecTimeAllocation; 
        }
    }

    }
    

    //lq:make t=1 to disappear the influence of the scale
    //update 11.18:whatever t is, the curve is ok. But t=1 will be easy to 
    //understand the variables
    //update 11.22:the time allocation decides the shape of the curve
    //but i don't know how.
    //  for(int i=0;i<corridorsize;i++)
    //  {
    //      corridor[i].t=1;
    //  }
    return;
    //corridor[0].t=1;

}

