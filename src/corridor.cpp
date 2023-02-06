#include"corridor.h"
pair<Square, bool> inflateSquare(Square square,vector<vector<GridNode>>const&MapWithObs){
    //Square MAXsquare=square;
    bool collide=false;
    if(MapWithObs[square.p1.first*2][square.p1.second*2].index==1){//传进去的是点，判断一个顶点就可以了
        return make_pair(square,false);//如果传进去的点在障碍物内，返回false。
    }
    int iter = 0;
    int x_iter=0;
    double id_x,id_y;
    while(iter < _max_inflate_iter) // 迭代次数也就是最大扩展距离
    {
        Square presquare=square;
        double x_lo=max(0.0,square.p4.first-step_length);
        double x_up=min(X_size-1.0,square.p3.first+step_length);
        
        // Y Axis
        double y_lo = max(0.0,square.p4.second-step_length);
        //最大值\最小值决定了只走一圈，Y_size是处理边界条件用的
        double y_up = min(Y_size-1.0,square.p1.second+step_length);
        // Y+ p1.p2
        // ############################################################################################################
        collide = false;
        for(id_y = square.p1.second; id_y <= y_up; id_y+=step_length )//这里不用for循环会更好理解
        {   
            if( collide == true) 
                break;
            
            for( id_x = square.p2.first; id_x >= square.p1.first; id_x-=step_length ) // p1p2
            {
                int occupy=MapWithObs[id_x*2][id_y*2].index;
                    if(occupy==1) // the voxel is occupied
                    {   
                        collide = true;
                        break;
                    }
            }
        }
        

        if(collide)
        {
            square.p1.second = max(id_y-2*step_length, square.p1.second); 
            square.p2.second = max(id_y-2*step_length, square.p2.second);   // _step_length = 1, 若有障碍说明之前cube就已经到达边界 
        }
        else
            square.p1.second=square.p2.second= id_y - step_length;  // 没有碰撞，膨胀成功

        //Y- p4p3
        // ############################################################################################################
        collide = false;
        for(id_y = square.p4.second; id_y>=y_lo; id_y-=step_length )//这里不用for循环会更好理解
        {   
            if( collide == true) 
                break;
            
            for( id_x = square.p4.first; id_x <= square.p3.first; id_x+=step_length ) // p4p3
            {
                int occupy=MapWithObs[id_x*2][id_y*2].index;
                    if(occupy==1) // the voxel is occupied
                    {   
                        collide = true;
                        break;
                    }
            }
            
        }

        if(collide)
        {
            square.p4.second = max(id_y+2*step_length, square.p4.second); 
            square.p3.second = max(id_y+2*step_length, square.p3.second);   // _step_length = 1, 若有障碍说明之前cube就已经到达边界 
        }
        else
            square.p4.second=square.p3.second= id_y + step_length;  // 没有碰撞，膨胀成功
        
        //X+ p2p3
        // ############################################################################################################
        if(x_iter<_max_inflate_iter_x){
        collide = false;
        for(id_x = square.p2.first; id_x<=x_up; id_x+=step_length )//这里不用for循环会更好理解
        {   
            if( collide == true) 
                break;
            
            for( id_y = square.p3.second; id_y <= square.p2.second; id_y+=step_length ) // p4p3
            {
                int occupy=MapWithObs[id_x*2][id_y*2].index;
                    if(occupy==1) // the voxel is occupied
                    {   
                        collide = true;
                        break;
                    }
            }
            
        }

        if(collide)
        {
            square.p2.first = max(id_x-2*step_length, square.p2.first); 
            square.p3.first = max(id_x-2*step_length, square.p3.first);   // _step_length = 1, 若有障碍说明之前cube就已经到达边界 
        }
        else
            square.p2.first=square.p3.first= id_x - step_length;  // 没有碰撞，膨胀成功
        //X- p1p4
        // ############################################################################################################
        collide = false;
        for(id_x = square.p1.first; id_x>=x_lo; id_x-=step_length )//这里不用for循环会更好理解
        {   
            if( collide == true) 
                break;
            
            for( id_y = square.p4.second; id_y <= square.p1.second; id_y+=step_length ) // p4p3
            {
                int occupy=MapWithObs[id_x*2][id_y*2].index;
                    if(occupy==1) // the voxel is occupied
                    {   
                        collide = true;
                        break;
                    }
            }
            
        }

        if(collide)
        {
            square.p4.first = max(id_x+2*step_length, square.p4.first); 
            square.p1.first = max(id_x+2*step_length, square.p1.first);   // _step_length = 1, 若有障碍说明之前cube就已经到达边界 
        }
        else
            square.p4.first=square.p1.first= id_x + step_length;  // 没有碰撞，膨胀成功
        }

        if(presquare==square){//膨胀到了最大边界，square不变化了
            break;
        }

        iter++;
        x_iter++;     
    }
    return make_pair(square,true);
}
vector<Square>CorridorGenerate(vector<GridNode*>AstarPath,vector<vector<GridNode>>const&MapWithObs){
    vector<Square> SquareList;
    GridNode* gridnode;
    Square lstSquare;
    int counttmp=0;
    //这里控制步长，减少运算量
    for (int i = 0; i < (int)AstarPath.size(); i++)
    {
        gridnode = AstarPath[i];//类型是girdnode*
        double x=gridnode->coord[0];
        double y=gridnode->coord[1];
         //如果Astar上的点在上一个square的内部，不要,节约时间
        if(counttmp>0&&IsInSquare(gridnode,SquareList[counttmp-1])) continue;

        Square square(gridnode);//使用构造函数直接由gridnode生产square       
        auto result = inflateSquare(square, MapWithObs);
        square = result.first;
        square.SetBox();

        if(result.second == false) continue;
        //如果两个square一样，就不要
        if(lstSquare==square) continue;      
        //如果被包含,或者v重合区域过高，也不要
        if(counttmp>0&&IsContain(SquareList[counttmp-1],square)) continue;

        lstSquare = square;
        SquareList.push_back(square);        
        counttmp++;
    }
    std::cout<<"squarenum_beforetrim="<<counttmp<<endl;
    SquareList=CorridorTrim(SquareList);
    OutputCorridorToCSV(SquareList);
    std::cout<<"squarenum_aftertrim="<<SquareList.size()<<endl;

    return SquareList;
}
void PrintCorridor(const vector<Square>&corridor){
    int n=corridor.size();
    for(int i=0;i<n;i++)
    {
      cout<<corridor[i].p1.first<<"  "<<corridor[i].p2.first<<"  "<<corridor[i].p3.first<<"  "<<corridor[i].p4.first<<endl;
      cout<<corridor[i].p1.second<<"  "<<corridor[i].p2.second<<"  "<<corridor[i].p3.second<<"  "<<corridor[i].p4.second<<endl;
      cout<<corridor[i].center.first<<"  "<<corridor[i].center.second<<"  "<<corridor[i].t<<endl;
      cout<<endl;
        
    }
    return ;

}
bool IsContain(const Square &fsq,const Square &bsq){
    if(fsq.p4.first<=bsq.p4.first&&fsq.p4.second<=bsq.p4.second)
    {
        if(fsq.p2.first>=bsq.p2.first&&fsq.p2.second>=bsq.p2.second)
        {
            return true;
        }
    }
    double center_distance=0.0;
    center_distance=(fsq.center.first-bsq.center.first)*(fsq.center.first-bsq.center.first)+
                    (fsq.center.second-bsq.center.second)*(fsq.center.second-bsq.center.second);
    center_distance=sqrt(center_distance);  
    // 可能需要改成面积重复率，不过先这样吧  
    if(center_distance<min_center_distance) return true;
    return false;
}
bool IsOverlap(const Square &fsq,const Square &bsq){
    double flength=fsq.box[0].second-fsq.box[0].first;
    double fwidth=fsq.box[1].second-fsq.box[1].first;
    double blength=bsq.box[0].second-bsq.box[0].first;
    double bwidth=bsq.box[1].second-bsq.box[1].first;
    double CenterDistance_x=abs(fsq.center.first-bsq.center.first);
    double CenterDistance_y=abs(fsq.center.second-bsq.center.second);
    if(2*CenterDistance_x>=(flength+blength)||2*CenterDistance_y>=(fwidth+bwidth)) return false;
    else return true;
}
bool IsInSquare(const GridNode *gridnode,const Square&square){
    if(gridnode->coord[0]<square.p2.first&&gridnode->coord[0]>square.p4.first)
    {
        if(gridnode->coord[1]<square.p2.second&&gridnode->coord[1]>square.p4.second)
        {
            return true;
        }
    }
    return false;
}
vector<Square> CorridorTrim(vector<Square>&corridor){
    int n=corridor.size();
    if(n==1||n==0) return corridor;
    vector<Square> ret;
    Square SquareNow=corridor[0];
    ret.push_back(corridor[0]);
    int i=1;
    while(i<n){
        if(IsOverlap(SquareNow,corridor[i])) i++;
        else{
            SquareNow=corridor[i-1];
            ret.push_back(SquareNow);
        }
    }
    if(IsContain(SquareNow,corridor[n-1])==false)    ret.push_back(corridor[n-1]);
    return ret;
}
void OutputCorridorToCSV(const vector<Square>&corridor){
        //记录corridor边框
    
    ofstream outfile;
    for(int k=0;k<corridor.size();k++){
        Square square=corridor[k];
        vector<double>OutputCorridor;
        for(double i=square.p4.first;i<=square.p2.first;i+=square.p2.first-square.p4.first)
            {
             for(double j=square.p4.second;j<=square.p2.second;j+=1)
                {
                        OutputCorridor.push_back(i);
                        OutputCorridor.push_back(j);
                }
            }
        for(double i=square.p4.second;i<=square.p2.second;i+=square.p2.second-square.p4.second)
            {
             for(double j=square.p4.first;j<=square.p2.first;j+=0.1)
                {
                        OutputCorridor.push_back(j);
                        OutputCorridor.push_back(i);
                }
            }
                string s=to_string(k);
                outfile.open("/home/liuqiao/project/Btraj/MyBtraj 1-12/csv/corridor"+s+".csv");                
                for(int i=0;i<OutputCorridor.size()-1;i+=2) outfile<<OutputCorridor[i]<<","<<OutputCorridor[i+1]<<endl;
                outfile.close();
    }
        
            
        

}