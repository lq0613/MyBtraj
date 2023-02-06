#ifndef CORRIDOR_H
#define CORRIDOR_H
#include"data_type.h"
using namespace std;
pair<Square, bool> inflateSquare(Square square,vector<vector<GridNode>>const&MapWithObs);
vector<Square> CorridorGenerate(vector<GridNode*>AstarPath,vector<vector<GridNode>>const&MapWithObs);
//void Output(vector<Square> CorridorGenerate);
void PrintCorridor(const vector<Square>&corridor);
void OutputCorridorToCSV(const vector<Square>&corridor);
bool IsContain(const Square &fsq,const Square &bsq);
bool IsOverlap(const Square &fsq,const Square &bsq);
bool IsInSquare(const GridNode*gridnode,const Square&square);
vector<Square> CorridorTrim(vector<Square>&corridor);
#endif