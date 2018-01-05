#ifndef APPROACH2_H
#define APPROACH2_H
#include<fstream>
#include<iostream>
#include<map>
#include<set>
#include<algorithm>
#include<queue>
#include<cmath>
#include "Gmatrix.h"
#include "Approach1.h"
#define fout cout

class Approach2{
public:
	Approach1 app1;
	int K,P,w0; double C;
	struct sketch;

	map<pair<long long,long long>,long long> freqCount; //[u,v].freq
	set<long long> vertices;

	vector<int> intersect(vector<int> s1, vector<int> s2);

	map<long long,vector<int>> G1,G2; // V->Z
	int numberofgroups;
	map<int,Gmatrix> Gmatrices; //Z.Gmatrix;
	Gmatrix outlierSketch;

	void lookup(string file);
	pair<double,int> sourcePartitioning(sketch &s);
	pair<double,int> destinationPartitioning(sketch &s);

	bool terminate(sketch &s);
	void construct(sketch &s);

	double divide(double x, double y);

	Approach2(string file, int and_rows, int and_cols, int or_rows, int or_cols, int outlier_rows, int outlier_cols, int depth, int modulo, int w0, double C);

	void add(long long u, long long v, long long freq);
	long long query(long long u, long long v);
};

#endif
