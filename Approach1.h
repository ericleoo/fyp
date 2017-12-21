#ifndef APPROACH1_H
#define APPROACH1_H

#include<fstream>
#include<iostream>
#include<map>
#include<set>
#include<algorithm>
#include<queue>
#include<cmath>
#include "Gmatrix.h"

using namespace std;
#define fout cout

class Approach1{
public:

	struct outCmp;
	struct inCmp;

	int K,P,w0,mode;
	double C;

	map<long long,set<long long>> outNeighbour, inNeighbour;
	map<pair<long long,long long>,long long> freqCount; //[u,v].freq
	map<long long,long long> inTotalFreq, outTotalFreq; // V.total freq
	set<long long> vertices;
	
	vector<long long> sumDistinctEdges;

	vector<int> intersect(set<int> s1, set<int> s2);

	double divide(double one, double two);

	void lookup(string file);

	struct sketch{
		int rows,cols,l,r;
		sketch(){}
		sketch(int a, int b, int c, int d){
			rows = a; cols = b;
			l = c; r = d;
		}
	};

	int partition1(int l, int r, vector<long long> &sortedVertices, bool sourceNodeGrouping);

	bool terminate(int l, int r, vector<long long> &v, int rows, int cols, bool sourceNodeGrouping);

	map<long long,int> G; // V.Z
	int numberofgroups;
	map<int,Gmatrix> Gmatrices; //Z.Gmatrix;
	Gmatrix outlierSketch;

	void construct(int l, int r, vector<long long> &v, int rows, int cols);

	void reset();
	
	long long getDistinctEdges(int l, int r);

	Approach1(string data_sample_file, int rows, int cols, int outlier_rows, int outlier_cols, int depth, int modulo, int w0, double C);
	Approach1();

	void add(long long u, long long v, long long freq);
	long long query(long long u, long long v);
};

#endif
