#ifndef GSKETCH_H
#define GSKETCH_H

#include<fstream>
#include<iostream>
#include<map>
#include<set>
#include<algorithm>
#include<queue>
#include<cmath>
#include "CountMin.h"

using namespace std;

#define fout cout

class Gsketch{
public:
	struct outCmp;
	struct inCmp;

	int K,P,w0,mode;
	double C;

	map<int,set<int>> outNeighbour, inNeighbour;
	map<pair<int,int>,long long> freqCount; //[u,v].freq
	map<int,long long> inTotalFreq, outTotalFreq; // V.total freq
	set<int> vertices;
	map<int,int> G; // V.Z
	int numberofgroups;
	map<int,CountMin> Sketches; //Z.Gmatrix;
	CountMin outlierSketch;

	double divide(double one, double two);

	void lookup(string file);

	struct sketch{
		int width,l,r;
		sketch(){}
		sketch(int a, int c, int d){
			width = a;
			l = c; r = d;
		}
	};

	int partition1(int l, int r, vector<int> &sortedVertices, bool sourceNodeGrouping);
	bool terminate(int l, int r, vector<int> &v, int width, bool sourceNodeGrouping);

	void construct(int l, int r, vector<int> &v, int width);

	void reset();

	Gsketch(string data_sample_file, int width, int outlier_width, int depth, int modulo,int w0, double C);

	void add(int u, int v, long long freq);

	long long query(int u, int v);

	Gsketch();
};

#endif
