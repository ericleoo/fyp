#ifndef APPROACH1_H
#define APPROACH1_H

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <algorithm>
#include <queue>
#include <cmath>
#include <string>
#include <limits>
#include <utility>
#include <unordered_set>
#include <unordered_map>
#include "Gmatrix.h"
#include "HASH.h"
using namespace std;
#define fout cout

class Approach1
{
  public:
	struct outCmp;
	struct inCmp;

	int K, P, w0, mode;
	double C;
	double outlierPercentage1, outlierPercentage2;

	unordered_map<long long, set<pair<long long, long long>>> outNeighbour, inNeighbour;
	unordered_map<long long, long long> inTotalFreq, outTotalFreq; // V.total freq
	unordered_set<long long> vertices1, vertices2;
	vector<long long> sortedVertices;

	vector<long long> sumDistinctEdges;

	vector<int> intersect(set<int> s1, set<int> s2);

	double divide(double one, double two);

	void lookup(string file);

	struct sketch
	{
		int rows, cols, l, r;
		sketch() {}
		sketch(int a, int b, int c, int d)
		{
			rows = a;
			cols = b;
			l = c;
			r = d;
		}
	};
	double getVars(double &outVar, double &inVar);

	int partition1(int l, int r, vector<long long> &sortedVertices);

	bool terminate(int l, int r, int rows, int cols);

	unordered_map<long long, int> G; // V.Z
	int numberofgroups;
	unordered_map<int, Gmatrix> Gmatrices; //Z.Gmatrix;
	Gmatrix outlierSketch;

	void construct(int l, int r, vector<long long> &v, int rows, int cols);

	void reset();

	long long getDistinctEdges(int l, int r);

	Approach1(string data_sample_file, int rows, int cols, int depth, int modulo, int w0, double C);
	Approach1();

	void add(long long u, long long v, long long freq);
	long long query(long long u, long long v);

	void clearAll();
	bool sorting(string &s);

	void setup(string data_sample_file, int rows, int cols, int depth, int modulo, int w0, double C);

	double getPercentage(string s);
	int getMode(string s);

	long long queryNodeAggFreq(int x);

	unordered_set<pair<int, int>,HASH> getHeavyHitterEdges(long long F);
};

#endif
