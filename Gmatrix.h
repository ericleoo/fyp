#ifndef GMATRIX_H
#define GMATRIX_H
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <random>
#include <iostream>
#include <utility>
#include <string>
#include "HASH.h"
using namespace std;

#define check(x, y) (((x %= y) += y) %= y)

class Gmatrix
{
  private:
	mt19937 generator;
	uniform_int_distribution<> uniform;
	void generateRandomHashConstants();
	vector<vector<vector<long long>>> count;
	bool isEmpty;

  public:
	int rows, cols, depth, P;
	vector<pair<int, int>> hashConstants;
	int g(int k, long long v, int MOD);
	vector<int> gi(int depth, int g, int MOD);
	void add(long long i, long long j, long long val);
	long long query(long long i, long long j);
	Gmatrix(int rows, int cols, int depth);
	Gmatrix(int rows, int cols, int depth, int P);
	Gmatrix();
	long long modInverse(long long b, long long m);
	long long gcdExtended(long long a, long long b, long long *x, long long *y);

	void recurse(vector<int> U, vector<int> V, int k, unordered_set<pair<int, int>,HASH> &E, vector<vector<pair<vector<int>, vector<int>>>> &Q);
	unordered_set<pair<int, int>,HASH> getHeavyHitterEdges(long long F);

	void logging(string s);

	long long queryNodeOutgoingFreq(int u);
	long long queryNodeIncomingFreq(int v);

};

#endif
