#ifndef GMATRIX_H
#define GMATRIX_H

#include<vector>
#include<set>
#include<fstream>
#include<random>
using namespace std;

class Gmatrix{
private:
	mt19937 generator;
	uniform_int_distribution<> uniform;
	void generateRandomHashConstants();
	vector<vector<vector<long long>>> count;
	bool isEmpty;
public:
	int rows,cols,depth,P;
	vector<pair<int,int>> hashConstants;
	int g(int k, long long v, int MOD);
	vector<int> Gmatrix::gi(int depth, int g, int MOD);
	void add(long long i, long long j, long long val);
	long long query(long long i,long long j);
	Gmatrix(int rows, int cols, int depth);
	Gmatrix(int rows, int cols, int depth, int P);
	Gmatrix();
	long long modInverse(long long b, long long m);
	long long gcdExtended(long long a, long long b, long long *x, long long *y);
	long long check(long long &x,long long &y);

	void gMatrix::recurse(vector<int> U, vector<int> V, int k, set<pair<int,int>> &E, vector<set<int>> &S, vector<set<int>> &D);
	set<pair<int,int>> getHeavyHitterEdges(long long F);
};

#endif

