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
	void add(long long i, long long j, long long val);
	long long query(long long i,long long j);
	Gmatrix(int rows, int cols, int depth);
	Gmatrix(int rows, int cols, int depth, int P);
	Gmatrix();
};

#endif

