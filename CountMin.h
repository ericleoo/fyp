#ifndef CM_H
#define CM_H

#include<vector>
#include<set>
#include<iostream>
#include<fstream>
#include<random>
using namespace std;

class CountMin{
private:
	mt19937 generator;
	uniform_int_distribution<> uniform;
	void generateRandomHashConstants();
	vector<vector<long long>> count;
	bool isEmpty;
public:
	int width,depth,P;
	vector<pair<int,int>> hashConstants;
	int g(int k, long long v, int MOD);
	void add(int u, int v, long long val);
	long long query(int u, int v);

	void add(long long idx, long long val);
	long long query(long long idx);

	CountMin(int width, int depth);
	CountMin(int width, int depth, int P);
	CountMin();
};

#endif
