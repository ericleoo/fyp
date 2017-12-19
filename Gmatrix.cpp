#include "Gmatrix.h"
//random_device random;

Gmatrix::Gmatrix(int rows, int cols, int depth) : rows(rows),cols(cols),depth(depth),P(1000000007){
	if(!rows || !cols){ isEmpty = true; return;}
	isEmpty = false;
	count = vector<vector<vector<long long>>>(rows,vector<vector<long long>>(cols,vector<long long>(depth)));
	hashConstants.resize(depth);
	generateRandomHashConstants();
}

Gmatrix::Gmatrix(int rows, int cols, int depth, int P) : rows(rows),cols(cols),depth(depth),P(P){
	if(!rows || !cols){ isEmpty = true; return;}
	isEmpty = false;
	count = vector<vector<vector<long long>>>(rows,vector<vector<long long>>(cols,vector<long long>(depth)));
	hashConstants.resize(depth);
	generateRandomHashConstants();
}

Gmatrix::Gmatrix(){}

void Gmatrix::generateRandomHashConstants(){
	/*
	Generate pairwise distinct pairs (a,b) for the hash function
	*/
	random_device random;
	set<pair<int,int>> abpairs;
	generator = mt19937(random());
	uniform = uniform_int_distribution<>(1,P-1);
	while((int)abpairs.size() < depth){
		int a = uniform(generator);
		int b = uniform(generator);
		abpairs.insert({a,b});
	}
	for(int i=0;i<depth;i++){
		hashConstants[i] = *abpairs.begin();
		abpairs.erase(abpairs.begin());
	}
}

int Gmatrix::g(int k, long long v, int MOD){
	/*
	k-th hash function that maps vertex v into 0..colsOD-1
	*/
	long long a = hashConstants[k].first;
	long long b = hashConstants[k].second;
	long long res = (a * v) + b;
	return (res % P) % MOD;
}

void Gmatrix::add(long long i, long long j, long long val){
	/*
	Add val to the current frequency of edge (i,j)
	*/
	if(isEmpty) return;
	for(int k=0;k<depth;k++){
		count[g(k,i,rows)][g(k,j,cols)][k] += val;
	}
}

long long Gmatrix::query(long long i, long long j){
	/*
	Returns the frequency of edge (i,j)
	*/
	if(isEmpty) return 0;
	long long ret;
	for(int k=0;k<depth;k++){
		if(k==0) ret = count[g(k,i,rows)][g(k,j,cols)][k];
		else ret = min(ret,count[g(k,i,rows)][g(k,j,cols)][k]);
	}
	return ret;
}
