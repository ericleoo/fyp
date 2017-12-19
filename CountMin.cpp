#include "CountMin.h"

CountMin::CountMin(int width, int depth) : width(width),depth(depth),P(1000000007){
	if(!width){ isEmpty = true; return;}
	isEmpty = false;
	count = vector<vector<long long>>(width,vector<long long>(depth));
	hashConstants.resize(depth);
	generateRandomHashConstants();
}

CountMin::CountMin(int width, int depth, int P) : width(width),depth(depth),P(P){
	if(!width){ isEmpty = true; return;}
	isEmpty = false;
	count = vector<vector<long long>>(width,vector<long long>(depth));
	hashConstants.resize(depth);
	generateRandomHashConstants();
}

CountMin::CountMin(){}

void CountMin::generateRandomHashConstants(){
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

int CountMin::g(int k, long long v, int MOD){
	/*
	k-th hash function that maps vertex v into 0..colsOD-1
	*/
	long long a = hashConstants[k].first % P;
	long long b = hashConstants[k].second % P;
	v %= P;
	long long res = a;
	res %= P;
	res *= v;
	res %= P;
	res += b;
	res %= P;
	return res % MOD;
}

void CountMin::add(int u, int v, long long val){
	/*
	Add val to the current frequency of edge (i,j)
	*/
	// cout << u << " " << v << " " << val << '\n';

	if(isEmpty) return;

	long long idx = u;

	int cnt = 0, temp = v;
	while(temp) cnt++,temp/=10;

	long long multiplier = 1;
	while(cnt--) multiplier *= 10;

	idx *= multiplier;
	idx += v;

	add(idx,val);
}

void CountMin::add(long long idx, long long val){
	/*
	Add val to the current frequency of edge (i,j)
	*/
	// cout << idx << " " << val << '\n';
	if(isEmpty) return;
	for(int k=0;k<depth;k++){
		count[g(k,idx,width)][k] += val;
	}
}

long long CountMin::query(int u, int v){
	/*
	Add val to the current frequency of edge (i,j)
	*/

	if(isEmpty) return 0;

	long long idx = u;

	int cnt = 0, temp = v;
	while(temp) cnt++,temp/=10;

	long long multiplier = 1;
	while(cnt--) multiplier *= 10;

	idx *= multiplier;
	idx += v;

	return query(idx);
}

long long CountMin::query(long long idx){
	/*
	Returns the frequency of edge (i,j)
	*/

	if(isEmpty) return 0;

	long long ret;
	for(int k=0;k<depth;k++){
		if(k==0) ret = count[g(k,idx,width)][k];
		else ret = min(ret,count[g(k,idx,width)][k]);
	}
	return ret;
}
