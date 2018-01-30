#include "Gmatrix.h"
//random_device random;

Gmatrix::Gmatrix(int rows, int cols, int depth) : rows(rows), cols(cols), depth(depth), P(1000000007)
{
	if (!rows || !cols)
	{
		isEmpty = true;
		return;
	}
	isEmpty = false;
	count = vector<vector<vector<long long>>>(rows, vector<vector<long long>>(cols, vector<long long>(depth)));
	hashConstants.resize(depth);
	generateRandomHashConstants();
}

Gmatrix::Gmatrix(int rows, int cols, int depth, int P) : rows(rows), cols(cols), depth(depth), P(P)
{
	if (!rows || !cols)
	{
		isEmpty = true;
		return;
	}
	isEmpty = false;
	count = vector<vector<vector<long long>>>(rows, vector<vector<long long>>(cols, vector<long long>(depth)));
	hashConstants.resize(depth);
	generateRandomHashConstants();
}

Gmatrix::Gmatrix() {}

void Gmatrix::generateRandomHashConstants()
{
	/*
	Generate pairwise distinct pairs (a,b) for the hash function
	*/
	random_device random;
	set<pair<int, int>> abpairs;
	generator = mt19937(random());
	uniform = uniform_int_distribution<>(1, P - 1);
	while ((int)abpairs.size() < depth)
	{
		int a = uniform(generator);
		int b = uniform(generator);
		abpairs.insert({a, b});
	}
	for (int i = 0; i < depth; i++)
	{
		hashConstants[i] = *abpairs.begin();
		abpairs.erase(abpairs.begin());
	}
}

int Gmatrix::g(int k, long long v, int MOD)
{
	/*
	k-th hash function that maps vertex v into 0..colsOD-1
	*/
	long long a = hashConstants[k].first;
	long long b = hashConstants[k].second;
	long long res = a;
	check(res,P);
	res *= v;
	check(res,P);
	res += b;
	check(res,P);
	check(res,MOD);
	return res;
}

long long Gmatrix::check(long long &x,long long &y){
    (((x%=y)+=y)%=y);
}

vector<int> Gmatrix::gi(int w, int g, int MOD){
	long long a = hashConstants[w].first;
	long long b = hashConstants[w].second;
	vector<int> ret;
	for(int k=-(g/MOD);k<=(P-1-g)/MOD;k++){
		long long temp = k % P;
		temp *= MOD;
		check(temp,P);
		temp += g;
		check(temp,P);
		temp -= b;
		check(temp,P);
		temp *= modInverse(a,P);
		check(temp,P);
		ret.push_back(temp);
	}
	sort(ret.begin(),ret.end());
	return ret;
}

void Gmatrix::add(long long i, long long j, long long val)
{
	/*
	Add val to the current frequency of edge (i,j)
	*/
	if (isEmpty)
		return;
	for (int k = 0; k < depth; k++)
	{
		count[g(k, i, rows)][g(k, j, cols)][k] += val;
	}
}

long long Gmatrix::query(long long i, long long j)
{
	/*
	Returns the frequency of edge (i,j)
	*/
	if (isEmpty)
		return 0;
	long long ret;
	for (int k = 0; k < depth; k++)
	{
		if (k == 0)
			ret = count[g(k, i, rows)][g(k, j, cols)][k];
		else
			ret = min(ret, count[g(k, i, rows)][g(k, j, cols)][k]);
	}
	return ret;
}

// Function to find modulo inverse of b. It returns
// -1 when inverse doesn't
long long Gmatrix::modInverse(long long b, long long m)
{
	long long x, y; // used in extended GCD algorithm
	long long g = gcdExtended(b, m, &x, &y);
	// Return -1 if b and m are not co-prime
	if (g != 1)
		return -1;
	// m is added to handle negative x
	return (x % m + m) % m;
}

// C function for extended Euclidean Algorithm (used to
// find modular inverse.
long long Gmatrix::gcdExtended(long long a, long long b, long long *x, long long *y)
{
	// Base Case
	if (a == 0)
	{
		*x = 0, *y = 1;
		return b;
	}
	long long x1, y1; // To store results of recursive calong long
	long long gcd = gcdExtended(b % a, a, &x1, &y1);
	// Update x and y using results of recursive
	// calong long
	*x = y1 - (b / a) * x1;
	*y = x1;

	return gcd;
}

void Gmatrix::recurse(vector<int> U, vector<int> V, int k, set<pair<int,int>> &E, vector<set<int>> &S, vector<set<int>> &D){
	if(k >= depth){
		vector<pair<int,int>> temp(U.size() * V.size());
		int i = 0;
		for(auto it:U) for(auto it2:V){
			temp[i++] = {it,it2};
		}
		E.insert(temp.begin(),temp.end());
		return;
	}
	for(auto it:Q[k]){
		vector<int> X(it.first.size());

		int idx = 0;

		for(auto u:it.first){
			bool ok = true;
			for(int i=0;i<depth;i++){
				if(!S[i].count(u)){
					ok = false;
					break;
				}
			}
			if(ok) X[idx++] = u;
		}

		X.resize(idx);

		vector<int> Y(it.second.size());
		idx = 0;
		for(auto v:it.second){
			bool ok = true;
			for(int i=0;i<depth;i++){
				if(!D[i].count(v)){
					ok = false;
					break;
				}
			}
			if(ok) Y[idx++] = v;
		}

		Y.resize(idx);
		
		vector<int> nU(U.size() + X.size());
		vector<int> nV(V.size() + Y.size());

		nU.resize(set_intersection (U.begin(), U.end(), X.begin(), X.end(), nU.begin()) - nU.begin());
		if((int)nU.size() == 0) continue;

		nV.resize(set_intersection (V.begin(), V.end(), Y.begin(), Y.end(), nV.begin()) - nV.begin());
		if((int)nV.size() == 0) continue;

		recurse(nU,nV,k+1);
	}
}

set<pair<int,int>> Gmatrix::getHeavyHitterEdges(long long F){
	vector<pair<vector<int>,vector<int>>> Q[depth];
	vector<set<int>> S(depth), D(depth);

	for(int i=0;i<rows;i++) for(int j=0;j<cols;j++) for (int k = 0; k < depth; k++){
		if(count[i][j][k] >= F){
			vector<int> U = gi(i,rows), V = gi(j,cols);
	
			S[k].insert(U.begin(),U.end());
			D[k].insert(V.begin(),V.end());

			Q[k].push_back({U,V});
		}
	}

	set<pair<int,int>> ret;

	for(auto it:Q[0]){
		vector<int> X(it.first.size());
		int idx = 0;

		for(auto u:it.first){
			bool ok = true;
			for(int i=0;i<depth;i++){
				if(!S[i].count(u)){
					ok = false;
					break;
				}
			}
			if(ok) X[idx++] = u;
		}

		X.resize(idx);

		vector<int> Y(it.second.size());
		idx = 0;
		for(auto v:it.second){
			bool ok = true;
			for(int i=0;i<depth;i++){
				if(!D[i].count(v)){
					ok = false;
					break;
				}
			}
			if(ok) Y[idx++] = v;
		}
		Y.resize(idx);
		recurse(X,Y,1,ret,S,D);
	}

	return ret;
}