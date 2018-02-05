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
	count = vector<vector<vector<long long>>>(depth, vector<vector<long long>>(rows, vector<long long>(cols)));
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
	count = vector<vector<vector<long long>>>(depth, vector<vector<long long>>(rows, vector<long long>(cols)));
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
	check(res, P);
	res *= v;
	check(res, P);
	res += b;
	check(res, P);
	check(res, MOD);
	return res;
}

vector<int> Gmatrix::gi(int w, int g, int MOD)
{
	long long a = hashConstants[w].first;
	long long b = hashConstants[w].second;
	vector<int> ret;
	for (int k = -(g / MOD); k <= (P - 1 - g) / MOD; k++)
	{
		long long temp = k % P;
		temp *= MOD;
		check(temp, P);
		temp += g;
		check(temp, P);
		temp -= b;
		check(temp, P);
		temp *= modInverse(a, P);
		check(temp, P);
		ret.push_back(temp);
	}
	sort(ret.begin(), ret.end());
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
		count[k][g(k, i, rows)][g(k, j, cols)] += val;
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
			ret = count[k][g(k, i, rows)][g(k, j, cols)];
		else
			ret = min(ret, count[k][g(k, i, rows)][g(k, j, cols)]);
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


void Gmatrix::recurse(vector<int> U, vector<int> V, int k, unordered_set<pair<int, int>, HASH> &E, vector<vector<pair<vector<int>, vector<int>>>> &Q)
{
	if (k >= depth)
	{
		for (auto it : U)
			for (auto it2 : V)
				E.insert({it, it2});
		return;
	}

	for (auto it : Q[k])
	{
		vector<int> nU;
		vector<int> nV;
        
        set_intersection(U.begin(), U.end(), it.first.begin(), it.first.end(), back_inserter(nU));
        
        if ((int)nU.size() == 0)
			continue;
        
        set_intersection(V.begin(), V.end(), it.second.begin(), it.second.end(), back_inserter(nV));
        
		if ((int)nV.size() == 0)
			continue;

		recurse(nU, nV, k + 1, E, Q);
	}
}

unordered_set<pair<int, int>, HASH> Gmatrix::getHeavyHitterEdges(long long F)
{
	vector<vector<pair<vector<int>, vector<int>>>> Q(depth);
    vector<int> S, D;
    
    cout << "\nGetting source and destination set...\n";
    
    ProgressBar bar(5);
    
    bar.SetNIter(depth * rows * cols);
    bar.Reset();
    
	for (int k = 0; k < depth; k++)
	{
        bar.Update();
        set<int> tS, tD;
        
		for (int i = 0; i < rows; i++){
			for (int j = 0; j < cols; j++){
                bar.Update();
				if (count[k][i][j] >= F)
				{
					vector<int> U = gi(k, i, rows);
					vector<int> V = gi(k, j, cols);

					tS.insert(U.begin(), U.end());
					tD.insert(V.begin(), V.end());
                    
				}
            }
        }
        
        if(k == 0){
            S = vector<int>(tS.begin(),tS.end());
            D = vector<int>(tD.begin(),tD.end());
        }
        else{
            S.resize((int)(
                set_intersection(S.begin(),S.end(),tS.begin(),tS.end(),S.begin()) - S.begin()
            ));
            
            D.resize((int)(
                set_intersection(D.begin(),D.end(),tD.begin(),tD.end(),D.begin()) - D.begin()
            ));
        }
	}
	
	cout << "\nGetting pairs set..\n";
    
    bar.SetNIter(depth*rows*cols);
    bar.Reset();
	
	for (int k=0;k<depth;k++){
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                bar.Update();
                if(count[k][i][j] >= F)
                {
                    vector<int> U = gi(k, i, rows);
					vector<int> V = gi(k, j, cols);
                    
                    U.resize((int)(
                        set_intersection(U.begin(),U.end(),S.begin(),S.end(),U.begin()) - U.begin()
                    ));
                    
                    V.resize((int)(
                        set_intersection(V.begin(),V.end(),D.begin(),D.end(),V.begin()) - V.begin()
                    ));
                    
                    Q[k].push_back(pair<vector<int>, vector<int>>(U, V));
                }
            }
        }
    }    
    
    for (int k=0;k<depth;k++) count[k].clear();
    
    S.clear();
    D.clear();

	unordered_set<pair<int, int>, HASH> ret;
    
    cout << "\nRecursing..\n";
    int it = 0;
    
    bar.SetNIter((int)(Q[0].size()));
    bar.Reset();
    
	for (auto it : Q[0]){
		bar.Update();
        recurse(it.first, it.second, 1, ret, Q);
    }
	return ret;
}
