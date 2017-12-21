#include "Approach1.h"

struct Approach1::outCmp{
	Approach1* m;
	outCmp(Approach1* p) : m(p) {};

	bool operator() ( long long i, long long j )
	{
		if(!m->outNeighbour.count(i) && !m->outNeighbour.count(j)) return i < j;
		if(!m->outNeighbour.count(i)) return true;
		if(!m->outNeighbour.count(j)) return false;
		return m->outTotalFreq[i] * m->outNeighbour[j].size() < m->outTotalFreq[j] * m->outNeighbour[i].size();
	}
};

struct Approach1::inCmp{
	Approach1* m;
	inCmp(Approach1* p) : m(p) {};

	bool operator() ( long long i, long long j )
	{
		if(!m->inNeighbour.count(i) && !m->inNeighbour.count(j)) return i < j;
		if(!m->inNeighbour.count(i)) return true;
		if(!m->inNeighbour.count(j)) return false;
		return m->inTotalFreq[i] * m->inNeighbour[j].size() < m->inTotalFreq[j] * m->inNeighbour[i].size();
	}
};

vector<int> Approach1::intersect(set<int> s1, set<int> s2){
	vector<int> v(s1.size() + s2.size());
	auto it = set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),v.begin());
	v.resize(it-v.begin());
	return v;
}

double Approach1::divide(double one, double two){
	if(two < 1e-9) return 0;
	else return one/two;
}

void Approach1::lookup(string file){
	/*
	Reads data sample. Gathers statistics
	*/
	ifstream fin(file);
	long long u,v; long long freq;
	double temp;

	fout << "LOOKINGUP\n";

	while(fin >> u >> v >> temp){
		freq = temp;

		if(!freqCount.count({u,v}))
			freqCount[{u,v}] = 0;
		freqCount[{u,v}] += freq;

		outNeighbour[u].insert(v);
		outTotalFreq[u] += freq;

		inNeighbour[v].insert(u);
		inTotalFreq[v] += freq;

		vertices.insert(u);
		vertices.insert(v);
	}

	fout << "FINISHED\n";
}

// int Approach1::partition1(int l, int r, vector<int> &sortedVertices, bool sourceNodeGrouping=true){
// 	int tot = 0;
// 	for(int i=l;i<=r;i++){
// 		if(sourceNodeGrouping) tot += outNeighbour[sortedVertices[i]].size();
// 		else tot += inNeighbour[sortedVertices[i]].size();
// 	}
// 	int pivot = l;
// 	int minz = 2000000000;
// 	int cur = 0;
// 	for(int i=l;i<=r;i++){
// 		if(sourceNodeGrouping) cur += outNeighbour[sortedVertices[i]].size();
// 		else cur += inNeighbour[sortedVertices[i]].size();
//
// 		int curr = tot - cur;
// 		if(cur >= curr && cur-curr < minz) minz = cur-curr, pivot = i;
// 		else if(cur <= curr && curr-cur < minz) minz = curr-cur, pivot = i;
// 	}
// 	return pivot;
// }
int Approach1::partition1(int l, int r, vector<long long> &sortedVertices, bool sourceNodeGrouping=true){
	long long Ftot = 0;
	double Ctot = 0;

	for(int i=l;i<=r;i++){
		if(sourceNodeGrouping){
			if(!outTotalFreq.count(sortedVertices[i])) continue;
			Ftot += outTotalFreq[sortedVertices[i]];
			Ctot += divide(outNeighbour[sortedVertices[i]].size() * outNeighbour[sortedVertices[i]].size(),outTotalFreq[sortedVertices[i]]);
		}		
		else{
			if(!inTotalFreq.count(sortedVertices[i])) continue;
			Ftot += inTotalFreq[sortedVertices[i]];
			Ctot += divide(inNeighbour[sortedVertices[i]].size() * inNeighbour[sortedVertices[i]].size(),inTotalFreq[sortedVertices[i]]);
		}
	}

	long long F1 = 0;
	double C1 = 0;

	double minE = 1e15;
	int ret = l;

	for(int pivot = l;pivot < r; pivot++){

		if(sourceNodeGrouping && outTotalFreq.count(sortedVertices[pivot])){
			F1 += outTotalFreq[sortedVertices[pivot]];
			C1 += divide(outNeighbour[sortedVertices[pivot]].size() * outNeighbour[sortedVertices[pivot]].size(),outTotalFreq[sortedVertices[pivot]]);
		}
		else if(inTotalFreq.count(sortedVertices[pivot])){
			F1 += inTotalFreq[sortedVertices[pivot]];
			C1 += divide(inNeighbour[sortedVertices[pivot]].size() * inNeighbour[sortedVertices[pivot]].size(),inTotalFreq[sortedVertices[pivot]]);
		}

		if(F1 * C1 + (Ftot - F1) * (Ctot - C1) < minE){
			ret = pivot;
			minE = F1 * C1 + (Ftot - F1) * (Ctot - C1);
		}
	}

	return ret;
}

long long Approach1::getDistinctEdges(int l, int r){
	if(l == 0) return sumDistinctEdges[r];
	else return sumDistinctEdges[r] - sumDistinctEdges[l-1];
}

bool Approach1::terminate(int l, int r, int rows, int cols, bool sourceNodeGrouping=true){
	if(rows < w0 || cols < w0) return true;
	return divide(getDistinctEdges(l,r),rows*cols) <= C;
}

void Approach1::construct(int l, int r, vector<long long> &v, int rows, int cols){
	//
	// double totFreq = 0, totDeg = 0;
	// for(int i=l;i<=r;i++){
	// 	if(mode == 0) totFreq += outTotalFreq[v[i]], totDeg += outNeighbour[v[i]].size();
	// 	else totFreq += inTotalFreq[v[i]], totDeg += inNeighbour[v[i]].size();
	// }
	// cout << "PARTITION " << numberofgroups << ": " << divide(totFreq,totDeg) << '\n';

	for(int i=l;i<=r;i++) G[v[i]] = numberofgroups;
	Gmatrices[numberofgroups] = Gmatrix(rows,cols,K,P);
	numberofgroups++;
}

void Approach1::reset(){
	G.clear();
	numberofgroups = 0;
	Gmatrices.clear();
}

/*
Suppose you want to divide rows
Currently, there are S.rows rows

Suppose you want to partition such that:
	S1 gets (pivot - S.l + 1) / (S.r-S.l+1)
	S2 gets (S.r - (pivot+1) + 1) / (S.r-S.l+1)

S1.rows = ((pivot - S.l + 1) * S.rows) / (S.r-S.l+1) + ((((pivot - S.l + 1) * S.rows) % (S.r-S.l+1) != 0)?(1):(0))
S2.rows = ((S.r - (pivot+1) + 1) * S.rows) / (S.r-S.l+1)
*/

Approach1::Approach1(string data_sample_file, int rows, int cols, int outlier_rows, int outlier_cols, int depth, int modulo,int w0, double C):K(depth),P(modulo),w0(w0),C(C){
	cout << "ROWS,COLS: " << rows << " " << cols << '\n';
	if(rows!=0 && cols != 0){
		lookup(data_sample_file);
		fout << "BEGIN\n";
		reset();
		double outVar = 0, inVar = 0;
		for(auto it:freqCount){
			long long u = it.first.first, v = it.first.second;
			long long freq = it.second;
			outVar += ((double)freq - divide(outTotalFreq[u],outNeighbour[u].size()))
					* ((double)freq - divide(outTotalFreq[u],outNeighbour[u].size()));
			inVar += ((double)freq - divide(inTotalFreq[u],inNeighbour[u].size()))
					* ((double)freq - divide(inTotalFreq[u],inNeighbour[u].size()));
		}
		fout << "SORTING\n";
		vector<long long> sortedVertices(vertices.begin(),vertices.end());
		if(outVar <= inVar) sort(sortedVertices.begin(),sortedVertices.end(),outCmp(this));
		else sort(sortedVertices.begin(),sortedVertices.end(),inCmp(this));
		fout << "SORTED\n";
		
		sumDistinctEdges.resize(sortedVertices.size(),0);
		
		for(int i=0;i<sortedVertices.size();i++){
			if(i==0){
				if(outVar <= inVar) sumDistinctEdges[i] = outNeighbour[sortedVertices[i]].size();
				else sumDistinctEdges[i] = inNeighbour[sortedVertices[i]].size();
			}
			else{
				if(outVar <= inVar) sumDistinctEdges[i] = sumDistinctEdges[i-1] + outNeighbour[sortedVertices[i]].size();
				else sumDistinctEdges[i] = sumDistinctEdges[i-1] + inNeighbour[sortedVertices[i]].size();
			}
		}

		if(outVar <= inVar) mode = 0;
		else mode = 1;

		queue<sketch> q;
		q.push(sketch(rows,cols,0,sortedVertices.size()-1));

		while(!q.empty()){
			sketch S = q.front(); q.pop();
			if(S.l > S.r) continue;
			int pivot = partition1(S.l,S.r,sortedVertices,(outVar < inVar));
			fout << S.l << " " << S.r << " " << pivot << '\n';
			sketch S1,S2;

			if(outVar <= inVar){
				S1 = sketch(
					//S.rows/2,
					//((pivot - S.l + 1) * S.rows) / (S.r-S.l+1) + ((((pivot - S.l + 1) * S.rows) % (S.r-S.l+1) != 0)?(1):(0)),
					(int)((S.rows * getDistinctEdges(S.l,pivot))/getDistinctEdges(S.l,S.r)),
					S.cols,
					S.l,
					pivot
				);
				S2 = sketch(
					//1*S.rows/2 + (((1*S.rows) % 2 != 0)?(1):(0)),
					// ((S.r - (pivot+1) + 1) * S.rows) / (S.r-S.l+1),
					S.rows - (int)((S.rows * getDistinctEdges(S.l,pivot))/getDistinctEdges(S.l,S.r)),
					S.cols,
					pivot+1,
					S.r
				);
			}
			else{
				S1 = sketch(
					S.rows,
					(int)((S.cols * getDistinctEdges(S.l,pivot))/getDistinctEdges(S.l,S.r)),
					// ((pivot - S.l + 1) * S.cols) / (S.r-S.l+1) + ((((pivot - S.l + 1) * S.cols) % (S.r-S.l+1) != 0)?(1):(0)),
					S.l,
					pivot
				);
				S2 = sketch(
					S.rows,
					S.cols - (int)((S.cols * getDistinctEdges(S.l,pivot))/getDistinctEdges(S.l,S.r)),
					//1*S.cols/2 + (((1*S.cols) % 2 != 0)?(1):(0)),
					// ((S.r - (pivot+1) + 1) * S.cols) / (S.r-S.l+1),
					pivot+1,
					S.r
				);
			}

			if(!terminate(S1.l,S1.r,S1.rows,S1.cols,outVar<=inVar)) q.push(S1);
			else construct(S1.l,S1.r,sortedVertices,S1.rows,S1.cols);
			if(!terminate(S2.l,S2.r,S2.rows,S2.cols,outVar<=inVar)) q.push(S2);
			else construct(S2.l,S2.r,sortedVertices,S2.rows,S2.cols);
		}
		fout << "NUMBER OF PARTITIONS: " << numberofgroups << '\n';
		fout << G.size() << "\n";
		for(int i=0;i<=numberofgroups-1;i++){
			long long cnt = 0;
			for(auto it:G) if(it.second == i) cnt++;
			fout << i << " " << Gmatrices[i].rows << " " << Gmatrices[i].cols << ": " << cnt << "\n";
		}
	}
	outlierSketch = Gmatrix(outlier_rows,outlier_cols,K,P);
}

void Approach1::add(long long u, long long v, long long freq){
	if(mode == 0){
		if(G.count(u)){
			Gmatrices[G[u]].add(u,v,freq);
		}else{
			outlierSketch.add(u,v,freq);
		}
	}
	else{
		if(G.count(v)){
			Gmatrices[G[v]].add(u,v,freq);
		}
		else{
			outlierSketch.add(u,v,freq);
		}
	}
}

long long Approach1::query(long long u, long long v){
	if(mode == 0){
		if(G.count(u)) return Gmatrices[G[u]].query(u,v);
		else return outlierSketch.query(u,v);
	}
	else{
		if(G.count(v)) return Gmatrices[G[v]].query(u,v);
		else return outlierSketch.query(u,v);
	}
}

Approach1::Approach1(){}
