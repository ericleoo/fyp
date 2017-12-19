#include "Gsketch.h"

struct Gsketch::outCmp{
	Gsketch* m;
	outCmp(Gsketch* p) : m(p) {};

	bool operator() ( int i, int j )
	{
		if(!m->outNeighbour.count(i) && !m->outNeighbour.count(j)) return i < j;
		if(!m->outNeighbour.count(i)) return true;
		if(!m->outNeighbour.count(j)) return false;
		return m->outTotalFreq[i] * m->outNeighbour[j].size() < m->outTotalFreq[j] * m->outNeighbour[i].size();
	}
};

struct Gsketch::inCmp{
	Gsketch* m;
	inCmp(Gsketch* p) : m(p) {};

	bool operator() ( int i, int j )
	{
		if(!m->inNeighbour.count(i) && !m->inNeighbour.count(j)) return i < j;
		if(!m->inNeighbour.count(i)) return true;
		if(!m->inNeighbour.count(j)) return false;
		return m->inTotalFreq[i] * m->inNeighbour[j].size() < m->inTotalFreq[j] * m->inNeighbour[i].size();
	}
};

double Gsketch::divide(double one, double two){
	if(two < 0.0001) return 0;
	else return one/two;
}

void Gsketch::lookup(string file){
	/*
	Reads data sample. Gathers statistics
	*/
	ifstream fin(file);
	int u,v; long long freq;
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

int Gsketch::partition1(int l, int r, vector<int> &sortedVertices, bool sourceNodeGrouping=true){
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

bool Gsketch::terminate(int l, int r, vector<int> &v, int width, bool sourceNodeGrouping=true){
	if(width < w0) return true;
	int distinctedges = 0;
	for(int i=l;i<=r;i++){
		if(sourceNodeGrouping) distinctedges += outNeighbour[v[i]].size();
		else distinctedges += inNeighbour[v[i]].size();
	}
	if(divide(distinctedges,width) <= C) fout << "WOW\n";
	return divide(distinctedges,width) <= C;
}

void Gsketch::construct(int l, int r, vector<int> &v, int width){
	//
	// double totFreq = 0, totDeg = 0;
	// for(int i=l;i<=r;i++){
	// 	if(mode == 0) totFreq += outTotalFreq[v[i]], totDeg += outNeighbour[v[i]].size();
	// 	else totFreq += inTotalFreq[v[i]], totDeg += inNeighbour[v[i]].size();
	// }
	// cout << "PARTITION " << numberofgroups << ": " << divide(totFreq,totDeg) << '\n';

	for(int i=l;i<=r;i++) G[v[i]] = numberofgroups;
	Sketches[numberofgroups] = CountMin(width,K,P);
	numberofgroups++;
}

void Gsketch::reset(){
	G.clear();
	numberofgroups = 0;
	Sketches.clear();
}

Gsketch::Gsketch(string data_sample_file, int width, int outlier_width, int depth, int modulo,int w0, double C):K(depth),P(modulo),w0(w0),C(C){
	lookup(data_sample_file);
	fout << "BEGIN\n";
	reset();
	double outVar = 0, inVar = 0;
	for(auto it:freqCount){
		int u = it.first.first, v = it.first.second;
		long long freq = it.second;
		outVar += ((double)freq - divide(outTotalFreq[u],outNeighbour[u].size()))
				* ((double)freq - divide(outTotalFreq[u],outNeighbour[u].size()));
		inVar += ((double)freq - divide(inTotalFreq[u],inNeighbour[u].size()))
				* ((double)freq - divide(inTotalFreq[u],inNeighbour[u].size()));
	}
	fout << "SORTING\n";
	vector<int> sortedVertices(vertices.begin(),vertices.end());
	if(outVar <= inVar) sort(sortedVertices.begin(),sortedVertices.end(),outCmp(this));
	else sort(sortedVertices.begin(),sortedVertices.end(),inCmp(this));
	fout << "SORTED\n";

	if(outVar <= inVar) mode = 0;
	else mode = 1;

	queue<sketch> q;
	q.push(sketch(width,0,sortedVertices.size()-1));

	while(!q.empty()){
		sketch S = q.front(); q.pop();
		if(S.l > S.r) continue;
		int pivot = partition1(S.l,S.r,sortedVertices,(outVar < inVar));
		//fout << S.l << " " << S.r << " " << pivot << '\n';
		sketch S1,S2;

		S1 = sketch(
			S.width/2,
			S.l,
			pivot
		);
		S2 = sketch(
			S.width/2 + (S.width % 2),
			pivot+1,
			S.r
		);

		if(!terminate(S1.l,S1.r,sortedVertices,S1.width,outVar<=inVar)) q.push(S1);
		else construct(S1.l,S1.r,sortedVertices,S1.width);
		if(!terminate(S2.l,S2.r,sortedVertices,S2.width,outVar<=inVar)) q.push(S2);
		else construct(S2.l,S2.r,sortedVertices,S2.width);
	}
	outlierSketch = CountMin(outlier_width,K,P);
	//
	// fout << "NUMBER OF PARTITIONS: " << numberofgroups << '\n';
	// fout << G.size() << "\n";
	// for(int i=0;i<numberofgroups-1;i++){
	// 	int cnt = 0;
	// 	for(auto it:G) if(it.second == i) cnt++;
	// 	fout << i << " " << Sketches[i].rows << " " << Sketches[i].cols << ": " << cnt << "\n";
	// }
}

void Gsketch::add(int u, int v, long long freq){
	if(mode == 0){
		if(G.count(u)) Sketches[G[u]].add(u,v,freq);
		else outlierSketch.add(u,v,freq);
	}
	else{
		if(G.count(v)) Sketches[G[v]].add(u,v,freq);
		else outlierSketch.add(u,v,freq);
	}
}

long long Gsketch::query(int u, int v){
	if(mode == 0){
		if(G.count(u)) return Sketches[G[u]].query(u,v);
		else return outlierSketch.query(u,v);
	}
	else{
		if(G.count(v)) return Sketches[G[v]].query(u,v);
		else return outlierSketch.query(u,v);
	}
}

Gsketch::Gsketch(){}
