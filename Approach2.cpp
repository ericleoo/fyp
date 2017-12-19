#include "Approach2.h"

struct Approach2::sketch{
	int rows,cols;
	vector<long long> vertices1;
	vector<long long> vertices2;

	set<pair<long long,pair<long long,long long>>> edgeset;

	sketch(){}

	sketch(int a, int b, vector<long long> v, vector<long long> w, set<pair<long long,pair<long long,long long>>> edges){
		rows = a; cols = b;
		vertices1 = v;
		vertices2 = w;
		edgeset = edges;
	}

	sketch(int a, int b, vector<long long> v, vector<long long> w, set<pair<long long,pair<long long,long long>>> edges, bool source){
		rows = a; cols = b;
		vertices1 = v;
		vertices2 = w;

		set<long long> vertexset1 = set<long long>(v.begin(),v.end());
		set<long long> vertexset2 = set<long long>(w.begin(),w.end());

		if(source){
			for(auto it:edges) if(vertexset1.count(it.second.first)){
				edgeset.insert(it);
			}
		}
		else{
			for(auto it:edges) if(vertexset2.count(it.second.second)){
				edgeset.insert(it);
			}
		}
	}
};

vector<int> Approach2::intersect(vector<int> s1, vector<int> s2){
	vector<int> v(s1.size() + s2.size());
	auto it = set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),v.begin());
	v.resize(it-v.begin());
	return v;
}

double Approach2::divide(double one, double two){
	if(two < 0.0001) return 0;
	else return one/two;
}

void Approach2::lookup(string file){
	ifstream fin(file);
	long long u,v; long long freq;
	double temp;

	fout << "LOOKINGUP\n";
	while(fin >> u >> v >> temp){
		freq = temp;
		if(!freqCount.count({u,v}))
			freqCount[{u,v}] = 0;
		freqCount[{u,v}] += freq;
		vertices.insert(u);
		vertices.insert(v);
	}
	fout << "FINISHED\n";
}

pair<double,int> Approach2::sourcePartitioning(sketch &s){
	map<long long,long long> cnt;
	for(auto it:s.edgeset){
		cnt[it.second.first] += it.first;
	}

	vector<pair<long long,long long>> temp(s.vertices1.size());
	for(int i=0;i<s.vertices1.size();i++) temp[i] = {cnt[s.vertices1[i]],s.vertices1[i]};
	sort(temp.begin(),temp.end());

	for(int i=0;i<s.vertices1.size();i++) s.vertices1[i] = temp[i].second;

	double totalFreq = 0;
	double curTotal = 0;

	for(auto it:s.vertices1){
		if(cnt.count(it)){
			curTotal += divide(1,cnt[it]);
			totalFreq += cnt[it];
		}
	}

	double freq = 0;
	double cur = 0;

	double minz = 1e18;
	int ret = 0;

	for(int pivot = 0; pivot < s.vertices1.size()-1; pivot++){
		if(cnt.count(s.vertices1[pivot])){
			freq += cnt[s.vertices1[pivot]];
			cur  += divide(1,cnt[s.vertices1[pivot]]);
		}

		double E1 = freq * cur, E2 = (totalFreq - freq) * (curTotal - cur);
		if(E1 + E2 < minz){
			minz = E1 + E2;
			ret = pivot;
		}
	}

	return {minz,ret};
}

pair<double,int> Approach2::destinationPartitioning(sketch &s){
	map<long long,long long> cnt;
	for(auto it:s.edgeset){
		cnt[it.second.second] += it.first;
	}

	vector<pair<long long,long long>> temp(s.vertices2.size());
	for(int i=0;i<s.vertices2.size();i++) temp[i] = {cnt[s.vertices2[i]],s.vertices2[i]};
	sort(temp.begin(),temp.end());

	for(int i=0;i<s.vertices2.size();i++) s.vertices2[i] = temp[i].second;

	double totalFreq = 0;
	double curTotal = 0;

	for(auto it:s.vertices2){
		if(cnt.count(it)){
			curTotal += divide(1,cnt[it]);
			totalFreq += cnt[it];
		}
	}

	double freq = 0;
	double cur = 0;

	double minz = 1e18;
	int ret = 0;


	for(int pivot = 0; pivot < s.vertices2.size()-1; pivot++){
		if(cnt.count(s.vertices2[pivot])){
			freq += cnt[s.vertices2[pivot]];
			cur  += divide(1,cnt[s.vertices2[pivot]]);
		}

		double E1 = freq * cur, E2 = (totalFreq - freq) * (curTotal - cur);

		//cout << E1 + E2 << " ";
		if(E1 + E2 < minz){
			minz = E1 + E2;
			ret = pivot;
		}
	}

	return {minz,ret};
}

bool Approach2::terminate(sketch &s){
	if(s.rows < w0 || s.cols < w0) return true;
	if(divide(s.edgeset.size(),min(s.rows,s.cols)) <= C) cout << "WOW\n";
	return divide(s.edgeset.size(),min(s.rows,s.cols)) <= C;
}

void Approach2::construct(sketch &s){
	cout << "PARTITION #" << numberofgroups << ": " << s.rows << " " << s.cols << " " << s.edgeset.size() << '\n';

	for(auto it:s.vertices1){
		G1[it].push_back(numberofgroups);
	}
	for(auto it:s.vertices2){
		G2[it].push_back(numberofgroups);
	}
	Gmatrices[numberofgroups] = Gmatrix(s.rows,s.cols,K,P);
	numberofgroups++;
}

Approach2::Approach2(string file, int and_rows, int and_cols, int or_rows, int or_cols, int outlier_rows, int outlier_cols, int depth, int modulo, int w0, double C):K(depth),P(modulo),w0(w0),C(C){

	/*
	allocate p<1 of gmatrix for (u,v) such that u and v have been seen before
	allocate q<1 of gmatrix for (u,v) such that either u or v have been seen before
	allocate (1-p-q) of gmatrix for outliers
	*/

	lookup(file);

	set<pair<long long,pair<long long,long long>>> edges;
	for(auto it:freqCount) edges.insert({it.second,it.first});

	numberofgroups = 0;

	queue<sketch> q;
	q.push(sketch(
		and_rows,
		and_cols,
		vector<long long>(vertices.begin(),vertices.end()),
		vector<long long>(vertices.begin(),vertices.end()),
		edges
	));

	while(!q.empty()){
		sketch S = q.front(); q.pop();
		pair<double,int> one = sourcePartitioning(S);
		pair<double,int> two = destinationPartitioning(S);

		if(one.first <= two.first){
			sketch S1(
				S.rows/2,
				S.cols,
				vector<long long>(S.vertices1.begin(),S.vertices1.begin() + one.second + 1),
				vector<long long>(S.vertices2),
				S.edgeset,
				true
			);
			sketch S2(
				S.rows/2+(S.rows%2),
				S.cols,
				vector<long long>(S.vertices1.begin() + one.second + 1,S.vertices1.end()),
				vector<long long>(S.vertices2),
				S.edgeset,
				true
			);

			if(!terminate(S1)) q.push(S1);
			else construct(S1);

			if(!terminate(S2)) q.push(S2);
			else construct(S2);
		}
		else{
			sketch S1(
				S.rows,
				S.cols/2,
				vector<long long>(S.vertices1),
				vector<long long>(S.vertices2.begin(),S.vertices2.begin() + two.second + 1),
				S.edgeset,
				false
			);
			sketch S2(
				S.rows,
				S.cols/2+(S.cols%2),
				vector<long long>(S.vertices1),
				vector<long long>(S.vertices2.begin() + two.second + 1,S.vertices2.end()),
				S.edgeset,
				false
			);

			if(!terminate(S1)) q.push(S1);
			else construct(S1);

			if(!terminate(S2)) q.push(S2);
			else construct(S2);
		}
	}
	app1 = Approach1(file,or_rows,or_cols,outlier_rows,outlier_cols,K,P,w0,C);
}

void Approach2::add(long long u, long long v, long long freq){
	vector<int> intersection = intersect(G1[u],G2[v]);
	if(intersection.size()){
		Gmatrices[intersection[0]].add(u,v,freq);
	}
	else{
		//cout << "WEW\n";
		app1.add(u,v,freq);
	}
}

long long Approach2::query(long long u, long long v){
	vector<int> intersection = intersect(G1[u],G2[v]);
	if(intersection.size()){
		return Gmatrices[intersection[0]].query(u,v);
	}
	else{
		//cout << "WEW\n";
		return app1.query(u,v);
	}
}
