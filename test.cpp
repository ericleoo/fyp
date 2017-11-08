#include<fstream>
#include<iostream>
#include<map>
#include<set>
#include<algorithm>
#include<queue>
#include<cmath>
#include "Approach1.h"
#include "Approach2.h"
#include "CountMin.h"
#include "Gsketch.h"
using namespace std;
//
#define GRAPH_STREAM_FILE "tweet_stream_hashed_refined"
#define DATA_SAMPLE_FILE "tweet_stream_hashed_refined_reservoir.txt"
#define QUERY_FILE "tweet_stream_hashed_refined_reservoir_2.txt"

// #define GRAPH_STREAM_FILE "graph_freq_comp1.txt"
// #define DATA_SAMPLE_FILE "graph_freq_comp1_reservoir.txt"
// #define QUERY_FILE "ip_graph_refined_reservoir_3M.txt"

// #define GRAPH_STREAM_FILE 

#define w0 400
#define C 0.09
#define N 3500
#define M 3500
#define K 50
#define P 1000000007
// #define Outlier_Percentage 0.462
// #define Outlier_Percentage 0.358
// #define Outlier_Percentage 0.158
// #define Outlier_Percentage 0.358
#define Outlier_Percentage 0.214
#define USE_OUTLIER_SKETCH (Outlier_Percentage > 1e-7)
#define ROWS ((USE_OUTLIER_SKETCH)?((int)(ceil(sqrt(1.0-Outlier_Percentage) * N))):(N))
#define COLS ((USE_OUTLIER_SKETCH)?((int)(ceil(sqrt(1.0-Outlier_Percentage) * M))):(M))
#define OUTLIER_ROWS ((USE_OUTLIER_SKETCH)?(((int)sqrt(N*M-ROWS*COLS))):(0))
#define OUTLIER_COLS ((USE_OUTLIER_SKETCH)?(((int)sqrt(N*M-ROWS*COLS))):(0))

// #define AND_Percentage 0.10
#define AND_Percentage 0.05
#define OR_ROWS ((int)ceil(sqrt(1.0-AND_Percentage) * ROWS))
#define OR_COLS ((int)ceil(sqrt(1.0-AND_Percentage) * COLS))
#define AND_ROWS ((int)sqrt(ROWS*COLS-OR_ROWS*OR_COLS))
#define AND_COLS ((int)sqrt(ROWS*COLS-OR_ROWS*OR_COLS))

#define APPROACH 2

int main(){
	//fout = ofstream("output.out");

	if(APPROACH == 1){
		Approach1 app = Approach1(string(DATA_SAMPLE_FILE),ROWS,COLS,OUTLIER_ROWS,OUTLIER_COLS,K,P,w0,C);
		ifstream fin(GRAPH_STREAM_FILE);
		int one = 0, two = 0;
		double are = 0,are2 = 0,are3=0, are4 = 0;
		double tot = 0;
		Gmatrix control = Gmatrix(N,M,K,P);
		// CountMin control2 = CountMin(N*M,K,P);

		// Gsketch control3 = Gsketch(string(DATA_SAMPLE_FILE),(int)((1.0-Outlier_Percentage)*N*M),(int)((Outlier_Percentage)*N*M),K,P,w0*w0,C);
		set<pair<long long,pair<long long,long long>>> pq;
		//priority_queue<pair<long long,pair<long long,long long>>> pq;

		long long u,v; long long freq;
		double temp;


		for(int tc=0;/*(tc < 12000000) &&*/ (fin >> u >> v >> temp);tc++){
			if(tc % 100000 == 0) cout << tc << " " << pq.size() << '\n';
			freq = temp;
			pq.insert({freq,{u,v}});
			while(pq.size() > 500) pq.erase(pq.begin());
			
			
			// cout << "F1\n";
			app.add(u,v,freq);
			// cout << "F2\n";
			control.add(u,v,freq);
			// cout << "F3\n";
			// control2.add(u,v,freq);
			//
			// control3.add(u,v,freq);
			// cout << "F4\n";
			if(app.query(u,v) < freq){
				cout << app.query(u,v) << " " << freq << " " << u << " " << v << '\n';
			}
		}

		int tc = 500;
		long long minz = 1000000000000, maxz = 0, avg = 0;

		// fin = ifstream(QUERY_FILE);
		// while(tc--){
		// 	fin >> u >> v >> temp;

		for(auto it:pq){
			//auto it = *pq.begin(); pq.erase(pq.begin());
			long long u,v,freq;
			freq = it.first;
			u = it.second.first;
			v = it.second.second;

			if((app.mode == 0 && app.G.count(u)) || (app.mode == 1 && app.G.count(v)))
				one++;
			else two++;
			are  += (app.query(u,v) - freq);
			are2 += (control.query(u,v) - freq);
			// are3 += (control2.query(u,v) - temp);
			// are4 += (control3.query(u,v) - temp);
			minz = min(minz,freq);
			maxz = max(maxz,freq);
			avg  += freq;
			tot += freq;
		}
		cout << "HIT: " << one/(double)(one+two) << '\n';
		cout << "MISS: " << two/(double)(one+two) << '\n';
		cout << "OBSERVED ERROR: " << are/tot << '\n';
		cout << "OBSERVED ERROR CONTROL: " << are2/tot << '\n';
		// cout << "OBSERVED ERROR CONTROL2: " << are3/tot << '\n';
		// cout << "OBSERVED ERROR CONTROL3: " << are4/tot << '\n';
		cout << "MIN: " << minz << '\n';
		cout << "MAX: " << maxz << '\n';
		cout << "AVG: " << avg/500 << '\n';
	}
	else{
		Approach2 app = Approach2(string(DATA_SAMPLE_FILE),AND_ROWS,AND_COLS,OR_ROWS,OR_COLS,OUTLIER_ROWS,OUTLIER_COLS,K,P,w0,C);
		ifstream fin(GRAPH_STREAM_FILE);
		int one = 0, two = 0, three = 0;
		double are = 0,are2 = 0,are3=0, are4 = 0;
		double tot = 0;
		Gmatrix control = Gmatrix(N,M,K,P);
		// CountMin control2 = CountMin(N*M,K,P);

		// Gsketch control3 = Gsketch(string(DATA_SAMPLE_FILE),(int)((1.0-Outlier_Percentage)*N*M),(int)((Outlier_Percentage)*N*M),K,P,w0*w0,C);
		// priority_queue<pair<long long,pair<long long,long long>>> pq;
		set<pair<long long,pair<long long,long long>>> pq;
		long long u,v,freq;
		double temp;


		for(int tc=0;/*(tc < 12000000) &&*/ (fin >> u >> v >> temp);tc++){
			if(tc % 100000 == 0) cout << tc << " " << pq.size() << '\n';
			freq = temp;
			pq.insert({freq,{u,v}});
			while(pq.size() > 500) pq.erase(pq.begin());
			// cout << "F1\n";

			app.add(u,v,freq);
			// cout << "F2\n";
			control.add(u,v,freq);
			// cout << "F3\n";
			// control2.add(u,v,freq);
			//
			// control3.add(u,v,freq);
			// cout << "F4\n";

			if(app.query(u,v) < freq){
				cout << app.query(u,v) << " " << freq << " " << u << " " << v << '\n';
			}

		}

		int tc = 500;
		long long minz = 1000000000000, maxz = 0, avg = 0;

		// fin = ifstream(QUERY_FILE);
		// while(tc--){
		// 	fin >> u >> v >> temp;

		ofstream fout("log.txt");

		for(auto it:pq){
			//auto it = pq.top(); pq.pop();
			long long u,v,freq;
			freq = it.first;
			u = it.second.first;
			v = it.second.second;

			if(app.G1.count(u) && app.G2.count(v)) three++;
			else if(app.G1.count(u) || app.G2.count(v))
				one++;
			else two++;

			if(app.query(u,v) < freq) fout << app.query(u,v) << " " << freq << '\n';

			are  += (app.query(u,v) - freq);
			are2 += (control.query(u,v) - freq);
			// are3 += (control2.query(u,v) - temp);
			// are4 += (control3.query(u,v) - temp);
			minz = min(minz,freq);
			maxz = max(maxz,freq);
			avg  += freq;
			tot += freq;
		}
		cout << "HIT2: " << three/(double)(one+two+three) << '\n';
		cout << "HIT: " << one/(double)(one+two+three) << '\n';
		cout << "MISS: " << two/(double)(one+two+three) << '\n';
		cout << "OBSERVED ERROR: " << are/tot << '\n';
		cout << "OBSERVED ERROR CONTROL: " << are2/tot << '\n';
		// cout << "OBSERVED ERROR CONTROL2: " << are3/tot << '\n';
		// cout << "OBSERVED ERROR CONTROL3: " << are4/tot << '\n';
		cout << "MIN: " << minz << '\n';
		cout << "MAX: " << maxz << '\n';
		cout << "AVG: " << avg/500 << '\n';
	}

}
