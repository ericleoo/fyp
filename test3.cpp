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
//#define GRAPH_STREAM_FILE "tweet_stream_hashed_refined"
//#define DATA_SAMPLE_FILE "tweet_stream_hashed_refined_reservoir.txt"
//#define QUERY_FILE "tweet_stream_hashed_refined_reservoir_2.txt"

#define GRAPH_STREAM_FILE "graph_freq_comp1.txt"
#define DATA_SAMPLE_FILE "graph_freq_comp1_reservoir2.txt"
#define QUERY_FILE "graph_freq_comp1_reservoir.txt"

//#define GRAPH_STREAM_FILE "ip_graph_refined"
//#define DATA_SAMPLE_FILE "ip_graph_refined_reservoir.txt"
//#define QUERY_FILE "ip_graph_refined_reservoir2.txt"

#define w0 300
#define C 0.3
#define N 2000
#define M 2000
#define K 50
#define P 1000000007
// #define Outlier_Percentage 0.462
//#define Outlier_Percentage 0.05
#define Outlier_Percentage 0.196
//#define Outlier_Percentage 0.158
// #define Outlier_Percentage 0.358
//#define Outlier_Percentage 0.214
//#define Outlier_Percentage 0.05
#define USE_OUTLIER_SKETCH (Outlier_Percentage > 1e-7)
#define ROWS ((USE_OUTLIER_SKETCH)?((int)(ceil(sqrt(1.0-Outlier_Percentage) * N))):(N))
#define COLS ((USE_OUTLIER_SKETCH)?((int)(ceil(sqrt(1.0-Outlier_Percentage) * M))):(M))
#define OUTLIER_ROWS ((USE_OUTLIER_SKETCH)?(((int)sqrt(N*M-ROWS*COLS))):(0))
#define OUTLIER_COLS ((USE_OUTLIER_SKETCH)?(((int)sqrt(N*M-ROWS*COLS))):(0))

// #define AND_Percentage 0.10
#define AND_Percentage 0.362
#define OR_ROWS ((int)ceil(sqrt(1.0-AND_Percentage) * ROWS))
#define OR_COLS ((int)ceil(sqrt(1.0-AND_Percentage) * COLS))
#define AND_ROWS ((int)sqrt(ROWS*COLS-OR_ROWS*OR_COLS))
#define AND_COLS ((int)sqrt(ROWS*COLS-OR_ROWS*OR_COLS))

#define APPROACH 1

//#define fout cout
#define ffout std::cout
// ofstream fout;
//ofstream ffout;

void evaluate1(set<pair<long long,pair<long long,long long>>> &pq, Gsketch &app, CountMin &control){
	double are = 0,are2 = 0,are3=0, are4 = 0;
	long long one = 0, two = 0;
	double tot = 0;
	long long minz = 1000000000000000000LL, maxz = 0, avg = 0;

	// fin = ifstream(QUERY_FILE);
	// while(tc--){
	// 	fin >> u >> v >> temp;

	for(auto it:pq){
		//auto it = *pq.begin(); pq.erase(pq.begin());
		long long u,v,freq;
		freq = it.first;
		u = it.second.first;
		v = it.second.second;
		are  += (app.query(u,v) - freq);
		are2 += (control.query(u,v) - freq);
		minz = min(minz,freq);
		maxz = max(maxz,freq);
		avg  += freq;
		tot += freq;
	}
	ffout << "OBSERVED ERROR: " << are/tot << '\n';
	ffout << "OBSERVED ERROR CONTROL: " << are2/tot << '\n';
	ffout << "MIN: " << minz << '\n';
	ffout << "MAX: " << maxz << '\n';
	ffout << "AVG: " << avg/500 << '\n';	
}

int main(){
	//ffout = ofstream("x_graph_1_0.196.out",ofstream::app);
	
	//Approach1 app = Approach1(string(DATA_SAMPLE_FILE),ROWS,COLS,OUTLIER_ROWS,OUTLIER_COLS,K,P,w0,C);
	ifstream fin(GRAPH_STREAM_FILE);
	
	CountMin control2 = CountMin(N*M,K,P);
	Gsketch control3 = Gsketch(string(DATA_SAMPLE_FILE),N*M-(int)((Outlier_Percentage)*N*M),(int)((Outlier_Percentage)*N*M),K,P,w0*w0,C);
	
	set<pair<long long,pair<long long,long long>>> pq;
	
	long long u,v; long long freq;
	double temp;

	for(int tc=0;/*(tc < 12000000) &&*/ (fin >> u >> v >> temp);tc++){
		if(tc % 100000 == 0) ffout << tc << " " << pq.size() << '\n';
		if((tc % 1000000 == 0)) evaluate1(pq,control3,control2);
		
		freq = temp;
		pq.insert({freq,{u,v}});
		while(pq.size() > 500) pq.erase(pq.begin());
		
		// cout << "F1\n";
		// app.add(u,v,freq);
		// cout << "F2\n";
		// control.add(u,v,freq);
		// cout << "F3\n";
		control2.add(u,v,freq);
		//
		control3.add(u,v,freq);
		// cout << "F4\n";
		if(control3.query(u,v) < freq){
			ffout << "PROBLEM: " << control3.query(u,v) << " " << freq << " " << u << " " << v << '\n';
		}
	}

	evaluate1(pq,control3,control2);
}
