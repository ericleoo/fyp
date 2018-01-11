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

#define GRAPH_STREAM_FILE "tweet_stream_hashed_refined"
#define DATA_SAMPLE_FILE "tweet_stream_hashed_reservoir5.txt"
#define QUERY_FILE "tweet_stream_hashed_reservoir2.txt"
/*
#define GRAPH_STREAM_FILE "graph_shuffled.txt"
#define DATA_SAMPLE_FILE "graph_freq_comp1_reservoir5_SORTED_0.171396_0_.txt"
#define QUERY_FILE "graph_freq_comp1_reservoir.txt"
*/
/*
#define GRAPH_STREAM_FILE "ip_graph_refined"
#define DATA_SAMPLE_FILE "ip_graph_refined_reservoir5.txt"
#define QUERY_FILE "ip_graph_refined_reservoir2.txt"
*/

#define w0 300
#define C 0.1
#define N 4700
#define M 4700
#define K 9
#define P 1000000007
#define APPROACH 1

#define ffout std::cout

void evaluate1(set<pair<long long,pair<long long,long long>>> &pq, Approach1 &app, Gmatrix &control){
    double are = 0,are2 = 0,are3=0, are4 = 0, avg = 0;
	long long one = 0, two = 0;
	long long num = 0;
    
	ifstream fin = ifstream(QUERY_FILE);
	long long u,v,freq;
	double temp;
    
    long long cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0;
    
    for(auto it:pq){
        long long u,v,freq;
        freq = it.first;
        u = it.second.first;
        v = it.second.second;
    
        if((app.mode == 0 && app.G.count(u)) || (app.mode == 1 && app.G.count(v)))
			one++;
		else two++;
        
        double e1 = (app.query(u,v) - freq)/freq;
        double e2 = (control.query(u,v) - freq)/freq;
        
        if(e1 <= 5) cnt1++;
        else{
            bool fnd = ((app.mode == 0 && app.G.count(u)) || (app.mode == 1 && app.G.count(v)));
            //cout << "Bad: " << fnd << " " << e1 << '\n';
            if(fnd) cnt3++;
            else cnt4++;
        }
        if(e2 <= 5) cnt2++;
        
		are  += (app.query(u,v) - freq);
		are2 += (control.query(u,v) - freq);
        
        are3 += e1;
        are4 += e2;
        
        num++;
        
        avg += freq;
	}
	ffout << "HIT: " << one/(double)(one+two) << '\n';
	ffout << "MISS: " << two/(double)(one+two) << '\n';
	ffout << "OBSERVED ERROR: " << are/avg << '\n';
	ffout << "OBSERVED ERROR CONTROL: " << are2/avg << '\n';
    ffout << "AVERAGE RELATIVE ERROR: " << are3/num << '\n';
	ffout << "AVERAGE RELATIVE ERROR CONTROL: " << are4/num << '\n';
    ffout << "EFFECTIVE QUERIES: " << cnt1 << "/" << num << " = " << (cnt1/(double)num) * 100.0 << '\n';
    ffout << "EFFECTIVE QUERIES CONTROL: " << cnt2 << "/" << num << " = " << (cnt2/(double)num) * 100.0 << '\n';
    ffout << "AVERAGE: " << avg/num << '\n';
    
    ffout << "BAD PARTITIONING: " << cnt3 << "/" << num << " = " << (cnt3/(double)num) * 100.0 << '\n';
    ffout << "BAD OUTLIER: " << cnt4 << "/" << num << " = " << (cnt4/(double)num) * 100.0 << '\n';
    ffout << "TOTAL BAD: " << cnt3+cnt4 << "/" << num << " = " << ((cnt3+cnt4)/(double)num) * 100.0 << '\n';
}

void evaluate2(set<pair<long long,pair<long long,long long>>> &pq, Approach2 &app, Gmatrix &control){
	int tc = 500;
	long long one = 0, two = 0, three = 0;
	double are = 0,are2 = 0,are3=0, are4 = 0;
	double tot = 0;
	long long minz = 1000000000000000000LL, maxz = 0, avg = 0;

	// fin = ifstream(QUERY_FILE);
	// while(tc--){
	// 	fin >> u >> v >> temp;

	//ofstream fout("log.txt");

	for(auto it:pq){
		//auto it = pq.top(); pq.pop();
		long long u,v,freq;
		freq = it.first;
		u = it.second.first;
		v = it.second.second;
		
		//if(app.intersercti(app)) three++;
		if(app.G1.count(u) && app.G2.count(v) && app.intersect(app.G1[u],app.G2[v]).size()) three++;
		else if(app.G1.count(u) || app.G2.count(v)) one++;
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
	ffout << "HIT2: " << three/(double)(one+two+three) << '\n';
	ffout << "HIT: " << one/(double)(one+two+three) << '\n';
	ffout << "MISS: " << two/(double)(one+two+three) << '\n';
	ffout << "OBSERVED ERROR: " << are/tot << '\n';
	ffout << "OBSERVED ERROR CONTROL: " << are2/tot << '\n';
	// cout << "OBSERVED ERROR CONTROL2: " << are3/tot << '\n';
	// cout << "OBSERVED ERROR CONTROL3: " << are4/tot << '\n';
	ffout << "MIN: " << minz << '\n';
	ffout << "MAX: " << maxz << '\n';
	ffout << "AVG: " << avg/500 << '\n';	
}

int main(){
	//ffout = ofstream("x_graph_1_0.196.out",ofstream::app);
	//fout = ofstream("output.out");
	//freopen("output.out","w",stdout);	

	if(APPROACH == 1){
		Approach1 app = Approach1(string(DATA_SAMPLE_FILE),N,M,K,P,w0,C);
		ifstream fin(GRAPH_STREAM_FILE);
		
		Gmatrix control = Gmatrix(N,M,K,P);
		set<pair<long long,pair<long long,long long>>> pq;
		
		long long u,v; long long freq;
		double temp;


		for(int tc=0;/*(tc < 12000000) &&*/ (fin >> u >> v >> temp);tc++){
			if(tc % 1000000 == 0) ffout << tc << " " << pq.size() << '\n';
			if((tc % 5000000 == 0)) evaluate1(pq,app,control);
			
			freq = temp;
			pq.insert({freq,{u,v}});
			while(pq.size() > 500) pq.erase(pq.begin());

			app.add(u,v,freq);
			control.add(u,v,freq);
            
			if(app.query(u,v) < freq){
				ffout << "PROBLEM: " << app.query(u,v) << " " << freq << " " << u << " " << v << '\n';
			}
		}
		evaluate1(pq,app,control);
	}
	else{
		Approach2 app = Approach2();//string(DATA_SAMPLE_FILE),AND_ROWS,AND_COLS,OR_ROWS,OR_COLS,OUTLIER_ROWS,OUTLIER_COLS,K,P,w0,C);
		ifstream fin(GRAPH_STREAM_FILE);
		Gmatrix control = Gmatrix(N,M,K,P);
		// CountMin control2 = CountMin(N*M,K,P);

		// Gsketch control3 = Gsketch(string(DATA_SAMPLE_FILE),(int)((1.0-Outlier_Percentage)*N*M),(int)((Outlier_Percentage)*N*M),K,P,w0*w0,C);
		set<pair<long long,pair<long long,long long>>> pq;
		long long u,v,freq;
		double temp;


		for(int tc=0;/*(tc < 12000000) &&*/ (fin >> u >> v >> temp);tc++){
			if(tc % 100000 == 0) ffout << tc << " " << pq.size() << '\n';
			if((tc % 1000000 == 0)) evaluate2(pq,app,control);
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
				ffout << app.query(u,v) << " " << freq << " " << u << " " << v << '\n';
			}
		}

		evaluate2(pq,app,control);
	}

}
