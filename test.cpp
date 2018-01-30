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

/*
#define GRAPH_STREAM_FILE "tweet_stream_hashed_refined"
#define DATA_SAMPLE_FILE "tweet_stream_hashed_reservoir5.txt"
#define QUERY_FILE "tweet_stream_hashed_reservoir2.txt"
#define P 78508967
#define total_freq 146039643LL
*/

/*
#define GRAPH_STREAM_FILE "graph_shuffled.txt"
//#define DATA_SAMPLE_FILE "graph_freq_comp1_reservoir5.txt"
#define DATA_SAMPLE_FILE "graph_freq_comp1_reservoir5_SORTED_0.171396_0_.txt"
#define QUERY_FILE "graph_freq_comp1_reservoir.txt"
#define P 1372146661
#define total_freq 9812800185LL
*/


#define GRAPH_STREAM_FILE "ip_graph_refined"
#define DATA_SAMPLE_FILE "ip_graph_refined_reservoir5.txt"
#define QUERY_FILE "ip_graph_refined_reservoir2.txt"
#define P 12714851
#define total_freq 436186619LL

#define w0 100
#define C 0.1
#define N 1000
#define M 1000
#define K 10


#define APPROACH 1

#define ffout std::cout

void evaluate1x(Approach1 &app, Gmatrix &control){
	double are = 0,are2 = 0,are3=0, are4 = 0, avg = 0;
	long long one = 0, two = 0;
	long long num = 0;
    
	ifstream fin = ifstream(QUERY_FILE);
	long long u,v,freq;
	double temp;
    
    long long cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0;
    
	while(fin >> u >> v >> temp){
		freq = temp;
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

void evaluate2x(Approach2 &app, Gmatrix &control){
	int tc = 500;
	long long one = 0, two = 0, three = 0;
	double are = 0,are2 = 0,are3=0, are4 = 0;
	double tot = 0;
	long long minz = 1000000000000, maxz = 0, avg = 0;

	long long u,v,freq;
	double temp;
	ifstream fin = ifstream(QUERY_FILE);
	while(fin >> u >> v >> temp){
		freq = temp;
		if(app.G1.count(u) && app.G2.count(v) && app.intersect(app.G1[u],app.G2[v]).size()) three++;
		else if(app.G1.count(u) || app.G2.count(v)) one++;
		else two++;

		//if(app.query(u,v) < freq) fout << app.query(u,v) << " " << freq << '\n';

		are  += (app.query(u,v) - freq);
		are2 += (control.query(u,v) - freq);
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
	ffout << "MIN: " << minz << '\n';
	ffout << "MAX: " << maxz << '\n';
	ffout << "AVG: " << avg/1000000 << '\n';	
}

int main(){
	//ffout = ofstream("o2_graph_0.494.out",ofstream::app);
	//fout = ofstream("output.out");
	//freopen("output.out","w",stdout);	

	if(APPROACH == 1){
        Approach1 app = Approach1(string(DATA_SAMPLE_FILE),N,M,K,P,w0,C);
		ifstream fin(GRAPH_STREAM_FILE);
		
		Gmatrix control = Gmatrix(N,M,K,P);
		
		long long u,v; long long freq;
		double temp;

		set<pair<int,int>> heavy1, heavy2, heavy3;

		for(int tc=0;/*(tc < 12000000) &&*/ (fin >> u >> v >> temp);tc++){
			if(tc % 1000000 == 0) ffout << tc << '\n';
			if((tc % 5000000 == 0)) evaluate1x(app,control);
			freq = temp;
			app.add(u,v,freq);
			control.add(u,v,freq);
			if(app.query(u,v) < freq){
				ffout << app.query(u,v) << " " << freq << " " << u << " " << v << '\n';
			}

			if(freq >= 0.01 * total_freq){
				heavy1.insert({u,v});
			}

			if(freq >= 0.001 * total_freq){
				heavy2.insert({u,v});
			}

			if(freq >= 0.0001 * total_freq){
				heavy3.insert({u,v});
			}
		}

		evaluate1x(app,control);

		set<pair<int,int>> &hh1 = app.getHeavyHitterEdges(0.01 * total_freq);
		set<pair<int,int>> &hh2 = app.getHeavyHitterEdges(0.001 * total_freq);
		set<pair<int,int>> &hh3 = app.getHeavyHitterEdges(0.0001 * total_freq);

		int fp1 = 0, fp2 = 0, fp3 = 0;

		for(auto it:heavy1){
			if(!hh1.count(it)){
				cout << "WRONG\n";
				break;
			}
		}

		for(auto it:heavy2){
			if(!hh2.count(it)){
				cout << "WRONG\n";
				break;
			}
		}

		for(auto it:heavy3){
			if(!hh3.count(it)){
				cout << "WRONG\n";
				break;
			}
		}

		for(auto it:hh1){
			if(!heavy1.count(it)) fp1++;
		}

		for(auto it:hh2){
			if(!heavy2.count(it)) fp2++;
		}

		for(auto it:hh3){
			if(!heavy3.count(it)) fp3++;
		}
		
		cout << "False positive rate 1%: " << fp/((double)heavy1.size()) << '\n';
		cout << "False positive rate 0.1%: " << fp/((double)heavy2.size()) << '\n';
		cout << "False positive rate 0.01%: " << fp/((double)heavy3.size()) << '\n';
	}
    /*
	else{
		Approach2 app = Approach2();//string(DATA_SAMPLE_FILE),AND_ROWS,AND_COLS,OR_ROWS,OR_COLS,OUTLIER_ROWS,OUTLIER_COLS,K,P,w0,C);
		ifstream fin(GRAPH_STREAM_FILE);
		Gmatrix control = Gmatrix(N,M,K,P);
		long long u,v,freq;
		double temp;


		for(int tc=0;(fin >> u >> v >> temp);tc++){
			if(tc % 100000 == 0) ffout << tc << '\n';
			if((tc % 1000000 == 0)) evaluate2x(app,control);
			freq = temp;
			app.add(u,v,freq);
			control.add(u,v,freq);
			if(app.query(u,v) < freq){
				ffout << app.query(u,v) << " " << freq << " " << u << " " << v << '\n';
			}
		}

		evaluate2x(app,control);
	}
    */

}
