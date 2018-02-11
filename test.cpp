#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <algorithm>
#include <queue>
#include <cmath>
#include <chrono>
#include <string>
#include "HASH.h"
#include "Approach1.h"
#include "Approach2.h"
#include "CountMin.h"
#include "Gsketch.h"
using namespace std;

/*
#define GRAPH_STREAM_FILE "tweet_stream_hashed_refined"
#define DATA_SAMPLE_FILE "tweet_stream_hashed_reservoir5.txt"
#define QUERY_FILE "tweet_stream_hashed_reservoir2.txt"
#define P 17813333
#define total_freq 146039643LL
#define NUM_OF_LINES 78508963
*/

#define GRAPH_STREAM_FILE "graph_shuffled.txt"
//#define DATA_SAMPLE_FILE "graph_freq_comp1_reservoir5.txt"
#define DATA_SAMPLE_FILE "graph_freq_comp1_reservoir5_SORTED_0.171396_0_.txt"
#define QUERY_FILE "graph_freq_comp1_reservoir.txt"
#define P 56175521
#define total_freq 9812800185LL
#define NUM_OF_LINES 1372146644

/*
#define GRAPH_STREAM_FILE "ip_graph_refined"
#define DATA_SAMPLE_FILE "ip_graph_refined_reservoir5.txt"
#define QUERY_FILE "ip_graph_refined_reservoir2.txt"
#define P 4213103
#define total_freq 436186619LL
#define NUM_OF_LINES 12714850
*/

#define w0 400
#define C 0.1
#define N 1000
#define M 1000
#define K 10

//#define HH_CONSTANT 100
// #define HH_CONSTANT 1000
#define HH_CONSTANT 10000

#define APPROACH 1

#define ffout std::cout

void logging(string s){
	fstream logg = fstream("logg.txt",fstream::app);
	logg << s;
	logg.close();
}

void evaluate1x(Approach1 &app, Gmatrix &control)
{
	double are = 0, are2 = 0, are3 = 0, are4 = 0, avg = 0;
	long long one = 0, two = 0;
	long long num = 0;

	ifstream fin = ifstream(QUERY_FILE);
	long long u, v, freq;
	double temp;

	long long cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0;

	while (fin >> u >> v >> temp)
	{
		freq = temp;
		if ((app.mode == 0 && app.G.count(u)) || (app.mode == 1 && app.G.count(v)))
			one++;
		else
			two++;

		double e1 = (app.query(u, v) - freq) / freq;
		double e2 = (control.query(u, v) - freq) / freq;

		if (e1 <= 5)
			cnt1++;
		else
		{
			bool fnd = ((app.mode == 0 && app.G.count(u)) || (app.mode == 1 && app.G.count(v)));
			//cout << "Bad: " << fnd << " " << e1 << '\n';
			if (fnd)
				cnt3++;
			else
				cnt4++;
		}
		if (e2 <= 5)
			cnt2++;

		are += (app.query(u, v) - freq);
		are2 += (control.query(u, v) - freq);

		are3 += e1;
		are4 += e2;

		num++;

		avg += freq;
	}

	string z;

	z = "HIT: " + to_string(one / (double)(one + two)) + "\n"
		+ "MISS: " + to_string(two / (double)(one + two)) + "\n"
		+ "OBSERVED ERROR: " + to_string(are / avg) + "\n"
		+ "OBSERVED ERROR CONTROL: " + to_string(are2 / avg) + "\n"
		+ "AVERAGE RELATIVE ERROR: " + to_string(are3 / num) + "\n"
		+ "AVERAGE RELATIVE ERROR CONTROL: " + to_string(are4 / num) + "\n"
		+ "EFFECTIVE QUERIES: " + to_string(cnt1) + "/" + to_string(num) + " = " + to_string((cnt1 / (double)num) * 100.0) + "\n"
		+ "EFFECTIVE QUERIES CONTROL: " + to_string(cnt2) + "/" + to_string(num) + " = " + to_string((cnt2 / (double)num) * 100.0) + "\n"
		+ "AVERAGE: " + to_string(avg / num) + "\n"
		+ "BAD PARTITIONING: " + to_string(cnt3) + "/" + to_string(num) + " = " + to_string((cnt3 / (double)num) * 100.0) + "\n"
		+ "BAD OUTLIER: " + to_string(cnt4) + "/" + to_string(num) + " = " + to_string((cnt4 / (double)num) * 100.0) + "\n"
		+ "TOTAL BAD: " + to_string(cnt3 + cnt4) + "/" + to_string(num) + " = " + to_string(((cnt3 + cnt4) / (double)num) * 100.0) + "\n";
	cout << z; logging(z);
	
}

void evaluate2x(Approach2 &app, Gmatrix &control)
{
	int tc = 500;
	long long one = 0, two = 0, three = 0;
	double are = 0, are2 = 0, are3 = 0, are4 = 0;
	double tot = 0;
	long long minz = 1000000000000, maxz = 0, avg = 0;

	long long u, v, freq;
	double temp;
	ifstream fin = ifstream(QUERY_FILE);
	while (fin >> u >> v >> temp)
	{
		freq = temp;
		if (app.G1.count(u) && app.G2.count(v) && app.intersect(app.G1[u], app.G2[v]).size())
			three++;
		else if (app.G1.count(u) || app.G2.count(v))
			one++;
		else
			two++;

		//if(app.query(u,v) < freq) fout << app.query(u,v) << " " << freq << '\n';

		are += (app.query(u, v) - freq);
		are2 += (control.query(u, v) - freq);
		minz = min(minz, freq);
		maxz = max(maxz, freq);
		avg += freq;
		tot += freq;
	}
	ffout << "HIT2: " << three / (double)(one + two + three) << '\n';
	ffout << "HIT: " << one / (double)(one + two + three) << '\n';
	ffout << "MISS: " << two / (double)(one + two + three) << '\n';
	ffout << "OBSERVED ERROR: " << are / tot << '\n';
	ffout << "OBSERVED ERROR CONTROL: " << are2 / tot << '\n';
	ffout << "MIN: " << minz << '\n';
	ffout << "MAX: " << maxz << '\n';
	ffout << "AVG: " << avg / 1000000 << '\n';
}

void heavyHitter(int divisor, unordered_set<pair<int, int>,HASH> &heavy1, Approach1 &app)
{
	fstream logg;
	string z = "\nGetting HH for app\n";
	
	auto start = std::chrono::high_resolution_clock::now();	
	unordered_set<pair<int, int>,HASH> hh1 = app.getHeavyHitterEdges(total_freq / divisor);
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	
	z = "\nElapsed time: " + to_string(elapsed.count()) + " s\nDone.\nGetting false positives\n";
	cout << z; logging(z);	  
	
	int fp1 = 0;

	for (auto it : hh1)
	{
		if (!heavy1.count(it))
			fp1++;
	}
	
	
	z = "False positive rate " + to_string((1.0 / (double)divisor) * 100.0) + "%: " + to_string(fp1 / ((double)heavy1.size())) + "\n";
	cout << z; logging(z);
	
	for (auto it : heavy1)
	{
		if (!hh1.count(it))
		{
			cout << "WRONG\n";
			break;
		}
	}

	z = "Edge set sizes:\nActual: " + to_string(heavy1.size()) + "\ngMatrix with partitioning: " + to_string(hh1.size()) + "\n";

	cout << z; logging(z);

	hh1.clear();
	hh1.rehash(0);
}

void heavyHitterControl(int divisor, unordered_set<pair<int, int>,HASH> &heavy1, Gmatrix &control)
{
	fstream logg;
	
	string z = "\nGetting HH for control\n";
	cout << z; logging(z);
	
	auto start = std::chrono::high_resolution_clock::now();	
	unordered_set<pair<int, int>,HASH> hhc = control.getHeavyHitterEdges(total_freq / divisor);
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;

	z = "\nElapsed time: " + to_string(elapsed.count()) + " s\n";
	cout << z; logging(z);

	z = "Done.\nGetting false positives\n";
	cout << z; logging(z);
	
	int fpc = 0;
	for (auto it : hhc)
	{
		if (!heavy1.count(it))
			fpc++;
	}

	z = "False positive control " + to_string((1.0 / (double)divisor) * 100.0) + "%: " + to_string(fpc / ((double)heavy1.size())) + "\n";
	cout << z; logging(z);
	
	for (auto it : heavy1)
	{
		if (!hhc.count(it))
		{
			cout << "WRONG\n";
			break;
		}
	}

	z = "Edge set sizes:\nActual: " + to_string(heavy1.size()) + "\ngMatrix: " + to_string(hhc.size()) + "\n";
	cout << z; logging(z);
	
	hhc.clear();
	hhc.rehash(0);
}

unordered_map<int,long long> aggFreq;

void evalNodeAgg(Approach1 &app, Gmatrix &control){
	set<pair<long long,int>> vertices;

	int idx = 0;
	long long perc = 0;
	for(auto it:aggFreq){

		long long cur = (idx * (long long)100) / (long long) aggFreq.size();
		if(cur > perc){
			string z = to_string(cur) + "%\n";
			cout << z; logging(z);
			perc = cur;
		}

		vertices.insert({it.second,it.first});
		while((int)vertices.size() > 500)
			vertices.erase(vertices.begin());
		idx++;
	}

	double error1 = 0;
	double error2 = 0;
	double error3 = 0;
	double error4 = 0;
	double tot = 0;

	double cnt1 = 0, cnt2 = 0;

	for(auto it:vertices){
		int u = it.second;
		long long freq = it.first;

		long long one = app.queryNodeAggFreq(u);
		long long two = control.queryNodeOutgoingFreq(u);		

		if(one < freq) cout << "WEIRD1\n", logging("WEIRD1\n");
		if(two < freq) cout << "WEIRD2\n", logging("WEIRD2\n"); 

		error1 += one - freq;
		error2 += two - freq;

		double e1 = (double)(one - freq) / (double)freq;
		double e2 = (double)(two - freq) / (double)freq;

		if(e1 <= 5) cnt1+=1;
		if(e2 <= 5) cnt2+=1;

		error3 += e1;
		error4 += e2;

		tot += freq;
	}

	string z = "Observed error app: " + to_string(error1 / tot) 
				+ "\nObserved error control: " + to_string(error2/tot) + "\n" 
				+ "Average relative error app: " + to_string(error3/500) + "\n"
				+ "Average relative error control: " + to_string(error4/500) + "\n"
				+ "Effective queries: " + to_string(cnt1) + "/" + to_string(500) + " = " + to_string((cnt1 / (double)500) * 100.0) + "\n"
				+ "Effective queries control: " + to_string(cnt2) + "/" + to_string(500) + " = " + to_string((cnt2 / (double)500) * 100.0) + "\n";
	cout << z; logging(z);
}

int main()
{	
	if (APPROACH == 1)
	{
		Approach1 app = Approach1(string(DATA_SAMPLE_FILE), N, M, K, P, w0, C);
		ifstream fin(GRAPH_STREAM_FILE);

		Gmatrix control = Gmatrix(N, M, K, P);

		long long u, v;
		long long freq;
		double temp;

		//unordered_set<pair<int, int>,HASH> heavy1;		
		aggFreq.clear();

		long long perc = 0;
		
		string z = to_string(perc) + "%" +"\n";
		cout << z; logging(z);
		
		for (long long tc = 0; /*(tc < 12000000) &&*/ (fin >> u >> v >> temp); tc++)
		{

			long long cur = (tc * 100) / NUM_OF_LINES;
			if(cur > perc){
				z = to_string(cur) + "%" + "\n";
				cout << z; logging(z);
				perc = cur;
			}

			if((tc % 10000000 == 0)) evaluate1x(app,control);
			freq = temp;
			app.add(u, v, freq);
			control.add(u, v, freq);

			if (app.query(u, v) < freq)
			{
				ffout << app.query(u, v) << " " << freq << " " << u << " " << v << '\n';
			}
			//if (freq >= total_freq / HH_CONSTANT)
			//	heavy1.insert({u, v});

			if(aggFreq.count(u)) aggFreq[u] += freq;
			else aggFreq[u] = freq; 
		}

		evalNodeAgg(app,control);
		evaluate1x(app,control);
		//heavyHitter(HH_CONSTANT,heavy1,app);
		//heavyHitterControl(HH_CONSTANT,heavy1,control);
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
