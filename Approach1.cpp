#include "Approach1.h"

struct Approach1::outCmp
{
	Approach1 *m;
	outCmp() {}
	outCmp(Approach1 *p) : m(p){};

	bool operator()(long long i, long long j)
	{
		if (!m->outNeighbour.count(i) && !m->outNeighbour.count(j))
			return i < j;
		if (!m->outNeighbour.count(i))
			return true;
		if (!m->outNeighbour.count(j))
			return false;
		return m->outTotalFreq[i] * m->outNeighbour[j].size() < m->outTotalFreq[j] * m->outNeighbour[i].size();
	}
};

struct Approach1::inCmp
{
	Approach1 *m;
	inCmp() {}
	inCmp(Approach1 *p) : m(p){};

	bool operator()(long long i, long long j)
	{
		if (!m->inNeighbour.count(i) && !m->inNeighbour.count(j))
			return i < j;
		if (!m->inNeighbour.count(i))
			return true;
		if (!m->inNeighbour.count(j))
			return false;
		return m->inTotalFreq[i] * m->inNeighbour[j].size() < m->inTotalFreq[j] * m->inNeighbour[i].size();
	}
};

vector<int> Approach1::intersect(set<int> s1, set<int> s2)
{
	vector<int> v(s1.size() + s2.size());
	auto it = set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), v.begin());
	v.resize(it - v.begin());
	return v;
}

double Approach1::divide(double one, double two)
{
	if (two < 1e-9)
		return 0;
	else
		return one / two;
}

double Approach1::getVars(double &outVar, double &inVar)
{
	outVar = 0;
	inVar = 0;

	double one = 0, two = 0;

	for (auto it : outNeighbour)
	{
		long long u = it.first;
		double avg = divide(outTotalFreq[u], outNeighbour[u].size());

		double var = 0;
		for (auto it2 : outNeighbour[u])
		{
			var += (((double)it2.second) - avg) * (((double)it2.second) - avg);
		}

		if (var / outNeighbour[u].size() > 100)
		{
			vertices1.erase(u);
		}
		else
		{
			outVar += var;
			one += outNeighbour[u].size();
		}
	}

	outVar /= one;

	/*
	for (auto it : inNeighbour)
	{
		long long u = it.first;
		double avg = divide(inTotalFreq[u], inNeighbour[u].size());

		double var = 0;
		for (auto it2 : inNeighbour[u])
		{
			var += (((double)it2.second) - avg) * (((double)it2.second) - avg);
		}

		if (var / inNeighbour[u].size() > 100)
			vertices2.erase(u);
		else
		{
			inVar += var;
			two += inNeighbour[u].size();
		}
	}

	inVar /= two;
	*/
	////////////////////////////////////////////////////////
}

void Approach1::lookup(string file)
{
	ifstream fin(file);

	long long u, v, freq;
	double temp;

	fout << "LOOKINGUP\n";

	int num_lines = 0;
	while (fin >> u >> v >> temp)
		num_lines++;

	fin = ifstream(file);

	sortedVertices.resize(num_lines);
	int idx = 0;

	while (fin >> u >> v >> temp)
	{
		freq = temp;

		if (outNeighbour.count(u))
		{
			auto it = outNeighbour[u].lower_bound({v, 0});
			if (it != outNeighbour[u].end() && it->first == v)
			{
				auto cur = *it;
				outNeighbour[u].erase(it);
				cur = {v, cur.second + freq};
				outNeighbour[u].insert(cur);
			}
			else
				outNeighbour[u].insert({v, freq});
		}
		else
			outNeighbour[u].insert({v, freq});

		outTotalFreq[u] += freq;
		/*
		/////////////////////////////////////////////////////////
		if (inNeighbour.count(v))
		{
			auto it = inNeighbour[v].lower_bound({u, 0});
			if (it != inNeighbour[v].end() && it->first == u)
			{
				auto cur = *it;
				inNeighbour[v].erase(it);
				cur = {u, cur.second + freq};
				inNeighbour[v].insert(cur);
			}
			else
				inNeighbour[v].insert({u, freq});
		}
		else
			inNeighbour[v].insert({u, freq});

		inTotalFreq[v] += freq;
		*/
		////////////////////////////////////////////////////////

		vertices1.insert(u);
		//vertices2.insert(v); ///////////////////////////////////

		if (file.find("SORTED") != std::string::npos)
		{
			mode = getMode(file);
			auto &x = ((mode == 0) ? (u) : (v));

			if (idx == 0)
				sortedVertices[idx++] = x;
			else if (sortedVertices[idx - 1] != x)
				sortedVertices[idx++] = x;
		}
	}

	double outVar = 0, inVar = 0;
	getVars(outVar, inVar);

	fin = ifstream(file);

	//=============================================================
	// Calculate outlier percentage
	//=============================================================
	if (file.find("SORTED") == std::string::npos)
	{
		unordered_set<long long> A_1, A_2;
		int cnt1 = 0, cnt2 = 0;
		int it = 0;

		while (fin >> u >> v >> temp)
		{
			freq = temp;
			if (it >= (9 * num_lines) / 10)
			{

				if (!A_1.count(u))
					cnt1++;
				/*
				///////////////////////
				if (!A_2.count(v))
					cnt2++;
				///////////////////////
				*/
			}
			else
			{

				if (vertices1.count(u))
					A_1.insert(u);
				/*
				////////////////////////
				if (vertices2.count(v))
					A_2.insert(v);
				////////////////////////
				*/
			}
			it++;
		}

		outlierPercentage1 = cnt1 / (double)(num_lines - 1 - ((9 * num_lines) / 10) + 1);
		outlierPercentage2 = cnt2 / (double)(num_lines - 1 - ((9 * num_lines) / 10) + 1); /////////////////////////

		A_1.clear();
		A_2.clear();
		A_1.rehash(0);
		A_2.rehash(0);
	}

	//mode = ((outVar > inVar && outlierPercentage2 < outlierPercentage1) ? (1) : (0));
	mode = 0;
	fout << "FINISHED\n";
}

int Approach1::partition1(int l, int r, vector<long long> &sortedVertices)
{

	long long Ftot = 0;
	double Ctot = 0;

	for (int i = l; i <= r; i++)
	{

		auto &Adj = ((mode == 0) ? (outNeighbour) : (inNeighbour));
		auto &totFreq = ((mode == 0) ? (outTotalFreq) : (inTotalFreq));

		if (!totFreq.count(sortedVertices[i]))
			continue;

		Ftot += totFreq[sortedVertices[i]];

		if (Adj[sortedVertices[i]].size() == 0)
		{
			cout << "REEEE\n";
		}
		if (totFreq[sortedVertices[i]] == 0)
		{
			cout << "???\n";
		}

		Ctot += divide(Adj[sortedVertices[i]].size() * Adj[sortedVertices[i]].size(), totFreq[sortedVertices[i]]);
	}

	long long F1 = 0;
	double C1 = 0;

	double minE = std::numeric_limits<double>::max();
	int ret = l;

	for (int pivot = l; pivot < r; pivot++)
	{

		auto &Adj = ((mode == 0) ? (outNeighbour) : (inNeighbour));
		auto &totFreq = ((mode == 0) ? (outTotalFreq) : (inTotalFreq));

		if (totFreq.count(sortedVertices[pivot]))
		{
			F1 += totFreq[sortedVertices[pivot]];
			C1 += divide(Adj[sortedVertices[pivot]].size() * Adj[sortedVertices[pivot]].size(), totFreq[sortedVertices[pivot]]);
		}

		if (pivot == l || pivot == r - 1)
			cout << pivot << ": " << F1 * C1 + (Ftot - F1) * (Ctot - C1) << '\n';

		if (F1 * C1 + (Ftot - F1) * (Ctot - C1) <= minE)
		{
			ret = pivot;
			minE = F1 * C1 + (Ftot - F1) * (Ctot - C1);
		}
	}

	cout << "ret: " << ret << " " << minE << '\n';

	return ret;
}

long long Approach1::getDistinctEdges(int l, int r)
{
	if (l == 0)
		return sumDistinctEdges[r];
	else
		return sumDistinctEdges[r] - sumDistinctEdges[l - 1];
}

bool Approach1::terminate(int l, int r, int rows, int cols)
{
	if (rows < w0 || cols < w0)
		return true;
	return divide(getDistinctEdges(l, r), rows * cols) <= C;
}

void Approach1::construct(int l, int r, vector<long long> &v, int rows, int cols)
{
	for (int i = l; i <= r; i++)
		G[v[i]] = numberofgroups;
	Gmatrices[numberofgroups] = Gmatrix(rows, cols, K, P);
	numberofgroups++;
}

void Approach1::reset()
{
	G.clear();
	numberofgroups = 0;
	Gmatrices.clear();
}

bool Approach1::sorting(string &s)
{
	// returns true iff file hasn't been sorted yet

	if (s.find("SORTED") != std::string::npos)
	{
		return false;
	}

	s = s.substr(0, s.size() - 4);
	s += "_SORTED_";
	if (mode == 0)
		s += to_string(outlierPercentage1);
	else
		s += to_string(outlierPercentage2);
	s += '_';
	s += mode + '0';
	s += "_.txt";

	auto &V = ((mode == 0) ? (vertices1) : (vertices2));

	sortedVertices = vector<long long>(V.begin(), V.end());
	fout << "SORTING\n";
	if (mode == 0)
		sort(sortedVertices.begin(), sortedVertices.end(), outCmp(this));
	else
		sort(sortedVertices.begin(), sortedVertices.end(), inCmp(this));
	fout << "SORTED\n";
	ofstream sOut(s);

	fout << "CREATING SORTED FILE\n";

	for (int i = 0; i < sortedVertices.size(); i++)
	{
		if (mode == 0)
		{
			long long u = sortedVertices[i];
			for (auto it : outNeighbour[u])
			{
				long long v = it.first, freq = it.second;
				sOut << u << " " << v << " " << freq << '\n';
			}
		}
		else
		{
			long long v = sortedVertices[i];
			for (auto it : inNeighbour[v])
			{
				long long u = it.first, freq = it.second;
				sOut << u << " " << v << " " << freq << '\n';
			}
		}
	}

	sOut.close();
	return true;
}

double Approach1::getPercentage(string s)
{
	auto idx = s.find("SORTED");
	while (s[idx] != '_')
		idx++;
	idx++;
	int j = idx;
	while (s[j] != '_')
		j++;
	return stod(s.substr(idx, j - 1 - idx + 1));
}

int Approach1::getMode(string s)
{
	auto idx = s.find("SORTED");
	while (s[idx] != '_')
		idx++;
	idx++;
	int j = idx;
	while (s[j] != '_')
		j++;
	j++;
	return s[j] - '0';
}

void Approach1::clearAll()
{
	outNeighbour.clear();
	inNeighbour.clear();
	inTotalFreq.clear();
	outTotalFreq.clear();
	vertices1.clear();
	vertices2.clear();
	sortedVertices.clear();

	outNeighbour.rehash(0);
	inNeighbour.rehash(0);
	inTotalFreq.rehash(0);
	outTotalFreq.rehash(0);
	vertices1.rehash(0);
	vertices2.rehash(0);

}

void Approach1::setup(string data_sample_file, int rows, int cols, int depth, int modulo, int w0, double C)
{
	int N = rows, M = cols;

	clearAll();

	lookup(data_sample_file);

	reset();

	int outlier_rows, outlier_cols;

	if (sorting(data_sample_file))
	{
		if (!mode)
		{
			cout << "Mode 0.\n";
			cout << "outlier: " << outlierPercentage1 << '\n';
			rows = (int)(ceil(sqrt(1.0 - outlierPercentage1) * N));
			cols = (int)(ceil(sqrt(1.0 - outlierPercentage1) * M));
		}
		else
		{
			cout << "Mode 1.\n";
			cout << "outlier: " << outlierPercentage2 << '\n';
			rows = (int)(ceil(sqrt(1.0 - outlierPercentage2) * N));
			cols = (int)(ceil(sqrt(1.0 - outlierPercentage2) * M));
		}
		outlier_rows = (int)sqrt(N * M - rows * cols);
		outlier_cols = (int)sqrt(N * M - rows * cols);

		cout << "ROWS,COLS: " << rows << " " << cols << '\n';

		clearAll();
		setup(data_sample_file, N, M, depth, modulo, w0, C);
		return;
	}
	else
	{
		mode = getMode(data_sample_file);
		if (!mode)
		{
			outlierPercentage1 = getPercentage(data_sample_file);
			cout << "Mode 0.\n";
			cout << "outlier: " << outlierPercentage1 << '\n';
			rows = (int)(ceil(sqrt(1.0 - outlierPercentage1) * N));
			cols = (int)(ceil(sqrt(1.0 - outlierPercentage1) * M));
		}
		else
		{
			outlierPercentage2 = getPercentage(data_sample_file);
			cout << "Mode 1.\n";
			cout << "outlier: " << outlierPercentage2 << '\n';
			rows = (int)(ceil(sqrt(1.0 - outlierPercentage2) * N));
			cols = (int)(ceil(sqrt(1.0 - outlierPercentage2) * M));
		}
		outlier_rows = (int)sqrt(N * M - rows * cols);
		outlier_cols = (int)sqrt(N * M - rows * cols);

		cout << "ROWS,COLS: " << rows << " " << cols << '\n';
	}

	if (rows != 0 && cols != 0)
	{

		sumDistinctEdges.resize(sortedVertices.size(), 0);

		for (int i = 0; i < sortedVertices.size(); i++)
		{
			auto &Adj = ((mode == 0) ? (outNeighbour) : (inNeighbour));
			if (i == 0)
				sumDistinctEdges[i] = Adj[sortedVertices[i]].size();
			else
				sumDistinctEdges[i] = sumDistinctEdges[i - 1] + Adj[sortedVertices[i]].size();
		}

		queue<sketch> q;
		q.push(sketch(rows, cols, 0, sortedVertices.size() - 1));

		while (!q.empty())
		{
			sketch S = q.front();
			q.pop();
			if (S.l > S.r)
				continue;

			if (terminate(S.l, S.r, S.rows, S.cols))
			{
				construct(S.l, S.r, sortedVertices, S.rows, S.cols);
				continue;
			}
			int pivot = partition1(S.l, S.r, sortedVertices);
			fout << S.l << " " << S.r << " " << pivot << '\n';
			sketch S1, S2;

			if (!mode)
			{
				int one = S.rows / 2;
				int two = S.rows - one;

				S1 = sketch(
					one,
					S.cols,
					S.l,
					pivot);
				S2 = sketch(
					two,
					S.cols,
					pivot + 1,
					S.r);
			}
			else
			{
				int one = S.cols / 2;
				int two = S.cols - one;
				S1 = sketch(
					S.rows,
					one,
					S.l,
					pivot);
				S2 = sketch(
					S.rows,
					two,
					pivot + 1,
					S.r);
			}
			q.push(S1);
			q.push(S2);
		}
		fout << "NUMBER OF PARTITIONS: " << numberofgroups << '\n';

		fout << G.size() << "\n";
		for (int i = 0; i <= numberofgroups - 1; i++)
		{
			long long cnt = 0;
			for (auto it : G)
				if (it.second == i)
					cnt++;
			fout << i << " " << Gmatrices[i].rows << " " << Gmatrices[i].cols << ": " << cnt << "\n";
		}
	}
	outlierSketch = Gmatrix(outlier_rows, outlier_cols, K, P);

	clearAll();
}

Approach1::Approach1(string data_sample_file, int rows, int cols, int depth, int modulo, int w0, double C) : K(depth), P(modulo), w0(w0), C(C)
{
	setup(data_sample_file, rows, cols, depth, modulo, w0, C);
}

void Approach1::add(long long u, long long v, long long freq)
{
	if (mode == 0)
	{
		if (G.count(u))
		{
			Gmatrices[G[u]].add(u, v, freq);
		}
		else
		{
			outlierSketch.add(u, v, freq);
		}
	}
	else
	{
		if (G.count(v))
		{
			Gmatrices[G[v]].add(u, v, freq);
		}
		else
		{
			outlierSketch.add(u, v, freq);
		}
	}
}

long long Approach1::query(long long u, long long v)
{
	if (mode == 0)
	{
		if (G.count(u))
			return Gmatrices[G[u]].query(u, v);
		else
			return outlierSketch.query(u, v);
	}
	else
	{
		if (G.count(v))
			return Gmatrices[G[v]].query(u, v);
		else
			return outlierSketch.query(u, v);
	}
}

Approach1::Approach1() {}

unordered_set<pair<int, int>,HASH> Approach1::getHeavyHitterEdges(long long F)
{
	unordered_set<pair<int, int>,HASH> ret;
	for (auto it : Gmatrices)
	{
		auto cur = it.second.getHeavyHitterEdges(F);
		ret.insert(cur.begin(), cur.end());
	}
	auto cur = outlierSketch.getHeavyHitterEdges(F);
	ret.insert(cur.begin(), cur.end());
	return ret;
}


long long Approach1::queryNodeAggFreq(int x)
{
	if (mode == 0)
	{
		if (G.count(x))
			return Gmatrices[G[x]].queryNodeOutgoingFreq(x);
		else
			return outlierSketch.queryNodeOutgoingFreq(x);
	}
	else
	{
		if (G.count(x))
			return Gmatrices[G[x]].queryNodeIncomingFreq(x);
		else
			return outlierSketch.queryNodeIncomingFreq(x);
	}
}