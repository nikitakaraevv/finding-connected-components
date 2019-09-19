#include <iostream> 
#include <list> 
#include <stack> 
#include <random>
#include <fstream>
#include <map>
#include <set>
#include "scc.hpp"
using namespace std; 


Graph::Graph(int n) {
	this->n = n;
	this->adj = new list<int>[n];
}


Graph::Graph(int n, float p) {
	this->n = n;
	this->adj = new list<int>[n];
	std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(0.0,1.0);

	for (int v = 0; v < n; v++) { 
		for (int u = v+1; u < n; u++) { 
			double gen = distribution(generator);
			if (gen <= p) addEdge(u, v);
		}
	}
}


int Graph::getBiggestCompSize() {
	return biggestCompSize;
}

Graph Graph::getTranspose() {
	Graph gT(n);
	for (int v = 0; v < n; v++) {
		 list<int>::iterator it;
		 for (it = adj[v].begin(); it != adj[v].end(); it++)
			gT.adj[*it].push_back(v);
	}
	return gT;
}

void Graph::addEdge(int v1, int v2) {
	adj[v1].push_back(v2);
}

void Graph::DFS(int v, bool visited[]) { 
	visited[v] = true;
	this->counter++;
	components[currComponent].push_back(v);
	cout << v << " ";
	list<int>::iterator i;
	for (i = adj[v].begin(); i != adj[v].end(); i++) {
		if (!visited[*i]) {
			DFS(*i, visited);
		}
	}
}

void Graph::firstStage(stack<int> &s, int v, bool visited[]) { 
	visited[v] = true;
	list<int>::iterator i;
	for (i = adj[v].begin(); i != adj[v].end(); i++) {
		if (!visited[*i]) {
			firstStage(s, *i, visited);
		}
	}
	s.push(v);

}

map<int, list<int> > Graph::SCC() {
	stack<int> s;
	bool *visited = new bool[n];
	for (int i = 0; i < n; i++)
		visited[i] = false;

	for (int i = 0; i < n; i++) {
		if (!visited[i]) 
			firstStage(s, i, visited);
	}

	Graph gT = getTranspose();

	for (int i = 0; i < n; i++)
		visited[i] = false;

	gT.currComponent = 0;
	biggestCompSize = 0;
	while(!s.empty()) {
		int v = s.top();
		s.pop();
		if (!visited[v]) {
			gT.counter = 0;
			cout << endl << gT.currComponent++ << " component: " << endl;
			gT.DFS(v, visited);
			cout << endl << " number of nodes: " << gT.counter << endl;
			if (gT.counter > biggestCompSize) biggestCompSize = gT.counter;
			
		}
	}
	delete [] visited;
	return gT.components;
}
	

void Graph::timeTest() {

		int numOfPoints = 10;
		int power = 6;
		float p = 0.2;
		for (int i = 0; i < power; i++) {
		  	std::fstream fs;
		  	fs.open ("timeTestSCC.txt", std::fstream::in | std::fstream::out | std::fstream::app);
		  	numOfPoints *= 10;
		  	double avTime = 0.0f;
		  	for (int j = 0; j < 10; j++) {

		  		
				Graph g(numOfPoints,  p + 0.06 * j );
				auto start = chrono::system_clock::now();

				g.SCC();

				auto end = chrono::system_clock::now();
		 
		    	chrono::duration<double> elapsed_seconds = end - start;
		    	time_t end_time = chrono::system_clock::to_time_t(end);
		    	avTime += elapsed_seconds.count();
				cout << "finished computation at " << ctime(&end_time)
		              << "elapsed time: " << elapsed_seconds.count() << "s\n";
		  	}
		  	avTime /= 10.0f;
		  	fs << numOfPoints << " " << avTime << endl;
		  	fs.close();
		  }

	}

int main(int argc, char *argv[]) {
		cout << "Enter the file name" << endl;
		int x, source, target;
		ifstream inFile;
		inFile.open(argv[1]);
		if (!inFile) {
    		cerr << "Unable to open file " << argv[1] << endl;
    		exit(1);  
		}
		
		map<int, int> dict;  
		set<int> setOfIds;
		while (inFile >> x) {
			setOfIds.insert(x);
		}
		int n = setOfIds.size();
		cout << "n = " << n << endl;
		int i = 0;
		for (auto it = setOfIds.begin(); it != setOfIds.end(); it++)
			dict[*it] = i++;
  		
		Graph g(n); 
		
		
		inFile.close();
		inFile.open(argv[1]);
		
		while (inFile >> source >> target)  
			g.addEdge(dict[source], dict[target]); 

		inFile.close();
		map<int, list<int> > result;  
    	result = g.SCC();
    	cout << "Biggest component size: " << g.getBiggestCompSize() << endl;

    	Graph::timeTest();
		return 0;
	}