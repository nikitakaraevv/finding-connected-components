#include <iostream> 
#include <list> 
#include <stack> 
#include <random>
#include <fstream>
#include <map>
#include <set>
using namespace std; 

class Graph {

	int n, counter, currComponent, biggestCompSize;
	list<int> *adj;
	map<int, list<int> > components;  

	public:

	Graph(int n);

	Graph(int n, float p);
    
	~Graph(){
		delete[] adj;
	}

	int getBiggestCompSize();

	Graph getTranspose();

	void addEdge(int v1, int v2);

	void DFS(int v, bool visited[]);

	void firstStage(stack<int> &s, int v, bool visited[]);

	static void timeTest();

	map<int, list<int> > SCC();
	
};
