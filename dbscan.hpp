#include <iostream>
#include <cmath>
#include <list>
#include <fstream>
#include <random>
#include <ctime> 
#include <stdio.h>
#include <float.h>
using namespace std;



struct Node {

	double *coord;
	int label;
	int visited;
};

class Points {

	double *averageKnnDistances;
	list<int> *adj;
	Node *nodes;
	double eps;
	int minNeigh;

public:
	int n;

	int dim;

	Points(int _n, int _dim, int numOfClusters, double _eps, int _minNeigh);

	void addEdge(int v1, int v2);

	Points(double _eps, int _minNeigh, string dsname);

	~Points() {
		delete[] nodes;
		delete[] adj;
	}


	void print();

	list<int> getNumNeighbors(int ndNum);

	bool compare(Node nd1, Node nd2);

	int noiseCount();

	int dbscan();

	double eudclidDist(Node Node1, Node Node2);

	void knnDist(int k);

	std::fstream openFile(string fileName);

	void writeClustering();


	void writeEdges();


	void writeKnnDists();


	void timeTest(double eps, int neigh);
};

