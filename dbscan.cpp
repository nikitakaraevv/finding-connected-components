#include <iostream>
#include <cmath>
#include <list>
#include <fstream>
#include <random>
#include <ctime> 
#include <stdio.h>
#include <float.h>
using namespace std;
#include "dbscan.hpp"


Points::Points(int _n, int _dim, int numOfClusters, double _eps, int _minNeigh) {
	n = _n;
	eps = _eps;
	dim = _dim;
	minNeigh = _minNeigh;
	double clusterSize = 15.0;
	double interClusterDist = 10.0;

	adj = new list<int>[n];
	std::default_random_engine generator(time(0));
	nodes = (Node *)calloc(n, sizeof(Node));
	for (int c = 0; c < numOfClusters; c++) {
  		std::uniform_real_distribution<double> distribution(0.0 + clusterSize * c, interClusterDist + clusterSize * c);
		for (int i = c * n / (double)numOfClusters; i < (c + 1) * n / (double)numOfClusters; i++) {
			nodes[i].coord = (double *)calloc(dim, sizeof(double));
			nodes[i].label = -1; 
			nodes[i].visited = 0; 
			for (int j = 0; j < dim; j++) {
				nodes[i].coord[j] = distribution(generator);
			}

		}

	}

	std::default_random_engine prob_generator(time(0));
		std::uniform_real_distribution<double> distribution(0.0, 1.0);
		double mindist = DBL_MAX;
		double currdist;
		for (int v = 0; v < n; v++) { 
		for (int u = v+1; u < n; u++) { 
			currdist = eudclidDist(nodes[u], nodes[v]);
			if (currdist < mindist)  mindist = currdist;
		}
	}

	for (int v = 0; v < n; v++) { 
		for (int u = 0; u < n; u++) { 
			if (u!=v){
				double gen = distribution(prob_generator);
				currdist = eudclidDist(nodes[u], nodes[v]);
				if (gen <= pow(mindist / pow(currdist, 4.0), 0.5)) addEdge(u, v);
			}
		}
	}
}


void Points::addEdge(int v1, int v2) {
	
	std::default_random_engine prob_generator(time(0));
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	double gen = distribution(prob_generator);

	if (gen <= 0.5) adj[v1].push_back(v2);
	else  adj[v2].push_back(v1);
}



Points::Points(double _eps, int _minNeigh, string dsname) {
	
   ifstream inFile;
   
   inFile.open(dsname);
  
   double one, two, three, four, eps;
   string label;
   
   inFile >> n >> dim;
   cout << n << endl;
   eps = _eps;
   minNeigh = _minNeigh;
   nodes = (Node *)calloc(n, sizeof(Node));
  
   int i = 0;
   while (inFile >> one >> two) {
        nodes[i].coord = (double *)calloc(dim, sizeof(double));
        nodes[i].coord[0] = one;
        nodes[i].coord[1] = two;
        nodes[i].label = -1; 
        nodes[i].visited = 0; 
        i++;
   } 

   inFile.close();


}

void Points::print() {
	for (int i = 0; i < n; i++) {
		cout << i << " point, coords " << nodes[i].coord[0] << 
		" " << nodes[i].coord[1] <<  " label " << nodes[i].label
		 << " visited " << nodes[i].visited << endl;
	}

}


list<int> Points::getNumNeighbors(int ndNum) {
	list<int> neighborsList;

	for (int i = 0; i < n; i++) {
		if (i != ndNum) {
			double dist = 0;
			for (int j = 0; j < dim; j++) {
				dist += pow(nodes[i].coord[j] - nodes[ndNum].coord[j] , 2);

			}
			dist = pow(dist, 0.5);

			if (dist <= eps) {
				neighborsList.push_back(i);
			}

		}
	}
	return neighborsList;

}

bool Points::compare(Node nd1, Node nd2) {
	if ( nd1.coord == nd2.coord && nd1.visited == nd2.visited)
		return true;
	return false;
}

int Points::noiseCount() {
	int count = 0;
	for (int i = 0; i < n; i++) {
		if (nodes[i].label < 0)
			count++;
	}
	return count;
}

int Points::dbscan() {
	int c = 0;
	int lst_size = 0;
	for (int i = 0; i < n; i++) {
		if (nodes[i].visited == 1) {
			continue;
		} else {
			nodes[i].visited = 1;
			list<int> neighbors = getNumNeighbors(i);
			lst_size = neighbors.size();
			if (lst_size < minNeigh)
				nodes[i].label = -1; //noise
			else {
				c++; 
				nodes[i].label = c;
				while (!neighbors.empty()) {
					if (nodes[neighbors.front()].visited == 0) {
						nodes[neighbors.front()].visited = 1;
						
						list <int> ngh = getNumNeighbors(neighbors.front());
						int nghsz = ngh.size();
						if (nghsz >= minNeigh) {
							for (int j = 0; j < nghsz; j++) {
								neighbors.push_back(ngh.front());
								ngh.pop_front();
							}
						}
					}
					if (nodes[neighbors.front()].label  < 0)
						nodes[neighbors.front()].label = c;
					neighbors.pop_front();
				}
			}
		}
	}
	return c;
}

double Points::eudclidDist(Node Node1, Node Node2) {
	double dist = 0.0f;
	for (int j = 0; j < dim; j++)
		dist += pow(Node1.coord[j] - Node2.coord[j], 2);
	return pow(dist, 0.5);
}

void Points::knnDist(int k) {
	double *distances;
	averageKnnDistances = new double[n];
	double currDist;
	for (int i = 0; i < n; i++) {
		distances = new double[n];
		currDist = 0.0f;
		for (int j = 0; j < n; j++) 
			distances[j] = eudclidDist(nodes[i], nodes[j]);
		sort(distances, distances + sizeof(distances) / sizeof(distances[0]));
		for (int neigh = 0; neigh < k+1; neigh++)
			currDist += distances[neigh];
		
		currDist /= k;
		averageKnnDistances[i] = currDist;
	}

}


void Points::writeClustering() {
	  string fileName = "dbscan_clustering.txt";
	  if (remove(fileName.c_str( )) != 0)
	       cout << "Remove " << fileName << " operation failed" << endl;
	  else
	       cout << fileName << " has been removed." << endl;
	  std::fstream fs;
	  fs.open (fileName, std::fstream::in | std::fstream::out | std::fstream::app);
	  

	  for (int i = 0; i < n; i++) {
	  	for (int j = 0; j < dim; j++) {
	  		fs <<nodes[i].coord[j];
	  		fs << " ";
	  	}
	  	fs << nodes[i].label << endl;
	  }
	  fs.close();

	  cout << fileName << " has been created." << endl;
  
}


void Points::writeEdges() {
	  string fileName = "edges.txt";
	  if (remove(fileName.c_str( )) != 0)
	       cout << "Remove " << fileName << " operation failed" << endl;
	  else
	       cout << fileName << " has been removed." << endl;
	  std::fstream fs;
	  fs.open (fileName, std::fstream::in | std::fstream::out | std::fstream::app);

	  for (int i = 0; i < n; i++) {
	  	for (list<int>::iterator it = adj[i].begin(); it != adj[i].end(); it++)
	  		fs << i << " " << *it << endl;
	  }

	  fs.close();
	  cout << fileName << " has been created." << endl;
}


void Points::writeKnnDists() {

	  string fileName = "knnDistances.txt";
	  if (remove(fileName.c_str( )) != 0)
	       cout << "Remove " << fileName << " operation failed" << endl;
	  else
	       cout << fileName << " has been removed." << endl;
	  std::fstream fs;
	  fs.open (fileName, std::fstream::in | std::fstream::out | std::fstream::app);

	  for (int i = 0; i < n; i++) {
	  	for (int j = 0; j < dim; j++) {
	  		fs << nodes[i].coord[j] << " ";
	  	}
	  	fs << averageKnnDistances[i] << endl;
	  }
	  fs.close();
	  cout << fileName << " has been created." << endl;
}


void Points::timeTest(double eps, int neigh){

	int numOfPoints = 10;
	int power = 5;
	for (int i = 0; i < power; i++) {
	  	std::fstream fs;
	  	fs.open ("timetest.txt", std::fstream::in | std::fstream::out | std::fstream::app);
	  	numOfPoints *= 10;
	  	double av_time = 0.0f;
	  	for (int j = 0; j < 10; j++) {

	  		Points pnts(numOfPoints, 2, 3, eps + (eps / 10.0) * (j - 5), neigh);
		
			auto start = chrono::system_clock::now();

			int c = pnts.dbscan();

			auto end = chrono::system_clock::now();
	 
	    	chrono::duration<double> elapsed_seconds = end - start;
	    	time_t end_time = chrono::system_clock::to_time_t(end);
	    	av_time+=elapsed_seconds.count();
			cout << "finished computation at " << ctime(&end_time)
	              << "elapsed time: " << elapsed_seconds.count() << "s\n";
	  	}
	  	av_time /= 10.0f;
	  	fs << numOfPoints << " " << av_time << endl;
	  	fs.close();
	  }

}


int main(int argc, char *argv[]) {
	if (argc < 3)
		cout << "Enter eps end min number of neighbors" << endl;
	double eps = 0;
	eps = atof(argv[1]);
	int neigh = 0;
	neigh = atoi(argv[2]);
	if (argc > 3)
		string filename = argv[3];
	//Points::time_test(eps, neigh);
	//eps = 1;
	//neigh = 6000;

	Points pnts(100, 2, 3, eps, neigh);
	//pnts.print();
	//pnts.knnDist(5);
	//pnts.write_knn_dists();
	//auto start = chrono::system_clock::now();
    // Some computation here
	//Points pnts(eps, neigh, filename);
	pnts.print();
	pnts.knnDist(5);
	pnts.writeKnnDists();
	pnts.writeEdges();
	int c = pnts.dbscan();
	pnts.writeClustering();
	cout << "number of clusters " << c << endl;
	//int c = pnts.dbscan(); 
}