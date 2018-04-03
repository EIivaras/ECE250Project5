/*****************************************
* Instructions
*  - Replace 'uwuserid' with your uWaterloo User ID
*  - Select the current calendar term and enter the year
*  - List students with whom you had discussions and who helped you
*
* uWaterloo User ID:  zwalford @uwaterloo.ca
* Submitted for ECE 250
* Department of Electrical and Computer Engineering
* University of Waterloo
* Calender Term of Submission:  (Winter|Spring|Fall) 201N
*
* By submitting this file, I affirm that
* I am the author of all modifications to
* the provided code.
*
* The following is a list of uWaterloo User IDs of those students
* I had discussions with in preparing this project:
*    -
*
* The following is a list of uWaterloo User IDs of those students
* who helped me with this project (describe their help; e.g., debugging):
*    -enruizno (test files)
*    -gndaryee (test files)
*****************************************/

#ifndef WEIGHTED_GRAPH_H
#define WEIGHTED_GRAPH_H

#ifndef nullptr
#define nullptr 0
#endif

#include <iostream>
#include <limits>
#include "Exception.h"

// include whatever classes you want

class Weighted_graph {
private:
	// your implementation here
	//  you can add both private member variables and private member functions
	double **adjacency;
	int *degree_array;
	int edges;
	int vertices;

	static const double INF;

public:

	// Member Fxns

	Weighted_graph(int = 50);
	~Weighted_graph();

	int degree(int) const;
	int edge_count() const;
	double adjacent(int, int) const;
	double distance(int, int);

	void insert(int, int, double);

	// Friends

	friend std::ostream &operator<<(std::ostream &, Weighted_graph const &);
};

const double Weighted_graph::INF = std::numeric_limits<double>::infinity();

// Your implementation here

//////////////////////////////////////////////////////////////////////////////
// PUBLIC MEMBER FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

Weighted_graph::Weighted_graph(int n):
	adjacency(new double*[n]),
	degree_array(new int[n]),
	edges(0),
	vertices(n){
	for (int i = 0; i < n; i++) {
		adjacency[i] = new double[n];
		// Will set visited and cost inside distance, since they need to be reset on every distance calculation
		degree_array[i] = 0;
		for (int j = 0; j < n; j++) {
			if (i == j)
			{
				adjacency[i][j] = 0; // A node is a distance of 0 from itself
			}
			else {
				adjacency[i][j] = INF; // Each node starts with a distance of infinity between itself and any other node
			}
		}
	}
	return;
}

Weighted_graph::~Weighted_graph() {
	for (int i = 0; i < vertices; i++) {
		delete[] adjacency[i];
	}
	delete[] adjacency;
	delete[] degree_array;
	edges = 0;
	vertices = 0;
}

int Weighted_graph::degree(int node) const {
	if ((node < 0) || (node >= vertices)) { // Only criteria is that the node that is being asked for exists
		throw illegal_argument();
	}
	return degree_array[node];
}

int Weighted_graph::edge_count() const {
	return edges;
}

double Weighted_graph::adjacent(int node1, int node2) const {
	if ((node1 < 0 || node2 < 0) || (node1 >= vertices || node2 >= vertices)) { // Make sure the two nodes exist
		throw illegal_argument();
	}
	return adjacency[node1][node2];
}

double Weighted_graph::distance(int node1, int node2) {
	if ((node1 < 0 || node2 < 0) || (node1 >= vertices || node2 >= vertices)) { // Make sure the nodes we're finding the distance between exist
		throw illegal_argument();
	}
	bool *visited = new bool[vertices]; // Array that will track whether a node has been visited in the algorithm
	double *cost = new double[vertices]; // Array that will track the current shortest distance from the starting node to any other node
	double shortest_distance = INF; // The current shortest distance for a SINGLE iteration
	bool finished = false; // We're going to use a while loop to run the algorithm, and this is the variable that will be the continuing condition

	for (int i = 0; i < vertices; i++) {
		visited[i] = false; // Every node starts unvisited by the algorithm
		// The cost array is initialized below
	}

	// STEP 1: Initialize node1's visited status to true (cost = 0 is taken care of at formation of adjacency matrix)
	visited[node1] = true;
	// STEP 2: Look through all nodes node1 (the node we're basing all connections off) is adjacent to, update the cost array
	for (int i = 0; i < vertices; i++) {
		cost[i] = adjacency[node1][i]; // Inherently takes care of all distances of INF and 0, as this is done for the adjacency matrix in the constructor
	}

	double current_shortest_path;
	int node_csp;
	while (!finished) {
		current_shortest_path = INF;
		node_csp = -1; // node current shortest path - we initialize to -1, so that we can check if it's -1 after to determine whether or not any new shortest path to an unvisited node has been located
		// STEP 3: Step through the cost array and find the node that is the shortest path from the starting node, provided it has not been visited already
		for (int i = 0; i < vertices; i++) {
			if (visited[i] != true && cost[i] < current_shortest_path) {
				current_shortest_path = cost[i];
				node_csp = i;
			}
		}

		// STEP 4: Check that we've found a shortest distance not visited; if we've visited everything, let's break the loop
		if (current_shortest_path == INF) {
			finished = true;
			continue;
		}

		// STEP 5: For the node with the shortest distance, mark it as visited, and look through all the nodes it is adjacent to and update the cost array accordingly
		visited[node_csp] = true;

		for (int i = 0; i < vertices; i++) {
			if ((adjacency[node_csp][i] + current_shortest_path) < cost[i]) {
				cost[i] = adjacency[node_csp][i] + current_shortest_path;
			}
		}
	}

	shortest_distance = cost[node2]; // The shortest distance should be the cost stored to get to node2, as the shortest distance is being updated with every node found
	std::cout << "Shortest Distance: " << shortest_distance << std::endl;

	// Cleanup
	delete[] visited;
	delete[] cost;
	return shortest_distance;
}

void Weighted_graph::insert(int node1, int node2, double weight) {
	if ((node1 < 0 || node2 < 0) || (node1 >= vertices || node2 >= vertices)) { // Make sure the vertices exist
		throw illegal_argument();
	}
	else if (node1 == node2){ // Make sure we're not inserting an edge between the same vertice
		throw illegal_argument();
	}
	else if (weight <= 0) { // Make sure we're not selecting a weight of zero or less
		throw illegal_argument();
	}
	if (adjacency[node1][node2] == INF || adjacency[node1][node2] == 0) { // As long as an edge doesn't already exist between the two nodes, we update the edge count and degree
		degree_array[node1]++;
		degree_array[node2]++;
		edges++;
	}
	adjacency[node1][node2] = weight;
	adjacency[node2][node1] = weight;
	return;
}

// You can modify this function however you want:  it will not be tested

std::ostream &operator<<(std::ostream &out, Weighted_graph const &graph) {
	for (int i = 0; i < graph.vertices; i++) {
		out << "||";
		for (int j = 0; j < graph.vertices; j++) {
			out << "     " << graph.adjacency[i][j];
		}
		out << "      ||" << std::endl;
	}
	return out;
}

// Is an error showing up in ece250.h or elsewhere?
// Did you forget a closing '}' ?

#endif