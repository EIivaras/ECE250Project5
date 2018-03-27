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
*    -
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
				adjacency[i][j] = 0;
			}
			else {
				adjacency[i][j] = INF;
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
	if ((node < 0) || (node >= vertices)) {
		throw illegal_argument();
	}
	return degree_array[node];
}

int Weighted_graph::edge_count() const {
	return edges;
}

double Weighted_graph::adjacent(int node1, int node2) const {
	if ((node1 < 0 || node2 < 0) || (node1 >= vertices || node2 >= vertices)) {
		throw illegal_argument();
	}
	return adjacency[node1][node2];
}

double Weighted_graph::distance(int node1, int node2) {
	if ((node1 < 0 || node2 < 0) || (node1 >= vertices || node2 >= vertices)) {
		throw illegal_argument();
	}
	bool *visited = new bool[vertices];
	double *cost = new double[vertices];
	double shortest_distance = INF;
	bool finished = false;

	for (int i = 0; i < vertices; i++) {
		visited[i] = false;
		cost[i] = INF;
	}

	// STEP 1: Initialize node1's visited status to true (cost = 0 is taken care of at formation of adjacency matrix)
	visited[node1] = true;
	// STEP 2: Look through all nodes node1 is adjacent to, update the cost array
	for (int i = 0; i < vertices; i++) {
		cost[i] = adjacency[node1][i];
		std::cout << cost[i] << std::endl;
	}

	double current_shortest_path;
	int node_csp;
	while (!finished) {
		current_shortest_path = INF;
		node_csp = -1; // node current shortest path
		// STEP 3: Step through the cost array and find the node that is the shortest path from the starting node, provided it has not been visited already
		for (int i = 0; i < vertices; i++) {
			if (visited[i] != true && cost[i] < INF) {
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

		// NOTE: Costs need to be updated if shorter paths to already discovered nodes are found
		// NOTE: When looking at the third node in a string (if you got 3 shortest distance's in a row), you should save a sum of the distance you've travelled so far
		// NOTE: Every time you start at an existing node and check its adjacents, make the saved distance equal to that node's distance in the cost array, and THEN continue
			// - Ideally, if you've been updating the cost array each time a new path is found, this should be the current shortest distance to that node
	}

	shortest_distance = cost[node2];
	std::cout << "Shortest Distance: " << shortest_distance << std::endl;

	// Cleanup
	delete visited;
	delete cost;
	return shortest_distance;
}

void Weighted_graph::insert(int node1, int node2, double weight) {
	if ((node1 < 0 || node2 < 0) || (node1 >= vertices || node2 >= vertices)) {
		throw illegal_argument();
	}
	else if (node1 == node2){
		throw illegal_argument();
	}
	else if (weight <= 0) {
		throw illegal_argument();
	}
	adjacency[node1][node2] = weight;
	adjacency[node2][node1] = weight;
	degree_array[node1]++;
	degree_array[node2]++;
	edges++;
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