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

// include whatever classes you want

class Weighted_graph {
private:
	// your implementation here
	//  you can add both private member variables and private member functions

	static const double INF;

public:
	// Member Variables

	double **adjacency;
	bool *visited;
	double *cost;

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
	visited(new bool[n]),
	cost(new double[n]) {
	for (int i = 0; i < n; i++) {
		adjacency[i] = new double[n];
		// Will set visited and cost inside distance, since they need to be reset on every distance calculation
	}
	return;
}

Weighted_graph::~Weighted_graph() {
	// IMPLEMENTATION REQUIRED
}

int Weighted_graph::degree(int placeholder) const {
	// IMPLEMENTATION REQUIRED
	return 0;
}

int Weighted_graph::edge_count() const {
	// IMPLEMENTATION REQUIRED
	return 0;
}

double Weighted_graph::adjacent(int placeholder1, int placeholder2) const {
	// IMPLEMENTATION REQUIRED
	return 0;
}

double Weighted_graph::distance(int placeholder1, int placeholder2) {
	// IMPLEMENTATION REQUIRED
	return 0;
}

void Weighted_graph::insert(int placeholder1, int placeholder2, double placeholder3) {
	// IMPLEMENTATION REQUIRED
	return;
}

// You can modify this function however you want:  it will not be tested

std::ostream &operator<<(std::ostream &out, Weighted_graph const &graph) {
	return out;
}

// Is an error showing up in ece250.h or elsewhere?
// Did you forget a closing '}' ?

#endif