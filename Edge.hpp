#ifndef EDGE_HPP
#define EDGE_HPP

#include "include/spline/SplineConstruct.hpp"

class Edge{
public:
	int from, to; //nodes indices
	int out, in; //the numbers of the connections in the nodes from and to.
	SplineConstruct S;
	Edge(int from_, int to_, int out_, int in_, SplineConstruct S_){
		from = from_;
		to = to_;
		out = out_;
		in = in_;
		S = S_;
	};
	bool equal(Edge*);
	std::string dbg(bool spline_dbg = false);
	virtual ~Edge(){};
};


#endif /* end of include guard: EDGE_HPP */
