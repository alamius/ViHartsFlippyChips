#ifndef EDGE_CPP
#define EDGE_CPP

#include "Edge.hpp"
#include <sstream>

std::string Edge::dbg(bool spline_dbg){
	stringstream result;
	result << "E(" << char(80 + from) << out << " -- " << char(80 + to) << in << ")";
	if(spline_dbg) result << ": " << S.dbg();
	return result.str();
}
bool Edge::equal(Edge* other){
	return (
		from == other->from &&
		to   == other->to   &&
		in   == other->in   &&
		out  == other->out
	) || (
		from == other->to   &&
		to   == other->from &&
		in   == other->out  &&
		out  == other->in
	);
}


#endif /* end of include guard: EDGE_CPP */
