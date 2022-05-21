#ifndef NODE_CPP
#define NODE_CPP

#include "Node.hpp"

Node::Node(Vector value_, int p_, float t_, int q_, float u_, const std::vector<Point>& points){
	value = value_;
	p = p_;
	t = t_;
	q = q_;
	u = u_;
	Vector a = Spline(
		points[p],
		points[(p + 1) % points.size()]
	).dL(t);
	Vector b = Spline(
		points[q],
		points[(q + 1) % points.size()]
	).dL(u);
	flipped = a.x*b.y - a.y*b.x < 0;
}
Node::Node(Vector value_, int n_, int e_, float t_, float u_){
	//only for intersections Spline - straight Line (Spline/SplineConstruct::intersect_linear) !!
	value = value_;
	p = n_;
	q = e_;
	t = t_;
	u = u_;
}
std::string Node::dbg(){
	stringstream ss;
	ss  << "N("
		<< char(65 + p) << ": " << t << " x "
		<< char(65 + q) << ": " << u
	<< ")";
	return ss.str();
}
bool Node::operator>(Node other){
	return (
		(p > other.p) ||
		(
			p == other.p &&
			t  > other.t
		)
	);
}
bool Node::operator<(Node other){
	return (
		(p < other.p) ||
		(
			p == other.p &&
			t  < other.t
		)
	);
}


#endif /* end of include guard: NODE_CPP */
