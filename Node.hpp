#ifndef NODE_HPP
#define NODE_HPP

#include "include/canvas/Vector.hpp"
#include "Edge.hpp"

class Node {
	//a node represents an intersection within the Line of a Chip:
	//  it contains the Vector of the intersection,
	//  but also the corresponding Splines in the Chip, by their index in Chip::points,
	//  the parameters on these Splines and
	//  the edges that go out to other nodes.
	//  a node also saves the information, whether the two axes that go through it have to be flipped to draw a correct image
	//  it also contains an index that is set when the nodes along the Line of the Chip is sorted
	//a special version of this object is constructed to save the intersections with a straight line, where q is always -1 and u refers to the line instead of a Spline.
public:
	Vector value; //the vector of the intersection
	int p; //index of the first spline's begin point in the line (Chip::points[p])
	int q; //index of the second spline's begin point  (Chip::points[q])
	float t; //parameter of the intersection along the first spline (value = Spline(Chip::points[p] ~~ Chip::points[p+1])(t))
	float u; //parameter along the second spline
	Edge* edges[4];
	bool edges_walked[4];
	bool flipped; //tells the Chip whether the Spline on points[p] is left or right of the Spline on points[q]
	int index = -1; //tells the Chip, whether the Line at points[p] is dark or not (and points[q] bright respectivly), by looking at index % 2
	Node(){};
	Node(Vector value_, int p_, float t_, int q_, float u_, const std::vector<Point>& points);
	//only for intersections between a Spline and a straight line!!
	Node(Vector value_, int n_, int e_, float t_, float u_);
	std::string dbg();
	bool operator>(Node other);
	bool operator<(Node other);
	virtual ~Node(){};
};

#endif /* end of include guard: NODE_HPP */
