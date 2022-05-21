#ifndef FACE_HPP
#define FACE_HPP

#include <string>
#include <vector>
#include "Edge.hpp"

class Face : public std::vector<Edge*>{
	//extension of Edge* std::vector.
	//added: knows whether it is an outside or an inside face once that is set in Chip::make_faces
public:
	int inside = -1; //-1: not decided yet, 0: outside, 1: inside
	Face(){}; //edges is an empty array
	Face(std::vector<Edge*> edges_);
	std::string to_str();
	std::string dbg(string indent);
	virtual ~Face(){};
};


#endif /* end of include guard: FACE_HPP */
