#ifndef FACE_CPP
#define FACE_CPP

#include "Face.hpp"

Face::Face(std::vector<Edge*> edges_){
	for(int e = 0; e < edges_.size(); e++){
		push_back(edges_[e]);
	}
}
std::string Face::to_str(){
	std::string result = "(";
	for(int e = 0; e < size(); e++){
		result += char(80 + operator[](e)->from);
	}
	result += ")";
	return result;
}
std::string Face::dbg(string indent = ""){
	stringstream result;
	result << indent << "F({\n";
	for(int e = 0; e < size(); e++){
		result << indent << "  " << operator[](e)->dbg() << ", \n";
	}
	result << indent << "})";
	return result.str();
}


#endif /* end of include guard: FACE_CPP */
