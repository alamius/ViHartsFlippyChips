#ifndef FACE_HPP
#define FACE_HPP

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


#endif /* end of include guard: FACE_HPP */
