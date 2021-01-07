#ifndef EDGE_HPP
#define EDGE_HPP

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
    std::string dbg(bool spline_dbg);
    virtual ~Edge(){};
};
std::string Edge::dbg(bool spline_dbg = false){
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


#endif /* end of include guard: EDGE_HPP */
