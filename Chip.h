#ifndef CHIP_H
#define CHIP_H

#include "Edge.hpp"
#include "Node.hpp"
#include "Face.hpp"

colorint* std_color_func(colorint result[COLOR_LEN], float t, float v){
    make_color(
        result,
        255.0f*v,
        255.0f*(1.0f - t)*(1.0f - t)
    );
    return result;
}
colorint std_alpha_map(colorint alpha){
    if(alpha > 0){
        return 255;
    }else{
        return 0;
    }
}

Vector interpolate(Vector p, Vector q, float t){
    return p * (1.0f - t) + q * t;
}

class Chip{
private:
    //points (Vector of position and Vector of direction combined, defined in Spline.h) of the Chip. the line of the chip is the combination of the splines from every point to the next
    std::vector<Point> points;
    //will be created at the intersections of the line defined by points, containing the intersection Vector and the indices of the points and correcponding parameters (see Nodes)
    std::vector<Node> nodes;
    //the faces of the chip with their edges
    std::vector<Face> faces;
public:
    //copies the points
    Chip(std::vector<Point> points_);
    //colors all faces that are Face::inside with t_prec*v_prec small quadrilaterals with color given by color_func(t, v)
    void color(
        int t_prec,
        int v_prec,
        colorint* (*color_func)(colorint result[COLOR_LEN], float t, float v),
        bool apply_gauss
    );
    //colors a area defined by four Splines in t_prec*v_prec small quadrilaterals with color given by color_func(t, v)
    template <typename CanvasT>
    void color_stripe(
        CanvasT* C,
        int t_prec,
        int v_prec,
        colorint* (*color_func)(colorint result[COLOR_LEN], float t, float v),
        const Vector& P, const Vector& p1, const Vector& p2,
        const Vector& Q, const Vector& q1, const Vector& q2,
        const Vector& R, const Vector& r1, const Vector& r2,
        const Vector& S, const Vector& s1, const Vector& s2
    );
    //returns all intersections of the Line with itself. mistakes might be caused by Spline::intersect which has faults
    std::vector<Vector> intersect();
    //returns all intersections of the Line with a straight line A + a*t for t in R
    std::vector<Vector> intersect_linear(Vector A, Vector a);
    //subroutine of make_edges: takes intersection and follows the Line forewards of backwards
    void follow_edge(int p_curr, float t_curr, int sign, int d, int from, bool SplineConstruct_approximate);
    //looks through all Nodes that are set by Chip::intersect and finds the Edges that connect them (every Node has four Edges)
    void make_edges(bool SplineConstruct_approximate);
    //goes through all Edges, following them and connecting them into faces
    void make_faces();
    //draw the Line defined by Chip::points onto Canvas C
    template <typename CanvasT>
    void draw(CanvasT* C, int samples = 30);
    //draw a net of the Line onto Canvas C
    template <typename CanvasT>
    void draw_net(
        CanvasT* C,
        int t_prec = 5,
        int v_prec = 5,
        int face_from = 0,
        int face_to = -1, // -1 is turned into faces.size() - 1
        bool draw_E = true,
        bool draw_F = false
    );
    //marking the private Chip::points on Canvas C
    template <typename CanvasT>
    void mark_points(CanvasT* C);
    //transform all Points in Chip::points according to matrix
    void transform(float a, float b, float c, float d, float e, float f); // [x, y][[a, b], [c, d]] + [e, f]
    //return a latex compatible text version of the Line
    std::string latex();
    std::string dbg(std::string indent);
    virtual ~Chip(){};
};
Chip::Chip(std::vector<Point> points_){
    for(int p = 0; p < points_.size(); p++){
        points.push_back(points_[p]);
    }
}
void Chip::color(
    int t_prec = 30, //number of separate layers of quadrilaterals from the corner outwards
    int v_prec = 30, //number of separate layers of quadrilaterals across one corner
    colorint* (*color_func)(colorint result[COLOR_LEN], float t, float v) = &std_color_func, //tells the color for quadrilateral at (t, v)
    bool apply_gauss = true //whether the resulting image should be passed through a gaussian filter
){
    intersect();
    //catching unintersecting Line (error is printed by intersect)
    if(nodes.size() == 0) return;
    //the boolean argument is SplineConstruct_approximate: whether the full SplineConstruct along one edge should be approximated by a Spline (currently needed because coloring uses Splines defined from dL(0) and dL(1))
    make_edges(true);
    make_faces();

    //the canvas LC is used several times to paint separate layers that are written to hard drive in between and then all loaded onto a BasicCanvas in the end
    LC = new LayeredCanvas();
    //there will be as many layers as the maximum number of edges one face has (because for every edge there is a corner and corners overlap and are therefore written to separate layers)
    int max_edges = 0;
    for(int f = 0; f < faces.size(); f++){
        if(faces[f].size() > max_edges && faces[f].inside == 1){
            max_edges = faces[f].size();
        }
    }
    //variables for face coloring
    Vector P, Q, R;
    Vector p1, p2, q1, q2, r1, r2;
    //the Splines that actually define the quadrilaterals (stored to be reusable for a whole row)
    Spline* F[v_prec];
    //running variables: to from P (center of the corner) out to Q/R and v from PQ other PR
    float t, v;
    float dt = 1.0f/t_prec;
    float dv = 1.0f/v_prec;
    //color_func writes to this
    colorint fillcolor[COLOR_LEN];
    //readability variables
    Face* face;
    Edge* ePQ,* eQ_,* e_R,* eRP,* eQP,* ePP; //Edges from P to Q, Q to other, other to R, R to P, Q to P and P to P
    //going through the layers, writing the corners and then saving and clearing the canvas
    for(int layer = 0; layer < max_edges; layer++){
        for(int f = 0; f < faces.size(); f++){
            face = &(faces[f]);
            if(!face->inside) continue;
            if(face->size() == 1){
                if(layer != 0) continue; //single-edges faces (loops) go to layer 0
                ePP = face->at(0);
                if(dbg_file_lvl >= 2) std::cout << "color face " << face->to_str() << " edge " << char(80 + ePP->from) << char(80 + ePP->to) << '\n';
                P = nodes[ePP->from].value;
                Q = ePP->S(.5);
                p1 = ePP->S.dL(0);
                p2 = ePP->S.dL(1);
                q1 = ePP->S.dL(.5);
                q2 = q1.mult_complex(Vector(0, -1)); //turn right 90°
                F[0] = new Spline(P, P, Vector(0, 0), Vector(0, 0));
                for(int v_ = 1; v_ <= v_prec; v_++){
                    v = dv*v_;
                    F[v_] = new Spline(P, P, p1*v, p2*pow(v, 1.3));
                    for(float t = 0; t < 1; t += dt){
                        LC->setcolor(color_func(fillcolor, t, t));
                        LC->quadrilateral_unchecked(
                            (*(F[v_-1]))(t+dt),
                            (*(F[v_-1]))(t),
                            (*(F[v_  ]))(t+dt),
                            (*(F[v_  ]))(t)
                        );
                    }
                }
                if(dbg_file_lvl >= 4)
                    LC->write(
                        filename+"_f"+_to_str(f)+"="+face->to_str()+"_e"+char(80+ePP->from)+char(80+ePP->to)
                    );
                continue;
            }else if(face->size() == 2){
                if(dbg_file_lvl >= 2) std::cout << "color face " << face->to_str() << '\n';
                for(int e = 0; e < face->size(); e++){
                    if(e != layer) continue; //only edges of the same index as the current canvas layer are drawn
                    if(dbg_file_lvl >= 2) std::cout << "  color by edge " << char(80 + face->at(e)->from) << char(80 + face->at(e)->to) << '\n';
                    int g = (e + 1) % 2;
                    ePQ = face->at(e);
                    eQP = face->at(g);
                    P = nodes[ePQ->from].value;
                    Q = nodes[eQP->from].value;
                    color_stripe(
                        LC, t_prec, v_prec, color_func,
                        P,  ePQ->S.dL(0), Vector(0, 0),
                        Q,  ePQ->S.dL(1), Vector(0, 0),
                        P, -eQP->S.dL(1), Vector(0, 0),
                        Q, -eQP->S.dL(0), Vector(0, 0)
                    );
                    if(dbg_file_lvl >= 4)
                        LC->write(
                            filename+"_f"+_to_str(f)+"="+face->to_str()+"_e"+char(80+face->at(e)->from)+char(80+face->at(e)->to)
                        );
                }
                continue;
            }else{
                if(dbg_file_lvl >= 2) std::cout << "color face " << face->to_str() << '\n';
                for(int e = 0; e < face->size(); e++){
                    if(e != layer) continue; //only edges of the same index as the current canvas layer are drawn
                    if(dbg_file_lvl >= 2) std::cout << "  color by edge " << char(80 + face->at(e)->from) << char(80 + face->at(e)->to) << '\n';
                    //edge from R to P
                    ePQ = face->at(e);
                    eQ_ = faces[f][(e + 2*face->size() + 1) % face->size()];
                    e_R = face->at((e + 2*face->size() - 2) % face->size());
                    eRP = face->at((e + 2*face->size() - 1) % face->size());
                    P = nodes[face->at(e)->from].value;
                    Q = nodes[face->at(e)->to  ].value;
                    R = nodes[    eRP    ->from].value;
                    color_stripe(
                        LC, t_prec, v_prec, color_func,
                        P,  ePQ->S.dL(0),   Vector(0, 0),
                        Q,  ePQ->S.dL(1),   eQ_->S.dL(0),
                        P, -eRP->S.dL(1),   Vector(0, 0),
                        R, -eRP->S.dL(0),   e_R->S.dL(1)
                    );
                    if(dbg_file_lvl >= 4)
                        LC->write(filename+"_f"+_to_str(f)+"="+face->to_str()+"_e"+char(80+face->at(e)->from)+char(80+face->at(e)->to));
                }
            }
        }
        if(dbg_file_lvl >= 3){
            std::cout << "written layer " << layer << " of " << max_edges << "; size == " << sizeof(BasicCanvas) << '\n';
            LC->write(filename+"_layer"+_to_str(layer), *write_bg_color);
            if(dbg_file_lvl >= 4) LC->dump(filename+"_layer"+_to_str(layer));
        }
        LC->save(filename+"_layer"+_to_str(layer));
        LC->clear();
    }
    LC->clear();
    if(dbg_file_lvl >= 2) std::cout << "adding image layers together:" << '\n';
    if(dbg_file_lvl >= 4) LC->dump(filename+"0_bg");
    for(int bc = 0; bc < max_edges; bc++){
        LC->load(filename+"_layer"+_to_str(bc), 0, true);
        if(dbg_file_lvl >= 4) LC->dump(filename+_to_str(bc+1)+"_loaded");
    }
    if(apply_gauss){
        LC->filter(kernel_gauss, true);
    }
    LC->save(filename+"_added", std_alpha_map);
    delete LC;
    BC = new BasicCanvas();
    BC->background(*write_bg_color);
    BC->load(filename+"_added");
    BC->write(filename);
    if(dbg_file_lvl >= 4) BC->dump(filename);
    delete BC;
}
template <typename CanvasT>
void Chip::color_stripe(
    CanvasT* C,
    int t_prec,
    int v_prec,
    colorint* (*color_func)(colorint result[COLOR_LEN], float t, float v),
    const Vector& P, const Vector& p1, const Vector& p2,
    const Vector& Q, const Vector& q1, const Vector& q2,
    const Vector& R, const Vector& r1, const Vector& r2,
    const Vector& S, const Vector& s1, const Vector& s2
){
    // P 1 ~~ Q 1
    // 2      2
    // |      |
    // R 1 ~~ S 1
    // 2      2
    Spline* PR = new Spline(P, R, p2, r2);
    Spline* QS = new Spline(Q, S, q2, s2);
    float v;
    float dv = 1.0f/v_prec;
    float dt = 1.0f/t_prec;
    Spline* F[v_prec+1];
    F[0] = new Spline(P, Q, p1, q1);
    colorint fillcolor[COLOR_LEN];
    for(int v_ = 1; v_ <= v_prec; v_++){
        v = v_*dv;
        F[v_] = new Spline(
            PR->L(v),
            QS->L(v),
            interpolate(p1, r1, v),
            interpolate(q1, s1, v)
        );
        for(float t = 0; t <= 1.0f-dt; t += dt){
            C->setcolor(color_func(fillcolor, t, v));
            C->quadrilateral_unchecked(
                (*(F[v_-1]))(t+dt),
                (*(F[v_-1]))(t),
                (*(F[v_  ]))(t+dt),
                (*(F[v_  ]))(t)
            );
        }
    }
}
std::vector<Vector> Chip::intersect(){
    if(dbg_file_lvl >= 4) std::cout << "intersections of:" << '\n';
    Spline* P;
    Spline* Q;
    Intersection intersections[3];
    int intersections_len;
    std::vector<Node> nodes_ = std::vector<Node>(); //unsorted
    std::vector<Vector> result;
    bool neighbors, identical;
    for(int p = 0; p < points.size(); p++){
        if(dbg_file_lvl >= 4) std::cout << "  spline following " << char(65 + p) << '\n';
        for(int q = p; q < points.size(); q++){
            if(dbg_file_lvl >= 4) std::cout << "    with spline following " << char(65 + q) << '\n';
            intersections_len = 0;
            P = new Spline(points[p], points[(p+1) % points.size()]);
            Q = new Spline(points[q], points[(q+1) % points.size()]);
            neighbors = (
                (q-p == 1) ||
                (p == 0) && (q == points.size()-1)
            );
            identical = (p == q);
            P->intersect(*Q, intersections, &intersections_len, neighbors, identical);
            for(int i = 0; i < intersections_len; i++){
                nodes_.push_back(Node(intersections[i].value, p, intersections[i].t, q, intersections[i].u, points));
                result.push_back(intersections[i].value);
                if(dbg_file_lvl >= 4) std::cout << "      found: " << nodes_.back().dbg() << '\n';
            }
        }
    }
    //sorting nodes_ to nodes: every node twice, sorted in as {value, p, t, q, u} and {value, q, u, p, t}
    //using p/q as first sorting order and t/u as second
    if(nodes_.size() == 0){
        std::cerr << "no nodes, no intersections, no Chip!" << std::endl;
        nodes = std::vector<Node>();
        return result;
    }
    int count;
    for(int n = 0; n < nodes_.size(); n++){
        count = 0;
        for(int m = 0; m < nodes_.size(); m++){
            if(nodes_[m].p < nodes_[n].p || nodes_[m].p == nodes_[n].p && nodes_[m].t < nodes_[n].t){
                count++;
            }
            if(nodes_[m].q < nodes_[n].p || nodes_[m].q == nodes_[n].p && nodes_[m].u < nodes_[n].t){
                count++;
            }
        }
        nodes_[n].index = count;
        if(dbg_file_lvl >= 3) cout << "node " << n << ": index == " << count << std::endl;
    }
    nodes = std::vector<Node>({nodes_[0]});
    bool inserted;
    Node to_insert = nodes_[0];
    for(int n = 1; n < nodes_.size(); n++){
        inserted = false;
        to_insert = nodes_[n];
        //inserting nodes_[n] once as the original and once in the flipped form
        for(int m = 0; m < nodes.size(); m++){
            if(nodes[m] > to_insert){
                nodes.insert(nodes.begin()+m, to_insert);
                inserted = true;
                break;
            }
        }
        if(!inserted){
            nodes.push_back(to_insert);
        }
    }
    if(dbg_file_lvl > 2){
        std::cout << "nodes:" << '\n';
        for(int n = 0; n < nodes.size(); n++){
            std::cout << "  nodes[" << n << "] == " << nodes[n].dbg() << '\n';
        }
    }
    return result;
}
std::vector<Vector> Chip::intersect_linear(Vector A = Vector(0, 0), Vector a = Vector(1, 1)){
    if(dbg_file_lvl >= 4) std::cout << "intersections of: diagonal (0, 0)->(1, 1) with" << '\n';
    Spline* P;
    Intersection intersections[3];
    int intersections_len;
    std::vector<Node> nodes_ = std::vector<Node>(); //unsorted
    std::vector<Vector> result;
    for(int p = 0; p < points.size(); p++){
        if(dbg_file_lvl >= 4) std::cout << "  spline following " << char(65 + p) << '\n';
        intersections_len = 0;
        P = new Spline(points[p], points[(p+1) % points.size()]);
        P->intersect_linear(A, a, intersections, &intersections_len);
        for(int i = 0; i < intersections_len; i++){
            nodes_.push_back(Node(intersections[i].value, p, -1, intersections[i].t, intersections[i].u));
            result.push_back(intersections[i].value);
            if(dbg_file_lvl >= 4) std::cout << "      found: " << nodes_.back().dbg() << '\n';
        }
    }
    return result;
}
void Chip::follow_edge(int p_curr, float t_curr, int sign, int d, int from, bool SplineConstruct_approximate = true){
    SplineConstruct edge_curr;
    int m;
    // int dbg_file_lvl = 5;
    bool done;
    int m_final;
    float t_final = (sign == +1) ? 1 : 0;
    bool axis_q;
    if(dbg_file_lvl >= 3) std::cout << "follow_edge " << char(80 + from) << char(48 + d) << "(" << char(65 + p_curr) << ": " << t_curr << ", " << (sign > 0 ? '+' : '-') << ")\n";
    while(true){
        done = false;
        m = 0;
        for(
            //going through the nodes forward ? if following edge forwards : else backwards
            m = (sign == +1) ? 0 : (nodes.size()-1);
            (
                sign == +1 && m < nodes.size() ||
                sign == -1 && m >= 0
            );
            m += sign
        ){
            //if on same axis #TODO: describe well
            if(
                p_curr == nodes[m].p && (
                    sign == +1 && t_curr < nodes[m].t && nodes[m].t < t_final || //look for immediate after if moving forward
                    sign == -1 && t_curr > nodes[m].t && nodes[m].t > t_final  //look for immediate before else
                )
            ){
                done = true;
                m_final = m;
                t_final = nodes[m].t;
                axis_q = false;
            }
            if( //if on other axis
                p_curr == nodes[m].q && (
                    sign == +1 && t_curr < nodes[m].u && nodes[m].u < t_final ||
                    sign == -1 && t_curr > nodes[m].u && nodes[m].u > t_final
                )
            ){
                done = true;
                m_final = m;
                t_final = nodes[m].u;
                axis_q = true;
            }
        }
        if(done){
            Spline to_add = Spline(
                points[p_curr],
                points[(p_curr + 1) % points.size()]
            ).subspline(t_curr, t_final);
            edge_curr.add(
                to_add,
                (t_final - t_curr) * sign
            );
            int in = (
                (sign == +1) ? 2 : 0
            ) + (!axis_q + nodes[m_final].flipped) % 2;
            if(SplineConstruct_approximate && edge_curr.splines.size() > 1){
                SplineConstruct edge_curr_copy = edge_curr;
                edge_curr = SplineConstruct();
                edge_curr += edge_curr_copy.approximate();
                if(dbg_file_lvl >= 2){
                    std::cout << "spline approx for " << char(80 + from) << d << " -- " << char(80 + m_final) << in << "\n";
                    if(dbg_file_lvl >= 4){
                        std::cout << "  original: " << edge_curr_copy.dbg("  ") << '\n';
                        std::cout << "  approx:   " << edge_curr.dbg("  ") << '\n';
                    }
                }
            }
            nodes[from].edges[d] = new Edge(
                from, //from node index
                m_final, //to node index
                d, //out on direction d
                in, //into ->to along this direction
                edge_curr
            );
            // if(sign == -1) result->S.flip();
            if(dbg_file_lvl >= 4) std::cout
                << "  edge_curr.add( Spline " << char(65+p_curr) << char(65+(p_curr+1)%points.size())
                << "(" << t_curr << " .. " << t_final << ")); return Edge: " << nodes[from].edges[d]->dbg() << ";"
            << '\n';
            return;
        }else{
            if(dbg_file_lvl >= 4) std::cout
                << "  edge_curr.add( Spline " << char(65+p_curr) << char(65+(p_curr+1)%points.size())
                << "(" << t_curr << " .. " << 1 * (sign == +1) + 0 * (sign == -1) << ")); p: "
                << char(65 + p_curr) <<  "; next p: "
                << char(65 + (p_curr + sign + points.size()) % points.size()) << "; next t: "
                << 0 * (sign == +1) + 1 * (sign == -1) << ";"
            << '\n';
            edge_curr.add(
                Spline(
                    points[p_curr],
                    points[(p_curr + 1) % points.size()]
                ).subspline(t_curr, 1 * (sign == +1) + 0 * (sign == -1)),
                (1 - t_curr) * (sign == +1) + t_curr * (sign == -1)
            );
            p_curr += sign + points.size(); //next along edge is previous if moving backwards
            p_curr %= points.size();
            t_curr = 0 * (sign == +1) + 1 * (sign == -1); //start of the spline if moving forewards, end of the spline if moving backwards
        }
    }
}
void Chip::make_edges(bool SplineConstruct_approximate = true){
    for(int n = 0; n < nodes.size(); n++){
        if(
            nodes[n].flipped
        ){
            if(dbg_file_lvl >= 2) std::cout << "swapped axes for node " << char(80 + n) << '\n';
            //d == 0: following .p forward
            follow_edge(nodes[n].p, nodes[n].t, +1, 0, n, SplineConstruct_approximate);
            //d == 1: following .q forward
            follow_edge(nodes[n].q, nodes[n].u, +1, 1, n, SplineConstruct_approximate);
            //d == 2: following .p backwards
            follow_edge(nodes[n].p, nodes[n].t, -1, 2, n, SplineConstruct_approximate);
            //d == 3: following .q backwards
            follow_edge(nodes[n].q, nodes[n].u, -1, 3, n, SplineConstruct_approximate);
        }else{
            //d == 0: following .q forward
            follow_edge(nodes[n].q, nodes[n].u, +1, 0, n, SplineConstruct_approximate);
            //d == 1: following .p forward
            follow_edge(nodes[n].p, nodes[n].t, +1, 1, n, SplineConstruct_approximate);
            //d == 2: following .q backwards
            follow_edge(nodes[n].q, nodes[n].u, -1, 2, n, SplineConstruct_approximate);
            //d == 3: following .p backwards
            follow_edge(nodes[n].p, nodes[n].t, -1, 3, n, SplineConstruct_approximate);
        }
        if(dbg_file_lvl >= 3) std::cout << '\n';
    }
}
void Chip::make_faces(){
    if(dbg_file_lvl >= 4) std::cout << "Chip/make_faces" << '\n';
    //reset faces and edges_walked
    faces = std::vector<Face>();
    for(int n = 0; n < nodes.size(); n++){
        for(int d = 0; d < 4; d++){
            nodes[n].edges_walked[d] = false;
        }
    }
    //the current node index and current direction index
    int n_curr;
    int d_curr;
    //current face that the edges, that are walked, are added to
    Face face_curr;
    Edge* edge; //temp var for storing ->to and ->in while moving along the edge
    for(int n = 0; n < nodes.size(); n++){
        for(int d = 0; d < 4; d++){
            if(nodes[n].edges_walked[d]){
                continue;
            }
            //start with new, empty face at node n, direction d
            face_curr = Face();
            n_curr = n;
            d_curr = d;
            while(
                //back at the start node and direction
                !(n_curr == n && d_curr == d) ||
                face_curr.size() == 0
            ){
                //add edge to face and mark in edges_walked
                face_curr.push_back(nodes[n_curr].edges[d_curr]);
                nodes[n_curr].edges_walked[d_curr] = true;
                //copy the state of the current edge so changing n_curr is no problem
                edge = nodes[n_curr].edges[d_curr];
                if(dbg_file_lvl >= 3) std::cout << "  " << char(80 + n_curr) << d_curr << " -- " << char(80 + edge->to) << edge->in;
                //follow edge->to next node and into the right direction and then turn that direction
                n_curr = edge->to;
                d_curr = edge->in;
                d_curr = (d_curr + 1) % 4;
                if(dbg_file_lvl >= 3) std::cout << " -> " << char(80 + n_curr) << d_curr << "\n";
            }
            if(dbg_file_lvl >= 3){
                std::cout << "end while: " << char(80 + n_curr) << d_curr << " == " << char(80 + n) << d << "; push: " << face_curr.to_str() << '\n';
            }
            faces.push_back(face_curr);
        }
    }
    //outer and inner faces
    //finding outer face: first intersection with line (0, 0)->(1, 1)
    std::vector<Node> diagonal_intersections;
    std::vector<Intersection> intersections;
    float min_u = 2, min_t;
    int min_n, min_e;
    for(int n = 0; n < nodes.size(); n++){
        for(int e = 2; e < 4; e++){
            intersections = nodes[n].edges[e]->S.intersect_linear(Vector(0, 0), Vector(1, 1));
            for(int i = 0; i < intersections.size(); i++){
                if(intersections[i].u < min_u){
                    min_u = intersections[i].u;
                    min_t = intersections[i].t;
                    min_n = n;
                    min_e = e;
                }
            }
        }
    }
    if(min_u == 2){
        std::cerr << "Chip::make_faces(): there were no intersections with line x == y, please check the spline and, if needed, change the algorithm to search on other lines!" << '\n';
        return;
    }
    if(dbg_file_lvl >= 2) std::cout << "Chip::make_faces(): closest edge: " << nodes[min_n].edges[min_e]->dbg() << '\n';
    //if the detected edge comes from below the straight line and goes above it,
    //    it must be the outside facing edge, the same edge going the other way is part of an inside face then.
    Vector before_intersect = nodes[min_n].edges[min_e]->S(min_t - .02); //Vector along the edge, just before the intersection
    bool outside = before_intersect.y < before_intersect.x; //is below Line x == y

    int f_done;
    bool done = false;
    for(int f = 0; f < faces.size(); f++){
        for(int e = 0; e < faces[f].size(); e++){
            if(faces[f][e] == nodes[min_n].edges[min_e]){
                faces[f].inside = !outside;
                f_done = f;
                done = true;
                break;
            }
        }
        if(done) break;
    }
    if(!done){
        std::cerr << "Chip::make_faces(): first face not found! this is a bug!" << '\n';
        return;
    }

    std::vector<int> all_faces;
    for(int f = 0; f < faces.size(); f++){
        if(f != f_done){ //only undecided go into all_faces
            all_faces.push_back(f);
        }
    }
    std::vector<int> inside_faces = {};
    std::vector<int> outside_faces = {};
    if(outside){
        outside_faces.push_back(f_done);
    }else{
        inside_faces.push_back(f_done);
    }
    // bool done;
    while(all_faces.size() > 0){
        for(int f : all_faces){
            for(int e = 0; e < faces[f].size(); e++){
                done = false;
                for(int f_in : inside_faces){
                    for(int e_in = 0; e_in < faces[f_in].size(); e_in++){
                        if(dbg_file_lvl >= 5) std::cout << "  faces[" << f << "][" << e << "] == faces[" << f_in << "][" << e_in << "]: " << (faces[f][e]->equal(faces[f_in][e_in])) << '\n';
                        if(faces[f][e]->equal(faces[f_in][e_in])){
                            for(int i = all_faces.size() - 1; i >= 0; i--){
                                if(all_faces[i] == f){
                                    all_faces.erase(all_faces.begin() + i);
                                }
                            }
                            outside_faces.push_back(f);
                            faces[f].inside = 0;
                            if(dbg_file_lvl >= 2) std::cout << "added face " << faces[f].to_str() << " to outside_faces because of edge " << faces[f][e]->dbg() << " common with face " << faces[f_in].to_str() << " from inside_faces!" << '\n';
                            done = true;
                            break;
                        }
                    }
                    if(done) break;
                }
                if(done) break;
                for(int f_out : outside_faces){
                    for(int e_out = 0; e_out < faces[f_out].size(); e_out++){
                        if(dbg_file_lvl >= 5) std::cout << "  faces[" << f << "][" << e << "] == faces[" << f_out << "][" << e_out << "]: " << faces[f][e]->dbg() << " == " << faces[f_out][e_out]->dbg() << ": " << (faces[f][e]->equal(faces[f_out][e_out])) << '\n';
                        if(faces[f][e]->equal(faces[f_out][e_out])){
                            for(int i = all_faces.size() - 1; i >= 0; i--){
                                if(all_faces[i] == f){
                                    all_faces.erase(all_faces.begin() + i);
                                }
                            }
                            inside_faces.push_back(f);
                            faces[f].inside = 1;
                            if(dbg_file_lvl >= 2) std::cout << "added face " << faces[f].to_str() << " to inside_faces because of edge " << faces[f][e]->dbg() << " common with face " << faces[f_out].to_str() << " from outside_faces!" << '\n';
                            done = true;
                            break;
                        }
                    }
                    if(done) break;
                }
                if(done) break;
            }
            if(done) break;
        }
    }
}
template <typename CanvasT>
void Chip::draw(CanvasT* C, int samples){
    for(int p = 0; p < points.size(); p++){
        Spline(points[p], points[(p+1) % points.size()]).draw(C, samples);
    }
}
template <typename CanvasT>
void Chip::draw_net(
    CanvasT* C,
    int t_prec,
    int v_prec,
    int face_from,
    int face_to,
    bool draw_E,
    bool draw_F
){
    if(C == NULL){
        C = new BasicCanvas();
    }
    if(face_to == -1){
        face_to = faces.size();
    }
    SplineConstruct PQ, PR;
    Spline QR = Spline(Vector(), Vector(), Vector(), Vector());
    Edge* ePQ,* eQ_,* e_R,* eRP;
    Vector P, Q, R;
    Vector p1, p2, q1, q2, r1, r2;
    Spline* E[t_prec+1];
    Spline* F[v_prec];
    float t, v;
    float dt = 1.0f/t_prec;
    float dv = 1.0f/v_prec;
    for(int f = face_from; f < face_to; f++){
        if(!faces[f].inside) continue;
        if(faces[f].size() == 1){
            if(dbg_file_lvl >= 2) std::cout << "draw net of face " << faces[f].to_str() << '\n';
            P = nodes[faces[f][0]->from].value;
            Q = faces[f][0]->S(.5);
            p1 = faces[f][0]->S.dL(0);
            p2 = faces[f][0]->S.dL(1);
            q1 = faces[f][0]->S.dL(.5);
            q2 = q1.mult_complex(Vector(0, -1)); //turn right 90°
            PQ += Spline(P, Q, p1, q1);
            PR += Spline(P, Q, -p2, -q1);
            if(draw_E){
                if(drawing) BC->setcolor(GREEN);
                E[0] = new Spline(P, P, Vector(0, 0), Vector(0, 0));
                for(int t_ = 1; t_ <= t_prec; t_++){
                    t = dt*t_;
                    E[t_] = new Spline(
                        PQ(t),
                        PR(t),
                        (PQ.dL(t).mult_complex(Vector(0, -1))),
                        (PR.dL(t).mult_complex(Vector(0,  1)))
                    );
                    if(drawing){
                        E[t_]->draw(C);
                        //     for(float v = 0; v < 1; v += dv){
                        //         BC->line((*(E[t_-1]))(v), (*(E[t_]))(v));
                        //     }
                    }
                }
            }
            if(draw_F){
                F[0] = new Spline(P, P, Vector(0, 0), Vector(0, 0));
                if(drawing){
                    BC->setcolor(BLUE);
                    F[0]->draw(C);
                }
                for(int v_ = 1; v_ <= v_prec; v_++){
                    v = dv*v_;
                    F[v_] = new Spline(P, P, p1*v, p2*v);
                    if(drawing){
                        F[v_]->draw(C);
                        //     for(float v = 0; v < 1; v += dv){
                        //         BC->line((*(E[t_-1]))(v), (*(E[t_]))(v));
                        //     }
                    }
                }
            }
            continue;
        }else if(faces[f].size() == 2){
            if(dbg_file_lvl >= 2) std::cout << "draw net of face " << faces[f].to_str() << '\n';
            P = nodes[faces[f][0]->from].value;
            Q = nodes[faces[f][1]->from].value;
            p1 = faces[f][0]->S.dL(0);
            p2 = faces[f][1]->S.dL(1);
            q1 = faces[f][1]->S.dL(0);
            q2 = faces[f][0]->S.dL(1);
            PQ = faces[f][0]->S;
            PR = faces[f][1]->S;
            if(draw_E){
                if(drawing) BC->setcolor(GREEN);
                E[0] = new Spline(P, P, Vector(0, 0), Vector(0, 0));
                for(int t_ = 1; t_ <= t_prec; t_++){
                    t = dt*t_;
                    E[t_] = new Spline(
                        PQ(t),
                        PR(t),
                        (PQ.dL(t).mult_complex(Vector(0, -1))),
                        (PR.dL(t).mult_complex(Vector(0,  1)))
                    );
                    if(drawing){
                        E[t_]->draw(C);
                        //     for(float v = 0; v < 1; v += dv){
                        //         BC->line((*(E[t_-1]))(v), (*(E[t_]))(v));
                        //     }
                    }
                }
            }
            if(draw_F){
                F[0] = new Spline(P, Q, p1, q2);
                if(drawing){
                    BC->setcolor(BLUE);
                    F[0]->draw(C);
                }
                for(int v_ = 1; v_ <= v_prec; v_++){
                    v = dv*v_;
                    F[v_] = new Spline(P, Q, p1*(1.0f-v) - p2*v, q2*(1.0f-v) - q1*v);
                    if(drawing){
                        F[v_]->draw(C);
                        //     for(float v = 0; v < 1; v += dv){
                        //         BC->line((*(E[t_-1]))(v), (*(E[t_]))(v));
                        //     }
                    }
                }
            }
            continue;
        }
        if(dbg_file_lvl >= 2) std::cout << "draw net of face " << faces[f].to_str() << '\n';
        for(int e = 0; e < faces[f].size(); e++){
            std::cout << "  draw net on edge " << char(80 + faces[f][e]->from) << '\n';
            //edge from R to P
            ePQ = faces[f][e];
            eQ_ = faces[f][(e + 2*faces[f].size() + 1) % faces[f].size()];
            e_R = faces[f][(e + 2*faces[f].size() - 2) % faces[f].size()];
            eRP = faces[f][(e + 2*faces[f].size() - 1) % faces[f].size()];
            P = nodes[faces[f][e]->from].value;
            Q = nodes[faces[f][e]->to  ].value;
            R = nodes[    eRP    ->from].value;
            p1 = ePQ->S.dL(0);
            p2 = eRP->S.dL(1);
            q1 = eQ_->S.dL(0);
            q2 = ePQ->S.dL(1);
            r1 = eRP->S.dL(0);
            r2 = e_R->S.dL(1);
            if(dbg_file_lvl >= 5){
                std::cout << "  P  == " << P << '\n';
                std::cout << "  Q  == " << Q << '\n';
                std::cout << "  R  == " << R << '\n';
                std::cout << "  p1 == " << p1 << '\n';
                std::cout << "  p2 == " << p2 << '\n';
                std::cout << "  q1 == " << q1 << '\n';
                std::cout << "  q2 == " << q2 << '\n';
                std::cout << "  r1 == " << r1 << '\n';
                std::cout << "  r2 == " << r2 << '\n';
            }
            PQ = faces[f][e]->S;
            PR = eRP->S;
            PR.flip();
            QR = Spline(Q, R, q1, r2);
            //dbg
            if(dbg_file_lvl >= 5){
                std::cout << "  A(t) == " << PQ.dbg() << '\n';
                std::cout << "  B(t) == " << PR.dbg() << '\n';
                std::cout << "  C(t) == " << QR.dbg() << '\n';
            }
            if(draw_F){
                F[0] = new Spline(P, QR(0), p1, q2);
                if(drawing){
                    BC->setcolor(BLUE);
                    F[0]->draw(C);
                }
                for(int v_ = 1; v_ <= v_prec; v_++){
                    v = dv*v_;
                    F[v_] = new Spline(
                        P,
                        QR(v),
                        -p2*v + p1*(1.0f-v),
                        -r1*v + q2*(1.0f-v)
                    );
                    if(drawing){
                        F[v_]->draw(C);
                        // for(float t = 0; t < 1; t += dt){
                        //     BC->line((*(F[v_-1]))(t), (*(F[v_]))(t));
                        // }
                    }
                }
            }
            if(draw_E){
                if(drawing) BC->setcolor(GREEN);
                E[0] = new Spline(P, P, Vector(0, 0), Vector(0, 0));
                for(int t_ = 1; t_ <= t_prec; t_++){
                    t = dt*t_;
                    E[t_] = new Spline(
                        PQ(t),
                        PR(t),
                        (q1*t - p2*(1.0f-t))*t,
                        (r2*t - p1*(1.0f-t))*t
                    );
                    if(drawing){
                        E[t_]->draw(C);
                    //     for(float v = 0; v < 1; v += dv){
                    //         BC->line((*(E[t_-1]))(v), (*(E[t_]))(v));
                    //     }
                    }
                }
            }
        }
    }
}
template <typename CanvasT>
void Chip::mark_points(CanvasT* C){
    for(int p = 0; p < points.size(); p++){
        _plus(C, points[p].V);
    }
}
void Chip::transform(float a, float b, float c, float d, float e, float f){
    for(int p = 0; p < points.size(); p++){
        points[p].transform(a, b, c, d, e, f);
    }
}
std::string Chip::latex(){
    stringstream result;
    Point* point_next;
    for(int p = 0; p < points.size(); p++){
        result << "        \\node ("
            << (char)(65+p) << ") at ("
            << points[ p ].V.x << ", "
            << points[ p ].V.y << ") {$\\bullet$};\n"
        ;
    }
    for(int p = 0; p < points.size(); p++){
        point_next = &(points[(p+1) % points.size()]);
        result << "        \\draw ("
            << (char)(65+p) << ") node[below] {$"
            << (char)(65+p) << "$} .. controls +({"
            <<  points[ p ].v.x << "/3}, {" <<  points[ p ].v.y << "/3}) and +({"
            << -point_next->v.x << "/3}, {" << -point_next->v.y << "/3}) .. ("
            << (char)(65 + (p+1) % points.size()) << ");\n"
        ;
    }
    return result.str();
}
std::string Chip::dbg(std::string indent = ""){
    stringstream result;
    result << indent << "C: \n";
    result << indent << "  points: \n";
    for(int p = 0; p < points.size(); p++){
        result << indent << "    " << points[p].dbg() << "\n";
    }
    result << indent << "  nodes: \n";
    for(int n = 0; n < nodes.size(); n++){
        result << indent << "    " << nodes[n].dbg() << "\n";
        if(nodes[n].edges[0] == NULL){
            continue;
        }
        result << indent << "    edges: \n";
        for(int d = 0; d < 4; d++){
            result << indent << "      " << nodes[n].edges[d]->dbg() << "\n";
        }
    }
    result << indent << "  faces: \n";
    for(int f = 0; f < faces.size(); f++){
        result << indent << "    F(";
        for(int e = 0; e < faces[f].size(); e++){
            result << char(80 + faces[f][e]->from);
        }
        result << ") " << (faces[f].inside ? "inside" : "outside");
        result << "\n";
    }
    return result.str();
}

#endif /* end of include guard: CHIP_H */
