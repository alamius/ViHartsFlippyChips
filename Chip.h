#ifndef CHIP_H
#define CHIP_H

class Edge{
public:
    int from, to; //nodes indices
    int out, in; //the numbers of the connections in the nodes from and to.
    SplineConstruct S;
    Edge(){};
    Edge(int from_, int to_, int out_, int in_, SplineConstruct S_){
        from = from_;
        to = to_;
        out = out_;
        in = in_;
        S = S_;
    };
    bool operator==(Edge);
    bool equal(Edge*);
    std::string dbg(bool spline_dbg);
    virtual ~Edge(){};
};
std::string Edge::dbg(bool spline_dbg=false){
    stringstream result;
    result << "E(" << char(80 + from) << out << " -- " << char(80 + to) << in << ")";
    if(spline_dbg) result << ": " << S.dbg();
    return result.str();
}
bool Edge::operator==(Edge other){
    return (
        from == other.from &&
        to   == other.to   &&
        in   == other.in   &&
        out  == other.out
    ) || (
        from == other.to   &&
        to   == other.from &&
        in   == other.out  &&
        out  == other.in
    );
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

class Node{
public:
    Vector value; //the vector of the intersection
    int p; //index of the first spline's begin point in the line (Chip::points[p])
    int q; //index of the second spline's begin point  (Chip::points[q])
    float t; //parameter of the intersection along the first spline (value = Spline(Chip::points[p] ~~ Chip::points[p+1])(t))
    float u; //parameter along the second spline
    Edge* edges[4];
    bool edges_walked[4];
    bool flipped;
    int index = -1;
    Node(){
        for(int d = 0; d < 4; d++){
            edges[d] = new Edge();
        }
    };
    Node(int index_, Vector value_, int p_, float t_, int q_, float u_, std::vector<Point>& points){
        index = index_;
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
    };
    Node(int index_, Vector value_, int p_, float t_, int q_, float u_){
        index = index_;
        value = value_;
        p = p_;
        t = t_;
        q = q_;
        u = u_;
    };
    Node(Vector value_, int n_, int e_, float t_, float u_){ //only for intersections Spline - straight Line (Spline/SplineConstruct::intersect_linear)
        value = value_;
        p = n_;
        q = e_;
        t = t_;
        u = u_;
    };
    std::string dbg();
    bool operator>(Node other);
    bool operator<(Node other);
    virtual ~Node(){};
};
std::string Node::dbg(){
    stringstream ss;
    ss  << "N("
        // << value << "; {"
        << char(65 + p) << ": " << t << " x "
        << char(65 + q) << ": " << u
    // << "}";
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

class Face{
public:
    std::vector<Edge*> edges;
    int inside = -1; //not decided yet. 0 for outside, 1 for inside
    Face(std::vector<Edge*> edges_ = std::vector<Edge*>()){
        for(int e = 0; e < edges_.size(); e++){
            edges.push_back(edges_[e]);
        }
    }
    Edge* operator[](int);
    int size();
    void push_back(Edge*);
    std::string to_str();
    std::string dbg(string indent);
    virtual ~Face(){};
};
Edge* Face::operator[](int i){
    return edges[i];
}
int Face::size(){
    return edges.size();
}
void Face::push_back(Edge* e){
    edges.push_back(e);
}
std::string Face::to_str(){
    std::string result = "(";
    for(int e = 0; e < size(); e++){
        result += char(80 + edges[e]->from);
    }
    result += ")";
    return result;
}
std::string Face::dbg(string indent = ""){
    stringstream result;
    result << indent << "F({\n";
    for(int e = 0; e < size(); e++){
        result << indent << "  " << edges[e]->dbg() << ", \n";
    }
    result << indent << "})";
    return result.str();
}

colorint* std_color_func(colorint result[COLOR_LEN], float t, float v){
    make_color(
        result,
        255.0f*v,
        //255.0f*(t*t*t*(1.0f-v) + v),
        255.0f*(1.0f - t)*(1.0f - t)
        //255.0f*(1.0f-.99*(t*t*t*(1.0f-v) + v))
        // 255.0f*pow(.01, t)
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

class Chip{
private:
    //points (Vector of position and Vector of direction combined, defined in Spline.h) of the Chip. the line of the chip is the combination of the splines from every point to the next
    std::vector<Point> points;
    //will be created at the intersections of the line defined by points, containing the intersection Vector and the indices of the points and correcponding parameters (see Nodes)
    std::vector<Node> nodes;
    //the faces of the chip with their edges
    std::vector<Face> faces;
public:
    //only copies the points
    Chip(std::vector<Point> points_);
    //
    void color(int, int, bool, bool, int, colorint* (*)(colorint color[COLOR_LEN], float t, float v), bool apply_gauss);
    std::vector<Vector> intersect();
    std::vector<Vector> intersect_linear(Vector A, Vector a);
    void follow_edge(int p_curr, float t_curr, int sign, int d, int from, bool SplineConstruct_approximate);
    void make_edges(bool SplineConstruct_approximate);
    void test_edges();
    void make_faces();
    template <typename CanvasT>
    void draw(CanvasT* C, int samples = 30);
    template <typename CanvasT>
    void draw_net(
        CanvasT* C, // = (CanvasT*)NULL,
        int t_prec = 5,
        int v_prec = 5,
        int face_from = 0,
        int face_to = -1, // -1 is turned into faces.size() - 1
        bool draw_E = true,
        bool draw_F = false
    );
    template <typename CanvasT>
    void mark_points(CanvasT* C);
    void transform(float a, float b, float c, float d, float e, float f); // [x, y][[a, b], [c, d]] + [e, f]
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
    int t_prec = 30,
    int v_prec = 30,
    bool draw_E = false,
    bool draw_F = true,
    int pensize = 1,
    colorint* (*color_func)(colorint result[COLOR_LEN], float t, float v) = &std_color_func,
    bool apply_gauss = true
){
    intersect();
    make_edges(true);
    make_faces();
    LC = new LayeredCanvas();
    LC->setpensize(pensize);
    int max_edges = 0;
    for(int f = 0; f < faces.size(); f++){
        if(faces[f].size() > max_edges && faces[f].inside == 1){
            max_edges = faces[f].size();
        }
    }
    const int LC_len = max_edges;
    //variables for face coloring
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
    colorint color[COLOR_LEN];
    for(int e_LC = 0; e_LC < LC_len; e_LC++){
        for(int f = 0; f < faces.size(); f++){
            if(!faces[f].inside) continue;
            if(faces[f].size() == 1){
                if(e_LC != 0) continue; //only single loops go to layer 0
                if(dbg_file_lvl >= 2) std::cout << "color face " << faces[f].to_str() << " edge " << char(80 + faces[f][0]->from) << char(80 + faces[f][0]->to) << '\n';
                P = nodes[faces[f][0]->from].value;
                Q = faces[f][0]->S(.5);
                p1 = faces[f][0]->S.dL(0);
                p2 = faces[f][0]->S.dL(1);
                q1 = faces[f][0]->S.dL(.5);
                q2 = q1.mult_complex(Vector(0, -1)); //turn right 90°
                PQ += Spline(P, Q, p1, q1);
                PR += Spline(P, Q, -p2, -q1);
                if(draw_E){
                    // LC->setcolor(GREEN);
                    E[0] = new Spline(P, P, Vector(0, 0), Vector(0, 0));
                    for(int t_ = 1; t_ <= t_prec; t_++){
                        t = dt*t_;
                        E[t_] = new Spline(
                            PQ(t),
                            PR(t),
                            (PQ.dL(t).mult_complex(Vector(0, -1))),
                            (PR.dL(t).mult_complex(Vector(0,  1)))
                        );
                        // E[t_]->draw(LC);
                        for(float v = 0; v < 1; v += dv){
                            LC->setcolor(color_func(color, t, v));
                            LC->quadrilateral_unchecked(
                                (*(E[t_-1]))(v+dv),
                                (*(E[t_-1]))(v),
                                (*(E[t_  ]))(v+dv),
                                (*(E[t_  ]))(v)
                            );
                        }
                    }
                }
                if(draw_F){
                    F[0] = new Spline(P, P, Vector(0, 0), Vector(0, 0));
                    // LC->setcolor(BLUE);
                    // F[0]->draw(LC);
                    for(int v_ = 1; v_ <= v_prec; v_++){
                        v = dv*v_;
                        if(0 && nodes[faces[f][0]->from].index % 2 == 0){
                            v = 1.0f - v;
                        }
                        F[v_] = new Spline(P, P, p1*v, p2*v);
                        // F[v_]->draw(LC);
                        for(float t = 0; t < 1; t += dt){
                            LC->setcolor(color_func(color, t, t));
                            LC->quadrilateral_unchecked(
                                (*(F[v_-1]))(t+dt),
                                (*(F[v_-1]))(t),
                                (*(F[v_  ]))(t+dt),
                                (*(F[v_  ]))(t)
                            );
                        }
                    }
                }
                if(dbg_file_lvl >= 4)
                    LC->write(
                        filename+"_f"+_to_str(f)+"="+faces[f].to_str()+"_e"+char(80+faces[f][0]->from)+char(80+faces[f][0]->to)
                    );
                continue;
            }else if(faces[f].size() == 2){
                if(dbg_file_lvl >= 2) std::cout << "color face " << faces[f].to_str() << '\n';
                for(int e = 0; e < faces[f].size(); e++){
                    if(e != e_LC) continue; //only edges of the same index as the current canvas layer are drawn
                    if(dbg_file_lvl >= 2) std::cout << "  color by edge " << char(80 + faces[f][e]->from) << char(80 + faces[f][e]->to) << '\n';
                    int g = (e + 1) % 2;
                    if(1){
                        P = nodes[faces[f][e]->from].value;
                        Q = nodes[faces[f][g]->from].value;
                        p1 = faces[f][e]->S.dL(0);
                        p2 = faces[f][g]->S.dL(1);
                        q1 = faces[f][g]->S.dL(0);
                        q2 = faces[f][e]->S.dL(1);
                        PQ = faces[f][e]->S;
                        PR = faces[f][g]->S;
                    }
                    if(draw_E){
                        // LC->setcolor(GREEN);
                        E[0] = new Spline(P, P, Vector(0, 0), Vector(0, 0));
                        for(int t_ = 1; t_ <= t_prec; t_++){
                            t = dt*t_;
                            E[t_] = new Spline(
                                PQ(t),
                                PR(t),
                                (PQ.dL(t).mult_complex(Vector(0, -1))),
                                (PR.dL(t).mult_complex(Vector(0,  1)))
                            );
                            // E[t_]->draw(LC);
                            for(float v = 0; v < 1; v += dv){
                                LC->setcolor(color_func(color, t, v));
                                LC->quadrilateral_unchecked(
                                    (*(E[t_-1]))(v+dv),
                                    (*(E[t_-1]))(v),
                                    (*(E[t_  ]))(v+dv),
                                    (*(E[t_  ]))(v)
                                );
                            }
                        }
                    }
                    if(draw_F){
                        F[0] = new Spline(P, Q, p1, q2);
                        // LC->setcolor(BLUE);
                        // F[0]->draw(LC);
                        for(int v_ = 1; v_ <= v_prec; v_++){
                            v = dv*v_;
                            if(0 && nodes[faces[f][e]->from].index % 2 == 0){
                                v = 1.0f - v;
                            }
                            F[v_] = new Spline(P, Q, p1*(1.0f-v) - p2*v, q2*(1.0f-v) - q1*v);
                            // F[v_]->draw(LC);
                            for(float t = 0; t < .99; t += dt){
                                LC->setcolor(color_func(color, t, v));
                                LC->quadrilateral_unchecked(
                                    (*(F[v_-1]))(t+dt),
                                    (*(F[v_-1]))(t),
                                    (*(F[v_  ]))(t+dt),
                                    (*(F[v_  ]))(t)
                                );
                            }
                        }
                    }
                    if(dbg_file_lvl >= 4)
                        LC->write(
                            filename+"_f"+_to_str(f)+"="+faces[f].to_str()+"_e"+char(80+faces[f][e]->from)+char(80+faces[f][e]->to)
                        );
                }
                continue;
            }else{
                if(dbg_file_lvl >= 2) std::cout << "color face " << faces[f].to_str() << '\n';
                for(int e = 0; e < faces[f].size(); e++){
                    if(e != e_LC) continue; //only edges of the same index as the current canvas layer are drawn
                    if(dbg_file_lvl >= 2) std::cout << "  color by edge " << char(80 + faces[f][e]->from) << char(80 + faces[f][e]->to) << '\n';
                    //edge from R to P
                    if(1){
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
                        PQ = faces[f][e]->S;
                        PR = eRP->S;
                        PR.flip();
                        QR = Spline(Q, R, q1, r2);
                    }
                    if(draw_E){
                        // LC->setcolor(GREEN);
                        E[0] = new Spline(P, P, Vector(0, 0), Vector(0, 0));
                        for(int t_ = 1; t_ <= t_prec; t_++){
                            t = dt*t_;
                            E[t_] = new Spline(
                                PQ(t),
                                PR(t),
                                (q1*t - p2*(1.0f-t))*t,
                                (r2*t - p1*(1.0f-t))*t
                            );
                            // if(drawing){
                            // E[t_]->draw(LC);
                            for(float v = 0; v < 1; v += dv){
                                LC->setcolor(color_func(color, t, v));
                                LC->quadrilateral_unchecked(
                                    (*(E[t_-1]))(v+dv),
                                    (*(E[t_-1]))(v),
                                    (*(E[t_  ]))(v+dv),
                                    (*(E[t_  ]))(v)
                                );
                            }
                        }
                    }
                    if(draw_F){
                        F[0] = new Spline(P, QR(0), p1, q2);
                        // LC->setcolor(BLUE);
                        // F[0]->draw(LC);
                        for(int v_ = 1; v_ <= v_prec; v_++){
                            v = dv*v_;
                            if(0 && nodes[faces[f][e]->from].index % 2 == 0){
                                v = 1.0f - v;
                            }
                            F[v_] = new Spline(
                                P,
                                QR(v),
                                -p2*v + p1*(1.0f-v),
                                -r1*v + q2*(1.0f-v)
                            );
                            // F[v_]->draw(LC);
                            for(float t = 0; t < .99; t += dt){
                                LC->setcolor(color_func(color, t, v));
                                LC->quadrilateral_unchecked(
                                    (*(F[v_-1]))(t+dt),
                                    (*(F[v_-1]))(t),
                                    (*(F[v_  ]))(t+dt),
                                    (*(F[v_  ]))(t)
                                );
                                // LC->line((*(F[v_-1]))(t), (*(F[v_]))(t));
                            }
                        }
                    }
                    if(dbg_file_lvl >= 4)
                        LC->write(filename+"_f"+_to_str(f)+"="+faces[f].to_str()+"_e"+char(80+faces[f][e]->from)+char(80+faces[f][e]->to));
                }
            }
        }
        if(dbg_file_lvl >= 3){
            std::cout << "written layer " << e_LC << " of " << LC_len << "; size == " << sizeof(BasicCanvas) << '\n';
            LC->write(filename+"_layer"+_to_str(e_LC), *write_bg_color);
            if(dbg_file_lvl >= 4) LC->dump(filename+"_layer"+_to_str(e_LC));
        }
        LC->save(filename+"_layer"+_to_str(e_LC));
        LC->clear();
    }
    LC->clear();
    if(dbg_file_lvl >= 2) std::cout << "adding image layers together:" << '\n';
    if(dbg_file_lvl >= 4) LC->dump(filename+"0_bg");
    for(int bc = 0; bc < LC_len; bc++){
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
                nodes_.push_back(Node(-1, intersections[i].value, p, intersections[i].t, q, intersections[i].u, points));
                result.push_back(intersections[i].value);
                if(dbg_file_lvl >= 4) std::cout << "      found: " << nodes_.back().dbg() << '\n';
            }
            // getch();
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
        // getch();
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
        // getch();
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
            //going through the nodes forward if following edge forwards else backwards
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
            // if(sign == -1) edge_curr.splines.back().flip();
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
        if(dbg_file_lvl >= 2) std::cout << '\n';
    }
}
void Chip::test_edges(){
    if(dbg_file_lvl >= 2) std::cout << "test_edges" << '\n';
    Vector V, v;
    Vector rotate_left = Vector(0, 1); //read as complex number i
    for(int n = 0; n < nodes.size(); n++){
        if(dbg_file_lvl >= 3) std::cout << "  node " << char(80 + n) << '\n';
        for(int d = 0; d < 4; d++){
            // int d = 3;
            if(dbg_file_lvl >= 3) std::cout << "    edge " << d << ": " << nodes[n].edges[d]->dbg() << ": \n" << nodes[n].edges[d]->S.dbg("      ") << "\n";
            if(drawing)
                for(float t = 0; t <= 1; t += .1){
                    V = nodes[n].edges[d]->S(t);
                    v = nodes[n].edges[d]->S.dL(t).norm(.05 * (1.0f-t));
                    BC->line(V, V + v.mult_complex(rotate_left));
                }
        }
        if(dbg_file_lvl >= 2) std::cout << '\n';
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
    //finding outer face: most edges
    // TODO !!!
    // int max_edges = 0, max_f = 0;
    // for(int f = 1; f < faces.size(); f++){
    //     if(faces[f].size() > max_edges){
    //         max_f = f;
    //         max_edges = faces[f].size();
    //     }
    // }
    // if(dbg_file_lvl >= 2) std::cout << "chosen face " << faces[max_f].to_str() << " as outer face." << '\n';
    // faces[max_f].inside = 0;
    // int f_done = max_f;
    // bool done = false;
    //
    //finding outer face: first intersection with line (0, 0)->(1, 1)
    std::vector<Node> diagonal_intersections;
    std::vector<Intersection> intersections;
    float min_u = 2;
    int min_n, min_e;
    for(int n = 0; n < nodes.size(); n++){
        for(int e = 2; e < 4; e++){
            // if( //edge goes forward
            //     (nodes[n].edges[e].p + 1) % points.size() ==
            //     nodes[n].edges[e].p
            // ) continue;
            intersections = nodes[n].edges[e]->S.intersect_linear(Vector(0, 0), Vector(1, 1));
            for(int i = 0; i < intersections.size(); i++){
                // diagonal_intersections.push_back(Node(
                //     intersections[i].value,
                //     n, e,
                //     intersections[i].t,
                //     intersections[i].u
                // ));
                if(intersections[i].u < min_u){
                    min_u = intersections[i].u;
                    min_n = n;
                    min_e = e;
                }
            }
        }
    }
    if(min_u == 2){
        std::cerr << "there were no intersections with line x == y, please check the spline and, if needed, change the algorithm to search on other lines!" << '\n';
        return;
    }
    int f_done;
    bool done = false;
    for(int f = 0; f < faces.size(); f++){
        for(int e = 0; e < faces[f].size(); e++){
            if(faces[f][e] == nodes[min_n].edges[min_e]){
                faces[f].inside = 0;
                f_done = f;
                done = true;
                break;
            }
        }
        if(done) break;
    }
    if(!done){
        std::cerr << "first face not found again!" << '\n';
        return;
    }

    std::vector<int> all_faces;
    for(int f = 0; f < faces.size(); f++){
        if(f != f_done){ //not decided
            all_faces.push_back(f);
        }
    }
    std::vector<int> inside_faces = {};
    std::vector<int> outside_faces = {f_done};
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
            // getch();
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
