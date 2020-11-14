#ifndef CHIPS_TESTS_CPP
#define CHIPS_TESTS_CPP

void test_snake(){
    BC = new BasicCanvas();
    Spline A = Spline(Point(Vector(.2, .2), Vector(0.5, 2)), Point(Vector(.8, .8), Vector(.5, 2)));
    float step = .01;
    float dt = .5;
    A.draw(BC);
    Spline B = A;
    for(float t = 0; t <= 1.0f-dt; t += step){
        BC->setcolor(BLUE);
        B = A.subspline(t, t+dt);
        B.transform(1, 0, 0, 1, 0, .1);
        B.draw(BC);
    }
    BC->write("test_snake");
    delete BC;
}
void test_Splines(){
    BC = new BasicCanvas();
    Spline A = Spline(Point(Vector(.3, .3), Vector(0, 2)), Point(Vector(.7, .3), Vector(0, -2)));
    std::cout << "A:" << A.dbg() << '\n';
    if(drawing) A.draw(BC);
    A.transform(1, 0, 0, 1, 0, .1);
    A.flip();
    std::cout << "A:" << A.dbg() << '\n';
    if(drawing) A.draw(BC);
    // return;
    Spline B = Spline(Point(Vector(.3, .7), Vector(3,-2)), Point(Vector(.7, .7), Vector(3, 2)));
    if(drawing){
        A.draw(BC);
        B.draw(BC);
    }
    std::vector<Intersection> intersections;
    Intersection I[3];
    int I_len = 0;
    A.intersect(B, I, &I_len);
    std::cout << "intersections A x B:" << '\n';
    for(int i = 0; i < I_len; i++){
        intersections.push_back(I[i]);
        std::cout << "  " << intersections.back().value << '\n';
        std::cout << "  A(" << I[i].t << ") == " << A(I[i].t) << '\n';
        std::cout << "  B(" << I[i].u << ") == " << B(I[i].u) << '\n';
    }
    intersections = std::vector<Intersection>();
    I_len = 0;
    B.intersect(B, I, &I_len, false, true);
    std::cout << "intersections B x B:" << '\n';
    for(int i = 0; i < I_len; i++){
        intersections.push_back(I[i]);
        std::cout << "  " << intersections.back().value << '\n';
        std::cout << "  B(" << I[i].t << ") == " << B(I[i].t) << '\n';
        std::cout << "  B(" << I[i].u << ") == " << B(I[i].u) << '\n';
    }
    BC->setcolor(BLUE);
    float t1 = intersections[0].t;
    float t2 = intersections[0].u;
    Spline C = B.subspline(intersections[0].t, intersections[0].u);
    // C.transform(.5, 0, 0, .5, .2, .2);
    if(drawing){
        C.draw(BC);
    }
    BC->write("test_Splines");
    delete BC;
    std::cout << "test_Splines finished" << '\n';
}
void test_SplineConstructs(){
    BC = new BasicCanvas();
    Spline A = Spline(Point(Vector(.5, .6), Vector(-.1, .5)), Point(Vector(.25, .5), Vector( 0, -.5)));
    Spline B = Spline(Point(Vector(.25, .5), Vector( 0, -.5)), Point(Vector(.5, .1), Vector( .1, -.5)));
    SplineConstruct S = SplineConstruct({A, B}, {1, 2});
    std::cout << S.dbg() << '\n';
    for(float t = 0; t <= 1.01; t += .1){
        // BC->line(T(t), T(t) + Vector(.02, 0));
        std::cout << "S(" << t << ") == " << S(t).to_str("(%.2f, %.2f)") << '\n';
    }
    S.flip();
    std::cout << S.dbg() << '\n';
    SplineConstruct T = SplineConstruct({A.subspline(.5, 1), B.subspline(0, .7)}, {.5, .7});
    A.transform(-1, 0, 0, 1, 1, 0);
    B.transform(-1, 0, 0, 1, 1, 0);
    T.transform(1, 0, 0, 1, .1, 0);
    std::cout << T.dbg() << '\n';
    if(drawing){
        A.draw(BC);
        B.draw(BC);
        S.draw(BC);
        T.draw(BC);
    }
    for(float t = 0; t <= 1.01; t += .1){
        BC->line(T(t), T(t) + Vector(.02, 0));
        std::cout << "S.dbg_L(" << t << ") == " << S.dbg_L(t) << S(t).to_str("(%.2f, %.2f)") << '\n';
    }
    Vector a, b;
    for(float t = 0; t < 1; t += .1){
        if(t < .5){
            a = A(2*t);
        }else{
            a = B(2*t - 1);
        }
        b = S(t);
        BC->line(T(t), T(t) + Vector(.02, 0));
    }
    BC->write("test_SplineConstructs");
    delete BC;
    std::cout << "test_SplineConstructs finished" << '\n';
}
void test_SplineConstruct_approximate(){
    BC = new BasicCanvas();
    Spline A = Spline(Point(Vector(.5, .6), Vector(-.1, .5)), Point(Vector(.25, .5), Vector( 0, -.5)));
    Spline B = Spline(Point(Vector(.25, .5), Vector( 0, -.5)), Point(Vector(.5, .1), Vector( .1, -.5)));
    SplineConstruct S = SplineConstruct({A, B}, {1, 2});
    std::cout << S.dbg() << '\n';
    Spline AB1 = S.approximate(300, 0.01, 1);
    Spline AB2 = S.approximate(1000, 0.01, 1);
    std::cout << "S.approx(300): " << AB1.dbg() << '\n';
    std::cout << "S.approx(1000): " << AB2.dbg() << '\n';
    SplineConstruct AB1SC;
    AB1SC += AB1;
    std::cout << "dist: " << AB1SC.dist(AB2) << '\n';
    // A.transform(-1, 0, 0, 1, 1, 0);
    // B.transform(-1, 0, 0, 1, 1, 0);
    if(drawing){
        // A.draw(BC);
        // B.draw(BC);
        S.draw(BC);
        AB1.draw(BC);
        AB2.draw(BC);
    }
    BC->write("test_SplineConstruct_approximate");
    delete BC;
    std::cout << "test_SplineConstruct_approximate finished" << '\n';
}
void test_basic_canvas(){ //this will not look proper, the requested behavior is only functional in LayeredCanvas
    BC = new BasicCanvas();
    Vector a = Vector(.2, .2);
    Vector b = Vector(.5, .1);
    Vector c = Vector(.3, .9);
    Vector d = Vector(.9, .5);
    BC->setpensize(5);
    BC->line(a, b);
    BC->setcolor(0, 255, 0, 255);
    BC->line_variation(a+Vector(0, .1), b+Vector(0, .1));
    BC->write("test_basic_canvas");
    delete BC;
    std::cout << "test_basic_canvas finished" << '\n';
}
void test_BC_LC(){
    LC = new LayeredCanvas();
    Vector a = Vector(.2, .2);
    Vector b = Vector(.5, .1);
    Vector c = Vector(.3, .9);
    Vector d = Vector(.9, .5);
    // BC->setcolor(GREEN);
    // BC->quadrilateral_unchecked(a, b, c, d);
    Vector ab0, ab1, cd0, cd1;
    colorint color[COLOR_LEN];
    LC->setpensize(5);
    for(float t = 0; t < 1; t += .2){
        ab0 = a*(1.0f-t) + b*t;
        ab1 = a*(0.8f-t) + b*(t+.2);
        cd0 = c*(1.0f-t) + d*t;
        cd1 = c*(0.8f-t) + d*(t+.2);
        for(float v = 0; v < 1; v += .2){
            LC->setcolor(std_color_func(color, t, v));
            LC->quadrilateral_unchecked(
                ab0*(1.0f-v) + cd0*v,
                ab1*(1.0f-v) + cd1*v,
                ab0*(0.8f-v) + cd0*(v+.2),
                ab1*(0.8f-v) + cd1*(v+.2)
            );
        }
    }
    LC->quadrilateral_unchecked(Vector(0.23812, 0.64456), Vector(0.23465, 0.70827), Vector(0.23819, 0.64664), Vector(0.23465, 0.70827));
    // LC->setcolor(100, 0, 0, 127);
    // LC->quadrilateral_unchecked(a, b, c, d);
    // LC->setcolor(200, 0, 0, 127);
    // LC->quadrilateral_unchecked(a*.5, b*.5, c*.5, d*.5);
    LC->save("test");
    LC->dump("test0");
    LC->write("test_BC_LC_1");
    delete LC;
    BC = new BasicCanvas();
    // BC->load("test");
    // BC->dump("test1");
    // BC->write("test_BC_LC_2");
    BC->background(0, 60, 0, 200);
    BC->dump("test2");
    BC->load("test");
    BC->dump("test3");
    BC->write("test_BC_LC_3");
    delete BC;
    std::cout << "test_BC_LC finished" << '\n';
}
void test_knot_1(bool print = false){
    Point A = Point(Vector(.45, .7), Vector(-1, 0));
    Point B = Point(Vector(.15, .9), Vector(.5, .5));
    Point C = Point(Vector(.45, .1), Vector( 1, 0));
    Point D = Point(Vector(.75, .9), Vector(-1, 0));
    Point E = Point(Vector(.75, .1), Vector( 1, 0));

    Chip chip = Chip({A, B, C, D, E});
    // std::vector<Vector> intersections;
    // if(print) std::cout << "chip.latex():\n" << chip.latex() << '\n';
    // if(print){
    //     intersections = chip.intersect();
    //     for(int i = 0; i < intersections.size(); i++){
    //         std::cout << "        \\node ("
    //             << char(80+i) << ") at "
    //             << intersections[i] << " {$\\bullet$};\n";
    //         std::cout << "        \\draw ("
    //             << char(80+i) << ") node[below] {$"
    //             << char(80+i) << "$};\n";
    //     }
    // }
    // chip.transform(.9, 0, 0, .9, .2, 0);

    // if(drawing){
    //     BC->clear();
    //     chip.draw(BC);
    //     chip.mark_points(BC);
    // }
    // intersections = chip.intersect();
    // if(dbg_file_lvl >= 3){
    //     if(drawing) BC->setcolor(WHITE);
    //     std::cout << "intersections: " << '\n';
    //     for(int i = 0; i < intersections.size(); i++){
    //         if(drawing) _cross(intersections[i]);
    //         std::cout << "  intersection " << i << ": " << intersections[i] << '\n';
    //     }
    // }
    // chip.intersect();
    // chip.make_edges();
    // chip.test_edges();
    // chip.make_faces();
    // std::cout << "chip:\n" << chip.dbg() << '\n';
    chip.color(50*WIDTH/2000, 50*WIDTH/2000);
    // chip.draw_net(5, 5, 0, -1, false, true);
    // BC->setcolor(BLUE);
    // chip.draw(BC);
    std::cout << "test_knot_1 finished" << '\n';
}
void test_intersect_linear(){
    Point A = Point(Vector(.45, .7), Vector(-1, 0));
    Point B = Point(Vector(.15, .9), Vector(.5, .5));
    Point C = Point(Vector(.45, .1), Vector( 1, 0));
    Point D = Point(Vector(.75, .9), Vector(-1, 0));
    Point E = Point(Vector(.75, .1), Vector( 1, 0));

    Chip chip = Chip({A, B, C, D, E});
    BC = new BasicCanvas();
    if(drawing){
        // chip.draw(BC);
        chip.mark_points(BC);
    }
    //only the first that should be found
    // Intersection result[3];
    // int result_len = 0;
    // Spline(B, C).intersect_linear(Vector(0, 0), Vector(1, 1), result, &result_len);
    // for(int i = 0; i < result_len; i++){
    //     if(drawing) _cross(result[i].value);
    //     std::cout << "  intersection " << i << ": " << result[i].value.to_str() << '\n';
    // }
    std::vector<Vector> intersections = chip.intersect_linear();
    std::cout << "intersections: " << '\n';
    BC->setcolor(GREEN);
    for(int i = 0; i < intersections.size(); i++){
        if(drawing) _cross(intersections[i]);
        std::cout << "  intersection " << i << ": " << intersections[i] << '\n';
    }
    BC->dump("test_intersect_linear");
    BC->write("test_intersect_linear");
    delete BC;
    std::cout << "test_intersect_linear finished" << '\n';
}
void test_ABCD(int t_prec = 5, int v_prec = 5){
    Point A = Point(Vector(.15, .15), Vector(-.3,  .3));
    Point B = Point(Vector(.85, .15), Vector(-.3, -.3));
    Point C = Point(Vector(.85, .85), Vector( .3, -.3));
    Point D = Point(Vector(.15, .85), Vector( .3,  .3));
    Chip chip = Chip(std::vector<Point>({A, B, C, D}));
    if(drawing) chip.draw(BC);
    chip.intersect();
    chip.make_edges();
    chip.test_edges();
    // return;
    chip.make_faces();
    std::cout << "chip:\n" << chip.dbg() << '\n';
    chip.draw_net(BC, t_prec, v_prec);
    std::cout << "test_ABCD finished" << '\n';
}
void test_knot_2(){
    Point A = Point(Vector(.1, .3), Vector( .7, -2));
    Point B = Point(Vector(.3, .6), Vector(-.4, 0));
    Point C = Point(Vector(.4, .3), Vector( .7, .7));
    Point D = Point(Vector(.8, .6), Vector( 0,  .5));
    Point E = Point(Vector(.45,.5), Vector( 0, -.5));
    Point F = Point(Vector(.8, .1), Vector(-.7,-.3));
    Point G = Point(Vector(.5, .7), Vector(-1, .3));
    Chip chip = Chip({A, B, C, D, E, F, G});

    if(drawing){
        BC->clear();
        chip.draw(BC);
        chip.mark_points(BC);
    }
    chip.intersect();
    chip.make_edges();
    chip.test_edges();
    chip.make_faces();
    chip.draw_net(BC);
    std::cout << "chip:\n" << chip.dbg() << '\n';
    std::cout << "test_ABCDE finished" << '\n';
}
void test_P(){
    SplineConstruct O = SplineConstruct({
        Spline(
            Point(Vector(+0.23, +0.71), Vector(-0.50, +0.05)),
            Point(Vector(+0.15, +0.90), Vector(+0.37, +0.37))),
        Spline(
            Point(Vector(+0.15, +0.90), Vector(+0.19, +0.19)),
            Point(Vector(+0.23, +0.71), Vector(+0.02, -0.46)))
    }, {0.743479, 0.386428});
    Spline A = O.approximate();
    Spline B = O.flipped().approximate();
    if(drawing){
        O.draw(BC);
        BC->setcolor(BLUE);
        A.draw(BC);
        BC->setcolor(GREEN);
        B.draw(BC);
    }
}
void test_create(){
    #if 0
    if(!drawing){
        std::cerr << "No clicking without drawing!" << '\n';
        return;
    }
    int x, y, x2, y2;
    std::vector<Point> points;
    Chip chip = Chip(points);
    Vector A, a;

    for(int c = 0; c < 5; c++){
        // getMousePressRelease(&x, &y, &x2, &y2);
        A = screen_to_coordinate(Vector(x, y));
        a = screen_to_coordinate(Vector(x2, y2))-A;
        points.push_back(Point(A, a*3.0f));
        std::cout << "added Point to Chip: " << points[c].dbg() << std::endl;
        BC->clear();
        chip = Chip(points);
        chip.draw(BC);
    }
    BC->setcolor(BLUE);
    std::vector<Vector> intersections = chip.intersect();
    for(int i = 0; i < intersections.size(); i++){
        _cross(intersections[i], .01);
    }
    std::cout << "test_create finished" << '\n';
    #endif
}
void test_random(int n=5){
    std::vector<Point> points;
    Chip chip = Chip(points);
    Vector A, a;
    for(int c = 0; c < n; c++){
        A = Vector(random_float(.2, .8), random_float(.2, .8));
        a = Vector(random_float(-.5, .5), random_float(-.5, .5));
        points.push_back(Point(A, a));
        std::cout << "added Point to Chip: " << points[c].dbg() << std::endl;
        chip = Chip(points);
        if(drawing){
            BC->clear();
            chip.draw(BC);
        }
    }
    if(drawing) BC->setcolor(BLUE);
    std::vector<Vector> intersections = chip.intersect();
    if(drawing){
        for(int i = 0; i < intersections.size(); i++){
            _cross(intersections[i], .01);
        }
    }
    std::cout << "test_random finished" << '\n';
}
void test_PQR(
    int t_prec = 5, int v_prec = 5, bool print = false,
    Vector P = Vector(0.2, 0.25),
    Vector p1 = Vector(0.7, 0),
    Vector p2 = Vector(-0.05, -0.5),
    Vector Q = Vector(0.7, 0.15),
    Vector q1 = Vector(0.2, 1),
    Vector q2 = Vector(1.2, 0.25),
    Vector R = Vector(0.45, 0.75),
    Vector r1 = Vector(-0.3, -0.25),
    Vector r2 = Vector(-0.3, 1.25)
){
    #if 0
    std::vector<Point> points;
    points.push_back(Point(P, p2));
    points.push_back(Point(P, p1));
    points.push_back(Point(Q, q2));
    points.push_back(Point(Q, q1));
    points.push_back(Point(R, r2));
    points.push_back(Point(R, r1));

    Chip chip = Chip(points);
    if(drawing){
        BC->setcolor(BLUE);
        chip.draw(BC);
    }

    Spline A = Spline(P, Q,  p1,  q2);
    Spline B = Spline(P, R, -p2, -r1);
    Spline C = Spline(Q, R,  q1,  r2);
    Spline* E[t_prec+1];
    Spline* F[v_prec];
    F[0] = new Spline(P, C(0), p1, q2);
    if(drawing) F[0]->draw(BC);
    float t, v;
    for(int v_ = 1; v_ <= v_prec; v_++){
        v = (float)v_/v_prec;
        F[v_] = new Spline(
            P,
            C(v),
            -p2*v + p1*(1.0f-v),
            -r1*v + q2*(1.0f-v)
        );
        if(drawing){
            F[v_]->draw(BC);
            for(float t = 0; t < 1; t += 1.0f/t_prec){
                BC->line(T(t), T(t) + Vector(.02, 0));
            }
        }
        if(print) std::cout << "F[" << v_ << "]: " << F[v_]->latex() << '\n';
    }
    if(drawing) BC->setcolor(GREEN);
    E[0] = new Spline(P, P, Vector(0, 0), Vector(0, 0));
    for(int t_ = 1; t_ <= t_prec; t_++){
        t = (float)t_/t_prec;
        E[t_] = new Spline(
            A(t),
            B(t),
            (q1*t - p2*(1.0f-t))*t,
            (r2*t - p1*(1.0f-t))*t
        );
        if(drawing){
            E[t_]->draw(BC);
            for(float v = 0; v < 1; v += 1.0f/v_prec){
                BC->line(T(t), T(t) + Vector(.02, 0));
            }
        }
        if(print) std::cout << "E[" << t_ << "]: " << E[t_]->latex() << '\n';
    }
    std::cout << "test_PQR finished" << '\n';
    #endif
}

#endif /* end of include guard: CHIPS_TESTS_CPP */
