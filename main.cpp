#include <stdio.h>
#include <iostream>
// #include <graphics.h>
#include <chrono>
#include <thread>
#include <vector>
#include <random>

//own
#include <Vector.h>
static int dbg_file_lvl = 0;
#include <utils.random.h>
// #include <graphics.utils.h>
// #include <mouse.h>
const unsigned int WIDTH = 2000, HEIGHT = 1600;

#ifndef COLOR_LEN
#define COLOR_LEN 2
#endif

#include <basic_canvas.h>
#include <canvas_layered.h>
BasicCanvas* BC;
LayeredCanvas* LC;
#include "canvas.include.h"
const bool drawing = 1;

#include <Spline.h>
std::string filename;
colorint (*write_bg_color)[COLOR_LEN];
#include "Chip.h"
#include "test.cpp"

int main(int argc, char const *argv[]){
    colors_init();
    #if COLOR_LEN == 4
        write_bg_color = &DARKGREEN;
    #else
        write_bg_color = &TRANSPARENT;
    #endif
    make_kernel_gauss();
    if(dbg_file_lvl >= 2){
        std::cout << "size of BasicCanvas: " << sizeof(BasicCanvas)/1000 << "kB" << '\n';
        std::cout << "size of LayeredCanvas: " << sizeof(LayeredCanvas)/1000 << "kB" << '\n';
    }
    if(argc <= 1){
        filename = "image";
    }else{
        filename = argv[1];
    }
    srand(time(NULL));
    // initMouse();

    Point A = Point(Vector(.45, .7), Vector(-1, 0));
    Point B = Point(Vector(.15, .9), Vector(.5, .5));
    Point C = Point(Vector(.45, .1), Vector( 1, 0));
    Point D = Point(Vector(.75, .9), Vector(-1, 0));
    Point E = Point(Vector(.75, .1), Vector( 1, 0));

    Chip chip = Chip({A, B, C, D, E});
    chip.color(
        50 * WIDTH/2000,
        50 * WIDTH/2000,
        false, //draw_E
        true, //draw_F
        1, //pensize
        &std_color_func, //color_func
        false //apply_gauss
    );

    // test_PQR(10, 10);
    // test_snake();
    // test_Splines();
    // test_SplineConstructs();
    // test_SplineConstruct_approximate();
    // test_basic_canvas();
    // test_BC_LC();
    // test_ABCD();
    // test_intersect_linear();
    // test_knot_1(false);
    // test_knot_2();
    // test_P();
    // test_create();
    // for(int a = 0; a < 5; a++){
    //     test_random(7);
    //     getch();
    // }
    // if(drawing) getch();

    // char wait;
    // std::cin >> wait;
    std::cout << "finished" << '\n';

    // std::this_thread::sleep_for(std::chrono::milliseconds(50));
    // endMouse();
    // std::cout << "write image to " << filename << ".ppm" << '\n';
    // BC->write(filename);
    // delete BC;
    // delete LC;
    return 0;
}
