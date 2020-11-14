#include <stdio.h>
#include <iostream>
// #include <graphics.h>
#include <chrono>
#include <thread>
#include <vector>
#include <random>

//own
#include <Vector.h>
static int dbg_file_lvl = 3;
#include <utils.random.h>
// #include <graphics.utils.h>
// #include <mouse.h>
const int WIDTH = 2000, HEIGHT = 1600;

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
    make_kernel_gauss();
    write_bg_color = &DARKGREEN;
    std::cout << "size of BasicCanvas: " << sizeof(BasicCanvas)/1000 << "kB" << '\n';
    std::cout << "size of LayeredCanvas: " << sizeof(LayeredCanvas)/1000 << "kB" << '\n';
    if(argc <= 1){
        filename = "image";
    }else{
        filename = argv[1];
    }
    srand(time(NULL));
    // initMouse();

    // test_PQR(10, 10);
    // test_snake();
    // test_Splines();
    // test_SplineConstructs();
    // test_SplineConstruct_approximate();
    // test_basic_canvas();
    // test_BC_LC();
    // test_ABCD();
    // test_intersect_linear();
    test_knot_1(false);
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
