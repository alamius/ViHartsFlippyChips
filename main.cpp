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

//for include/canvas_color.h used by include/basic_canvas.h and include/canvas_layered.h
#ifndef COLOR_LEN
#define COLOR_LEN 4
#endif

//for include/canvas_kernel.h
#ifndef CANVAS_KERNEL_SIZE
#define CANVAS_KERNEL_SIZE 9
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

bool is_key(std::string arg, std::string key){ //--key* -> true
    if(arg.size() < 2 + key.size()){
        return false;
    }
    return (arg.substr(2, key.size()) == key);
}
int get_equal_sign_index(std::string arg){
    int result = -1; //first index of '='
    for(int i = 0; i < arg.size(); i++){
        if(arg[i] == '='){
            result = i;
            break;
        }
    }
    if(result == -1){
        if(dbg_file_lvl >= 1) std::cout << "get_equal_sign_index: no equal sign in arg '" << arg << "'. return -1.\n";
        return result;
    }
    return result + 1; //index after '='
}
std::vector<std::string> get_values(const std::string arg, const char sep = ','){ //--key=arg1<sep>arg2...
    std::vector<std::string> result = std::vector<std::string>();
    int index = get_equal_sign_index(arg);
    if(index == -1){
        if(dbg_file_lvl >= 1) std::cout << "get_values: no equal sign in arg '" << arg << "'. return empty result.\n";
        return result;
    }
    int i;
    for(i = 0; index + i + 1 < arg.size(); i++){
        if(arg[index + i] == sep){
            result.push_back(arg.substr(index, i));
            index += i + 1; //index after sep
            i = 0;
        }
    }
    if(index + i < arg.size()){
        result.push_back(arg.substr(index, arg.size() - index));
    }
    return result;
}
int get_int(const std::string arg){ //--key=int
    int result = 0;
    int index = get_equal_sign_index(arg);
    if(index == -1){
        if(dbg_file_lvl >= 1) std::cout << "get_int: no equal sign in arg '" << arg << "'. return -1.\n";
        return -1;
    }
    for(int i = index; i < arg.size(); i++){
        if('0' <= arg[i] && arg[i] <= '9'){
            result *= 10;
            result += arg[i] - '0';
        }else{
            std::cout << "found invalid character in int in arg '" << arg << "'. return -1\n";
            return -1;
        }
    }
    return result;
}
float get_float(const std::string arg){ //--key=float
    float result = 0;
    int index = get_equal_sign_index(arg);
    if(index == -1){
        if(dbg_file_lvl >= 1) std::cout << "get_float: no equal sign in arg '" << arg << "'. return 0.\n";
        return 0;
    }
    int i = index;
    int sign;
    if(arg[i] == '-'){
        sign = -1;
        i++;
    }else{
        sign = +1;
    }
    for(; i < arg.size(); i++){
        if('0' <= arg[i] && arg[i] <= '9'){
            result *= 10;
            result += sign * (arg[i] - '0');
        }else if(arg[i] == '.'){
            break;
        }else{
            std::cout << "found invalid character in int in arg '" << arg << "'. return -1\n";
            return -1;
        }
    }
    if(i < arg.size()){
        int base_index = i + 1;
        for(i = base_index; i < arg.size(); i++){
            result += sign * (arg[i] - '0') * pow(.1, i - base_index + 1);
        }
    }
    return result;
}
std::string get_string(const std::string arg){ //--key=string
    std::string result = "";
    int index = get_equal_sign_index(arg);
    if(index == -1){
        if(dbg_file_lvl >= 1) std::cout << "get_int: no equal sign in arg '" << arg << "'. return \"\".\n";
        return result;
    }
    return arg.substr(index, arg.size() - index);
}
std::vector<Vector> to_Vectors(std::vector<std::string> values){
    std::vector<Vector> result = std::vector<Vector>();
    int base_index;
    float* coordinate;
    for(int v = 0; v < values.size(); v++){
        result.push_back(Vector(0, 0));
        coordinate = &(result.back().x);
        base_index = 1; //after '.'
        for(int i = 1; i < values[v].size(); i++){
            if(values[v][i] == '.'){
                coordinate = &(result.back().y);
                base_index = i + 1;
            }else{
                if('0' <= values[v][i] && values[v][i] <= '9'){
                    *coordinate += (values[v][i] - '0') * pow(.1, i - base_index + 1);
                }
            }
        }
    }
    return result;
}
std::vector<Point> to_Points(std::vector<std::string> values){
    std::vector<Point> result = std::vector<Point>();
    int base_index;
    Vector* vector;
    float* coordinate;
    int sign = +1;
    float factor;
    for(int v = 0; v < values.size(); v++){
        //create new Vector to be filled with value
        result.push_back(Point(Vector(0, 0), Vector(0, 0)));
        //make pointer to place vector
        vector = &(result.back().V);
        //make pointer to place vector x-coordinate
        coordinate = &(vector->x);
        base_index = 1; //after '.'
        for(int i = 1; i < values[v].size(); i++){
            if( //detects "->" -> not place but direction vector
                values[v][ i ] == '-' && values[v].size() > i + 1 &&
                values[v][i+1] == '>'
            ){
                vector = &(result.back().v);
                i++;
            }else if(values[v][i] == '*'){ //a factor for direction vector
                if(vector == &(result.back().v)){
                    factor = 0;
                    i++;
                    if(i == values[v].size()){
                        std::cerr << "empty argument multiplication. ignored" << '\n';
                        break;
                    }
                    if(values[v][i] == '-'){
                        sign = -1;
                        i++;
                    }else{
                        sign = +1;
                    }
                    for(; i < values[v].size(); i++){
                        if(values[v][i] == '.'){
                            break;
                        }else if('0' <= values[v][i] && values[v][i] <= '9'){
                            factor *= 10;
                            factor += values[v][i] - '0';
                        }else{
                            std::cerr << "unknown symbol in float: " << values[v][i] << " in " << values[v] << '\n';
                        }
                    }
                    if(i < values[v].size()){
                        base_index = i + 1;
                        for(i = base_index; i < values[v].size(); i++){
                            if('0' <= values[v][i] && values[v][i] <= '9'){
                                factor += sign * (values[v][i] - '0') * pow(.1, i - base_index + 1);
                            }else{
                                std::cerr << "unknown symbol in float: " << values[v][i] << " in " << values[v] << '\n';
                            }
                        }
                    }
                    (*vector) *= factor;
                }else{
                    std::cerr << "vector argument multiplication only allowed on direction vector: '.5.5->.5.5*5'. ignoring in '" << values[v] << "'" << '\n';
                }
            }else if(values[v][i] == '.' || values[v][i] == '-'){ //float of format "[\.-][\d]+"
                if(coordinate == &(vector->x)){
                    coordinate = &(vector->y);
                }else{
                    coordinate = &(vector->x);
                }
                base_index = i + 1;
                sign = values[v][i] == '.' ? +1 : -1;
            }else if('0' <= values[v][i] && values[v][i] <= '9'){
                *coordinate += sign * (values[v][i] - '0') * pow(.1, i - base_index + 1);
            }else{
                std::cerr << "invalid character in '" << values[v] << "': << '" << values[v][i] << "'" << '\n';
            }
        }
    }
    return result;
}

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
