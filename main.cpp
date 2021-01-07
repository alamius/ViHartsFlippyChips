#include <stdio.h>
#include <iostream>
// #include <graphics.h>
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
#include "Commandline.hpp"

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
    filename = "image";
    srand(time(NULL));
    // initMouse();
    Chip* chip = NULL;

    std::string arg;
    std::vector<std::string> values;
    std::vector<Point> points;
    bool apply_gauss = false;
    int t_prec = 50;
    int v_prec = 50;
    int regular = 3;
    float regular_factor = 3;
    for(int a = 1; a < argc; a++){
        arg = argv[a];
        if(arg.substr(0, 2) == "--"){
            if(is_key(arg, "help")){
                if(is_key(arg, "help-chip-values")){
                    std::cout << help_chip_values << '\n';
                }else{
                    std::cout << help << '\n';
                }
            }else if(is_key(arg, "chip")){
                if(is_key(arg, "chip-file")){
                    if(is_key(arg, "chip-file=")){
                        std::string chip_file_name = get_string(arg);
                        std::cout << argument_message(arg) << "file name read as: '" << chip_file_name << "'\n";
                        std::ifstream chip_file(chip_file_name, std::ios::in);
                        if(chip_file.is_open()){
                            std::string line;
                            while(std::getline(chip_file, line)){
                                values.push_back(line);
                            }
                            chip_file.close();
                        }else{
                            std::cerr << argument_message(arg) << "Unable to open file '" << chip_file_name << "'\n";
                            continue;
                        }
                        points = to_Points(values);
                        for(int i = 0; i < values.size(); i++){
                            std::cout << argument_message(arg) << "values[" << i << "]: " << values[i] << " -> " << points[i].dbg() << '\n';
                        }
                    }else{
                        std::cerr << argument_message(arg) << "must get an argument: --chip-file=<file1>\n";
                    }
                }else{
                    if(is_key(arg, "chip=")){
                        values = get_values(arg, ',');
                        points = to_Points(values);
                        for(int i = 0; i < values.size(); i++){
                            std::cout << argument_message(arg) << "values[" << i << "] == " << values[i] << " -> " << points[i].dbg() << '\n';
                        }
                    }else{
                        points.push_back(Point(Vector(.45, .7), Vector(-1, 0)));
                        points.push_back(Point(Vector(.15, .9), Vector(.5, .5)));
                        points.push_back(Point(Vector(.45, .1), Vector( 1, 0)));
                        points.push_back(Point(Vector(.75, .9), Vector(-1, 0)));
                        points.push_back(Point(Vector(.75, .1), Vector( 1, 0)));
                        std::cout << argument_message(arg) << "using standard Chip Points" << '\n';
                    }
                }
                chip = new Chip(points);
                std::cout << argument_message(arg) << "defined Chip as: " << chip->dbg() << '\n';
            }else if(is_key(arg, "regular-factor=")){
                regular_factor = get_float(arg);
                std::cout << argument_message(arg) << "using regular factor " << regular_factor << " (must be set before --regular!)" << '\n';
            }else if(is_key(arg, "regular")){
                if(is_key(arg, "regular=")){
                    regular = get_int(arg);
                }else{
                    std::cout << argument_message(arg) << "using standard regular with 3 points, use --regular=<n> for others." << '\n';
                    regular = 3;
                }
                int P = regular;
                float t;
                for(int p = 0; p < P; p++){
                    t = M_PI * 2 * p / P;
                    points.push_back(Point(
                        Vector(
                            -sin(t*(P-1))/6 + sin(t)/4 + .5,
                            cos( t*(P-1))/6 + cos(t)/4 + .5
                        ), Vector(
                            -cos(t*(P-1))*4/6 + cos(t)/4,
                            -sin(t*(P-1))*4/6 - sin(t)/4
                        )*regular_factor
                    ));
                }
                std::cout << argument_message(arg) << "using regular Chip with " << regular << " points. (subsequent --chip* overwrites this again!)" << '\n';
                chip = new Chip(points);
            }else if(is_key(arg, "color")){
                if(chip == NULL){
                    std::cerr << argument_message(arg) << "called without Chip! create one with --chip* or --regular before this option." << '\n';
                    exit(1);
                }
                chip->color(
                    t_prec * WIDTH/2000,
                    50 * WIDTH/2000,
                    &std_color_func,
                    apply_gauss
                );
            }else if(is_key(arg, "draw")){
                if(chip == NULL){
                    std::cerr << argument_message(arg) << "called without Chip! create one with --chip* or --regular before this option." << '\n';
                    exit(1);
                }
                BC = new BasicCanvas();
                chip->draw(BC);
                BC->write(filename+".draw");
            }else if(is_key(arg, "gauss")){
                std::cout << argument_message(arg) << "using gauss." << '\n';
                apply_gauss = true;
            }else if(is_key(arg, "no-gauss")){
                std::cout << argument_message(arg) << "not using gauss." << '\n';
                apply_gauss = false;
            }else if(is_key(arg, "dbg")){
                if(is_key(arg, "dbg=")){
                    dbg_file_lvl = get_int(arg);
                }else{
                    dbg_file_lvl = 1;
                }
                std::cout << argument_message(arg) << "using dbg level " << dbg_file_lvl << "." << '\n';
            }else{
                std::cerr << argument_message(arg) << "the option '" << arg << "' was not understood!" << '\n';
            }
        }else{
            filename = arg;
            std::cout << argument_message(arg) << "output file set to '" << filename << "'\n";
        }
    }
    std::cout << "finished" << '\n';

    // endMouse();
    return 0;
}
