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

static const int argument_spacing = 20;
std::string argument_message(std::string arg, int arg_spacing = argument_spacing, std::string message = ""){
    for(int i = 0; i < arg.size(); i++){
        if(arg[i] == '='){
            arg = arg.substr(0, i);
        }
    }
    std::string result = arg + ":";
    while(result.size() < arg_spacing){
        result += " ";
    }
    result += message;
    return result;
}
const std::string help = "command line options of chip:\n"
    "<filename>          sets the filename under which the output is written\n"
    "--chip              uses a predefined Chip\n"
    "--chip=<values>     uses the values for the Points of the Chip, see --help-chip-values for details\n"
    "--chip-file=<file>  reads <file> as values\n"
    "--color             colors and exports the created Chip to <filename>.ppm\n"
    "--draw              draws and exports the outline of the created Chip to <filename>.draw.ppm\n"
    "--gauss             activates gaussian blurring for the subsequent images; deactivate with --no-gauss\n"
    "--help              prints this help\n"
    "--help-chip-values  explains the use of in-line chip values (also applicable to chip files)\n"
    "--regular[=<n>]     creates a regular <n>-faced Chip\n"
    "--regular-factor=<f>sets the factor for the direction vectors of the n-faced Chip to <f>, use before --regular!\n"
    "--dbg[=<n>]         sets the debug level to <n> (0 to 5, where 5 is super verbose)\n"
;

const std::string help_chip_values = "using chip values:\n"
    "general form:      .012.345->.-67.-89[*1[.23]] (.- means either . or -)\n"
    "as regex:          [.-]\\d+[.-]\\d+->[.+]\\d[.+]\\d<\\*\\d+<.\\d+>> (<...> is optional)\n"
    "example:           .15.85->.3-3*5.5\n"
    "meaning:           a 2d vector (.15, .85), followed by a 2d direction (.3, -.3), optionally scaled by an arbitrary number (5.5)\n"
    "negative numbers:  .75 is 0.75 while -75 is -0.75, so the - replaces the . \n"
    "                   as no values with an absolute greater than 1 are allowed \n"
    "                   (thats the reason for the weird factor)\n"
    "ranges:            the first vector (before ->) must always be in [0, 1]×[0, 1]\n"
    "                   the second is free\n"
;

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
                std::cout << "using dbg level " << dbg_file_lvl << "." << '\n';
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
