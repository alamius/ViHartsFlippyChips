#ifndef COMMANDLINE_HPP
#define COMMANDLINE_HPP

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


#endif /* end of include guard: COMMANDLINE_HPP */
