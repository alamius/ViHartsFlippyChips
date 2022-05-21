#ifndef CANVAS_INCLUDE_H
#define CANVAS_INCLUDE_H

#include "include/canvas/color.hpp"
#include "include/canvas/Vector.hpp"

extern colorint WHITE[COLOR_LEN];
extern colorint BLUE[COLOR_LEN];
extern colorint CYAN[COLOR_LEN];
extern colorint GREEN[COLOR_LEN];
extern colorint DARKGREEN[COLOR_LEN];
extern colorint YELLOW[COLOR_LEN];
extern colorint RED[COLOR_LEN];
extern colorint TRANSPARENT[COLOR_LEN];
void colors_init();

template <typename CanvasT> void std_line(CanvasT* C, Vector a, Vector b);
template <typename CanvasT> void _line (  CanvasT* C, Vector a, Vector b);
template <typename CanvasT> void _cross(  CanvasT* C, Vector a, float size = .03);
template <typename CanvasT> void _plus (  CanvasT* C, Vector a, float size = .03);


#endif /* end of include guard: CANVAS_INCLUDE_H */
