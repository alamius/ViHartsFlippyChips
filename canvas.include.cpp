#ifndef CANVAS_INCLUDE_H
#define CANVAS_INCLUDE_H

#include "include/canvas/color.hpp"
#include "include/canvas/Vector.hpp"

colorint BLACK[COLOR_LEN];
colorint WHITE[COLOR_LEN];
colorint BLUE[COLOR_LEN];
colorint CYAN[COLOR_LEN];
colorint GREEN[COLOR_LEN];
colorint DARKGREEN[COLOR_LEN];
colorint YELLOW[COLOR_LEN];
colorint RED[COLOR_LEN];
colorint TRANSPARENT[COLOR_LEN];
void colors_init(){
	make_color(BLACK, 0, 0, 0);
	make_color(BLUE, 0, 0, 255);
	make_color(WHITE, 255, 255, 255);
	make_color(CYAN, 0, 255, 255);
	make_color(GREEN, 0, 255, 0);
	make_color(DARKGREEN, 0, 64, 0);
	make_color(YELLOW, 128, 255, 0);
	make_color(RED, 255, 0, 0);
	make_color(TRANSPARENT, 0, 0, 0, 0);
}

template <typename CanvasT>
void std_line(CanvasT* C, Vector a, Vector b){
	C->line(a, b);
}
template <typename CanvasT>
void _line(CanvasT* C, Vector a, Vector b){
	C->line(a, b);
}
template <typename CanvasT>
void _cross(CanvasT* C, Vector a, float size = .03){
	C->line(
		a + Vector(-size, -size),
		a + Vector(+size, +size)
	);
	C->line(
		a + Vector(-size, +size),
		a + Vector(+size, -size)
	);
}
template <typename CanvasT>
void _plus(CanvasT* C, Vector a, float size = .03){
	C->line(
		a + Vector(-size, 0),
		a + Vector(+size, 0)
	);
	C->line(
		a + Vector(0, -size),
		a + Vector(0, +size)
	);
}

#include "include/canvas/BasicCanvas.hpp"
template void _plus<BasicCanvas>(BasicCanvas*, Vector, float);
template void _cross<BasicCanvas>(BasicCanvas* C, Vector a, float size = .03);

#endif /* end of include guard: CANVAS_INCLUDE_H */
