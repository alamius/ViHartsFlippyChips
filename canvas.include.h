#ifndef CANVAS_INCLUDE_H
#define CANVAS_INCLUDE_H

colorint WHITE[COLOR_LEN];
colorint BLUE[COLOR_LEN];
colorint CYAN[COLOR_LEN];
colorint GREEN[COLOR_LEN];
colorint DARKGREEN[COLOR_LEN];
colorint YELLOW[COLOR_LEN];
colorint RED[COLOR_LEN];
colorint TRANSPARENT[COLOR_LEN];
void colors_init(){
	make_color(BLUE, 0, 0, 255);
	make_color(WHITE, 255, 255, 255);
	make_color(CYAN, 0, 255, 255);
	make_color(GREEN, 0, 255, 0);
	make_color(DARKGREEN, 0, 64, 0);
	make_color(YELLOW, 128, 255, 0);
	make_color(RED, 255, 0, 0);
	make_color(TRANSPARENT, 0, 0, 0, 0);
}

void std_line(Vector a, Vector b){
	BC->line(a, b);
}
void (*_line)(Vector a, Vector b){ std_line };
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
void _cross(Vector a, float size = .03){
	_cross(BC, a, size);
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
void _plus(Vector a, float size = .03){
	_plus(BC, a, size);
}


#endif /* end of include guard: CANVAS_INCLUDE_H */
