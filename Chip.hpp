#ifndef CHIP_HPP
#define CHIP_HPP

#include <vector>

#include "include/canvas/BasicCanvas.hpp"
#include "include/canvas/LayeredCanvas.hpp"
#include "include/canvas/Vector.hpp"

#include "Edge.hpp"
#include "Node.hpp"
#include "Face.hpp"

colorint* std_color_func(colorint result[COLOR_LEN], float t, float v);
colorint std_alpha_map(colorint alpha);

Vector interpolate(Vector p, Vector q, float t);

class Chip {
private:
	//points (Vector of position and Vector of direction combined, defined in Spline.h) of the Chip. the line of the chip is the combination of the splines from every point to the next
	std::vector<Point> points;
	//will be created at the intersections of the line defined by points, containing the intersection Vector and the indices of the points and correcponding parameters (see Nodes)
	std::vector<Node> nodes;
	//the faces of the chip with their edges
	std::vector<Face> faces;
public:
	//copies the points
	Chip(std::vector<Point> points_);
	//colors all faces that are Face::inside with t_prec*v_prec small quadrilaterals with color given by color_func(t, v)
	void color(
		int t_prec,
		int v_prec,
		colorint* (*color_func)(colorint result[COLOR_LEN], float t, float v),
		bool apply_gauss
		bool apply_gauss,
		bool apply_gauss_after_bg
	);
	//colors a area defined by four Splines in t_prec*v_prec small quadrilaterals with color given by color_func(t, v)
	template <typename CanvasT>
	void color_stripe(
		CanvasT* C,
		int t_prec,
		int v_prec,
		colorint* (*color_func)(colorint result[COLOR_LEN], float t, float v),
		const Vector& P, const Vector& p1, const Vector& p2,
		const Vector& Q, const Vector& q1, const Vector& q2,
		const Vector& R, const Vector& r1, const Vector& r2,
		const Vector& S, const Vector& s1, const Vector& s2
	);
	//returns all intersections of the Line with itself. mistakes might be caused by Spline::intersect which has faults
	std::vector<Vector> intersect();
	//returns all intersections of the Line with a straight line A + a*t for t in R
	std::vector<Vector> intersect_linear(Vector A, Vector a);
	//subroutine of make_edges: takes intersection and follows the Line forewards of backwards
	void follow_edge(int p_curr, float t_curr, int sign, int d, int from, bool SplineConstruct_approximate);
	//looks through all Nodes that are set by Chip::intersect and finds the Edges that connect them (every Node has four Edges)
	void make_edges(bool SplineConstruct_approximate);
	//goes through all Edges, following them and connecting them into faces
	void make_faces();
	//draw the Line defined by Chip::points onto Canvas C
	template <typename CanvasT>
	void draw(CanvasT* C, int samples = 30);
	//draw a net of the Line onto Canvas C
	template <typename CanvasT>
	void draw_net(
		CanvasT* C,
		int t_prec = 5,
		int v_prec = 5,
		int face_from = 0,
		int face_to = -1, // -1 is turned into faces.size() - 1
		bool draw_E = true,
		bool draw_F = false,
		bool drawing = true
	);
	//marking the private Chip::points on Canvas C
	template <typename CanvasT>
	void mark_points(CanvasT* C);
	//transform all Points in Chip::points according to matrix
	void transform(float a, float b, float c, float d, float e, float f); // [x, y][[a, b], [c, d]] + [e, f]
	//return a latex compatible text version of the Line
	std::string latex();
	std::string dbg(std::string indent = "");
	virtual ~Chip(){};
};

#endif /* end of include guard: CHIP_HPP */
