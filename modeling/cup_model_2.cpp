#include "numerical/geometry.h"
#include "UI/stl_encoder.h"
#include "triangulate/octatree.h"
#include "triangulate/parametric_surface_adaptive_dist.h"
#include <vector>
#include <functional>


// in cm
const double h = 6.4;
const double r0 = 2.0;
const double r1 = 2.8;
const double h_b = 0.4;
const double r_b = r0 + (r1 - r0)*(h_b / h);
double ring_r = 0.1;
const vec2 n = normalize(vec2(h, r0 - r1));
const vec2 ring_pos = vec2(r1, h) + ring_r * n;
const double half_thickness = 0.015;



/* Implicit surface modeling */


// https://www.iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
double segment_dist(vec2 p, vec2 a, vec2 b) {
	vec2 pa = p - a, ba = b - a;
	double h = dot(pa, ba) / dot(ba, ba);
	return length(pa - ba * clamp(h, 0., 1.));
}

double map(vec3 p) {
	vec2 s = vec2(length(p.xy()), p.z);
	double side = segment_dist(s, vec2(r0, 0), vec2(r1, h));
	double bottom = segment_dist(s, vec2(0, h_b), vec2(r_b, h_b));
	double ring = length(s - ring_pos) - ring_r;
	return min(min(side, bottom), ring) - half_thickness;
}



/* Parametric surface modeling */

// point of intersection
vec2 POI(vec2 P0, vec2 d0, vec2 P1, vec2 d1) {
	double t = det(d1, P1 - P0) / det(d1, d0);
	return P0 + d0 * t;
}
vec2 POI_circ(vec2 p0, vec2 d, vec2 c, double r) {
	vec2 p = p0 - c;
	double c2 = d.sqr(), c1 = dot(p, d), c0 = p.sqr() - r * r;
	double t = (-sqrt(c1 * c1 - c0 * c2) - c1) / c2;
	return p0 + t * d;
}

// parametric surface
vec3 map_p(double u, double v) {
	const vec2 p0 = vec2(0, h_b - half_thickness);
	const vec2 p2 = vec2(r0, 0) - half_thickness * n;
	const vec2 p1 = POI(p0, vec2(1, 0), p2, n.rot());
	const vec2 p3 = vec2(r0, 0) + half_thickness * n;
	const vec2 p4 = POI_circ(p3, n.rot(), ring_pos, ring_r + half_thickness);
	const vec2 p5 = vec2(r1, h) - half_thickness * n;
	const vec2 p7 = vec2(0, h_b + half_thickness);
	const vec2 p6 = POI(p7, vec2(1, 0), p5, n.rotr());
	v *= 7.;
	vec2 c;
	if (v <= 1.) {
		c = mix(p0, p1, v);
	}
	else if (v <= 2.) {
		c = mix(p1, p2, v - 1.);
	}
	else if (v <= 3.) {
		double t = PI * (v - 2.);
		c = vec2(r0, 0) - half_thickness * n*cos(t) - half_thickness * n.rot()*sin(t);
	}
	else if (v <= 4.) {
		c = mix(p3, p4, v - 3.);
	}
	else if (v <= 5.) {
		double a0 = atan2(p4.y - ring_pos.y, p4.x - ring_pos.x);
		double a1 = atan2(p5.y - ring_pos.y, p5.x - ring_pos.x);
		double t = a0 + (a1 - a0)*(v - 4.);
		c = ring_pos + (ring_r + half_thickness)*cossin(t);
	}
	else if (v <= 6.) {
		c = mix(p5, p6, v - 5.);
	}
	else {
		c = mix(p6, p7, v - 6.);
	}
	return vec3(c.x*cossin(u), c.y);
}




int main(int argc, char* argv[]) {


	const int UN = 400;
	const int VN = 400;
	double tol = 1. / sqrt(UN*VN);

	const double thickness = 0.015;
	const double eps = 1e-5;

	std::vector<triangle_3d> trigs;
	//trigs = ScalarFieldTriangulator_octatree::octatree_cylindrical(map, 5., -2., 8., 60, 20, 25, 3);

	std::vector<triangle_3d> trigs_p;
	trigs_p = AdaptiveParametricSurfaceTriangulator_dist(map_p).triangulate_adaptive(0., 2.*PI, 0., 1., 12, 15, 12, 1. / 400., true, false);
	trigs.insert(trigs.begin(), trigs_p.begin(), trigs_p.end());


	// place it on the plane
	if (0) {
		double min_z = INFINITY;
		for (int i = 0; i < (int)trigs.size(); i++) {
			min_z = min(min(min_z, trigs[i][0].z), min(trigs[i][1].z, trigs[i][2].z));
		}
		for (int i = 0; i < (int)trigs.size(); i++) {
			trigs[i][0].z -= min_z, trigs[i][1].z -= min_z, trigs[i][2].z -= min_z;
		}
	}

	// output
	writeSTL(argv[1], &trigs[0], trigs.size(), nullptr, STL_CCW);
	return 0;
}

