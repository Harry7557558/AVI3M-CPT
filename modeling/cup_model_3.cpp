#include "numerical/geometry.h"
#include "UI/stl_encoder.h"
#include "triangulate/octatree.h"
#include "triangulate/parametric_surface_adaptive_dist.h"
#include <vector>
#include <functional>


// in cm
const vec2 s0 = vec2(2.4, 0);
const vec2 s1 = vec2(3.6, 8.5);
const vec2 s2 = vec2(3.8, 8.7);
const vec2 s3 = vec2(4.0, 10.0);
const vec2 sb = vec2(2.45, 0.7);
const vec2 sb0 = vec2(0, sb.y);
const double half_thickness = 0.05;



/* Implicit surface modeling */

// https://www.iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
double segment_dist(vec2 p, vec2 a, vec2 b) {
	vec2 pa = p - a, ba = b - a;
	double h = dot(pa, ba) / dot(ba, ba);
	return length(pa - ba * clamp(h, 0., 1.));
}

double map(vec3 p) {
	vec2 s = vec2(length(p.xy()), p.z);
	double side = min(
		min(segment_dist(s, s0, sb),
			segment_dist(s, sb, s1)
		),
		min(segment_dist(s, s1, s2),
			segment_dist(s, s2, s3)));
	double bottom = segment_dist(s, sb0, sb);
	return min(side, bottom) - half_thickness;
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
	const vec2 n0 = half_thickness * vec2(0, -1);
	const vec2 n1 = half_thickness * normalize(sb - s0).rotr();
	const vec2 n2 = half_thickness * normalize(s1 - sb).rotr();
	const vec2 n3 = half_thickness * normalize(s2 - s1).rotr();
	const vec2 n4 = half_thickness * normalize(s3 - s2).rotr();
	const vec2 p0 = sb0 + n0;
	const vec2 p2 = s0 - n1;
	const vec2 p1 = POI(p0, n0.rot(), p2, n1.rot());
	const vec2 p3 = s0 + n1;
	const vec2 p4 = POI(p3, n1.rot(), sb + n2, n2.rot());
	const vec2 p5 = POI(sb + n2, n2.rot(), s1 + n3, n3.rot());
	const vec2 p6 = s2 + n3;
	const vec2 p7 = s2 + n4;
	const vec2 p8 = s3 + n4;
	const vec2 p9 = s3 - n4;
	const vec2 p10 = POI(p9, n4.rot(), s2 - n3, n3.rot());
	const vec2 p11 = s1 - n3;
	const vec2 p12 = s1 - n2;
	const vec2 p13 = POI(p12, n2.rot(), sb - n0, n0.rot());
	const vec2 p14 = sb0 - n0;
	v *= 14.;
	vec2 c;
	if (v <= 1.) {
		c = mix(p0, p1, v);
	}
	else if (v <= 2.) {
		c = mix(p1, p2, v - 1.);
	}
	else if (v <= 3.) {
		double t = PI * (v - 2.);
		c = s0 - n1 * cos(t) + n1.rotr()*sin(t);
	}
	else if (v <= 4.) {
		c = mix(p3, p4, v - 3.);
	}
	else if (v <= 5.) {
		c = mix(p4, p5, v - 4.);
	}
	else if (v <= 6.) {
		c = mix(p5, p6, v - 5.);
	}
	else if (v <= 7.) {
		double t = mix(atan2(n3.y, n3.x), atan2(n4.y, n4.x), v - 6.);
		c = s2 + half_thickness * cossin(t);
	}
	else if (v <= 8.) {
		c = mix(p7, p8, v - 7.);
	}
	else if (v <= 9.) {
		double t = PI * (v - 8.);
		c = s3 + n4 * cos(t) + n4.rot()*sin(t);
	}
	else if (v <= 10.) {
		c = mix(p9, p10, v - 9.);
	}
	else if (v <= 11.) {
		c = mix(p10, p11, v - 10.);
	}
	else if (v <= 12.) {
		double t = mix(atan2(-n3.y, -n3.x), atan2(-n2.y, -n2.x), v - 11.);
		c = s1 + half_thickness * cossin(t);
	}
	else if (v <= 13.) {
		c = mix(p12, p13, v - 12.);
	}
	else if (v <= 14.) {
		c = mix(p13, p14, v - 13.);
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
	//trigs = ScalarFieldTriangulator_octatree::octatree_cylindrical(map, 5., -2., 12., 60, 20, 30, 3);

	std::vector<triangle_3d> trigs_p;
	trigs_p = AdaptiveParametricSurfaceTriangulator_dist(map_p).triangulate_adaptive(0., 2.*PI, 0., 1., 12, 31, 12, 1. / 400., true, false);
	trigs.insert(trigs.begin(), trigs_p.begin(), trigs_p.end());


	// place it on the plane
	if (1) {
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

