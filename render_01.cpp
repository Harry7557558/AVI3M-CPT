#include <stdio.h>
#include <vector>
#include <algorithm>
#include "numerical/geometry.h"
#include "raytracing/bvh.h"
#include "preview/brdf.h"

#include "preview/gui.h"
void Init() {}
void render() {}




BVH* BVH_R = 0;


texture nautilus_macromphalus_2("textures\\nautilus_macromphalus_2.png");
texture neptunea_arthritica_2("textures\\neptunea_arthritica_2.png");
texture haliotis_discus_2("textures\\haliotis_discus_2.png");
texture pitar_lupanaria_2("textures\\pitar_lupanaria_2.png");

texture table_texture("textures\\1f7dca9c22f324751f2a5a59c9b181dfe3b5564a04b724c657732d0bf09c99db.jpg");  // Shadertoy wood texture


struct cup {
	const vec2 cupdim = vec2(20., 9.);
	vec3 pos;
	double scale;
	texture* tex;
	vec2 tex_sc;
	mat2 tex_ang;
	cup(vec3 pos, double scale, texture* tex, double tex_scale, double tex_ang) {
		this->pos = pos, this->scale = scale;
		this->tex = tex;
		this->tex_sc = (1. / tex_scale) * (cupdim / vec2(tex->res)) * max(tex->res.x / cupdim.x, tex->res.y / cupdim.y);
		this->tex_ang = rotationMatrix2d(tex_ang);
	}
	bool intersect(vec3 ro, vec3 rd, double &t, vec3 &n) {
		bool k = intersectScene(BVH_R, ro - pos, rd, t, n);
		return k;
	}
	vec3 getTexture(vec3 p) {
		p -= pos;
		vec2 pxy = tex_ang * p.xy();
		vec2 uv = vec2(atan2(pxy.y, pxy.x) / (2.*PI) + 0.5, p.z / cupdim.y);
		uv = (uv - vec2(0, 0.5)) * tex_sc + vec2(0, 0.5);
		if (uv.x > 0. && uv.x < 1. && uv.y > 0. && uv.y < 1.)
			return tex->fetch(uv.x, uv.y);
		return vec3(1.);
	}
};

std::vector<cup> Cups({
	cup(vec3(180, 110, 80), 1.0, &neptunea_arthritica_2, 1.0, 1.5*PI),
	});




vec3 calcCol(vec3 ro, vec3 rd, uint32_t &seed) {

	vec3 m_col = vec3(1.), col;
	vec3 tot_col(0.);

	const bool DLS = true;  // direct light sample

							// "recursive" ray-tracing
	const int MAX_ITER = 24;  // set to 32 in final rendering
	for (int iter = 0; iter < MAX_ITER; iter++) {
		vec3 n, min_n;
		vec2 uv, min_uv;
		double min_t = INFINITY, t;
		ro += 1e-6*rd;  // alternate of t>1e-6

		int intersect_id = -1;  // which object the ray hits
		const int PLANE = 0x00, WALL = 0x01, CEILING = 0x02,
			TABLE = 0x10, CUP_0 = 0x100,
			CEILING_LIGHT = 0x20, WINDOW_LIGHT = 0x21;
		const bool CEILING_LIGHT_ON = false;
		const bool WINDOW_LIGHT_ON = true;

		const vec3 CUP_POS = vec3(180, 110, 80);

		auto intersectObjects = [&](vec3 ro, vec3 rd, int &intersect_id, double &min_t, vec3 &min_n, vec2 &min_uv, bool BreakWhenIntersect = false) {

			double min_t_0 = min_t;
			double t; vec3 n; vec2 uv;

			// intersect plane
			if ((t = -ro.z / rd.z) > 0. && t < min_t) {
				min_t = t;
				min_n = vec3(0, 0, ro.z > 0. ? 1 : -1.);
				intersect_id = PLANE;
				if (BreakWhenIntersect) return true;
			}

			// intersect wall
			if (1) {
				if (intersectParallelogram(vec3(0, 0, 0), vec3(360, 0, 0), vec3(0, 0, 220), ro, rd, t = min_t, n, uv)) {
					min_t = t, min_n = n, min_uv = uv, intersect_id = WALL;
					if (BreakWhenIntersect) return true;
				}
				if (intersectParallelogram(vec3(0, 0, 0), vec3(0, 220, 0), vec3(0, 0, 220), ro, rd, t = min_t, n, uv)) {
					min_t = t, min_n = n, min_uv = uv, intersect_id = WALL;
					if (BreakWhenIntersect) return true;
				}
				if (intersectParallelogram(vec3(0, 220, 0), vec3(360, 0, 0), vec3(0, 0, 220), ro, rd, t = min_t, n, uv)) {
					min_t = t, min_n = n, min_uv = uv, intersect_id = WALL;
					if (BreakWhenIntersect) return true;
				}
				if (intersectParallelogram(vec3(360, 0, 0), vec3(0, 220, 0), vec3(0, 0, 220), ro, rd, t = min_t, n, uv)) {
					min_t = t, min_n = n, min_uv = uv, intersect_id = WALL;
					if (BreakWhenIntersect) return true;
				}
			}

			// intersect ceiling
			if (1) {
				if (intersectParallelogram(vec3(0, 0, 220), vec3(360, 0, 0), vec3(0, 220, 0), ro, rd, t = min_t, n, uv)) {
					min_t = t, min_n = n, min_uv = uv;
					intersect_id = CEILING;
					if (BreakWhenIntersect) return true;
				}
			}

			// intersect table
			if (1) {
				const double TABLE_L = 110, TABLE_W = 50;
				if (intersectParallelogram(vec3(180 - TABLE_L / 2, 110 - TABLE_W / 2, 80), vec3(TABLE_L, 0, 0), vec3(0, TABLE_W, 0), ro, rd, t = min_t, n, uv)) {
					min_t = t, min_n = n, min_uv = uv;
					intersect_id = TABLE;
					if (BreakWhenIntersect) return true;
				}
			}

			// intersect cup
			if (1) {
				for (int i = 0; i < (int)Cups.size(); i++) {
					if (Cups[i].intersect(ro, rd, t = min_t, n)) {
						min_t = t, min_n = n;
						intersect_id = CUP_0 + i;
					}
				}
			}

			return min_t < min_t_0;

		};


		auto intersectLights = [&](vec3 ro, vec3 rd, int &intersect_id, double &min_t, vec3 &min_n, vec2 &min_uv) {

			double min_t_0 = min_t;
			double t; vec3 n; vec2 uv;

			// intersect ceiling light
			if (CEILING_LIGHT_ON) {
				const int SIZE = 30;
				if (intersectParallelogram(vec3(180 - SIZE, 110 - SIZE, 220 - 5), vec3(2 * SIZE, 0, 0), vec3(0, 2 * SIZE, 0), ro, rd, t = min_t, n, uv)) {
					min_t = t, min_n = n, min_uv = uv;
					intersect_id = CEILING_LIGHT;
				}
			}

			// intersect window light
			if (WINDOW_LIGHT_ON) {
				const double PAD = 40, BOTTOM = 80;
				if (intersectParallelogram(vec3(0.01, PAD, BOTTOM), vec3(0, 220 - 2 * PAD, 0), vec3(0, 0, 220 - PAD - BOTTOM), ro, rd, t = min_t, n, uv)) {
					int i = int(11 * uv.x), j = int(11 * uv.y);
					if ((i + j) & 1) {
						min_t = t, min_n = n, min_uv = uv;
						intersect_id = WINDOW_LIGHT;
					}
				}
			}

			return min_t < min_t_0;
		};


		// intersect objects and lights
		bool hitObjects = intersectObjects(ro, rd, intersect_id, min_t, min_n, min_uv);
		bool hitLights = intersectLights(ro, rd, intersect_id, min_t, min_n, min_uv);
		if (hitLights) hitObjects = false;


		// color calculation and update ray

		// background
		if (intersect_id == -1) {
			return vec3(0.);
			// calculate incoming light
			vec3 col = vec3(5.*pow(max(dot(rd, normalize(vec3(0.1, 0.3, 1))), 0.), 5.));
			return m_col * col;
		}

		//min_n = dot(rd, min_n) < 0. ? min_n : -min_n;

		{
			// plane
			if (intersect_id == PLANE) {
				// ray calculation
				vec3 wi = -rd;
				vec3 wo = randdir_cosWeighted(min_n, seed);
				// color calculation
				vec2 p = ro.xy() + min_t * rd.xy();
				col = int(floor(p.x) + floor(p.y)) & 1 ? vec3(0.8) : vec3(0.6);
				col = vec3(0.9, 0.8, 0.8);
				m_col *= col;
				// update ray
				ro = ro + rd * min_t;
				rd = wo;
			}

			// wall
			if (intersect_id == WALL) {
				// ray calculation
				vec3 wi = -rd;
				vec3 wo = randdir_cosWeighted(min_n, seed);
				// color calculation
				col = vec3(0.9);
				m_col *= col;
				// update ray
				ro = ro + rd * min_t;
				rd = wo;
			}

			// ceiling
			if (intersect_id == CEILING) {
				// ray calculation
				vec3 wi = -rd;
				vec3 wo = randdir_cosWeighted(min_n, seed);
				// color calculation
				col = vec3(0.9);
				m_col *= col;
				// update ray
				ro = ro + rd * min_t;
				rd = wo;
			}

			// table
			if (intersect_id == TABLE) {
				// ray calculation
				vec3 wi = -rd;
				//vec3 wo = randdir_cosWeighted(min_n, seed);
				vec3 wo = mix(randdir_cosWeighted(min_n, seed), rd - 2.*dot(rd, min_n)*min_n, 0.);
				// color calculation
				col = vec3(0.4);
				col *= vec3(.6) + .8*table_texture.fetch(2.*min_uv.x, min_uv.y);
				m_col *= col;
				// update ray
				ro = ro + rd * min_t;
				rd = wo;
			}

			// cup
			if (intersect_id >= CUP_0) {
				// ray calculation
				vec3 wi = -rd;
				vec3 wo = randdir_cosWeighted(min_n, seed);
				// color calculation
				{
					col = Cups[intersect_id - CUP_0].getTexture(ro + rd * min_t);
					col *= vec3(0.9, 0.9, 0.95);
					m_col *= col;
				}
				// update ray
				ro = ro + rd * min_t;
				rd = wo;
			}

		}

		{
			// ceiling light
			if (intersect_id == CEILING_LIGHT) {
				if (DLS && iter != 0) return tot_col;
				// calculate incoming light
				col = vec3(15.*abs(dot(rd, min_n)));
				return tot_col + m_col * col;
			}

			// ceiling light
			if (intersect_id == WINDOW_LIGHT) {
				if (DLS && iter != 0) return tot_col;
				// calculate incoming light
				col = vec3(10.*abs(dot(rd, min_n)));
				return tot_col + m_col * col;
			}
		}


		// (approximated) direct light sample
		if (DLS) {

			// ceiling light
			if (CEILING_LIGHT_ON) {
				double u = rand01(seed), v = rand01(seed);
				const int SIZE = 30;
				vec3 p = vec3(180 - SIZE, 110 - SIZE, 220 - 5) + 2 * SIZE * vec3(u, v, 0);
				vec3 sd = normalize(p - ro);
				double sdn = dot(sd, min_n);
				if (sdn > 0.) {
					if (!intersectLights(ro + 1e-6*sd, sd, intersect_id, t = INFINITY, min_n, uv)) fprintf(stderr, "N");  // should hit
					if (!intersectObjects(ro + 1e-6*sd, sd, intersect_id, t, n, uv, true)) {
						double sphang = 4. * SIZE * SIZE / (p - ro).sqr() * sdn;
						tot_col += m_col * vec3(15. * sphang * abs(dot(sd, min_n)) / 4.);
					}
				}
			}

			// window light, seems to be too biased, use for preview
			if (WINDOW_LIGHT_ON) {
				double u, v;
				do {
					u = rand01(seed), v = rand01(seed);
				} while (((int(11.*u) + int(11.*v)) & 1) == 0);
				const double PAD = 40, BOTTOM = 80;
				vec3 a = vec3(0, 220 - 2 * PAD, 0), b = vec3(0, 0, 220 - PAD - BOTTOM);
				vec3 p = vec3(0.01, PAD, BOTTOM) + a * u + b * v;
				vec3 sd = normalize(p - ro);
				double sdn = dot(sd, min_n);
				if (sdn > 0.) {
					if (!intersectLights(ro + 1e-6*sd, sd, intersect_id, t = INFINITY, min_n, uv)) printf("N");
					if (!intersectObjects(ro + 1e-6*sd, sd, intersect_id, t, n, uv, true)) {
						double sphang = length(a)*length(b) * (60. / 121.) / (p - ro).sqr() * sdn;
						tot_col += m_col * vec3(10. * sphang * abs(dot(sd, min_n))) / 4.;
					}
				}
			}
		}

		if (dot(m_col, vec3(1. / 3.)) < 1e-3) break;
	}
	return tot_col + m_col;
}






#define STB_IMAGE_WRITE_IMPLEMENTATION
#include ".libraries/stb_image_write.h"

const int W = 624 * 2, H = 361 * 2;


typedef unsigned char byte;
struct rgba {
	byte r, g, b, a;
} IMG[H][W];


BVH_Triangle *STL;
int STL_N;

int main() {

	std::vector<BVH_Triangle*> BT;
	readBinarySTL("modeling\\cup_model_3.stl", BT);
	BVH_R = new BVH;
	vec3 Min(INFINITY), Max(-INFINITY);
	constructBVH(BVH_R, BT, Min, Max);


	for (int Frame = 0; Frame < 5; Frame++) {

		if (Frame == 0) CamP = vec3(163.408657, 65.814623, 85.000000), ScrO = vec3(168.511802, 114.313749, 77.900684), ScrA = vec3(22.976396, -8.627498, 0.000000), ScrB = vec3(-0.000000, -0.000000, 14.198633);
		if (Frame == 1) CamP = vec3(183.658651, 157.046185, 85.943891), ScrO = vec3(192.245463, 109.190110, 77.902103), ScrA = vec3(-24.468910, 1.902879, 0.000000), ScrB = vec3(-0.022016, -0.283099, 14.195793);
		if (Frame == 2) CamP = vec3(130.004577, 128.159571, 86.596217), ScrO = vec3(184.497942, 123.086631, 76.999145), ScrA = vec3(-9.447228, -26.009323, 0.000000), ScrB = vec3(0.451344, -0.163939, 16.001710);
		if (Frame == 3) CamP = vec3(231.014517, 98.226568, 94.527099), ScrO = vec3(178.284971, 96.196159, 77.124866), ScrA = vec3(6.222721, 26.963174, 0.000000), ScrB = vec3(-2.792662, 0.644507, 15.750269);
		if (Frame == 4) CamP = vec3(163.902981, 154.207237, 88.771787), ScrO = vec3(191.336645, 114.731752, 77.923389), ScrA = vec3(-23.061521, -8.397307, 0.000000), ScrB = vec3(0.388230, -1.066196, 14.153221);

		printf("\n\n================= Frame #%d ==================\n", Frame);

		// rendering
		Render_Exec([](int beg, int end, int step, bool* sig) {
			const int WIN_SIZE = W * H;
			const int SPP = 256;

			for (int k = beg; k < end; k += step) {
				int i = k % W, j = k / W;
				vec3 col(0.);
				for (int u = 0; u < SPP; u++) {
					uint32_t seed = hashu((u*W + i)*H + j);
					vec3 CamD = ScrO + ((i + rand01(seed)) / W)*ScrA + ((j + rand01(seed)) / H)*ScrB;
					col += calcCol(CamP, normalize(CamD - CamP), seed);
				}
				col /= SPP;
				IMG[H - 1 - j][i] = rgba{
					byte(255.99*clamp(col.x, 0., 1.)),
					byte(255.99*clamp(col.y, 0., 1.)),
					byte(255.99*clamp(col.z, 0., 1.)),
					255 };
				if (beg == 0 && i == 0) printf("%.1lf%%\n", k * 100. / end);  // progress line
			}
			if (sig) *sig = true;
		}, W*H);


		// output
		stbi_write_png(&(std::string("render_01_") + char('0' + Frame) + ".png")[0], W, H, 4, &IMG[0][0], 4 * W);

	}

	return 0;
}

