

#include "gui.h"


#include <vector>
#include "raytracing/bvh.h"


int Trig_N;
BVH* BVH_R = 0;



#include "brdf.h"


texture nautilus_macromphalus_2("D:\\Homework\\AVI3M\\AVI3M-CPT\\textures\\nautilus_macromphalus_2.png");
texture neptunea_arthritica_2("D:\\Homework\\AVI3M\\AVI3M-CPT\\textures\\neptunea_arthritica_2.png");
texture haliotis_discus_2("D:\\Homework\\AVI3M\\AVI3M-CPT\\textures\\haliotis_discus_2.png");
texture pitar_lupanaria_2("D:\\Homework\\AVI3M\\AVI3M-CPT\\textures\\pitar_lupanaria_2.png");

texture wood_texture("D:\\Homework\\AVI3M\\AVI3M-CPT\\textures\\1f7dca9c22f324751f2a5a59c9b181dfe3b5564a04b724c657732d0bf09c99db.jpg");  // Shadertoy wood texture



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

	cup(vec3(50, 10, 0), 1.0, &nautilus_macromphalus_2, 1.0, 0.75*PI),
	cup(vec3(50, 20, 0), 1.0, &neptunea_arthritica_2, 1.0, 0.75*PI),
	cup(vec3(50, 30, 0), 1.0, &haliotis_discus_2, 1.0, 0.75*PI),
	cup(vec3(50, 40, 0), 1.0, &pitar_lupanaria_2, 1.0, 0.75*PI),

	});



// dome light
vec3 calcCol(vec3 ro, vec3 rd, uint32_t &seed) {

	vec3 m_col = vec3(1.2), col;
	vec3 tot_col(0.);

	const bool DLS = true;  // direct light sample

	// "recursive" ray-tracing
	const int MAX_ITER = 6;  // set to 32 in final rendering
	for (int iter = 0; iter < MAX_ITER; iter++) {
		vec3 n, min_n;
		vec2 uv, min_uv;
		double min_t = INFINITY, t;
		ro += 1e-6*rd;  // alternate of t>1e-6

		int intersect_id = -1;  // which object the ray hits
		const int PLANE = 0x00, WALL = 0x01, CEILING = 0x02,
			CEILING_LIGHT = 0x10, WINDOW_LIGHT = 0x11,
			CUP_0 = 0x100;
		const bool CEILING_LIGHT_ON = false;
		const bool WINDOW_LIGHT_ON = true;

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
				if (intersectParallelogram(vec3(359.99, PAD, BOTTOM), vec3(0, 220 - 2 * PAD, 0), vec3(0, 0, 220 - PAD - BOTTOM), ro, rd, t = min_t, n, uv)) {
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
				vec2 p = 0.2*(ro.xy() + min_t * rd.xy());
				col = int(floor(p.x) + floor(p.y)) & 1 ? vec3(1.0) : vec3(0.8);
				col = vec3(0.5) + 0.5*wood_texture.fetch(0.1*p.x, 0.1*p.y);
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
				col = vec3(0.95, 0.9, 0.9) * abs(min_n.x) + vec3(0.9, 0.95, 0.9) * abs(min_n.y);
				col = vec3(0.9, 0.92, 0.95);
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

			// window light
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
				vec3 p = vec3(359.99, PAD, BOTTOM) + a * u + b * v;
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



void render_RT_BVH() {
	static uint32_t call_time = 0;
	call_time = lcg_next(call_time);

	SPP++;

	Render_Exec([](int beg, int end, int step, bool* sig) {
		const int WIN_SIZE = _WIN_W * _WIN_H;
		for (int k = beg; k < end; k += step) {
			int i = k % _WIN_W, j = k / _WIN_W;
			// brute force Monte-Carlo sampling
			const int N = 1;
			vec3 col(0.);
			for (int u = 0; u < N; u++) {
				uint32_t seed = hashu(u*WIN_SIZE + k + call_time);
				vec3 ro = CamP;
				vec3 rd = scrDir(vec2(i + rand01(seed), j + rand01(seed)));
				if (0) {  // DOF
					double dist = length(ro - Center);
					vec3 fp = ro + rd * dist;
					double a = 2.*PI*rand01(seed), r = 0.01*dist*rand01(seed);
					ro = ro + (normalize(ScrA)*cos(a) + normalize(ScrB)*sin(a))*r;
					rd = normalize(fp - ro);
				}
				col += calcCol(ro, rd, seed);
			}
			if (SPP == 1) colorBuffer[i][j] = vec3(0.);
			colorBuffer[i][j] += col / N;
			Canvas(i, j) = toCOLORREF(colorBuffer[i][j] / SPP);
		}
		if (sig) *sig = true;
	}, _WIN_W*_WIN_H);
}



void render() {
	auto t0 = NTime::now();
	// initialize window
	for (int i = 0, l = _WIN_W * _WIN_H; i < l; i++) _WINIMG[i] = 0;
	for (int i = 0; i < _WIN_W; i++) for (int j = 0; j < _WIN_H; j++) _DEPTHBUF[i][j] = INFINITY;
	getScreen(CamP, ScrO, ScrA, ScrB);
	printf("W=%d,H=%d; CamP=vec3(%lf,%lf,%lf),ScrO=vec3(%lf,%lf,%lf),ScrA=vec3(%lf,%lf,%lf),ScrB=vec3(%lf,%lf,%lf);\n",
		_WIN_W, _WIN_H, CamP.x, CamP.y, CamP.z, ScrO.x, ScrO.y, ScrO.z, ScrA.x, ScrA.y, ScrA.z, ScrB.x, ScrB.y, ScrB.z);


	render_RT_BVH();


	double t = fsec(NTime::now() - t0).count();
	sprintf(text, "[%d×%d, %d]  %dspp  %.1fms (%.1ffps)\n", _WIN_W, _WIN_H, Trig_N, SPP, 1000.0*t, 1. / t);
	SetWindowTextA(_HWND, text);
}


// ============================================== User ==============================================



bool inited = false;
void Init() {
	if (inited) return; inited = true;

	std::vector<BVH_Triangle*> BT;
	readBinarySTL("D:\\Homework\\AVI3M\\AVI3M-CPT\\modeling\\cup_model_3.stl", BT);

	Center = vec3(0.);
	for (int i = 0; i < (int)Cups.size(); i++)
		Center += Cups[i].pos / Cups.size();
	Center += vec3(0, 0, 5);

	rz = 0.2*PI;

	BVH_R = new BVH;
	vec3 Min(INFINITY), Max(-INFINITY);
	constructBVH(BVH_R, BT, Min, Max);

}



