#include "modeling/parametric/surfaces.h"
#include "triangulate/parametric_surface_adaptive_dist.h"
#include "UI/stl_encoder.h"

int main(int argc, char* argv[]) {

	const int i = 40;

	std::vector<triangle_3d> temp;
	auto S = ParamSurfaces[i];
	auto info = ParamSurfaceInfo::info[i];
	printf("%d - %s\n", i, S.name);
	temp = AdaptiveParametricSurfaceTriangulator_dist(S.P).triangulate_adaptive(S.u0, S.u1, S.v0, S.v1, 64, 64, 16, 0.001*pow(determinant(info.InertiaTensor_u), 1. / 6.), false, false);
	int TN = temp.size();
	translateToCOM_shell(&temp[0], TN);
	scaleGyrationRadiusTo_shell(&temp[0], TN, 1.5);

	// balanced position
	mat3 M = axis_angle(cross(info.minGravPotential_vec, vec3(0, 1e-8, -1)), acos(-info.minGravPotential_vec.z));
	M = rotationMatrix_z(0.8*PI)*M;
	for (int i = 0; i < TN; i++) temp[i].applyMatrix(M);
	translateShape_onPlane(&temp[0], TN);

	writeSTL(argv[1], &temp[0], temp.size(), nullptr, STL_CCW);

}
