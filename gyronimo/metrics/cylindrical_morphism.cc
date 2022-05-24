#include <gyronimo/metrics/cylindrical_morphism.hh>

namespace gyronimo {

// Maps the curvilinear coordinates `q` into cartesian coordinates `x`.
IR3 cylindrical_morphism::operator()(const IR3 &q) const {
	return {q[IR3::u] * cos(q[IR3::v]), q[IR3::u] * sin(q[IR3::v]), q[IR3::w]};
}

// Inverse transform from cartesian coordinates `x` into curvilinear coordinates `q`.
IR3 cylindrical_morphism::inverse(const IR3 &x) const {
	return {sqrt(x[IR3::u]*x[IR3::u] + x[IR3::v]*x[IR3::v]),
			atan2(x[IR3::v], x[IR3::u]), x[IR3::w]};
}

// Returns the morphism's derivatives, correspondent to the covariant basis vectors in point `q`.
dIR3 cylindrical_morphism::del(const IR3 &q) const {
	double sn = sin(q[IR3::v]);
	double cn = cos(q[IR3::v]);
	double r  = q[IR3::u];
	return {
		cn, 	sn, 	0,
		- r * sn, r * cn, 0,
		0, 		0, 		1
	};
}

// Returns the jacobian of the transformation in point `q`.
double cylindrical_morphism::jacobian(const IR3 &q) const {
	return q[IR3::u];
}

// Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point `q`.
dIR3 cylindrical_morphism::del_inverse(const IR3 &q) const {
	double sn = sin(q[IR3::v]);
	double cn = cos(q[IR3::v]);
	double ir = 1 / q[IR3::u];
	return {
		cn, -sn * ir, 0,
		sn,  cn * ir, 0,
		0,		0   , 1
	};
}

}; // end namespace gyronimo