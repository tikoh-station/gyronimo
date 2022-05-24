
#ifndef GYRONIMO_CYLINDRICAL_MORPHISM
#define GYRONIMO_CYLINDRICAL_MORPHISM

#include <cmath>

#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

class cylindrical_morphism : public morphism {

public:

	//! Class Constructor
	cylindrical_morphism() : morphism() {};

	//! Class Destructor
	~cylindrical_morphism() {};

	//! Maps the curvilinear coordinates `q` into cartesian coordinates `x`.
	IR3 operator()(const IR3 &q) const override;

	//! Inverse transform from cartesian coordinates `x` into curvilinear coordinates `q`.
	IR3 inverse(const IR3 &x) const override;

	//! Returns the morphism's derivatives, correspondent to the covariant basis vectors in point `q`.
	dIR3 del(const IR3 &q) const override;

	//! Returns the jacobian of the transformation in point `q`.
	double jacobian(const IR3 &q) const override;

	//! Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point `q`.
	dIR3 del_inverse(const IR3 &q) const override;

private:

}; // end class cylindrical_morphism

} // end namespace gyronimo

#endif // GYRONIMO_CYLINDRICAL_MORPH