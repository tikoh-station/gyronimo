#ifndef GYRONIMO_MORPHISM
#define GYRONIMO_MORPHISM

#include <gyronimo/core/IR3algebra.hh>
#include <gyronimo/core/contraction.hh>

namespace gyronimo {

class morphism {

public:

	//! Class Constructor
	morphism() {};

	//! Class Destructor
	virtual ~morphism() {};

	//! Maps the curvilinear coordinates `q` into cartesian coordinates `x`.
	virtual IR3 operator()(const IR3 &q) const = 0;

	//! Inverse transform from cartesian coordinates `x` into curvilinear coordinates `q`.
	virtual IR3 inverse(const IR3 &x) const = 0;

	//! Returns the morphism's derivatives, correspondent to the covariant basis vectors in point `q`.
	virtual dIR3 del(const IR3 &q) const = 0;

	//! Returns the jacobian of the transformation in point `q`.
	virtual double jacobian(const IR3 &q) const;

	//! Returns the morphism's inverse derivatives, correspondent to the contravariant basis vectors in point `q`.
	virtual dIR3 del_inverse(const IR3 &q) const;

	//! Returns the covariant components of cartesian `A` in position `q`.
	virtual IR3 to_covariant(const IR3 &q, const IR3 &A) const;

	//! Returns the contravariant components of `A` in position `q`.
	virtual IR3 to_contravariant(const IR3 &q, const IR3 &A) const;

	//! Returns the cartesian vector from its covariant components `A` in position `q`.
	virtual IR3 from_covariant(const IR3 &q, const IR3 &A) const;

	//! Returns the cartesian vector from its contravariant components `A` in position `q`.
	virtual IR3 from_contravariant(const IR3 &q, const IR3 &A) const;

private:

}; // end class morphism

} // end namespace gyronimo

#endif // GYRONIMO_MORPHISM