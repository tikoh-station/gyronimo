
#ifndef IR3FIELD_ANALYTICAL
#define IR3FIELD_ANALYTICAL

#include <functional>
#include <utility>

#include <gyronimo/fields/IR3field.hh>

namespace gyronimo {

class IR3field_analytical : public IR3field {

public:

	typedef std::function<IR3(IR3,double)> field_function;

	// constructors
	IR3field_analytical(double m_factor, double t_factor, const metric_covariant* g, const field_function &field_contravariant)
		: IR3field(m_factor, t_factor, g), field_contravariant_(std::move(field_contravariant)) {};
	
	// destructor
	~IR3field_analytical() override {};

	IR3 contravariant(const IR3& position, double time) const override {
		return field_contravariant_(position, time);
	};

	const field_function* function() const {return &field_contravariant_;};

private:

	const field_function field_contravariant_;

}; // end class IR3field_analytical

}; // end namespace gyronimo

#endif