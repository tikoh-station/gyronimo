#ifndef GYRONIMO_ELECTROMAGNETIC_SYSTEM
#define GYRONIMO_ELECTROMAGNETIC_SYSTEM

#include <gyronimo/core/codata.hh>
#include <gyronimo/core/contraction.hh>
#include <gyronimo/fields/IR3field.hh>
#include <gyronimo/metrics/morphism.hh>

namespace gyronimo {

class electromagnetic_system {

public:

	typedef std::array<double, 6> state;

	electromagnetic_system(
		double Lref, double Vref, double qom,
		const IR3field *E, const IR3field *B,
		const morphism *morph);
	~electromagnetic_system() {};

	state operator()(const state& x, double t) const;

	IR3 get_position(const state &s) const;
	IR3 get_velocity(const state &s) const;

	state generate_state(const IR3& pos, const IR3 &vel) const;

	double Lref() const {return Lref_;};
	double Tref() const {return Tref_;};
	double Vref() const {return Vref_;};
	double Oref() const {return Oref_;};
	const IR3field* electric_field() const {return electric_field_;};
	const IR3field* magnetic_field() const {return magnetic_field_;};
	const morphism* morph() const {return field_morph_;};

private:

	const double Lref_, Vref_, Tref_;
	const double Oref_;
	const IR3field* electric_field_;
	const IR3field* magnetic_field_;
	const double iEfield_time_factor_, iBfield_time_factor_;
	const morphism *field_morph_;

}; // end class electromagnetic_system

} // end namespace gyronimo

#endif // GYRONIMO_ELECTROMAGNETIC_SYSTEM