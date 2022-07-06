
#ifndef GYRONIMO_CENTERED_STEPPER
#define GYRONIMO_CENTERED_STEPPER

#include <cmath>

#include <gyronimo/core/codata.hh>
#include <gyronimo/fields/IR3field.hh>
#include <gyronimo/metrics/morphism.hh>
#include <gyronimo/dynamics/electromagnetic_system.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>	

#include <boost/numeric/odeint.hpp>

namespace gyronimo {

//! Gyronimo implementation for classical boris stepper class.
class centered_stepper {

public:

	//! State variable of the stepper.
	typedef std::array<double, 9> state;

	//! Class Constructor.
	centered_stepper(const double Lref, const double Vref, const double qom, 
			const IR3field *E, const IR3field *B, 
			const morphism *morph = NULL, const metric_covariant *metric = NULL)
			: Lref_(Lref), Vref_(Vref), Tref_(Lref/Vref), qom_(qom),
			Oref_(	B ? qom * gyronimo::codata::e / gyronimo::codata::m_proton * B->m_factor() : 
					E ? qom * gyronimo::codata::e / gyronimo::codata::m_proton * E->m_factor() : 0), 
			electric_(E), magnetic_(B),
			iBfield_time_factor_(B ? Tref_ / B->t_factor() : 0.0),
			iEfield_time_factor_(E ? Tref_ / E->t_factor() : 0.0), 
			field_morph_(morph), metric_(metric) {
	
		//@todo test if fields are consistent (normalization and cartesian metric)
	};

	//! Class Destructor.
	~centered_stepper() {};

	//! Returns the order of the method.
	unsigned short order() const {return 2;};

	//! Performs a time step `dt` update to the state `in` and returns the result.
	state do_step(const state &in, const double t, const double dt) const;

	//! Returns the vector position of the state normalized to `Lref`.
	IR3 get_position(const state& s);
	//! Returns the vector velocity of the state normalized to `Vref`.
	IR3 get_velocity(const state& s);
	//! Returns the kinetic energy of the state, normalized to `Uref`.
	double get_kinetic_energy(const state& s);

	//! Returns the `centered_stepper::state` from a normalized point in phase-space.
	state generate_state(const IR3 &position, const IR3 &velocity, const IR3 &cartesian_velocity) const;

	//! Creates the first `classical_boris::state` from a normalized point in cartesian phase-space.
	state generate_initial_state(const IR3 &cartesian_position, const IR3 &cartesian_velocity, const double &t, const double &dt) const;

	
	//! Returns the reference gyration frequency scale `Oref`.
	double Oref() const {return Oref_;};
	//! Returns the pointer to the electric field.
	const IR3field* electric() const {return electric_;};
	//! Returns the pointer to the magnetic field.
	const IR3field* magnetic() const {return magnetic_;};
	//! Returns the pointer to the morphism.
	const morphism* field_morph() const {return field_morph_;};
	//! Returns the pointer to the morphism.
	const morphism* inner_morph() const {return field_morph_;};

	//! Returns the reference length scale `Lref`.
	double Lref() const {return Lref_;};
	//! Returns the reference time scale `Tref`.
	double Tref() const {return Tref_;};
	//! Returns the reference velocity scale `Vref`.
	double Vref() const {return Vref_;};

private:

	const double Lref_, Vref_, Tref_;
	const double qom_, Oref_;
	const IR3field *electric_;
	const IR3field *magnetic_;
	const double iBfield_time_factor_;
	const double iEfield_time_factor_;
	const morphism *field_morph_;
	const metric_covariant *metric_;

	//! Performs a boris step in cartesian coordinates.
	IR3 cartesian_boris(const IR3 &vmh, const double &Oref, const IR3 &Ek, const IR3 &Bk, const double &dt) const;

}; // end class centered_stepper

} // end namespace gyronimo

#endif // GYRONIMO_CENTERED_STEPPER