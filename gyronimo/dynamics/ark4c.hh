
#ifndef GYRONIMO_ARK4C
#define GYRONIMO_ARK4C

#include <cmath>
#include <utility>

#include <gyronimo/core/error.hh>
#include <gyronimo/core/codata.hh>
#include <gyronimo/fields/IR3field.hh>
#include <gyronimo/metrics/morphism.hh>
#include <gyronimo/dynamics/electromagnetic_system.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>

#include <boost/numeric/odeint.hpp>

namespace gyronimo {

//! Gyronimo implementation for classical boris stepper class.
class ark4c {

public:

	//! State variable of the stepper.
	typedef std::array<double,6> state;

	//! Class Constructor.
	ark4c(const double Lref, const double Vref, const double qom, 
		const IR3field *E, const IR3field *B, const morphism *morph = NULL);

	//! Class Destructor.
	~ark4c() {};

	//! Returns the order of the method.
	unsigned short order() const {return 4;};

	//! Performs a time step `dt` update to the state `in` and returns the result.
	state do_step(const state &in, const double t, const double dt) const;

	//! Returns the vector position of the state normalized to `Lref`.
	IR3 get_position(const state& s) const;
	//! Returns the vector velocity of the state normalized to `Vref`.
	IR3 get_velocity(const state& s) const;
	//! Returns the kinetic energy of the state, normalized to `Uref`.
	double get_kinetic_energy(const state& s) const;
	double get_parallel_energy(const state& s, double &time) const;
	double get_perpendicular_energy(const state& s, double &time) const;

	//! Returns the `exact_stepper::state` from a normalized point in phase-space.
	state generate_state(const IR3 &position, const IR3 &velocity) const;

	//! Creates the first `exact_stepper::state` from a normalized point in cartesian phase-space.
	state generate_initial_state(const IR3 &cartesian_position, 
		const IR3 &cartesian_velocity, const double &tinit, const double &dt) const;

	
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

	//! Interpolate at time t \element [0, 1]
	IR3 interpolate(IR3 q0, IR3 q1, IR3 u0, IR3 u1, IR3 du0, const IR3::index i, double dt) const;

	const double Lref_, Vref_, Tref_;
	const double qom_, Oref_, iOref_;
	const IR3field *electric_;
	const IR3field *magnetic_;
	const double iBfield_time_factor_;
	const double iEfield_time_factor_;
	const morphism *field_morph_;

	const electromagnetic_system em_;
	const odeint_adapter<electromagnetic_system> sys_;

	const double a0_, b0_, c0_;
	const double alpha_, beta_, gamma_;
	const IR3 k0_, k1_, k2_, k3_, k4_;


	class straight_line_system {

	public:

		typedef std::array<double, 3> state;

		typedef boost::numeric::odeint::runge_kutta4<state> odeint_stepper;

		straight_line_system(const morphism *m, const state &dx) : m_(m), dx_(dx) {};

		void operator()(const state &in, state &out, const double &t) {
			IR3 res = m_->to_contravariant(get_IR3(in), dx_);
			out = generate_state(res);
			return;
		}

		state generate_state(const IR3 &q) const {
			return {q[IR3::u], q[IR3::v], q[IR3::w]};
		};

		IR3 get_IR3(const state &qs) const {
			return {qs[0], qs[1], qs[2]};
		};

	private:

		const morphism *m_;
		const IR3 dx_;

	}; // end class straight_line_system



	//! Performs a `runge-kutta` step in curvilinear coordinates across the line set by the vector dx.
	IR3 straight_line_step(const IR3 &q_old, const IR3 &dx) const;
	
	//! Calculates the cartesian displacement vector.
	std::pair<IR3, IR3> displacement(const IR3 &q_field, IR3 v_old, const double &t_field, const double &dt) const;

}; // end class ark4c

} // end namespace gyronimo

#endif // GYRONIMO_ARK4C