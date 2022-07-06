#include <gyronimo/dynamics/cylindrical_boris.hh>

namespace gyronimo {

// Performs a time step `dt` update to the state `in` and returns the result.
cylindrical_boris::state cylindrical_boris::do_step(const cylindrical_boris::state &in, const double t, const double dt) const {
	
	// extract position and velocity from state
	IR3 q_old = cylindrical_boris::get_position(in);
	IR3 u_old = cylindrical_boris::get_velocity(in);

	// advance particle across phase-space
	IR3 Ek = {0, 0, 0};
	if(electric_) {
		Ek = electric_->contravariant(q_old, t/electric_->t_factor());
		Ek *= IR3({1, q_old[IR3::u], 1});
	}
	IR3 Bk = {0, 0, 0};
	if(magnetic_) {
		Bk = magnetic_->contravariant(q_old, t/magnetic_->t_factor());
		Bk *= IR3({1, q_old[IR3::u], 1});
	}

	IR3 x_old = {q_old[IR3::u], 0, q_old[IR3::w]};

	// advance with boris step
	auto [x_new, v_new] = classical_boris_update(x_old, u_old, Ek, Bk, dt);

	// std::cout << "v_new: [ " << v_new[IR3::u] << " " << v_new[IR3::v] << " " << v_new[IR3::w] << " ]" << std::endl;

	// convert back to cylindrical coordinates
	IR3 q_new = {
		sqrt(x_new[IR3::u] * x_new[IR3::u] + x_new[IR3::v] * x_new[IR3::v]),
		q_old[IR3::v] + atan(x_new[IR3::v] / x_new[IR3::u]),
		x_new[IR3::w]
	};

	double sa = x_new[IR3::v] / q_new[IR3::u];
	double ca = x_new[IR3::u] / q_new[IR3::u];
	IR3 u_new = {
		ca * v_new[IR3::u] + sa * v_new[IR3::v],
		ca * v_new[IR3::v] - sa * v_new[IR3::u],
		v_new[IR3::w]
	};

	return {q_new[IR3::u], q_new[IR3::v], q_new[IR3::w], 
			u_new[IR3::u], u_new[IR3::v], u_new[IR3::w]};
}

// Returns the vector position of the state normalized to `Lref`.
IR3 cylindrical_boris::get_position(const cylindrical_boris::state& s) {
	return {s[0], s[1], s[2]};
}

// Returns the vector velocity of the state normalized to `Vref`.
IR3 cylindrical_boris::get_velocity(const cylindrical_boris::state& s) {
	return {s[3], s[4], s[5]};
}

// Returns the kinetic energy of the state, normalized to `Uref`.
double cylindrical_boris::get_kinetic_energy(const cylindrical_boris::state& s) {
	return s[3]*s[3] + s[4]*s[4] + s[5]*s[5];
}

// Returns the `classical_boris::state` from a normalized point in phase-space.
cylindrical_boris::state cylindrical_boris::generate_state(const IR3 &pos, const IR3 &vel) const {
	return {pos[IR3::u], pos[IR3::v], pos[IR3::w],
			vel[IR3::u], vel[IR3::v], vel[IR3::w]};
}

// Creates the first `cylindrical_boris::state` from a normalized point in cartesian phase-space.
cylindrical_boris::state cylindrical_boris::generate_initial_state(
		const IR3 &cartesian_position, const IR3 &cartesian_velocity, const double &tinit, const double &dt) const {

	electromagnetic_system em(1, 1, qom_, electric_, magnetic_, field_morph_);
	odeint_adapter<electromagnetic_system> sys(&em);
	boost::numeric::odeint::runge_kutta4<electromagnetic_system::state> rk4;

	IR3 qk = cartesian_position;
	IR3 uk = cartesian_velocity;
	if(field_morph_) {
		IR3 qk = field_morph_->inverse(cartesian_position);
		IR3 uk = field_morph_->to_contravariant(qk, cartesian_velocity);
	}
	electromagnetic_system::state in = em.generate_state(qk, uk);
	electromagnetic_system::state out;

	rk4.do_step(sys, in, tinit, out, -0.5*dt);
	IR3 qmh = em.get_position(out);
	IR3 umh = em.get_velocity(out);

	if(field_morph_) {

		IR3 vmh = field_morph_->from_contravariant(qmh, umh);
		umh = field_morph_->to_contravariant(qk, vmh);
		return generate_state(qk, umh);

	} else return generate_state(cartesian_position, umh);

}

std::pair<IR3, IR3> cylindrical_boris::classical_boris_update(const IR3 &x_old, const IR3 &v_old, const IR3 &Ek, const IR3 &Bk, const double &dt) const {

	// advance particle across phase-space

	// step 1
	IR3 v_minus = v_old;
	IR3 half_E_impulse = {0, 0, 0};
	if(electric_) {
		half_E_impulse = 0.5 * Oref_ * dt * Ek;
		v_minus += half_E_impulse;
	}

	// step 2
	IR3 v_plus = v_minus;
	if(magnetic_) {
		double B = sqrt(inner_product(Bk, Bk));
		IR3 b = B!=0 ? (1/B) * Bk : Bk;

		double t = tan(0.5 * Oref_ * B * dt);
		double s = (2 * t) / (1 + t*t);

		IR3 v_prime = v_minus + t * cartesian_cross_product(v_minus, b);
		v_plus += s * cartesian_cross_product(v_prime, b);;
	}

	// step 3
	IR3 v_new = v_plus; 
	if(electric_) v_new += half_E_impulse;
	
	// step 4
	IR3 x_new = x_old + v_new * dt;

	return {x_new, v_new};
}

} // end namespace gyronimo