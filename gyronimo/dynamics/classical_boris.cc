#include <gyronimo/dynamics/classical_boris.hh>

namespace gyronimo {

// Performs a time step `dt` update to the state `in` and returns the result.
classical_boris::state classical_boris::do_step(const classical_boris::state &in, const double t, const double dt) {
	
	// extract position and velocity from state
	IR3 q_old = get_position(in);
	IR3 u_old = get_velocity(in);

	// solve in cartesian coordinates
	if(field_morph_) {

		// coordinate conversion to cartesian
		dIR3 ek = field_morph_->del(q_old);
		IR3 x_old = (*field_morph_)(q_old);
		IR3 v_old = contraction<second>(ek, u_old);

		// calculate fields
		IR3 Ek = {0, 0, 0};
		IR3 Bk = {0, 0, 0};
		if(electric_) {
			Ek = electric_->contravariant(q_old, t * iEfield_time_factor_);
			Ek = contraction<second>(ek, Ek);
		}
		if(magnetic_) {
			Bk = magnetic_->contravariant(q_old, t * iBfield_time_factor_);
			Bk = contraction<second>(ek, Bk);
		}

		// perform cartesian step
		auto [x_new, v_new] = cartesian_boris(x_old, v_old, Oref_, Ek, Bk, dt);

		// conversion back to curvilinear coordinates
		IR3 q_new = field_morph_->inverse(x_new);
		IR3 u_new = field_morph_->to_contravariant(q_new, v_new);

		return {q_new[IR3::u], q_new[IR3::v], q_new[IR3::w],
				u_new[IR3::u], u_new[IR3::v], u_new[IR3::w]};

	} else {

		// calculate fields
		IR3 Ek = {0, 0, 0};
		IR3 Bk = {0, 0, 0};
		if(electric_) Ek = electric_->contravariant(q_old, t * iEfield_time_factor_);
		if(magnetic_) Bk = magnetic_->contravariant(q_old, t * iBfield_time_factor_);

		// perform cartesian step
		auto [x_new, v_new] = cartesian_boris(q_old, u_old, Oref_, Ek, Bk, dt);

		return {x_new[IR3::u], x_new[IR3::v], x_new[IR3::w],
				v_new[IR3::u], v_new[IR3::v], v_new[IR3::w]};
	}

}

// Returns the vector position of the state normalized to `Lref`.
IR3 classical_boris::get_position(const classical_boris::state& s) const {
	return {s[0], s[1], s[2]};
}

// Returns the vector velocity of the state normalized to `Vref`.
IR3 classical_boris::get_velocity(const classical_boris::state& s) const {
	return {s[3], s[4], s[5]};
}

// Returns the kinetic energy of the state, normalized to `Uref`.
double classical_boris::get_kinetic_energy(const classical_boris::state& s) const {

	if(field_morph_) {

		IR3 q = {s[0], s[1], s[2]};
		IR3 u = {s[3], s[4], s[5]};
		IR3 v = contraction<second>(field_morph_->del(q), u);
		return v[IR3::u] * v[IR3::u] + v[IR3::v] * v[IR3::v] + v[IR3::w] * v[IR3::w];

	} else return s[3]*s[3] + s[4]*s[4] + s[5]*s[5];
}

// Returns the `classical_boris::state` from a normalized point in phase-space.
classical_boris::state classical_boris::generate_state(const IR3 &pos, const IR3 &vel) const {
	return {pos[IR3::u], pos[IR3::v], pos[IR3::w],
			vel[IR3::u], vel[IR3::v], vel[IR3::w]};
}

// Creates the first `classical_boris::state` from a normalized point in cartesian phase-space.
classical_boris::state classical_boris::generate_initial_state(const IR3 &cartesian_position, const IR3 &cartesian_velocity, const double &t, const double &dt) const {
	
	if(field_morph_) {
		
		// coordinate conversion to cartesian
		IR3 q_old = field_morph_->inverse(cartesian_position);
		dIR3 ek = field_morph_->del(q_old);

		// calculate fields
		IR3 Ek = {0, 0, 0};
		IR3 Bk = {0, 0, 0};
		if(electric_) {
			Ek = electric_->contravariant(q_old, t * iEfield_time_factor_);
			Ek = contraction<second>(ek, Ek);
		}
		if(magnetic_) {
			Bk = magnetic_->contravariant(q_old, t * iBfield_time_factor_);
			Bk = contraction<second>(ek, Bk);
		}

		auto [xmh, vmh] = cartesian_boris(cartesian_position, cartesian_velocity, Oref_, Ek, Bk, -0.5*dt);

		dIR3 eek = inverse(ek);

		return generate_state(q_old, contraction<second>(eek, vmh));

	} else {

		// calculate fields
		IR3 Ek = {0, 0, 0};
		IR3 Bk = {0, 0, 0};
		if(electric_) Ek = electric_->contravariant(cartesian_position, t * iEfield_time_factor_);
		if(magnetic_) Bk = magnetic_->contravariant(cartesian_position, t * iBfield_time_factor_);

		auto [xmh, vmh] = cartesian_boris(cartesian_position, cartesian_velocity, Oref_, Ek, Bk, -0.5*dt);

		return generate_state(cartesian_position, vmh);
	}

}

// Performs a boris step in cartesian coordinates.
std::pair<IR3, IR3> classical_boris::cartesian_boris(const IR3 &x_old, const IR3 &v_old, const double &Oref, const IR3 &Ek, const IR3 &Bk, const double &dt) const {

	
	// step 1
	IR3 half_E_impulse = {0, 0, 0};
	IR3 v_new = v_old;

	if(electric_) {
		half_E_impulse = (0.5 * Oref * dt) * Ek;
		v_new += half_E_impulse;
	}

	// step 2
	if(magnetic_) {

		double B = std::sqrt(inner_product(Bk, Bk));
		IR3 b = B!=0 ? (1/B) * Bk : Bk;

		double t = std::tan((0.5 * Oref * dt) * B);
		double s = (2 * t) / (1 + t*t);

		IR3 v_prime = v_new + t * cartesian_cross_product(v_new, b);
		v_new += s * cartesian_cross_product(v_prime, b);
	}

	// step 3
	if(electric_) v_new += half_E_impulse;

	// step 4
	IR3 x_new = x_old + v_new * dt;

	return {x_new, v_new};
}

} // end namespace gyronimo