#include <gyronimo/dynamics/boris_advanced.hh>

namespace gyronimo {

// Performs a time step `dt` update to the state `in` and returns the result.
boris_advanced::state boris_advanced::do_step(const boris_advanced::state &in, const double t, const double dt) const {

	// extract position and velocity from state
	IR3 qk = {in[0], in[1], in[2]};
	IR3 umh = {in[3], in[4], in[5]};

	if(field_morph_) {

		dIR3 ek  = field_morph_->del(qk);
		dIR3 eek = inverse(ek);
		IR3 vmh = contraction<second>(ek, umh);

		// calculate fields
		IR3 Ek = {0, 0, 0};
		IR3 Bk = {0, 0, 0};
		if(electric_) {
			Ek = electric_->contravariant(qk, t * iEfield_time_factor_);
			Ek = contraction<second>(ek, Ek);
		}
		if(magnetic_) {
			Bk = magnetic_->contravariant(qk, t * iBfield_time_factor_);
			Bk = contraction<second>(ek, Bk);
		}

		// perform cartesian step
		IR3 vph = cartesian_boris(vmh, Oref_, Ek, Bk, dt);
		
		// RK2 with minimum error coefficients
		IR3 uuph = contraction<second>(eek, vph);
		IR3 qpt = qk + (2./3. * dt) * uuph;
		dIR3 eept = field_morph_->del_inverse(qpt);
		IR3 upt = contraction<second>(eept, vph);

		// new state
		IR3 qp1 = qk + (0.25 * dt) * uuph + (0.75 * dt) * upt;
		IR3 uph = field_morph_->to_contravariant(qp1, vph);

		return {qp1[IR3::u], qp1[IR3::v], qp1[IR3::w],
				uph[IR3::u], uph[IR3::v], uph[IR3::w]};

	} else {

		// calculate fields
		IR3 Ek = {0, 0, 0};
		IR3 Bk = {0, 0, 0};
		if(electric_) Ek = electric_->contravariant(qk, t * iEfield_time_factor_);
		if(magnetic_) Bk = magnetic_->contravariant(qk, t * iBfield_time_factor_);

		// perform cartesian step
		IR3 vph = cartesian_boris(umh, Oref_, Ek, Bk, dt);

		IR3 qp1 = qk + dt * vph;

		return {qp1[IR3::u], qp1[IR3::v], qp1[IR3::w],
				vph[IR3::u], vph[IR3::v], vph[IR3::w]};

	}

}

// Returns the vector position of the state normalized to `Lref`.
IR3 boris_advanced::get_position(const boris_advanced::state& s) const {
	return {s[0], s[1], s[2]};
}

// Returns the vector velocity of the state normalized to `Vref`.
IR3 boris_advanced::get_velocity(const boris_advanced::state& s) const {
	return {s[3], s[4], s[5]};
}

// Returns the kinetic energy of the state, normalized to `Uref`.
double boris_advanced::get_kinetic_energy(const boris_advanced::state& s) const {
	if(field_morph_) {

		IR3 qk = {s[0], s[1], s[2]};
		IR3 umh = {s[3], s[4], s[5]};
		IR3 vmh = field_morph_->from_contravariant(qk, umh);

		return vmh[IR3::u] * vmh[IR3::u] + vmh[IR3::v] * vmh[IR3::v] + vmh[IR3::w] * vmh[IR3::w];

	} else return s[3]*s[3] + s[4]*s[4] + s[5]*s[5];
}

// Returns the parallel energy of the state, normalized to `Uref`.
double boris_advanced::get_parallel_energy(const boris_advanced::state& s, double &time) const {
	IR3 q = {s[0], s[1], s[2]};
	IR3 u = {s[3], s[4], s[5]};
	if(magnetic_) {
		IR3 b = magnetic_->contravariant_versor(q, time * iBfield_time_factor_);
		if(field_morph_) {
			IR3 v = field_morph_->from_contravariant(q, u);
			IR3 u_cov = field_morph_->to_covariant(q, v);
			double vpp = inner_product(u_cov, b);
			return vpp * vpp;
		} else {
			double vpp = inner_product(u, b);
			return vpp * vpp;
		}
	} else {
		error(__func__, __FILE__, __LINE__, "null magnetic field.", 1);
		return 0;
	}
}

// Returns the perpendicular energy of the state, normalized to `Uref`.
double boris_advanced::get_perpendicular_energy(const boris_advanced::state& s, double &time) const {
	IR3 q = {s[0], s[1], s[2]};
	IR3 u = {s[3], s[4], s[5]};
	if(magnetic_) {
		IR3 b = magnetic_->contravariant_versor(q, time * iBfield_time_factor_);
		if(field_morph_) {
			IR3 v = field_morph_->from_contravariant(q, u);
			IR3 b_cartesian = field_morph_->from_contravariant(q, b);
			IR3 vperp = cartesian_cross_product(v, b_cartesian);
			return inner_product(vperp, vperp);
		} else {
			IR3 vperp = cartesian_cross_product(u, b);
			return inner_product(vperp, vperp);
		}
	} else {
		error(__func__, __FILE__, __LINE__, "null magnetic field.", 1);
		return 0;
	}
}

// Returns the `boris_advanced::state` from a normalized point in phase-space.
boris_advanced::state boris_advanced::generate_state(const IR3 &pos, const IR3 &vel) const {
	return {pos[IR3::u], pos[IR3::v], pos[IR3::w],
			vel[IR3::u], vel[IR3::v], vel[IR3::w]};
}

// Creates the first `boris_advanced::state` from a normalized point in cartesian phase-space.
boris_advanced::state boris_advanced::generate_initial_state(
		const IR3 &cartesian_position, const IR3 &cartesian_velocity, const double &tinit, const double &dt) const {

	electromagnetic_system em(1, 1, qom_, electric_, magnetic_, field_morph_);
	odeint_adapter<electromagnetic_system> sys(&em);
	boost::numeric::odeint::runge_kutta4<electromagnetic_system::state> rk4;

	IR3 qk = cartesian_position;
	IR3 uk = cartesian_velocity;
	if(field_morph_) {
		qk = field_morph_->inverse(cartesian_position);
		uk = field_morph_->to_contravariant(qk, cartesian_velocity);
	}
	electromagnetic_system::state in = em.generate_state(qk, uk);
	electromagnetic_system::state out;

	rk4.do_step(sys, in, tinit, out, -0.5*dt);
	IR3 qmh = em.get_position(out);
	IR3 umh = em.get_velocity(out);

	if(field_morph_) {
		IR3 vmh = field_morph_->from_contravariant(qmh, umh);
		IR3 uumh = field_morph_->to_contravariant(qk, vmh);

		return generate_state(qk, uumh);

	} else return generate_state(qk, umh);
}

// Performs a boris step in cartesian coordinates.
IR3 boris_advanced::cartesian_boris(const IR3 &v_old, const double &Oref, const IR3 &Ek, const IR3 &Bk, const double &dt) const {

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

	return v_new;
}

} // end namespace gyronimo