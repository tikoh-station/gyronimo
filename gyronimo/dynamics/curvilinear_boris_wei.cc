#include <gyronimo/dynamics/curvilinear_boris_wei.hh>

namespace gyronimo {

// Performs a time step `dt` update to the state `in` and returns the result.
curvilinear_boris_wei::state curvilinear_boris_wei::do_step(const curvilinear_boris_wei::state &in, const double t, const double dt) const {

	// extract position and velocity from state
	IR3 qm1 = {in[0], in[1], in[2]};
	IR3 umh = {in[3], in[4], in[5]};
	IR3 qmh = {in[6], in[7], in[8]};

	if(field_morph_) {

		IR3 qk = qm1 + umh * dt;
		dIR3 ek  = field_morph_->del(qk);
		dIR3 eek = inverse(ek);

		IR3 vmh = field_morph_->from_contravariant(qmh, umh);

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

		IR3 uuph = contraction<second>(eek, vph);
		IR3 qph = qk + (0.5 * dt) * uuph;
		dIR3 eeph = field_morph_->del_inverse(qph);
		IR3 uph = contraction<second>(eeph, vph);

		return { qk[IR3::u],  qk[IR3::v],  qk[IR3::w],
				uph[IR3::u], uph[IR3::v], uph[IR3::w],
				qph[IR3::u], qph[IR3::v], qph[IR3::w]};

	} else {

		IR3 xk = qm1 + umh * dt;

		// calculate fields
		IR3 Ek = {0, 0, 0};
		IR3 Bk = {0, 0, 0};
		if(electric_) Ek = electric_->contravariant(xk, t * iEfield_time_factor_);
		if(magnetic_) Bk = magnetic_->contravariant(xk, t * iBfield_time_factor_);

		// perform cartesian step
		IR3 vph = cartesian_boris(umh, Oref_, Ek, Bk, dt);

		IR3 xph = xk + (0.5 * dt) * vph;

		return { xk[IR3::u],  xk[IR3::v],  xk[IR3::w],
				vph[IR3::u], vph[IR3::v], vph[IR3::w],
				xph[IR3::u], xph[IR3::v], xph[IR3::w]};

	}

}

// Returns the vector position of the state normalized to `Lref`.
IR3 curvilinear_boris_wei::get_position(const curvilinear_boris_wei::state& s) {
	return {s[0], s[1], s[2]};
}

// Returns the vector velocity of the state normalized to `Vref`.
IR3 curvilinear_boris_wei::get_velocity(const curvilinear_boris_wei::state& s) {
	return {s[3], s[4], s[5]};
}

// Returns the kinetic energy of the state, normalized to `Uref`.
double curvilinear_boris_wei::get_kinetic_energy(const curvilinear_boris_wei::state& s) {
	if(field_morph_) {

		IR3 qmh = {s[6], s[7], s[8]};
		IR3 umh = {s[3], s[4], s[5]};
		dIR3 emh = field_morph_->del(qmh);
		IR3 vmh = contraction<second>(emh, umh);

		return vmh[IR3::u] * vmh[IR3::u] + vmh[IR3::v] * vmh[IR3::v] + vmh[IR3::w] * vmh[IR3::w];

	} else return s[3]*s[3] + s[4]*s[4] + s[5]*s[5];
}

// Returns the `curvilinear_boris_wei::state` from a normalized point in phase-space.
curvilinear_boris_wei::state curvilinear_boris_wei::generate_state(const IR3 &pos, const IR3 &vel, const IR3 &pph) const {
	return {pos[IR3::u], pos[IR3::v], pos[IR3::w],
			vel[IR3::u], vel[IR3::v], vel[IR3::w],
			pph[IR3::u], pph[IR3::v], pph[IR3::w]};
}

// Creates the first `curvilinear_boris_wei::state` from a normalized point in cartesian phase-space.
curvilinear_boris_wei::state curvilinear_boris_wei::generate_initial_state(
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

	rk4.do_step(sys, in, tinit, out, -dt);
	IR3 qm1 = em.get_position(out);

	return generate_state(qm1, umh, qmh);
}

// Performs a boris step in cartesian coordinates.
IR3 curvilinear_boris_wei::cartesian_boris(const IR3 &v_old, const double &Oref, const IR3 &Ek, const IR3 &Bk, const double &dt) const {

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