#include <gyronimo/dynamics/centered_stepper.hh>

namespace gyronimo {

// Performs a time step `dt` update to the state `in` and returns the result.
centered_stepper::state centered_stepper::do_step(const centered_stepper::state &in, const double t, const double dt) const {

	// extract position and velocity from state
	IR3 q_old = {in[0], in[1], in[2]};
	IR3 u_old = {in[3], in[4], in[5]};
	IR3 v_old = {in[6], in[7], in[8]};

	// solve in cartesian coordinates
	if(field_morph_) {

		// coordinate conversion to cartesian
		dIR3 ek   = field_morph_->del(q_old);
		dIR3 em1  = field_morph_->del(q_old - u_old * dt);
		dIR3 eek  = inverse(ek);
		dIR3 eem1 = inverse(em1);

		// IR3 v_old = 0.5 * (contraction<second>(ek, u_old) + contraction<second>(em1, u_old));

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
		IR3 v_new = cartesian_boris(v_old, Oref_, Ek, Bk, dt);

		// std::cout << "v_old: [" << v_old[0] << " " << v_old[1] << " " << v_old[2] << " ]" << std::endl;
		// std::cout << "v_new: [" << v_new[0] << " " << v_new[1] << " " << v_new[2] << " ]" << std::endl;

		// conversion back to curvilinear coordinates
		IR3 u_new = 1.5 * contraction<second>(eek, v_new) + (-0.5) * contraction<second>(eem1, v_new);
		IR3 q_new = q_old + u_new * dt;

		return {q_new[IR3::u], q_new[IR3::v], q_new[IR3::w],
				u_new[IR3::u], u_new[IR3::v], u_new[IR3::w],
				v_new[IR3::u], v_new[IR3::v], v_new[IR3::w]};

	} else {

		// calculate fields
		IR3 Ek = {0, 0, 0};
		IR3 Bk = {0, 0, 0};
		if(electric_) Ek = electric_->contravariant(q_old, t * iEfield_time_factor_);
		if(magnetic_) Bk = magnetic_->contravariant(q_old, t * iBfield_time_factor_);

		// perform cartesian step
		IR3 v_new = cartesian_boris(u_old, Oref_, Ek, Bk, dt);
		IR3 x_new = q_old + v_new * dt;

		return {x_new[IR3::u], x_new[IR3::v], x_new[IR3::w],
				v_new[IR3::u], v_new[IR3::v], v_new[IR3::w]};
	}

}

// Returns the vector position of the state normalized to `Lref`.
IR3 centered_stepper::get_position(const centered_stepper::state& s) {
	return {s[0], s[1], s[2]};
}

// Returns the vector velocity of the state normalized to `Vref`.
IR3 centered_stepper::get_velocity(const centered_stepper::state& s) {
	return {s[3], s[4], s[5]};
}

// Returns the kinetic energy of the state, normalized to `Uref`.
double centered_stepper::get_kinetic_energy(const centered_stepper::state& s) {
	return s[6]*s[6] + s[7]*s[7] + s[8]*s[8];
}

// Returns the `centered_stepper::state` from a normalized point in phase-space.
centered_stepper::state centered_stepper::generate_state(const IR3 &pos, const IR3 &vel, const IR3 &vct) const {
	return {pos[IR3::u], pos[IR3::v], pos[IR3::w],
			vel[IR3::u], vel[IR3::v], vel[IR3::w],
			vct[IR3::u], vct[IR3::v], vct[IR3::w]};
}

// Creates the first `centered_stepper::state` from a normalized point in cartesian phase-space.
centered_stepper::state centered_stepper::generate_initial_state(const IR3 &cartesian_position, const IR3 &cartesian_velocity, const double &tinit, const double &dt) const {

	electromagnetic_system em(1, 1, qom_, electric_, magnetic_, field_morph_);
	odeint_adapter<electromagnetic_system> sys(&em);
	boost::numeric::odeint::runge_kutta4<electromagnetic_system::state> rk4;

	electromagnetic_system::state in = em.generate_state(cartesian_position, cartesian_velocity);
	electromagnetic_system::state out;

	rk4.do_step(sys, in, tinit, out, -0.5*dt);
	IR3 xmh = em.get_position(out);
	IR3 vmh = em.get_velocity(out);

	if(field_morph_) {

		IR3 qk = field_morph_->inverse(cartesian_position);
		IR3 qmh = field_morph_->inverse(xmh);
		IR3 umh = field_morph_->to_contravariant(qmh, vmh);

		return generate_state(qk, umh, vmh);

	} else return generate_state(cartesian_position, vmh, vmh);

	// if(field_morph_) {

	// 	// coordinate conversion to cartesian
	// 	IR3 q_old = field_morph_->inverse(cartesian_position);
	// 	dIR3 ek = field_morph_->del(q_old);

	// 	// calculate fields
	// 	IR3 Ek = {0, 0, 0};
	// 	IR3 Bk = {0, 0, 0};
	// 	if(electric_) {
	// 		Ek = electric_->contravariant(q_old, tinit * iEfield_time_factor_);
	// 		Ek = contraction<second>(ek, Ek);
	// 	}
	// 	if(magnetic_) {
	// 		Bk = magnetic_->contravariant(q_old, tinit * iBfield_time_factor_);
	// 		Bk = contraction<second>(ek, Bk);
	// 	}

	// 	IR3 vmh = cartesian_boris(cartesian_velocity, Oref_, Ek, Bk, (-0.5 * dt));

	// 	IR3 xm1 = cartesian_position - dt * vmh;
	// 	IR3 qm1 = field_morph_->inverse(xm1);
	// 	dIR3 eem1 = field_morph_->del_inverse(qm1);
	// 	dIR3 eek = inverse(ek);
	// 	IR3 u_old = 0.5 * (contraction<second>(eek, vmh) + contraction<second>(eem1, vmh));

	// 	return generate_state(q_old, u_old, vmh);

	// } else {

	// 	// calculate fields
	// 	IR3 Ek = {0, 0, 0};
	// 	IR3 Bk = {0, 0, 0};
	// 	if(electric_) Ek = electric_->contravariant(cartesian_position, tinit * iEfield_time_factor_);
	// 	if(magnetic_) Bk = magnetic_->contravariant(cartesian_position, tinit * iBfield_time_factor_);

	// 	IR3 vmh = cartesian_boris(cartesian_velocity, Oref_, Ek, Bk, -0.5*dt);

	// 	return generate_state(cartesian_position, vmh, vmh);
	// }

}

// Performs a boris step in cartesian coordinates.
IR3 centered_stepper::cartesian_boris(const IR3 &v_old, const double &Oref, const IR3 &Ek, const IR3 &Bk, const double &dt) const {


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