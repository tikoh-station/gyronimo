#include <gyronimo/dynamics/exact_stepper.hh>

namespace gyronimo {

// Performs a time step `dt` update to the state `in` and returns the result.
exact_stepper::state exact_stepper::do_step(const exact_stepper::state &in, const double t, const double dt) const {

	// extract position and velocity from state
	IR3 q_old = get_position(in);
	IR3 v_old = get_velocity(in);

	// advance particle across phase-space

	dIR3 ek = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	if(field_morph_) {
		ek = field_morph_->del(q_old);
		v_old = contraction<second>(ek, v_old);
	}

	auto [dx1, v1] = displacement(q_old, v_old, t, 0.5 * dt);
	double t_field = t + 0.5 * dt;
	IR3 q_field = straight_line_step(q_old, dx1);

	dIR3 e_field = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	if(field_morph_) e_field = field_morph_->del(q_field);

	IR3 v_new = v_old;
	IR3 dx = {0, 0, 0};

	IR3 Ek = {0, 0, 0};
	IR3 Bk = {0, 0, 0};

	double B2 = 0;
	if(magnetic_) {
		Bk = magnetic_->contravariant(q_field, t_field * iBfield_time_factor_);
		if(field_morph_) Bk = contraction<second>(e_field, Bk);
		B2 = Bk[IR3::u] * Bk[IR3::u] + Bk[IR3::v] * Bk[IR3::v] + Bk[IR3::w] * Bk[IR3::w];
	}

	if(B2 != 0) {

		IR3 vExB = {0, 0, 0};
		double Epar = 0;

		double B = sqrt(B2);
		double iB = 1 / B;

		double tn = tan(- 0.5 * Oref_ * B * dt);
		double idenominator = 2 / (1 + tn * tn);
		double sn = tn * idenominator;
		double cn = idenominator - 1;

		IR3 bk = iB * Bk;
		double vpar = inner_product(v_old, bk);

		if(electric_) {

			Ek = electric_->contravariant(q_field, t_field * iEfield_time_factor_);
			if(field_morph_) Ek = contraction<second>(e_field, Ek);
			vExB = iB * gyronimo::cartesian_cross_product(Ek, bk);
			Epar = inner_product(Ek, bk);
			v_old -= vExB;

		}

		double iOmega = iB * iOref_;
		IR3 vprime = gyronimo::cartesian_cross_product(bk, v_old);

		v_new = cn * v_old + sn * vprime + ((1 - cn) * vpar) * bk;
		dx = iOmega * ((cn - 1) * vprime + sn * (vpar * bk - v_old)) + (vpar * dt) * bk;

		if(electric_) {

			IR3 electric_impulse = (Oref_ * Epar * dt) * bk;
			v_new += vExB + electric_impulse;
			dx += dt * vExB + (0.5 * dt) * electric_impulse;

		}

	} else {

		if(electric_) {

			Ek = electric_->contravariant(q_field, t_field * iEfield_time_factor_);
			if(field_morph_) Ek = contraction<second>(e_field, Ek);
			IR3 electric_impulse = (Oref_ * dt) * Ek;
			v_new += electric_impulse;
			dx = (v_old + 0.5 * electric_impulse) * dt;

		} else dx = v_old * dt;

	}


	IR3 q_new = straight_line_step(q_old, dx);
	if(field_morph_) {
		IR3 u_new = field_morph_->to_contravariant(q_new, v_new);

		return {q_new[IR3::u], q_new[IR3::v], q_new[IR3::w], 
				u_new[IR3::u], u_new[IR3::v], u_new[IR3::w]};
	}

	return {q_new[IR3::u], q_new[IR3::v], q_new[IR3::w], 
			v_new[IR3::u], v_new[IR3::v], v_new[IR3::w]};
}

// Returns the vector position of the state normalized to `Lref`.
IR3 exact_stepper::get_position(const exact_stepper::state& s) const {
	return {s[0], s[1], s[2]};
}

// Returns the vector velocity of the state normalized to `Vref`.
IR3 exact_stepper::get_velocity(const exact_stepper::state& s) const {
	return {s[3], s[4], s[5]};
}

// Returns the kinetic energy of the state, normalized to `Uref`.
double exact_stepper::get_kinetic_energy(const exact_stepper::state& s) const {
	IR3 q = {s[0], s[1], s[2]};
	IR3 u = {s[3], s[4], s[5]};
	IR3 v = field_morph_->from_contravariant(q, u);
	return v[IR3::u] * v[IR3::u] + v[IR3::v] * v[IR3::v] + v[IR3::w] * v[IR3::w];
}

// Returns the `exact_stepper::state` from a normalized point in phase-space.
exact_stepper::state exact_stepper::generate_state(const IR3 &pos, const IR3 &vel) const {
	return {pos[IR3::u], pos[IR3::v], pos[IR3::w],
			vel[IR3::u], vel[IR3::v], vel[IR3::w]};
}

// Creates the first `exact_stepper::state` from a normalized point in cartesian phase-space.
exact_stepper::state exact_stepper::generate_initial_state(const IR3 &cartesian_position, 
		const IR3 &cartesian_velocity, const double &tinit, const double &dt) const {
	
	if(field_morph_) {

		gyronimo::IR3 q0 = field_morph_->inverse(cartesian_position);
		gyronimo::IR3 u0 = field_morph_->to_contravariant(q0, cartesian_velocity);
		return generate_state(q0, u0);

	} else return generate_state(cartesian_position, cartesian_velocity);
}

// Performs a `runge-kutta` step in curvilinear coordinates across the line set by the vector dx.
IR3 exact_stepper::straight_line_step(const IR3 &q_old, const IR3 &dx) const {

	IR3 q_new = {0, 0, 0};
	if(field_morph_) {

		exact_stepper::straight_line_system sys(field_morph_, dx);
		exact_stepper::straight_line_system::odeint_stepper rk;

		exact_stepper::straight_line_system::state st_old = sys.generate_state(q_old);
		exact_stepper::straight_line_system::state st_new;
		rk.do_step(sys, st_old, 0, st_new, 1);
		q_new = sys.get_IR3(st_new);

	} else q_new = q_old + dx;

	return q_new;
}

// Calculates the cartesian displacement vector.
std::pair<IR3, IR3> exact_stepper::displacement(const IR3 &q_field, IR3 v_old, const double &t_field, const double &dt) const {

	// advance particle across phase-space

	dIR3 e_field = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	if(field_morph_) e_field = field_morph_->del(q_field);

	IR3 v_new = v_old;
	IR3 dx = {0, 0, 0};

	IR3 Ek = {0, 0, 0};
	IR3 Bk = {0, 0, 0};

	double B2 = 0;
	if(magnetic_) {
		Bk = magnetic_->contravariant(q_field, t_field * iBfield_time_factor_);
		if(field_morph_) Bk = contraction<second>(e_field, Bk);
		B2 = Bk[IR3::u] * Bk[IR3::u] + Bk[IR3::v] * Bk[IR3::v] + Bk[IR3::w] * Bk[IR3::w];
	}

	if(B2 != 0) {

		IR3 vExB = {0, 0, 0};
		double Epar = 0;

		double B = sqrt(B2);
		double iB = 1 / B;

		double tn = tan(- 0.5 * Oref_ * B * dt);
		double idenominator = 2 / (1 + tn * tn);
		double sn = tn * idenominator;
		double cn = idenominator - 1;

		IR3 bk = iB * Bk;
		double vpar = inner_product(v_old, bk);

		if(electric_) {

			Ek = electric_->contravariant(q_field, t_field * iEfield_time_factor_);
			if(field_morph_) Ek = contraction<second>(e_field, Ek);
			vExB = iB * gyronimo::cartesian_cross_product(Ek, bk);
			Epar = inner_product(Ek, bk);
			v_old -= vExB;

		}

		double iOmega = iB * iOref_;
		IR3 vprime = gyronimo::cartesian_cross_product(bk, v_old);

		v_new = cn * v_old + sn * vprime + ((1 - cn) * vpar) * bk;
		dx = iOmega * ((cn - 1) * vprime + sn * (vpar * bk - v_old)) + (vpar * dt) * bk;

		if(electric_) {

			IR3 electric_impulse = (Oref_ * Epar * dt) * bk;
			v_new += vExB + electric_impulse;
			dx += dt * vExB + (0.5 * dt) * electric_impulse;

		}

	} else {

		if(electric_) {

			Ek = electric_->contravariant(q_field, t_field * iEfield_time_factor_);
			if(field_morph_) Ek = contraction<second>(e_field, Ek);
			IR3 electric_impulse = (Oref_ * dt) * Ek;
			v_new += electric_impulse;
			dx = (v_old + 0.5 * electric_impulse) * dt;

		} else dx = v_old * dt;

	}

	return {dx, v_new};
}

} // end namespace gyronimo