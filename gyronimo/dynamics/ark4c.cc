#include <gyronimo/dynamics/ark4c.hh>

namespace gyronimo {

// Class Constructor.
ark4c::ark4c(const double Lref, const double Vref, const double qom, const IR3field *E, const IR3field *B, const morphism *morph)
		: Lref_(Lref), Vref_(Vref), Tref_(Lref/Vref), qom_(qom),
		Oref_(	B ? qom * gyronimo::codata::e * Tref_ / gyronimo::codata::m_proton * B->m_factor() : 
				E ? qom * gyronimo::codata::e * Tref_ / gyronimo::codata::m_proton * E->m_factor() : 0), 
		iOref_(	B ? 1 / Oref_ : E ? 1 / Oref_ : 0), electric_(E), magnetic_(B),
		iBfield_time_factor_(B ? Tref_ / B->t_factor() : 0.0),
		iEfield_time_factor_(E ? Tref_ / E->t_factor() : 0.0), 
		field_morph_(morph), 
		em_(Lref, Vref, qom, E, B, morph), sys_(&em_),
		a0_((2. + std::pow(2., -1./3.) + std::pow(2., 1./3.)) / 3.), b0_(1. - 2. * a0_), c0_(a0_),
		alpha_(0.5 * a0_), beta_(0.5), gamma_(1. - alpha_),
		k0_({alpha_, beta_, gamma_}), 
		k1_((4 - 3 * k0_) * k0_ * k0_ * k0_),
		k2_((k0_ - 1) * k0_ * (2 * k0_ * k0_ - k0_ - 1)), 
		k3_((k0_ - 1) * k0_ * k0_ * k0_), 
		k4_((0.5 * (k0_ - 1) * (k0_ - 1) * k0_ * k0_)) {
		// a0_(7./24.), b0_(0.75), c0_(-1./24.),
		// alpha_(0), beta_(2./3.), gamma_(0) {

	//@todo test if fields are consistent (normalization and cartesian metric)

	// std::cout << "ark4c parameter list:" << std::endl;
	// std::cout << "a0: " << a0_ << "\tb0: " << b0_ << "\tc0: " << c0_ << std::endl;
	// std::cout << "alpha: " << alpha_ << "\tbeta: " << beta_ << "\tgamma: " << gamma_ << std::endl;
	// std::cout << std::endl;
}

// Performs a time step `dt` update to the state `in` and returns the result.
ark4c::state ark4c::do_step(const ark4c::state &in, const double t, const double dt) const {

	// extract position and velocity from state
	IR3 q0 = get_position(in);
	IR3 u0 = get_velocity(in);

	// estimate field positions
	boost::numeric::odeint::runge_kutta4<electromagnetic_system::state> rk4;

	electromagnetic_system::state st0 = em_.generate_state(q0, u0);
	electromagnetic_system::state dst0 = em_(st0, t);
	IR3 du0 = em_.get_velocity(dst0);

	electromagnetic_system::state st1;
	rk4.do_step(sys_, st0, t, st1, dt);
	IR3 q1 = em_.get_position(st1);
	IR3 u1 = em_.get_velocity(st1);

	// create the interpolator
	IR3 qA = interpolate(q0, q1, u0, u1, du0, IR3::u, dt);
	IR3 qB = interpolate(q0, q1, u0, u1, du0, IR3::v, dt);
	IR3 qC = interpolate(q0, q1, u0, u1, du0, IR3::w, dt);

	// advance particle across phase-space
	IR3 v0 = u0;
	if(field_morph_) v0 = field_morph_->from_contravariant(q0, u0);

	auto [dxA, vA] = displacement(qA, v0, t + alpha_ * dt, a0_ * dt);
	auto [dxB, vB] = displacement(qB, vA, t + beta_  * dt, b0_ * dt);
	auto [dxC, vC] = displacement(qC, vB, t + gamma_ * dt, c0_ * dt);

	if(field_morph_) {

		IR3 q_new = straight_line_step(q0, dxA + dxB + dxC);
		IR3 u_new = field_morph_->to_contravariant(q_new, vC);

		return {q_new[IR3::u], q_new[IR3::v], q_new[IR3::w], 
				u_new[IR3::u], u_new[IR3::v], u_new[IR3::w]};

	} else {

		IR3 x_new = q0 + dxA + dxB + dxC;
		IR3 v_new = vC;

		return {x_new[IR3::u], x_new[IR3::v], x_new[IR3::w], 
				v_new[IR3::u], v_new[IR3::v], v_new[IR3::w]};
	}

	// // THIRD ORDER CONFIGURATION WITH RUNGE KUTTA 4

	// electromagnetic_system::state st1;
	// electromagnetic_system::state st2;
	// electromagnetic_system::state st3;

	// rk4.do_step(sys, st0, t, st1, alpha_ * dt);
	// rk4.do_step(sys, st0, t, st2, beta_  * dt);
	// rk4.do_step(sys, st0, t, st3, gamma_ * dt);

	// IR3 q1 = em.get_position(st1);
	// IR3 q2 = em.get_position(st2);
	// IR3 q3 = em.get_position(st3);

	// // advance particle across phase-space

	// dIR3 ek = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	// IR3 v0 = u0;
	// if(field_morph_) v0 = field_morph_->from_contravariant(q0, u0);

	// auto [dx1, v1] = displacement(q1, v0, t + alpha_ * dt, a0_ * dt);
	// auto [dx2, v2] = displacement(q2, v1, t + beta_  * dt, b0_ * dt);
	// auto [dx3, v3] = displacement(q3, v2, t + gamma_ * dt, c0_ * dt);

	// if(field_morph_) {

	// 	IR3 q_new = straight_line_step(q0, dx1 + dx2 + dx3);
	// 	IR3 u_new = field_morph_->to_contravariant(q_new, v3);

	// 	return {q_new[IR3::u], q_new[IR3::v], q_new[IR3::w], 
	// 			u_new[IR3::u], u_new[IR3::v], u_new[IR3::w]};

	// } else {

	// 	IR3 x_new = q0 + dx1 + dx2 + dx3;
	// 	IR3 v_new = v3;

	// 	return {x_new[IR3::u], x_new[IR3::v], x_new[IR3::w], 
	// 			v_new[IR3::u], v_new[IR3::v], v_new[IR3::w]};
	// }
	// //
}

// Returns the vector position of the state normalized to `Lref`.
IR3 ark4c::get_position(const ark4c::state& s) const {
	return {s[0], s[1], s[2]};
}

// Returns the vector velocity of the state normalized to `Vref`.
IR3 ark4c::get_velocity(const ark4c::state& s) const {
	return {s[3], s[4], s[5]};
}

// Returns the kinetic energy of the state, normalized to `Uref`.
double ark4c::get_kinetic_energy(const ark4c::state& s) const {
	IR3 q = {s[0], s[1], s[2]};
	IR3 u = {s[3], s[4], s[5]};
	if(field_morph_) {
		IR3 v = field_morph_->from_contravariant(q, u);
		return v[IR3::u] * v[IR3::u] + v[IR3::v] * v[IR3::v] + v[IR3::w] * v[IR3::w];
	} else {
		return u[IR3::u] * u[IR3::u] + u[IR3::v] * u[IR3::v] + u[IR3::w] * u[IR3::w];
	}
}

// Returns the parallel energy of the state, normalized to `Uref`.
double ark4c::get_parallel_energy(const ark4c::state& s, double &time) const {
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
double ark4c::get_perpendicular_energy(const ark4c::state& s, double &time) const {
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

// Returns the `ark4c::state` from a normalized point in phase-space.
ark4c::state ark4c::generate_state(const IR3 &pos, const IR3 &vel) const {
	return {pos[IR3::u], pos[IR3::v], pos[IR3::w],
			vel[IR3::u], vel[IR3::v], vel[IR3::w]};
}

// Creates the first `ark4c::state` from a normalized point in cartesian phase-space.
ark4c::state ark4c::generate_initial_state(const IR3 &cartesian_position, 
		const IR3 &cartesian_velocity, const double &tinit, const double &dt) const {
	
	if(field_morph_) {

		gyronimo::IR3 q0 = field_morph_->inverse(cartesian_position);
		gyronimo::IR3 u0 = field_morph_->to_contravariant(q0, cartesian_velocity);
		return generate_state(q0, u0);

	} else return generate_state(cartesian_position, cartesian_velocity);
}

// Performs a `runge-kutta` step in curvilinear coordinates across the line set by the vector dx.
IR3 ark4c::straight_line_step(const IR3 &q_old, const IR3 &dx) const {

	IR3 q_new = {0, 0, 0};
	if(field_morph_) {

		ark4c::straight_line_system sys(field_morph_, dx);
		ark4c::straight_line_system::odeint_stepper rk;

		ark4c::straight_line_system::state st_old = sys.generate_state(q_old);
		ark4c::straight_line_system::state st_new;
		rk.do_step(sys, st_old, 0, st_new, 1);
		q_new = sys.get_IR3(st_new);

	} else q_new = q_old + dx;

	return q_new;
}

// Calculates the cartesian displacement vector.
std::pair<IR3, IR3> ark4c::displacement(const IR3 &q_field, IR3 v_old, const double &t_field, const double &dt) const {

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

// Interpolate at time t \element [0, 1]
IR3 ark4c::interpolate(IR3 q0, IR3 q1, IR3 u0, IR3 u1, IR3 du0, IR3::index i, double dt) const {
	// double a0 = (4 - 3 * k) * k * k * k;
	// double a1 = dt * (k - 1) * k;
	// double a2 = (2 * k * k - k - 1);
	// double a3 = 0.5 * a1 * a1;

	// IR3 q = (1 - a0) * q0 + a0 * q1 + (a1 * a2) * u0 + (a1 * k * k) * u1 + a3 * du0;
	// return q;

	IR3 q = (1 - k1_[i]) * q0 + k1_[i] * q1 + (dt * k2_[i]) * u0 + (dt * k3_[i]) * u1 + (dt * dt * k4_[i]) * du0;
	return q;
}

} // end namespace gyronimo