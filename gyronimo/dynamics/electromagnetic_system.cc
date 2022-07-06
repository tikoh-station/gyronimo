#include <gyronimo/dynamics/electromagnetic_system.hh>

namespace gyronimo {

electromagnetic_system::electromagnetic_system(double Lref, double Vref, 
			double qom, const IR3field* E, const IR3field* B, 
			const morphism *morph/*, const metric_covariant *metric*/) : 
		Lref_(Lref), Vref_(Vref), Tref_(Lref / Vref),
		Oref_(  B ? qom * codata::e / codata::m_proton * B->m_factor() :
				E ? qom * codata::e / codata::m_proton * E->m_factor() : 0),
		electric_field_(E), magnetic_field_(B), 
		iEfield_time_factor_(E ? Tref_ / E->t_factor() : 0), 
		iBfield_time_factor_(B ? Tref_ / B->t_factor() : 0),
		field_morph_(morph)/*, metric_(metric)*/ {

	// test if fields are consistent (normalization and metrics)
}

electromagnetic_system::state electromagnetic_system::operator()(const state &x, double t) const {
	/*
	IR3 pos = {x[0], x[1], x[2]};
	IR3 vel = {x[3], x[4], x[5]};

	IR3 E = {0, 0, 0};
	IR3 B = {0, 0, 0};

	if(field_morph_) {

		IR3 q = pos;
		dIR3 ek = {1, 0, 0, 0, 1, 0, 0, 0, 1};
		if(field_morph_) {
			q = field_morph_->inverse(q);
			ek = field_morph_->del(q);
		}

		if(electric_field_) {
			E = electric_field_->contravariant(q, t * iEfield_time_factor_);
			E = contraction<second>(ek, E);
		}
		if(magnetic_field_) {
			B = magnetic_field_->contravariant(q, t * iBfield_time_factor_);
			B = contraction<second>(ek, B);
		}

		IR3 acl = Oref_ * (E + cartesian_cross_product(vel, B));

		return {vel[IR3::u], vel[IR3::v], vel[IR3::w],
				acl[IR3::u], acl[IR3::v], acl[IR3::w]};

	} else {

		if(electric_field_) E = electric_field_->contravariant(pos, t * iEfield_time_factor_);
		if(magnetic_field_) B = magnetic_field_->contravariant(pos, t * iBfield_time_factor_);

		IR3 acl = Oref_ * (E + cartesian_cross_product(vel, B));

		return {vel[IR3::u], vel[IR3::v], vel[IR3::w],
				acl[IR3::u], acl[IR3::v], acl[IR3::w]};
	}
	*/

	IR3 q = {x[0], x[1], x[2]};
	IR3 u = {x[3], x[4], x[5]};

	IR3 E = {0, 0, 0};
	IR3 B = {0, 0, 0};

	if(field_morph_) {

		dIR3 e = field_morph_->del(q);
		IR3 e1 = {e[dIR3::uu], e[dIR3::vu], e[dIR3::wu]};
		IR3 e2 = {e[dIR3::uv], e[dIR3::vv], e[dIR3::wv]};
		IR3 e3 = {e[dIR3::uw], e[dIR3::vw], e[dIR3::ww]};
		SM3 g = {inner_product(e1, e1), inner_product(e1, e2), inner_product(e1, e3),
				 inner_product(e2, e2), inner_product(e2, e3), inner_product(e3, e3)};
		SM3 ig = inverse(g);
		dSM3 dg = field_morph_->g_del(q);

		if(electric_field_) E = electric_field_->contravariant(q, t * iEfield_time_factor_);
		if(magnetic_field_) B = magnetic_field_->contravariant(q, t * iBfield_time_factor_);

		/* INERTIAL FORCE TERM

			G^{m}_{jk} = 1/2 g^{mi} (d_{k} g_{ij} + d_{j} g_{ik} - d_{i} g_{jk})

			G^{m}_{jk} * u^{j} * u{k} = 
				= 1/2 g^{mi} (d_{k} g_{ij} * u^{j} * u{k} 
							+ d_{j} g_{ik} * u^{j} * u{k} 
							- d_{i} g_{jk} * u^{j} * u{k})
		*/

		dIR3 chr1 = contraction<second>(dg, u);
		dIR3 chr2 = contraction<third>(dg, u);
		dIR3 chr3 = contraction<first>(dg, u);

		IR3 chr_1 = contraction<second>(chr1, u);
		IR3 chr_2 = contraction<second>(chr2, u);
		IR3 chr_3 = contraction<first>(chr3, u);

		IR3 inertial_covariant = 0.5 * (chr_1 + chr_2 - chr_3);

		double jacobian = field_morph_->jacobian(q);
		IR3 vxB_covariant = (Oref_ * jacobian) * cartesian_cross_product(u, B);

		IR3 magnetic_forces_contravariant = contraction(ig, vxB_covariant - inertial_covariant);

		IR3 a = Oref_ * E + magnetic_forces_contravariant;

		return {u[IR3::u], u[IR3::v], u[IR3::w],
				a[IR3::u], a[IR3::v], a[IR3::w]};

	} else {

		if(electric_field_) E = electric_field_->contravariant(q, t * iEfield_time_factor_);
		if(magnetic_field_) B = magnetic_field_->contravariant(q, t * iBfield_time_factor_);

		IR3 a = Oref_ * (E + cartesian_cross_product(u, B));

		return {u[IR3::u], u[IR3::v], u[IR3::w],
				a[IR3::u], a[IR3::v], a[IR3::w]};
	}
}

IR3 electromagnetic_system::get_position(const state &s) const {
	return {s[0], s[1], s[2]};
}

IR3 electromagnetic_system::get_velocity(const state &s) const {
	return {s[3], s[4], s[5]};
}

electromagnetic_system::state electromagnetic_system::generate_state(const IR3 &pos, const IR3 &vel) const {
	return {pos[IR3::u], pos[IR3::v], pos[IR3::w],
			vel[IR3::u], vel[IR3::v], vel[IR3::w]};
}

} // end namespace gyronimo