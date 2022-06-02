#include <gyronimo/dynamics/electromagnetic_system.hh>

namespace gyronimo {

electromagnetic_system::electromagnetic_system(double Lref, double Vref, 
		double qom, const IR3field* E, const IR3field* B, const morphism *morph) : 
		Lref_(Lref), Vref_(Vref), Tref_(Lref / Vref),
		Oref_(  B ? qom * codata::e / codata::m_proton * B->m_factor() :
				E ? qom * codata::e / codata::m_proton * E->m_factor() : 0),
		electric_field_(E), magnetic_field_(B), 
		iEfield_time_factor_(E ? Tref_ / E->t_factor() : 0), 
		iBfield_time_factor_(B ? Tref_ / B->t_factor() : 0),
		field_morph_(morph) {

	// test if fields are consistent (normalization and metrics)
}

electromagnetic_system::state electromagnetic_system::operator()(const state &x, double t) const {

	IR3 pos = {x[0], x[1], x[2]};
	IR3 vel = {x[3], x[4], x[5]};

	IR3 E = {0, 0, 0};
	IR3 B = {0, 0, 0};

	if(field_morph_) {

		IR3 q = field_morph_->inverse(pos);
		dIR3 e = field_morph_->del(q);

		if(electric_field_) {
			E = electric_field_->contravariant(q, t * iEfield_time_factor_);
			E = contraction<second>(e, E);
		}
		if(magnetic_field_) {
			B = magnetic_field_->contravariant(q, t * iBfield_time_factor_);
			B = contraction<second>(e, B);
		}

	} else {

		if(electric_field_) E = electric_field_->contravariant(pos, t * iEfield_time_factor_);
		if(magnetic_field_) B = magnetic_field_->contravariant(pos, t * iBfield_time_factor_);
	}

	IR3 acl = Oref_ * (E + cartesian_cross_product(vel, B));

	return {vel[IR3::u], vel[IR3::v], vel[IR3::w],
			acl[IR3::u], acl[IR3::v], acl[IR3::w]};
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