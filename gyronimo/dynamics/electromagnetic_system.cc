#include <gyronimo/fields/electromagnetic_system.hh>

namespace gyronimo {

electromagnetic_system::electromagnetic_system(double Lref, double Vref, 
		double qom, const IR3field* E, const IR3field* B) : 
		Lref_(Lref), Vref_(Vref), Tref_(Lref / Vref),
		Oref_(  B ? qom * codata::e / codata::m_proton * B->m_factor() :
				E ? qom * codata::e / codata::m_proton * E->m_factor() : 0),
		electric_field_(E), magnetic_field_(B), 
		iEfield_time_factor_(E ? ), {

	// test if fields are consistent (normalization)
}

state electromagnetic_system::operator()(const state &x, double t) const {

	IR3 pos = {x[0], x[1], x[2]};
	IR3 vel = {x[3], x[4], x[5]};

	IR3 E = electric_field_->contravariant(pos, t * iEfield_time_factor_);
	IR3 B = magnetic_field_->contravariant(pos, t * iBfield_time_factor_);

	IR3 acl = Oref_ * (E + cartesian_cross_product(vel, B));

	return {vel[IR3::u], vel[IR3::v], vel[IR3::w],
			acl[IR3::u], acl[IR3::v], acl[IR3::w]};
}

IR3 electromagnetic_system::get_position(const state &s) const {
	return {x[0], x[1], x[2]};
}

IR3 electromagnetic_system::get_velocity(const state &s) const {
	return {x[3], x[4], x[5]};
}

state electromagnetic_system::generate_state(const IR3 &pos, const IR3 &vel) const {
	return {pos[IR3::u], pos[IR3::v], pos[IR3::w],
			vel[IR3::u], vel[IR3::v], vel[IR3::w]};
}

} // end namespace gyronimo