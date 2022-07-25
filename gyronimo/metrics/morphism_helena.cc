// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Paulo Rodrigues.

// ::gyronimo:: is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ::gyronimo:: is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with ::gyronimo::.  If not, see <https://www.gnu.org/licenses/>.

// @morphism_helena.cc, this file is part of ::gyronimo::

#include <gsl/gsl_multiroots.h>
#include <gyronimo/core/error.hh>
#include <gyronimo/metrics/morphism_helena.hh>

namespace gyronimo{

morphism_helena::morphism_helena(
    const parser_helena *p, const interpolator2d_factory *ifactory)
    : parser_(p), R_(nullptr), z_(nullptr), 
	R0_(p->rmag()), squaredR0_(p->rmag()*p->rmag()) {
  double Rgeo = p->rgeo();
  double a = p->eps()*Rgeo;
  dblock_adapter s_range(p->s()), chi_range(p->chi());
  R_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(
          parser_helena::narray_type(a*p->x() + Rgeo)));
  z_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(
          parser_helena::narray_type(a*p->y() + 0.0)));

  guu_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(p->covariant_g11()));
  guv_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(p->covariant_g12()));
  gvv_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(p->covariant_g22()));
  gww_ = ifactory->interpolate_data(
      s_range, chi_range, dblock_adapter(p->covariant_g33()));
}
morphism_helena::~morphism_helena() {
  if(R_) delete R_;
  if(z_) delete z_;
  if(guu_) delete guu_;
  if(guv_) delete guv_;
  if(gvv_) delete gvv_;
  if(gww_) delete gww_;
}
IR3 morphism_helena::operator()(const IR3& q) const {
  double s = q[IR3::u], phi = q[IR3::w];
  double chi = this->reduce_chi(q[IR3::v]);
  double R = (*R_)(s, chi);
  return {R*std::cos(phi), -R*std::sin(phi), (*z_)(s, chi)};
}
int morphism_helena::root_f(const gsl_vector *q, void *target, gsl_vector *f) {
//   double u = reduce(gsl_vector_get(q, 0), 1);
  reduce2(q->data[0], q->data[1]);
  double u = q->data[0];
  double v = reduce(gsl_vector_get(q, 1), 2*std::numbers::pi);
  double target_R = ((struct root_target*) target)->R;
  double target_z = ((struct root_target*) target)->z;
  interpolator2d* interpolator_R = ((struct root_target*) target)->i_R;
  interpolator2d* interpolator_z = ((struct root_target*) target)->i_z;
  gsl_vector_set(f, 0, target_R - (*interpolator_R)(u, v));
  gsl_vector_set(f, 1, target_z - (*interpolator_z)(u, v));
  return GSL_SUCCESS;
}
IR3 morphism_helena::inverse(const IR3& X) const {
  const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_dnewton;
  gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(T, 2);

  double x = X[IR3::u], y = X[IR3::v], z = X[IR3::w];
  double R = std::sqrt(x*x + y*y);
  struct root_target target = {R, z, R_, z_};
  gsl_multiroot_function f = {&root_f, 2, &target};

  gsl_vector *guess_q = gsl_vector_alloc(2);
  gsl_vector_set(guess_q, 0, 0.5);
  gsl_vector_set(guess_q, 1, 0.0);

  const double tolerance = 1.0e-08;
  gsl_multiroot_fsolver_set(s, &f, guess_q);
  for (auto i : std::views::iota(1, 1000)) {
    // reduce2(s->x->data[0], s->x->data[1]);
	// double u = s->x->data[0];
    // double v = reduce(gsl_vector_get(s->x, 1), 2*std::numbers::pi);
	// std::cout << "# its: " << u << " " << v << " " 
	//   << gsl_vector_get(s->x, 0) << " " << gsl_vector_get(s->x, 1) << std::endl;
    if (gsl_multiroot_fsolver_iterate(s)) break;  // breaks if stuck;
    if (gsl_multiroot_test_residual(s->f, tolerance) != GSL_CONTINUE) break;
  }
  if (std::max(gsl_vector_get(s->f, 0), gsl_vector_get(s->f, 1)) > tolerance) {
    //   reduce2(s->x->data[0], s->x->data[1]);
    //   double u = s->x->data[0];
    //   double v = reduce(gsl_vector_get(s->x, 1), 2*std::numbers::pi);
    //   std::cout << "# error log: " << std::max(gsl_vector_get(s->f, 0), gsl_vector_get(s->f, 1))
	//     << " " << x << " " << y << " " << z << " " << R << " " << u << " " << v << std::endl;
	  error(__func__, __FILE__, __LINE__,
          "above tolerance after max iterations.", 1);
  }
  reduce2(s->x->data[0], s->x->data[1]);
  double u = s->x->data[0];
  double v = reduce(gsl_vector_get(s->x, 1), 2*std::numbers::pi);
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(guess_q);
  return {u, v, std::atan2(-y, x)};
}
dIR3 morphism_helena::del(const IR3& q) const {
  double s = q[IR3::u], phi = q[IR3::w];
  double chi = this->reduce_chi(q[IR3::v]);
  double R = (*R_)(s, chi),
      Ru = (*R_).partial_u(s, chi), Rv = (*R_).partial_v(s, chi);
  double cos = std::cos(phi), sin = std::sin(phi);
  return {Ru*cos, Rv*cos, -R*sin, -Ru*sin, -Rv*sin, -R*cos,
      (*z_).partial_u(s, chi), (*z_).partial_v(s, chi), 0.0};
}
double morphism_helena::jacobian(const IR3 &q) const {
	double s = q[IR3::u];
	double chi = this->reduce_chi(q[IR3::v]);
	double guu = squaredR0_ * (*guu_)(s, chi);
	double guv = squaredR0_ * (*guv_)(s, chi);
	double gvv = squaredR0_ * (*gvv_)(s, chi);
	double gww = squaredR0_ * (*gww_)(s, chi);
	return std::sqrt((guu*gvv - guv*guv)*gww);
}
dSM3 morphism_helena::g_del(const IR3& q) const {
  double s = q[IR3::u];
  double chi = this->reduce_chi(q[IR3::v]);
  return {
      squaredR0_*(*guu_).partial_u(s, chi),
      squaredR0_*(*guu_).partial_v(s, chi), 0.0, // d_i g_uu
      squaredR0_*(*guv_).partial_u(s, chi),
      squaredR0_*(*guv_).partial_v(s, chi), 0.0, // d_i g_uv
      0.0, 0.0, 0.0, // d_i g_uw
      squaredR0_*(*gvv_).partial_u(s, chi),
      squaredR0_*(*gvv_).partial_v(s, chi), 0.0, //d_i g_vv
      0.0, 0.0, 0.0, // d_i g_vw
      squaredR0_*(*gww_).partial_u(s, chi),
      squaredR0_*(*gww_).partial_v(s, chi), 0.0}; // d_i g_ww
}

//! Reduces an arbitrary angle chi to the interval [0:pi].
double morphism_helena::reduce_chi(double chi) const {
  double rchi = reduce(chi, 2*std::numbers::pi);
  if(parser_->is_symmetric() && rchi > std::numbers::pi)
      rchi = 2*std::numbers::pi - rchi;
  return rchi;
}

} // end namespace gyronimo