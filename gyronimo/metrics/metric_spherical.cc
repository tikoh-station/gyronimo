// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Paulo Rodrigues and Manuel Assunção.

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

// @metric_spherical.cc, this file is part of ::gyronimo::

#include <cmath>
#include <gyronimo/metrics/metric_spherical.hh>

namespace gyronimo {

metric_spherical::metric_spherical(const morphism_spherical *morph)
    : metric_connected(morph),
	  Lref_(morph->Lref()),
      Lref_squared_(Lref_*Lref_),
      Lref_cube_(Lref_*Lref_*Lref_), 
	  iLref_squared_(1/Lref_squared_) {
}

//! General-purpose implementation of the covariant metric.
/*!
	Implements the covariant metric in spherical coordinates: 
	@f$ g_{\alpha\beta} = \left(\begin{matrix}
		L_{ref}^2 && 0 && 0 \\
		0 && L_{ref}^2 \, R^2 && 0 \\
		0 && 0 && L_{ref}^2 \, R^2 \, \sin^2 \theta
	\end{matrix}\right) @f$
*/
SM3 metric_spherical::operator()(const IR3& r) const {
  double factor = Lref_squared_*r[IR3::u]*r[IR3::u];
  double sinv = std::sin(r[IR3::v]);
  return {Lref_squared_, 0.0, 0.0, factor, 0.0, factor*sinv*sinv};
}

//! General-purpose implementation of the inverse (i.e., contravariant metric).
/*!
	Implements the contravariant metric in spherical coordinates: 
	@f$ g^{\alpha\beta} = \left(\begin{matrix}
		\frac{1}{L_{ref}^2} && 0 && 0 \\
		0 && \frac{1}{L_{ref}^2 \, R^2} && 0 \\
		0 && 0 && \frac{1}{L_{ref}^2 \, R^2 \, \sin^2 \theta}
	\end{matrix}\right) @f$
*/
SM3 metric_spherical::inverse(const IR3& r) const {
  double ifactor = iLref_squared_/(r[IR3::u]*r[IR3::u]);
  double isinv = 1/std::sin(r[IR3::v]);
  return {iLref_squared_, 0.0, 0.0, ifactor, 0.0, ifactor*isinv*isinv};
}
dSM3 metric_spherical::del(const IR3& r) const {
  double cosv = std::cos(r[IR3::v]), sinv = std::sin(r[IR3::v]);
  double factor = 2.0*Lref_squared_*r[IR3::u];
  return {
      0.0, 0.0, 0.0, // d_i g_uu (i=u,v,w)
      0.0, 0.0, 0.0, // d_i g_uv
      0.0, 0.0, 0.0, // d_i g_uw
      factor, 0.0, 0.0, //d_i g_vv
      0.0, 0.0, 0.0, // d_i g_vw
      factor*sinv*sinv, factor*r[IR3::u]*sinv*cosv, 0.0}; // d_i g_ww
}

//! General-purpose implementation of the Jacobian of the transformation in point @f$ q^\alpha @f$.
/*!
	Implements the Jacobian in spherical coordinates: 
	@f$ J = L_{ref}^3 \, R^2 \, \sin \theta @f$
*/
double metric_spherical::jacobian(const IR3& r) const {
  return Lref_cube_*r[IR3::u]*r[IR3::u]*std::sin(r[IR3::v]);
}

//! General-purpose implementation of the Jacobian gradient in point @f$ q^\alpha @f$.
/*!
	Implements the Jacobian gradient in spherical coordinates: 
	@f$ \nabla J = \left( L_{ref}\,r\,\cos\phi\,(1+\sin^2\theta), 
	L_{ref}\,r\,\sin\phi\,(1+\sin^2\theta), 
	L_{ref}\,r\,\sin\theta\,\cos\theta \right) @f$
*/
IR3 metric_spherical::del_jacobian(const IR3& r) const {
  double cosv = std::cos(r[IR3::v]), sinv = std::sin(r[IR3::v]);
  double cosw = std::cos(r[IR3::w]), sinw = std::sin(r[IR3::w]);
  double r0 = Lref_squared_*r[IR3::u];
  double factor = (1+sinv*sinv);
  return {r0*cosw*factor, r0*sinw*factor, r0*sinv*cosv};
}
IR3 metric_spherical::to_covariant(const IR3& B, const IR3& r) const {
  double factor = Lref_squared_*r[IR3::u]*r[IR3::u];
  double sinv = std::sin(r[IR3::v]);
  return {Lref_squared_*B[IR3::u],
      factor*B[IR3::v], factor*sinv*sinv*B[IR3::w]};
}
IR3 metric_spherical::to_contravariant(const IR3& B, const IR3& r) const {
  double factor = Lref_squared_*r[IR3::u]*r[IR3::u];
  double sinv = std::sin(r[IR3::v]);
  return {B[IR3::u]/Lref_squared_,
      B[IR3::v]/factor, B[IR3::w]/(factor*sinv*sinv)};
}
ddIR3 metric_spherical::christoffel_first_kind(const IR3& q) const {
	double r = Lref_squared_*q[IR3::u];
	double s = std::sin(q[IR3::v]);
	double c = std::cos(q[IR3::v]);
	double rss = r*s*s;
	double r2sc = r*q[IR3::u]*s*c;
	return {
		0, 0, 0, -r, 0, -rss,
		0, r, 0, 0, 0, -r2sc,
		0, 0, rss, 0, r2sc, 0
	};
}
ddIR3 metric_spherical::christoffel_second_kind(const IR3& q) const {
	double r = q[IR3::u];
	double ir = 1/q[IR3::u];
	double s = std::sin(q[IR3::v]);
	double c = std::cos(q[IR3::v]);
	return {
		0, 0, 0, -r, 0, -r*s*s,
		0, ir, 0, 0, 0, -s*c,
		0, 0, ir, 0, c/s, 0
	};
}

} // end namespace gyronimo.
