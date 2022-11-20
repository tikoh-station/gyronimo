// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Manuel Assunção.

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

// @metric_connected.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/metric_connected.hh>

namespace gyronimo {

metric_connected::metric_connected(const morphism *morph) 
		: morph_(morph) {
	
	if(morph_ == nullptr)
		error(__func__, __FILE__, __LINE__, " morphism cannot be nullptr.", 1);
}

//! General-purpose implementation of the local metric matrix `SM3`.
/*!
	Implemented the rule
	@f$ g_{\alpha\beta} = \textbf{e}_\alpha \cdot \textbf{e}_\beta @f$
*/
SM3 metric_connected::operator()(const IR3& r) const {
	dIR3 e = morph_->del(r);
	IR3 e1 = {e[dIR3::uu], e[dIR3::vu], e[dIR3::wu]};
	IR3 e2 = {e[dIR3::uv], e[dIR3::vv], e[dIR3::wv]};
	IR3 e3 = {e[dIR3::uw], e[dIR3::vw], e[dIR3::ww]};
	SM3 g = {inner_product(e1, e1), inner_product(e1, e2), inner_product(e1, e3),
			 inner_product(e2, e2), inner_product(e2, e3), inner_product(e3, e3)};
	return g;
}

//! General implementation of `metric_connected` derivatives
/*!
	Implemented the rule
	@f$ \partial_\gamma \, g_{\alpha\beta} = \Gamma_{\alpha\beta\gamma} + \Gamma_{\beta\alpha\gamma} @f$
*/
dSM3 metric_connected::del(const IR3& r) const {
	ddIR3 CF = christoffel_first_kind(r);
	return {
		CF[ddIR3::uuu] + CF[ddIR3::uuu], // uuu
		CF[ddIR3::uuv] + CF[ddIR3::uuv], // uuv
		CF[ddIR3::uuw] + CF[ddIR3::uuw], // uuw
		CF[ddIR3::uuv] + CF[ddIR3::vuu], // uvu
		CF[ddIR3::uvv] + CF[ddIR3::vuv], // uvv
		CF[ddIR3::uvw] + CF[ddIR3::vuw], // uvw
    	CF[ddIR3::uuw] + CF[ddIR3::wuu], // uwu
		CF[ddIR3::uvw] + CF[ddIR3::wuv], // uwv
		CF[ddIR3::uww] + CF[ddIR3::wuw], // uww
		CF[ddIR3::vuv] + CF[ddIR3::vuv], // vvu
		CF[ddIR3::vvv] + CF[ddIR3::vvv], // vvv
		CF[ddIR3::vvw] + CF[ddIR3::vvw], // vvw
    	CF[ddIR3::vuw] + CF[ddIR3::wuv], // vwu
		CF[ddIR3::vvw] + CF[ddIR3::wvv], // vwv
		CF[ddIR3::vww] + CF[ddIR3::wvw], // vww
		CF[ddIR3::wuw] + CF[ddIR3::wuw], // wwu
		CF[ddIR3::wvw] + CF[ddIR3::wvw], // wwv
		CF[ddIR3::www] + CF[ddIR3::www]  // www
	};
}

//! General-purpose implementation of the Jacobian.
double metric_connected::jacobian(const IR3& r) const {
	return morph_->jacobian(r);
}

//! General-purpose implementation of the Jacobian gradient.
/*!
    Implements the rule
    @f$ \partial_\alpha J = J \left(
		\Gamma^1_{\alpha 1} + \Gamma^2_{\alpha 2} + \Gamma^3_{\alpha 3}
	\right) @f$
*/
IR3 metric_connected::del_jacobian(const IR3& r) const {
	double J = jacobian(r);
	dIR3 ee = morph_->del_inverse(r);
	ddIR3 de = morph_->ddel(r);
	IR3 con = {
		ee[dIR3::uu] * de[ddIR3::uuu] + ee[dIR3::uv] * de[ddIR3::vuu] + ee[dIR3::uw] * de[ddIR3::wuu] + 
		ee[dIR3::vu] * de[ddIR3::uuv] + ee[dIR3::vv] * de[ddIR3::vuv] + ee[dIR3::vw] * de[ddIR3::wuv] + 
		ee[dIR3::wu] * de[ddIR3::uuw] + ee[dIR3::wv] * de[ddIR3::vuw] + ee[dIR3::ww] * de[ddIR3::wuw],
		ee[dIR3::uu] * de[ddIR3::uuv] + ee[dIR3::uv] * de[ddIR3::vuv] + ee[dIR3::uw] * de[ddIR3::wuv] + 
		ee[dIR3::vu] * de[ddIR3::uvv] + ee[dIR3::vv] * de[ddIR3::vvv] + ee[dIR3::vw] * de[ddIR3::wvv] + 
		ee[dIR3::wu] * de[ddIR3::uvw] + ee[dIR3::wv] * de[ddIR3::vvw] + ee[dIR3::ww] * de[ddIR3::wvw],
		ee[dIR3::uu] * de[ddIR3::uuw] + ee[dIR3::uv] * de[ddIR3::vuw] + ee[dIR3::uw] * de[ddIR3::wuw] + 
		ee[dIR3::vu] * de[ddIR3::uvw] + ee[dIR3::vv] * de[ddIR3::vvw] + ee[dIR3::vw] * de[ddIR3::wvw] + 
		ee[dIR3::wu] * de[ddIR3::uww] + ee[dIR3::wv] * de[ddIR3::vww] + ee[dIR3::ww] * de[ddIR3::www]
	};
	return (J * con);
}

//! Implementation of Christoffel symbols of the first kind @f$ \Gamma_{\alpha\beta\gamma} @f$
/*!
    Implements the rule
    @f$ \Gamma_{\alpha\beta\gamma} = \textbf{e}_\alpha \cdot 
	\frac{\partial^2 \textbf{x}}{\partial q^\beta \, \partial q^\gamma} @f$
*/
ddIR3 metric_connected::christoffel_first_kind(const IR3& r) const {
	dIR3 e = morph_->del(r);
	ddIR3 ddR = morph_->ddel(r);
	return contraction<first, first>(e, ddR);
}

//! Implementation of Christoffel symbols of the second kind @f$ \Gamma^\alpha_{\beta\gamma} @f$
/*!
    Implements the rule
    @f$ \Gamma^\alpha_{\beta\gamma} = \textbf{e}^\alpha \cdot 
	\frac{\partial^2 \textbf{x}}{\partial q^\beta \, \partial q^\gamma} @f$
*/
ddIR3 metric_connected::christoffel_second_kind(const IR3& r) const {
	dIR3 ee = morph_->del_inverse(r);
	ddIR3 ddR = morph_->ddel(r);
	return contraction<second, first>(ee, ddR);
}

} // end namespace gyronimo.
