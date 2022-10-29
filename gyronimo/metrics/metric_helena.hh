// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021-2022 Paulo Rodrigues and Manuel Assunção.

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

// @metric_helena.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_METRIC_HELENA
#define GYRONIMO_METRIC_HELENA

#include <gyronimo/metrics/metric_connected.hh>
#include <gyronimo/metrics/morphism_helena.hh>
#include <gyronimo/interpolators/interpolator2d.hh>

namespace gyronimo{

//! Covariant metric in `HELENA` curvilinear coordinates.
/*!
    Builds the metric from the information provided by a `parser_helena` object.
    The actual type of 2d interpolators to use is set by the specific
    `interpolator2d_factory` object pointer provided to the constructor.
    
    The right-handed coordinates are the square root of the poloidal flux per
    radian normalised to its value at the boundary (`u`, or `HELENA`
    @f$s=\sqrt{\Psi/\Psi_b}@f$), an angle on the poloidal cross section (`v`, or
    `HELENA` @f$\chi : B^\phi=q B^\chi@f$, in rads), and the toroidal angle
    (`w`, or `HELENA` @f$\phi@f$, in rads) measured **clockwise** when looking
    from the torus` top. The covariant metric units are m^2, as required by the
    parent class.

    @todo override to_contravariant to use the contravariant metric (it is
    available, right?) instead of the more costly default inversion.
*/
class metric_helena : public metric_connected {
 public:
  metric_helena(const morphism_helena *morph);
  virtual ~metric_helena() override;

  virtual SM3 operator()(const IR3& position) const override;
  virtual dSM3 del(const IR3& position) const override;

  virtual ddIR3 christoffel_first_kind (const IR3& r) const override;
  virtual ddIR3 christoffel_second_kind(const IR3& r) const override;

  const parser_helena* parser() const {return parser_;};

  virtual const morphism_helena* morph() const override {
    return static_cast<const morphism_helena*>(this->metric_connected::morph());
  };

 private:
  const parser_helena* parser_;
  const interpolator2d_factory *ifactory_;
  interpolator2d *guu_, *guv_, *gvv_, *gww_;
  double R0_, squaredR0_;
};

} // end namespace gyronimo

#endif // GYRONIMO_METRIC_HELENA
