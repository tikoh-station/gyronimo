// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022 Jorge Ferreira and Paulo Rodrigues.

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

// @metric_vmec.cc, this file is part of ::gyronimo::

#include <gyronimo/metrics/metric_vmec.hh>

namespace gyronimo{

metric_vmec::metric_vmec(const morphism_vmec* morph, 
        const interpolator1d_factory *ifactory) 
    : metric_connected(morph), parser_(morph->parser()),
      mnmax_(parser_->mnmax()), mnmax_nyq_(parser_->mnmax_nyq()),
      ns_(parser_->ns()), mpol_(parser_->mpol()), ntor_(parser_->ntor()), 
      signsgs_(parser_->signgs()), nfp_(parser_->nfp()),
      xm_(parser_->xm()), xn_(parser_->xn()), xm_nyq_(parser_->xm_nyq()), xn_nyq_(parser_->xn_nyq()),
      Rmnc_(nullptr), Zmns_(nullptr), gmnc_(nullptr)
      {
    // set radial grid block
    dblock_adapter s_range(parser_->radius());
    dblock_adapter s_half_range(parser_->radius_half());
    // set spectral components 
    Rmnc_ = new interpolator1d* [xm_.size()];
    Zmns_ = new interpolator1d* [xm_.size()];
    gmnc_ = new interpolator1d* [xm_.size()];
//@todo NEED TO FIX AXIS AND EDGE! TBI! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #pragma omp parallel for
    for(size_t i=0; i<xm_.size(); i++) {
      std::slice s_cut (i, s_range.size(), xm_.size());
      std::valarray<double> rmnc_i = (parser_->rmnc())[s_cut];
      Rmnc_[i] = ifactory->interpolate_data( s_range, dblock_adapter(rmnc_i));
      std::valarray<double> zmnc_i = (parser_->zmns())[s_cut];
      Zmns_[i] = ifactory->interpolate_data( s_range, dblock_adapter(zmnc_i));
      // note that gmnc is defined at half mesh
      std::slice s_h_cut (i+xm_nyq_.size(), s_half_range.size(), xm_nyq_.size());
      std::valarray<double> gmnc_i = (parser_->gmnc())[s_h_cut];
      gmnc_[i] = ifactory->interpolate_data( s_half_range, dblock_adapter(gmnc_i));
    };
}
metric_vmec::~metric_vmec() {
  if(Rmnc_) delete Rmnc_;
  if(Zmns_) delete Zmns_;
  if(gmnc_) delete gmnc_;
}
SM3 metric_vmec::operator()(const IR3& position) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];
  double R = 0.0, dR_ds = 0.0, dR_dtheta = 0.0, dR_dzeta = 0.0;
  double Z = 0.0, dZ_ds = 0.0, dZ_dtheta = 0.0, dZ_dzeta = 0.0;

   #pragma omp parallel for reduction(+: R, Z, dR_ds, dR_dtheta, dR_dzeta, dZ_ds, dZ_dtheta, dZ_dzeta)
  for (size_t i = 0; i<xm_.size(); i++) {  
    double m = xm_[i]; double n = xn_[i];
    double cosmn = std::cos( m*theta - n*zeta );
    double sinmn = std::sin( m*theta - n*zeta );
    double rmnc_i = (*Rmnc_[i])(s); 
    double zmns_i = (*Zmns_[i])(s);
    // assuming for now that vmec equilibrium has stellarator symmetry.
    R += rmnc_i * cosmn; 
    Z += zmns_i * sinmn;
    dR_ds += (*Rmnc_[i]).derivative(s) * cosmn; 
    dR_dtheta -= m * rmnc_i * sinmn; 
    dR_dzeta += n * rmnc_i * sinmn;
    dZ_ds += (*Zmns_[i]).derivative(s) * sinmn; 
    dZ_dtheta += m * zmns_i * cosmn; 
    dZ_dzeta -= n * zmns_i * cosmn; 
  };
  return {
    dR_ds * dR_ds + dZ_ds * dZ_ds,                      // g_uu
    dR_ds * dR_dzeta + dZ_ds * dZ_dzeta,                // g_uw
    dR_ds * dR_dtheta + dZ_ds * dZ_dtheta,              // g_uv
    R * R + dR_dzeta * dR_dzeta + dZ_dzeta * dZ_dzeta,  // g_vv
    dR_dtheta * dR_dzeta + dZ_dtheta * dZ_dzeta,        // g_vw
    dR_dtheta * dR_dtheta + dZ_dtheta * dZ_dtheta       // g_ww
  };
}
dSM3 metric_vmec::del(const IR3& position) const {
  double s = position[IR3::u];
  double zeta = position[IR3::v];
  double theta = position[IR3::w];
  double R = 0.0, Z = 0.0;
  double dR_ds = 0.0,        dR_dtheta = 0.0,       dR_dzeta = 0.0;
  double d2R_ds2 = 0.0,      d2R_dsdtheta = 0.0,    d2R_dsdzeta = 0.0;
  double d2R_dthetads = 0.0, d2R_dtheta2 = 0.0,     d2R_dthetadzeta = 0.0; 
  double d2R_dzetads = 0.0,  d2R_dzetadtheta = 0.0, d2R_dzeta2 = 0.0;
  double dZ_ds = 0.0,        dZ_dtheta = 0.0,       dZ_dzeta = 0.0;
  double d2Z_ds2 = 0.0,      d2Z_dsdtheta = 0.0,    d2Z_dsdzeta = 0.0;
  double d2Z_dthetads = 0.0, d2Z_dtheta2 = 0.0,     d2Z_dthetadzeta = 0.0; 
  double d2Z_dzetads = 0.0,  d2Z_dzetadtheta = 0.0, d2Z_dzeta2 = 0.0;

  #pragma omp parallel for reduction(+: R, dR_ds, dR_dtheta, dR_dzeta, d2R_ds2, d2R_dsdtheta, d2R_dsdzeta, d2R_dtheta2, d2R_dthetadzeta, d2R_dzeta2, Z, dZ_ds ,dZ_dtheta, dZ_dzeta, d2Z_ds2, d2Z_dsdtheta, d2Z_dsdzeta, d2Z_dtheta2, d2Z_dthetadzeta, d2Z_dzeta2)
  for (size_t i = 0; i<xm_.size(); i++) {  
    double m = xm_[i]; double n = xn_[i];
    double cosmn = std::cos( m*theta - n*zeta );
    double sinmn = std::sin( m*theta - n*zeta );
    double rmnc_i = (*Rmnc_[i])(s); 
    double zmns_i = (*Zmns_[i])(s);
    double d_rmnc_i = (*Rmnc_[i]).derivative(s); 
    double d_zmns_i = (*Zmns_[i]).derivative(s); 
    double d2_rmnc_i = (*Rmnc_[i]).derivative2(s);
    double d2_zmns_i = (*Zmns_[i]).derivative2(s);
    // assuming for now that vmec equilibrium has stellarator symmetry.
    R += rmnc_i * cosmn; Z += zmns_i * sinmn;
    dR_ds += d_rmnc_i * cosmn; 
    dR_dtheta -= m * rmnc_i * sinmn; 
    dR_dzeta += n * rmnc_i * sinmn;
    d2R_ds2 += d2_rmnc_i * cosmn; 
    d2R_dsdtheta -= m * d_rmnc_i * sinmn;
    d2R_dsdzeta += n * d_rmnc_i * sinmn;
    d2R_dtheta2 -= m * m * rmnc_i * cosmn;
    d2R_dthetadzeta += m * n * rmnc_i * cosmn;
    d2R_dzeta2 -= n * n * rmnc_i * cosmn;
    dZ_ds += d_zmns_i * sinmn; 
    dZ_dtheta += m * zmns_i * cosmn; 
    dZ_dzeta -= n * zmns_i * cosmn; 
    d2Z_ds2 += d2_zmns_i * sinmn;
    d2Z_dsdtheta += m * d_zmns_i * cosmn;
    d2Z_dsdzeta -= n * d_zmns_i * cosmn;
    d2Z_dtheta2 -= m * m * zmns_i * sinmn;
    d2Z_dthetadzeta += m * n * zmns_i * sinmn;
    d2Z_dzeta2 -= n * n * zmns_i * sinmn;
}
//@todo still need to test this carefully. Find a way to test d_g!
  return {
      2 * (dR_ds * d2R_ds2      + dZ_ds * d2Z_ds2), 
      2 * (dR_ds * d2R_dsdzeta  + dZ_ds * d2Z_dsdzeta), // d_i g_uu
      2 * (dR_ds * d2R_dsdtheta + dZ_ds * d2Z_dsdtheta), 
      dR_ds * d2R_dsdzeta       + dR_dzeta * d2R_ds2      + dZ_ds * d2Z_dsdzeta      + dZ_dzeta * d2Z_ds2,
      dR_ds * d2R_dzeta2        + dR_dzeta * d2R_dsdzeta  + dZ_ds * d2Z_dzeta2       + dZ_dzeta * d2Z_dsdzeta,// d_i g_uv
      dR_ds * d2R_dthetadzeta   + dR_dzeta * d2R_dsdzeta  + dZ_ds * d2Z_dthetadzeta  + dZ_dzeta * d2Z_dsdtheta, 
      dR_ds * d2R_dsdtheta      + dR_dtheta * d2R_ds2      + dZ_ds * d2Z_dsdtheta     + dZ_dtheta * d2Z_ds2,
      dR_ds * d2R_dthetadzeta   + dR_dtheta * d2R_dsdzeta  + dZ_ds * d2Z_dthetadzeta  + dZ_dtheta * d2Z_dsdzeta, // d_i g_uw
      dR_ds * d2R_dtheta2       + dR_dtheta * d2R_dsdtheta + dZ_ds * d2Z_dtheta2      + dZ_dtheta * d2Z_dsdtheta, 
      2 * (R * dR_ds     + dR_dzeta * d2R_dsdzeta     + dZ_dzeta * d2Z_dsdzeta), 
      2 * (R * dR_dzeta  + dR_dzeta * d2R_dzeta2      + dZ_dzeta * d2Z_dzeta2),  // d_i g_vv
      2 * (R * dR_dtheta + dR_dzeta * d2R_dthetadzeta + dZ_dzeta * d2Z_dthetadzeta),  
      dR_dtheta * d2R_dsdzeta     + dR_dzeta * d2R_dsdtheta     + dZ_dtheta * d2Z_dsdzeta      + dZ_dzeta * d2Z_dsdtheta,
      dR_dtheta * d2R_dzeta2      + dR_dzeta * d2R_dthetadzeta  + dZ_dtheta * d2Z_dzeta2       + dZ_dzeta * d2Z_dthetadzeta, // d_i g_vw
      dR_dtheta * d2R_dthetadzeta + dR_dzeta * d2R_dtheta2      + dZ_dtheta * d2Z_dthetadzeta  + dZ_dzeta * d2Z_dtheta2,  
      2 * (dR_dtheta * d2R_dsdtheta     + dZ_dtheta * d2Z_dsdtheta), 
      2 * (dR_dtheta * d2R_dthetadzeta  + dZ_dtheta * d2Z_dthetadzeta), // d_i g_ww
      2 * (dR_dtheta * d2R_dtheta2      + dZ_dtheta * d2Z_dtheta2),  
  };
}

} // end namespace gyronimo