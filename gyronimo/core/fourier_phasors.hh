// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2023 Manuel Assunção.

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

// @fourier_phasors.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_FOURIER_PHASORS
#define GYRONIMO_FOURIER_PHASORS

#include <complex>
#include <valarray>

namespace gyronimo {

class fourier_phasors {
 public:
  using cis_container_t = std::valarray<std::complex<double>>;

  fourier_phasors(const std::valarray<double>& modes);
  fourier_phasors(const std::valarray<int>& modes);
  ~fourier_phasors() {};

  cis_container_t operator()(const double angle) const;

 private:
  std::valarray<int> mode_list_;
  std::valarray<size_t> mode_map_;

  std::valarray<int> lrint(const std::valarray<double>& darray) const;
  void synthesize_modes(const std::valarray<int>& modes);
  cis_container_t multi_pow(const std::complex<double>& base, 
      const std::valarray<int>& exps) const;
};

}

#endif