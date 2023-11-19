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

// @fourier_phasors.cc, this file is part of ::gyronimo::

#include <gyronimo/core/fourier_phasors.hh>
#include <algorithm>
#include <vector>
#include <numeric>

namespace gyronimo {

fourier_phasors::fourier_phasors(const std::valarray<double>& modes)
    : mode_list_(modes.size()), mode_map_(modes.size()) {
  synthesize_modes(lrint(modes));
}

fourier_phasors::fourier_phasors(const std::valarray<int>& modes)
    : mode_list_(modes.size()), mode_map_(modes.size()) {
  synthesize_modes(modes);
}

fourier_phasors::cis_container_t fourier_phasors::operator()(
    const double angle) const {
  std::complex<double> e_iangle = {std::cos(angle), std::sin(angle)};
  cis_container_t powers = multi_pow(e_iangle, mode_list_);
  return powers[mode_map_];
}

std::valarray<int> fourier_phasors::lrint(
    const std::valarray<double>& darray) const {
  std::valarray<int> iarray(darray.size());
  std::transform(std::begin(darray), std::end(darray), std::begin(iarray),
      [](double d) -> int { return std::lrint(d); });
  return iarray;
}

void fourier_phasors::synthesize_modes(const std::valarray<int>& modes) {

  // create arrays to sort and unsort the modes
  std::valarray<size_t> sorter(modes.size());
  std::iota(std::begin(sorter), std::end(sorter), 0);
  std::stable_sort(std::begin(sorter), std::end(sorter), 
      [&modes](size_t s1, size_t s2) -> bool { return modes[s1] < modes[s2]; });

  std::valarray<size_t> unsorter(modes.size());
  for(size_t i = 0; i < modes.size(); ++i) unsorter[sorter[i]] = i;

  // create map
  std::valarray<int> temp = modes[sorter];
  std::vector<int> sorted_modes(std::begin(temp), std::end(temp));
  std::valarray<size_t> sorted_map(modes.size());
  sorted_map[0] = 0;
  for(size_t i = 1, s = 0; i < modes.size(); ++i) {
    if(sorted_modes[i-1] != sorted_modes[i]) ++s;
    sorted_map[i] = s;
  }
  mode_map_ = sorted_map[unsorter];

  // create synthesized mode list
  auto last = std::unique(std::begin(sorted_modes), std::end(sorted_modes));
  sorted_modes.erase(last, std::end(sorted_modes));
  mode_list_ = std::valarray<int>(sorted_modes.data(), sorted_modes.size());

  return;
}

fourier_phasors::cis_container_t fourier_phasors::multi_pow(
    const std::complex<double>& base, const std::valarray<int>& exps) const {
  cis_container_t base_powers(exps.size());
  std::transform(std::begin(exps), std::end(exps), std::begin(base_powers), 
      [&base](int e) -> std::complex<double> { 
        return e < 0 ? std::pow(1.0/base, -e) : std::pow(base, e); 
      });
  return std::move(base_powers);
}

}