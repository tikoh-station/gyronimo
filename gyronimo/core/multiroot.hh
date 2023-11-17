// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2022-2023 Paulo Rodrigues and Manuel Assunção.

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

// @multiroot.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_MULTIROOT
#define GYRONIMO_MULTIROOT

#include <ranges>
#include <algorithm>
#include <functional>
#include <gsl/gsl_multiroots.h>
#include <gyronimo/core/error.hh>
#include <gyronimo/core/generators.hh>

namespace gyronimo {

//! Interface to [GSL](https://www.gnu.org/software/gsl) multiroot solver.
/*!
    Intended usage:
    ```
    container_t guess = {1.0, ..., 2.2};
    container_t root = multiroot(1.0e-09, 75)(my_rooting_function, guess);
    ```
    where `container_t` is any type following `SizedContiguousRange` and
    `my_rooting_function` is a `std::function<container_t(container_t)>` object
    containing the function to find the root of.
*/
class multiroot {
 public:
  template<SizedContiguousRange UserArgs>
  using user_function_t = typename std::function<UserArgs(const UserArgs&)>;
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  using user_dfunction_t = typename std::function<UserDArgs(const UserArgs&)>;
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  using user_fdfunction_t = typename std::function<std::pair<UserArgs,UserDArgs>(const UserArgs&)>;
  multiroot(double tolerance, size_t max_iterations)
      : max_iterations_(max_iterations), tolerance_(tolerance) {};
  template<SizedContiguousRange UserArgs>
  UserArgs operator()(
      user_function_t<UserArgs>& f, const UserArgs& guess) const;
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  UserArgs operator()(
      user_function_t<UserArgs>& f, user_dfunction_t<UserArgs,UserDArgs>& df,
      user_fdfunction_t<UserArgs,UserDArgs>& fdf, const UserArgs& guess) const;
 private:
  size_t max_iterations_;
  double tolerance_;
  template<SizedContiguousRange UserArgs>
  static int translation_function(
      const gsl_vector *args_gsl, void *f_pointer, gsl_vector *eval_gsl);
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  static int translation_dfunction(
      const gsl_vector *args_gsl, void *f_pointer, gsl_matrix *deval_gsl);
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  static int translation_fdfunction(const gsl_vector *args_gsl, 
    void *f_pointer, gsl_vector *eval_gsl, gsl_matrix *deval_gsl);
  template<SizedContiguousRange UserArgs>
  auto allocate_gsl_objects(
      user_function_t<UserArgs>& f, const UserArgs& guess) const;
  template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
  auto allocate_gsl_objects(
      user_function_t<UserArgs>& f, user_dfunction_t<UserArgs,UserDArgs>& df,
      user_fdfunction_t<UserArgs,UserDArgs>& fdf, const UserArgs& guess) const;
  inline void deallocate_gsl_objects(
      gsl_multiroot_fsolver*, gsl_vector*, gsl_multiroot_function*, void**) const;
  inline void deallocate_gsl_objects(gsl_multiroot_fdfsolver*, 
      gsl_vector*, gsl_multiroot_function_fdf*, void**) const;
};

template<SizedContiguousRange UserArgs>
UserArgs multiroot::operator()(
    user_function_t<UserArgs>& f, const UserArgs& guess) const {
  auto [solver, guess_gsl, struct_f_gsl, par] = allocate_gsl_objects(f, guess);
  for(auto iteration : std::views::iota(1u, max_iterations_)) {
    int flag = gsl_multiroot_fsolver_iterate(solver);
    switch(flag) {
      case GSL_ENOPROG: error(__func__, __FILE__, __LINE__,
          "iteration is stuck.", 1);
      case GSL_ENOPROGJ: error(__func__, __FILE__, __LINE__,
          "jacobian not improving the solution.", 1);
      case GSL_EBADFUNC: error(__func__, __FILE__, __LINE__,
          "singular user-supplied function (Inf/NaN).", 1);
    }
    if(gsl_multiroot_test_residual(solver->f, tolerance_) == GSL_SUCCESS) break;
  }
  if(gsl_multiroot_test_residual(solver->f, tolerance_) == GSL_CONTINUE)
      error(__func__, __FILE__, __LINE__,
          "still above tolerance after max iterations.", 1);
  UserArgs root = generate_sized<UserArgs>(solver->x->size);
  std::copy(solver->x->data, solver->x->data + solver->x->size, root.begin());
  deallocate_gsl_objects(solver, guess_gsl, struct_f_gsl, par);
  return root;
}

template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
UserArgs multiroot::operator()(
    user_function_t<UserArgs>& f, user_dfunction_t<UserArgs,UserDArgs>& df,
      user_fdfunction_t<UserArgs,UserDArgs>& fdf, const UserArgs& guess) const {
  auto [solver, guess_gsl, struct_f_gsl, par] = allocate_gsl_objects(f, df, fdf, guess);
  for(auto iteration : std::views::iota(1u, max_iterations_)) {
    int flag = gsl_multiroot_fdfsolver_iterate(solver);
    switch(flag) {
      case GSL_ENOPROG: error(__func__, __FILE__, __LINE__,
          "iteration is stuck.", 1);
      case GSL_ENOPROGJ: error(__func__, __FILE__, __LINE__,
          "jacobian not improving the solution.", 1);
      case GSL_EBADFUNC: error(__func__, __FILE__, __LINE__,
          "singular user-supplied function (Inf/NaN).", 1);
    }
    if(gsl_multiroot_test_residual(solver->f, tolerance_) == GSL_SUCCESS) break;
  }
  if(gsl_multiroot_test_residual(solver->f, tolerance_) == GSL_CONTINUE)
      error(__func__, __FILE__, __LINE__,
          "still above tolerance after max iterations.", 1);
  UserArgs root = generate_sized<UserArgs>(solver->x->size);
  std::copy(solver->x->data, solver->x->data + solver->x->size, root.begin());
  deallocate_gsl_objects(solver, guess_gsl, struct_f_gsl, par);
  return root;
}

template<SizedContiguousRange UserArgs>
int multiroot::translation_function(
    const gsl_vector *args_gsl, void *f_pointer, gsl_vector *eval_gsl) {
  UserArgs args = generate_sized<UserArgs>(args_gsl->size);
  std::copy(args_gsl->data, args_gsl->data + args_gsl->size, args.begin());
  UserArgs eval = (*static_cast<user_function_t<UserArgs>*>(
    static_cast<void**>(f_pointer)[0]))(args);
  std::ranges::copy(eval, eval_gsl->data);
  return GSL_SUCCESS;
}

template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
int multiroot::translation_dfunction(
    const gsl_vector *args_gsl, void *f_pointer, gsl_matrix *deval_gsl) {
  UserArgs args = generate_sized<UserArgs>(args_gsl->size);
  std::copy(args_gsl->data, args_gsl->data + args_gsl->size, args.begin());
  UserDArgs deval = (*static_cast<user_dfunction_t<UserArgs,UserDArgs>*>(
    static_cast<void**>(f_pointer)[1]))(args);
  std::ranges::copy(deval, deval_gsl->data);
  return GSL_SUCCESS;
}

template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
int multiroot::translation_fdfunction(
    const gsl_vector *args_gsl, void *f_pointer, 
    gsl_vector *eval_gsl, gsl_matrix *deval_gsl) {
  UserArgs args = generate_sized<UserArgs>(args_gsl->size);
  std::copy(args_gsl->data, args_gsl->data + args_gsl->size, args.begin());
  std::pair<UserArgs,UserDArgs> eval = 
      (*static_cast<user_fdfunction_t<UserArgs,UserDArgs>*>(
        static_cast<void**>(f_pointer)[2]))(args);
  std::ranges::copy(eval.first, eval_gsl->data);
  std::ranges::copy(eval.second, deval_gsl->data);
  return GSL_SUCCESS;
}

template<SizedContiguousRange UserArgs>
auto multiroot::allocate_gsl_objects(
    user_function_t<UserArgs>& f, const UserArgs& guess) const {
  const size_t n = guess.size();
  const gsl_multiroot_fsolver_type *type = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver *solver = gsl_multiroot_fsolver_alloc(type, n);
  gsl_vector *guess_gsl = gsl_vector_alloc(n);
  void **par = new void*{&f};
  auto *struct_f_gsl =
      new gsl_multiroot_function {&translation_function<UserArgs>, n, par};
  if(!solver || !guess_gsl || !struct_f_gsl) error(
      __func__, __FILE__, __LINE__, "gsl object allocation failed.", 1);
  std::ranges::copy(guess, guess_gsl->data);
  gsl_multiroot_fsolver_set(solver, struct_f_gsl, guess_gsl);
  return std::tuple<
      gsl_multiroot_fsolver*, gsl_vector*, gsl_multiroot_function*, void**
          >({solver, guess_gsl, struct_f_gsl, par});
}

template<SizedContiguousRange UserArgs, SizedContiguousRange UserDArgs>
auto multiroot::allocate_gsl_objects(
    user_function_t<UserArgs>& f, user_dfunction_t<UserArgs,UserDArgs>& df,
    user_fdfunction_t<UserArgs,UserDArgs>& fdf, const UserArgs& guess) const {
  const size_t n = guess.size();
  const gsl_multiroot_fdfsolver_type *type = gsl_multiroot_fdfsolver_newton;
  gsl_multiroot_fdfsolver *solver = gsl_multiroot_fdfsolver_alloc(type, n);
  gsl_vector *guess_gsl = gsl_vector_alloc(n);
  void **par = new void*[3]{&f, &df, &fdf};
  auto *struct_fdf_gsl =
      new gsl_multiroot_function_fdf {&translation_function<UserArgs>,
          &translation_dfunction<UserArgs,UserDArgs>, 
          &translation_fdfunction<UserArgs,UserDArgs>, n, par};
  if(!solver || !guess_gsl || !struct_fdf_gsl) error(
      __func__, __FILE__, __LINE__, "gsl object allocation failed.", 1);
  std::ranges::copy(guess, guess_gsl->data);
  gsl_multiroot_fdfsolver_set(solver, struct_fdf_gsl, guess_gsl);
  return std::tuple<
      gsl_multiroot_fdfsolver*, gsl_vector*, gsl_multiroot_function_fdf*,
      void**>({solver, guess_gsl, struct_fdf_gsl, par});
}

inline void multiroot::deallocate_gsl_objects(
    gsl_multiroot_fsolver* solver,
    gsl_vector* guess_gsl, gsl_multiroot_function* struct_f_gsl, void** par) const {
  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(guess_gsl);
  delete struct_f_gsl;
  delete par;
}

inline void multiroot::deallocate_gsl_objects(
    gsl_multiroot_fdfsolver* solver, gsl_vector* guess_gsl, 
    gsl_multiroot_function_fdf* struct_fdf_gsl, void** par) const {
  gsl_multiroot_fdfsolver_free(solver);
  gsl_vector_free(guess_gsl);
  delete struct_fdf_gsl;
  delete par;
}

} // end namespace gyronimo.

#endif // GYRONIMO_MULTIROOT
