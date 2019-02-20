#include "incomplete_qr_builder.h"

template<> double SparseIncompleteQRBuilder<double>::sqnorm(const double &x) {
  return x*x;
}

template<> double SparseIncompleteQRBuilder<double>::inner(const double &x, const double &y) {
  return x*y;
}

template<> bool SparseIncompleteQRBuilder<double>::cmp_larger_abs_coef(
  const SparseIncompleteQRBuilder<double>::i_c &a,
  const SparseIncompleteQRBuilder<double>::i_c &b) {
  return std::abs(a.second) > std::abs(b.second);
}

template<> double SparseIncompleteQRBuilder<typename std::complex<double> >::sqnorm(const std::complex<double> &x) {
  return std::norm(x);
}

template<> std::complex<double> SparseIncompleteQRBuilder<std::complex<double> >::inner(const std::complex<double> &x, const std::complex<double> &y) {
  return x*std::conj(y);
}

template<> bool SparseIncompleteQRBuilder<std::complex<double> >::cmp_larger_abs_coef(
  const SparseIncompleteQRBuilder<std::complex<double> >::i_c &a,
  const SparseIncompleteQRBuilder<std::complex<double> >::i_c &b) {
  return std::norm(a.second) > std::norm(b.second);
}
