#include "Bar.hpp"

Bar::Bar(size_t i_, size_t j_) : i(i_), j(j_), l0(0.0), kn(0.0), fn(0.0) {}

void Bar::init(double kn_, vec2r &posi, vec2r &posj) {
  vec2r branch = posj - posi;
  l0 = branch.normalize();
  kn = kn_;
  fn = 0.0;
}