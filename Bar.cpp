#include "Bar.hpp"

Bar::Bar(size_t t_i, size_t t_j) : i(t_i), j(t_j), l0(0.0), kn(0.0), fn(0.0) {}

/**
 *  Initializes the Bar object with the given parameters.
 *
 *  @param t_kn   stiffness
 *  @param posi   position of the start point
 *  @param posj   position of the end point
 */
void Bar::init(double t_kn, vec2r &posi, vec2r &posj) {
  vec2r branch = posj - posi;
  l0 = branch.normalize();
  kn = t_kn;
  fn = 0.0;
}
