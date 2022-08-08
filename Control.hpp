#ifndef CONTROL_HPP
#define CONTROL_HPP

#define VELOCITY_CONTROL 0
#define FORCE_CONTROL 1

/**
 * @brief Un control contient des flags pour dire si une force ou une vitesse est
 *        imposée selon les direction x et y. Il concerne uniquement les noeuds
 *
 */
class Control {
public:
  int xmode;
  double xvalue;
  int ymode;
  double yvalue;

  // Ctor
  Control(int xmode_, double xvalue_, int ymode_, double yvalue_);
};

#endif /* end of include guard: CONTROL_HPP */
