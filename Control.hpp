#ifndef CONTROL_HPP
#define CONTROL_HPP

#define VELOCITY_CONTROL 0
#define FORCE_CONTROL 1

/**
 *  Un control contient des flags pour dire si une force ou une vitesse est
 *  imposée selon les direction x et y. Il concerne uniquement les noeuds
 *
 */
class Control {
public:
  int xmode;
  double xvalue;
  int ymode;
  double yvalue;

  // Ctor
  Control(int t_xmode, double t_xvalue, int t_ymode, double t_yvalue);
};

#endif /* end of include guard: CONTROL_HPP */
