#include "Neighbor.hpp"

Neighbor::Neighbor(size_t t_ic, size_t t_jc, size_t t_in, size_t t_jn)
    : ic(t_ic), jc(t_jc), in(t_in), jn(t_jn), n(), contactState(0), fn(0.0), ft(0.0), glueState(0), fn_coh(0.0),
      ft_coh(0.0), kn_coh(0.0), kt_coh(0.0), fn_coh_max(0.0), ft_coh_max(0.0), yieldPower(2.0), Gc(0.0) {}