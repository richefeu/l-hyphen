#include "Neighbor.hpp"

Neighbor::Neighbor(size_t ic_, size_t jc_, size_t in_, size_t jn_)
    : ic(ic_), jc(jc_), in(in_), jn(jn_), n(), contactState(0), fn(0.0), ft(0.0), glueState(0), fn_coh(0.0),
      ft_coh(0.0), kn_coh(0.0), kt_coh(0.0), fn_coh_max(0.0), ft_coh_max(0.0), yieldPower(2.0) {}