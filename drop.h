#include "bp.h"

std::tuple<int, real_t, real_t>
drop_sc(FactorGraph & fg, int t, int maxit_bp, real_t tol_bp, real_t damping_bp, int maxit_sc, real_t tol_sc, real_t damping_sc);
