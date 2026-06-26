from rings.dual_steenrod import DualSteenrod
from functools import cache

# Computes the series from Theorem 5.21 of Lawson - En ring spectra and Dyer--Lashof operations, up to degree k.
def dl_relation_series(k):
  lhs = DualSteenrod.one() + DualSteenrod.zeta(1) + DualSteenrod.zeta(1)**2
  for i in range(1, k - 1):
    lhs += dyer_lashof_xi_1(i)

  rhs = DualSteenrod.one()
  i = 1
  while 2**i - 1 <= k + 2:
    rhs += DualSteenrod.zeta(2**i - 1)
    i += 1
  
  return lhs*rhs

@cache
def dyer_lashof_xi_1(k):
  return dl_relation_series(k + 1).get_homogeneous_part(2 + k)
