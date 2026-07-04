from rings.dual_steenrod import DualSteenrod
from functools import cache
from dyer_lashof.adem import adem_sequence

def xi(i):
  return DualSteenrod.xi(i)

# Computes DL ops via the series from Theorem 5.21 of Lawson - En ring spectra and Dyer--Lashof operations, up to degree k.
# Uses lower indexing, extends so that supplying -1 gives xi_1 for ease of computation.
@cache
def dyer_lashof_xi_1(k):
  if k == -1:
    return xi(1)
  if k == 0:
    return xi(1)**2

  summands = []
  i = 1
  while 2**i - 1 < 2 + k:
    summands.append(xi(i)*dyer_lashof_xi_1(k - (2**i - 1)))
  return sum(summands, DualSteenrod.zero())