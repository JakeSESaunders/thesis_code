from power_series.power_series import PowerSeries
from rings.free_over_dl_ops import FreeOverDLAlg
from dyer_lashof.formal import DL
from dyer_lashof.free_over_dl_ops import dyer_lashof_free
from power_series.dyer_lashof import dl_power_series
from combinatorics.indices import is_power_of_two

# Code for computing spherical classes congruent to a given Dyer--Lashof operation mod decomposables.
# Works for both H_*Y^0 and S//2.
# Assume that the generator in degree 1 is a or x_1, and has coaction given by \xi_1 \otimes 1 + 1 \otimes a or x_1.

def Qx_n(n):
  x_n = FreeOverDLAlg.x(n)
  return dl_power_series(dyer_lashof_free, x_n)

def X_coeffs(d):
  i = is_power_of_two(d)
  if i is None:
    return FreeOverDLAlg.zero()
  if i == 0:
    return FreeOverDLAlg.one()
  if i == 1:
    return FreeOverDLAlg.x(1)
  operation = tuple([1 for i in range(i - 1)])
  return FreeOverDLAlg.Q(operation, 1)

X = PowerSeries(X_coeffs)

# The power series given by 1/X(t) + 1/t. You can work out with pen and paper that this bizarre formula gives the right thing without considering -ve coefficients.
X_t_minus_1_over_t = (X.shift(-1).multiplicative_inverse() - PowerSeries.constant(FreeOverDLAlg.one())).shift(-1)

X1 = PowerSeries.constant(FreeOverDLAlg.x(1))

u1 = Qx_n(1).circ(X) + X_t_minus_1_over_t + X1

# These are power series with spherical coefficients.
def u(i):
  if i == 1:
    return u1
  return Qx_n(i).circ(X)

def spherical(x: FreeOverDLAlg):
  # reduce to single summand # TODO doesn't work like that!!! don't even need to do this.
  if len(x.summands) != 1:
    return sum([spherical(y) for y in x.get_summands_as_polys()], FreeOverDLAlg.zero())
  
  # reduce to single factor
  summand = x.summands[0]
  if len(summand) != 1:
    result = FreeOverDLAlg.one()
    for factor in summand:
      result *= spherical(FreeOverDLAlg.generator(factor))
    return result
  
  # break into cases based on length of operation.
  dl_op: DL = summand[0]
  I, n = dl_op.operation, dl_op.generator

  if len(I) == 0:
    if n == 1:
      raise ValueError("No spherical class for x_1.")
    return FreeOverDLAlg.x(n)

  if len(I) == 1:
    if n == 1 and I[0] == 1:
      raise ValueError("No spherical class for Q_1x_1.")
    return u(n).coeff(n + I[0])

  spherical_remainder = spherical(FreeOverDLAlg.Q(I[1:], n))
  power_series = dl_power_series(
    dyer_lashof_free,
    spherical_remainder
  )
  return power_series.coeff(spherical_remainder.max_degree() + I[0])