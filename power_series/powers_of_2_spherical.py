from power_series.power_series import PowerSeries
from rings.free_over_dl_ops import PowersOfTwo
from combinatorics.indices import is_power_of_two
from power_series.dyer_lashof import dl_power_series
from dyer_lashof.free_over_dl_ops import dl_ops_on_generators_powers_of_2

# TODO a lot of this is copied from s_modmod_2_spherical.py. Would be good to have a common abstraction.

def Qx_2_pow_i(i):  
  x2_i = PowersOfTwo.x(i, specified_power=True)
  return dl_power_series(dl_ops_on_generators_powers_of_2, x2_i) # Returns Qx_{2^i}(t)

def X_coeffs(d):
  i = is_power_of_two(d)
  if i is None:
    return PowersOfTwo.zero()
  if i == 0:
    return PowersOfTwo.one()
  if i == 1:
    return PowersOfTwo.generator(
      (0, ())
    )
  operation = tuple([1 for i in range(i - 1)])
  return PowersOfTwo.generator(
    (0, operation)
  )

X = PowerSeries(X_coeffs)

# The power series given by 1/X(t) + 1/t. You can work out with pen and paper that this bizarre formula gives the right thing without considering -ve coefficients.
X_t_minus_1_over_t = (X.shift(-1).multiplicative_inverse() - PowerSeries.constant(PowersOfTwo.one())).shift(-1)

X1 = PowerSeries.constant(PowersOfTwo.generator(
  (0, ())
))

u1 = Qx_2_pow_i(0).circ(X) + X_t_minus_1_over_t + X1

def u(i):
  if i == 0:
    return u1
  return Qx_2_pow_i(i).circ(X)

def spherical(x):
  # reduce to single summand
  if len(x.summands) != 1:
    return sum([spherical(y) for y in x.get_summands_as_polys()], PowersOfTwo.zero())
  
  # reduce to single factor
  summand = x.summands[0]
  if len(summand) != 1:
    result = PowersOfTwo.one()
    for factor in summand:
      result *= spherical(PowersOfTwo.generator(factor))
    return result
  
  # break into cases based on length of I.
  x, I = summand[0]
  if len(I) == 0:
    if x == 0:
      raise ValueError("No spherical class for x_1.")
    return PowersOfTwo.x(x)

  if len(I) == 1:
    if x == 0 and I[0] == 1:
      raise ValueError("No spherical class for Q_1x_1.")
    return u(x).coeff(2**x + I[0])

  spherical_remainder = spherical(PowersOfTwo.Q(I[1:], x))
  power_series = dl_power_series(
    dl_ops_on_generators_powers_of_2,
    spherical_remainder
  )
  return power_series.coeff(spherical_remainder.max_degree() + I[0])