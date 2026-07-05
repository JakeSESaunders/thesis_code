from power_series.power_series import PowerSeries
from rings.free_over_dl_ops import PowersOfTwo
from combinatorics.indices import is_power_of_two

# TODO a lot of this is copied from s_modmod_2_spherical.py. Would be good to have a common abstraction.

def Qx2_pow_i_coeffs(i, d):
  if d < 2**i:
    return PowersOfTwo.zero()
  if d == 2**i:
    return PowersOfTwo.generator(
      (i, ())
    )**2
  
  P = (i, (d - 2**i, ))
  return PowersOfTwo.generator(P)

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

def Qx_2_pow_i(i):
  return PowerSeries(lambda d: Qx2_pow_i_coeffs(i, d))

# The power series given by 1/X(t) + 1/t. You can work out with pen and paper that this bizarre formula gives the right thing.
X_t_minus_1_over_t = (X.shift(-1).multiplicative_inverse() - PowerSeries.constant(PowersOfTwo.one())).shift(-1)

X1 = PowerSeries.constant(PowersOfTwo.generator(
  (0, ())
))

u1 = Qx_2_pow_i(0).circ(X) + X_t_minus_1_over_t + X1

def u(i):
  return Qx_2_pow_i(i).circ(X)


def spherical(k, i):
  if i == 0:
    if k < 1:
      raise ValueError("There is no spherical class for Q1a.")
    global u1
    return u1.coeff(2**i + k)

  return u(i).coeff(2**i + k)