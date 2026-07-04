from power_series.power_series import PowerSeries
from rings.free_over_dl_ops import SModMod2
from combinatorics.indices import is_power_of_two

def Qx1_coeffs(d):
  if d == 0:
    return SModMod2.zero()
  if d == 1:
    return SModMod2.generator(())**2
  
  return SModMod2.generator((d - 1, ))

Qx1 = PowerSeries(Qx1_coeffs)

def X_coeffs(d):
  i = is_power_of_two(d)
  if i is None:
    return SModMod2.zero()
  if i == 0:
    return SModMod2.one()
  if i == 1:
    return SModMod2.generator(())
  operation = tuple([1 for i in range(i - 1)])
  return SModMod2.generator(operation)

X = PowerSeries(X_coeffs)

# The power series given by 1/X(t) + 1/t. You can work out with pen and paper that this bizarre formula gives the right thing.
X_t_minus_1_over_t = (X.shift(-1).multiplicative_inverse() - PowerSeries.constant(SModMod2.one())).shift(-1)

X1 = PowerSeries.constant(SModMod2.generator(()))

u = Qx1.circ(X) + X_t_minus_1_over_t + X1

# Returns a spherical class congruent to Q_ia mod decomposables.
# Requires i > 1.
def spherical(i):
  if i < 1:
    raise ValueError("There is no spherical class for Q1a.")
  global u
  return u.coeff(1 + i)