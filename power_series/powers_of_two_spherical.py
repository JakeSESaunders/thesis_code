from power_series.power_series import PowerSeries
from rings.free_over_dl_ops import FreePowersOfTwo
from dyer_lashof.free_over_dl_ops import dyer_lashof_free_powers_of_two
from power_series.dyer_lashof import dl_power_series
from combinatorics.indices import is_power_of_two
from dyer_lashof.monomial import DyerLashofMonomial

# Code for computing spherical classes congruent to a given Dyer--Lashof operation mod decomposables using work of Hoyer.
# Assume that the generator in degree 1 is a or x_1, and has coaction given by \xi_1 \otimes 1 + 1 \otimes a or x_1.

def Qx(n: int):
  """Returns the power series whose t^0, ..., t^nn coefficients are 0 and whose t^{n + k} coefficient is Q_kx_n."""
  if is_power_of_two(n) is None: raise ValueError(f"Cannot define power series for x_{n} where {n} is not a power of 2.")
  x_n = FreePowersOfTwo.x(n)
  return dl_power_series(dyer_lashof_free_powers_of_two, x_n)

def X_coeffs(d):
  """Returns the power series whose t^1 coefficient is 1, whose t^{2**i} coefficient is Q_1^{i - 1} and all other coefficients zero.""" 
  i = is_power_of_two(d)
  if i is None:
    return FreePowersOfTwo.zero()
  if i == 0:
    return FreePowersOfTwo.one() # t^1.
  if i == 1:
    return FreePowersOfTwo.x(1)
  
  operation_tuple = tuple([1 for _ in range(i - 1)])
  return FreePowersOfTwo.QIx(operation_tuple, 1)

X = PowerSeries(X_coeffs)

# The power series given by 1/X(t) + 1/t. You can work out with pen and paper that this bizarre formula gives the right thing without considering -ve coefficients.
X_t_minus_1_over_t = (X.shift(-1).multiplicative_inverse() - PowerSeries.constant(FreePowersOfTwo.one())).shift(-1)

X1 = PowerSeries.constant(FreePowersOfTwo.x(1))

u1 = Qx(1).circ(X) + X_t_minus_1_over_t + X1 # Hoyer 3.4.1 says this is spherical.

# Hoyer Lemma 3.4.2 says these are spherical.
# TODO don't need these anymore.
def u(i):
  if i == 1:
    return u1
  return Qx(i).circ(X)

def spherical_from_spherical(QI: DyerLashofMonomial, x: FreePowersOfTwo):
  """If x is a homogeneous spherical class, return the spherical class corresponding to QIx."""
  if not x.is_homogeneous(): raise ValueError(f"Cannot compute spherical class corresponding to non-homogeneous class {x}.")
  
  if QI.length() == 0:
    return x # x is already spherical.

  if QI.length() > 1:
    Qk, QJ = QI.pop()
    return spherical_from_spherical(Qk, spherical_from_spherical(QJ, x))

  # Now QI is a single operation, so we can take the (|x| + k)-coefficient of the relevant power series.
  k = QI.operation[0]
  d = x.max_degree()
  return dl_power_series(dyer_lashof_free_powers_of_two, x).circ(X).coeff(d + k) # this has primitive coeffs by Hoyer 3.4.2.

def spherical(QI: DyerLashofMonomial, n: int):
  """Returns the spherical class corresponding to the generator QIx_n of H_*Y^0."""
  if is_power_of_two(n) is None: raise ValueError(f"Cannot find spherical class of generator x_{n} as {n} is not a power of 2.")
  if QI.upper: raise NotImplementedError("Not yet implemented spherical classes for upper indexed Dyer--Lashof monomials.")
  if not QI.is_reduced(): raise NotImplementedError("Not yet implemented spherical classes for non-reduced Dyer--Lashof monomials.")

  if n != 1:
    return spherical_from_spherical(QI, FreePowersOfTwo.x(n)) # x(n) is already spherical, so we can apply the previous result.
  
  # Now we know we are dealing with the separate case of x_1.
  if QI.length() == 0: raise ValueError("There is no spherical lift of x_1.")
  if QI.operation[-1] == 1: raise ValueError("There is no spherical lift of Q_1^ix_1.") # has to be all 1's as QI is reduced.

  J, k = QI.operation[:-1], QI.operation[-1]
  spherical_Qkx1 = u1.coeff(1 + k)
  return spherical_from_spherical(
    DyerLashofMonomial(J, upper=False), spherical_Qkx1
  )