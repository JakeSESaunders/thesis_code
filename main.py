# Your code goes here.
from functools import cache
from rings.dual_steenrod import DualSteenrod
from rings.bo import HomologyMO
from combinatorics.indices import is_power_of_two
from rings.free_over_dl_ops import PowersOfTwo, SModMod2
from homomorphisms.free_to_MO import free_to_MO
from power_series.hessenberg import hessenberg, hessenberg_det
from power_series.power_series import PowerSeries
from combinatorics.partitions import get_unordered_partitions
from power_series.powers_of_2_spherical import X, u1, u, Qx_2_pow_i

@cache
def xi_t_coeffs(d):
  if d == 1:
    return DualSteenrod.one()
  r = is_power_of_two(d)
  if r is not None:
    return DualSteenrod.xi(r)
  return DualSteenrod.zero()

@cache
def xi_t_coeffs_power(k, d):
  if k == 0:
    return DualSteenrod.one()
  if k == 1:
    return xi_t_coeffs(d)

  result = DualSteenrod.zero()
  for i in range(0, d + 1):
    lower_power = xi_t_coeffs_power(k - 1, i)
    next_power = xi_t_coeffs(d - i)
    result += lower_power*next_power

  return result

# TODO this is not quite right, there should be no x_2 appearing in u_1(t).
print("u_1(t) =", u1)
print()
for i in range(1, 4):
  print(f"u_{2**i}(t) =", u(i))
  print()

print(X)
print()
print(Qx_2_pow_i(0))
print(Qx_2_pow_i(0).circ(X))