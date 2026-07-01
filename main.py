# Your code goes here.
from functools import cache
from rings.dual_steenrod import DualSteenrod
from combinatorics.indices import is_power_of_two

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

for k in range(5, 10):
  print(k, 
    xi_t_coeffs_power(k, 6)
  )