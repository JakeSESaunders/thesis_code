# Your code goes here.
from rings.free_over_dl_ops import PowersOfTwo
from dyer_lashof.free_over_dl_ops import dyer_lashof_powers_of_two

x1 = PowersOfTwo.x(0)
x2 = PowersOfTwo.x(1)

print(dyer_lashof_powers_of_two(
  (4, 2), x1, upper=False
))

print(dyer_lashof_powers_of_two(
  (4, 2), x2, upper=False
))