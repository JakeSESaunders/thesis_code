# Your code goes here.
from rings.free_over_dl_ops import PowersOfTwo
from dyer_lashof.free_over_dl_ops import dyer_lashof_powers_of_two
from power_series.powers_of_2_spherical import spherical
from homomorphisms.free_to_MO import free_to_MO
from homomorphisms.manifold_to_e_basis import e_to_manifold_basis

x1 = PowersOfTwo.x(0)
x2 = PowersOfTwo.x(1)
x4 = PowersOfTwo.x(2)

Q1x4 = spherical(
  1, 2
)

ftMO = free_to_MO(Q1x4)

print(
  ftMO
)

print(
  e_to_manifold_basis(ftMO)
)

print(Q1x4)

""" for i in range(2, 16):
  print(i)
  ftMO = free_to_MO(spherical(i, 1))
  print(
    ftMO
  )
  print(
    e_to_manifold_basis(ftMO)
  )
  print() """