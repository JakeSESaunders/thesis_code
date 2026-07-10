from functools import cache
from math import floor, log2
from combinatorics.indices import is_power_of_two
from dyer_lashof.monomial import DyerLashofMonomial
from rings.free_over_dl_ops import FreePowersOfTwo, FormalDLOp
from power_series.powers_of_two_spherical import spherical
from homomorphisms.manifold_to_e_basis import e_to_manifold_basis
from homomorphisms.free_to_MO import free_to_MO

@cache
def get_classes_in_degree(n) -> list[FormalDLOp]:
  result = []
  k = is_power_of_two(n)
  if k is not None:
    result.append(
      FormalDLOp(DyerLashofMonomial.empty(upper=False), n)
    )
  for i in range(1, floor(n/2) + 1):
    difference = n - 2*i
    if difference == 0:
      continue
    for dl_op in get_classes_in_degree(i):
      if dl_op.QI.length() != 0 and dl_op.QI.operation[0] < difference:
        continue
      new_gen = dl_op.Q(difference)
      result.append(new_gen)

  return result

# Returns the classes of the form Q_{2^i_1} \cdots Q_{2^i_n}e_{2^{i_{n + 1}}} with i_1, \dots, i_{n + 1} a non-decreasing seq of positive integers
@cache
def get_good_classes_in_degree(n):
  result = []
  k = is_power_of_two(n)
  if k is not None:
    result.append(
      FormalDLOp(DyerLashofMonomial.empty(upper=False), n)
    )
  for i in range(1, floor(n/2) + 1):
    difference = n - 2*i
    if difference == 0:
      continue
    if is_power_of_two(difference) is None:
      continue
    for dl_op in get_good_classes_in_degree(i):
      if dl_op.QI.length() != 0 and dl_op.QI.operation[0] < difference:
        continue
      new_gen = dl_op.Q(difference)
      if new_gen.generator < difference:
        continue
      result.append(new_gen)

  return result

@cache
def get_classes_to_kill_in_degree(n):
  result = []
  # get the classes of the form Q_ne_{2^i} where n > 2^i.
  for i in range(0, floor(log2(n))):
    difference = n - 2*2**i
    if is_power_of_two(difference) is not None and difference > 2**i:
      result.append(
        FormalDLOp(
          DyerLashofMonomial.single(difference), 2**i
        )
      )
  # now get the classes which are prefaced by a non-power of 2.
  for i in range(1, floor(n/2) + 1):
    difference = n - 2*i
    if difference == 0:
      continue
    for dl_op in get_good_classes_in_degree(i):
      if dl_op.QI.length() != 0 and dl_op.QI.operation[0] < difference:
        continue
      if is_power_of_two(difference) is not None:
        continue
      new_gen = dl_op.Q(difference)
      result.append(new_gen)
  return result

@cache
def get_classes_already_killed_in_degree(n):
  result = []
  for i in range(1, floor(n/2) + 1):
    difference = n - 2*i
    if difference == 0:
      continue
    for dl_op in get_classes_to_kill_in_degree(i):
      if dl_op.QI.length() != 0 and dl_op.QI.operation[0] < difference:
        continue
      new_gen = dl_op.Q(difference)
      result.append(new_gen)

  return result

def print_spherical_classes(_degree):
  for i in range(1, _degree):
    print(i)
    print("All:       ", [str(x) for x in get_classes_in_degree(i)])
    print("Good:      ", [str(x) for x in get_good_classes_in_degree(i)])
    print("Spherical: ", [str(spherical(dl_op.QI, dl_op.x)) for dl_op in get_classes_in_degree(i) if not (dl_op.x == 1 and (dl_op.QI.length() == 0 or dl_op.QI.operation[-1] == 1))])
    print("To Kill:   ", [str(x) for x in get_classes_to_kill_in_degree(i)])
    print("Killed:    ", [str(x) for x in get_classes_already_killed_in_degree(i)])
    print()

def print_spherical_table(_degree):
  print("""\\begin{longtable}{|c|c|C{3cm}|C{3cm}|C{3cm}|}""")
  print("""\\hline""")
  print("""Degree & Class & Spherical & $f^0(\\text{Spherical})$ & Manifold Basis\\\\ \\hline""")
  last_degree = 0
  for i in range(1, _degree):
    for dl_op in get_classes_in_degree(i):
      
      if (dl_op.x == 1 and (dl_op.QI.length() == 0 or dl_op.QI.operation[-1] == 1)):
        continue

      x = FreePowersOfTwo.generator(dl_op)
      degree = x.max_degree()
      spherical_x = spherical(dl_op.QI, dl_op.x)
      image_x = free_to_MO(spherical_x)
      manifold_basis_image_x = e_to_manifold_basis(image_x)
      if degree != last_degree:
        print(f"""${degree}$ & ${x}$ & ${spherical_x}$ & ${image_x}$ & ${manifold_basis_image_x}$ \\\\ \\hline""")
      
  print("""\\end{longtable}""")