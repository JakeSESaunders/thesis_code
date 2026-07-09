from rings.polynomial import PolyRing
from combinatorics.admissible import is_tuple_of_nondecreasing_positive_integers
from combinatorics.indices import is_power_of_two
from dyer_lashof.formal import DL

# Currently only supports having at most one class in each degree.
class FreeOverDLAlg(PolyRing):
  def is_valid_generator(dl_op):
    if not isinstance(dl_op, DL):
      return False
    if not dl_op.lower:
      return False
    if not dl_op.is_reduced():
      return False
    return True # child classes can implement further restrictions.

  @classmethod
  def generator(cls, dl_op: DL):
    leading_zeros, shorter_op = dl_op.strip_leading_zeros() # Deal with any leading zeros.
    return cls([(shorter_op, )])**(2**leading_zeros)

  def x(n: int):
    return FreeOverDLAlg.generator(
      DL.from_generator(n, upper=False)
    )

  # See above for parameters.
  # I is lower indexing.
  def Q(I, n: int):
    if isinstance(I, list):
      I = tuple(I)
    if not isinstance(I, tuple):
      I = [I]
      
    return FreeOverDLAlg.generator(
      DL(list(I), n, upper=False)
    )

  @classmethod
  def degree(cls, dl_op: DL):
    return dl_op.degree()

  def symbol(dl_op: DL):
    return str(dl_op)

# The ring H_*\S//2.
# Generators are classes dl_op of type DL such that dl_op.generator = 1.
class SModMod2(FreeOverDLAlg):
  def is_valid_generator(dl_op):
    return dl_op.generator == 1 and FreeOverDLAlg.is_valid_generator(dl_op)

  def symbol(dl_op: DL):
    result = 'a'
    operation = ''
    for i in reversed(dl_op.operation):
        operation = f"{operation}Q_{i}"
    return f"{operation}{result}"

  def a():
    return SModMod2.generator(
      DL.from_generator(1)
    )

  def Q(I: tuple):
    return SModMod2.generator(
      DL(list(I), 1)
    )

""" # The ring H_*Y^0, free over the DL ops on generators in degree 1, 2, 4, 8, ... (powers of 2).
# Generators are classes dl_op of type DL such that dl_op.generator is a power of 2.
class PowersOfTwo(FreeOverDLAlg):
  def is_valid_generator(dl_op):
    return (is_power_of_two(dl_op.generator) is not None) and FreeOverDLAlg.is_valid_generator(dl_op) """