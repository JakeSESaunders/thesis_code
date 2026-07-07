from rings.polynomial import PolyRing
from combinatorics.admissible import is_tuple_of_nondecreasing_positive_integers
from combinatorics.indices import is_power_of_two

# The ring H_*\S//2.
# Generators are tuples of non-decreasing positive integers (lower index notation).
# e.g. (1, 2, 2, 5) represents Q_1Q_2^2Q_5a.
# Note the empty tuple (as a generator) represents the class a.
class SModMod2(PolyRing):
  def is_valid_generator(I):
    return is_tuple_of_nondecreasing_positive_integers(I)

  def a():
    return SModMod2([
        (
          ()
        )
      ])

  @classmethod
  def generator(cls, I):
    if I == 0:
      return cls.one()

    # deal with the case that some number of the leading terms is zero
    leading_zeros = 0
    while leading_zeros < len(I) and I[leading_zeros] == 0:
      leading_zeros += 1

    if leading_zeros == len(I):
      return SModMod2.a()**(2**leading_zeros)

    return cls([(
      I[leading_zeros:], 
    )])**(2**leading_zeros)

  @classmethod
  def degree(cls, I):
    if not cls.is_valid_generator(I):
      raise ValueError(f"Cannot determine degree of invalid generator {I}.")
    result = 1
    for i in reversed(range(len(I))):
      result = 2*result + I[i]
  
    return result

  def symbol(I):
    result = 'a'
    operation = ''
    for i in range(len(I)):
      operation = f"{operation}Q_{I[i]}"
    return f"{operation}{result}"

  def Q(I):
    return SModMod2.generator(I)

# The ring H_*Y_1, free over the DL ops on generators in degree 1, 2, 4, 8, ... (powers of 2).
# Generators are pairs (x, I), of a non-negative int x and a tuple I of non-decreasing positive integers, where e.g (3, (1, 5, 6)) denotes Q_1Q_5Q_6x_{2^3}.
class PowersOfTwo(PolyRing):
  def is_valid_generator(P):
    if not isinstance(P, tuple):
      return False
    if len(P) != 2:
      return False
    x, I = P

    if not isinstance(x, int):
      return False
    if x < 0:
      return False
    
    if not is_tuple_of_nondecreasing_positive_integers(I):
      return False

    return True

  @classmethod
  def generator(cls, P):
    n, I = P

    if len(I) == 0:
      return PowersOfTwo([
        (
          (n, ())
        )
      ])

    # deal with the case that some number of the leading terms is zero
    leading_zeros = 0
    while leading_zeros < len(I) and I[leading_zeros] == 0:
      leading_zeros += 1

    if leading_zeros == len(I):
      return PowersOfTwo.x(n)**(2**leading_zeros)

    operation = [(
      I[leading_zeros:], 
    )]
    P_new = (n, operation)
    return cls(P_new)**(2**leading_zeros)

  @classmethod
  def degree(cls, P):
    if not cls.is_valid_generator(P):
      raise ValueError(f"Cannot determine degree of invalid generator {P}.")
    x, I = P
    result = 2**x
    for i in reversed(range(len(I))):
      result = 2*result + I[i]
  
    return result

  def symbol(P):
    x, I = P
    result = f'x_{2**x}'
    operation = ''
    for i in range(len(I)):
      operation = f"{operation}Q_{I[i]}"
    return f"{operation}{result}"

  def x(i, specified_power=True):
    if specified_power:
      return PowersOfTwo.generator((i, ()))
    else:
      j = is_power_of_two(i)
      if j is None:
        raise ValueError(f"Specified generator {j} is not a power of 2.")
      return PowersOfTwo.x(j, specified_power=True)

  def Q(I, x, specified_power=True):
    if specified_power:
      return PowersOfTwo.generator((x, I))
    else:
      j = is_power_of_two(x)
      if j is None:
        raise ValueError(f"Specified generator {j} is not a power of 2.")
      return PowersOfTwo.Q(I, j, specified_power=True)
    