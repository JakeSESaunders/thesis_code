from functools import cache
from collections import Counter
from strings.polynomial_string import polynomial_to_string

# Base class for a graded polynomial ring with addition mod 2.
# Objects of this class are elements of the polynomial ring.
# "summands" is a list of sorted tuples representing summands of a polynomial.
# The default implementation sorts the list of summands and each tuple of factors, but does not check for duplicates.
# We instead check this when performing arithmetic operations.
# Each element of the tuple should be a "generator", a hashable type with equality and a comparison (le/ge) dunder method, such that is_valid_generator returns True. 
# By convention, zero has summands [] and one has summands [()].
class PolyRing:
  """An element of a polynomial ring over F2, the field with 2 elements.
  
  Keyword arguments:
  summands -- A list of tuples (containing no duplicates), each tuple consisting of generators of the polynomial ring.
  """
  def __init__(self, summands: list[tuple]):
    R = self.__class__

    new_summands = [] # a list of sorted summands
    counter = Counter(summands)
    for summand in counter:
      # check the input for duplicates
      if counter[summand] != 1: raise ValueError(f"List of summands may not contain duplicates: summand {summand} appeared {counter[summand]} times.")
      new_summand = []

      # check for invalid generators
      for factor in summand:
        if not R.is_valid_generator(factor): raise ValueError(f"Summand {summand} contains invalid generator {factor}.")
        new_summand.append(factor)
      
      new_summands.append(tuple(sorted(new_summand)))
    self.summands = sorted(new_summands)

  @classmethod
  def is_valid_generator(cls, x) -> bool:
    """Returns a bool indicating if the parameter x can be interpreted as a generator of this polynomial ring."""
    return True

  @classmethod
  @cache
  def get_basis(cls, degree: int) -> list:
    """Returns an additive basis of this PolyRing in the given degree. The basis elements must consist of classes with a single summand."""
    raise NotImplementedError("Base class PolyRing does not implement get_basis.") 

  @classmethod
  def generator(cls, x):
    """Return the element of this polynomial ring with a single summand having single factor x."""
    return cls([(x, )])

  @classmethod
  def degree(cls, x) -> int:
    """Given a generator x, return the degree of x."""
    if not cls.is_valid_generator(x): raise ValueError(f"Cannot determine degree of invalid generator {x}.")
    raise NotImplementedError("Base class PolyRing does not implement degree.")

  @classmethod
  def degree_of_summand(cls, summand) -> int:
    """Return the degree of summand: the sum of the degrees of its factors.
    If the summand is empty, i.e. is 0, then this returns 0 rather than -\infty.
    
    Keyword arguments:
    summand: a tuple of generators of this polynomial ring.
    """
    degree_of_summand = 0
    if not isinstance(summand, tuple): raise ValueError(f"Cannot obtain degree of summand {summand} which is not a tuple of generators.")
    for factor in summand:
      degree_of_summand += cls.degree(factor)
    return degree_of_summand

  @classmethod
  def one(cls):
    """Return the multiplicative identity of this polynomial ring."""
    return cls([()])

  @classmethod
  def zero(cls):
    """Return the additive identity of this polynomial ring."""
    return cls([])

  # TODO this could probably be faster.
  def __eq__(self, other):
    if not isinstance(other, self.__class__):
      return False
    
    for summand in self.summands + other.summands:
      if (summand not in other.summands) or (summand not in self.summands):
        return False

    return True

  def __add__(self, other):
    R = self.__class__
    if not isinstance(other, R): raise ValueError("{self} and {other} are not elements of the same ring.")

    summands = self.summands.copy()
    for summand in other.summands:
      if summand in summands:
        summands.remove(summand)
      else:
        summands.append(summand)
    return R(summands)

  def __sub__(self, other):
    return self + other # we're mod 2 so + and - are the same.

  def __mul__(self, other):
    R = self.__class__
    if not isinstance(other, R): raise ValueError("{self} and {other} are not elements of the same ring.")
    # If there are no summands, the polynomial is zcopy meero.
    if self.summands == [] or other.summands == []:
      return R.zero()

    # If self has multiple summands, multiply each summand by other and add the result.
    if len(self.summands) > 1:
      polys = self.get_summands_as_polys()
      return sum([poly * other for poly in polys], R.zero())
    
    # Do the same for other.
    if len(other.summands) > 1:
      polys = other.get_summands_as_polys()
      return sum([self * poly for poly in polys], R.zero())
    
    # Now self and other both have only a single summand.
    summand_1 = self.summands[0]
    summand_2 = other.summands[0]
    
    product = tuple(sorted(summand_1 + summand_2))
    return R([product])

  def __pow__(self, other):
    if not isinstance(other, int): raise ValueError(f"Cannot raise polynomial to the power of non-integer value {other}.")
    if other < 0: raise ValueError(f"Cannot raise polynomial to the power of negative integer {other}.")
    result = self.__class__.one()
    for i in range(0, other):
      result *= self
    return result

  @classmethod
  def symbol(cls, x) -> str:
    """Return the string symbol to use for the generator x."""
    raise NotImplementedError("Function symbol not implemented for base PolyRing class.")

  def __str__(self):
    R = self.__class__
    return polynomial_to_string(R.symbol, self.summands)

  def __hash__(self):
    return f"{self.__class__.__name__}, {self.__str__()}".__hash__()

  def max_degree(self) -> int:
    """Returns the maximum degree of all summands of a polynomial."""
    R = self.__class__
    current_max = 0
    for summand in self.summands:
      degree_of_summand = R.degree_of_summand(summand)
      current_max = max(current_max, degree_of_summand)
    return current_max

  def is_homogeneous(self) -> bool:
    """Returns true if all summands of this polynomial have the same degree."""
    R = self.__class__

    if len(self.summands) == 0:
      return True

    zeroth_summand = self.summands[0]
    degree_of_zeroth_summand = R.degree_of_summand(zeroth_summand)
    for summand in self.summands:
      if degree_of_zeroth_summand != R.degree_of_summand(summand):
        return False

    return True

  def get_homogeneous_part(self, d: int):
    """Returns the homogeneous polynomial consisting of the degree d summands."""
    R = self.__class__

    result = R.one()
    for summand in self.summands:
      degree = R.degree_of_summand(summand)
      if degree == d:
        result += R([summand])
    
    return result

  # Return a list of homogeneous summands.
  # Quicker than repeating get_homogeneous_part several times as we only iterate once.
  def get_homogeneous(self) -> list:
    """Returns a list of polynomials, each of which is a homogeneous summand of the specified polynomial.
    This is quicker than calling get_homogeneous_part several times as we only iterate once.
    """
    R = self.__class__

    summands_by_degree = {}
    for summand in self.summands:
      degree = R.degree_of_summand(summand)
      if degree not in summands_by_degree:
        summands_by_degree[degree] = [summand]
      else:
        summands_by_degree[degree].append(summand)
    
    homogeneous_polys = []
    for degree in summands_by_degree:
      homogeneous_part_for_degree = R(summands_by_degree[degree])
      homogeneous_polys.append(homogeneous_part_for_degree)
    return homogeneous_polys

  def get_summands_as_polys(self) -> list:
    """Returns a list of monomials which are the summands of the specified polynomial."""
    R = self.__class__
    return [R([summand]) for summand in self.summands]