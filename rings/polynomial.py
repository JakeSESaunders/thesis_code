from combinatorics.partitions import get_ordered_partitions
from functools import cache

# Base class for a graded polynomial ring with addition mod 2.
# Objects of this class are elements of the polynomial ring.
# "summands" is a list of sorted tuples representing summands of a polynomial.
# It should contain no duplicates as we are working mod 2, and this could break equality!
# The default implementation sorts the list of summands and each tuple of factors, but does not check for duplicates.
# We instead check this when performing arithmetic operations.
# Each element of the tuple should be a "generator", a hashable type with equality and a comparison (le/ge) dunder method, such that is_valid_generator returns True. 
# By convention, zero has summands [] and one has summands [()].
class PolyRing:
  def __init__(self, summands):
    new_summands = []
    for summand in summands:
      new_summand = []
      for factor in summand:
        if not self.__class__.is_valid_generator(factor):
          raise ValueError(f"Generator with index {factor} is not valid.")
        new_summand.append(factor)
      new_summands.append(tuple(sorted(new_summand)))
    self.summands = sorted(new_summands)

  def is_valid_generator(x):
    return True

  @classmethod
  @cache
  def get_basis(cls, degree):
    raise NotImplementedError("Base class PolyRing does not implement get_basis.") 

  @classmethod
  def generator(cls, i):
    if i == 0:
      return cls.one()
    return cls([(i, )])

  # given generator x, return the degree of x.
  @classmethod
  def degree(cls, x):
    if not cls.is_valid_generator(x):
      raise ValueError(f"Cannot determine degree of invalid generator {x}.")
    raise NotImplementedError("Base class PolyRing does not implement degree.")

  # when supplied with a tuple of generators, return the sum of their degrees, giving the degree of the summand.
  # NOTE if the summand is empty, i.e. is 0, then this returns degree 0 rather than -\infty.
  @classmethod
  def degree_of_summand(cls, summand):
    degree_of_summand = 0
    for factor in summand:
      degree_of_summand += cls.degree(factor)
    return degree_of_summand

  @classmethod
  def one(cls):
    return cls([()])

  @classmethod
  def zero(cls):
    return cls([])

  # NOTE vulnerable to poorly formed objects with e.g. summands that appear twice or have extra zeros.
  # TODO this could probably be faster.
  def __eq__(self, other):
    if not isinstance(other, self.__class__):
      return False
    
    for summand in self.summands + other.summands:
      if (summand not in other.summands) or (summand not in self.summands):
        return False

    return True

  def __add__(self, other):
    if not isinstance(other, self.__class__):
      raise ValueError("{self} and {other} are not elements of the same ring.")
    summands = self.summands.copy()
    for summand in other.summands:
      if summand in summands:
        summands.remove(summand)
      else:
        summands.append(summand)
    return self.__class__(summands)

  def __sub__(self, other):
    return self + other # we're mod 2 so + and - are the same.

  def __mul__(self, other):
    cls = self.__class__
    if not isinstance(other, cls):
      raise ValueError("{self} and {other} are not elements of the same ring.")
    # If there are no summands, the polynomial is zero.
    if self.summands == [] or other.summands == []:
      return cls.zero()

    # If self has multiple summands, multiply each summand by other and add the result.
    if len(self.summands) > 1:
      polys = self.get_summands_as_polys()
      return sum([poly * other for poly in polys], cls.zero())
    
    # Do the same for other.
    if len(other.summands) > 1:
      polys = other.get_summands_as_polys()
      return sum([self * poly for poly in polys], cls.zero())
    
    # Now self and other both have only a single summand.
    summand_1 = self.summands[0]
    summand_2 = other.summands[0]
    
    product = tuple(sorted(summand_1 + summand_2))
    return cls([product])

  def __pow__(self, other):
    if not isinstance(other, int):
      raise ValueError(f"Cannot raise polynomial to the power of non-integer value {other}.")
    result = self.__class__.one()
    for i in range(0, other):
      result *= self
    return result

  # Returns the string symbol to use for the generator x.
  def symbol(x):
    raise NotImplementedError("Function symbol not implemented for base PolyRing class.")

  def __str__(self):
    if len(self.summands) == 0:
      return "0"

    result = ""
    # Add a string part for each summand.
    for summand in self.summands:
      result_summand = ""
      for factor in summand:
        result_summand += self.__class__.symbol(factor) # Add a string part for each factor.
      if len(summand) == 0:
        result_summand = "1"
      if len(result) == 0:
        result = result_summand
      else:
        result += f" + {result_summand}"
    return result

  def __hash__(self):
    return f"{self.__class__.__name__}, {self.__str__()}".__hash__()

  def max_degree(self):
    cls = self.__class__
    current_max = 0
    for summand in self.summands:
      degree_of_summand = cls.degree_of_summand(summand)
      current_max = max(current_max, degree_of_summand)
    return current_max

  def is_homogeneous(self):
    cls = self.__class__

    if len(self.summands) == 0:
      return True

    zeroth_summand = self.summands[0]
    degree_of_zeroth_summand = cls.degree_of_summand(zeroth_summand)
    for summand in self.summands:
      if degree_of_zeroth_summand != cls.degree_of_summand(summand):
        return False

    return True

  def get_homogeneous_part(self, d):
    cls = self.__class__

    result = cls.one()
    for summand in self.summands:
      degree = cls.degree_of_summand(summand)
      if degree == d:
        result += cls([summand])
    
    return result

  # Return a list of homogeneous summands.
  # Quicker than repeating get_homogeneous_part several times as we only iterate once.
  def get_homogeneous(self):
    cls = self.__class__

    summands_by_degree = {}
    for summand in self.summands:
      degree = cls.degree_of_summand(summand)
      if degree not in summands_by_degree:
        summands_by_degree[degree] = [summand]
      else:
        summands_by_degree[degree].append(summand)
    
    homogeneous_polys = []
    for degree in summands_by_degree:
      homogeneous_part_for_degree = cls(summands_by_degree[degree])
      homogeneous_polys.append(homogeneous_part_for_degree)
    return homogeneous_polys

  def get_summands_as_polys(self):
    cls = self.__class__
    return [cls([summand]) for summand in self.summands]

# Base class for polynomial rings with at most one generator in each positive integer degree.
class NatPolyRing(PolyRing):
  def __init__(self, summands):
    new_summands = []
    for summand in summands:
      new_summand = []
      for factor in summand:
        if factor == 0:
          continue # NOTE If every factor is 0, then we add an empty tuple to the summands, and this represents 1, as expected, since e.g. e_0 = 1.
        if factor < 0:
          raise ValueError("Negative indices are not supported.")
        if not self.__class__.is_valid_generator(factor):
          raise ValueError(f"Generator with index {factor} is not valid.")
        new_summand.append(factor)
      new_summands.append(tuple(sorted(new_summand)))
    self.summands = sorted(new_summands)

  def is_valid_generator(i):
    if not isinstance(i, int):
      return False
    if i < 0:
      return False
    return True # NOTE we deal with generator with index 0 appropriately

  @classmethod
  def degree(cls, x):
    if not cls.is_valid_generator(x):
      raise ValueError("Cannot get degree of non-integer value {x}.")
    return x

  @classmethod
  @cache
  def get_basis(cls, degree):
    basis = []
    for partition in get_ordered_partitions(degree):
      partition_is_valid = True
      for part in partition:
        if not cls.is_valid_generator(part):
          partition_is_valid = False
          break
      if partition_is_valid:
        basis.append(
          cls.from_partition(partition)
        )
    return basis

  @classmethod
  def from_partition(cls, partition):
    return cls([partition])
