from combinatorics.partitions import get_ordered_partitions
from functools import cache

# Base class for a graded polynomial ring with addition mod 2.
# Objects of this class are elements of the polynomial ring.
# Current limitation means only classes of 
class PolyRing:
  def __init__(self, summands):
    new_summands = []
    for summand in summands:
      new_summand = []
      for factor in summand:
        if factor == 0:
          continue
        if factor < 0:
          raise ValueError("Negative indices are not supported.")
        if not self.__class__.is_valid_generator(factor):
          raise ValueError(f"Generator with index {factor} is not valid.")
        new_summand.append(factor)
      new_summands.append(tuple(sorted(new_summand)))
    self.summands = sorted(new_summands)

  def is_valid_generator(i):
    return True

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
  def generator(cls, i):
    if i == 0:
      return cls.one()
    return cls([(i, )])

  @classmethod
  def from_partition(cls, partition):
    return cls([partition])

  @classmethod
  def one(cls):
    return cls([()])

  @classmethod
  def zero(cls):
    return cls([])

  # NOTE vulnerable to poorly formed objects with e.g. summands that appear twice or have extra zeros.
  def __eq__(self, other):
    if not isinstance(other, self.__class__):
      return False
    for summand in self.summands:
      if summand not in other.summands:
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
    if not isinstance(other, self.__class__):
      raise ValueError("{self} and {other} are not elements of the same ring.")
    # If there are no summands, the polynomial is zero.
    if self.summands == [] or other.summands == []:
      return self.__class__.zero()

    # If self has multiple summands, multiply each summand by other and add the result.
    if len(self.summands) > 1:
      polys = self.get_summands_as_polys()
      return sum([poly * other for poly in polys], self.__class__.zero())
    
    # Do the same for other.
    if len(other.summands) > 1:
      polys = other.get_summands_as_polys()
      return sum([self * poly for poly in polys], self.__class__.zero())
    
    # Now self and other both have only a single summand.
    summand_1 = self.summands[0]
    summand_2 = other.summands[0]
    
    product = tuple(sorted(summand_1 + summand_2))
    return self.__class__([product])

  def __pow__(self, other):
    if not isinstance(other, int):
      raise ValueError(f"Cannot raise polynomial to the power of non-integer value {other}.")
    result = self.__class__.one()
    for i in range(0, other):
      result *= self
    return result

  # Returns the string symbol to use for the generator of degree i
  def symbol(i):
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
    current_max = 0
    for summand in self.summands:
      candidate_max = sum(summand)
      current_max = max(current_max, candidate_max)
    return current_max

  def is_homogeneous(self):
    if len(self.summands) == 0:
      return True
    degree = sum(self.summands[0])
    for summand in self.summands:
      if degree != sum(summand):
        return False
    return True

  def get_homogeneous_part(self, d):
    result = self.__class__.one()
    for summand in self.summands:
      degree = sum(summand)
      if degree == d:
        result += self.__class__([summand])
    
    return result

  # Return a list of homogeneous summands.
  def get_homogeneous(self):
    summands_by_degree = {}
    for summand in self.summands:
      degree = sum(summand)
      if degree not in summands_by_degree:
        summands_by_degree[degree] = [summand]
      else:
        summands_by_degree[degree].append(summand)
    
    homogeneous_polys = []
    for degree in summands_by_degree:
      homogeneous_polys.append(self.__class__(summands_by_degree[degree]))
    return homogeneous_polys

  def get_summands_as_polys(self):
    return [self.__class__([summand]) for summand in self.summands]