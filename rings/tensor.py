from typing import TypeVar, Generic
from rings.polynomial import PolyRing

# An object that represents an element in the tensor product of two rings R, S.
# Summands is a list of summands, each summand a pure tensor given by a basis element of R and a basis element of S.
class Tensor():
  def __init__(self, summands):
    self.summands = summands

  # Given a left and a right tensor factor, return the sum of pure tensors.
  @classmethod
  def from_factors(cls, left: PolyRing, right: PolyRing):
    result_summands = []
    for l in left.get_summands_as_polys():
      for r in right.get_summands_as_polys():
        result_summands.append(
          (l, r)
        )
    return cls(result_summands)

  @classmethod
  def get_basis(cls, R, S, degree):
    basis = []
    for i in range(0, degree + 1):
      j = degree - i
    
      R_basis = R.get_basis(i)
      S_basis = S.get_basis(j)

      for v in R_basis:
        for w in S_basis:
          basis.append(
            cls.from_factors(v, w)
          )

    return basis

  @classmethod
  def one(cls, R, S):
    return cls.from_factors(R.one(), S.one())

  @classmethod
  def zero(cls):
    return cls([])

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
    result = self.__class__(summands)
    return result

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
    left_1, right_1 = self.summands[0]
    left_2, right_2 = other.summands[0]
    
    product_left = left_1*left_2
    product_right = right_1*right_2

    return self.__class__([(product_left, product_right)])

  def __pow__(self, other):
    if not isinstance(other, int):
      raise ValueError(f"Cannot raise polynomial to the power of non-integer value {other}.")
    if len(self.summands) == 0:
      return self.__class__.zero()
    left, right = self.summands[0]
    R, S = left.__class__, right.__class__
    result = self.__class__.one(R, S)
    for i in range(0, other):
      result *= self
    return result

  def __str__(self):
    if len(self.summands) == 0:
      return "0"
    string_result = ""
    for summand in self.summands:
      if string_result != "":
        string_result += " + "
      left, right = summand
      string_result += f"{left.__str__()}⊗{right.__str__()}"
    return string_result

  def __hash__(self):
    return f"{self.__class__.__name__}, {self.__str__()}".__hash__()

  def max_degree(self):
    current_max = 0
    for summand in self.summands:
      left, right = summand
      candidate_max = left.max_degree() + right.max_degree()
      current_max = max(current_max, candidate_max)
    return current_max

  def is_homogeneous(self):
    if len(self.summands) == 0:
      return True
    left_0, right_0 = self.summands[0]
    degree = left_0.max_degree() + right_0.max_degree()
    for summand in self.summands:
      left, right = summand
      if degree != left.max_degree() + right.max_degree():
        return False
    return True

  def get_homogeneous_part(self, d):
    if len(self.summands) == 0:
      return self.__class__.zero()
    left, right = self.summands[0]
    R, S = left.__class__, right.__class__
    result = self.__class__.one(R, S)
    for summand in self.summands:
      left, right = summand
      degree = left.max_degree() + right.max_degree()
      if degree == d:
        result += self.__class__([summand])
    
    return result

  # Return a list of homogeneous summands.
  def get_homogeneous(self):
    summands_by_degree = {}
    for summand in self.summands:
      left, right = summand
      degree = left.max_degree() + right.max_degree()
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