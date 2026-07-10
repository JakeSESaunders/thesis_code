from rings.polynomial import PolyRing
from collections import Counter

class Tensor():
  """An object that represents an element in the tensor product of two rings R, S."""
  def __init__(self, summands: list[tuple]):
    """
    Keyword arguments:
    summands: a list of tuples, each tuple a "pure tensor" given by a pair of a basis element of R and a basis element of S.
    """
    if len(summands) == 0:
      self.summands = summands
      return

    R = summands[0][0].__class__ # the left tensor factor
    S = summands[0][1].__class__ # the right tensor factor

    counter = Counter(summands)
    for summand in counter:
      # check the input for duplicates
      if counter[summand] != 1: raise ValueError(f"List of summands may not contain duplicates: summand {summand} appeared {counter[summand]} times.")
      
      r, s = summand

      # check for tensor factors from rings other than R and S.
      if not isinstance(r, R): raise ValueError(f"Left tensor factor {r} of summand {summand} is not in ring {R}.")
      if not isinstance(s, S): raise ValueError(f"Right tensor factor {s} of summand {summand} is not in ring {S}.")

      # check that the provided left/right tensor factors are indeed basis elements
      if not len(r.summands) in [0, 1]: raise ValueError(f"Left tensor factor {r} of summand {summand} is not a product of generators of {R}.")
      if not len(s.summands) in [0, 1]: raise ValueError(f"Right tensor factor {s} of summand {summand} is not a product of generators of {S}.")

    self.summands = summands

  @classmethod
  def from_factors(cls, left: PolyRing, right: PolyRing):
    """Given a left and a right tensor factor, return their tensor product decomposed as a sum of pure tensors."""
    result_summands = []
    for l in left.get_summands_as_polys():
      for r in right.get_summands_as_polys():
        result_summands.append(
          (l, r)
        )
    return cls(result_summands)

  @classmethod
  def get_basis(cls, R: type[PolyRing], S: type[PolyRing], degree: int) -> list:
    """Returns an additive basis of the tensor product of R and S in the given degree."""
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
    """Return the pure tensor of the multiplicative identities of R and S."""
    return cls.from_factors(R.one(), S.one())

  @classmethod
  def zero(cls):
    """Return the pure tensor of the additive identities of R and S."""
    return cls([])

  def __eq__(self, other):
    T = self.__class__
    if not isinstance(other, T):
      return False
    for summand in self.summands:
      if summand not in other.summands:
        return False
    return True

  def __add__(self, other):
    T = self.__class__
    if not isinstance(other, T): raise ValueError("{self} and {other} are not elements of the same ring.")
    summands = self.summands.copy()
    for summand in other.summands:
      if summand in summands:
        summands.remove(summand)
      else:
        summands.append(summand)
    result = T(summands)
    return result

  def __sub__(self, other):
    return self + other # we're mod 2 so + and - are the same.

  def __mul__(self, other):
    T = self.__class__
    if not isinstance(other, self.__class__): raise ValueError("{self} and {other} are not elements of the same ring.")
    # If there are no summands, the polynomial is zero.
    if self.summands == [] or other.summands == []:
      return T.zero()

    # If self has multiple summands, multiply each summand by other and add the result.
    if len(self.summands) > 1:
      polys = self.get_summands_as_polys()
      return sum([poly * other for poly in polys], T.zero())
    
    # Do the same for other.
    if len(other.summands) > 1:
      polys = other.get_summands_as_polys()
      return sum([self * poly for poly in polys], T.zero())
    
    # Now self and other both have only a single summand.
    left_1, right_1 = self.summands[0]
    left_2, right_2 = other.summands[0]
    
    product_left = left_1*left_2
    product_right = right_1*right_2

    return T.from_factors(product_left, product_right)

  def __pow__(self, other):
    T = self.__class__
    if not isinstance(other, int): raise ValueError(f"Cannot raise polynomial to the power of non-integer value {other}.")
    if other < 0: raise ValueError(f"Cannot raise polynomial to negative power {other}.")
    if len(self.summands) == 0:
      return T.zero()
    left, right = self.summands[0]
    R, S = left.__class__, right.__class__
    result = T.one(R, S)
    for _ in range(0, other):
      result *= self
    return result

  # TODO upgrade to use polyring string code.
  # TODO include option to replace unicode with latex \otimes.
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

  def max_degree(self) -> int:
    """Returns the maximum degree of all summands of a polynomial, where the degree is the sum of the degrees of the tensor factors."""
    current_max = 0
    for summand in self.summands:
      left, right = summand
      candidate_max = left.max_degree() + right.max_degree()
      current_max = max(current_max, candidate_max)
    return current_max

  def is_homogeneous(self) -> bool:
    """Returns true if all summands of this polynomial have the same degree."""
    if len(self.summands) == 0:
      return True
    left_0, right_0 = self.summands[0]
    degree = left_0.max_degree() + right_0.max_degree()
    for summand in self.summands:
      left, right = summand
      if degree != left.max_degree() + right.max_degree():
        return False
    return True

  def get_homogeneous_part(self, d: int):
    """Returns the homogeneous polynomial consisting of the degree d summands."""
    T = self.__class__
    if len(self.summands) == 0:
      return T.zero()
    left, right = self.summands[0]
    R, S = left.__class__, right.__class__
    result = T.zero(R, S)
    for summand in self.summands:
      left, right = summand
      degree = left.max_degree() + right.max_degree()
      if degree == d:
        result += T.from_factors(left, right)
    
    return result

  # Return a list of homogeneous summands.
  def get_homogeneous(self) -> list:
    """Returns a list of polynomials, each of which is a homogeneous summand of the specified polynomial.
    This is quicker than calling get_homogeneous_part several times as we only iterate once.
    """
    T = self.__class__
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
      homogeneous_polys.append(T(summands_by_degree[degree]))
    return homogeneous_polys

  def get_summands_as_polys(self) -> list:
    """Returns a list of monomials which are the summands of the specified polynomial."""
    T = self.__class__
    return [T([summand]) for summand in self.summands]