from collections import Counter
from itertools import product
from config import DEBUG

class TruncPolySignature:
  # A truncated polynomial ring has n generators X_0, ..., X_{n - 1}.
  # The list degree_of_generator has length n, and the ith entry is the degree |X_i|.
  # The list order_of_generator has length n, and the ith entry is the minimum j such that X_i**j = 0.
  # We call this the _signature_ of the truncated polynomial ring
  # Aliases is a list of strings. The ith string determines the label of the ith generator in string conversions.
  def __init__(self, degree_of_generator, order_of_generator, aliases=None):
    if len(degree_of_generator) != len(order_of_generator):
      raise ValueError("Supplied lists are of different lengths, so cannot specify a truncated polynomial ring.")
    if len(aliases) != len(degree_of_generator):
      raise ValueError("List of aliases is incorrect length.")
    self.degree_of_generator = degree_of_generator
    self.order_of_generator = order_of_generator
    self.aliases = aliases

  def number_of_generators(self):
    return len(self.degree_of_generator)

  def dimension(self):
    total = 0
    for i in range(self.number_of_generators()):
      total += self.degree_of_generator[i] * (self.order_of_generator[i] - 1)
    return total

  # Return a list of tuples corresponding to an additive basis of the truncated poly ring.
  # Note the list is sorted such that, if j > 0, then l[i + j] is never a factor of l[i].
  def basis(self):
    return product(*[range(self.order_of_generator[i]) for i in range(self.number_of_generators())])

  # Returns true only if the signature data is exactly equal.
  # E.g. if |X_0| = 1, |X_1| = 2 and both have order 2, this is not equal to
  # a signature with |X_1| = 1, |X_0| = 2 for both generators of order 2.
  def __eq__(self, other):
    if len(self.degree_of_generator) != len(other.degree_of_generator):
      return False
    for i in range(len(self.degree_of_generator)):
      if self.degree_of_generator[i] != other.degree_of_generator[i]:
        return False
      if self.order_of_generator[i] != other.order_of_generator[i]:
        return False
    return True
  
  def is_summand_zero(self, summand):
    power_in_summand = Counter(summand)
    for i in power_in_summand.keys():
      if self.order_of_generator[i] <= power_in_summand[i]:
        return True
    return False

  def degree_of_summand(self, summand):
    total_degree = 0
    for generator in summand:
      total_degree += self.degree_of_generator[generator]
    return total_degree

# Represents an element of a truncated polynomial ring over Z/2.
# The argument _summands_ is a list of tuples representing the summands of the polynomial.
# It should be sorted, and each tuple element of the list should also be sorted.
# A list element of the form (0, 1, 1, 2, 2, 2) represents a summand X_0X_1^2X_2^3.
class TruncPoly:
  def __init__(self, signature: TruncPolySignature, summands: list[tuple[int]]):
    for summand in summands:
      if signature.is_summand_zero(summand):
        raise ValueError("Cannot instantiate truncated polynomial with zero summand {summand}.")
    self.signature = signature
    self.summands = summands

  def copy(self):
    return TruncPoly(self.signature, self.summands.copy())

  def is_zero(self):
    return len(self.summands) == 0

  def zero(signature: TruncPolySignature):
    return TruncPoly(signature, [])

  def one(signature: TruncPolySignature):
    return TruncPoly(signature, [()])

  def __str__(self):
    result = ""
    if self.is_zero():
      return "0"
    for summand in self.summands:
      result_summand = ""
      for factor in summand:
        if self.signature.aliases is not None:
          result_summand += f"{self.signature.aliases[factor]}"
        else:
          result_summand += f"X_{factor}"
      if len(summand) == 0:
        result_summand = "1"
      if len(result) == 0:
        result += result_summand
      else:
        result += f" + {result_summand}"
    return result

  def __add__(self, other):
    if self.signature != other.signature:
      raise ValueError("Cannot add polynomials in different rings.")
    result_summands = list(self.summands)
    for summand in other.summands:
      if summand in result_summands:
        result_summands.remove(summand)
      else:
        result_summands.append(summand)
    result_summands = tuple(sorted(result_summands))
    return TruncPoly(self.signature, result_summands)

  def __sub__(self, other):
    return self + other # we're working mod 2 so they're the same.

  def __mul__(self, other):
    if self.signature != other.signature:
      raise ValueError("Cannot multiply polynomials in different rings.")
    
    if self.is_zero():
      return TruncPoly.zero(self.signature)
    if other.is_zero():
      return TruncPoly.zero(self.signature)

    if len(self.summands) > 1:
      return sum([poly * other for poly in self.get_summands_as_polys()], TruncPoly.zero(self.signature))
    if len(other.summands) > 1:
      return sum([self * poly for poly in other.get_summands_as_polys()], TruncPoly.zero(self.signature))

    summand_1 = self.summands[0]
    summand_2 = other.summands[0]
    product = tuple(sorted(summand_1 + summand_2))
    if self.signature.is_summand_zero(product):
      return TruncPoly.zero(self.signature)
    return TruncPoly(self.signature, [product])

  def __pow__(self, other):
    if not isinstance(other, int):
      raise ValueError(f"Cannot raise truncated polynomial to non-integer value {other}.")
    if other < 0:
      return (self.inverse())**(-int)
    if other == 0:
      return TruncPoly.one(self.signature)
    if other == 1:
      return self
    return self * (self)**(other - 1)

  def get_summands_as_polys(self):
    return [TruncPoly(self.signature, [summand]) for summand in self.summands]

  def get_summands_of_degree(self, d):
    summands_of_degree_d = [summand for summand in self.summands if self.signature.degree_of_summand(summand) == d]
    return TruncPoly(self.signature, summands_of_degree_d)

  def inverse(self):
    # print(f"Computing inverse of {self}...")
    current_inverse = TruncPoly.zero(self.signature)
    # We want to find a truncated poly which multiplies with self to give 1. 1 is the "desired product".
    desired_product = TruncPoly.one(self.signature)
    current_product = TruncPoly.zero(self.signature)

    indexing_set = self.signature.basis()
    # Convert the indices to a monomial, e.g. (1, 3, 2) -> (0, 1, 1, 1, 2, 2) for X_0X_1^3X_2^2.
    for indices in indexing_set:
      monomial = []
      for i in range(len(indices)):
        monomial += [i for j in range(indices[i])]
      monomial = tuple(sorted(monomial))
      if monomial in (current_product - desired_product).summands:
        monomial_poly = TruncPoly(self.signature, [monomial])
        current_inverse += monomial_poly
        current_product += self * monomial_poly

    # as a sanity check, determine whether self*current_inverse is the desired product.
    if not (self*current_inverse - desired_product).is_zero():
      raise RuntimeError("Calculated inverse is not actual inverse.")

    if DEBUG:
      print("Computed!")

    return current_inverse
    
