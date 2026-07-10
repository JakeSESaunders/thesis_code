from combinatorics.admissible import is_tuple_of_nondecreasing_nonnegative_integers, is_tuple_of_nondecreasing_positive_integers
from dyer_lashof.adem import adem_sequence

# Operation is a list of integers, generator is an integer representing the degree of a class.
# x_4 is DL
class DyerLashofMonomial:
  """A monomial in the Dyer--Lashof operations, represented by a tuple of integers and a bool upper indicating whether we are using upper or lower indexing notation.
  For instance, DyerLashofMonomial((1, 2, 2, 3), upper=False) represents Q_1Q_2Q_2Q_3
  """
  def __init__(self, operation: tuple[int] | int, upper: bool):
    """Keyword arguments:
    operation: a tuple of integers corresponding to the indices of the operations. For example, (1, 2, 2, 3) represents Q_1Q_2Q_2Q_3.
    """
    if not isinstance(operation, tuple):
      operation = (operation, )
    self.operation = operation
    self.upper = upper
    self.lower = not upper # shorthand.

  # See Curtis --- The Dyer--Lashof algebra and the \Lambda-algebra p233.
  def length(self) -> int:
    """Get the length of the operation (the number of factors of the monomial)."""
    return len(self.operation)

  # See Curtis --- The Dyer--Lashof algebra and the \Lambda-algebra p233.
  # Note that DL ops in his paper act from the right instead of from the left, like they do here.
  def excess(self) -> int:
    if len(self.operation) == 0:
      return 0
    if len(self.operation) == 1:
      return self.operation[0]
    return self.operation[0] - (sum(self.operation[1:]))

  def empty(upper: bool):
    """The empty Dyer--Lashof monomial"""
    return DyerLashofMonomial((), upper)

  def single(k: int, upper: bool):
    return DyerLashofMonomial((k, ), upper)

  def to_tuple(self):
    """Returns a tuple uniquely identifying the Dyer--Lashof monomial."""
    return (
      self.operation,
      self.upper
    )

  def to_upper(self, d: int):
    """Convert the monomial to upper indexing notation, based on the degree d of the generator to which the operation is being applied."""
    if self.upper:
      return self
    current_degree = d
    new_operation = []
    for i in reversed(self.operation):
      new_operation.append(i + current_degree)
      current_degree = 2*current_degree + i
    return DyerLashofMonomial(reversed(new_operation), upper=True)

  def to_lower(self, d: int):
    """Convert the monomial to lower indexing notation, based on the degree d of the generator to which the operation is being applied."""
    if self.lower:
      return self
    current_degree = d
    new_operation = []
    for i in reversed(self.operation):
      new_operation.append(i - current_degree)
      current_degree = current_degree + i
    return DyerLashofMonomial(reversed(new_operation), upper=False)

  def Q(self, k: int):
    """Add a new operation to the front of this one. For example, (Q_3Q_4^2).Q(2) gives Q_2Q_3Q_4^2."""
    new_operation = (k, ) + self.operation
    return DyerLashofMonomial(new_operation, self.upper)

  def degree(self, d: int) -> int:
    """If x is a class of degree d, returns the degree of the operation applied to x.
    For example, Q_1Q_2.degree(5) = 2*(2*5 + 2) + 1.
    """
    if self.upper:
      return self.to_lower().degree()
    
    for i in reversed(self.operation):
      d = 2*d + i

    return d

  def pop(self) -> tuple:
    """Take the first operation off of the front, return the remaining DL op. Returns a pair (Qk, Qj) of operations."""
    if len(self.operation) == 0:
      raise ValueError("Cannot pop the empty Dyer--Lashof monomial.")
    return (DyerLashofMonomial.single(self.operation[0], self.upper), DyerLashofMonomial(self.operation[1:], self.upper))

  # Remove squares from a Dyer--Lashof operation.
  def strip_leading_zeros(self, d: int = -1):
    """Remove any leading zeros from the operation, or leading "squares" if the monomial uses upper indexing notation.
    
    Keyword arguments:
    d: the degree of the class on which the operator is operating. This must be supplied if the monomial is in upper indexing notation.
    """
    if self.upper:
      if d == -1: raise ValueError(f"Must specify degree to strip leading squares from Dyer--Lashof monomial {self} in upper indexing notation.")
      leading_zeros, shorter_op_lower = self.to_lower(d).strip_leading_zeros()
      return (leading_zeros, shorter_op_lower.to_upper(d))

    number_of_leading_zeros = 0
    for i in range(len(self.operation)):
      if self.operation[i] == 0:
        number_of_leading_zeros += 1
    return number_of_leading_zeros, DyerLashofMonomial(self.operation[number_of_leading_zeros:], self.upper)

  def adem(self) -> list:
    if self.upper:
      excess = self.excess()
      return [dl_op.to_upper(excess) for dl_op in self.to_lower(excess).adem()]

    admissible_operations = adem_sequence(tuple(self.operation))
    return [DyerLashofMonomial(admissible_operation, upper=False) for admissible_operation in admissible_operations]

  # See Curtis --- The Dyer--Lashof algebra and the \Lambda-algebra p233.
  # Note that DL ops in his paper act from the right instead of from the left, like they do here.
  def is_admissible(self, d: int = -1) -> bool:
    """If lower, checks if the monomial is a tuple of nondecreasing nonnegative integers.
    If upper, checks if 2i_{j + 1} >= i_j for Q^{i_1}Q^{i_2}...Q^{i_n}.
    Keyword arguments:
    d: the degree of the class on which the operator is operating. This must be supplied if the monomial is in upper indexing notation."""
    if self.upper:
      if d == -1: raise ValueError(f"Must specify degree to check if upper-indexed Dyer--Lashof monomial {self} in upper indexing notation is admissible.")
      return self.to_lower(d).is_admissible()
    return is_tuple_of_nondecreasing_nonnegative_integers(tuple(self.operation))

  def is_reduced(self, d: int = -1) -> bool:
    """If lower, checks if the monomial is a tuple of nondecreasing positive integers.
    
    Keyword arguments:
    d: the degree of the class on which the operator is operating. This must be supplied if the monomial is in upper indexing notation."""
    if self.upper:
      if d == -1: raise ValueError(f"Must specify degree to check if upper-indexed Dyer--Lashof monomial {self} in upper indexing notation is reduced.")
      return self.to_lower().is_reduced()
    return is_tuple_of_nondecreasing_positive_integers(tuple(self.operation))

  def __eq__(self, other):
    if not isinstance(other, DyerLashofMonomial): return False
    return self.to_tuple() == other.to_tuple()

  def __lt__(self, other):
    if not isinstance(other, DyerLashofMonomial): raise ValueError(f"Cannot compare Dyer--Lashof monomial {self} with non-monomial {other}.")
    if self.upper != other.upper: raise ValueError(f"Cannot compare Dyer--Lashof monomials {self} and {other} which use different indexing notation.")
    return self.to_tuple() < other.to_tuple()

  def __hash__(self):
    return self.to_tuple().__hash__()

  def __str__(self): # TODO implement config.
    result = ""
    length = len(self.operation)
    for i in range(0, length):
      if self.upper:
        result += f"Q^{{{self.operation[i]}}}"
      else:
        result += f"Q_{{{self.operation[i]}}}"
    return result

  def __add__(self, other):
    """Concatenate Dyer--Lashof monomials."""
    if not isinstance(other, DyerLashofMonomial):
      raise ValueError(f"Cannot concatenate Dyer--Lashof monomial {self} to non-Dyer--Lashof monomial {other}.")
    if not (self.upper == other.upper):
      raise ValueError(f"Cannot concatenate Dyer--Lashof monomials using differing indexing notations.")
    return DyerLashofMonomial(self.operation + other.operation, self.upper)