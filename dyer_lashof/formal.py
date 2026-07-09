from combinatorics.admissible import is_tuple_of_nondecreasing_nonnegative_integers, is_tuple_of_nondecreasing_positive_integers
from dyer_lashof.adem import adem_sequence

# Operation is a list of integers, generator is an integer representing the degree of a class.
# Q_1Q_2Q_3x_4 is DL([1, 2, 3], 4)
class DL:
  def __init__(self, operation, generator, upper):
    if isinstance(operation, tuple):
      operation = list(operation)
    if not isinstance(operation, list):
      operation = [operation]
    self.operation = operation
    self.generator = generator
    self.upper = upper
    self.lower = not upper

  def from_generator(generator, upper):
    return DL([], generator, upper)

  def to_tuple(self):
    return (
      tuple(self.operation),
      self.generator,
      self.upper
    )

  def adem(self):
    if self.upper:
      return [dl_op.to_upper() for dl_op in self.to_lower().adem()]
    admissible_operations = adem_sequence(tuple(self.operation))
    return [DL(admissible_operation, self.generator, upper=False) for admissible_operation in admissible_operations]

  # If lower, checks if is a tuple of nondecreasing nonnegative integers.
  def is_admissible(self):
    if self.upper:
      return self.to_lower().is_admissible()
    return is_tuple_of_nondecreasing_nonnegative_integers(tuple(self.operation))

  # If lower, checks if is a tuple of nondecreasing positive integers
  def is_reduced(self):
    if self.upper:
      return self.to_lower().is_reduced()
    return is_tuple_of_nondecreasing_positive_integers(tuple(self.operation))

  def __eq__(self, other):
    if not isinstance(other, DL):
      return False
    return self.to_tuple() == other.to_tuple()

  def __lt__(self, other):
    return self.to_tuple() < other.to_tuple()

  def __hash__(self):
    return self.to_tuple().__hash__()

  def __str__(self):
    result = ""
    length = len(self.operation)
    for i in range(0, length):
      if self.upper:
        result += f"Q^{{{self.operation[i]}}}"
      else:
        result += f"Q_{{{self.operation[i]}}}"
    result += f"x_{{{self.generator}}}"
    return result

  # Add an operation to the front.
  def Q(self, i):
    new_operation = self.operation.copy()
    new_operation.insert(0, i)
    return DL(new_operation, self.generator, self.upper)

  def degree(self):
    if self.upper:
      return self.to_lower().degree()
    
    d = self.generator
    for i in reversed(self.operation):
      d = 2*d + i

    return d

  # Take the first operation off of the front, return the remaining DL op.
  def pop(self):
    if len(self.operation) == 0:
      raise ValueError("Cannot pop DL with no operation.")
    return (self.operation[0], DL(self.operation[1:], self.generator, self.upper))

  def to_upper(self):
    if self.upper:
      return self
    current_degree = self.generator
    new_operation = []
    for i in reversed(self.operation):
      new_operation.append(i + current_degree)
      current_degree = 2*current_degree + i
    return DL(reversed(new_operation), self.generator, upper=True)

  def to_lower(self):
    if self.lower:
      return self
    current_degree = self.generator
    new_operation = []
    for i in reversed(self.operation):
      new_operation.append(i - current_degree)
      current_degree = current_degree + i
    return DL(reversed(new_operation), self.generator, upper=False)

  # Remove squares from a Dyer--Lashof operation.
  def strip_leading_zeros(self):
    if self.upper:
      leading_zeros, shorter_op_lower = self.to_lower().strip_leading_zeros()
      return (leading_zeros, shorter_op_lower.to_upper())

    number_of_leading_zeros = 0
    for i in range(len(self.operation)):
      if self.operation[i] == 0:
        number_of_leading_zeros += 1
    return number_of_leading_zeros, DL(self.operation[number_of_leading_zeros:], self.generator, self.upper)