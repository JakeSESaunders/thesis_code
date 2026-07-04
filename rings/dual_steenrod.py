from rings.polynomial import NatPolyRing
from combinatorics.indices import is_one_less_than_power_of_two
from functools import cache
from rings.tensor import Tensor

class BaseDualSteenrod(NatPolyRing):
  def is_valid_generator(i):
    return (is_one_less_than_power_of_two(i) is not None)

class DualSteenrod(BaseDualSteenrod):
  def symbol(i):
    index = is_one_less_than_power_of_two(i)
    return f"ξ_{index}"

  def xi(i):
    return DualSteenrod([(2**i - 1, )])

class DualSteenrodTensorSquare(Tensor):
  @classmethod
  def get_basis(cls, degree):
    return super().get_basis(DualSteenrod, DualSteenrod, degree)

  @classmethod
  def one(cls, R=DualSteenrod, S=DualSteenrod):
    return DualSteenrodTensorSquare([(DualSteenrod.one(), DualSteenrod.one())])

class DualSteenrodConjugate(BaseDualSteenrod):
  def symbol(i):
    index = is_one_less_than_power_of_two(i)
    return f"ζ_{index}"

  def zeta(i):
    return DualSteenrodConjugate([(2**i - 1, )])

# Returns the polynomial in the xi basis given by zeta n, according to the inductive formula given in the proof of
# Milnor - The Steenrod algebra and its dual, Lemma 10
@cache
def xi_to_zeta(n):
  if n == 0 or n == 1:
    return DualSteenrodConjugate.zeta(n)
  summands = [
    xi_to_zeta(i) * DualSteenrodConjugate.zeta(n - i)**(2**i) for i in range(1, n)
  ]
  summands.append(DualSteenrodConjugate.zeta(n))
  return sum(summands, DualSteenrodConjugate.zero())


class DualSteenrodDLOp(BaseDualSteenrod):
  def symbol(i):
    index = is_one_less_than_power_of_two(i)
    return f"Q_1^{index}xi_1"

  def Q(i):
    return DualSteenrodDLOp([(i, )])