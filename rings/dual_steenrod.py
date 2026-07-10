from rings.nat_polynomial import NatPolyRing
from combinatorics.indices import is_one_less_than_power_of_two
from functools import cache
from rings.tensor import Tensor
from strings.dual_steenrod_string import get_dual_steenrod_symbol

class BaseDualSteenrod(NatPolyRing):
  @classmethod
  def is_valid_generator(cls, i: int) -> bool:
    return (is_one_less_than_power_of_two(i) is not None)

class DualSteenrod(BaseDualSteenrod):
  @classmethod
  def symbol(cls, i: int) -> str:
    return get_dual_steenrod_symbol(True, i)

  def xi(i: int):
    return DualSteenrod.generator(2**i - 1)

class DualSteenrodConjugate(BaseDualSteenrod):
  @classmethod
  def symbol(cls, i: int) -> str:
    return get_dual_steenrod_symbol(False, i)

  def zeta(i: int):
    return DualSteenrodConjugate.generator(2**i - 1)

class DualSteenrodTensorSquare(Tensor):
  @classmethod
  def get_basis(cls, degree: int) -> list:
    return super().get_basis(DualSteenrod, DualSteenrod, degree)

  @classmethod
  def one(cls, R=DualSteenrod, S=DualSteenrod):
    return DualSteenrodTensorSquare.from_factors(DualSteenrod.one(), DualSteenrod.one())

# 
@cache
def xi_to_zeta(n: int):
  """Returns the polynomial in the xi basis given by zeta n, according to the inductive formula given in the proof of Milnor - The Steenrod algebra and its dual, Lemma 10."""
  if n == 0 or n == 1:
    return DualSteenrodConjugate.zeta(n)
  summands = [
    xi_to_zeta(i) * DualSteenrodConjugate.zeta(n - i)**(2**i) for i in range(1, n)
  ]
  summands.append(DualSteenrodConjugate.zeta(n))
  return sum(summands, DualSteenrodConjugate.zero())

class DualSteenrodDLOp(BaseDualSteenrod):
  @classmethod
  def symbol(cls, i: int) -> str:
    index = is_one_less_than_power_of_two(i)
    return f"Q_1^{index}xi_1" # TODO replace with DL monomial string generation code.

  def Q(i: int):
    """Return the generator Q_1^i\xi_1."""
    return DualSteenrodDLOp.generator(2**(i + 1) - 1)