from rings.nat_polynomial import NatPolyRing
from rings.tensor import Tensor
from combinatorics.indices import is_one_less_than_power_of_two
from strings.manifold_string import get_manifold_basis_symbol
from config import SURROUND_INDICES_WITH_BRACES

class HomologyMO(NatPolyRing):
  @classmethod
  def symbol(cls, i: int) -> str:
    if SURROUND_INDICES_WITH_BRACES:
      return f"e_{{{i}}}"
    return f"e_{i}"

  def e(i: int):
    return HomologyMO.generator(i)

class CohomologyBO(NatPolyRing):
  @classmethod
  def symbol(cls, i: int) -> str:
    if SURROUND_INDICES_WITH_BRACES:
      return f"w_{{{i}}}"
    return f"w_{i}"

  def w(i: int):
    return CohomologyBO.generator(i)

class TensorCohomologyBO(Tensor):
  @classmethod
  def get_basis(cls, degree: int) -> list:
    return super().get_basis(CohomologyBO, CohomologyBO, degree)

  @classmethod
  def one(cls, R=CohomologyBO, S=CohomologyBO):
    return TensorCohomologyBO.from_factors(CohomologyBO.one(), CohomologyBO.one())

class HomotopyMO(NatPolyRing):
  @classmethod
  def symbol(cls, i: int) -> str:
    if not HomotopyMO.is_valid_generator(i): raise ValueError(f"Cannot find symbol for invalid generator {i}.")
    return get_manifold_basis_symbol(i)

  @classmethod
  def is_valid_generator(cls, i: int) -> bool:
    return (is_one_less_than_power_of_two(i) is None)

class HomologyMOManifoldBasis(NatPolyRing):
  @classmethod
  def symbol(cls, i: int) -> str:
    return get_manifold_basis_symbol(i)

  def M(i: int):
    return HomologyMOManifoldBasis.generator(i)

class EulerHomology(NatPolyRing):
  @classmethod
  def symbol(cls, i: int) -> str:
    if i != 2: raise ValueError(f"Cannot find symbol for invalid generator {i}.")
    return "t"

  @classmethod
  def is_valid_generator(cls, i: int) -> bool:
    return (i == 2)

  def t():
    return EulerHomology.generator(2)

# TODO
class HomologyMOBinaryBasis(NatPolyRing):
  @classmethod
  def symbol(cls, i: int) -> str:
    return f"[B_{i}]" # TODO

  def B(i: int):
    return HomologyMOBinaryBasis.generator(i)