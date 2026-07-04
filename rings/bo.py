from rings.polynomial import NatPolyRing
from rings.tensor import Tensor
from combinatorics.indices import is_one_less_than_power_of_two, get_dold_manifold_indices_for_degree
from config import DISPLAY_DOLD_MANIFOLD_INDICES_MANIFOLD_BASIS

class HomologyMO(NatPolyRing):
  def symbol(i):
    return f"e_{i}"

  def e(i):
    return HomologyMO.generator(i)

class EulerHomology(NatPolyRing):
  def symbol(i):
    if i != 2:
      raise ValueError(f"Cannot find symbol for invalid generator {i}.")
    return "t"

  def is_valid_generator(i):
    return (i == 2)

  def t():
    return EulerHomology.generator(2)

class CohomologyBO(NatPolyRing):
  def symbol(i):
    return f"w_{i}"

  def w(i):
    return CohomologyBO.generator(i)

class TensorCohomologyBO(Tensor):
  @classmethod
  def get_basis(cls, degree):
    return super().get_basis(CohomologyBO, CohomologyBO, degree)

  @classmethod
  def one(cls, R=CohomologyBO, S=CohomologyBO):
    return TensorCohomologyBO([(CohomologyBO.one(), CohomologyBO.one())])

def get_manifold_basis_symbol(i):
  if i % 2 == 0:
    return f"[RP^{i}]"
  index = is_one_less_than_power_of_two(i)
  if index is None:
    if DISPLAY_DOLD_MANIFOLD_INDICES_MANIFOLD_BASIS:
      r, s = get_dold_manifold_indices_for_degree(i)
      return f"[P({r}, {s})]"
    return f"[M_{i}]"
  if i == 1:
    return "(e_1)"
  if i == 3:
    return "(Q_1e_1)"
  return f"(Q_1^{index - 1}e_1)"

class HomotopyMO(NatPolyRing):
  def symbol(i):
    return get_manifold_basis_symbol(i)

  def is_valid_generator(i):
    return (is_one_less_than_power_of_two(i) is None)

class HomologyMOManifoldBasis(NatPolyRing):
  def symbol(i):
    return get_manifold_basis_symbol(i)

  def M(i):
    return HomologyMOManifoldBasis.generator(i)

class HomologyMOBinaryBasis(NatPolyRing):
  def symbol(i):
    return f"[B_{i}]" # TODO

  def B(i):
    return HomologyMOBinaryBasis.generator(i)