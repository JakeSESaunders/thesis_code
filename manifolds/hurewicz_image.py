from dyer_lashof.bo import dyer_lashof_HMO
from rings.truncated_poly import TruncPoly
from combinatorics.partitions import get_ordered_partitions, get_weighted_partitions, to_partition
from combinatorics.multinomial import multinomial_mod_2
from combinatorics.indices import get_dold_manifold_indices_for_degree, is_one_less_than_power_of_two
from manifolds.sw_classes import tangential_sw_classes_dold_manifold, tangential_sw_classes_projective_space
from rings.bo import CohomologyBO, HomologyMO
from homomorphisms.w_to_e_basis import change_basis_w_to_e
from config import USE_LAGRANGE_INVERSION_FOR_HUREWICZ_IMAGE_OF_PROJECTIVE_SPACES
from dyer_lashof.monomial import DyerLashofMonomial

def hurewicz_image_projective_space_lagrange_inversion(n: int) -> HomologyMO:
  """Uses the multinomial formula for Lagrange inversion to compute the Hurewicz image of RP^n in the e-basis of H_*MO."""
  if n % 2 != 0: raise ValueError("{n} is not even, so this calculation gives an unintended result.")
  if n <= 0: raise ValueError("{n} is not positive, so this calculation gives an unintended result.")

  result = HomologyMO.zero()
  for partition in get_weighted_partitions(n):
    ll = partition + (n, )
    if multinomial_mod_2(ll) == 1:
      result += HomologyMO.from_partition(to_partition(partition, n))
  return result

def from_total_sw_class_to_sw_number_polynomial(total_sw_class: TruncPoly) -> CohomologyBO:
  """Given a class in a truncated polynomial ring representing the total Stiefel--Whitney class of a manifold, return the dual of the corresponding polynomial in the w-basis of H^*BO."""
  result = CohomologyBO.zero()
  dimension = total_sw_class.signature.dimension()
  summands_of_degree = [total_sw_class.get_summands_of_degree(i) for i in range(0, dimension + 1)]
  # Each partition of the dimension of the manifold encodes an SW class
  for partition in get_ordered_partitions(dimension):
    product = TruncPoly.one(total_sw_class.signature)
    for part in partition:
      product *= summands_of_degree[part]
    if not product.is_zero():
      result += CohomologyBO.from_partition(partition)
  return result

def hurewicz_image_via_BO_coaction(tangential_sw_class: TruncPoly) -> HomologyMO:
  """Compute the Hurewicz image of a manifold with the given tangential Stiefel--Whitney class using the comultiplication in BO (this is slow)."""
  normal_sw_class = tangential_sw_class.inverse()
  normal_sw_poly = from_total_sw_class_to_sw_number_polynomial(normal_sw_class)
  return change_basis_w_to_e(normal_sw_poly)

def hurewicz_image_projective_space_via_BO_coaction(n: int) -> HomologyMO:
  """Compute the Hurewicz image of a projective space using the comultiplication in BO (this is slow)."""
  if n % 2 != 0: raise ValueError("Dimension {n} is odd and therefore the calculation gives an unintended result.")
  tangential_sw_class = tangential_sw_classes_projective_space(n)
  return hurewicz_image_via_BO_coaction(tangential_sw_class)

def hurewicz_image_dold_manifold_via_BO_coaction(n: int) -> HomologyMO:
  """Compute the Hurewicz image of a Dold manifold space using the comultiplication in BO (this is slow)."""
  if n % 2 == 0:
    raise ValueError("Dimension {n} is even and therefore invalid.")
  r, s = get_dold_manifold_indices_for_degree(n)
  tangential_sw_class = tangential_sw_classes_dold_manifold(r, s) 
  return hurewicz_image_via_BO_coaction(tangential_sw_class)

def hurewicz_image_even(n: int) -> HomologyMO:
  """Compute the Hurewicz image of an even dimensional projective space."""
  if USE_LAGRANGE_INVERSION_FOR_HUREWICZ_IMAGE_OF_PROJECTIVE_SPACES:
    return hurewicz_image_projective_space_lagrange_inversion(n)
  return hurewicz_image_projective_space_via_BO_coaction(n)

def hurewicz_image_odd(n: int) -> HomologyMO:
  """Compute the Hurewicz image of an odd dimensional manifold, or Q_1^i in the e-basis if applicable."""
  index = is_one_less_than_power_of_two(n)
  if index is None:
    return hurewicz_image_dold_manifold_via_BO_coaction(n)
  Q1index = DyerLashofMonomial(tuple([1 for i in range(index - 1)]), upper=False)
  return dyer_lashof_HMO(Q1index, HomologyMO.e(1))

def hurewicz_image(n: int) -> HomologyMO:
  """Return the Hurewicz image of the generating manifold in degree n in the e-basis of H_*MO."""
  return hurewicz_image_even(n) if (n % 2 == 0) else hurewicz_image_odd(n)