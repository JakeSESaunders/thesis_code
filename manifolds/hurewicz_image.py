from dyer_lashof.bo import dyer_lashof
from rings.truncated_poly import TruncPoly
from combinatorics.partitions import get_ordered_partitions, get_weighted_partitions, to_partition
from combinatorics.multinomial import multinomial_mod_2
from combinatorics.indices import get_dold_manifold_indices_for_degree, is_one_less_than_power_of_two
from manifolds.sw_classes import tangential_sw_classes_dold_manifold, tangential_sw_classes_projective_space
from rings.bo import CohomologyBO, HomologyMO
from homomorphisms.w_to_e_basis import change_basis_w_to_e
from config import USE_LAGRANGE_INVERSION_FOR_HUREWICZ_IMAGE_OF_PROJECTIVE_SPACES

# Uses the multinomial formula for Lagrange inversion to compute the H
def hurewicz_image_projective_space_lagrange_inversion(n):
  if n % 2 != 0:
    raise ValueError("{n} is not even, so this calculation is erroneous.")
  if n <= 0:
    raise ValueError("{n} is not positive, so this calculation is erroneous.")
  result = HomologyMO.zero()
  for partition in get_weighted_partitions(n):
    ll = partition + (n, )
    if multinomial_mod_2(ll) == 1:
      result += HomologyMO.from_partition(to_partition(partition, n))
  return result

# Uses the conjectural observation that Q_1h[\RP^n] = h[M_{2n + 1}] (where M is the relevant Dold manifold) to compute the Hurewicz image of the Dold manifold 
# NOTE This observation is definitely dodgy.
def hurewicz_image_dold_manifold_dyer_lashof(n):
  if n % 2 == 0:
    raise ValueError("{n} is not odd.")
  k = (n - 1)//2
  if k < 2:
    raise ValueError(f"{2} = (n - 1)/2 < 2 is insufficiently large to make this calculation.")
  hurewicz_image_RPn = hurewicz_image_projective_space_lagrange_inversion(k)
  return dyer_lashof(hurewicz_image_RPn, 1)

def from_total_sw_class_to_sw_number_polynomial(total_sw_class: TruncPoly):
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

def hurewicz_image_via_BO_coaction(tangential_sw_class: TruncPoly):
  normal_sw_class = tangential_sw_class.inverse()
  normal_sw_poly = from_total_sw_class_to_sw_number_polynomial(normal_sw_class)
  return change_basis_w_to_e(normal_sw_poly)

def hurewicz_image_projective_space_via_BO_coaction(n):
  if n % 2 != 0:
    raise ValueError("Dimension {n} is odd and therefore invalid.")
  tangential_sw_class = tangential_sw_classes_projective_space(n)
  return hurewicz_image_via_BO_coaction(tangential_sw_class)

# Computes the normal Stiefel--Whitney number polynomial of the relevant dold manifold and converts to the e_i basis.
def hurewicz_image_dold_manifold_via_BO_coaction(n):
  if n % 2 == 0:
    raise ValueError("Dimension {n} is even and therefore invalid.")
  r, s = get_dold_manifold_indices_for_degree(n)
  tangential_sw_class = tangential_sw_classes_dold_manifold(r, s) 
  return hurewicz_image_via_BO_coaction(tangential_sw_class)

def hurewicz_image_even(n):
  if USE_LAGRANGE_INVERSION_FOR_HUREWICZ_IMAGE_OF_PROJECTIVE_SPACES:
    return hurewicz_image_projective_space_lagrange_inversion(n)
  return hurewicz_image_projective_space_via_BO_coaction(n)

def hurewicz_image_odd(n):
  index = is_one_less_than_power_of_two(n)
  if index is None:
    return hurewicz_image_dold_manifold_via_BO_coaction(n)
  return dyer_lashof(HomologyMO.e(1), tuple([1 for i in range(index - 1)]))

def hurewicz_image(n):
  return hurewicz_image_even(n) if (n % 2 == 0) else hurewicz_image_odd(n)