from functools import cache
from homomorphisms.matrix import dual, from_matrix
from coaction.w_basis_coproduct import coproduct
from rings.bo import CohomologyBO, HomologyMO, TensorCohomologyBO
import numpy as np
from sympy import Matrix
from config import SEE_CHANGE_OF_BASIS_MATRIX_STATUS_FOR_E_TO_W_BASIS, SEE_CHANGE_OF_BASIS_MATRIX_STATUS_FOR_E_TO_W_BASIS_TELL_AFTER_EACH_ROW

def w(i):
  return CohomologyBO.w(i)

@cache
def multiply_two_monomials(left: tuple[int], right: tuple[int]) -> CohomologyBO:
  """Given two monomials left and right, return the sum of the monomials in H^*BO that include left \otimes right as a summand."""
  f = lambda i: coproduct(w(i))
  s = TensorCohomologyBO.from_factors(left, right)
  return dual(f, s, CohomologyBO)

@cache
def multiply_two_polynomials(left: CohomologyBO, right: CohomologyBO) -> CohomologyBO:
  """Given two polynomials left and right, return the sum of the monomials in H^*BO that include left \otimes right as a summand."""
  result = CohomologyBO.zero()
  for summand in left.get_summands_as_polys():
    for other_summand in right.get_summands_as_polys():
      product_of_summands = multiply_two_monomials(summand, other_summand)
      result += product_of_summands

  return result

@cache
def change_of_basis_matrix_e_to_w(degree):
  """Get the numpy array representing the duality linear map from the e-basis of H_*MO to the w-basis of H^*BO."""
  homology_basis = HomologyMO.get_basis(degree)
  homology_dim = len(homology_basis)
  cohomology_basis = CohomologyBO.get_basis(degree)
  cohomology_dim = len(cohomology_basis)

  matrix = np.zeros((homology_dim, cohomology_dim))

  if SEE_CHANGE_OF_BASIS_MATRIX_STATUS_FOR_E_TO_W_BASIS:
    print(f"Computing change of basis matrix for coproduct in degree {degree}.")
    print(f"There are {homology_dim} coproducts to be computed.")
  for i in range(homology_dim):
    p = homology_basis[i]
    # p is a basis element of H_*MO, so a monomial in the e_i.
    # We need to figure out where the duality sends it.
    # We know e_i = (w_1^i)*, so we need to determine the product of the (w_1^i)^*
    partition = p.summands[0]

    result = CohomologyBO.one()
    for factor in partition:
      result = multiply_two_polynomials(result, CohomologyBO.from_partition([1 for i in range(factor)]))

    for result_summand in result.get_summands_as_polys():
      j = cohomology_basis.index(result_summand)
      matrix[i][j] = 1
    if SEE_CHANGE_OF_BASIS_MATRIX_STATUS_FOR_E_TO_W_BASIS_TELL_AFTER_EACH_ROW:
      print(f"Completed row {i} of matrix...")

  return matrix

@cache
def change_basis_e_to_w(r: HomologyMO) -> CohomologyBO:
  """Given r \in H_*MO, return the linear dual r* \in H^*BO."""
  R = HomologyMO
  S = CohomologyBO
  if len(r.summands) == 0:
    return S.zero()
  if not len(r.summands) == 1:
    return sum([change_basis_e_to_w(x) for x in r.get_summands_as_polys()], S.zero())

  degree = r.max_degree()
  M = change_of_basis_matrix_e_to_w(degree)
  return from_matrix(M, R, S, degree)(r)

@cache
def change_basis_w_to_e(s: CohomologyBO) -> HomologyMO:
  """Given r \in H^*BO, return the linear dual r* \in H_*MO."""
  R = HomologyMO
  S = CohomologyBO
  if len(s.summands) == 0:
    return R.zero()
  if not len(s.summands) == 1:
    return sum([change_basis_w_to_e(x) for x in s.get_summands_as_polys()], R.zero())

  degree = s.max_degree()
  M = Matrix(change_of_basis_matrix_e_to_w(degree)).inv_mod(2)
  M = np.array(M)
  return from_matrix(M, S, R, degree)(s)