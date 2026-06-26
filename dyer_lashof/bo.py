from rings.bo import HomologyMO
from combinatorics.partitions import get_weighted_partitions, to_partition
from combinatorics.multinomial import multinomial_mod_2, binomial_mod_2

# Implements the calculation of Priddy, Corollary 2.3.
# Uses upper indexing, though this agrees with lower indexing as the degree of e_0^{-1} is 0.
def dyer_lashof_e_0_inverse(n):
  if n < 0:
    raise ValueError("Dyer--Lashof operations with negative index are undefined.")

  if n == 0:
    return HomologyMO.one()

  Qne0_inverse = HomologyMO.zero()
  for ll in get_weighted_partitions(n):
    if multinomial_mod_2(ll) != 0:
      Qne0_inverse += HomologyMO.from_partition(to_partition(ll, n))
  return Qne0_inverse

# Implements the calculation of Priddy, Theorem 1.1
# Uses upper indexing.
def dyer_lashof_e_k_one_component(n, k):
  if n < k:
    return HomologyMO.zero()
  if n == k:
    return HomologyMO.e(k)**2
  
  Qnek = HomologyMO.zero()
  for u in range(0, k + 1):
    if binomial_mod_2(n - k + u - 1, u) % 2 == 0:
      continue
    summand = []
    if n + u != 0:
      summand.append(n + u)
    if k - u != 0:
      summand.append(k - u)
    summand.sort()
    Qnek += HomologyMO.from_partition(tuple(summand))
  return Qnek

# Uses the Cartan formula to compute Q^ke_0^{-1}e_n = \sum_{i + j = n}Q^ie_0^{-1}Q^je_n
def dyer_lashof_single_indeterminate(n, k, upper=True):
  if not upper:
    k = k + n
  result = HomologyMO.zero()
  for i in range(0, n + 1):
    j = n - i
    dl_op = dyer_lashof_e_0_inverse(i)
    dl_op_2 = dyer_lashof_e_k_one_component(j, k)
    result += dyer_lashof_e_0_inverse(i) * dyer_lashof_e_k_one_component(j, k)
  return result

# Uses the Cartan formula to compute Q^{monomial}P(e_1, e_2, \dots) for a polynomial P in H_*(BO).
# Can use upper or lower indexing.
def dyer_lashof(poly: HomologyMO, monomial, upper=False):
  # Correct for if the user accidentally provided an integer instead of a list/tuple
  if isinstance(monomial, int):
    monomial = [monomial]

  if not upper:
    # convert the monomial to upper indexing notation
    # only do this if the class is homogeneous, else error.
    if not poly.is_homogeneous():
      raise NotImplementedError("Not yet implemented lower indexing notation for Dyer--Lashof operations on non-homogeneous classes.")
    degree = poly.max_degree()

    new_monomial = []
    for i in range(len(monomial)):
      new_monomial.insert(0, sum(new_monomial) + degree + monomial[-(i + 1)])
    return dyer_lashof(poly, tuple(new_monomial), upper=True)
  
  # Use additivity to reduce to the calculation of the DL op on each summand.
  if len(poly.summands) == 0:
    return poly
  if len(poly.summands) > 1:
    dl_op_summands = [dyer_lashof(summand, monomial, upper) for summand in poly.get_summands_as_polys()]
    return sum(dl_op_summands, HomologyMO.zero())
  
  # Now we have a single summand, and we reduce to the case that the monomial has length 1.
  if len(monomial) == 0:
    return poly
  if len(monomial) > 1:
    # Given monomial Q_{i_0}...Q_{i_n}, first compute Q_{i_1}...Q_{i_n}, then apply Q_{i_0} to it.
    return dyer_lashof(dyer_lashof(poly, monomial[1:], upper=True), (monomial[0], ), upper=True)

  # Now the length of the monomial is 1,
  operation = monomial[0]
  # and we have only 1 summand.
  summand = poly.summands[0]
  
  # If the summand has length 0, then the class is given by 1, which has vanishing operation Q^i unless i = 0.
  if len(summand) == 0:
    if operation == 0:
      return HomologyMO.one()
    else:
      return HomologyMO.zero()
  
  # We use the Cartan formula to reduce to the case that the summand has only 1 factor.
  if len(summand) > 1:
    first_factor = summand[0]
    partial_summand = HomologyMO.from_partition(tuple(summand[1:]))
    result_summands = []
    for i in range(0, operation + 1):
      operation_first_factor = dyer_lashof_single_indeterminate(i, first_factor, upper=True)
      operation_partial_summand = dyer_lashof(partial_summand, (operation - i,), upper=True)
      result_summands.append(operation_first_factor*operation_partial_summand)
    return sum(result_summands, HomologyMO.zero())
  
  # We now have a single operation and a single summand with a single factor.
  # This is what dyer_lashof_single_indeterminate can handle.
  factor = summand[0]
  return dyer_lashof_single_indeterminate(operation, factor)