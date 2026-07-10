from rings.bo import HomologyMO
from combinatorics.partitions import get_weighted_partitions, to_partition
from combinatorics.multinomial import multinomial_mod_2, binomial_mod_2
from dyer_lashof.monomial import DyerLashofMonomial
from dyer_lashof.extend import extend_dl

def dyer_lashof_e_0_inverse(n):
  """Implements the calculation of Q^ne_0^{-1} due to Priddy, Corollary 2.3.

  Uses upper indexing, though this agrees with lower indexing as the degree of e_0^{-1} is 0."""
  if n < 0:
    raise ValueError("Dyer--Lashof operations with negative index are undefined.")

  if n == 0:
    return HomologyMO.one()

  Qne0_inverse = HomologyMO.zero()
  for ll in get_weighted_partitions(n):
    if multinomial_mod_2(ll) != 0:
      Qne0_inverse += HomologyMO.from_partition(to_partition(ll, n))
  return Qne0_inverse

def dyer_lashof_e_k_one_component(n, k):
  """Implements the calculation of Q^ne_k due to Priddy, Theorem 1.1.

  Uses upper indexing."""
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

def dyer_lashof_single_indeterminate(n: int, upper: bool, k: int) -> HomologyMO:
  """Uses the Cartan formula to compute Q^ne_0^{-1}e_k = \sum_{i + j = k}Q^ie_0^{-1}Q^je_n.

  NOTE: Usually we use n as a generator and k as an operation, but here we use Priddy's notation for ease of comparison."""
  if not upper:
    n = n + k
  result = HomologyMO.zero()
  for i in range(0, n + 1):
    j = n - i
    result += dyer_lashof_e_0_inverse(i) * dyer_lashof_e_k_one_component(j, k)
  return result

def dyer_lashof_HMO(QI: DyerLashofMonomial, x: HomologyMO) -> HomologyMO:
  """Compute the result of QIx in H_*MO."""
  return extend_dl(dyer_lashof_single_indeterminate, QI, x)