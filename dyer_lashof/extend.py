from dyer_lashof.monomial import DyerLashofMonomial
from collections.abc import Callable
from rings.polynomial import PolyRing

def extend_dl(dl_ops_for_generators: Callable[[int, bool, object], object], QI: DyerLashofMonomial, poly: PolyRing) -> PolyRing:
  """Apply the Dyer--Lashof operation Q_I to polynomial poly in ring R.
  
  Keyword arguments:
  dl_ops_for_generators: a function sending (k, upper, x) where x is a generator of R to Q^kx (if upper) or Q_kx (if not upper).
  """

  # lower indexing notation only makes sense for homogeneous classes poly, so check this
  if QI.lower and not poly.is_homogeneous():
    raise ValueError(f"Cannot apply lower-indexed operation {QI} to non-homogeneous polynomial {poly}.")

  R = poly.__class__

  # Use additivity to reduce to the calculation of the DL op on each summand.
  if len(poly.summands) == 0:
    return poly
  if len(poly.summands) > 1:
    dl_op_summands = [extend_dl(dl_ops_for_generators, QI, summand) for summand in poly.get_summands_as_polys()]
    return sum(dl_op_summands, R.zero())

  # Now we have a single summand, and we reduce to the case that the monomial has length 1.
  if QI.length() == 0:
    return poly
  if QI.length() > 1:
    Qk, QI_shorter = QI.pop()
    QI_shorter_poly = extend_dl(dl_ops_for_generators, QI_shorter, poly)
    return extend_dl(dl_ops_for_generators, Qk, QI_shorter_poly)

  # We now have only one operation
  k = QI.operation[0]
  # and only one summand
  summand = poly.summands[0]

  # If the summand has length 0, then the class is given by 1 iff k == 0, and 0 otherwise, regardless of upper/lower indexing notation.
  if len(summand) == 0:
    if k == 0:
      return R.one()
    else:
      return R.zero()

  # We use the Cartan formula to reduce to the case that the summand has only 1 factor.
  if len(summand) > 1:
    result_summands = []
    first_factor = summand[0]
    remaining_factors = R([summand[1:]])
    
    for j in range(0, k + 1):
      operation_first_factor = dl_ops_for_generators(j, QI.upper, first_factor)
      operation_remaining_factors = extend_dl(dl_ops_for_generators, DyerLashofMonomial(k - j, QI.upper), remaining_factors)
      result_summands.append(operation_first_factor*operation_remaining_factors)
    return sum(result_summands, R.zero())
  
  # We now have a single operation and a single summand with a single factor.
  # This is what dl_ops_for_generators can handle.
  factor = summand[0]
  return dl_ops_for_generators(k, QI.upper, factor)