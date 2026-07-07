# Let x be an element of a PolyRing R.
# Let dl_ops_for_generators be a class that encodes the DL ops Q_ix_r where x_r is a generator.
# Use the Cartan formula to extend this to the rest of the poly ring.
def extend_dl(dl_ops_for_generators, I, poly, upper=False):
  # correct for the user accidentally not wrapping a single element in a tuple
  if (not isinstance(I, tuple)) and (not isinstance(I, list)):
    I = (I, ) # NOTE does not correct for invalid inputs.

  if not upper:
    # convert the monomial to upper indexing notation
    # only do this if the class is homogeneous, else error.
    if not poly.is_homogeneous():
      raise NotImplementedError("Not yet implemented lower indexing notation for Dyer--Lashof operations on non-homogeneous classes.")
    degree = poly.max_degree()

    new_I = []
    for i in range(len(I)):
      new_I.insert(0, sum(new_I) + degree + I[-(i + 1)])
    return extend_dl(dl_ops_for_generators, tuple(new_I), poly, True)

  extend = lambda _I, _poly: extend_dl(dl_ops_for_generators, _I, _poly, upper)
  R = poly.__class__

  # Use additivity to reduce to the calculation of the DL op on each summand.
  if len(poly.summands) == 0:
    return poly
  if len(poly.summands) > 1:
    dl_op_summands = [extend(I, summand) for summand in poly.get_summands_as_polys()]
    return sum(dl_op_summands, R.zero())

  # Now we have a single summand, and we reduce to the case that the monomial has length 1.
  if len(I) == 0:
    return poly
  if len(I) > 1:
    Qi_poly = extend(I[-1], poly)
    return extend(I[:-1], Qi_poly)

  # We now have only one operation
  k = I[0]
  # and only one summand
  summand = poly.summands[0]

  # If the summand has length 0, then the class is given by 1, which has vanishing operation Q^i unless i = 0.
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
      operation_first_factor = dl_ops_for_generators(j, first_factor, upper)
      operation_remaining_factors = extend(k - j, remaining_factors)
      result_summands.append(operation_first_factor*operation_remaining_factors)
    return sum(result_summands, R.zero())
  
  # We now have a single operation and a single summand with a single factor.
  # This is what dyer_lashof_single_indeterminate can handle.
  factor = summand[0]
  return dl_ops_for_generators(k, factor, upper)