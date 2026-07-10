from rings.polynomial import PolyRing
from functools import cache

# 
@cache
def extend_multiplicatively(f, r: PolyRing, S):
  """Let R, S be polynomial rings (S can also be a Tensor). Let f be a function that takes a generator of R and returns an element s of S.
  image returns the extension of f to the entire ring R, not just the generators."""
  no_summands = len(r.summands)
  if no_summands == 0:
    return S.zero()
  if no_summands > 1:
    summands_as_polys = r.get_summands_as_polys()
    return sum([extend_multiplicatively(f, summand, S) for summand in summands_as_polys], S.zero())
  
  # We have a single summand only
  summand = r.summands[0]
  result = S.one()
  for factor in summand:
    result *= f(factor)

  return result