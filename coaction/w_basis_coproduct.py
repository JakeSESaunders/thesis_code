from rings.bo import CohomologyBO, TensorCohomologyBO
from homomorphisms.image import extend_multiplicatively
from functools import cache

def w(i):
  return CohomologyBO.w(i)

@cache
def coproduct_BO_on_generators(n):
  """The comultiplication of H^*BO, on generators."""
  if n == 0:
    return TensorCohomologyBO.one()
  coproduct_of_generator_summands = []
  for i in range(n + 1):
    coproduct_of_generator_summands.append(
      TensorCohomologyBO.from_factors(w(i), w(n - i)) # NOTE w(0) automatically evaluates to 1
    )
  result = sum(coproduct_of_generator_summands, TensorCohomologyBO.zero())
  return result

# Returns the coproduct as a Tensor.
@cache
def coproduct(w_poly: CohomologyBO):
  """The comultiplication of H^*BO."""
  return extend_multiplicatively(coproduct_BO_on_generators, w_poly, TensorCohomologyBO)