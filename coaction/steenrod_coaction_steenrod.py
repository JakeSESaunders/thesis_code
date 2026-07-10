from rings.dual_steenrod import DualSteenrod, DualSteenrodTensorSquare
from homomorphisms.image import extend_multiplicatively

# Implements the coaction of the dual Steenrod algebra on itself.
# Formula appears in Milnor - The Steenrod algebra and its dual, Theorem 3.

def xi(k):
  return DualSteenrod.xi(k)

def A_coaction_A_on_generators(k):
  """The coaction of the dual Steenrod algebra on itself, on generators."""
  summands = []
  for i in range(0, k + 1):
    left = xi(k - i)**(2**i)
    right = xi(i)
    summands.append(
      DualSteenrodTensorSquare.from_factors(left, right)
    )
  return sum(summands, DualSteenrodTensorSquare.zero())

def A_coaction_A(r: DualSteenrod):
  """The coaction of the dual Steenrod algebra on itself."""
  return extend_multiplicatively(A_coaction_A_on_generators, r, DualSteenrodTensorSquare)