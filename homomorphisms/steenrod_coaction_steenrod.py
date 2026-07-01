from rings.dual_steenrod import DualSteenrod, DualSteenrodTensorSquare

# Implements the coaction of the dual Steenrod algebra on itself.
# Formula appears in Milnor - The Steenrod algebra and its dual, Theorem 3.

def xi(k):
  return DualSteenrod.xi(k)

def coaction_xi(k):
  summands = []
  for i in range(0, k + 1):
    left = xi(k - i)**(2**i)
    right = xi(i)
    summands.append(
      DualSteenrodTensorSquare.from_factors(left, right)
    )
  return sum(summands, DualSteenrodTensorSquare.zero())