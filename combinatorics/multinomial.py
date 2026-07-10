from math import factorial, isclose

# Given a list of coefficients ll, return the multinomial (generalised binomial)
# coefficient given by sum(ll)!/(ll[0]! x ... x ll[-1]!).
# NOTE responds badly to negative inputs.
def multinomial(ll: list[int]) -> int:
  numerator = factorial(sum(ll))
  denominator = 1
  for l in ll:
    denominator *= factorial(l)
  result = numerator / denominator
  error = result - int(result)

  if not isclose(error, 0):
    raise RuntimeError(f"Multinomial coefficient for partition {ll} is {result}, which is not an integer.")
  return int(numerator / denominator)

# Uses a generalisation of Lucas' theorem to compute multinomial coefficients mod 2 by xoring binary digits.
# the xor and the sum are equal iff the numbers have disjoint binary digits iff the multinomial coeff is odd.
def multinomial_mod_2(ll: list[int]) -> int:
  xored = 0
  for l in ll:
    xored = xored ^ l
  return 1 if xored == sum(ll) else 0

def binomial_mod_2(n: int, k: int) -> int:
  return multinomial_mod_2([k, n - k])
