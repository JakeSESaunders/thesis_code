from rings.truncated_poly import TruncPoly, TruncPolySignature
from combinatorics.multinomial import binomial_mod_2
from config import DEBUG

def signature_cohomology_projective_space(d):
  return TruncPolySignature(
    [1],
    [d + 1],
    ["a"]
  )

def signature_cohomology_dold_manifold(m, n):
  return TruncPolySignature(
    [1, 2],
    [m + 1, n + 1],
    ["c", "d"]
  )

def tangential_sw_classes_projective_space(d):
  if DEBUG:
    print("Computing tangential SW classes for RP^{d}...")
  signature = signature_cohomology_projective_space(d)
  result = [
      tuple([0 for i in range(j)]) for j in range(d + 1) if binomial_mod_2(d + 1, j) == 1
    ]
  if DEBUG:
    print("Computed!")
  return TruncPoly(
    signature,
    result
  )

# Here, generator 0 is c and generator 1 is d.
def tangential_sw_classes_dold_manifold(m, n):
  signature = signature_cohomology_dold_manifold(m, n)
  left_factor = TruncPoly(
    signature,
    [(), (0, ), (1, )]
  )**(n + 1)
  right_factor = TruncPoly(
    signature,
    [(), (0, )]
  )**m
  return left_factor * right_factor