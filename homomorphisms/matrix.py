from functools import cache
from homomorphisms.image import extend_multiplicatively
import numpy as np
from sympy import Matrix
from collections.abc import Callable

# NOTE these functions take PolyRings and Tensors as arguments, so 
# NOTE the type of the Matrix argument is unclear.

@cache
def to_matrix(f, R: type, S: type, degree: int):
  """Let R and S be rings. If f is a function encoding the image of a linear map R -> S on the generators of R, then extend f to a monomials basis of R in degree d, and obtain the corresponding matrix in the canonical bases of R and S.
  
  Keyword arguments:
  f: a function encoding the image of a linear map R -> S on the generators of R.
  R: a type implementing get_basis, addition and multiplication (a PolyRing or Tensor).
  S: a type implementing get_basis, addition and multiplication (a PolyRing or Tensor).
  degree: the total degree of monomials to use as bases for R and S. 
  """
  R_basis = R.get_basis(degree)
  dim_R = len(R_basis)
  S_basis = S.get_basis(degree)
  dim_S = len(S_basis)

  matrix = np.zeros((dim_R, dim_S))
  
  for i in range(dim_R):
    r = R_basis[i]
    f_r = extend_multiplicatively(f, r, S)
    for j in range(dim_S):
      s = S_basis[j].summands[0]
      if s in f_r.summands:
        matrix[i][j] = 1
  return matrix

def from_matrix(M, R: type, S: type, degree: int) -> Callable:
  """Suppose given a matrix representing any linear map from R -> S
  Return the corresponding function, defined on arbitrary polynomials.
  
  Keyword arguments:
  M: a matrix. Must have dimensions given by the sizes of the bases of R and S in degree d.
  f: a function encoding the image of a linear map R -> S on the generators of R.
  R: a type implementing get_basis, addition, multiplication and get_summands_as_polys (a PolyRing or Tensor).
  S: a type implementing get_basis, addition and multiplication (a PolyRing or Tensor).
  degree: the total degree of monomials to use as bases for R and S.
  """
  def f(r):
    if len(r.summands) == 0:
      return S.zero()
    if not len(r.summands) == 1:
      return sum([f(x) for x in r.get_summands_as_polys()], S.zero())

    # Now we know s has a single summand, and it's in the basis by assumption.
    r_degree = r.max_degree()
    if r_degree != degree:
      raise ValueError(f"Degree of specified element {r_degree} is not expected degree {degree}.")

    R_basis = R.get_basis(degree)
    i = R_basis.index(r)
    S_basis = S.get_basis(degree)
    dim_S = len(S_basis)
    
    result = S.zero()
    for j in range(dim_S):
      if M[i][j] == 1:
        result += S_basis[j]

    return result

  return f

@cache
def dual(f: Callable, s, R: type):
  """Suppose f is a homomorphism from R to S
  Apply the dual, sending s* to the sum of those r* s.t. <s*, f(r)> = 1"""
  S = s.__class__
  if len(s.summands) == 0:
    return R.zero()
  if not len(s.summands) == 1:
    return sum([dual(f, x, R) for x in s.get_summands_as_polys()], R.zero())
  
  degree = s.max_degree()
  M = to_matrix(f, R, S, degree).transpose()
  return from_matrix(M, S, R, degree)(s)

@cache
def inverse(f: Callable, s, R: type):
  """Given an isomorphism f, use matrix inversion to get f^{-1}(s) \in R. Errors badly if f isn't an isomorphism."""
  S = s.__class__
  if len(s.summands) == 0:
    return R.zero()
  if not len(s.summands) == 1:
    return sum([inverse(f, x, R) for x in s.get_summands_as_polys()], R.zero())

  degree = s.max_degree()
  M = Matrix(to_matrix(f, R, S, degree))
  M_inv = M.inv_mod(2)
  M_inv = np.array(M_inv)
  return from_matrix(M_inv, S, R, degree)(s)