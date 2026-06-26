from functools import cache
from homomorphisms.image import image
import numpy as np
from sympy import Matrix

# Let R be a poly ring with generators x_i, and S some vector space.
# Let f be a function from int to S with f(i) encoding the image of the generator R.
# Return the matrix encoding the function f. 
@cache
def to_matrix(f, R, S, degree):
  R_basis = R.get_basis(degree)
  dim_R = len(R_basis)
  S_basis = S.get_basis(degree)
  dim_S = len(S_basis)

  matrix = np.zeros((dim_R, dim_S))
  
  for i in range(dim_R):
    r = R_basis[i]
    f_r = image(f, r, S)
    for j in range(dim_S):
      s = S_basis[j].summands[0] # NOTE for now let's say by convention the basis is just products of generators.
      if s in f_r.summands:
        matrix[i][j] = 1
  return matrix

# Suppose given a matrix representing any linear map from R -> S
# Return the corresponding function, defined on arbitrary polynomials.
def from_matrix(M, R, S, degree):
  def f(r):
    if len(r.summands) == 0:
      return S.zero()
    if not len(r.summands) == 1:
      return sum([f(x) for x in r.get_summands_as_polys], S.zero())

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

# Suppose f is a homomorphism from R to S
# Apply the dual, sending s* to the sum of r* s.t. <s*, f(r)> = 1
@cache
def dual(f, s, R):
  S = s.__class__
  if len(s.summands) == 0:
    return R.zero()
  if not len(s.summands) == 1:
    return sum([dual(f, x, R) for x in s.get_summands_as_polys()], R.zero())
  
  degree = s.max_degree()
  M = to_matrix(f, R, S, degree).transpose()
  return from_matrix(M, S, R, degree)(s)

@cache
def inverse(f, s, R):
  S = s.__class__
  if len(s.summands) == 0:
    return R.zero()
  if not len(s.summands) == 1:
    return sum([dual(f, x, R) for x in s.get_summands_as_polys()], R.zero())

  degree = s.max_degree()
  M = Matrix(to_matrix(f, R, S, degree))
  M_inv = M.inv_mod(2)
  return from_matrix(M, S, R, degree)(s)