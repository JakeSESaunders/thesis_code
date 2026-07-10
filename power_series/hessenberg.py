# NOTE this function assumes everything is mod 2 and doesn't check that the matrix is square.
def cofactor(M: list[list], i: int) -> list[list]:
  """Given an nxn matrix M and 0 <= i < n, return the cofactor expanded along the top row."""
  n = len(M)
  C = [
    M[j][:i] + M[j][(i + 1):] for j in range(1, n)
  ]
  return C

# NOTE
# Don't bother with signs --- assume everything is mod 2.
# Check for zeros before computing.
# Useful for Hessenberg matrices as most of these entries will be zero.
# Doesn't check if the matrix is square.
def determinant(M: list[list]):
  """Given an length n list of length n lists, compute the determinant recursively by Laplace expansion along the top row."""
  n = len(M)
  if n == 1:
    return M[0][0]

  top_row = M[0]
  cls = top_row[0].__class__

  result = cls.zero()

  for i in range(n):
    if top_row[i] == cls.zero():
      continue
    result += top_row[i] * determinant(cofactor(M, i))

  return result

# given a list L = [a_1, a_2, \dots, a_n], return the Hessenberg matrix as below.
"""
[
  a_1   1     0   ... 0
  a_2   a_1   1   ... 0
  a_3   a_2   a_1 ... 0
  .     .             .
  .           .       .
  .               .   1
  a_n   ...   ... ... a_1
]
"""
def hessenberg(L: list):
  """Return the Hessenberg matrix corresponding to the list L = [a_1, a_2, ..., a_n]."""
  R = L[0].__class__
  n = len(L)
  if n == 1:
    return [[L[0]]]
  smaller_hessenberg = hessenberg(L[:-1])
  big_hessenberg = []
  for i in range(n - 1):
    trailing_entry = R.zero()
    if i == n - 2:
      trailing_entry = R.one()
    big_hessenberg.append(
      smaller_hessenberg[i] + [trailing_entry]
    )
  big_hessenberg.append(
    [x for x in reversed(L)]
  )
  return big_hessenberg

def hessenberg_det(L: list):
  """Uses Hessenberg determinants to compute the coefficients of the multiplicative inverse of a power series
  See e.g. Inselberg - On determinants of Toeplitz--Hessenberg matrices arising in power series."""
  H = hessenberg(L)
  return determinant(H)