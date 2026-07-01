# If n = 2^i - 1, return n. Otherwise, return None.
def is_one_less_than_power_of_two(n):
  return is_power_of_two(n + 1)

def is_power_of_two(n):
  i = 0
  while 2**i <= n:
    if 2**i == n:
      return i
    i += 1
  return None

def get_dold_manifold_indices_for_degree(n):
  r = 1
  while 2**r < n + 1:
    if (n + 1) % (2**r) == 0:
      r += 1
    else:
      r -= 1
      s = (((n + 1) // 2**r) - 1) // 2
      return (2**r - 1, 2**r * s)
  raise ValueError("{n} is either even or one less than a power of 2.")