from collections import Counter
from config import DEBUG
from itertools import permutations

# TODO change terminology in this file to use order_depedent/order_indepednent.

# TODO update this file to use functools cache
PARTITIONS = {}
ENUMERATED = {}

def enumerate_partition(partition: tuple[int]) -> int:
  """Returns the index of this partition in the list of partitions. Useful for ordering bases of PolyRings."""
  global ENUMERATED
  n = sum(partition)
  
  if n not in ENUMERATED:
    if DEBUG:
      print(f"Enumerating partitions of size {n}")
    ENUMERATED[n] = {}
    partitions = get_ordered_partitions(n)
    for i in range(len(partitions)):
      ENUMERATED[n][tuple(partitions[i])] = i
    if DEBUG:
      print(f"Enumerated partitions!")

  return ENUMERATED[n][partition]

# Populate the global dictionary PARTITIONS with the ordered partitions of integers.
# NOTE This can only get into the 30s with any speed.
def load_ordered_partitions(n: int) -> list[tuple[int]]:
  global PARTITIONS
  if n in PARTITIONS:
    return
  
  if n < 0:
    raise ValueError("Cannot get partitions of a negative integer.")
  if n == 0:
    PARTITIONS[0] = [
      ()
    ] # The only partition of no elements is the empty partition
    return
  
  if (n - 1) not in PARTITIONS:
    load_ordered_partitions(n - 1)
  
  result = []
  for i in range(1, n):
    for start in PARTITIONS[i]:
      for end in PARTITIONS[n - i]:
        if start[-1] <= end[0]:
          if start + end not in result:
            result.append(tuple(start + end))
  result.append((n, )) # NOTE this is weird notation for a single element tuple
  # print(f"Loaded partitions of {n}.")
  PARTITIONS[n] = result

def get_ordered_partitions(n: int) -> list[tuple[int]]:
  global PARTITIONS
  if n not in PARTITIONS:
    load_ordered_partitions(n)
  return PARTITIONS[n]

def get_unordered_partitions(n: int) -> list[tuple[int]]:
  """Return a list of unordered partitions of n.
  
  For example, (1, 2) and (2, 1) are distinct as unordered partitions.
  """
  ordered_partitions = get_ordered_partitions(n)
  unordered_partitions = []
  for partition in ordered_partitions:
    perms = permutations(partition)
    for perm in perms:
      tuple_perm = tuple(perm)
      if tuple_perm not in unordered_partitions:
        unordered_partitions.append(tuple_perm)
  return unordered_partitions

def to_weighted_partition(partition: tuple[int], n: int) -> tuple[int]:
  weights = Counter(partition)
  weighted_partition = []
  for i in range(1, n + 1):
    if i in weights:
      weighted_partition.append(weights[i])
    else:
      weighted_partition.append(0)
  return tuple(weighted_partition)

# TODO is this the same as getting the conjugate???
def to_partition(weighted_partition: tuple[int], n: int=-1) -> tuple[int]:
  """Convert the partition (a_1, a_2, \dots, a_n) to a partition containing 1 a_1 times, 2 a_2 times, etc."""
  if n == -1:
    n = weighted_total(weighted_partition)
  partition = []
  for i in range(n):
    partition += [i + 1 for j in range(weighted_partition[i])]
  return tuple(partition)

def get_conjugate(partition: tuple[int]) -> tuple[int]:
  """Return the conjugate partition."""
  dual_partition = []
  for n in range(1, max(partition) + 1):
    count_greater_than_or_equal_to_n = 0
    for i in partition:
      if i >= n:
        count_greater_than_or_equal_to_n += 1
    dual_partition.append(count_greater_than_or_equal_to_n)

  dual_partition.sort()
  return tuple(dual_partition)

def get_weighted_partitions(n: int) -> list[tuple[int]]:
  """
  Returns a list of tuples (a_1, \dots, a_k) such that the weighted sum 1a_1 + 2a_2 + ... + ka_k = n.
  Note that a_1 is the 0th element of the tuple.
  """
  if n == 0:
    return [()]
  weighted_partitions = []
  
  for partition in get_ordered_partitions(n):
    weighted_partitions.append(to_weighted_partition(partition, n))
  return weighted_partitions

def total(partition: tuple[int]) -> int:
  """Return the sum of the parts of partition."""
  return sum(partition)

def weighted_total(weighted_partition: tuple[int]) -> int:
  """Given a partition (a_1, a_2, a_3, \dots, a_n), return 1a_1 + 2a_2 + ... + na_n."""
  total = 0
  for i in range(len(weighted_partition)):
    total += (i + 1) * weighted_partition[i]

  return total