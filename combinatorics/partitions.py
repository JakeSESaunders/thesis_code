from collections import Counter
from config import DEBUG
from itertools import permutations

# TODO update this file to use functools cache
PARTITIONS = {}
ENUMERATED = {}

def enumerate_partition(partition):
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
def load_ordered_partitions(n):
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

def get_ordered_partitions(n):
  global PARTITIONS
  if n not in PARTITIONS:
    load_ordered_partitions(n)
  return PARTITIONS[n]

def get_unordered_partitions(n):
  ordered_partitions = get_ordered_partitions(n)
  unordered_partitions = []
  for partition in ordered_partitions:
    perms = permutations(partition)
    for perm in perms:
      tuple_perm = tuple(perm)
      if tuple_perm not in unordered_partitions:
        unordered_partitions.append(tuple_perm)
  return unordered_partitions

def to_weighted_partition(partition, n):
  weights = Counter(partition)
  weighted_partition = []
  for i in range(1, n + 1):
    if i in weights:
      weighted_partition.append(weights[i])
    else:
      weighted_partition.append(0)
  return tuple(weighted_partition)

def to_partition(weighted_partition, n):
  partition = []
  for i in range(n):
    partition += [i + 1 for j in range(weighted_partition[i])]
  return tuple(partition)

def get_dual(partition):
  dual_partition = []
  for n in range(1, max(partition) + 1):
    count_greater_than_or_equal_to_n = 0
    for i in partition:
      if i >= n:
        count_greater_than_or_equal_to_n += 1
    dual_partition.append(count_greater_than_or_equal_to_n)

  dual_partition.sort()
  return tuple(dual_partition)

# Returns a list of lists [a_1, \dots, a_k] such that the weighted sum 1a_1 + 2a_2 + ... + ka_k = n NOTE this has 0 indexing
def get_weighted_partitions(n):
  if n == 0:
    return [()]
  weighted_partitions = []
  
  for partition in get_ordered_partitions(n):
    weighted_partitions.append(to_weighted_partition(partition, n))
  return weighted_partitions

def total(partition):
  return sum(partition)

def weighted_total(weighted_partition):
  total = 0
  for i in range(len(weighted_partition)):
    total += (i + 1) * weighted_partition[i]

  return total