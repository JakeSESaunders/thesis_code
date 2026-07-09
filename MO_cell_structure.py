from functools import cache
from math import floor, log2
from combinatorics.indices import is_power_of_two
from dyer_lashof.formal import DL

@cache
def get_classes_in_degree(n):
  result = []
  k = is_power_of_two(n)
  if k is not None:
    result.append(DL([], n, upper=False))
  for i in range(1, floor(n/2) + 1):
    difference = n - 2*i
    if difference == 0:
      continue
    for gen in get_classes_in_degree(i):
      if len(gen.operation) != 0 and gen.operation[0] < difference:
        continue
      new_gen = gen.Q(difference)
      result.append(new_gen)

  return result

# Returns the classes of the form Q_{2^i_1} \cdots Q_{2^i_n}e_{2^{i_{n + 1}}} with i_1, \dots, i_{n + 1} a non-decreasing seq of positive integers
@cache
def get_good_classes_in_degree(n):
  result = []
  k = is_power_of_two(n)
  if k is not None:
    result.append(DL([], n))
  for i in range(1, floor(n/2) + 1):
    difference = n - 2*i
    if difference == 0:
      continue
    if is_power_of_two(difference) is None:
      continue
    for gen in get_good_classes_in_degree(i):
      if len(gen.operation) != 0 and gen.operation[0] < difference:
        continue
      new_gen = gen.Q(difference)
      if new_gen.generator < difference:
        continue
      result.append(new_gen)

  return result

@cache
def get_classes_to_kill_in_degree(n):
  result = []
  # get the classes of the form Q_ne_{2^i} where n > 2^i.
  for i in range(0, floor(log2(n))):
    difference = n - 2*2**i
    if is_power_of_two(difference) is not None and difference > 2**i:
      result.append(
        DL(
          [difference], 2**i
        )
      )
  # now get the classes which are prefaced by a non-power of 2.
  for i in range(1, floor(n/2) + 1):
    difference = n - 2*i
    if difference == 0:
      continue
    for gen in get_good_classes_in_degree(i):
      if len(gen.operation) != 0 and gen.operation[0] < difference:
        continue # sequence has to be non-decreasing
      if is_power_of_two(difference) is not None:
        continue
      new_gen = gen.Q(difference)
      result.append(new_gen)
  return result

@cache
def get_classes_already_killed_in_degree(n):
  result = []
  for i in range(1, floor(n/2) + 1):
    difference = n - 2*i
    if difference == 0:
      continue
    for gen in get_classes_to_kill_in_degree(i):
      if len(gen.operation) != 0 and gen.operation[0] < difference:
        continue
      new_gen = gen.Q(difference)
      result.append(new_gen)

  return result