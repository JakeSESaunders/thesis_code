from rings.polynomial import PolyRing
from combinatorics.partitions import get_ordered_partitions
from functools import cache
from collections import Counter

class NatPolyRing(PolyRing):
  """Base class for polynomial rings with at most one generator in each positive integer degree."""

  # NOTE we override __init__ so we don't have to iterate over the summands twice to correct zeros.
  def __init__(self, summands: list[tuple[int]]):
    R = self.__class__

    new_summands = [] # a list of sorted summands
    counter = Counter(summands)
    for summand in counter:
      # check the input for duplicates
      if counter[summand] != 1: raise ValueError(f"List of summands may not contain duplicates: summand {summand} appeared {counter[summand]} times.")
      new_summand = []

      # check for invalid generators
      for factor in summand:
        if factor == 0:
          continue # NOTE If every factor is 0, then we add an empty tuple to the summands, and this represents 1, as expected, since e.g. e_0 = 1.
        if not R.is_valid_generator(factor): raise ValueError(f"Summand {summand} contains invalid generator {factor}.")
        new_summand.append(factor)

      new_summands.append(tuple(sorted(new_summand)))
    self.summands = sorted(new_summands)

  @classmethod
  def is_valid_generator(cls, i) -> bool:
    if not isinstance(i, int):
      return False
    if i < 0:
      return False
    return True # NOTE we deal with generator with index 0 appropriately

  @classmethod
  def degree(cls, x) -> int:
    if not cls.is_valid_generator(x): raise ValueError("Cannot get degree of non-integer value {x}.")
    return x

  @classmethod
  @cache
  def get_basis(cls, degree: int) -> list:
    basis = []
    for partition in get_ordered_partitions(degree):
      partition_is_valid = True
      for part in partition:
        if not cls.is_valid_generator(part):
          partition_is_valid = False
          break
      if partition_is_valid:
        basis.append(
          cls.from_partition(partition)
        )
    return basis

  @classmethod
  def from_partition(cls, partition: tuple[int]):
    """Given a partition (i.e. a tuple of ints) return the corresponding polynomial."""
    return cls([partition])
