from math import floor, ceil
from combinatorics.multinomial import binomial_mod_2

# NOTE This file uses lower indexing.
# TODO could modify to make this agnostic about the indexing.

def is_admissible(sequence: tuple[int]) -> int | None:
    """Return the first index at which a sequence is not admissible. Return None otherwise.
    For instance, (1, 3, 2) is not admissible at index 1, because 3 > 2."""
    for i in range(len(sequence) - 1):
        if sequence[i] > sequence[i + 1]:
            return i
    return None

# Given a pair (i, j), return the sequence of summands involved in the Adem relation for QiQj if i > j, or (i, j) itself if i <= j.
# See Lawson's paper Thm 5.2 for Adem relation formula using lower indexing notation.
def adem(r: int, s: int) -> list[tuple[int]]:
    if r < s:
        return [(r, s)]
    result = []
    lower_bound = max(s + 1, ceil((r + s)/2)) 
    upper_bound = floor((r + 2*s)/2)
    for j in range(lower_bound, upper_bound + 1):
        coefficient = binomial_mod_2(j - s - 1, 2*j - r - s)
        if coefficient % 2 == 0:
            continue
        result.append((r + 2*s - 2*j, j))
    return result

# Given a sequence, apply the Adem relations at the given index to obtain a (slightly more) admissible sequence.
def adem_sequence_at(sequence: tuple[int], index: int) -> list[tuple[int]]:
  if index >= len(sequence):
      raise ValueError(f'Index {index} not present in sequence {sequence} of length {len(sequence)}.')
  if index + 1 >= len(sequence):
      raise ValueError(f'Index {index + 1} not present in sequence {sequence} of length {len(sequence)}.')
  
  start = sequence[:index]
  end = sequence[(index + 2):]
  r = sequence[index]
  s = sequence[index + 1]

  swapped = adem(r, s)
  result = []
  for pair in swapped:
      result.append(start + pair + end)
  return result

# Given a sequence of integers, recursively apply the Adem relations to obtain a list of admissible sequences.
def adem_sequence(sequence: tuple[int]) -> list[tuple[int]]:
  if len(sequence) < 2:
      return sequence
  
  finished = (is_admissible(sequence) == None)
  result = [sequence]
  
  while not finished:
      finished = True
      new_result = []
      adem_indices = [is_admissible(s) for s in result]    
      for i in range(len(result)):
          sequence = result[i]
          adem_index = adem_indices[i]
          if adem_index == None:
              new_result.append(sequence)
              continue
          finished = False
          new_result += adem_sequence_at(sequence, adem_index)

      result = new_result
  
  return result