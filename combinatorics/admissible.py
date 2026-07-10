def is_tuple_of_nondecreasing_integers_greater_than_or_equal_to(n: int, I: tuple[int]) -> bool:
  if not isinstance(I, tuple):
      return False

  if len(I) == 0:
    return True # the condition is superflous

  previous = n
  for i in range(len(I)):
    if not isinstance(I[i], int):
      return False
    if I[i] < previous:
      return False
    previous = I[i]
  return True

def is_tuple_of_nondecreasing_nonnegative_integers(I: tuple[int]) -> bool:
  return is_tuple_of_nondecreasing_integers_greater_than_or_equal_to(0, I)

def is_tuple_of_nondecreasing_positive_integers(I: tuple[int]) -> bool:
  return is_tuple_of_nondecreasing_integers_greater_than_or_equal_to(1, I)