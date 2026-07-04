def is_tuple_of_nondecreasing_positive_integers(I):
  if not isinstance(I, tuple):
      return False

  if len(I) == 0:
    return True # the condition is superflous

  previous = 1
  for i in range(len(I)):
    if not isinstance(I[i], int):
      return False
    if I[i] < previous:
      return False
    previous = I[i]
  return True