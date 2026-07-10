from collections.abc import Callable
from collections import Counter
from config import USE_POWERS, SURROUND_INDICES_WITH_BRACES

def summand_to_string(get_symbol: Callable[[object], str], summand: tuple) -> str:
  """Return the string representing a summand of a polynomial ring.
  
  Keyword arguments:
  get_symbol: A function taking a generator and returning the symbol representing the generator.
  summand: A tuple of generators of the polynomial ring.
  """
  if len(summand) == 0:
    return "1"
  
  result = ""

  if not USE_POWERS:
    for factor in summand:
      result += get_symbol(factor)
    return result
  
  counter = Counter(summand)
  for factor in counter:
    if counter[factor] == 1:
      result += get_symbol(factor)
      continue
    if SURROUND_INDICES_WITH_BRACES:
      result += f"{get_symbol(factor)}^{{{counter[factor]}}}"
    else:
      result += f"{get_symbol(factor)}^{counter[factor]}"

  return result

def polynomial_to_string(get_symbol: Callable[[object], str], summands: list[tuple]) -> str:
  """Return the string representing an element of a polynomial ring.
  
  Keyword arguments:
  get_symbol: A function taking a generator and returning the symbol representing the generator.
  summand: A list of tuples of generators of the polynomial ring.
  """
  if len(summands) == 0:
    return "0"

  result = ""

  for summand in summands:
    summand_string = summand_to_string(get_symbol, summand)
    if len(result) == 0:
      result = summand_string
      continue
    result += f" + {summand_string}"

  return result