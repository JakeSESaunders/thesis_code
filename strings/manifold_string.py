from config import DISPLAY_DOLD_MANIFOLD_INDICES_MANIFOLD_BASIS, SURROUND_INDICES_WITH_BRACES, USE_BB_R_FOR_PROJECTIVE_SPACE
from combinatorics.indices import is_one_less_than_power_of_two, get_dold_manifold_indices_for_degree

def get_manifold_basis_symbol(i: int) -> str:
  """Return the symbol for the generator of the manifold basis in degree i."""
  if i % 2 == 0:
    return get_projective_space_symbol(i)

  index = is_one_less_than_power_of_two(i)
  if index is None:
    return get_dold_manifold_symbol(i)

  return get_Q1_symbol(index)

def get_Q1_symbol(i: int) -> str:
  """Return the symbol representing Q_1^{i}e_1."""
  result = "e_"
  if SURROUND_INDICES_WITH_BRACES:
    result += "{1}"
  else:
    result += "1"
  
  if i == 0:
    return result
  if i == 1:
    if SURROUND_INDICES_WITH_BRACES:
      result = "Q_{1}" + result
    else:
      result = "Q_1" + result
    return result
    
  if SURROUND_INDICES_WITH_BRACES:
    result = f"Q_{{{1}}}^{{{i}}}" + result
  else:
    result = f"Q_1^{i}" + result
  return result

def get_projective_space_symbol(i: int) -> str:
  """Return the symbol representing \RP^i."""
  result = "["
  if USE_BB_R_FOR_PROJECTIVE_SPACE:
    result += "\RP^"
  else:
    result += "RP^"
  if SURROUND_INDICES_WITH_BRACES:
    result += f"{{{i}}}"
  else:
    result += f"{i}"
  result += "]"
  return result

def get_dold_manifold_symbol(i: int) -> str:
  """Return the symbol representing the Dold manifold generator in degree i."""
  if DISPLAY_DOLD_MANIFOLD_INDICES_MANIFOLD_BASIS:
    r, s = get_dold_manifold_indices_for_degree(i)
    if SURROUND_INDICES_WITH_BRACES:
      return f"[P({{{r}}}, {{{s}}})]"
    else:
      return f"[P({r}, {s})]"
  
  if SURROUND_INDICES_WITH_BRACES:
    return f"[M_{{{i}}}]"
  return f"[M_{i}]"