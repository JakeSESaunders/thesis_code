from config import SURROUND_INDICES_WITH_BRACES, LATEX_FRIENDLY

def get_dual_steenrod_symbol(is_xi: bool, i: int) -> str:
  """Return the symbol for the generator of the dual Steenrod algebra in degree 2**i - 1.
  
  Keyword arguments:
  is_xi: bool which determines if the symbol is xi or zeta.
  """
  if LATEX_FRIENDLY:
    result = "\xi" if is_xi else "\zeta"
  else:
    result = "ξ" if is_xi else "ζ"
  if SURROUND_INDICES_WITH_BRACES:
    return f"{result}_{{{i}}}"
  return f"{result}_{i}"
