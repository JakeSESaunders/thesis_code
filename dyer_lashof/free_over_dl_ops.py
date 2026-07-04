from rings.free_over_dl_ops import SModMod2, PowersOfTwo
from dyer_lashof.adem import adem_sequence
from dyer_lashof.extend import extend_dl

def dl_ops_on_generators_s_modmod_2(k, I, upper=False):
  if upper == True:
    k -= SModMod2.degree(I)
  inadmissible_operation = tuple([k] + list(I))
  if k <= I[0]:
    return SModMod2.generator(inadmissible_operation)
  
  admissible_operations = adem_sequence(inadmissible_operation)
  result = SModMod2.zero()
  for operation in admissible_operations:
    result += SModMod2.generator(operation)

  return result

dyer_lashof_s_modmod_2 = lambda I, x, upper: extend_dl(dl_ops_on_generators_s_modmod_2, I, x, upper)

def dl_ops_on_generators_powers_of_2(k, P, upper=False):
  if upper == True:
    k -= PowersOfTwo.degree(P)
  
  x, I = P
  inadmissible_operation = tuple([k] + list(I))
  if len(I) == 0 or k <= I[0]:
    return PowersOfTwo.generator((x, inadmissible_operation))
  
  admissible_operations = adem_sequence(inadmissible_operation)
  result = PowersOfTwo.zero()
  for operation in admissible_operations:
    result += PowersOfTwo.generator((x, operation))

  return result

dyer_lashof_powers_of_two = lambda I, x, upper: extend_dl(dl_ops_on_generators_powers_of_2, I, x, upper)