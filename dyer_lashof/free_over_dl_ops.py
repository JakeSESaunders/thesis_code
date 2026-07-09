from rings.free_over_dl_ops import FreeOverDLAlg
from dyer_lashof.formal import DL
from rings.free_over_dl_ops import SModMod2
from dyer_lashof.adem import adem_sequence
from dyer_lashof.extend import extend_dl

def dl_ops_on_generators_free_over_dl_ops(R, k, dl_op: DL, upper):
  dl_op = dl_op.to_lower()
  if upper:
    k -= dl_op.degree()
  if k < 0:
    return R.zero()
    
  adem = dl_op.Q(k).adem()
  return sum(
    [R.generator(admissible_op) for admissible_op in adem],
    R.zero()
  )

def dl_ops_free(R, I, x, upper):
  return extend_dl(
    lambda _k, _dl_op, _upper: dl_ops_on_generators_free_over_dl_ops(R, _k, _dl_op, _upper),
    I,
    x,
    upper
  )

dyer_lashof_free = lambda I, x, upper: dl_ops_free(FreeOverDLAlg, I, x, upper)
dyer_lashof_s_modmod_2 = lambda I, x, upper: dl_ops_free(SModMod2, I, x, upper)
dyer_lashof_powers_of_two = lambda I, x, upper: dl_ops_free(PowersOfTwo, I, x, upper)