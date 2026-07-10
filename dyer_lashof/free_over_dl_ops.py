from rings.free_over_dl_ops import FreeOverDLAlg, FormalDLOp, SModMod2, FreePowersOfTwo
from dyer_lashof.monomial import DyerLashofMonomial
from dyer_lashof.extend import extend_dl

def dl_ops_on_generators_free_over_dl_ops(R: type[FreeOverDLAlg], k: int, upper: bool, dl_op: FormalDLOp):
  # convert to lower indexing
  if upper:
    d = R.degree(dl_op)
    k -= d

  if k < 0:
    return R.zero()

  Qk = DyerLashofMonomial.single(k, upper=False)
  QI = dl_op.QI
  x = dl_op.x

  QkI: DyerLashofMonomial = Qk + QI
  if len(QI.operation) == 0 or k <= QI.operation[0]:
    return R.QIx(QkI, x)

  adem : list[DyerLashofMonomial] = QkI.adem()
  return sum(
    [R.generator(FormalDLOp(QJ, x)) for QJ in adem],
    R.zero()
  )

def dyer_lashof_free(QI, x):
  return extend_dl(
  lambda _k, _upper, _x: dl_ops_on_generators_free_over_dl_ops(FreeOverDLAlg, _k, _upper, _x),
  QI,
  x
)

def dyer_lashof_s_modmod_2(QI, x):
  return extend_dl(
  lambda _k, _upper, _x: dl_ops_on_generators_free_over_dl_ops(SModMod2, _k, _upper, _x),
  QI,
  x
)

def dyer_lashof_free_powers_of_two(QI, x):
  return extend_dl(
  lambda _k, _upper, _x: dl_ops_on_generators_free_over_dl_ops(FreePowersOfTwo, _k, _upper, _x),
  QI,
  x
)
