from power_series.power_series import PowerSeries
from dyer_lashof.extend import extend_dl

# Let x be a homogeneous element of a PolyRing R.
# Suppose dl_fn is a function sending (i, r, upper) to Q_ir or Q^ir (depending on upper).
# This returns the power series x^2t^{\abs{x}} + Q_1xt^{\abs{x} + 1} + \cdots. 
def dl_power_series(dl_ops_for_generators, x):
  if not x.is_homogeneous():
    raise ValueError(f"Cannot produce Dyer--Lashof power series for non-homogeneous element {x}.")
  
  R = x.__class__
  degree = x.max_degree()

  def Qx_coeff(coeff):
    if coeff < degree:
      return R.zero()
    if coeff == degree:
      return x**2

    operation = coeff - degree
    return extend_dl(
      dl_ops_for_generators,
      operation,
      x,
      upper=False
    )

  return PowerSeries(Qx_coeff)