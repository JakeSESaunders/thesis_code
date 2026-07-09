from power_series.power_series import PowerSeries

# Let x be a homogeneous element of a PolyRing R.
# Suppose dyer_lashof sends (I, x, upper) to Q_Ir or Q^Ir (depending on upper).
# This returns the power series x^2t^{\abs{x}} + Q_1xt^{\abs{x} + 1} + \cdots. 
def dl_power_series(dyer_lashof, x):
  if not x.is_homogeneous():
    raise ValueError(f"Cannot produce Dyer--Lashof power series for non-homogeneous element {x}.")
  
  R = x.__class__
  degree = x.max_degree()

  def Qx_coeff(coeff):
    if coeff < degree:
      return R.zero()
    if coeff == degree:
      return x**2

    operation = coeff - degree # TODO should this be a plus?
    return dyer_lashof(
      operation,
      x,
      upper=False
    )

  return PowerSeries(Qx_coeff)