from power_series.power_series import PowerSeries
from collections.abc import Callable
from dyer_lashof.monomial import DyerLashofMonomial
from rings.polynomial import PolyRing

def dl_power_series(dyer_lashof: Callable[[DyerLashofMonomial, PolyRing], PolyRing], x: PolyRing):
  """Given a homogeneous polynomial x (of degree d), return the power series whose t^0, ... t^{d - 1} coefficients are 0, and t^{d + k} coefficient is Q_kx."""
  if not x.is_homogeneous():
    raise ValueError(f"Cannot produce Dyer--Lashof power series for non-homogeneous element {x}.")
  
  R = x.__class__
  degree = x.max_degree()

  def Qx_coeff(coeff):
    if coeff < degree:
      return R.zero()
    if coeff == degree:
      return x**2

    operation = DyerLashofMonomial.single(coeff - degree, upper=False)
    return dyer_lashof(
      operation,
      x
    )

  return PowerSeries(Qx_coeff)