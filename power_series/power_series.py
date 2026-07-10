from combinatorics.partitions import get_unordered_partitions, get_weighted_partitions
from combinatorics.multinomial import multinomial
from power_series.hessenberg import hessenberg_det
from config import NON_CONSTANT_TERMS_OF_POWER_SERIES_TO_DISPLAY
from collections.abc import Callable
from rings.polynomial import PolyRing

class PowerSeries:
  """A power series with coefficients in a polynomial ring"""
  def __init__(self, coeff: Callable[[int], PolyRing]):
    """Keyword arguments:
    coeff: function int -> R which sends i to the i-coefficient of the power series."""
    self.coeff = coeff
    self.R = coeff(0).__class__ # NOTE vulnerable to the series having mixed coefficients!

  def __add__(self, other):
    if not isinstance(other, PowerSeries):
      raise ValueError(f"Cannot add power series {self} to non-power series {other}.")
    if other.R != self.R:
      raise ValueError(f"Cannot add power series with coefficients {self.R} to power series with coefficients {other.R}.")

    return PowerSeries(lambda i: self.coeff(i) + other.coeff(i))

  def __sub__(self, other):
    return self + other # mod 2

  def __mul__(self, other):
    if not isinstance(other, PowerSeries):
      raise ValueError(f"Cannot multiply power series {self} with non-power series {other}.")
    if other.R != self.R:
      raise ValueError(f"Cannot multiply power series with coefficients {self.R} with power series with coefficients {other.R}.")

    def new_coeff(n):
      return sum([self.coeff(i)*other.coeff(n - i) for i in range(0, n + 1)], self.R.zero())

    return PowerSeries(new_coeff)

  def __str__(self):
    R = self.R
    c = self.coeff(0)
    if c == R.zero():
      result = "0"
    elif c == R.one():
      result = "1"
    else:
      result = f"({self.coeff(0)})"
    for i in range(1, 1+ NON_CONSTANT_TERMS_OF_POWER_SERIES_TO_DISPLAY):
      c = self.coeff(i)
      if c == R.zero():
        result = f"{result} + 0t^{i}"
      elif c == R.one():
        result = f"{result} + t^{i}"
      else:
        result = f"{result} + ({self.coeff(i)})t^{i}"
      
    result = f"{result} + ..."
    return result
    
  def constant(r: PolyRing):
    """Return the power series with r as constant coefficient and all other coefficients zero."""
    R = r.__class__
    return PowerSeries(lambda i: r if i == 0 else R.zero())

  def shift(self, n: int):
    """Shift the power series, so that the constant coefficient becomes the t^n coefficient.
    Works for negative n too, provided the 0, 1, ..., n - 1 coefficients are all zero.
    """
    R = self.R
    if n == 0:
      return self
    
    if n < 0:
      for i in range(n):
        if self.coeff(i) != R.zero():
          raise ValueError("Shifting power series to have nonzero negative coefficients is not supported.")

    return PowerSeries(lambda i: self.coeff(i - n) if i > n else R.zero())

  def __pow__(self, other):
    R = self.R
    if not isinstance(other, int):
      return False

    if other == 0:
      return PowerSeries.constant(R.one())

    n = abs(other)
    ratio = self if n > 0 else self.multiplicative_inverse()
    result = R.one()
    for i in range(n):
      result *= ratio

    return result

  # 
  def composition_inverse(self):
    """Obtain the composition inverse using Lagrange inversion.
    Only valid for power series with constant coeff 0 and linear coeff 1"""
    R = self.R
    if self.coeff(0) != R.zero():
      raise ValueError("Cannot find composition inverse of power series with non-zero constant coefficient.")
    if self.coeff(1) != R.one():
      raise ValueError("Cannot find composition inverse of power series with non-unit linear coefficient.")

    def coeffs_of_comp_inverse(d):
      if d == 0:
        return R.zero()
      if d == 1:
        return R.one()
      
      n = d - 1

      result = R.zero()
      for partition in get_weighted_partitions(n):
        if (multinomial(list(partition) + [n])//(n + 1)) % 2 == 0:
          continue
        factor = R.one()
        for i in range(len(partition)):
          lambda_i = partition[i]
          factor *= (self.coeff(i + 2)**lambda_i)
        result += factor
      return result

    return PowerSeries(coeffs_of_comp_inverse)

  # NOTE could redo this to use enumerative formula from Aguiar--Ardila
  def multiplicative_inverse(self):
    """Obtain the multiplicative inverse of this power series using Hessenberg determinants.
    Only valid for power series with constant coefficient 1.
    See e.g. Inselberg - On determinants of Toeplitz--Hessenberg matrices arising in power series."""
    R = self.R
    if self.coeff(0) != R.one():
      raise ValueError("Cannot find multiplicative inverse of power series with non-unit constant coefficient.")
    def coeffs_of_mult_inverse(n):
      if n == 0:
        return R.one()
      return hessenberg_det(
        [self.coeff(i) for i in range(1, n + 1)]
      )
    return PowerSeries(coeffs_of_mult_inverse)

  def circ(self, other):
    """If self = f(t) and other = g(t), returns the composite (f \circ g)(t)."""
    """Uses Faa di Bruno (see e.g https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula#Formal_power_series_version)"""
    if not isinstance(other, PowerSeries):
      raise ValueError(f"Cannot compose power series {self} with non-power series {other}.")
    if other.R != self.R:
      raise ValueError(f"Cannot compose power series with coefficients {self.R} with power series with coefficients {other.R}.")

    R = self.R
    def coeffs_of_circ(n):
      coeff = R.zero()

      for partition in get_unordered_partitions(n):
        product = self.coeff(len(partition))
        for part in partition:
          product *= other.coeff(part)
        coeff += product

      return coeff

    return PowerSeries(coeffs_of_circ)