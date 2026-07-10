from rings.polynomial import PolyRing
from dyer_lashof.monomial import DyerLashofMonomial
from combinatorics.admissible import is_tuple_of_nondecreasing_positive_integers
from combinatorics.indices import is_power_of_two

class FormalDLOp:
  """A pair consisting of an operation and a generator. Used for ease of type checking for free algebras over the Dyer--Lashof operations."""
  def __init__(self, QI: DyerLashofMonomial, x):
    """Keyword arguments:
    QI: A Dyer--Lashof monomial. Must be lower indexed and reduced (no leading zeros).
    x: An object representing a generator of the ring under the Dyer--Lashof algebra. 
    """
    if not QI.lower: raise ValueError("Upper indexing notation is not supported.")
    if not QI.is_reduced(): raise ValueError("Unreduced operations are not supported.")
    self.QI = QI
    self.x = x

  def Q(self, k: int):
    return FormalDLOp(self.QI.Q(k), self.x)

  def __lt__(self, other):
    if not isinstance(other, FormalDLOp): raise ValueError(f"Cannot compare formal Dyer--Lashof operation {self} with non-operation {other}.")
    return (self.x < other.x) or (self.QI < other.QI)

class FreeOverDLAlg(PolyRing):
  def is_valid_dl_generator(x) -> bool:
    """Return if x is a valid generator for this ring under the Dyer--Lashof operations, i.e. is Q_Ix a valid generator for this ring?"""
    return True

  @classmethod
  def is_valid_generator(cls, dl_op: FormalDLOp) -> bool:
    if not isinstance(dl_op, FormalDLOp):
      return False
    return cls.is_valid_dl_generator(dl_op.x)

  @classmethod
  def dl_generator(cls, x):
    return cls.generator(
      FormalDLOp(DyerLashofMonomial.empty(upper=False), x)
    )

  @classmethod
  def generator(cls, dl_op: FormalDLOp):
    leading_zeros, QI_shorter = dl_op.QI.strip_leading_zeros() # Deal with any leading zeros.
    return cls([(FormalDLOp(QI_shorter, dl_op.x), )])**(2**leading_zeros)

  @classmethod
  def degree_of_dl_generator(cls, x) -> int:
    """Degree of a generator for the ring under the Dyer--Lashof operations."""
    raise NotImplementedError("Base class FreeOverDLAlg does not implement degree_of_dl_generator.")

  @classmethod
  def degree(cls, dl_op: FormalDLOp) -> int:
    d = cls.degree_of_dl_generator(dl_op.x)
    return dl_op.QI.degree(d)

  def symbol_dl_generator(x) -> str:
    raise NotImplementedError("Base class FreeOverDLAlg does not implement symbol_dl_generator.")

  @classmethod
  def symbol(cls, dl_op: FormalDLOp) -> str:
    return f"{dl_op.QI}{cls.symbol_dl_generator(dl_op.x)}" # TODO implement config.

  @classmethod
  def QIx(cls, I: tuple[int], x):
    return cls.generator(
      FormalDLOp(DyerLashofMonomial(I, upper=False), x)
    )

class NatFreeOverDLAlg(FreeOverDLAlg):
  """Subclass of FreeOverDLAlg, representing a ring having at most one Dyer--Lashof generator in each positive degree."""
  def is_valid_dl_generator(n: int):
    return isinstance(n, int) and n > 0

  @classmethod
  def degree_of_dl_generator(cls, n: int) -> int:
    return n

class SModMod2(NatFreeOverDLAlg):
  """The ring H_*\S//2. The only Dyer--Lashof generator is a, in degree 1."""
  def is_valid_dl_generator(x):
    return (x == 1)

  def symbol_dl_generator(x):
    return "a"

  def a():
    return SModMod2.dl_generator(1)

  def Qa(I):
    if isinstance(I, int):
      I = (I, )
    if not isinstance(I, tuple[int]):
      raise ValueError("Operation must consist of an integer/a tuple of integers.")
    return SModMod2.QIx(I, SModMod2.a())

class FreePowersOfTwo(NatFreeOverDLAlg):
  """The ring H_*Y^0, free over the DL ops on generators in degree 1, 2, 4, 8, ... (powers of 2)."""
  def is_valid_dl_generator(x: int):
    if not isinstance(x, int):
      return False

    i = is_power_of_two(x)
    return (i is not None)

  def symbol_dl_generator(n: int):
    return f"x_{n}" # TODO make compatible with config.

  def x(n):
    return FreePowersOfTwo.dl_generator(n)