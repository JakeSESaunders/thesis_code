from rings.free_over_dl_ops import FreePowersOfTwo, FormalDLOp
from rings.bo import HomologyMO
from manifolds.hurewicz_image import hurewicz_image
from homomorphisms.image import extend_multiplicatively
from dyer_lashof.bo import dyer_lashof_HMO

def free_to_MO_on_generators(dl_op: FormalDLOp) -> HomologyMO:
  """The homomorphism from H_*Y0 to H_*MO sending each Dyer--Lashof generator to the projective space in the relevant degree, on generators."""
  if not FreePowersOfTwo.is_valid_generator(dl_op):
    raise ValueError(f"{dl_op} does not describe a valid generator, so cannot give image.")

  QI, x = dl_op.QI, dl_op.x
  image_of_x = hurewicz_image(x)
  return dyer_lashof_HMO(QI, image_of_x)

def free_to_MO(x : FreePowersOfTwo) -> HomologyMO:
  """The homomorphism from H_*Y0 to H_*MO sending each Dyer--Lashof generator to the projective space in the relevant degree."""
  return extend_multiplicatively(free_to_MO_on_generators, x, HomologyMO)