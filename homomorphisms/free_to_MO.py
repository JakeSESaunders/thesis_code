from rings.free_over_dl_ops import FreeOverDLAlg
from rings.bo import HomologyMO
from manifolds.hurewicz_image import hurewicz_image
from homomorphisms.image import image
from dyer_lashof.bo import dyer_lashof
from dyer_lashof.formal import DL

def free_to_MO_on_generators(dl_op: DL):
  if not FreeOverDLAlg.is_valid_generator(dl_op):
    raise ValueError(f"{dl_op} does not describe a valid generator, so cannot give image.")

  I, n = dl_op.operation, dl_op.generator
  image_of_x = hurewicz_image(n)
  return dyer_lashof(I, image_of_x, upper=False)

def free_to_MO(x : FreeOverDLAlg):
  return image(free_to_MO_on_generators, x, HomologyMO)