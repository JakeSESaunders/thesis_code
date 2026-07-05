from rings.free_over_dl_ops import PowersOfTwo
from rings.bo import HomologyMO
from manifolds.hurewicz_image import hurewicz_image
from homomorphisms.image import image
from dyer_lashof.bo import dyer_lashof

def free_to_MO_on_generators(P):
  if not PowersOfTwo.is_valid_generator(P):
    raise ValueError(f"{P} does not describe a valid generator, so cannot give image.")
  x, I = P
  image_of_x = hurewicz_image(2**x)
  return dyer_lashof(I, image_of_x, upper=False)

def free_to_MO(x : PowersOfTwo):
  return image(free_to_MO_on_generators, x, HomologyMO)