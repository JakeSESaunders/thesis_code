from manifolds.hurewicz_image import hurewicz_image
from rings.bo import HomologyMO, HomologyMOManifoldBasis
from homomorphisms.image import image
from homomorphisms.matrix import inverse

def manifold_to_e_basis(x : HomologyMOManifoldBasis):
  return image(hurewicz_image, x, HomologyMO)

def e_to_manifold_basis(x : HomologyMO):
  return inverse(hurewicz_image, x, HomologyMOManifoldBasis)