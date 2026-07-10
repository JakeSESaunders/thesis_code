from manifolds.hurewicz_image import hurewicz_image
from rings.bo import HomologyMO, HomologyMOManifoldBasis
from homomorphisms.image import extend_multiplicatively
from homomorphisms.matrix import inverse

def manifold_to_e_basis(x : HomologyMOManifoldBasis):
  """The homomorphism sending a bordism class of manifold to its Hurewicz image written in the e-basis of H_*MO, on generators."""
  return extend_multiplicatively(hurewicz_image, x, HomologyMO)

def e_to_manifold_basis(x : HomologyMO):
  """The homomorphism sending a bordism class of manifold to its Hurewicz image written in the e-basis of H_*MO."""
  return inverse(hurewicz_image, x, HomologyMOManifoldBasis)