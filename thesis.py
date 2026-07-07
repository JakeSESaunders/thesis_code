from rings.bo import HomologyMO
from dyer_lashof.bo import dyer_lashof
from manifolds.hurewicz_image import hurewicz_image
from homomorphisms.manifold_to_e_basis import e_to_manifold_basis 
# This file contains the automated calculations used in my thesis.

# Remark 3.3.5
M5 = hurewicz_image(5)
print(M5)
Q1M5 = dyer_lashof(1, M5, upper=False)
print(Q1M5)
Q1M5_manifold_basis = e_to_manifold_basis(Q1M5) # NOTE this takes a while
print(Q1M5_manifold_basis)

e2 = HomologyMO.e(2)
# Lemma 4.3.8
Q4e2 = dyer_lashof(4, e2, upper=True)
print(Q4e2)

# Lemma 4.3.9
Q6e2 = dyer_lashof(6, e2, upper=True)
print(Q6e2)

# Lemma 4.3.10
Q10e2 = dyer_lashof(10, e2, upper=True)
print(Q10e2)

# Table TODO
