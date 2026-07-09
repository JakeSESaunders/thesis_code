# Your code goes here.
from rings.free_over_dl_ops import FreeOverDLAlg
from MO_cell_structure import get_classes_in_degree, get_good_classes_in_degree, get_classes_to_kill_in_degree, get_classes_already_killed_in_degree
from homomorphisms.free_to_MO import free_to_MO
from dyer_lashof.free_over_dl_ops import dyer_lashof_free
from power_series.free_over_dl_ops_spherical import spherical, X, u1, X_t_minus_1_over_t
from power_series.dyer_lashof import dl_power_series
from homomorphisms.manifold_to_e_basis import e_to_manifold_basis

def print_spherical_classes(_degree):
  for i in range(1, _degree):
    print(i)
    print("All:       ", [str(x) for x in get_classes_in_degree(i)])
    print("Good:      ", [str(x) for x in get_good_classes_in_degree(i)])
    print("Spherical: ", [str(spherical(FreeOverDLAlg.Q(tuple(x.operation), x.generator, specified_power=False))) for x in get_classes_in_degree(i) if not (x.generator == 1 and (len(x.operation) == 0 or x.operation[-1] == 1))])
    print("To Kill:   ", [str(x) for x in get_classes_to_kill_in_degree(i)])
    print("Killed:    ", [str(x) for x in get_classes_already_killed_in_degree(i)])
    print()

def print_spherical_table(_degree):
  print("""\\begin{longtable}{|c|c|C{3cm}|C{3cm}|C{3cm}|}""")
  print("""\\hline""")
  print("""Degree & Class & Spherical & $f^0(\\text{Spherical})$ & Manifold Basis\\\\ \\hline""")
  last_degree = 0
  for i in range(1, _degree):
    for x_data in get_classes_in_degree(i):
      # print(x_data)
      
      if (x_data.generator == 1 and (len(x_data.operation) == 0 or x_data.operation[-1] == 1)):
        continue

      x = FreeOverDLAlg.Q(tuple(x_data.operation), x_data.generator)
      degree = x.max_degree()
      spherical_x = spherical(x)
      image_x = free_to_MO(spherical_x)
      manifold_basis_image_x = e_to_manifold_basis(image_x)
      if degree != last_degree:
        print(f"""${degree}$ & ${x}$ & ${spherical_x}$ & ${image_x}$ & ${manifold_basis_image_x}$ \\\\ \\hline""")
      
  print("""\\end{longtable}""")

# print_spherical_table(12)

spherical_Q3x1 = spherical(FreeOverDLAlg.Q(3, 1))

print(
  e_to_manifold_basis(free_to_MO(dl_power_series(dyer_lashof_free, spherical_Q3x1).coeff(6)))
)