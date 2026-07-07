# Your code goes here.
from rings.free_over_dl_ops import PowersOfTwo
from MO_cell_structure import get_classes_in_degree, get_good_classes_in_degree, get_classes_to_kill_in_degree, get_classes_already_killed_in_degree
from homomorphisms.free_to_MO import free_to_MO
from power_series.powers_of_2_spherical import spherical, u, Qx_2_pow_i
from power_series.dyer_lashof import dl_power_series
from dyer_lashof.free_over_dl_ops import dl_ops_on_generators_powers_of_2, dyer_lashof_powers_of_two

x1 = PowersOfTwo.x(0)
x2 = PowersOfTwo.x(1)
x4 = PowersOfTwo.x(2)

spherical_Q2x1 = spherical(PowersOfTwo.Q((2, ), 0, specified_power=True))
print(
  spherical_Q2x1
)

for i in range(1, 20):
  print(dyer_lashof_powers_of_two(i, spherical_Q2x1, upper=False))


""" print(
  dl_power_series(dl_ops_on_generators_powers_of_2, spherical_Q2x1)
) """

""" print(
  spherical(PowersOfTwo.Q((1, 2), 0, specified_power=True))
) """

""" for i in range(1, 10):
  print(i)
  print("All:       ", [str(x) for x in get_classes_in_degree(i)])
  print("Good:      ", [str(x) for x in get_good_classes_in_degree(i)])
  print("Spherical: ", [str(spherical(PowersOfTwo.Q(tuple(x.operation), x.generator, specified_power=False))) for x in get_classes_in_degree(i) if not (x.generator == 1 and (len(x.operation) == 0 or x.operation[-1] == 1))])
  print("To Kill:   ", [str(x) for x in get_classes_to_kill_in_degree(i)])
  print("Killed:    ", [str(x) for x in get_classes_already_killed_in_degree(i)])
  print() """

''' print("""\\begin{table}[]""")
print("""\\begin{tabular}{|c|c|c|c|}""")
print("""\\hline""")
print("""Degree & Class & Spherical & $f(\\text{Spherical})$ \\\\ \\hline""")
last_degree = 0
for i in range(1, 12):
  for x_data in get_classes_in_degree(i):
    # print(x_data)
    
    if (x_data.generator == 1 and (len(x_data.operation) == 0 or x_data.operation[-1] == 1)):
      continue

    x = PowersOfTwo.Q(tuple(x_data.operation), x_data.generator, specified_power=False)
    degree = x.max_degree()
    spherical_x = spherical(x)
    image_x = free_to_MO(spherical_x)
    if degree != last_degree:
      print(f"""${degree}$ & ${x}$ & ${spherical_x}$ & ${image_x}$ \\\\ \\hline""")
    
print("""\\end{tabular}""")
print("\\end{table}") '''