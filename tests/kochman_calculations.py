from rings.homology_bo import PolySW
from rings.bo import HomologyMO
from dyer_lashof.bo import dyer_lashof_single_indeterminate
from config import PRINT_FAILED_TESTS

def E(i):
  return HomologyMO.e(i)

# Test against Kochman's results

# We populate the dictionary in the order based on Kochman's calculations on p133-134 of his paper.
KOCHMAN_CALCULATIONS = {}
KOCHMAN_CALCULATIONS[1] = {}
KOCHMAN_CALCULATIONS[2] = {}
KOCHMAN_CALCULATIONS[3] = {}
KOCHMAN_CALCULATIONS[4] = {}
KOCHMAN_CALCULATIONS[5] = {}

KOCHMAN_CALCULATIONS[1][1] = E(1)**2
KOCHMAN_CALCULATIONS[1][2] = E(3) + E(2)*E(1) + E(1)**3
KOCHMAN_CALCULATIONS[2][2] = E(2)**2
KOCHMAN_CALCULATIONS[1][3] = E(1)**4
KOCHMAN_CALCULATIONS[2][3] = E(5) + E(4)*E(1) + E(3)*E(2) + E(2)**2*E(1)
KOCHMAN_CALCULATIONS[1][4] = E(5) + E(4)*E(1) + E(3)*E(2) + E(3)*E(1)**2 + E(2)**2*E(1) + E(2)*E(1)**3 + E(1)**5
KOCHMAN_CALCULATIONS[3][3] = E(3)**2
KOCHMAN_CALCULATIONS[2][4] = E(6) + E(5)*E(1) + E(4)*E(2) + E(4)*E(1)**2 + E(3)*E(2)*E(1) + E(2)**3 + E(2)**2*E(1)**2
KOCHMAN_CALCULATIONS[1][5] = E(3)**2 + E(2)**2*E(1)**2 + E(1)**6
KOCHMAN_CALCULATIONS[3][4] = E(7) + E(6)*E(1) + E(5)*E(2)  + E(4)*E(3) + E(3)**2*E(1)
KOCHMAN_CALCULATIONS[2][5] = E(5)*E(1)**2 + E(4)*E(1)**3 + E(3)*E(2)*E(1)**2 + E(2)**2*E(1)**3
KOCHMAN_CALCULATIONS[1][6] = E(7) + E(6)*E(1) + E(5)*E(2) + E(4)*E(3) + E(5)*E(1)**2 + E(3)**2*E(1) + E(3)*E(2)**2 + E(4)*E(1)**3 + E(3)*E(2)*E(1)**2 + E(2)**3*E(1) + E(3)*E(1)**4 + E(2)*E(1)**5 + E(1)**7
KOCHMAN_CALCULATIONS[4][4] = E(4)**2
KOCHMAN_CALCULATIONS[3][5] = E(5)*E(3) + E(6)*E(1)**2 + E(5)*E(2)*E(1) + E(4)*E(3)*E(1) + E(3)**2*E(2) + E(3)**2*E(1)**2
KOCHMAN_CALCULATIONS[2][6] = E(5)*E(3) + E(5)*E(2)*E(1) + E(4)*E(3)*E(1) + E(3)**2*E(2) + E(5)*E(1)**3 + E(4)*E(2)*E(1)**2 + E(2)**4 + E(4)*E(1)**4 + E(3)*E(2)*E(1)**3 + E(2)**3*E(1)**2 + E(2)**2*E(1)**4
KOCHMAN_CALCULATIONS[1][7] = E(1)**8
KOCHMAN_CALCULATIONS[4][5] = E(9) + E(8)*E(1) + E(7)*E(2) + E(6)*E(3) + E(5)*E(4) + E(4)**2*E(1)
KOCHMAN_CALCULATIONS[3][6] = E(6)*E(3) + E(6)*E(2)*E(1) + E(5)*E(3)*E(1) + E(5)*E(2)**2 + E(4)*E(3)*E(2) + E(3)**3 + E(6)*E(1)**3 + E(5)*E(2)*E(1)**2 + E(4)*E(3)*E(1)**2 + E(3)**2*E(1)**3
KOCHMAN_CALCULATIONS[2][7] = E(9) + E(8)*E(1) + E(7)*E(2) + E(6)*E(3) + E(5)*E(4) + E(5)*E(2)**2 + E(4)**2*E(1) + E(4)*E(2)**2*E(1) + E(3)*E(2)**3 + E(5)*E(1)**4 + E(2)**4*E(1) + E(4)*E(1)**5 + E(3)*E(2)*E(1)**4 + E(2)**2*E(1)**5
KOCHMAN_CALCULATIONS[1][8] = E(9) + E(8)*E(1) + E(7)*E(2) + E(6)*E(3) + E(5)*E(4) + E(7)*E(1)**2 + E(5)*E(2)**2 + E(4)**2*E(1) + E(3)**3 + E(6)*E(1)**3 + E(5)*E(2)*E(1)**2 + E(4)*E(3)*E(1)**2 + E(4)*E(2)**2*E(1) + E(3)*E(2)**3 + E(3)**2*E(2)*E(1) + E(5)*E(1)**4 + E(2)**4*E(1) + E(3)*E(2)*E(1)**4 + E(4)*E(1)**5 + E(3)*E(1)**6 + E(2)**2*E(1)**5 + E(2)*E(1)**7 + E(1)**9
KOCHMAN_CALCULATIONS[5][5] = E(5)**2
KOCHMAN_CALCULATIONS[4][6] = E(10) + E(9)*E(1) + E(8)*E(2) + E(6)*E(4) + E(8)*E(1)**2 + E(7)*E(2)*E(1) + E(6)*E(3)*E(1) + E(5)*E(4)*E(1) + E(4)**2*E(2) + E(4)**2*E(1)**2
# There is a typo in Kochman's calculation below: the term a_4a_3a_1^2 is degree 9, when the summand should have degree 10.
KOCHMAN_CALCULATIONS[3][7] = E(5)*E(3)*E(1)**2 + E(3)**2*E(2)**2 + E(6)*E(1)**4 + E(5)*E(2)*E(1)**3 + E(4)*E(3)*E(1)**3 + E(3)**2*E(2)*E(1)**2 + E(3)**2*E(1)**4
# There is a typo: the 3rd summand on the 4th line should be two summands: a_3^2a_1^4 and a_3a_2^2a_1^3.
KOCHMAN_CALCULATIONS[2][8] = E(10) + E(9)*E(1) + E(8)*E(2) + E(6)*E(4) + E(8)*E(1)**2 + E(7)*E(2)*E(1) + E(6)*E(2)**2 + E(6)*E(3)*E(1) + E(5)*E(4)*E(1) + E(4)**2*E(2) + E(5)*E(3)*E(1)**2 + E(5)*E(2)**2*E(1) + E(4)**2*E(1)**2 + E(4)*E(2)**3 + E(5)*E(2)*E(1)**3 + E(4)*E(3)*E(1)**3 + E(4)*E(2)**2*E(1)**2 + E(3)**2*E(2)*E(1)**2 + E(3)*E(2)**3*E(1) + E(2)**5 + E(5)*E(1)**5 + E(4)*E(2)*E(1)**4 + E(3)**2*E(1)**4 + E(3)*E(2)**2*E(1)**3 + E(2)**4*E(1)**2 + E(4)*E(1)**6 + E(3)*E(2)*E(1)**5 + E(3)*E(2)*E(1)**5 + E(2)**2*E(1)**6
# Another typo: a_10^2 should be a_5^2, otherwise it has the wrong degree.
# In the 4th summand, the computer says a_2^3a_1^4 should be a_3^2a_1^4. It seems reasonable to believe that this was a typo.
KOCHMAN_CALCULATIONS[1][9] = E(5)**2 + E(4)**2*E(1)**2 + E(3)**2*E(2)**2 + E(2)**3*E(1)**4 + E(2)**4*E(1)**2 + E(2)**2*E(1)**6 + E(1)**10

def check_against_kochman_calculations():
  for k in KOCHMAN_CALCULATIONS.keys():
    for n in KOCHMAN_CALCULATIONS[k].keys():
      our_calculation = dyer_lashof_single_indeterminate(n, k)
      kochman_calculation = KOCHMAN_CALCULATIONS[k][n]
      comparison = (our_calculation == kochman_calculation)
      message = f"Dyer--Lashof operation Q^{n}e_{k} as computed by the program does not agree with Kochman's computation."
      if PRINT_FAILED_TESTS and (comparison == False):
        print(message)
        print(f"Kochman's calculation: {kochman_calculation}")
        print(f"Our calculation: {our_calculation}")
        print(f"Difference: {kochman_calculation - our_calculation}")
        print()
      else:
        assert comparison, message 
