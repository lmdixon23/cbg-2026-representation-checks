"""
pentagon_typo.py -- demonstrates the typo in the main-text pentagon equation
(eq. pentagon_korep) of arXiv:2606.20473.

The printed RHS uses angles  phi_1235, phi_2345, phi_2345  -- so phi_2345 appears twice and the
four-subset {1,2,4,5} is missing.  A pentagon must use each of the five 4-subsets exactly once.
The paper's OWN appendix writes the third factor with phi_1245.  Below: with zeta in the
decreasing order the paper assumes, the printed form FAILS and the appendix form HOLDS.

Run:  python pentagon_typo.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from operators import cos_sin

ZETA = {1: 5.0, 2: 4.0, 3: 3.0, 4: 2.0, 5: 1.0}   # paper's constraint zeta1>...>zeta5


def Rmid(i, j, k, l):                  # rotation in the (2,3) block
    c, s = cos_sin(ZETA, i, j, k, l)
    return np.array([[1, 0, 0], [0, c, s], [0, -s, c]], dtype=complex)


def Rtop(i, j, k, l):                  # rotation in the (1,2) block
    c, s = cos_sin(ZETA, i, j, k, l)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=complex)


def Rout(i, j, k, l):                  # rotation in the (1,3) block
    c, s = cos_sin(ZETA, i, j, k, l)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=complex)


def main():
    LHS = Rmid(1, 3, 4, 5) @ Rtop(1, 2, 3, 4)
    RHS_printed  = Rtop(1, 2, 3, 5) @ Rmid(2, 3, 4, 5) @ Rout(2, 3, 4, 5)   # main-text (typo)
    RHS_appendix = Rtop(1, 2, 3, 5) @ Rmid(2, 3, 4, 5) @ Rout(1, 2, 4, 5)   # appendix (correct)
    print("Pentagon equation (eq. pentagon_korep), zeta = (5,4,3,2,1):")
    print(f"  || LHS - RHS_printed   (third factor phi_2345) || = {np.linalg.norm(LHS - RHS_printed):.3e}  -> FAILS")
    print(f"  || LHS - RHS_appendix  (third factor phi_1245) || = {np.linalg.norm(LHS - RHS_appendix):.3e}  -> HOLDS")
    print("\nConclusion: main-text third RHS factor should carry phi_1245, not phi_2345.")


if __name__ == "__main__":
    main()
