import numpy as np


# Unitary gate U1
def u1(circ, g, qa, qb):
    gfac = 1 / np.sqrt(1 + abs(g))
    thetaR = np.arcsin(gfac) * 2

    circ.h(qa)
    circ.cnot(qa, qb)
    circ.u(thetaR, 0, np.pi, qb)
    if g > 0:
        circ.h(qb)
        circ.cnot(qa, qb)
        circ.h(qb)
    return circ


# Unitary gate U
def u(circ, g, qa, qb):
    gfac = 1 / np.sqrt(1 + abs(g))
    thetaV = np.arcsin(gfac * np.sqrt(abs(g)))
    thetaW = np.arccos(gfac * np.sqrt(abs(g)) * np.sign(g))

    circ.x(qa)
    circ.u(thetaW, 0, 0, qb)
    circ.cnot(qa, qb)
    circ.x(qa)
    circ.u(-thetaW, 0, 0, qb)
    circ.u(thetaV, 0, 0, qb)
    circ.cnot(qa, qb)
    circ.x(qa)
    circ.u(-thetaV, 0, 0, qb)
    return circ
