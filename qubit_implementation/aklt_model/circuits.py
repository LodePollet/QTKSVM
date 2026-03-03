import os
from typing import Literal

import numpy as np
from qiskit import ClassicalRegister, QuantumRegister, QuantumCircuit

from unitaries import UfromData, M

GValue = Literal[
    -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, -0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1
]


def build_circuit(n: int, nphys: int, g: GValue,
                  meas_basis: list[Literal[0, 1, 2, 3]]) -> QuantumCircuit:
    # Create a circuit with qubit and classical registers
    qbits = QuantumRegister(n, 'q')
    clbits = ClassicalRegister(nphys, 'cl')
    circuit = QuantumCircuit(qbits, clbits)

    # State preparation
    U1data = np.genfromtxt(os.path.join('circuit_defs', f'U1_circ_g{g}.txt'))
    Udata = np.genfromtxt(os.path.join('circuit_defs', f'U_circ_g{g}.txt'))

    print(U1data)
    print(Udata)

    circuit = UfromData(U1data, circuit, [qbits[1], qbits[0]])
    for i in range(1, n - 2, 2):
        circuit = UfromData(Udata, circuit, [qbits[i + 2], qbits[i + 1], qbits[i]])

    # Basis rotation. Note as defined above we have to apply M[i]^\dagger,
    # because the columns of M[i] are the MUB-states, and not the rows.
    for i in range(1, n - 1, 2):
        circuit.unitary(np.conj(M[meas_basis[(i - 1) // 2]]).T, [qbits[i], qbits[i + 1]])

    # Projective measurement
    for i in range(1, n - 1):
        circuit.measure(qbits[i], clbits[i - 1])

    return circuit


if __name__ == '__main__':
    qc = build_circuit(8, 6, -1, [0, 0, 0])
    print(qc.draw())
