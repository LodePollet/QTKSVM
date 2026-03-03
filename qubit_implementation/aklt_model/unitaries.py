import numpy as np


sq3 = 1 / np.sqrt(3)
qp = np.exp(np.pi * 2j / 3) / np.sqrt(3)
qm = np.exp(np.pi * -2j / 3) / np.sqrt(3)

# Set of unitaries rotating into mutually unbiased bases
M = [np.array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0, 1]
               ]),
     np.array([[qm, sq3, qp, 0],
               [sq3, qm, qp, 0],
               [sq3, sq3, sq3, 0],
               [0, 0, 0, 1]
               ]),
     np.array([[sq3, qp, qm, 0],
               [qp, sq3, qm, 0],
               [sq3, sq3, sq3, 0],
               [0, 0, 0, 1]
               ]),
     np.array([[qp, qm, sq3, 0],
               [qm, qp, sq3, 0],
               [sq3, sq3, sq3, 0],
               [0, 0, 0, 1]
               ]),
     ]


# opdata: list containing instructions
# circ: the circuit to apply the instructions
# qbs: tuple of qubits, that the unitary acts on
def UfromData(opdata, circ, qbs):
    for op in opdata:
        optype = op[0]
        if optype == 0:
            circ.cx(qbs[int(op[1] - 1)], qbs[int(op[2] - 1)])
        elif optype == 1:
            circ.rx(op[1], qbs[int(op[2] - 1)])
        elif optype == 2:
            circ.ry(-op[1], qbs[int(op[2] - 1)])
        elif optype == 3:
            circ.rz(op[1], qbs[int(op[2] - 1)])
        elif optype == 5:
            continue
        else:
            raise TypeError(f'bad operation type {optype}')
    return circ