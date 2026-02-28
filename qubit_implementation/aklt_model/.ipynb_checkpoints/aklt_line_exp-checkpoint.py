import numpy as np
from tqdm import tqdm
import yaml

import qiskit
from qiskit import transpile, ClassicalRegister, QuantumRegister, QuantumCircuit
from qiskit.providers.aer import AerSimulator

from qiskit.quantum_info import hellinger_fidelity
from qiskit_aqt_provider import AQTProvider, aqt_pass_manager

date = '2023_07_28'

# total number of qubits, edge qubits are unphysical and should not be measured.
# total number of measured qutrits is Nphys/2 == N/2 - 1
N = 8
nshots_in_basis = 2000
runs = [2]

assert N % 2 == 0
Nphys = N - 2
# Number of projective measurements
# nshots = 1000000
nshots = nshots_in_basis * 4 ** (Nphys // 2)

gs = [
    -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
    -0.01,
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99,
]

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
            circ.cnot(qbs[int(op[1] - 1)], qbs[int(op[2] - 1)])
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


if __name__ == '__main__':

    aqt = AQTProvider('MY_TOKEN')
    aqt_backend = aqt.backends.aqtion_ms_minus
    aqt_pass_mng = aqt_pass_manager(lowest_angle=np.pi / 15)
    backend_sim = qiskit.Aer.get_backend('qasm_simulator')
    simulator = AerSimulator()

    meas_basis = np.random.randint(4, size=(nshots, Nphys // 2))
    meas_basis_unique = np.unique(meas_basis, axis=0)
    for run in runs:
        for g in gs:
            result_list = []
            nshots_final = 0  # can become slightly different than nshots, due to rounding
            for mb in tqdm(meas_basis_unique, desc=f'g = {g:.2f}'):
                # Create a circuit with qubit and classical registers
                qbits = QuantumRegister(N, 'q')
                clbits = ClassicalRegister(Nphys, 'cl')
                circuit = QuantumCircuit(qbits, clbits)

                # State preparation
                U1data = np.genfromtxt(f'./circuit-defs-9cnot/U1_circ_g{g}.txt')
                Udata = np.genfromtxt(f'./circuit-defs-9cnot/U_circ_g{g}.txt')
                circuit = UfromData(U1data, circuit, [qbits[1], qbits[0]])
                for i in range(1, N - 2, 2):
                    circuit = UfromData(Udata, circuit, [qbits[i + 2], qbits[i + 1], qbits[i]])

                # Basis rotation. Note as defined above we have to apply M[i]^\dagger,
                # because the columns of M[i] are the MUB-states, and not the rows.
                for i in range(1, N - 1, 2):
                    circuit.unitary(np.conj(M[mb[(i - 1) // 2]]).T, [qbits[i], qbits[i + 1]])

                # Projective measurement
                for i in range(1, N - 1):
                    circuit.measure(qbits[i], clbits[i - 1])

                shots_in_basis = int(np.ceil(nshots / meas_basis_unique.shape[0]))

                # Sim
                compiled_circuit = transpile(circuit, simulator)
                result_sim = simulator.run(compiled_circuit, shots=shots_in_basis).result()

                # Exp
                circuit_1 = transpile(circuit, aqt_backend)
                trans_qc = aqt_pass_mng.run(circuit_1)
                result_exp = aqt_backend.run(trans_qc, shots=shots_in_basis).result()
                # result_exp = simulator.run(compiled_circuit, shots=shots_in_basis).result()

                result_entry = {'Basis': [int(b) for b in mb],
                                'Qubits': int(N),
                                'g': float(g),
                                'Shots': shots_in_basis,
                                'Sim': dict(result_sim.get_counts()),
                                'Exp': dict(result_exp.get_counts()),
                                'Hellinger fidelity': float(
                                    hellinger_fidelity(result_sim.get_counts(),
                                                       result_exp.get_counts()))
                                }
                result_list.append(result_entry)

                nshots_final += shots_in_basis

                file_name = f"data-9cnot/{date}_aklt_exp_{N}qubits_g{str(g).replace('.', '')}_{nshots}shots_{run}.yaml"

                with open(file_name, 'w') as file:
                    t = yaml.dump(result_list, file)
