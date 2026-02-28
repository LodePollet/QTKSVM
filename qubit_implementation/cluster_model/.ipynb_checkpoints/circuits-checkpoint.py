from qiskit import ClassicalRegister, QuantumRegister, QuantumCircuit
from unitaries import u, u1


def build_circuit(n, nphys, g, meas_basis):
    # Create a circuit with qubit and classical registers
    qbits = QuantumRegister(n, 'q')
    clbits = ClassicalRegister(nphys, 'cl')
    circuit = QuantumCircuit(qbits, clbits)

    # State preparation
    circuit = u1(circuit, g, qbits[0], qbits[1])
    for i in range(1, n - 1):
        circuit = u(circuit, g, qbits[i], qbits[i + 1])

    # Apply random single qubit transformation to x basis (Hadamard),
    # y basis (phase and Hadamard), or z basis (identity)
    for i in range(1, n - 1):
        if meas_basis[i - 1] == 0:
            circuit.h(qbits[i])
        elif meas_basis[i - 1] == 1:
            circuit.sdg(qbits[i])
            circuit.h(qbits[i])
        elif meas_basis[i - 1] == 2:
            pass

    # Projective measurement
    #     circuit.measure_all()
    for i in range(1, n - 1):
        circuit.measure(qbits[i], clbits[i - 1])

    return circuit
