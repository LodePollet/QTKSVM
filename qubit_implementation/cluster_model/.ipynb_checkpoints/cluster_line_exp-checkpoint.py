from tqdm import tqdm
import qiskit
import numpy as np
import yaml
from qiskit import transpile
from qiskit.providers.aer import AerSimulator
from qiskit.quantum_info import hellinger_fidelity
from qiskit_aqt_provider import AQTProvider, aqt_pass_manager
from circuits import build_circuit

date = '2023_07_30'

N = 7
nshots_in_basis = 500
runs = [9, 10]

gs = [
    -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
    -0.01,
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
]

Nphys = N - 2
nshots = nshots_in_basis * 3 ** Nphys

if __name__ == '__main__':

    aqt = AQTProvider('MY_TOKEN')
    aqt_backend = aqt.backends.aqtion_ms_minus
    aqt_pass_mng = aqt_pass_manager(lowest_angle=np.pi / 15)
    backend_sim = qiskit.Aer.get_backend('qasm_simulator')
    simulator = AerSimulator()
    for run in runs:
        for g in gs:
            try:
                # random measurement basis array. 0: x-basis, 1: y-basis, 2:z-basis
                meas_basis = np.random.randint(3, size=(nshots, Nphys))
                meas_basis_unique = np.unique(meas_basis, axis=0)
                # meas_basis_unique = meas_basis_unique[:]

                # idx = np.random.randint(len(meas_basis_unique), size=300)
                # meas_basis_unique = meas_basis_unique[idx]

                result_list = []
                nshots_final = 0  # can become slightly different than nshots, due to rounding
                for mb in tqdm(meas_basis_unique):
                    circuit = build_circuit(N, Nphys, g, mb)

                    # Execute
                    shots_in_basis = int(np.ceil(nshots / meas_basis_unique.shape[0]))

                    # Sim
                    compiled_circuit = transpile(circuit, simulator)
                    result_sim = simulator.run(compiled_circuit, shots=shots_in_basis).result()

                    # Exp
                    trans_qc = aqt_pass_mng.run(circuit)
                    result_exp = aqt_backend.run(trans_qc, shots=shots_in_basis).result()
                    # result_exp = simulator.run(compiled_circuit, shots=shots_in_basis).result()

                    # compiled_circuit = transpile(trans_qc, simulator)
                    # result_exp = simulator.run(trans_qc, shots=shots_in_basis).result()

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

                    file_name = f"data/{date}_cluster_exp_{N}qubits_g{str(g).replace('.', '')}_{nshots}shots_{run}.yaml"

                    with open(file_name, 'w') as file:
                        t = yaml.dump(result_list, file)
            except:
                pass
