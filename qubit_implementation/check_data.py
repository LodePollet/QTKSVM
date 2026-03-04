import os
import numpy as np
import yaml

if __name__ == '__main__':

    for model in ['cluster_model', 'aklt_model']:
        file_names = os.listdir(os.path.join(model, 'data'))

        print(f'{model} data:')
        print(f"{'File':55} {'Mean':5} | {'Min':5} | {'Max':5}")

        for name in file_names:
            with open(os.path.join(model, 'data', name), 'r') as file:
                # Loading the data from the .yaml file
                result = yaml.load(file, Loader=yaml.FullLoader)

                # Each file contains measurement outcomes for the trapped-ion experiment and the Qiskit simulation
                # for each basis for one g value

            fid_mean = np.mean([f['Hellinger fidelity'] for f in result])
            fid_min = np.min([f['Hellinger fidelity'] for f in result])
            fid_max = np.max([f['Hellinger fidelity'] for f in result])
            print(f'{name:55} {fid_mean:4.3f} | {fid_min:4.3f} | {fid_max:4.3f}')

        print()
