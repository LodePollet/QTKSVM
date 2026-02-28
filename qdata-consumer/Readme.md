# Overview
This repository is a general client code for [TK-SVM] adapted for the use with for quantum models. See this [paper] for a theoretical description of the method. In this guide, only the technical differences between quantum and classical models will be covered. Make sure to be familiar with [TK-SVM] before reading this guide.

# Configuration and Installation
All the necessary dependencies required for installation, most importantly [ALPSCore], can be found in the [TK-SVM] repository. The client code contains [TK-SVM] as submodule which in turn contains two other submodules. To get all the submodules in one pull, use the flag `--recurse-submodules` when downloading.
```sh
git clone git@gitlab.physik.uni-muenchen.de:Nicolas.Sadoune/qdata-consumer.git --recurse-submodules
```
Once downloaded, the [TK-SVM] submodule is already on the correct branch adapted for quantum models.
For configuration you will eventually need to modify the source code in several places. These locations are often marked by a comment `//client-specific` in the source code. If you just want to run the example, no modifications to the source code are needed.

## Lattice and site type
First of all, specify which lattice you want to use in the file `inlcude/client/sim.hpp` by including the correct header for the lattice. Several lattices such as the chain, square, cubic, kagome, honeycomb and pyrochlore lattice are already implemented in the location  in the directory `include/frustmag/lattice`.
Next specify the site type, which is a basic container whose capacity must match the size of the operator basis. The operator basis is the Language in which local observables are expressed. In the follow-along guide, we consider a spin-1/2 chain and naturally choose the pauli operators X,Y,Z as operator basis. Since in this example the operator basis has size 3, we can repurpose the site type `spin_O3` originally used for the classical models. We include the header `orhto.hpp` which contains the definition for chain, square and cubic lattice.
```cpp
#include <frustmag/lattice/ortho.hpp>
#include <frustmag/site/spin_O3.hpp>
```
The size of the lattice is determined during runtime in the function `function update_phasepoint()` in the file `inlcude/client/sim.hpp`. The variable `n_line` counts the number of values in a single line of the input data, which corresponds to the number of physical sites in the system. The lattice is initialized calling the constructor of the lattice class
```cpp
lattice = {static_cast<size_t>(n_line), true, [&rng] {
    return site_type::random(rng);
}};
```
The constructor takes the number of unitcells as first argument. In case of a simple chain we have the number of unicells equals `n_line`. For illustration, the square-link latte (toric code) has `sqrt(n_line/2)` unitcells because there are two sites per unit cell and the lattice has 2 spatial dimensions. The corresponding code for the square-link lattice would be
 ```cpp
 lattice = {static_cast<size_t>(sqrt(n_line/2)), true, [&rng] {
    return site_type::random(rng);
}};
```
For more detail, see the constructor for lattice class in `include/frustmag/lattice/bravais.hpp`, which is the base class for all lattices. Almost at the end of the file we specify the `sim_base`
```cpp
using sim_base = client::sim<frustmag::lattice::chain>;
```
which is an alias for `ortho<1>`.


## POVM and data format
The data format are plain `.txt` containing one sample per line. Each sample consists of integers separated by a space character. In the sim class `inlcude/client/sim.hpp`, POVM outcomes are encoded in a one-to-one correspondence with integer values. We use the Pauli-6 POVM encoded as `0: x+, 1: x-, 2: y+, 3: y-, 4: z+, 5: z-` where `x+` is  the outcome correpsonding to the projector on the eigenvector of X with eigenvalue +1.
The sample provided under `example/samples` were obtained from ion traps realizing the the spin-1/2 family of states described in [paper]. Integer labels distinguish different data-sets. The correspondence between dataset label and g-value (Hamiltonian parameter) is
```
Run_0: g=-1, Run_1: g=-0.9, Run_2: g=-0.8, ... , Run_19: g=0.9, Run_20: g=0.99
```
If necessary, more labels can be introduced in `include/client/phasepoint.hpp`.
In the function `update()` in `include/client/sim.hpp` several POVM are already coded. To use another POVM than the Pauli-6, simply select it by commenting out other POVM definitions. For better understanding of the way that POVM are encoded in the sim class, have a look at the mathematica scripts under `POVM_definitions_mathematica`. In those mathematice notebooks you will find the construction of one SIC-POVM and one MUB-POVM for spin-1/2 and spin-1, based on the references [Decker03], [Renes03] and [Wootters89]. The choice of the POVM must be accounted for in the site_type `include/client/site/spin_O3`. In there, the function `random()` must be adapt. It is called when classifying against a set of random samples. For classical models this class is the infinite temperature class. For quantum models, we need to generate POVM outcomes uniformly. Currently selected in the `spin_O3` class is the Pauli-6 POVM selected.

## Compilation
Compile by running the following commands
```sh
mkdir build
cd build
cmake ..
make -j6
```
The result are four executables: `sample, learn, coeffs,` and `segregate_phases`. The latter three executables are very similiar to classical client codes, and can be used as described in [TK-SVM]. Instead of performing a Monte Carlo simulation, the `sample` executable reads in the the provided data, but cannot produce any data. The data has to be obtained by other means, e.g. DMRG simulation or experiments such as in the example.

In addition to the usual [TK-SVM] parameters, one new parameter needs to be specified. The new parameter, called `Nc` controlls the *sample average*. During the computation of feature vectors, the average over many clusters is taken. In case of system size restrictions, the number of clusters within a single sample is too small to get a somewhat accurate estimate of the feature vector. For example, the trapped ion data in the example consists of only 5 sites, which is merely one cluster if we are interested in the rank 5 feature vector. Therefore the average must be taken over several samples. `Nc` determines over how many samples the sample average is taken. In case of `Nc=1` the average is taken over clusters only. The cluster average is always taken automatically, and has no cotrolling parameter.

## Phase classification
Now we are ready to run the phase classification on the example data. Go to the directory `example/phase_diagram`. The provided file `phasediagram.ini` specifies all the necessary parameters.
```
## phasediagram.ini
datapath = "<path>/example/samples"
rank = 1
cluster = "1cell"
nu = 0.3
SEED = 22
[classifier]
policy = "fixed_from_sweep"
[sweep]
policy = "sweep_grid"
samples = 121500
Nc = 500
[sweep.sweep_grid.sweep1]
policy = "line_scan"
[sweep.sweep_grid.sweep1.line_scan]
N = 21
a.T = 0
b.T = 20
```
Make sure to replace `<path>` with the location of where you downloaded this repository on your local device. The label `T` is used for dataset distinction. Run the phase classification task with the following commands
```sh
<path>/build/sample phasediagram.ini
<path>/build/learn phasediagram.clone.h5 --noInfT
<path>/build/segregate-phases phasediagram.out.h5 --weight=lorentzian
python3 plot_phasediagram.py
```
The `sample` program will read in all the files starting from `Run_0.txt` to `Run_20.txt`. The flag `--noInfT` excludes the *infintie temperature* class (random samples). Succeeding, the learn program performs all possible 21*(21-1)/2 = 210 binary classification tasks. Finally, the executable `segregate-phases` constructs and partitions the graph. Note that each command takes the output of the previous command as input. At last the resulting graph and Fiedler values can be plotted with a simple python script. The resulting image `ClusterPhases.pdf` should display a clear bi-partite graph.


## Feature extraction
To determine the local features characteristic of each phase, all datasets that belong to the same phase according to the previous phase classification are pooled and shuffled. This results in two datasets, one for g < 0 which is called `Run_100.txt` and one for g > 0 called `Run_200.txt`.
Changing to the directory `example/feature_extraction`, you will find the intitialization file `SPTfeatures.ini` specifying runtime parameters to extract the local features from the topological phase (g < 0)
```
## SPTfeatures.ini
datapath = "<path>/example/samples"
rank = 4
cluster = "5cell"
nu = 0.3
SEED = 22
[classifier]
policy = "fixed_from_sweep"
[sweep]
policy = "sweep_grid"
samples = 1336000
Nc = 1000
#temperature is used for dataset distinction
[sweep.sweep_grid.sweep1]
policy = "line_scan"
[sweep.sweep_grid.sweep1.line_scan]
N = 1
a.T = 100
```
Run the feature extraction task by executing the following commands
```sh
<path>/build/sample SPTfeatures.ini
<path>/build/learn SPTfeatures.clone.h5 -i
<path>/build/coeffs SPTfeatures.out.h5
python3 plot_features.py SPTfeatures-P1-InfT.coeffs.txt 5 4
```
Contrary to before, this time only a single dataset is read by the `sample` executable. Next, a single binary classification task it performed by the `learn` executable, with the flag `-i` specifying the second class (random samples) against which the dataset `Run_100.txt` shall be classified. Finally the `coeffs` program extracts the coefficient vector which is subsequently ploted with a python script. Note that the python script also prints the analytical expression of the dominant features; to do so it takes the cluster size 5 and the rank 4 as command line arguments. To obtain the features characteristic of the trivial phase (g > 0) repeat the above steps using the parameter file `TrivialFeatures.ini`.


[//]: # (external links)

   [TK-SVM]: <https://gitlab.physik.uni-muenchen.de/tk-svm/tk-svm>
   [ALPSCore]: <https://github.com/ALPSCore/ALPSCore>
   [paper]: <https://arxiv.org/abs/2408.05017>
   [Wootters89]: <https://doi.org/10.1016/0003-4916(89)90322-9>
   [Renes03]: <https://doi.org/10.1063/1.1737053>
   [Decker03]: <https://doi.org/10.48550/arXiv.quant-ph/0308098>








