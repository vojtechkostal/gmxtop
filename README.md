# gmxtop

## Description
`gmxtop` is a Python package for parsing Gromacs topology files

## Installation
1. Clone this repository
2. OPTIONAL: create a new `conda`/`mamba` environment
3. install the package using `pip`
```bash
pip install .
```

## Usage
#### see `example.ipynb` file in *examples*
```python
from gmxtop import Topology

# load topology from file
fn_top = "./topol.top"
top = Topology(fn_top)

# modifying topology parameters
for atom in top.molecules['POPC'][0].atoms:
    atom.update(charge=0.0)

for atomtype in top.atomtypes:
    atomtype.update(sigma=0.31)

# remove virtual sites from a molecule
top.molecules['SOL'][0].remove_vsites()

# write topology into a file
top.write("./topol-processed.top", overwrite=True)
```