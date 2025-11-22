# rdfadf Package

## 1. Summary

`rdfadf` is for computing  
**RDF (Radial Distribution Function)** and  
**ADF (Angular Distribution Function)**  
for crystalline and amorphous materials with high numerical accuracy.

This script was developed for educational purposes and has not been optimized for computational speed. As a result, it may become significantly slower when applied to very large systems with a high number of particles. 
This script is primarily built using the functionality provided by the Atomic Simulation Environment (ASE) https://ase-lib.org/.

Key features:

- Fully **PBC-aware RDF / ADF** computation
- Native support for ASE (`ase.Atoms`)
- Flexible selection of central and neighbor atoms using:
  **element symbols / block (s,p,d,f) / group (G#) / period (P#) / All**
- RDF:
  - Normalized per central atom
  - Optional 4πr² normalization
  - Cumulative counts → direct coordination number (CN)
- ADF:
  - Distance filtering (r_min, r_max)
  - Normalized per central atom
  - Optional Gaussian smoothing
- Optimized for research, MD trajectory analysis, and materials informatics workflows

---

## 2. Installation

Install the package from the project root directory:

```bash
pip install .
```

Uninstall:

```bash
pip uninstall rdfadf -y
```

---

## 3. Inputs

The package accepts **ASE `Atoms` objects**:

```python
from ase.io import read

atoms = read("POSCAR")
atoms = read("structure.cif")
atoms = read("output.traj")
```

---

## 4. Label Specification

You can classify central and neighbor atoms using:

| Category | Example | Meaning |
|---------|----------|---------|
| Element | "Li", "O", "Co" | Specific element |
| All | "All" | All elements in the structure |
| Block | "s", "p", "d", "f" | s/p/d/f block |
| Group | "G16" | Group 16 |
| Period | "P3" | Period 3 |

Examples:

- `"Li"`–`"O"` → Li–O RDF  
- `"Co"`–`"All"` → Co around all elements  
- `"d"`–`"p"` → d-block vs p-block atoms  
- `"G1"`–`"O"` → Group 1 elements around O  

---

## 5. RDF Function

```
compute_general_rdf(atoms, center_label, neighbor_label, ...)
```

Output includes:

```
r,
g_raw_norm,
g_smooth,
g_raw2,
g_smooth2,
cumulative_raw,
cumulative_smooth
```

`cumulative_raw` is the coordination number (CN) curve.

---

## 6. ADF Function

```
compute_general_adf(atoms, center_label, neigh1_label, neigh2_label, ...)
```

Supports:

- Distance filtering (r_min, r_max)
- Normalization per central atom
- Gaussian smoothing

---

## 7. Example (RDF)

```python
from rdfadf import compute_general_rdf
from ase.io import read

atoms = read("POSCAR")

r, g_raw, g_smooth, g_raw2, g_smooth2, cn_raw, cn_smooth =     compute_general_rdf(
        atoms,
        center_label="Co",
        neighbor_label="O",
        rcut=6.0,
        bins=200,
        sigma=2.0
    )
```

---

## 8. Example (ADF)

```python
from rdfadf import compute_general_adf
from ase.io import read

atoms = read("POSCAR")

theta, adf_raw, adf_smooth = compute_general_adf(
    atoms,
    center_label="Co",
    neigh1_label="O",
    neigh2_label="O",
    rcut=4.0,
    r_min=1.0,
    r_max=3.0,
    bins=180,
    sigma=2.0
)
```

---

## 9. Notes

- ADF computation may be more expensive than RDF (O(N × n_neighbors²))
- `atoms.pbc = True` is required for proper PBC handling
- Large MD trajectories require frame-by-frame analysis
- Adjust `bins` and `rcut` depending on the system

---

## 10. License

MIT License is recommended.
