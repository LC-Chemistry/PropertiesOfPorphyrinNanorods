# PropertiesOfPorphyrinNanorods
Python code that created and diagonalized Hamiltonian matrices and, from that, is able to find the absorption intensity and wavelength corresponding to the transition to the first excited state. This code was developed Spring 2026 by Lillian Cassell.

To run the code you must define the arguments that are listed at the beginning of the code. Some of the arguments have default values and do not need to be defined every time you wish to run the code. All of the string arguments are case sensitive. An example argument to run the code is as such: python SpectraModel.py --N 24 --Shape "linear" > output.log

# SpectralModel — Molecular Aggregate Absorption Simulator

Simulates the absorption spectrum of molecular aggregates by constructing a Hamiltonian matrix, solving for its eigenvalues/eigenvectors, and computing the absorption cross section across a range of wavelengths. Parallel processing is used to accelerate the wavelength sweep.

---

## Requirements

- Python 3.8+
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)

Install dependencies if needed:

```bash
pip install numpy scipy
```

---

## Running the Script

```bash
python SpectralModel_parallell_updated.py --Shape <shape> [options]
```

`--Shape` is required. All other arguments are optional and will fall back to their defaults.

### Minimal example (linear aggregate, all defaults)

```bash
python SpectralModel_parallell_updated.py --Shape linear
```

### Full example with all arguments specified

```bash
python SpectralModel_parallell_updated.py \
  --Shape cylindrical\ vertical-dipole \
  --N 24 \
  --M 6 \
  --Umon 3.0 \
  --Vmon 23000 \
  --Beta -474.3 \
  --Gamma -670.8 \
  --Zeta -474.3 \
  --Xi -474.3 \
  --Workers 4
```

> **Login node tip:** If running on a shared login node, limit core usage to be considerate of other users (e.g. `--Workers 4`).

---

## Arguments

### Aggregate structure

| Argument | Type | Default | Description |
|---|---|---|---|
| `--Shape` | string | *(required)* | Aggregate geometry. See [Aggregate Shapes](#aggregate-shapes) below for valid values. |
| `--N` | int | `24` | Total number of monomers in the aggregate. |
| `--M` | int | `6` | Number of monomers per cross-sectional ring. **Only used for cylindrical shapes.** |

### Physical parameters

| Argument | Type | Default | Units | Description |
|---|---|---|---|---|
| `--Umon` | float | `3.0` | Debye | Transition dipole moment of a single monomer. |
| `--Vmon` | float | `23000.0` | cm⁻¹ | Monomer transition energy. |
| `--Beta` | float | `-474.3` | cm⁻¹ | Nearest-neighbor coupling energy (β). |
| `--Gamma` | float | `-670.8` | cm⁻¹ | Next-nearest-neighbor coupling energy (γ). |
| `--Zeta` | float | `-474.3` | cm⁻¹ | Intra-ring diagonal coupling energy (ζ). Used in cylindrical horizontal-dipole only. |
| `--Xi` | float | `-474.3` | cm⁻¹ | Inter-ring diagonal coupling energy (ξ). Used in cylindrical horizontal-dipole only. |

### Parallelism

| Argument | Type | Default | Description |
|---|---|---|---|
| `--Workers` | int | `None` | Number of CPU cores to use for the parallel absorption spectrum calculation. If not set, falls back to the `NUM_WORKERS` constant at the top of the script, then to all available cores. |

You can also set the worker count directly in the script by editing the constant near the top of the file:

```python
NUM_WORKERS = 8  # or set to None to use all available cores
```

Priority order: `--Workers` flag → `NUM_WORKERS` constant → all available cores.

---

## Aggregate Shapes

Pass one of the following strings to `--Shape`:

| `--Shape` value | Description | Extra args |
|---|---|---|
| `monomer` | Single monomer (no coupling). | — |
| `linear` | 1D chain of monomers with nearest-neighbor β coupling. | — |
| `ring` | Circular chain with periodic boundary (first and last monomer also coupled). | — |
| `double-stranded` | Two parallel strands; alternating γ (intra-strand) and β (inter-strand) coupling. | — |
| `cylindrical vertical-dipole` | Cylindrical aggregate with vertically oriented dipoles; β and γ couplings along and between rings. | `--M` |
| `cylindrical horizontal-dipole` | Cylindrical aggregate with horizontally oriented dipoles; adds ζ/ξ intra-ring couplings. | `--M` |

> **Note:** Shape strings with spaces must be quoted on the command line, e.g. `--Shape "cylindrical vertical-dipole"`.

---

## Output

The script prints the following to the terminal:

1. **Hamiltonian matrix** — the sparse matrix built from the aggregate geometry and coupling parameters.
2. **Eigenvalues and eigenvectors** — from diagonalizing the Hamiltonian.
3. **Transition energy** — ground to first excited state, in nanometers.
4. **Minimum absorption cross section** — the minimum value of *C*_abs over the scanned wavelength range, in arbitrary units.

---

## How Parallelism Works

The absorption cross section is computed independently for each wavelength from the monomer transition energy (~434 nm for the default `--Vmon`) up to 580 nm. Because each wavelength calculation is fully independent, they are distributed across CPU cores using Python's `multiprocessing.Pool.map`. The matrix construction and eigenvalue steps remain sequential.