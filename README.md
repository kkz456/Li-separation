# Li Isotope Separation: Structures, Vibrational Data, and Analysis Workflows

This repository contains the structural data, vibrational calculations, representative VASP inputs/outputs, and analysis scripts used in our computational study of **lithium isotope separation in rigid and flexible host materials**.

The workflow combines two complementary descriptors:

- **Equilibrium isotope fractionation**, characterized by the reduced partition function ratio, \(\beta^{7/6}\), calculated from vibrational frequencies.
- **Lithium transport accessibility**, characterized by the **CAVD bottleneck radius**.

Together, these quantities provide a framework for screening candidate materials by separating **thermodynamic isotope selectivity** from **geometric/transport accessibility**.

---

## Repository structure

```text
Li-separation/
├── code_ of_beta_cavd/
│   ├── Beta.py
│   └── cavd_bottleneck.py
│
├── Rigid host/
│   ├── 晶体结构/
│   ├── 振动数据/
│   │   ├── FREQ_MASS=6/
│   │   └── FREQ_MASS=7/
│   └── vesta文件/
│
├── Flexible host/
│   ├── structure/
│   └── 晶体结构/
│       ├── freq_mass=6/
│       ├── freq_mass=7/
│       └── structure/
│
├── Li_H2O_4_8_12/
│   ├── 结构/
│   └── 振动信息/
│       ├── freq=6/
│       └── freq=7/
│
└── VASP/
    ├── Rigid/
    ├── Flexible/
    ├── .gitignore
    └── POTCAR_README.md
```

### `Rigid host/`

Contains crystal structures and isotope-dependent vibrational data for the rigid-host screening set.

The vibrational calculations are separated into:

- `FREQ_MASS=6/`: calculations using the \(^{6}\mathrm{Li}\) isotope mass.
- `FREQ_MASS=7/`: calculations using the \(^{7}\mathrm{Li}\) isotope mass.

Representative structures used for visualization are also provided in VASP/VESTA formats.

### `Flexible host/`

Contains candidate flexible-host structures and corresponding isotope-dependent vibrational data.

The current dataset contains several structural classes, including layered/chalcogenide hosts, oxide tunnel or slab structures, hydroxide-type local environments, and MXene-like structures.

For selected systems, multiple Li sites are included to investigate how the **local coordination environment** modifies the isotope fractionation factor.

### `Li_H2O_4_8_12/`

Contains structures and vibrational data for Li-containing water clusters with 4, 8, and 12 water molecules.

These calculations are included as molecular/reference systems for comparison with solid-state hosts.

### `VASP/`

Contains representative VASP calculation directories, including geometry optimization and isotope-dependent vibrational calculations.

Typical files include:

```text
INCAR
KPOINTS
POSCAR
OUTCAR
```

The original VASP `POTCAR` files are **not distributed** in this repository. See [`VASP/POTCAR_README.md`](VASP/POTCAR_README.md) for details.

---

## Isotope fractionation calculation

The main script for calculating the Li isotope beta factor is:

```text
code_ of_beta_cavd/Beta.py
```

The script evaluates the harmonic reduced partition function ratio

\[
\beta^{7/6}
\]

from paired \(^{6}\mathrm{Li}\) and \(^{7}\mathrm{Li}\) vibrational modes.

The reported quantity used for screening is

\[
1000\ln\beta^{7/6}.
\]

For unit cells containing different numbers of Li atoms, the repository also uses a **per-Li normalized quantity** for cross-material comparison.

### Treatment of vibrational modes

The current implementation:

1. pairs the \(^{6}\mathrm{Li}\) and \(^{7}\mathrm{Li}\) vibrational modes using their eigenvectors;
2. excludes modes identified as imaginary (`fi` or `f/i`);
3. retains only mode pairs that are real and positive in both isotope calculations;
4. removes very-low-frequency acoustic/zero modes below a configurable cutoff;
5. evaluates the beta factor under the harmonic approximation.

If SciPy is available, mode matching uses the Hungarian algorithm; otherwise, a greedy matching scheme is used.

### Python dependencies

```bash
pip install numpy openpyxl
```

Optional:

```bash
pip install scipy
```

### Example

For a directory organized as

```text
calculation/
├── FREQ_MASS=6/
├── FREQ_MASS=7/
├── crystal_structure/
└── bottleneck.xlsx
```

run:

```bash
python Beta.py \
  --root . \
  --mass6-dir "FREQ_MASS=6" \
  --mass7-dir "FREQ_MASS=7" \
  --crystal-dir "crystal_structure" \
  --bottleneck "bottleneck.xlsx" \
  --temperature 300
```

The default acoustic-mode cutoff is:

```text
5 cm^-1
```

and can be changed with:

```bash
--acoustic-cutoff-cm1
```

> Note: \(\beta\) is a single-phase equilibrium isotope fractionation descriptor.  
> An isotope fractionation factor \(\alpha\) requires a second reference phase and is therefore not identified with the single-phase beta factor in this workflow.

---

## CAVD bottleneck calculation

The bottleneck analysis script is:

```text
code_ of_beta_cavd/cavd_bottleneck.py
```

It uses CAVD together with pymatgen-based structure processing to characterize the geometric accessibility of Li migration pathways.

The workflow includes:

1. reading optimized structures;
2. assigning oxidation states;
3. constructing CAVD-accessible void/channel information;
4. extracting Li-related bottleneck descriptors;
5. writing material-level and site-level results to CSV/Excel files.

The script contains optional manual overrides for oxidation states and site radii, which can be used when automatic chemical assignment fails for a specific structure.

Important dependencies include:

```text
cavd
numpy
pandas
pymatgen
```

Because CAVD installation can be environment-dependent, users should reproduce the analysis in a compatible CAVD environment.

---

## Representative systems

The repository contains representative calculations for materials used to illustrate different regions of the screening space, including structures associated with:

- Li-containing titanate/oxide hosts,
- LLTO-type structures,
- LATP-related structures,
- NiHO2 local Li environments,
- V2O5-related flexible hosts.

The exact structures used in the calculations are provided in the corresponding structure directories.

---

## VASP pseudopotentials

VASP `POTCAR` files are intentionally excluded because the PAW potential files are distributed under the VASP license and should not be redistributed publicly.

For reproducibility, the VASP calculation directories use metadata/specification files where appropriate, and isotope masses are controlled through the corresponding calculation inputs.

For the Li isotope phonon calculations:

- \(^{6}\mathrm{Li}\) calculations use the Li isotope mass specified in the corresponding `INCAR` where applicable;
- \(^{7}\mathrm{Li}\) calculations use the corresponding reference/default Li mass unless otherwise specified.

Users with a valid VASP license should reconstruct the required `POTCAR` files locally.

See:

```text
VASP/POTCAR_README.md
```

for additional information.

---

## Reproducibility notes

When reproducing the calculations, please keep the following points consistent:

- crystal structure;
- exchange-correlation functional;
- PAW potential choice;
- plane-wave cutoff;
- k-point sampling;
- electronic and ionic convergence criteria;
- Hubbard \(U\), if used;
- magnetic initialization, where relevant;
- isotope mass;
- vibrational calculation settings;
- treatment of imaginary and near-zero-frequency modes.

For isotope effects, the \(^{6}\mathrm{Li}\) and \(^{7}\mathrm{Li}\) calculations should differ only in the isotope mass unless a specific test explicitly requires otherwise.

---

## Data usage

The repository is intended to provide the computational data and scripts needed to inspect and reproduce the main analysis pipeline:

```text
Crystal structure
        ↓
VASP vibrational calculations
        ↓
6Li / 7Li mode matching
        ↓
1000 ln(beta7/6)
        ↓
CAVD bottleneck calculation
        ↓
Combined thermodynamic–transport screening
```

The raw and processed data are organized by host type to facilitate both reproduction of individual calculations and comparison across materials.

---

## Citation

If you use this repository, please cite the associated publication.

```text
[Publication information will be added after acceptance/publication.]
```

---

## License and third-party software

The scripts in this repository are provided for academic reproducibility.

VASP itself, VASP source code, and VASP PAW/POTCAR datasets are **not** included and remain subject to the VASP license.

Other third-party software, including CAVD and pymatgen, is subject to its own respective license.

---

## Contact

For questions regarding the calculations or data organization, please open an issue in this repository or contact the authors of the associated publication.
