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

The current dataset contains several structural classes, including:

- layered/chalcogenide hosts,
- oxide tunnel or slab structures,
- hydroxide-type local environments,
- MXene-like structures.

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

The original VASP `POTCAR` files are **not distributed** in this public repository. The full pseudopotential metadata and reconstruction instructions are documented below.

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

The script contains optional manual overrides for:

- oxidation states,
- site radii,

which can be used when automatic chemical assignment fails for a specific structure.

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

# VASP pseudopotentials and POTCAR metadata

## Why the actual POTCAR files are not included

The actual VASP `POTCAR` files used in this project are **not included in this public repository**.

VASP pseudopotential files are distributed as part of the licensed VASP package. To avoid redistribution of licensed PAW data, this repository records only the metadata needed to identify and reconstruct the pseudopotentials used in the calculations.

For reproducibility, the repository provides:

- the exact PAW dataset label;
- the corresponding release date;
- the element ordering through each `POSCAR`;
- isotope-mass modifications through `INCAR`, where relevant;
- the remaining VASP input and selected output files.

Users who wish to reproduce the calculations should obtain the corresponding PAW datasets through their own valid VASP license and reconstruct the required `POTCAR` files locally.

---

## PAW datasets used in this project

The following PAW-PBE datasets and release dates were used across the structures included in this repository:

```text
PAW_PBE H       15Jun2001
PAW_PBE Li      17Jan2003
PAW_PBE C       08Apr2002
PAW_PBE N       08Apr2002
PAW_PBE O       08Apr2002
PAW_PBE F       08Apr2002
PAW_PBE Al      04Jan2001
PAW_PBE Si      05Jan2001
PAW_PBE P       06Sep2000
PAW_PBE S       06Sep2000
PAW_PBE Sc      04Feb2005
PAW_PBE Ti      08Apr2002
PAW_PBE V       08Apr2002
PAW_PBE Mn      06Sep2000
PAW_PBE Fe      06Sep2000
PAW_PBE Co      02Aug2007
PAW_PBE Ni      02Aug2007
PAW_PBE Ga      08Apr2002
PAW_PBE Ge      05Jan2001
PAW_PBE Se      06Sep2000
PAW_PBE Zr_sv   04Jan2005
PAW_PBE Nb_pv   08Apr2002
PAW_PBE Mo      08Apr2002
PAW_PBE Pd      04Jan2005
PAW_PBE In      08Apr2002
PAW_PBE Sn      08Apr2002
PAW_PBE Te      08Apr2002
PAW_PBE Ba_sv   06Sep2000
PAW_PBE La      06Sep2000
PAW_PBE Ce      23Dec2003
PAW_PBE Nd      23Dec2003
PAW_PBE Hf      20Jan2003
PAW_PBE Ta      17Jan2003
PAW_PBE W_sv    04Sep2015
```

The pseudopotential label and release date should both be matched when reproducing the calculations.

---

## POTCAR metadata table

| Element | PAW dataset | Release date |
|---|---|---|
| H | PAW_PBE H | 15Jun2001 |
| Li | PAW_PBE Li | 17Jan2003 |
| C | PAW_PBE C | 08Apr2002 |
| N | PAW_PBE N | 08Apr2002 |
| O | PAW_PBE O | 08Apr2002 |
| F | PAW_PBE F | 08Apr2002 |
| Al | PAW_PBE Al | 04Jan2001 |
| Si | PAW_PBE Si | 05Jan2001 |
| P | PAW_PBE P | 06Sep2000 |
| S | PAW_PBE S | 06Sep2000 |
| Sc | PAW_PBE Sc | 04Feb2005 |
| Ti | PAW_PBE Ti | 08Apr2002 |
| V | PAW_PBE V | 08Apr2002 |
| Mn | PAW_PBE Mn | 06Sep2000 |
| Fe | PAW_PBE Fe | 06Sep2000 |
| Co | PAW_PBE Co | 02Aug2007 |
| Ni | PAW_PBE Ni | 02Aug2007 |
| Ga | PAW_PBE Ga | 08Apr2002 |
| Ge | PAW_PBE Ge | 05Jan2001 |
| Se | PAW_PBE Se | 06Sep2000 |
| Zr | PAW_PBE Zr_sv | 04Jan2005 |
| Nb | PAW_PBE Nb_pv | 08Apr2002 |
| Mo | PAW_PBE Mo | 08Apr2002 |
| Pd | PAW_PBE Pd | 04Jan2005 |
| In | PAW_PBE In | 08Apr2002 |
| Sn | PAW_PBE Sn | 08Apr2002 |
| Te | PAW_PBE Te | 08Apr2002 |
| Ba | PAW_PBE Ba_sv | 06Sep2000 |
| La | PAW_PBE La | 06Sep2000 |
| Ce | PAW_PBE Ce | 23Dec2003 |
| Nd | PAW_PBE Nd | 23Dec2003 |
| Hf | PAW_PBE Hf | 20Jan2003 |
| Ta | PAW_PBE Ta | 17Jan2003 |
| W | PAW_PBE W_sv | 04Sep2015 |


---

## Reconstructing POTCAR files locally

The order of concatenated PAW datasets for an individual VASP calculation must follow the element order in the corresponding `POSCAR`.

For example, if the `POSCAR` contains:

```text
Li Ti P O
```

then the locally reconstructed `POTCAR` should concatenate the corresponding PAW datasets in exactly the same order:

```text
PAW_PBE Li  17Jan2003
PAW_PBE Ti  08Apr2002
PAW_PBE P   06Sep2000
PAW_PBE O   08Apr2002
```

Similarly, for a structure with:

```text
Li La Ti O
```

the required sequence is:

```text
PAW_PBE Li  17Jan2003
PAW_PBE La  06Sep2000
PAW_PBE Ti  08Apr2002
PAW_PBE O   08Apr2002
```

The general reconstruction procedure is:

```text
POSCAR element order
        ↓
match each element to the PAW dataset listed above
        ↓
concatenate the corresponding licensed POTCAR components
        ↓
place the reconstructed POTCAR in the calculation directory
```

The reconstructed local `POTCAR` should **not** be committed back to the public repository.

---

## Li isotope treatment

The \(^{6}\mathrm{Li}\) and \(^{7}\mathrm{Li}\) calculations use the same Li electronic PAW dataset:

```text
PAW_PBE Li  17Jan2003
```

The isotope effect is introduced through the **nuclear mass**, rather than through a different electronic pseudopotential.

For \(^{6}\mathrm{Li}\) calculations, the Li atomic mass is explicitly changed in the corresponding `INCAR` using `POMASS` where applicable.

For \(^{7}\mathrm{Li}\) calculations, the corresponding reference/default Li mass is used unless explicitly overwritten in the calculation directory.

Therefore:

```text
Electronic PAW potential: unchanged
Nuclear isotope mass:     changed
```

This is the intended setup for comparing isotope-dependent vibrational frequencies and calculating

\[
1000\ln\beta^{7/6}.
\]

---

## Why `_sv` and `_pv` labels matter

Some elements in this project use PAW datasets with explicit semicore treatment, including:

```text
Zr_sv
Nb_pv
Ba_sv
W_sv
```

These labels are part of the pseudopotential definition and should not be replaced silently by the corresponding default potentials.

For reproducibility, both the **dataset label** and the **release date** should be matched.

For example:

```text
PAW_PBE Zr_sv  04Jan2005
```

should not automatically be replaced by a generic `Zr` PAW dataset.

The same applies to:

```text
PAW_PBE Nb_pv  08Apr2002
PAW_PBE Ba_sv  06Sep2000
PAW_PBE W_sv   04Sep2015
```

---

## Preventing accidental POTCAR uploads

The repository `.gitignore` should contain rules that exclude actual POTCAR files:

```gitignore
POTCAR
POTCAR.*
**/POTCAR
**/POTCAR.*
```

Before pushing new calculation directories to GitHub, it is recommended to verify that no actual POTCAR file is being tracked.

On Linux/macOS:

```bash
git ls-files | grep POTCAR
```

On Windows PowerShell:

```powershell
git ls-files | Select-String "POTCAR"
```

A properly cleaned public repository should list only metadata-related files such as:

```text
POTCAR_README.md
POTCAR.spec
```

and should not contain the actual licensed `POTCAR` content.

---

## Recommended POTCAR.spec format

For individual calculation directories, an optional `POTCAR.spec` file can be used to document the required pseudopotential sequence without redistributing the pseudopotential itself.

Example:

```text
# POTCAR specification
# Order follows POSCAR

PAW_PBE Li  17Jan2003
PAW_PBE Ti  08Apr2002
PAW_PBE P   06Sep2000
PAW_PBE O   08Apr2002
```

This allows a licensed VASP user to reconstruct the exact `POTCAR` locally.

---

## Reproducibility notes

When reproducing the calculations, the following settings should be kept consistent:

- crystal structure;
- element ordering;
- exchange-correlation functional;
- PAW potential label;
- PAW release date;
- plane-wave cutoff;
- k-point sampling;
- electronic convergence criteria;
- ionic convergence criteria;
- Hubbard \(U\), if used;
- magnetic initialization, where relevant;
- isotope mass;
- vibrational calculation settings;
- treatment of imaginary modes;
- treatment of near-zero-frequency acoustic modes.

For isotope calculations, the \(^{6}\mathrm{Li}\) and \(^{7}\mathrm{Li}\) calculations should differ only in isotope mass unless a specific methodological test explicitly requires otherwise.

---

## Data usage

The repository is intended to provide the computational data and scripts needed to inspect and reproduce the main analysis pipeline:

```text
Crystal structure
        ↓
VASP geometry optimization / vibrational calculations
        ↓
6Li / 7Li vibrational mode extraction and matching
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

