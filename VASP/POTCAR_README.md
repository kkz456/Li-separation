POTCAR files are intentionally omitted from this archive.

Reason
------
VASP PAW potential files are distributed under the VASP license and should not be redistributed publicly.

Reproducibility
---------------
Each calculation directory that originally contained a POTCAR now contains a POTCAR.spec file recording the PAW_PBE potential labels and release dates used, in the same order as the original POTCAR datasets.

For Li isotope phonon calculations:
- freq=6 directories override the Li mass through POMASS in INCAR (6.015).
- freq=7 directories use the default Li mass from the referenced Li PAW dataset unless otherwise specified in INCAR.

Users with a valid VASP license should reconstruct POTCAR locally from their licensed PAW dataset using the listed labels.
