# Least-Squares Fitting of Experimental Spectra

Fortran and Python implementation of least-squares Lorentzian fitting for five experimental spectra measured at 1, 3, 6, 10, and 15 atm. Compares unweighted and intensity-weighted fits and examines pressure-dependent line shifts and broadening.

## Files
- `main.f90` — Fortran fitting and output generation
- `read_spectra.py` — plotting and linear regression
- `NumericalSimulations_Practical1_LeastSquareFits.pdf` — project report

## Run
Requires `gfortran`, Python, NumPy, and Matplotlib. Place the input data (`spectre_1.txt` through `spectre_5.txt`) beside the code, then run:

```bash
gfortran -o fit main.f90
./fit
python3 read_spectra.py
```

**Note:** The input spectra are not included in this repository.