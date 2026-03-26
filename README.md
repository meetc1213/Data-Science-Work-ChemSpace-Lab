# Molecular Energy Prediction with Alchemical Perturbation Theory

A data science pipeline for predicting molecular energies across chemical space without running full DFT calculations for every molecule. Built on Alchemical Perturbation DFT (APDFT), the project frames energy prediction as a regression problem and uses hyperparameter tuning and residual error modeling to push prediction accuracy above 90%.

The `g_modules/` layer was written by Giorgio Domenichini (University of Vienna). The modeling, analysis, and optimization in `m_modules/` and `main.py` is my own work.

## Results

- Built and validated regression models to predict molecular energies, achieving 90%+ accuracy across isoelectronic diatomic series.
- Tuned model hyperparameters (power-law exponent `n`) via grid search across Z=5 to Z=16, evaluated on six loss metrics (RMSE, MAE, chi-squared, max error, std, integral error).
- Designed modular, reusable code with regression analysis, hyperparameter tuning, and hypothesis testing to reduce setup time for future research cycles.
- Identified statistically significant trends in how prediction error varies with atomic number and path parameterization, informing the direction of ongoing research.

---

## Problem

Running quantum chemistry calculations (DFT) for every molecule in a large chemical space is expensive. APDFT sidesteps this by treating a new molecule as a perturbation of a reference molecule: compute derivatives at the reference once, then predict energies elsewhere analytically.

This project asks: **which parameterization of the perturbation path minimizes prediction error?**

Two directions are studied independently:
1. **Alchemical:** changing the nuclear charges (atomic numbers) of both atoms
2. **Separation:** changing the bond length

---

## Model

### Alchemical Direction

The perturbation parameter lambda goes from 0 (reference molecule) to 1 (target molecule). The nuclear charge along the path is either:

**Linear (baseline model):**
```
Z(lambda) = Z_i + lambda * (Z_f - Z_i)
```

**Non-linear power-law (tuned model):**
```
Z(lambda) = (Z_i^n + lambda * (Z_f^n - Z_i^n))^(1/n)
```

The derivative `dZ/dlambda` feeds into the first-order energy prediction. The exponent `n` is the key hyperparameter. For example, starting from N2, both paths reach CO, BF, or BeNe at lambda=1, but the non-linear path changes how the gradient is weighted.

### Separation Direction

Three parameterizations of the bond length change are compared:

- **Linear:** `dS = d_f - d_i`
- **Power-law:** `dS = (-1/n) * d_i^(1+n) * (1/d_f^n - 1/d_i^n)`
- **Morse-based:** uses a Morse potential decay parameter `a`

---

## Hyperparameter Tuning

The exponent `n` is optimized via grid search over `[1.5, 2.5]` in steps of 0.001 for each isoelectronic series (Z=5 through Z=16). All six loss metrics are recorded at each value of `n`, and the minimizing exponent is identified per metric per series. Results are visualized as 3D scatter plots of optimal `n` vs. atomic number.

---

## Residual Modeling

After the first-order prediction, the residual error has systematic structure. A hierarchy of regression models is fit to the residuals and subtracted:

1. **Quadratic regression** over lambda (dominant second-order term)
2. **Quartic regression** (next leading term)
3. **Damped sinusoidal fit** (`A * exp(b*lambda) * sin(f*lambda + p)`) with up to 5 terms
4. **Beat pattern fit** (`A * cos(f1*x + p1) * cos(f2*x + p2) + b`)

Dominant frequencies for step 3 are identified by applying FFT to the (reflected) residual signal and fitting a damped harmonic to the spectrum.

---

## Loss Metrics

| Metric | Description |
|---|---|
| RMSE | Root mean squared error |
| MAE | Mean absolute error |
| Chi-squared | Sum of squared errors normalized by target |
| Max error | Worst-case absolute error |
| Std | Standard deviation of errors |
| Integral error | Trapezoidal integral of absolute error over lambda |

---

## Project Structure

```
.
├── main.py                  # Entry points: init_alc, init_sep, atomic_bomb (hyperparameter scan)
├── g_modules/
│   ├── FcMole.py            # Fractional nuclear charge molecule construction (PySCF wrapper)
│   ├── alch_deriv.py        # CPHF-based energy derivatives (1st, 2nd, 3rd order)
│   ├── AP_class.py          # APDFT_perturbator: generates predictions from derivatives
│   └── aaff.py              # Mixed charge/geometry derivative (alchemical force field)
├── m_modules/
│   ├── energy.py            # Data generation, linear/NL-Z and separation parameterizations
│   ├── fits.py              # Residual modeling (quadratic, quartic, Fourier, beat fits)
│   ├── stat_funcs.py        # Loss metrics, fit functions, FFT helpers
│   ├── plots.py             # Visualization, 3D PE surface, hyperparameter optimization plots
│   ├── E1.py                # Archived analysis
│   ├── E2.py                # Archived analysis
│   └── E3.py                # Archived analysis
├── data/
│   ├── alc/step=0.1/        # Cached DFT energies for alchemical scans
│   └── sep/step=0.1/        # Cached DFT energies for separation scans
└── data/figs/               # Output figures
```

---

## Dependencies

- [PySCF](https://pyscf.org/) (DFT calculations, CPHF)
- NumPy, SciPy (numerics, curve fitting, FFT)
- Matplotlib (plotting)
- pyperclip (LaTeX table output)

Basis set: uncontracted cc-pVDZ. Functional: PBE0 (RKS).

---

## Usage

```python
from main import init_alc, init_sep, atomic_bomb

# Predict energies for the 14-electron isoelectronic series at bond length d=2.1 Bohr
init_alc(atomic_number=7, step=0.1, d=2.1, n=2.0)

# Predict energies across bond lengths for N2
init_sep(Z1=7, Z2=7, d_i=2.1, d_f=3.0)

# Grid search over n for a full isoelectronic series, save all loss metrics
atomic_bomb(atomic_number=7)
```

DFT energies are cached as CSV files so repeated runs do not recompute them.
