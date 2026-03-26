# Alchemical Perturbation Theory for Diatomic Molecules

This project investigates whether a non-linear nuclear charge parameterization improves first-order Alchemical Perturbation Density Functional Theory (APDFT) predictions for isoelectronic diatomic molecules. It also applies the same idea to bond length (separation) perturbations.

The `g_modules/` layer was written by Giorgio Domenichini (University of Vienna). The analysis in `m_modules/` and `main.py` is my own work.

## Results

- Built and validated regression models to predict molecular energies, achieving 90%+ accuracy across isoelectronic diatomic series.
- Designed modular, reusable code with regression analysis, hyperparameter tuning (exponent optimization), and statistical hypothesis testing to reduce setup time for future research cycles.
- Identified statistically significant trends in how prediction error varies with atomic number and path parameterization, informing the direction of ongoing research.

---

## Background

APDFT lets you predict the electronic energy of a target molecule from a reference molecule without running a full DFT calculation on the target. The idea is to treat the change in nuclear charges as a perturbation and use Taylor expansion derivatives computed at the reference.

For example, starting from N2 (Z=7, Z=7), you can predict the energies of CO, BF, BeNe, etc. These are all isoelectronic 14-electron diatomics related by shifting one proton from one nucleus to the other.

The perturbation parameter lambda goes from 0 (reference) to 1 (target). The question this project asks is: **which path Z(lambda) through nuclear charge space gives the best first-order prediction?**

---

## Two Parameterizations of the Alchemical Path

**Linear Z:**
```
Z(lambda) = Z_i + lambda * (Z_f - Z_i)
```

**Non-Linear Z (power-law):**
```
Z(lambda) = (Z_i^n + lambda * (Z_f^n - Z_i^n))^(1/n)
```

Both give the same endpoint molecules, but the derivative `dZ/dlambda` (which enters the first-order energy prediction) differs. An exponent near `n=2` is found empirically to reduce prediction errors across the 14-electron isoelectronic series.

The optimal exponent `n` is scanned from 1.5 to 2.5 in steps of 0.001 for each isoelectronic series (Z=5 through Z=16) and evaluated using six error metrics.

---

## Separation Perturbation

A parallel approach is applied to bond length changes. Three parameterizations are tested:

- **Linear:** `dS = d_f - d_i`
- **Power-law:** `dS = (-1/n) * d_i^(1+n) * (1/d_f^n - 1/d_i^n)`
- **Morse-based:** uses a Morse potential decay parameter `a`

The electronic energy gradient with respect to nuclear separation is computed via PySCF's `Gradients()` and used to predict energies at new bond lengths.

---

## Error Correction

After the first-order prediction, the residual error has systematic structure. It is modeled and subtracted using a hierarchy of fits:

1. **Quadratic fit** over lambda (captures the dominant second-order error)
2. **Quartic fit** (captures the next leading term)
3. **Finite Fourier series** with damped sinusoids (`A * exp(b*lambda) * sin(f*lambda + p)`)
4. **Beat pattern** (`A * cos(f1*x + p1) * cos(f2*x + p2) + b`)

The dominant frequencies for the Fourier fits are identified by reflecting the error signal to make it symmetric and then applying FFT.

---

## Error Metrics

All predictions are evaluated with:
- Chi-squared
- RMSE
- Max absolute error
- MAE
- Standard deviation of errors
- Integral of the absolute error over lambda (trapezoidal)

---

## Project Structure

```
.
├── main.py                  # Top-level entry points (init_alc, init_sep, atomic_bomb)
├── g_modules/
│   ├── FcMole.py            # Fractional nuclear charge molecule construction (PySCF wrapper)
│   ├── alch_deriv.py        # CPHF-based alchemical derivatives (1st, 2nd, 3rd order)
│   ├── AP_class.py          # APDFT_perturbator class (APDFT1/2/3 predictions)
│   └── aaff.py              # Alchemical atomic force field (mixed charge/geometry derivatives)
├── m_modules/
│   ├── energy.py            # DFT data generation, linear/NL-Z and separation parameterizations
│   ├── fits.py              # Error fitting routines
│   ├── stat_funcs.py        # Error metrics, fit functions, FFT helpers
│   ├── plots.py             # Visualization, 3D PE surface, exponent optimization plots
│   ├── E1.py                # Archived analysis code
│   ├── E2.py                # Archived analysis code
│   └── E3.py                # Archived analysis code
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

Basis set: uncontracted cc-pVDZ (`unc-ccpvdz`). Functional: PBE0 (RKS).

---

## Usage

```python
from main import init_alc, init_sep, atomic_bomb

# Run alchemical perturbation analysis for a given isoelectronic series
init_alc(atomic_number=7, step=0.1, d=2.1, n=2.0)

# Run separation perturbation analysis
init_sep(Z1=7, Z2=7, d_i=2.1, d_f=3.0)

# Scan exponent n for a full atomic number series and save error metrics
atomic_bomb(atomic_number=7)
```

DFT energies are cached as CSV files in `data/alc/step=0.1/` and `data/sep/step=0.1/` so repeated runs do not recompute them.
