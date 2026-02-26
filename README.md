# Thin Airfoil Coding Assignment

This project contains a detailed aerodynamic analysis of the Clark-Y Airfoil leveraging classical Thin Airfoil Theory (TAT) extended to ground-effect aerodynamics.

## Files
- `thin_airfoil_report.md`: A comprehensive report of the Clark-Y airfoil analysis, including the methodology, formulations, computational results, and verifications.
- `thin_airfoil_clarky.jl`: The companion Julia script executing adaptive Gauss-Kronrod quadrature integrals to uncover the image-interaction matrix and bound circulation distributions.
- `assets/`: Contains generated plots of Lift Coefficient ($C_L$), Pressure Coefficient Distribution ($\Delta C_p$), and analytical cross-verifications.

## Usage
To execute the Julia script:
```bash
julia thin_airfoil_clarky.jl
```
It utilizes packages like `Plots` (with headless GR rendering backend), `QuadGK`, and `LaTeXStrings`. The outputs are saved into `.png` plots within your active directory.
