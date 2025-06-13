
# Lotka-Volterra System – Numerical Methods and Parameter Estimation

## Description

This MATLAB project simulates and analyzes the Lotka-Volterra predator-prey model using several numerical integration techniques and real-world data fitting. The system of differential equations is solved using:

- `ode45` (MATLAB built-in solver)
- Explicit Euler Method
- Implicit Euler Method
- Adams-Bashforth Method (2nd order)
- Two-step Trapezoidal Method

The final part estimates model parameters using experimental data from the file `MNUB_24L_P3_dane17.csv`.

---

## Files

- `main.m` – Contains the complete numerical analysis and plotting
- `MNUB_24L_P1_dane17.mat` – Experimental data (optional analysis)
- `MNUB_24L_P2_dane17.mat` – Experimental data (optional analysis)
- `MNUB_24L_P3_dane17.csv` – Real data used for parameter estimation
- `GUI.py` and `main.py` – Optional Python GUI interface (not required for MATLAB execution)

---

## Tasks Implemented

### Task 1 – Solving with `ode45`

Numerically solves the Lotka-Volterra system using high-precision ODE solver with default tolerances. The resulting solution is used as a reference for later error estimation.

---

### Task 2 – Custom Integration Methods

Implemented four methods:
- **(a) Explicit Euler** – first-order approximation
- **(b) Implicit Euler** – solved with `fsolve` for each timestep
- **(c) Adams-Bashforth (2nd order)** – multi-step method
- **(d) Two-step Trapezoidal Method** – predictor-corrector scheme

Each method outputs its own `x(t)` and `y(t)` values over the range [0, 1] with step size `h = 0.005`.

---

### Task 3 – Error Estimation

Compares each numerical method to the reference solution (`ode45`) and computes the relative aggregated error using:

```
error = sqrt(sum((y_interp - y_ref)^2)) / sqrt(sum(y_ref^2))
```

---

### Task 4 – Error vs Step Size Plot

Evaluates the dependency of the numerical error on step size `h` from 1/10000 to 1/100, only for valid steps (integer number of timesteps). Plots log-log error curves for all four methods.

---

### Task 5 – Parameter Estimation (Data Fitting)

Fits the model parameters (p1, p2, p3, p4) to real data using `fminsearch`, minimizing squared error between model output and experimental data. The estimated parameters are printed after optimization.

---

## Requirements

- MATLAB R2021 or newer
- Optimization Toolbox (`fminsearch`)
- Curve fitting optional, depending on tasks

---

## Output

The program generates:
- Multiple figures showing `x(t)` and `y(t)` for each method
- Error plots (Task 4)
- Optimized parameter values (Task 5)

---

## Notes

- Task 5 is self-contained and will automatically load `MNUB_24L_P3_dane17.csv`.
- Tasks 1–4 are entirely numerical and do not require external files.

---

## Author

Wojciech Malesiński – 5th semester, Biomedical Engineering, Warsaw University of Technology
