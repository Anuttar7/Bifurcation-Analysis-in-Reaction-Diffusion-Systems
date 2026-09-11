# Bifurcation Analysis in a Predator–Prey Model

> **MM6103 — Individual Term Project | IIT Bombay**

An implementation of the numerical methodology presented in *Numerical simulations of multicomponent ecological models with adaptive methods* by K. M. Owolabi and K. C. Patidar.

The project studies a predator–prey model with reaction–diffusion characteristics and uses a finite-difference approach based on **Exponential Time Differencing (ETD)**. The complete implementation is written in **MATLAB**.

---

## 📌 Project Overview

The project focuses on the numerical simulation and bifurcation analysis of a predator–prey reaction–diffusion system.

### Key aspects

- 🐇 **Model:** Predator–prey reaction–diffusion system
- 🔢 **Spatial discretisation:** Fourth-order centred finite differences
- ⏱️ **Time integration:** Exponential Time Differencing
- 🚀 **ETD scheme:** ETDRK4
- 🧮 **Matrix approximation:** Carathéodory–Fejér (CF) approximations
- 💻 **Implementation:** MATLAB

---

## 📁 Repository Structure

```text
.
├── Presentation.pdf
└── Code/
    ├── RUN_THIS.m
    └── ...
```

### `Presentation.pdf`

Contains the project report, including:

- Motivation
- Solution strategy
- Numerical methodology
- Simulation results

### `Code/`

Contains the MATLAB source files and dependencies required to reproduce the numerical simulations.

**`RUN_THIS.m` is the main entry point.** Run this file to generate the results; the remaining files serve as dependencies.

---

## 🧩 Model Used

The primary model is a **non-dimensionalised reaction–diffusion system**.

```math
\frac{\partial u}{\partial t} = \frac{\partial^2 u}{\partial x^2} + u(1-u) - \frac{\mu uv}{u + \phi} = f(u, v) \tag{1}
```
```math
\frac{\partial v}{\partial t} = D \frac{\partial^2 v}{\partial x^2} + \psi v - \frac{\psi v^2}{u} = g(u, v) \tag{2}
```

The variables \(u\) and \(v\) represent the **non-dimensional prey and predator densities**, respectively, with the densities non-dimensionalised using the carrying capacity.

The model contains three bifurcation parameters:

$$
\mu,\quad \psi,\quad \phi
$$

along with a constant diffusion-related parameter \(D\).

---

## ⚙️ Solution Strategy

The numerical solution procedure consists of the following steps:

1. **Spatial discretisation**  
   Represent the second-order spatial derivatives using a **fourth-order centred finite-difference scheme**.

2. **Stiffness treatment**  
   Integrate the linear component using an **exponential integrator**, reducing the stiffness associated with the linear terms.

3. **Time integration**  
   Approximate the resulting integral using a **fourth-order Runge–Kutta method**, giving the **ETDRK4** scheme.

4. **Non-diagonal matrix treatment**  
   Use **Carathéodory–Fejér (CF) approximations** to enable the ETDRK4 method to be applied to non-diagonal matrices.

---

## 📚 References

1. K. M. Owolabi and K. C. Patidar, “Numerical simulations of multicomponent ecological models with adaptive methods,” *Theoretical Biology and Medical Modelling*, vol. 13, no. 1, 2016.  
   DOI: `10.1186/s12976-016-0027-4`

2. S. M. Cox and P. C. Matthews, “Exponential Time Differencing for Stiff Systems,” *Journal of Computational Physics*, vol. 176, pp. 430–455, 2002.  
   DOI: `10.1006/jcph.2002.6995`

3. T. Schmelzer and L. Trefethen, “Evaluating matrix functions for exponential integrators via Caratheodory-Fejer approximation and contour integrals,” *Electronic Transactions on Numerical Analysis*, vol. 29, pp. 1–18, 2007.  
   [Publication link](https://www.researchgate.net/publication/251168138)
