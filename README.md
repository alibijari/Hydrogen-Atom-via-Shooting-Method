# 🪐 Hydrogen Atom: Shooting Method Solver

> **Numerically solving the radial Schrödinger equation for the hydrogen atom using Python, the shooting method, and Runge-Kutta integration.**  
> _Educational, well-commented, and ready for research or teaching!_

---

## 🚀 Features at a Glance

- 🎯 **Full numerical solution** for the hydrogen atom radial Schrödinger equation
- 🏹 **Shooting method:** Finds quantum energy eigenvalues (bound states)
- 🧮 **Runge-Kutta (RK4) integration** for accurate and stable propagation
- 📈 **Automatic plotting** of all eigenfunctions and energy levels
- 🧪 **Physical constants** and all parameters in SI units—fully transparent
- 💡 **Detailed code comments & physics explanations** in every section

---

## 🧑‍🔬 Theoretical Background

The script solves the **radial time-independent Schrödinger equation**:

\[
\left[ -\frac{\hbar^2}{2m}\frac{d^2}{dr^2} + \frac{l(l+1)\hbar^2}{2mr^2} - \frac{e^2}{4\pi\epsilon_0 r} \right] u(r) = E u(r)
\]

- \( u(r) = r R(r) \) is the reduced radial wavefunction.
- \( l \) is the orbital angular momentum quantum number.
- \( E \) is the energy eigenvalue.

**Goal:** Find physically normalizable eigenfunctions and quantized energy levels, just as in real hydrogen!

---

### 🏹 Shooting Method

- **How it works:**  
  - Guess the energy $E$.
  - Numerically integrate the equation from small $r$ to large $r$.
  - Adjust $E$ until the solution meets boundary conditions (vanishing at infinity).

- **Root-finding:**  
  Detects sign changes in the wavefunction at the end of integration to pinpoint true eigenvalues.

---

### 🔄 Runge-Kutta Integration (RK4)

- Ensures **accuracy and stability** for integrating ODEs (wavefunction propagation).
- Handles both oscillatory and decaying behaviors as required for quantum bound states.

---

## ⚡️ Usage

1. **Clone or download this repo**.
2. Make sure you have Python 3 and these libraries:
   ```bash
   pip install numpy matplotlib scipy
