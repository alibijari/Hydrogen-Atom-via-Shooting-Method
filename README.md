# Hydrogen Atom: Shooting Method Solver

This repository contains a **complete numerical solution of the radial Schrödinger equation for the hydrogen atom**, using the **shooting method** and **Runge-Kutta integration**. The code is written in Python and well-documented for educational and research purposes.

## Features

- **Numerical solution of the radial Schrödinger equation** for the hydrogen atom.
- **Shooting method** implementation to find bound state (energy eigenvalue) solutions.
- **Runge-Kutta (RK4) method** for stable integration of the wavefunction.
- Plots of **energy eigenvalues** and **radial wavefunctions** for different quantum numbers.
- Full set of **physical constants** and initial conditions defined in SI units.
- **Thorough code comments** and physical background explanations for each section.

## Theoretical Background

The code solves the radial part of the time-independent Schrödinger equation for the hydrogen atom:

\[
\left[ -\frac{\hbar^2}{2m}\frac{d^2}{dr^2} + \frac{l(l+1)\hbar^2}{2mr^2} - \frac{e^2}{4\pi\epsilon_0 r} \right] u(r) = E u(r)
\]

where \(u(r) = r R(r)\) is the scaled radial wavefunction, \(l\) is the orbital angular momentum quantum number, and \(E\) is the energy eigenvalue. The equation is solved for physically acceptable, normalizable solutions.

### Shooting Method

- The **shooting method** is a root-finding algorithm where the differential equation is integrated for a guessed energy. The guess is adjusted until the boundary conditions (vanishing at infinity) are met.
- The method tracks the sign change of the wavefunction at large \(r\) to identify eigenvalues.

### Runge-Kutta Integration

- A **4th-order Runge-Kutta (RK4)** method is used for numerical stability and accuracy in integrating the wavefunction from small \(r\) to large \(r\).

## Usage

1. **Clone the repository** and run the code in Python 3.
2. Required libraries: `numpy`, `matplotlib`, `scipy`
3. **Edit parameters** (such as number of eigenvalues, quantum number \(l\), grid size, etc.) in the script if needed.
4. The code will output:
    - Energy eigenvalues (in eV)
    - Plots of radial wavefunctions for the hydrogen atom

## Example Output

![Example plot of hydrogen atom wavefunctions](example_wavefunction.png)

## File Structure

- `Shooting method.py` : Main code for solving the Schrödinger equation using shooting method and Runge-Kutta integration.
- `Shooting method.pdf` : Educational lesson note and theoretical background.
- `README.md` : Project description and user guide.

## References

- Griffiths, D. J. *Introduction to Quantum Mechanics* (2nd Edition), Pearson, 2004.
- Bransden, B. H. & Joachain, C. J., *Quantum Mechanics*, 2nd Edition, Prentice Hall, 2000.
- Numerical Recipes in C: The Art of Scientific Computing (Ch. 17, Shooting Method)

## License

This project is for academic and educational use.

---

پ
