import matplotlib.pyplot as plt
import numpy as np
from scipy.special import genlaguerre, factorial

# ============================================
# 1. PHYSICAL CONSTANTS AND PARAMETERS
# ============================================
a_0 = 5.29e-11         # Bohr radius (meters)
a = 1.00054 * a_0      # (Approximate) reduced mass radius in the hydrogen atom (meters)
e = 1.602176634e-19    # Elementary charge (Coulombs)
m = 9.10938356e-31     # Electron mass (kg)
h_bar = 1.0545718e-34  # Reduced Planck constant (J·s)
E_0 = 13.606           # Ground state energy (eV, for reference)
k_e = 8.98755e9        # Coulomb constant (N·m²·C⁻²)
z = 1                  # Atomic number for hydrogen atom
psi_min = 0.0          # Initial value of wavefunction at r_min
dpsi_dr_guess = 1.0    # Initial guess for dψ/dr at r_min
num_points = 1000      # Number of points in spatial grid
num_eigenvalues = 10   # Number of energy levels to find
l = 0                  # Orbital angular momentum quantum number

# ============================================
# 2. RADIAL SCHRÖDINGER EQUATION (Physics)
# ============================================
def schrodinger_equation(r, y, E, l):
    """
    Returns the derivatives for the coupled first-order ODEs corresponding
    to the radial Schrödinger equation:
        dψ/dr = y[1]
        d²ψ/dr² = [l(l+1)/r² - 2m/(ħ²)(V(r) - E)] * ψ
    where:
        - ψ is the radial wave function
        - l is the angular quantum number
        - V(r) is the Coulomb potential
    The function normalizes ψ to reduce numeric overflow.
    """
    psi, dpsi_dr = y
    V = - (k_e * z * e**2) / r     # Hydrogen Coulomb potential, -Ze²/r

    # Normalize wavefunction for numeric stability
    max_psi = np.max(np.abs(psi))
    if max_psi != 0:
        psi_normalized = psi / max_psi
    else:
        psi_normalized = psi

    d2psi_dr2 = (l * (l + 1)) / r ** 2 - 2 * (m / h_bar ** 2) * (V - E) * psi_normalized
    return np.array([dpsi_dr, d2psi_dr2])

# ============================================
# 3. RUNGE-KUTTA 4TH-ORDER INTEGRATOR (Programming)
# ============================================
def runge_kutta_step(r, y, h, E, l):
    """
    Advances the solution of the ODE by one step using
    the classical fourth-order Runge-Kutta method.
    """
    k1 = h * schrodinger_equation(r, y, E, l)
    k2 = h * schrodinger_equation(r + h / 2, y + k1 / 2, E, l)
    k3 = h * schrodinger_equation(r + h / 2, y + k2 / 2, E, l)
    k4 = h * schrodinger_equation(r + h, y + k3, E, l)
    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6

# ============================================
# 4. SCHRÖDINGER SOLVER WITH SHOOTING METHOD (Physics & Programming)
# ============================================
def solve_schrodinger(E, l, num_points):
    """
    Integrates the radial Schrödinger equation from r_min to r_max
    for a trial energy E and quantum number l.
    Uses the shooting method: stops at the first sign change of ψ (node).
    """
    r_min = a_0                  # Start integration just above zero to avoid singularity
    r_max = 100 * a_0            # Integrate out to a large distance
    h = (r_max - r_min) / num_points
    r_values = np.linspace(r_min, r_max, num_points)
    psi_values = np.zeros(num_points)
    dpsi_dr_values = np.zeros(num_points)
    psi_values[0] = psi_min
    dpsi_dr_values[0] = dpsi_dr_guess

    psi_sign_changed = False
    for i in range(1, num_points):
        r = r_values[i]
        y = np.array([psi_values[i - 1], dpsi_dr_values[i - 1]])
        y = runge_kutta_step(r, y, h, E, l)
        psi_values[i] = y[0]
        dpsi_dr_values[i] = y[1]

        # If the wavefunction changes sign, this node marks a solution candidate
        if psi_values[i] * psi_values[i - 1] < 0:
            psi_sign_changed = True
        if psi_sign_changed:
            break

    return r_values, psi_values

# ============================================
# 5. SHOOTING METHOD: ENERGY SCAN LOOP
# ============================================
eigenvalues = []
for n in range(1, num_eigenvalues + 1):
    # Use analytical hydrogen-like energies as first guess (in eV)
    E_guess = -13.606 / n ** 2
    dE = 0.1
    convergence_threshold = 1e-50

    while True:
        r_values, psi_values = solve_schrodinger(E_guess, l, num_points)
        psi_at_infinity = psi_values[-1]
        # If wavefunction at large r is close to zero, accept energy
        if np.abs(psi_at_infinity) < convergence_threshold:
            eigenvalues.append(E_guess)
            break
        # Otherwise, adjust energy (coarse-to-fine) and retry
        dE *= 0.05
        E_guess -= dE
        if dE < 1e-50:
            break

print("Energy eigenvalues for l = 0 and n :")
for i, eigenvalue in enumerate(eigenvalues):
    print(f"Eigenvalue {i + 1}: {eigenvalue} eV")

# ============================================
# 6. ANALYTICAL HYDROGENIC RADIAL WAVE FUNCTIONS (Physics)
# ============================================
def hydrogen_wavefunction(n, L, r):
    """
    Computes the normalized analytical radial wavefunction for hydrogen atom.
    Uses generalized Laguerre polynomials.
    """
    R_nl = np.sqrt(((2 / (n * a)) ** 3) * (factorial(n - 1 - L) / (2 * n * factorial(n + L))))
    R_nl *= ((2 * r / (n * a)) ** L) * np.exp(-r / (n * a))
    if n - 1 - L >= 0:
        R_nl *= genlaguerre(n - 1 - L, 2 * L + 1)(2 * r / (n * a))
    return R_nl

# ============================================
# 7. PLOT ANALYTICAL RADIAL WAVE FUNCTIONS (Visualization)
# ============================================
r_min = a_0
r_max = 100 * a_0
r_values = np.linspace(r_min, r_max, num_points)

plt.figure(figsize=(12, 6))
for n in range(1, 11):  # Plot for n=1 to n=10
    psi = hydrogen_wavefunction(n, 2, r_values)
    plt.plot(r_values, psi, label=f'n={n}')

plt.xlabel('r (m)')
plt.ylabel('Radial wave function')
plt.title('Hydrogen Atom: Analytical Radial Wavefunctions')
plt.legend()
plt.grid(True)
plt.show()
