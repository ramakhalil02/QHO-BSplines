import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh


# Constants
k = 7  
r_min = 0.0
r_max = 30.0
N_pts = 20
N_quad = 10
ell = 0
delta = 1e-12
alpha = -1  
Z = 1
n_max = 5

def bsplgen(x_vals, tknot, k):
    N_basis = len(tknot) - k
    n_x = len(x_vals)
    Bval = np.zeros((n_x, N_basis))
    dBval = np.zeros((n_x, N_basis))

    for xi, x in enumerate(x_vals):
        B = np.zeros((len(tknot), k))
        for i in range(N_basis + k - 1):
            if tknot[i] <= x < tknot[i + 1] or (x == tknot[-1] and tknot[i + 1] == tknot[-1]):
                B[i, 0] = 1.0
        for j in range(1, k):
            for i in range(N_basis + k - j - 1):
                t1, t2 = tknot[i], tknot[i + j]
                t3, t4 = tknot[i + 1], tknot[i + j + 1]
                left = (x - t1) / (t2 - t1) * B[i, j - 1] if t2 != t1 else 0.0
                right = (t4 - x) / (t4 - t3) * B[i + 1, j - 1] if t4 != t3 else 0.0
                B[i, j] = left + right

        for i in range(N_basis):
            Bval[xi, i] = B[i, k - 1]

        for i in range(N_basis):
            if tknot[i + k - 1] != tknot[i]:
                dBval[xi, i] += (k - 1) * B[i, k - 2] / (tknot[i + k - 1] - tknot[i])
            if tknot[i + k] != tknot[i + 1]:
                dBval[xi, i] -= (k - 1) * B[i + 1, k - 2] / (tknot[i + k] - tknot[i + 1])

    return Bval, dBval


def gauss_legendre(n):
    return np.polynomial.legendre.leggauss(n)


def generate_linear_knots(r_min, r_max, N_pts, k):
    knots = np.linspace(r_min+delta, r_max, N_pts, endpoint=False)
    knots = np.concatenate((np.full(k, r_min), knots, np.full(k, r_max)))
    return knots

def generate_exponential_knots(r_min, r_max, N_pts, k, alpha):
    xi = np.linspace(0+delta, 1, N_pts)
    exp_knots = r_max * ( 1 - np.exp(alpha * xi) )
    tknot = np.concatenate((np.full(k, r_min), exp_knots, np.full(k, r_max)))
    return tknot

def generate_quadratic_knots(r_min, r_max, N_pts, k):
    xi = np.linspace(0 + 1e-6 , 1, N_pts)
    r_internal = (r_max - r_min) * xi**2
    knots = np.concatenate((np.full(k, r_min), r_internal, np.full(k, r_max)))
    return knots

def spherical_charge_potential(r, R0, Z):
    V = np.zeros_like(r)
    inside = r <= R0
    outside = ~inside
    V[inside] = -Z / (2 * R0) * (3 - (r[inside]**2) / R0**2)
    V[outside] = -Z / r[outside]
    #print(V)
    return V



def build_hydrogen_H(tknot, k, N_quad, ell, potential_type, R0, Z):
    N_basis = len(tknot) - k
    H = np.zeros((N_basis, N_basis))
    S = np.zeros((N_basis, N_basis))
    xi, wi = gauss_legendre(N_quad)

    for i in range(N_basis):
        for j in range(N_basis):
            left = max(tknot[i], tknot[j])
            right = min(tknot[i + k], tknot[j + k])
            if left >= right:
                continue

            m_start = np.searchsorted(tknot, left, side='left') - 1
            m_end = np.searchsorted(tknot, right, side='right') - 1

            for m in range(m_start, m_end):
                a, b = tknot[m], tknot[m+1]
                if b <= a:
                    continue
                xq = 0.5 * (b - a) * xi + 0.5 * (b + a)
                wq = 0.5 * (b - a) * wi

                Bval, dBval = bsplgen(xq, tknot, k)

                if potential_type == "coulomb":
                    V = -Z / xq
                elif potential_type == "sphere":
                    V = spherical_charge_potential(xq, R0, Z)

                integrand_H = (
                     0.5 * dBval[:, i] * dBval[:, j]
                    + 0.5 *(ell * (ell + 1)) * Bval[:, i] * Bval[:, j] / (xq ** 2)
                    + V * Bval[:, i] * Bval[:, j]
                )
                integrand_S = Bval[:, i] * Bval[:, j]

                H[i, j] += np.sum(wq * integrand_H)
                S[i, j] += np.sum(wq * integrand_S)

    return H, S



def analytical_eigenvalues(n_max, Z):
    return np.array([- (0.5 * Z**2) / (n ** 2) for n in range(1, n_max + 1)])

def analytical_radial_wavefunction(n, r, Z, ell):
    if n == 1:
        return 2 * (Z)**(3/2) * np.exp(-Z*r)
    elif n == 2:
        if ell == 0:
            return  (Z/2)**(3/2) * (2 - Z*r) * np.exp(-(Z*r)/2)
        elif ell == 1:
            return (1/np.sqrt(3)) * (Z/2)**(3/2) * (Z*r) * np.exp(-(Z*r)/2)
    elif n == 3:
        if ell == 0:
            return (2/27) * (Z/3)**(3/2) * ( 27 - 18*Z*r + 2*(Z*r)**2 ) * np.exp(-(Z*r)/3)
        elif ell == 1:
            return (1/27) * ( (2*Z) / 3)**(3/2) * (6*(Z*r) - (Z*r)**2) * np.exp(-(Z*r)/3)
        elif ell == 2:
            return (4/(27*np.sqrt(10))) * (Z/3)**(3/2) * (Z*r)**2 * np.exp(-(Z*r)/3)



def solve_hydrogen(ell, n_max, knots_type, potential_type, R0, Z):
    if knots_type == "linear":
        knots = generate_linear_knots(r_min, r_max, N_pts, k)
    elif knots_type == "exponential":
        knots = generate_exponential_knots(r_min, r_max, N_pts, k, alpha)
    elif knots_type == "quadratic":
        knots = generate_quadratic_knots(r_min, r_max, N_pts, k)

    H, S = build_hydrogen_H(knots, k, N_quad, ell, potential_type, R0, Z)

    H_reduced = H[1:-1, 1:-1]
    S_reduced = S[1:-1, 1:-1]
    
    eigvals_analytical = analytical_eigenvalues(n_max, Z)

    eigvals, eigvecs = eigh(H_reduced, S_reduced, eigvals_only=False)

    print(f"Numerical Eigenvalues ({knots_type} knots, {potential_type}):")
    print(np.round(eigvals[:n_max], 5))
    
    print("Analytical Eigenvalues: ")
    print(np.round(eigvals_analytical, 5))

    return eigvals, eigvecs, knots



def plot_wavefunction(knots, eigvecs, Z, ell, n, knots_type="linear"):
    r_plot = np.linspace(r_min, r_max, 500)
    Bval, _ = bsplgen(r_plot, knots, k)
    Bval = Bval[:, 1:-1]  

    psi_num = Bval @ eigvecs[:, n - 1]

    norm = np.sqrt(np.trapz(psi_num**2, r_plot))
    psi_num /= norm

    psi_ana = analytical_radial_wavefunction(n, r_plot, Z, ell)
    
    if np.sign(psi_num[1]) != np.sign(psi_ana[1]):
        psi_num *= -1

    plt.figure(figsize=(8,5))
    plt.plot(r_plot, psi_num, label=f'Numerical n={n}')
    plt.plot(r_plot, r_plot*psi_ana, '--', label=f'Analytical n={n}')
    plt.xlabel('$r$')
    plt.ylabel('$\psi(r)$')
    plt.title(f'Numerical vs Analytical Wavefunction for n={n} ({knots_type} knots)')
    plt.legend()
    plt.grid()
    plt.show()

def compare_probability_distributions(knots1, eigvecs1, Z1, knots2, eigvecs2, Z2, n):
    r_plot = np.linspace(r_min + delta, r_max, 500)
    Bval1, _ = bsplgen(r_plot, knots1, k)
    Bval1 = Bval1[:, 1:-1]

    Bval2, _ = bsplgen(r_plot, knots2, k)
    Bval2 = Bval2[:, 1:-1]

    psi1 = Bval1 @ eigvecs1[:, n-1]
    
    norm1 = np.sqrt(np.trapz(psi1**2, r_plot))
    psi1 /= norm1

    psi2 = Bval2 @ eigvecs2[:, n-1]
    
    norm2 = np.sqrt(np.trapz(psi2**2, r_plot))
    psi2 /= norm2
    
    if np.sign(psi1[10]) != np.sign(psi2[10]):
        psi1 *= -1

    plt.figure(figsize=(8, 5))
    plt.plot(r_plot, psi1, label='Hydrogen')
    plt.plot(r_plot, psi2, label='Hydrogen-like Uranium')
    plt.xlabel('$r$')
    plt.ylabel('$P_n(r)$')
    plt.title('Radial Wavefunction')
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    #eigvals, eigvecs, knots = solve_hydrogen(ell, n_max, knots_type="linear", potential_type="coulomb", R0=1, Z=1)
    
    #plot_wavefunction(knots, eigvecs, Z, ell, n=2, knots_type="linear") 
    #plot_wavefunction(knots, eigvecs, Z, ell, n=3, knots_type="exponential") 
    #plot_wavefunction(knots, eigvecs, Z, ell, n=1, knots_type="quadratic")  
    
    # Z1 = 1
    # R0_H = 1.5  
    # eigvals_C, eigvecs_H, knots_H = solve_hydrogen(ell, n_max,  knots_type="linear", potential_type="sphere", R0=R0_H, Z=Z1)

    # Z2 = 92
    # R0_U = 7.45 
    # eigvals_U, eigvecs_U, knots_U = solve_hydrogen(ell, n_max, knots_type="linear", potential_type="sphere", R0=R0_U, Z=Z2)
    
    # compare_probability_distributions(knots_H, eigvecs_H, Z1, knots_U, eigvecs_U, Z2, n=1)
    
    
    




    

