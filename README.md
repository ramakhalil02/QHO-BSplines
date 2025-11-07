# B-spline Collocation Method for Eigenvalue Problems

This project is from independent research and development conducted during my **Master of Science program in Computational Physics**. The entire methodology, analysis, and implementation of all numerical solvers were executed solely by me, and the full work is documented in the attached report and code.

---

## Table of Contents

* [About The Project](#about-the-project)
* [Core Objectives](#core-objectives)
* [Languages and Libraries](#languages-and-libraries)
* [Methods Implemented](#methods-implemented)
* [Key Findings](#key-findings)
* [Getting Started](#getting-started)
* [Full Project Report](#full-project-report)
* [Contact](#contact)

---

## About The Project

This work applies the **B-spline Collocation Method** to solve the **radial Schrödinger equation** for central potentials. The primary objective is to accurately compute the energy eigenvalues ($E_{nl}$) and radial wave functions ($R_{nl}(r)$) for the Hydrogen atom.

The project investigates how the method's accuracy and convergence depend on the chosen **knot sequence** and explores complex systems, including a comparison between a standard Coulomb potential and a more realistic **spherical charge distribution** potential in Hydrogen-like systems (like Uranium).

## Core Objectives

1.  **Implement B-spline Collocation for Eigenvalues:** Develop the B-spline basis and collocation matrix used to transform the radial Schrödinger equation into a generalized eigenvalue problem.

2.  **Optimize Knot Sequence:** Investigate different knot distributions (e.g., exponential) to optimize the convergence rate and accuracy, especially for deeply bound states.

3.  **Solve Hydrogen Atom:** Compute the first few energy levels of the Hydrogen atom and compare them to analytical solutions for validation.

4.  **Model Complex Potentials:** Analyze the difference in energy eigenvalues and radial probability distributions when replacing the point-charge Coulomb potential with a finite **spherical charge distribution** (modeling the finite size of the nucleus).

---

## Languages and Libraries

| Category | Tools & Libraries | Competency Demonstrated |
| :--- | :--- | :--- |
| **Language** | Python | Efficient development and numerical handling of systems of linear equations. |
| **Numerical** | NumPy, SciPy | Advanced array manipulation, linear algebra for generalized eigenvalue solution, and B-spline generation. |
| **Visualization** | Matplotlib | Generating high-quality plots of radial wave functions and probability distributions. |

---

## Methods Implemented

The computational core of the project implements the following techniques:

| Method | Role in Project | Key Implementation Detail |
| :--- | :--- | :--- |
| **B-spline Collocation** | Core numerical solver. | Uses B-spline basis functions and collocation points to transform the second-order differential equation into a linear generalized eigenvalue problem ($HC = \lambda SC$). |
| **Radial Schrödinger Eq.** | Governing equation. | Solved using the potential $V(r)$ and the centrifugal term $l(l+1)/r^2$. |
| **Exponential Knot Sequence** | Accuracy enhancement. | Used to improve resolution near the origin ($r=0$), where the wave function changes rapidly due to the strong potential. |
| **Spherical Charge Potential** | High-Z modeling. | Implemented to compare against the ideal Coulomb potential, demonstrating its effect on deeply bound states in Hydrogen-like ions. |

## Key Findings

* **High Accuracy:** The B-spline method accurately reproduced the analytical energy eigenvalues of the Hydrogen atom, confirming its validity.
* **Knot Sequence Importance:** The exponential knot sequence significantly improved the convergence rate and accuracy for low-$r$ phenomena compared to a uniform sequence.
* **Finite Nuclear Size Effect:** The analysis of the spherical charge potential showed a noticeable shift in the ground state radial probability distribution ($l=0$) away from the origin compared to the Coulomb case.
* **Potential Approximation:** The effect of the finite nuclear size potential diminished rapidly for higher angular momentum states ($l \ge 1$), as the centrifugal barrier dominates the behavior near the nucleus.

---

## Getting Started

### Execution

To run the simulation and generate the results and visualizations, execute the core solver script:

```bash
python QHO_Bsplines.py
```
---

## Full Project Report

For a complete breakdown of the theoretical derivations and full results, please see the final project report:

[**Full Project Report (PDF)**](Eigenvalue_Bspline.pdf)

---

## Contact

I'm happy to hear your feedback or answer any questions about this project!

**Author** Rama Khalil 

**Email**  rama.khalil.990@gmail.com
