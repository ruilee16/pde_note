# Lecture Notes: Solving the 1D Wave Equation via Fourier Transform

**Date:** Jan 27, 2026  
**Subject:** Partial Differential Equations  
**Topic:** The Cauchy Problem, d'Alembert's Formula, and Duhamel's Principle.

---

We seek to solve the initial value problem for the wave equation on the unbounded domain $\mathbb{R}$:

$$
\begin{cases}
u_{tt} - c^2 u_{xx} = 0, & x \in \mathbb{R}, t > 0 \\
u(x,0) = \phi(x) & (\text{Initial Displacement}) \\
u_t(x,0) = \psi(x) & (\text{Initial Velocity})
\end{cases}
$$

### Step 1: Spatial Fourier Transform
We apply the Fourier transform with respect to the spatial variable $x$. Let $\mathcal{F}[u] = \hat{u}$. The transform converts spatial derivatives into algebraic multiplication:

$$ \mathcal{F}[u_{xx}] = (2\pi i \xi)^2 \hat{u} = -4\pi^2 \xi^2 \hat{u} $$

The PDE transforms into an Ordinary Differential Equation (ODE) in time $t$:

$$ \hat{u}_{tt} + (2\pi c \xi)^2 \hat{u} = 0 $$

Let $\omega = 2\pi c \xi$. The general solution to this harmonic oscillator equation is:

$$ \hat{u}(\xi, t) = A(\xi) \cos(\omega t) + B(\xi) \sin(\omega t) $$

### Step 2: Applying Initial Conditions
We determine the coefficients $A(\xi)$ and $B(\xi)$ using the transformed initial conditions:
1.  **At $t=0$:**
   
    $$ \hat{u}(\xi, 0) = A(\xi) = \hat{\phi}(\xi) $$
3.  **Derivative at $t=0$:**

    $$ \hat{u}_t(\xi, 0) = \omega B(\xi) = \hat{\psi}(\xi) \implies B(\xi) = \frac{\hat{\psi}(\xi)}{\omega} = \frac{\hat{\psi}(\xi)}{2\pi c \xi} $$

Substituting these back, the solution in the frequency domain is:

$$ \hat{u}(\xi, t) = \hat{\phi}(\xi) \cos(2\pi c \xi t) + \hat{\psi}(\xi) \frac{\sin(2\pi c \xi t)}{2\pi c \xi} $$

### Step 3: Inverse Transform (Deriving d'Alembert's Formula)
We invert the solution term by term.

**Term 1: The Cosine Term (Displacement)**
Using Euler's formula $\cos(\theta) = \frac{e^{i\theta} + e^{-i\theta}}{2}$:

$$ \hat{\phi}(\xi) \cos(2\pi c \xi t) = \frac{1}{2} \hat{\phi}(\xi) \left( e^{2\pi i \xi (ct)} + e^{-2\pi i \xi (ct)} \right) $$
By the **Shift Theorem** $\mathcal{F}^{-1}[e^{2\pi i \xi h} \hat{f}(\xi)] = f(x+h)$, the inverse transform is:

$$ \frac{1}{2} [\phi(x+ct) + \phi(x-ct)] $$

**Term 2: The Sine Term (Velocity)**
Recall the Fourier transform of the indicator (box) function on the interval $(-a, a)$:

$$ \mathcal{F}[\mathbb{1}_{(-a,a)}](\xi) = \frac{\sin(2\pi a \xi)}{\pi \xi} $$
We manipulate the sine term to match this form. Let $a = ct$:

$$ \hat{\psi}(\xi) \frac{\sin(2\pi c t \xi)}{2\pi c \xi} = \hat{\psi}(\xi) \frac{1}{2c} \underbrace{\frac{\sin(2\pi (ct) \xi)}{\pi \xi}}_{\text{FT of } \mathbb{1}_{(-ct, ct)}} $$
Using the **Convolution Theorem** $\mathcal{F}^{-1}[\hat{f} \cdot \hat{g}] = f * g$:

$$ \mathcal{F}^{-1} \left[ \hat{\psi} \cdot \frac{1}{2c} \widehat{\mathbb{1}_{(-ct, ct)}} \right] = \psi * \left( \frac{1}{2c} \mathbb{1}_{(-ct, ct)} \right) $$
Writing out the convolution integral:

$$ \frac{1}{2c} \int_{-\infty}^{\infty} \mathbb{1}_{(-ct, ct)}(y) \psi(x-y) \, dy $$
Let $z = x-y$. The limits of integration become $x-ct$ to $x+ct$:

$$ = \frac{1}{2c} \int_{x-ct}^{x+ct} \psi(z) \, dz $$

### Final Result: d'Alembert's Formula
Combining both parts:

$$ u(x,t) = \frac{\phi(x+ct) + \phi(x-ct)}{2} + \frac{1}{2c} \int_{x-ct}^{x+ct} \psi(z) \, dz $$

---

## 2. The Non-Homogeneous Problem

We now solve the forced wave equation with zero initial conditions:

$$
\begin{cases}
u_{tt} - c^2 u_{xx} = F(x,t) \\
u(x,0) = 0, \quad u_t(x,0) = 0
\end{cases}
$$

### Step 1: Frequency Domain ODE
Transforming the PDE with respect to $x$:

$$ \hat{u}_{tt} + \omega^2 \hat{u} = \hat{F}(\xi, t) $$
where $\omega = 2\pi c \xi$. This is a non-homogeneous harmonic oscillator.

### Step 2: Variation of Parameters
We seek a particular solution of the form:

$$ \hat{u}(t) = C_1(t) \cos(\omega t) + C_2(t) \sin(\omega t) $$
The derivatives of the parameters must satisfy the system:
1. $C_1' \cos(\omega t) + C_2' \sin(\omega t) = 0$
2. $-C_1' \omega \sin(\omega t) + C_2' \omega \cos(\omega t) = \hat{F}(\xi, t)$

Solving for $C_1'$ and $C_2'$:

$$ C_1'(t) = -\frac{\hat{F}}{\omega} \sin(\omega t), \quad C_2'(t) = \frac{\hat{F}}{\omega} \cos(\omega t) $$
Integrating from $0$ to $t$ (using dummy variable $s$):

$$ C_1(t) = -\int_0^t \frac{\hat{F}(\xi, s)}{\omega} \sin(\omega s) \, ds $$

$$ C_2(t) = \int_0^t \frac{\hat{F}(\xi, s)}{\omega} \cos(\omega s) \, ds $$

Substitute back into the ansatz and use the identity $\sin(A)\cos(B) - \cos(A)\sin(B) = \sin(A-B)$:
$$ \hat{u}(\xi, t) = \int_0^t \hat{F}(\xi, s) \frac{\sin(\omega(t-s))}{\omega} \, ds $$
Substituting $\omega = 2\pi c \xi$:
$$ \hat{u}(\xi, t) = \int_0^t \hat{F}(\xi, s) \frac{\sin(2\pi c \xi (t-s))}{2\pi c \xi} \, ds $$

### Step 3: Inverse Transform (Duhamel's Principle)
We recognize the term multiplying $\hat{F}$ inside the integral as the transform of a box function with width $c(t-s)$:
$$ \frac{\sin(2\pi c \xi (t-s))}{2\pi c \xi} = \frac{1}{2c} \mathcal{F}\left[ \mathbb{1}_{(-(t-s)c, (t-s)c)} \right] $$

Apply the Convolution Theorem to the spatial variables:
$$ \mathcal{F}^{-1} \left[ \hat{F}(\cdot, s) \frac{\sin(\dots)}{\dots} \right] = F(\cdot, s) * \left( \frac{1}{2c} \mathbb{1}_{(-c(t-s), c(t-s))} \right) $$

This yields the integral over the spatial interval $[x-c(t-s), x+c(t-s)]$:
$$ \frac{1}{2c} \int_{x-c(t-s)}^{x+c(t-s)} F(y, s) \, dy $$

### Final Result: Non-Homogeneous Solution
Integrating over the time history $s$ from $0$ to $t$:
$$ u(x,t) = \frac{1}{2c} \int_0^t \int_{x-c(t-s)}^{x+c(t-s)} F(y, s) \, dy \, ds $$
This represents the integral of the source $F$ over the **characteristic triangle** in spacetime.
```
