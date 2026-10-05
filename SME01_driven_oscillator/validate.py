"""Independent validation of the SME01 gold solution and tests.

1. Closed-form solutions for the simple test cases.
2. Direct ODE integration (scipy DOP853) with the interpolant evaluated by explicit sums (no FFT).
3. Continuum constants from CCP PSI 2003 Maths 1 (a1(phi) = -pi/4, int e^-t phi = coth(pi/2)/4,
   unbounded response to phi).
4. Negative controls: plausible wrong implementations must fail at least one hidden test.

Run from this folder:  python validate.py
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp

from scicode.parse.parse import process_hdf5_to_tuple

HERE = Path(__file__).parent
H5 = str(HERE / "test_data.h5")
problem = json.loads((HERE / "problem.json").read_text(encoding="utf-8"))
gold = {"np": np}
for s in problem["sub_steps"]:
    exec(s["ground_truth_code"], gold)
modes, periodic, driven = gold["forcing_resonant_modes"], gold["periodic_response"], gold["driven_response"]

failures = []


def check(name, err, tol):
    ok = err <= tol
    print(f"  [{'ok' if ok else 'FAIL'}] {name}: {err:.2e} (tol {tol:.0e})")
    if not ok:
        failures.append(name)


def grid(N):
    return 2 * np.pi * np.arange(N) / N


def interpolant(mu, drop_resonant=False):
    """Trigonometric interpolant as an explicit cos/sin sum (independent of the FFT path)."""
    N = mu.size
    tg = grid(N)
    M = (N - 1) // 2
    a = [2 / N * np.dot(mu, np.cos(k * tg)) for k in range(M + 1)]
    b = [2 / N * np.dot(mu, np.sin(k * tg)) for k in range(M + 1)]
    if drop_resonant:
        a[1] = b[1] = 0.0
    return lambda t: a[0] / 2 + sum(a[k] * np.cos(k * t) + b[k] * np.sin(k * t) for k in range(1, M + 1))


def ode(mu, t, drop_resonant=False):
    f = interpolant(mu, drop_resonant)
    sol = solve_ivp(lambda s, y: [y[1], f(s) - y[0]], [0, t.max()], [0.0, 0.0], t_eval=t,
                    method="DOP853", rtol=1e-12, atol=1e-12)
    return sol.y[0]


def phi_closed(t):
    s = np.mod(t, 2 * np.pi)
    u = np.minimum(s, 2 * np.pi - s)
    return 0.5 * (np.sin(u) - u * np.cos(u))


print("1. Closed forms")
tg = grid(33)
check("step1 cos t + 0.5 cos 3t -> [1, 0]",
      np.abs(modes(np.cos(tg) + 0.5 * np.cos(3 * tg)) - [1, 0]).max(), 1e-13)
tg = grid(65)
check("step1 2 sin(t+0.3) -> [2 sin 0.3, 2 cos 0.3]",
      np.abs(modes(2 * np.sin(tg + 0.3) + 0.4 * np.sin(5 * tg)) - [2 * np.sin(.3), 2 * np.cos(.3)]).max(), 1e-13)
tg = grid(33)
exact = -np.sin(2 * tg) / 3 - 0.0375 * np.cos(3 * tg) + 0.0375 * np.cos(tg) + 2 / 3 * np.sin(tg)
check("step2 sin 2t + 0.3 cos 3t (hand-derived)",
      np.abs(periodic(np.sin(2 * tg) + 0.3 * np.cos(3 * tg), 1e-8) - exact).max(), 1e-13)
t = np.linspace(0, 20 * np.pi, 400)
check("step3 cos t -> (t/2) sin t", np.abs(driven(np.cos(grid(49)), t) - t / 2 * np.sin(t)).max(), 1e-11)

print("2. Direct ODE integration of the interpolant (DOP853, rtol 1e-12)")
cases = [
    ("step2 sin 2t + 0.3 cos 3t", 33, lambda g: np.sin(2 * g) + 0.3 * np.cos(3 * g), None),
    ("step2 |sin t|", 101, lambda g: np.abs(np.sin(g)), None),
    ("step3 2 sin(t+0.3) + 0.4 sin 5t", 65, lambda g: 2 * np.sin(g + 0.3) + 0.4 * np.sin(5 * g),
     np.linspace(0, 30, 500)),
    ("step3 forcing = phi (III.5)", 101, phi_closed, np.linspace(0, 40 * np.pi, 800)),
]
for name, N, f, t in cases:
    mu = f(grid(N))
    if t is None:
        ref = ode(mu, grid(N), drop_resonant=True)
        got = periodic(mu, 1e-3)
    else:
        ref = ode(mu, t)
        got = driven(mu, t)
    check(name, np.abs(got - ref).max() / max(1, np.abs(ref).max()), 1e-8)

print("2b. Variation of constants from the official worked solution (I.1-I.3): phi = sin x G(x) - cos x H(x)")


def coefficients(mu):
    N = mu.size
    tg = grid(N)
    M = (N - 1) // 2
    a = np.array([2 / N * np.dot(mu, np.cos(k * tg)) for k in range(M + 1)])
    b = np.array([2 / N * np.dot(mu, np.sin(k * tg)) for k in range(M + 1)])
    a[0] /= 2  # constant term a0/2
    return a, b


def int_cos(m, x):  # int_0^x cos(m t) dt
    return x if m == 0 else np.sin(m * x) / m


def int_sin(m, x):  # int_0^x sin(m t) dt
    return 0 * x if m == 0 else (1 - np.cos(m * x)) / m


def variation_of_constants(mu, x, drop_resonant=False):
    """F(x) = sin x G(x) - cos x H(x), with G, H integrated exactly mode by mode."""
    a, b = coefficients(mu)
    if drop_resonant:
        a[1] = b[1] = 0.0
    G = np.zeros_like(x)
    H = np.zeros_like(x)
    for k in range(len(a)):
        # cos kt cos t = [cos(k-1)t + cos(k+1)t]/2 ; sin kt cos t = [sin(k+1)t + sin(k-1)t]/2
        G += a[k] * (int_cos(k - 1, x) + int_cos(k + 1, x)) / 2 + b[k] * (int_sin(k + 1, x) + int_sin(k - 1, x)) / 2
        # cos kt sin t = [sin(k+1)t - sin(k-1)t]/2 ; sin kt sin t = [cos(k-1)t - cos(k+1)t]/2
        H += a[k] * (int_sin(k + 1, x) - int_sin(k - 1, x)) / 2 + b[k] * (int_cos(k - 1, x) - int_cos(k + 1, x)) / 2
    return np.sin(x) * G - np.cos(x) * H


x = np.linspace(0, 20 * np.pi, 400)
check("I.4.6 mu = sin: phi = (sin x - x cos x)/2",
      np.abs(driven(np.sin(grid(49)), x) - (np.sin(x) - x * np.cos(x)) / 2).max(), 1e-11)
check("I.4.6 mu = cos: phi = x sin x / 2",
      np.abs(driven(np.cos(grid(49)), x) - x * np.sin(x) / 2).max(), 1e-11)
for name, N, f, t in cases:
    mu = f(grid(N))
    if t is None:
        ref, got = variation_of_constants(mu, grid(N), drop_resonant=True), periodic(mu, 1e-3)
    else:
        ref, got = variation_of_constants(mu, t), driven(mu, t)
    check(name + " vs worked-solution formula", np.abs(got - ref).max() / max(1, np.abs(ref).max()), 1e-10)

print("3. Continuum results of CCP PSI 2003 (convergence as N grows)")
for N in (101, 401, 1601):
    g = grid(N)
    print(f"  N={N:5d}: a1(|sin t| samples) = {modes(np.abs(np.sin(g)))[0]: .3e}"
          f" | max|phi_N - phi_exact| = {np.abs(periodic(np.abs(np.sin(g)), 1e-3) - phi_closed(g)).max():.3e}"
          f" | a1(phi) + pi/4 = {modes(phi_closed(g))[0] + np.pi / 4: .3e}")
g = grid(1601)
check("a1(phi) = -pi/4 (III.3.4)", abs(modes(phi_closed(g))[0] + np.pi / 4), 1e-5)
a_mu, _ = coefficients(np.abs(np.sin(g)))
a_phi, b_phi = coefficients(phi_closed(g))
p = np.arange(1, 6)
check("III.2.1 a_2p(mu) = -4/(pi(4p^2-1))", np.abs(a_mu[2 * p] + 4 / (np.pi * (4 * p**2 - 1))).max(), 1e-5)
check("III.3.3 a_2p(phi) = 4/(pi(4p^2-1)^2), odd a_n = 0, b_n = 0",
      max(np.abs(a_phi[2 * p] - 4 / (np.pi * (4 * p**2 - 1) ** 2)).max(),
          np.abs(a_phi[2 * p + 1]).max(), np.abs(b_phi).max()), 1e-5)
S = sum(1 / ((4 * q**2 - 1) * (16 * q**4 - 1)) for q in range(1, 200000))
check("III.4 sum 1/((4p^2-1)(16p^4-1)) = (pi/4)[coth(pi/2)/4 - 2/pi + pi/8]",
      abs(S - np.pi / 4 * (0.25 / np.tanh(np.pi / 2) - 2 / np.pi + np.pi / 8)), 1e-12)
I_ref = 0.25 / np.tanh(np.pi / 2)
I_num = quad(lambda s: np.exp(-s) * phi_closed(np.array([s]))[0], 0, 60, limit=400)[0]
check("int_0^inf e^-t phi = coth(pi/2)/4 (II.3.2), closed-form phi", abs(I_num - I_ref), 1e-9)
mu_phi = phi_closed(g)
T = np.linspace(0, 200 * np.pi, 4001)
y = driven(mu_phi, T)
rate = np.abs(y[T > 150 * np.pi]).max() / (200 * np.pi)
check("III.5: response to phi grows like (pi/8) t", abs(rate - np.pi / 8), 5e-3)

print("4. Negative controls (each must fail at least one hidden test)")


def run_tests(step_index, funcs):
    step = problem["sub_steps"][step_index]
    targets = process_hdf5_to_tuple(step["step_number"], len(step["test_cases"]), h5py_file=H5)
    fails = 0
    for tc, target in zip(step["test_cases"], targets):
        ns = {"np": np, **funcs, "target": target}
        try:
            exec(tc, ns)
        except Exception:
            fails += 1
    return fails


def variant(step_index, old, new):
    """Gold functions with one textual change in the given step."""
    ns = {"np": np}
    for i, s in enumerate(problem["sub_steps"]):
        code = s["ground_truth_code"]
        if i == step_index:
            assert old in code, old
            code = code.replace(old, new)
        exec(code, ns)
    return {k: v for k, v in ns.items() if callable(v)}


controls = [
    ("step1: 1/N instead of 2/N normalisation", 0,
     "a1 = 2.0 / N * np.dot(mu, np.cos(t))", "a1 = 1.0 / N * np.dot(mu, np.cos(t))"),
    ("step1: sign error on b1 (FFT convention)", 0,
     "b1 = 2.0 / N * np.dot(mu, np.sin(t))", "b1 = -2.0 / N * np.dot(mu, np.sin(t))"),
    ("step2: forgets y'(0) = 0 (no sin t correction)", 1,
     "phi = y_p - y0 * np.cos(t) - dy0 * np.sin(t)", "phi = y_p - y0 * np.cos(t)"),
    ("step2: returns the periodic particular solution only", 1,
     "phi = y_p - y0 * np.cos(t) - dy0 * np.sin(t)", "phi = y_p"),
    ("step3: regularised division instead of secular term", 2,
     "y_s = 0.5 * t * (a1 * np.sin(t) - b1 * np.cos(t))", "y_s = 0.0 * t"),
    ("step3: float wave numbers from fftfreq (k == 1 test misses for N = 49)", 2,
     "k = np.fft.ifftshift(np.arange(-(N // 2), N // 2 + 1))", "k = np.fft.fftfreq(N, d=1.0 / N)"),
    ("step3: forgets y_s'(0) = -b1/2 in the initial slope", 2,
     "dy0 = np.real(np.sum(1j * k * y_hat)) - 0.5 * b1", "dy0 = np.real(np.sum(1j * k * y_hat))"),
]
print(f"  gold: step fails = {[run_tests(i, {k: v for k, v in gold.items() if callable(v)}) for i in range(3)]}")
for name, idx, old, new in controls:
    n = run_tests(idx, variant(idx, old, new))
    print(f"  [{'ok' if n else 'FAIL'}] {name}: fails {n}/3 tests")
    if not n:
        failures.append(name)

print("\nVALIDATION", "PASSED" if not failures else f"FAILED: {failures}")
