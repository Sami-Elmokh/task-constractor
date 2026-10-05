"""Independent validation of the SME02 gold solution and tests.

1. Closed forms (harmonic oscillator, pendulum turning points and barriers, tilted pendulum).
2. Exact pendulum period T = 4 K((1 + h)/2) (scipy.special.ellipk).
3. Direct ODE integration (scipy DOP853) of x'' = -f'(x), period measured by events.
4. Results from the exam (Centrale 2007 PSI Maths I): I.A.1, III.E.3, IV.A, IV.D, IV.E.3.
5. Negative controls: plausible wrong implementations must fail at least one hidden test.

Run from this folder:  python validate.py
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipk

from scicode.parse.parse import process_hdf5_to_tuple

HERE = Path(__file__).parent
H5 = str(HERE / "test_data.h5")
problem = json.loads((HERE / "problem.json").read_text(encoding="utf-8"))


def load(replace=None):
    """Gold functions, optionally with one textual change (step index, old, new)."""
    ns = {}
    exec(problem["required_dependencies"], ns)
    for i, s in enumerate(problem["sub_steps"]):
        code = s["ground_truth_code"]
        if replace and replace[0] == i:
            assert replace[1] in code, replace[1]
            code = code.replace(replace[1], replace[2])
        exec(code, ns)
    return ns


gold = load()
tp, period, barrier = gold["turning_points"], gold["oscillation_period"], gold["barrier_energy"]
failures = []


def check(name, err, tol):
    ok = err <= tol
    print(f"  [{'ok' if ok else 'FAIL'}] {name}: {err:.2e} (tol {tol:.0e})")
    if not ok:
        failures.append(name)


def ode_period(fp, a):
    """Time for the orbit started at rest at the left turning point a to return there."""
    ev = lambda t, y: y[1]
    ev.direction = 1
    s = solve_ivp(lambda t, y: [y[1], -fp(y[0])], [0, 1e4], [a, 0.0], events=ev,
                  method="DOP853", rtol=1e-13, atol=1e-14)
    t = s.t_events[0]
    return t[t > 1e-6][0]


pend, pend_p = (lambda x: -np.cos(x)), np.sin
pert, pert_p = (lambda x: -np.cos(x) + x**2 / 20), (lambda x: np.sin(x) + x / 10)
tilt, tilt_p = (lambda x: -np.cos(x) + 0.2 * x), (lambda x: np.sin(x) + 0.2)
x_tilt = -np.arcsin(0.2)

print("1. Closed forms")
check("harmonic f = 2x^2, h = 3: turning points = -+sqrt(h/2)",
      np.abs(tp(lambda x: 2 * x**2, lambda x: 4 * x, 3.0, 0.1) - [-np.sqrt(1.5), np.sqrt(1.5)]).max(), 1e-13)
for h in (0.5, 1 - 1e-6):
    check(f"pendulum h = {h}: turning points = -+arccos(-h)",
          np.abs(tp(pend, pend_p, h, 0.0) - [-np.arccos(-h), np.arccos(-h)]).max(), 1e-12)
check("pendulum barriers = [-pi, pi, 1]", np.abs(barrier(pend, pend_p, 0.0) - [-np.pi, np.pi, 1]).max(), 1e-13)
el, er = -np.pi + np.arcsin(0.2), np.pi + np.arcsin(0.2)
check("tilted barriers = [-pi + asin 0.2, pi + asin 0.2, f(left)]",
      np.abs(barrier(tilt, tilt_p, x_tilt) - [el, er, tilt(el)]).max(), 1e-13)

print("2. Exact pendulum period 4 K((1+h)/2)")
for h in (0.5, 0.99, 1 - 1e-6, 1 - 1e-9):
    ex = 4 * ellipk((1 + h) / 2)
    check(f"h = {h}", abs(period(pend, pend_p, h, 0.0) / ex - 1), 1e-8)

print("3. Direct ODE integration (DOP853, rtol 1e-13)")
for name, f, fp, h, x0 in [
    ("perturbed pendulum h = 0.8", pert, pert_p, 0.8, 0.0),
    ("perturbed pendulum h = 1.6 (orbit over 3 wells)", pert, pert_p, 1.6, 0.0),
    ("tilted pendulum h = 0.3", tilt, tilt_p, 0.3, x_tilt),
    ("tilted pendulum h = h* - 1e-4 (general test)", tilt, tilt_p, tilt(el) - 1e-4, x_tilt),
]:
    a, b = tp(f, fp, h, x0)
    check(name, abs(period(f, fp, h, x0) / ode_period(fp, a) - 1), 1e-9)
    # the ODE started at rest at a must reach b exactly (no hidden barrier in between)
    s = solve_ivp(lambda t, y: [y[1], -fp(y[0])], [0, ode_period(fp, a) / 2], [a, 0.0],
                  method="DOP853", rtol=1e-13, atol=1e-14)
    check(name + ": ODE reaches b", abs(s.y[0, -1] - b), 1e-7)

print("4. Results from the exam")
for h in (0.1, 7.0):
    check(f"I.A.1 harmonic period 2 pi/omega independent of h (h = {h})",
          abs(period(lambda x: 2 * x**2, lambda x: 4 * x, h, 0.0) - np.pi), 1e-12)
h = pert(0.0) + 1e-8
check("III.E.3 small oscillations: T -> 2 pi / sqrt(f''(0)) = 2 pi / sqrt(1.1)",
      abs(period(pert, pert_p, h, 0.0) / (2 * np.pi / np.sqrt(1.1)) - 1), 1e-7)
for h in (1 - 1e-6, 1 - 1e-9):
    check(f"IV.D separatrix: T ~ 2 ln(32/(1-h)) as h -> 1 (h = {h})",
          abs(period(pend, pend_p, h, 0.0) - 2 * np.log(32 / (1 - h))), 10 * (1 - h) * np.log(32 / (1 - h)))
e = 3.5  # IV.E.3: Newton iteration on sin e + e/10 = 0, independent of brentq
for _ in range(50):
    e -= (np.sin(e) + e / 10) / (np.cos(e) + 0.1)
check("IV.E.3 nearest equilibrium right of 0 (Newton)", abs(barrier(pert, pert_p, 0.0)[1] - e), 1e-13)
check("IV.E.1 it is unstable: f''(e) = cos e + 1/10 < 0", max(0.0, np.cos(e) + 0.1), 0.0)

print("5. Negative controls (each must fail at least one hidden test)")


def run_tests(funcs, step_index=None):
    """Number of failed hidden tests for one step (or for the general test if step_index is None)."""
    if step_index is None:
        step_id, tests = problem["problem_id"], problem["general_tests"]
    else:
        step_id = problem["sub_steps"][step_index]["step_number"]
        tests = problem["sub_steps"][step_index]["test_cases"]
    fails = 0
    for tc, target in zip(tests, process_hdf5_to_tuple(step_id, len(tests), h5py_file=H5)):
        try:
            exec(tc, {"np": np, **funcs, "target": target})
        except Exception:
            fails += 1
    return fails


controls = [
    ("step1: fixed-step march, no check for narrow barriers", 0,
     "if sign * fp(x) > 0 and sign * fp(y) < 0:", "if False:"),
    ("step1: stops at the first barrier even when h is above it", 0,
     "if f(xm) >= h:\n                    return brentq(lambda u: f(u) - h, x, xm, xtol=1e-15, rtol=4 * eps)",
     "return xm"),
    ("step1: assumes a symmetric well (a = 2 x0 - b)", 0,
     "ab = np.array([edge(-1), edge(1)])", "b = edge(1)\n    ab = np.array([2 * x0 - b, b])"),
    ("step2: half period only (a -> b)", 1, "T = 2 * (half(a) + half(b))", "T = half(a) + half(b)"),
    ("step2: small-amplitude formula 2 pi / sqrt(f''(x0))", 1,
     "T = 2 * (half(a) + half(b))", "T = 2 * np.pi / np.sqrt((fp(x0 + 1e-5) - fp(x0 - 1e-5)) / 2e-5)"),
    ("step2: midpoint rule with 2000 points in u (singular endpoints)", 1,
     "T = 2 * (half(a) + half(b))",
     "u = a + (np.arange(2000) + 0.5) * (b - a) / 2000\n"
     "    T = 2 * (b - a) / 2000 * sum(1 / np.sqrt(2 * (h - f(v))) for v in u)"),
    ("step3: assumes symmetric barriers (h* = f(e_right))", 2,
     "min(f(e_left), f(e_right))", "f(e_right)"),
    ("step3: takes the higher barrier (max instead of min)", 2,
     "min(f(e_left), f(e_right))", "max(f(e_left), f(e_right))"),
]
gold_funcs = {k: v for k, v in gold.items() if callable(v)}
gold_fails = [run_tests(gold_funcs, i) for i in range(3)] + [run_tests(gold_funcs)]
print(f"  gold: fails per step = {gold_fails[:3]}, general = {gold_fails[3]}")
if any(gold_fails):  # the controls only mean something if the gold passes every test
    failures.append("gold fails its own tests in this harness")
for name, idx, old, new in controls:
    funcs = {k: v for k, v in load((idx, old, new)).items() if callable(v)}
    n = run_tests(funcs, idx)
    total = len(problem["sub_steps"][idx]["test_cases"])
    print(f"  [{'ok' if n else 'FAIL'}] {name}: fails {n}/{total} tests")
    if not n:
        failures.append(name)

# informative: plain adaptive quad in u, a legitimate alternative that should pass
alt = {k: v for k, v in load((1, "T = 2 * (half(a) + half(b))",
                             "T = 2 * quad(lambda u: 1 / np.sqrt(2 * (h - f(u))), a, b, limit=500)[0]")).items()
       if callable(v)}
print(f"  [info] step2 alternative: plain scipy quad in u fails {run_tests(alt, 1)}/4 tests")

print("\nVALIDATION", "PASSED" if not failures else f"FAILED: {failures}")
