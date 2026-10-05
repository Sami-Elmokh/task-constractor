"""Build problem.json (SciCode schema) and test_data.h5 (SciCode HDF5 layout) for SME01.

Run from this folder:  python build_task.py
"""
import json
import re
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).parent
PROBLEM_ID = "SME01"
DEPENDENCIES = "import numpy as np"

GRID = "t_grid = 2*np.pi*np.arange(N)/N\n"
# closed form of the 2*pi-periodic, even response to |sin t| (CCP PSI 2003, parts I and III):
# on [0, pi], phi(u) = (sin u - u cos u) / 2, extended evenly and 2*pi-periodically
PHI_ABS_SIN = (
    "s = np.mod(t_grid, 2*np.pi)\n"
    "u = np.minimum(s, 2*np.pi - s)\n"
    "mu = 0.5*(np.sin(u) - u*np.cos(u))\n"
)

# (setup code, call expression, comparison) per test; comparison is "allclose" or "is_none"
TESTS = {
    1: [
        ("N = 33\n" + GRID + "mu = np.cos(t_grid) + 0.5*np.cos(3*t_grid)\n",
         "forcing_resonant_modes(mu)", "allclose"),
        ("N = 101\n" + GRID + "mu = np.abs(np.sin(t_grid))\n",
         "forcing_resonant_modes(mu)", "allclose"),
        ("N = 65\n" + GRID + "mu = 2*np.sin(t_grid + 0.3) + 0.4*np.sin(5*t_grid)\n",
         "forcing_resonant_modes(mu)", "allclose"),
    ],
    2: [
        ("N = 33\n" + GRID + "mu = np.sin(2*t_grid) + 0.3*np.cos(3*t_grid)\n",
         "periodic_response(mu, 1e-8)", "allclose"),
        ("N = 101\n" + GRID + "mu = np.abs(np.sin(t_grid))\n",
         "periodic_response(mu, 1e-3)", "allclose"),
        ("N = 65\n" + GRID + "mu = np.sin(t_grid) + np.cos(2*t_grid)\n",
         "periodic_response(mu, 1e-8)", "is_none"),
    ],
    3: [
        ("N = 49\n" + GRID + "mu = np.cos(t_grid)\nt = np.linspace(0, 20*np.pi, 400)\n",
         "driven_response(mu, t)", "allclose"),
        ("N = 65\n" + GRID + "mu = 2*np.sin(t_grid + 0.3) + 0.4*np.sin(5*t_grid)\n"
         "t = np.linspace(0, 30, 500)\n",
         "driven_response(mu, t)", "allclose"),
        ("N = 101\n" + GRID + PHI_ABS_SIN + "t = np.linspace(0, 40*np.pi, 800)\n",
         "driven_response(mu, t)", "allclose"),
    ],
}

GENERAL_TESTS = [
    ("N = 101\n" + GRID + "phi = periodic_response(np.abs(np.sin(t_grid)), 1e-3)\n"
     "t = np.linspace(0, 40*np.pi, 800)\n",
     "driven_response(phi, t)", "allclose"),
]

PROMPTS = {
    1: (
        "Write a function that extracts the component of the sampled forcing at the natural "
        "frequency of the oscillator. Identify mu with its trigonometric interpolant of degree "
        "M = (N-1)/2 through the N samples, mu(t) = a0/2 + sum_{k=1}^{M} [a_k cos(kt) + b_k sin(kt)], "
        "and return the array [a1, b1] of its cos(t) and sin(t) coefficients."
    ),
    2: (
        "Write a function that decides whether the response y(t) of y'' + y = mu(t), "
        "y(0) = y'(0) = 0, is 2*pi-periodic, and returns it on the sampling grid when it is. "
        "With [a1, b1] from the previous step, if sqrt(a1^2 + b1^2) > tol return None. "
        "Otherwise remove the cos(t) and sin(t) modes from the interpolant (treat them as "
        "numerical noise), solve the initial value problem exactly for the remaining "
        "trigonometric-polynomial forcing, and return y(t_j) for j = 0, ..., N-1 as a 1D array."
    ),
    3: (
        "Write a function that evaluates the exact solution y(t) of y'' + y = mu(t), "
        "y(0) = y'(0) = 0, at arbitrary times t, where mu is the full trigonometric interpolant "
        "of the samples, including its cos(t) and sin(t) modes. The function must remain valid "
        "at resonance, where the response is not bounded, and must return y at every requested "
        "time as a 1D array."
    ),
}

BACKGROUNDS = {
    1: (
        "Background\n"
        "For odd N the interpolant is unique and its complex coefficients are the scaled DFT "
        "c_k = (1/N) sum_j mu_j exp(-i k t_j), |k| <= M. Then a_k = 2 Re c_k and b_k = -2 Im c_k, "
        "i.e. a1 = (2/N) sum_j mu_j cos(t_j) and b1 = (2/N) sum_j mu_j sin(t_j). "
        "These are the only Fourier modes in resonance with the free oscillation cos(t), sin(t)."
    ),
    2: (
        "Background\n"
        "In Fourier space y'' + y = mu becomes (1 - k^2) y_k = c_k. For k != +-1 the periodic "
        "particular solution is y_p = sum_k c_k/(1-k^2) exp(ikt). The k = +-1 equations read "
        "0 * y_{+-1} = c_{+-1}: they are solvable only if a1 = b1 = 0 (Fredholm alternative), "
        "which is the necessary and sufficient condition for the response to be 2*pi-periodic. "
        "The initial conditions are enforced by adding the free oscillation A cos t + B sin t "
        "with A = -y_p(0) and B = -y_p'(0), where y_p'(0) = sum_k i k y_k."
    ),
    3: (
        "Background\n"
        "The resonant forcing a1 cos t + b1 sin t has the secular particular solution "
        "y_s(t) = (t/2) (a1 sin t - b1 cos t), which grows linearly with amplitude rate "
        "sqrt(a1^2 + b1^2)/2. Then y = y_p + y_s + A cos t + B sin t with A = -y_p(0) and "
        "B = -(y_p'(0) + y_s'(0)) = -y_p'(0) + b1/2."
    ),
}

PROBLEM_DESCRIPTION = (
    "Consider an undamped harmonic oscillator with natural frequency omega0 = 1 driven by a "
    "2*pi-periodic force mu(t), starting from rest:\n"
    "    y''(t) + y(t) = mu(t),   y(0) = y'(0) = 0.\n"
    "The forcing is known only through N equally spaced samples mu_j = mu(t_j), "
    "t_j = 2*pi*j/N, j = 0, ..., N-1, with N odd. Throughout, mu is identified with the unique "
    "trigonometric interpolant of degree (N-1)/2 through these samples, so every answer below is "
    "exact for that interpolant (no discretization error is allowed or expected). "
    "Implement functions that (1) extract the component of the forcing at the natural frequency, "
    "(2) decide whether the response is 2*pi-periodic and, if so, return it on the sampling grid, "
    "and (3) evaluate the exact response at arbitrary times, including the linearly growing "
    "response that appears at resonance."
)

PROBLEM_IO = (
    "Input\n"
    "mu:  1D float array of odd length N >= 3, samples of the forcing at t_j = 2*pi*j/N\n"
    "tol: float, threshold on sqrt(a1^2 + b1^2) used to declare resonance (step 2)\n"
    "t:   1D float array of evaluation times (step 3)\n"
    "Output\n"
    "Step 1: 1D array [a1, b1]; step 2: 1D array of length N or None; "
    "step 3: 1D array y(t) of the same length as t"
)

PROBLEM_BACKGROUND = (
    "Background\n"
    "Writing mu as a trigonometric polynomial reduces the ODE to independent scalar equations "
    "(1 - k^2) y_k = c_k for each Fourier mode. All modes except k = +-1 give bounded periodic "
    "responses; the k = +-1 modes are resonant with the free oscillation and, when present, "
    "force a response growing linearly in time. Hence the response is 2*pi-periodic if and only "
    "if the first Fourier coefficients a1, b1 of mu vanish. For mu = |sin t| they do in the "
    "continuum (on an odd grid, aliasing of the even harmonics leaves a residual a1 = O(N^-2), "
    "which the tolerance tol is meant to absorb); the "
    "periodic response is even with a1(phi) = -pi/4 != 0, so driving the oscillator with phi "
    "itself produces an unbounded response."
)


def load_gold(step):
    return (HERE / "gold" / f"step{step}.py").read_text(encoding="utf-8")


def split_header(code):
    """Function header = def line + docstring, as in SciCode's function_header field."""
    m = re.match(r"(def .*?'''.*?''')\n", code, flags=re.S)
    return m.group(1)


def render_test(setup, call, kind):
    if kind == "allclose":
        return setup + f"assert np.allclose({call}, target)"
    return setup + f"assert ({call} is None) == target"


def evaluate(setup, call, kind, namespace):
    ns = dict(namespace)
    exec(setup, ns)
    value = eval(call, ns)
    return (value is None) if kind == "is_none" else np.asarray(value)


def main():
    gold_ns = {"np": np}
    sub_steps = []
    targets = {}
    for step in (1, 2, 3):
        code = load_gold(step)
        exec(code, gold_ns)
        step_id = f"{PROBLEM_ID}.{step}"
        sub_steps.append({
            "step_number": step_id,
            "step_description_prompt": PROMPTS[step],
            "step_background": BACKGROUNDS[step],
            "ground_truth_code": code,
            "function_header": split_header(code),
            "test_cases": [render_test(*tc) for tc in TESTS[step]],
            "return_line": code.rstrip().splitlines()[-1],
        })
        targets[step_id] = [evaluate(*tc, gold_ns) for tc in TESTS[step]]

    problem = {
        "problem_name": "Driven_Oscillator_Resonance",
        "problem_id": PROBLEM_ID,
        "problem_description_main": PROBLEM_DESCRIPTION,
        "problem_background_main": PROBLEM_BACKGROUND,
        "problem_io": PROBLEM_IO,
        "required_dependencies": DEPENDENCIES,
        "sub_steps": sub_steps,
        "general_solution": "\n\n".join(load_gold(s) for s in (1, 2, 3)),
        "general_tests": [render_test(*tc) for tc in GENERAL_TESTS],
    }
    targets[PROBLEM_ID] = [evaluate(*tc, gold_ns) for tc in GENERAL_TESTS]

    (HERE / "problem.json").write_text(json.dumps(problem, indent=2), encoding="utf-8")
    with h5py.File(HERE / "test_data.h5", "w") as f:
        for step_id, values in targets.items():
            for i, value in enumerate(values, start=1):
                f.create_dataset(f"{step_id}/test{i}/var1", data=value)
    print("Wrote problem.json and test_data.h5")
    for step_id, values in targets.items():
        print(step_id, [np.shape(v) for v in values])


if __name__ == "__main__":
    main()
