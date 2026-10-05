"""Build problem.json (SciCode schema) and test_data.h5 (SciCode HDF5 layout) for SME02.

Run from this folder:  python build_task.py
"""
import json
import re
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).parent
PROBLEM_ID = "SME02"
DEPENDENCIES = "import numpy as np\nfrom scipy.integrate import quad\nfrom scipy.optimize import brentq"

# potentials used in the tests (Centrale 2007 PSI Maths I: parts I, IV and IV.E, plus a tilted pendulum)
HARMONIC = "f = lambda x: 2.0*x**2\nfp = lambda x: 4.0*x\n"
PENDULUM = "f = lambda x: -np.cos(x)\nfp = lambda x: np.sin(x)\n"
PERTURBED = "f = lambda x: -np.cos(x) + x**2/20\nfp = lambda x: np.sin(x) + x/10\n"
TILTED = "f = lambda x: -np.cos(x) + 0.2*x\nfp = lambda x: np.sin(x) + 0.2\nx0 = -np.arcsin(0.2)\n"

# (setup code, call expression) per test; all comparisons use np.allclose
TESTS = {
    1: [
        (HARMONIC, "turning_points(f, fp, 3.0, 0.1)"),
        (TILTED, "turning_points(f, fp, 0.3, x0)"),
        (PENDULUM, "turning_points(f, fp, 1 - 1e-6, 0.0)"),
        (PERTURBED, "turning_points(f, fp, 1.6, 0.0)"),
    ],
    2: [
        (HARMONIC, "oscillation_period(f, fp, 3.0, 0.1)"),
        (PENDULUM, "oscillation_period(f, fp, 0.5, 0.0)"),
        (PENDULUM, "oscillation_period(f, fp, 1 - 1e-6, 0.0)"),
        (PERTURBED, "oscillation_period(f, fp, 1.6, 0.0)"),
    ],
    3: [
        (PENDULUM, "barrier_energy(f, fp, 0.0)"),
        (PERTURBED, "barrier_energy(f, fp, 0.0)"),
        (TILTED, "barrier_energy(f, fp, x0)"),
    ],
}

GENERAL_TESTS = [
    (TILTED + "h_star = barrier_energy(f, fp, x0)[2]\n",
     "oscillation_period(f, fp, h_star - 1e-4, x0)"),
]

PROMPTS = {
    1: (
        "Write a function that returns the turning points of the orbit of energy h passing through "
        "x0, i.e. the ends a < x0 < b of the maximal open interval containing x0 on which f(x) < h. "
        "You may assume f(x0) < h, that a and b are finite with f(a) = f(b) = h, and that f'(a) != 0 "
        "and f'(b) != 0. Note that the region where f >= h separating two wells can be very narrow "
        "when h is close to the top of a barrier, and that an orbit with h above a barrier extends "
        "past it. Return np.array([a, b]) accurate to about 1e-12."
    ),
    2: (
        "Write a function that returns the period T of the closed orbit of energy h passing through "
        "x0, under the same assumptions as the previous step. The motion goes from a to b and back, "
        "so T = 2 * integral_a^b du / sqrt(2 (h - f(u))). The integrand is infinite at both turning "
        "points, and h may be very close to the energy of a barrier top, where the period becomes "
        "very long. Return T as a float with a relative accuracy of at least 1e-8."
    ),
    3: (
        "Write a function that, given a strict local minimum x0 of f (a stable equilibrium), returns "
        "the nearest unstable equilibria e_left < x0 < e_right (the nearest local maxima of f on each "
        "side) and the escape energy h_star = min(f(e_left), f(e_right)), the largest energy for which "
        "orbits around x0 stay confined to the well. The potential need not be symmetric. Return "
        "np.array([e_left, e_right, h_star])."
    ),
}

BACKGROUNDS = {
    1: (
        "Background\n"
        "Energy conservation y^2/2 + f(x) = h confines the motion to {f(x) <= h}. The orbit through x0 "
        "lives in the connected component ]a, b[ of {f < h} that contains x0; a and b are the turning "
        "points where the velocity vanishes. A robust search marches away from x0 and stops at the first "
        "point where f reaches h. Between two sample points, a barrier can be missed if it is narrower "
        "than the step; this is detected by a sign change of f' from + to - (resp. - to +), in which case "
        "the local maximum is located and compared with h."
    ),
    2: (
        "Background\n"
        "Solving y^2/2 + f(x) = h for y > 0 gives dt = dx / sqrt(2 (h - f(x))), so the half period is the "
        "time from a to b. Near a simple turning point h - f(u) ~ |f'(end)| |u - end|, an integrable "
        "1/sqrt singularity; the substitution u = end - L s^2 (L = end - x0) turns it into a bounded "
        "integrand. When h approaches the energy of a barrier top, f'(end) -> 0 and the integrand develops "
        "a narrow peak; the period diverges logarithmically (for the pendulum f = -cos x, "
        "T = 4 K(m) with m = (1 + h)/2, K the complete elliptic integral of the first kind)."
    ),
    3: (
        "Background\n"
        "Equilibria are the zeros of f'. An equilibrium e is stable when f''(e) > 0 (a minimum of the "
        "potential) and unstable when f''(e) < 0 (a maximum). Starting from a minimum and moving away, "
        "f' keeps the sign of the displacement until the first maximum, where it changes sign. An orbit "
        "around x0 is closed as long as h stays below both barriers, so the escape energy is the lower "
        "of the two barrier heights."
    ),
}

PROBLEM_DESCRIPTION = (
    "Consider a particle of unit mass moving on a line in a smooth potential f, governed by the "
    "conservative equation\n"
    "    x''(t) = -f'(x(t)).\n"
    "Along every trajectory the energy y^2/2 + f(x), with y = x', is constant; call it h. "
    "When the orbit of energy h through a point x0 is trapped between two simple turning points "
    "a < x0 < b, the motion is periodic. The potential is given as two Python callables, f and its "
    "derivative fp. Implement functions that (1) find the turning points a and b, (2) compute the "
    "period of the oscillation accurately, including energies extremely close to a barrier top "
    "(near a separatrix), and (3) find the barriers bounding a potential well and its escape energy."
)

PROBLEM_IO = (
    "Input\n"
    "f:  callable, potential energy f(x) (vectorization not required)\n"
    "fp: callable, derivative f'(x)\n"
    "h:  float, energy of the orbit (steps 1 and 2)\n"
    "x0: float, a point of the orbit (steps 1 and 2) or the bottom of a well (step 3)\n"
    "Output\n"
    "Step 1: np.array([a, b]); step 2: float period T; "
    "step 3: np.array([e_left, e_right, h_star])"
)

PROBLEM_BACKGROUND = (
    "Background\n"
    "Writing y = x' gives the planar system X' = (y, -f'(x)), whose trajectories lie on the level "
    "curves y^2/2 + f(x) = h of the energy. Near a stable equilibrium e (f'(e) = 0, f''(e) > 0) the "
    "motion is a small oscillation of period 2 pi / sqrt(f''(e)); as h grows the period changes, "
    "and it diverges when the orbit approaches an unstable equilibrium (the separatrix). For the "
    "pendulum f = -cos x the separatrix is h = 1; for the perturbed pendulum f = -cos x + x^2/20 "
    "the barriers are lowered and orbits above them encircle several wells; for a tilted pendulum "
    "f = -cos x + c x the two barriers of a well have different heights."
)


def load_gold(step):
    return (HERE / "gold" / f"step{step}.py").read_text(encoding="utf-8")


def split_header(code):
    """Function header = def line + docstring, as in SciCode's function_header field."""
    return re.match(r"(def .*?'''.*?''')\n", code, flags=re.S).group(1)


def render_test(setup, call):
    return setup + f"assert np.allclose({call}, target)"


def evaluate(setup, call, namespace):
    ns = dict(namespace)
    exec(setup, ns)
    return np.asarray(eval(call, ns))


def main():
    gold_ns = {}
    exec(DEPENDENCIES, gold_ns)
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
        "problem_name": "Conservative_Oscillator_Period",
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
        print(step_id, [np.round(v, 10).tolist() for v in values])


if __name__ == "__main__":
    main()
