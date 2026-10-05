# SME02 — Conservative_Oscillator_Period

**Status:** draft, pending SME sign-off  
**Domain:** Physics → Computational Physics  

## 1. What the task is about

A particle of unit mass moves on a line in a potential f:

    x''(t) = −f'(x(t))

Its energy h = x'²/2 + f(x) is conserved (exam preliminary question 2). When the particle is trapped between two turning points a < b, where the velocity vanishes and f(a) = f(b) = h, the motion is periodic with period (exam question II.E.1.d)

    T(h) = 2 ∫ₐᵇ du / √(2(h − f(u)))

The physics question: **how long does one oscillation take, for any potential and any energy, including energies just below the top of a barrier?** Close to small oscillations, T → 2π/√f''(e) (III.E.3). As the energy approaches an unstable equilibrium (the separatrix of the pendulum, IV.D), the period diverges logarithmically.

The model receives the potential f and its derivative f' as Python functions. The tests use the potentials from the exam: harmonic (part I), pendulum f = −cos x (part IV) and perturbed pendulum f = −cos x + x²/20 (IV.E). They also use a tilted pendulum f = −cos x + 0.2x, whose two barriers have different heights.

## 2. What the model must write

| Step | Function | What it must do |
|---|---|---|
| SME02.1 | `turning_points(f, fp, h, x0)` → `[a, b]` | Find the ends of the maximal interval containing x0 on which f < h (II.A.1) |
| SME02.2 | `oscillation_period(f, fp, h, x0)` → `T` | Compute the period integral accurately, even very close to a barrier top |
| SME02.3 | `barrier_energy(f, fp, x0)` → `[e_left, e_right, h*]` | From the bottom x0 of a well, find the nearest unstable equilibria on each side (IV.E.3) and the escape energy h* = the lower of the two barrier heights |
| general test | step 3 feeding step 2 | Period of the tilted pendulum at h* − 10⁻⁴, just below its escape energy (T ≈ 16.29, about 2.6 times the small-oscillation period) |

The exact prompts, function headers and test code are in `problem.json` (SciCode format). The expected values are in `test_data.h5`, and the reference (gold) solution is in `gold/`.

**What makes it hard:**
- the period integrand is **infinite at both turning points** (an integrable 1/√ singularity);
- near a separatrix it develops a **narrow, tall peak**, and the period diverges logarithmically;
- **the turning-point search can step over a narrow barrier.** At h = 1 − 10⁻⁶ the region where f ≥ h around π is only 2.8·10⁻³ wide. Conversely, when h is above a barrier (perturbed pendulum, h = 1.6), the orbit must extend past it and enclose three wells (IV.E.2).

## 3. Errors and why they matter

### 3a. Errors the hidden tests catch

Each error below is a plausible mistake. Each one fails at least one hidden test, while the gold solution passes them all (checked by `validate.py`). This shows that a failure means a real scientific or coding error, not an ambiguity in the task.

| Error | Tests failed | Why it matters |
|---|---|---|
| Fixed-step march without checking for narrow barriers | 1/4 (step 1) | Near the separatrix the march steps over the barrier at π and continues into the next well. The "orbit" then crosses an unstable equilibrium, which is physically impossible: a particle at energy h < 1 can never go over the top. The resulting period was off by a factor of 7 (241 instead of 34.6) |
| Stops at the first barrier even when h is above it | 2/4 (step 1) | It misses that an orbit with enough energy goes over the barrier. Physically, the particle is wrongly confined to one well (IV.E.2) |
| Assumes a symmetric well (a = 2x0 − b) | 2/4 (step 1) | It only holds for even potentials with x0 at the bottom. A tilted potential or an off-center starting point breaks it |
| Half period only (∫ from a to b) | 4/4 (step 2) | It forgets the return trip from b to a. The motion is periodic only after both halves |
| Small-amplitude formula 2π/√f''(x0) | 3/4 (step 2) | It's only the limit h → f(e) (III.E.3). For a pendulum at large amplitude, the period is longer and diverges at the separatrix. This is exactly the small-oscillation approximation that part IV of the exam explicitly asks students not to make |
| Midpoint rule with 2000 points in u | 4/4 (step 2) | Ignoring the endpoint singularity makes the error decrease only like N^−1/2. Even 2000 points miss the 10⁻⁵ tolerance |
| Assumes symmetric barriers (h* = f(e_right)) | 1/3 (step 3) | In the tilted potential, the left barrier (0.39) is much lower than the right one (1.65). The particle escapes over the lower barrier, so taking the right one overestimates the escape energy by a factor of 4 |
| Takes the higher barrier (max instead of min) | 1/3 (step 3) | Same physics error: confinement fails as soon as the lowest barrier is crossed |

**Accepted alternative:** a plain adaptive `scipy.integrate.quad` directly in u, using correct turning points, passes all 4 step-2 tests. That's expected: it's a correct method. The tests don't impose the gold solution's substitution.

### 3b. Errors found while building the task

| What happened | Why it matters |
|---|---|
| **Bug in the first gold prototype: skipped barrier.** The first `turning_points` marched with a fixed 10⁻² step. At h = 1 − 10⁻⁶ the barrier around π (2.8·10⁻³ wide) fell between two steps, and the period came out as 241 instead of 4K = 34.56. The exact elliptic-integral value exposed it. Fixed by detecting sign changes of f' inside each step. That case is now a hidden test. | Without an **independent** oracle (the exact pendulum period), this wrong value would have been stored as the expected answer, and correct code would have been graded as wrong (deck, slide 8). It's the same defect we found in SciCode task 9.1. |
| **Loss of precision at the turning points.** h − f(u) is a difference of two nearly equal numbers close to a turning point. That cost about 10⁻⁸ in relative accuracy and made `quad` emit roundoff warnings near the separatrix. Fixed by computing f(end) − f(u) with Simpson's rule on f' when \|u − end\| < 10⁻³. | The final accuracy is about 10⁻¹⁵ in normal cases and 7·10⁻¹² at h = 1 − 10⁻⁶, which leaves a safe margin against the 10⁻⁵ test tolerance. |
| **Bug in my first ODE oracle.** Starting at b with y = +10⁻³⁰⁰, it detected a spurious event at t ≈ 0 and returned half the period. Fixed by starting at rest at a and ignoring events at t ≈ 0. | Oracles need checking too. A factor of exactly 2 is a typical sign of a counting error, not a physics error. |
| **Bug in the first error-control harness.** The test namespace didn't contain `np`, so every test raised an error, **including the gold solution's**, and all the controls looked like they "failed". Fixed, and the script now refuses to report success unless the gold solution passes every test in the same harness. | A negative control only proves something if the gold solution passes the same test. Otherwise "it fails" means nothing. |

## 4. Why the gold solution is correct

`validate.py` compares the gold solution with sources that share none of its code:

1. **Closed forms.** Harmonic turning points ±√(h/2), pendulum turning points ±arccos(−h), pendulum barriers [−π, π, 1] (IV.A), and tilted-pendulum barriers [−π + asin 0.2, π + asin 0.2]. Agreement is about 10⁻¹⁵.
2. **Exact pendulum period** T = 4K((1 + h)/2), from `scipy.special.ellipk`. The relative error is 2·10⁻¹⁶ at h = 0.5, 7·10⁻¹² at h = 1 − 10⁻⁶ and 6·10⁻⁹ at h = 1 − 10⁻⁹.
3. **Direct ODE integration** (DOP853, rtol 10⁻¹³) of x'' = −f'(x), started at rest at a. The return time matches the period to 10⁻¹¹–10⁻¹⁵, and the orbit reaches b to 10⁻¹³. This confirms there is no hidden barrier in the interval.
4. **Exam results.**
   - The harmonic period 2π/ω does not depend on h (I.A.1).
   - The small-oscillation limit T → 2π/√1.1 for the perturbed pendulum is reproduced to 10⁻⁹ (III.E.3).
   - Near the separatrix, T ≈ 2 ln(32/(1 − h)) (IV.D).
   - The nearest equilibrium to the right of 0 for the perturbed pendulum, from an independent Newton iteration, is e = 3.49906 and it is unstable, with f''(e) = −0.84 (IV.E.1 and IV.E.3).

**Tolerances:** the default `np.allclose` (rtol 10⁻⁵, atol 10⁻⁸), which leaves a wide margin for any correct method.

**Note from the exam's jury report:** more than 50% of candidates wrote the stability condition as f''(e) ≥ 0, when it must be f''(e) > 0. The task avoids that degenerate case: every equilibrium in the tests has f''(e) ≠ 0.

## 5. Points for the SME reviewer

1. Is h = 1 − 10⁻⁶ a fair "hard regime", or is it too close to the separatrix for a benchmark?
2. Should f' be given as an input (as now), or should the model have to differentiate f itself?
3. Step 1 assumes the turning points are simple (f'(a), f'(b) ≠ 0). Should that assumption be stated even more explicitly in the prompt?
4. Domain: Physics → Computational Physics, or Mathematics → Computational Mechanics?

## 6. Files and how to reproduce

| File | Role |
|---|---|
| `TASK.md` | This explanation |
| `problem.json` | The task in SciCode format: prompts, function headers, test code |
| `gold/` | Reference solution, one file per step |
| `test_data.h5` | Expected values for the hidden tests |
| `build_task.py` | Regenerates `problem.json` and `test_data.h5` from `gold/` |
| `run_tests.py` | Runs the tests the way SciCode's test script does |
| `validate.py` | Independent checks and error controls |

```bash
source ../../.venv/Scripts/activate      # SciCode venv (needs numpy, scipy, h5py, scicode)
python build_task.py                     # regenerates problem.json and test_data.h5 from gold/
python run_tests.py --gold               # gold passes 4/4 (3 steps + general test)
python validate.py                       # independent checks + error controls
python run_tests.py --code-dir <dir>     # test candidate files SME02.1.py, SME02.2.py, SME02.3.py
```

Note: the official `eval/scripts/test_generated_code.py` loads tasks from Hugging Face, so it can't see this new task. `run_tests.py` assembles the files the same way (code + `process_hdf5_to_tuple` + test cases) and runs them with `sys.executable`.
