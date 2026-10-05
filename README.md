# task-constractor

Original SciCode-style scientific coding tasks. Each task is a folder in the SME task format.

| Task | Domain | Topic |
|---|---|---|
| [driven_oscillator_resonance](driven_oscillator_resonance/problem.md) | Physics → Computational Physics | Exact response of an undamped oscillator to a sampled periodic force, including resonance |
| [oscillator_period_energy](oscillator_period_energy/problem.md) | Physics → Computational Physics | Period–energy curve of a 1D potential well, from small oscillations up to escape |

Each task folder contains:

| File | Role |
|---|---|
| `problem.md` | The complete written problem (no solution, tests or answers) |
| `solution.py` | Complete working code |
| `main.py` | Runs the solution on one example |
| `solution_test.py` | pytest tests with independently verified expected answers |
| `requirements.txt` | Libraries used |
| `verification.md` | How the science and expected answers were checked, and why the tolerances were chosen |

## Running a task

```bash
cd driven_oscillator_resonance
python -m pip install -r requirements.txt
python main.py
python -m pytest -q
```
