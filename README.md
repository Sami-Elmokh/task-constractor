# task-constractor

Original SciCode-style scientific coding tasks. Each task is a folder in the SME task format.

| Task | Domain | Topic |
|---|---|---|
| [driven_oscillator_resonance](driven_oscillator_resonance/problem.md) | Physics → Computational Physics | Response of an undamped oscillator of any natural frequency to a sampled periodic force, at and near resonance |
| [oscillator_period_energy](oscillator_period_energy/problem.md) | Physics → Computational Physics | Period–energy curve of a 1D potential well, from small oscillations up to escape |
| [rectangular_pair_normal_form](rectangular_pair_normal_form/problem.md) | Mathematics → Numerical Linear Algebra | Normal form of a pair of integer matrices under simultaneous equivalence (chains and Jordan blocks), with exact ranks |

Each task folder contains:

| File | Role |
|---|---|
| `problem.md` | The complete written problem (no solution, tests or answers) |
| `solution.py` | Complete working code |
| `main.py` | Runs the solution on one example |
| `solution_test.py` | pytest tests with independently verified expected answers |
| `requirements.txt` | Libraries used |
| `verification.md` | How the science and expected answers were checked, and why the tolerances were chosen |

Where ChatGPT fails on each task: [HISTORIQUE.md](HISTORIQUE.md) (in French). Detailed analysis for the last task: [ERREURS_CHATGPT_rectangular_pair.md](ERREURS_CHATGPT_rectangular_pair.md).

## Running a task

```bash
cd driven_oscillator_resonance
python -m pip install -r requirements.txt
python main.py
python -m pytest -q
```
