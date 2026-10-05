# task-constractor

New SciCode-style tasks (research-coding problems with hidden numerical tests), built and validated in the [SciCode](https://github.com/scicode-bench/SciCode) format.

| Task | Domain | Topic | Status |
|---|---|---|---|
| [SME01_driven_oscillator](SME01_driven_oscillator/TASK.md) | Physics → Computational Physics | Driven oscillator y'' + y = μ(t): periodic response vs resonance | draft, pending SME sign-off |
| [SME02_conservative_oscillator](SME02_conservative_oscillator/TASK.md) | Physics → Computational Physics | Period of x'' = −f'(x) for any potential, up to the separatrix | draft, pending SME sign-off |

Each task folder contains:

| File | Role |
|---|---|
| `TASK.md` | Explanation of the task, the errors it catches and why they matter, and the validation evidence |
| `problem.json` | The task in SciCode format: prompts, function headers, test code |
| `gold/` | Reference solution, one file per step |
| `test_data.h5` | Expected values for the hidden tests |
| `build_task.py` | Regenerates `problem.json` and `test_data.h5` from `gold/` |
| `run_tests.py` | Runs the tests the way SciCode's test script does |
| `validate.py` | Independent checks of the gold solution and negative controls |

## Running

The scripts need a Python environment with the SciCode package installed (`pip install -e .` in a clone of SciCode), plus numpy, scipy and h5py.

```bash
cd SME01_driven_oscillator
python build_task.py           # regenerate problem.json and test_data.h5
python run_tests.py --gold     # the gold solution passes every test
python validate.py             # independent oracles + negative controls
```
