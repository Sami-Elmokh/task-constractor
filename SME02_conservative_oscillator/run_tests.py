"""Run SME02 tests the same way eval/scripts/test_generated_code.py does.

Each step file <code-dir>/SME02.<k>.py must contain everything needed for step k
(imports + functions of steps 1..k), exactly like files written by eval/scripts/gencode.py.

    python run_tests.py --gold                      # test the gold solution from problem.json
    python run_tests.py --code-dir path/to/files    # test candidate files (results also saved to --log-dir)
"""
import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
H5 = (HERE / "test_data.h5").resolve().as_posix()


def write_gold_files(problem, out_dir):
    out_dir.mkdir(exist_ok=True)
    code = problem["required_dependencies"] + "\n\n"
    for step in problem["sub_steps"]:
        code += step["ground_truth_code"] + "\n\n"
        (out_dir / f"{step['step_number']}.py").write_text(code, encoding="utf-8")


def assemble(code, step_id, tests):
    text = code + "\n\nfrom scicode.parse.parse import process_hdf5_to_tuple\n\n"
    text += f"targets = process_hdf5_to_tuple('{step_id}', {len(tests)}, h5py_file='{H5}')\n"
    for i, tc in enumerate(tests):
        text += f"target = targets[{i}]\n\n{tc}\n"
    return text


def run(path):
    r = subprocess.run([sys.executable, str(path)], capture_output=True, text=True, timeout=1800)
    return r.returncode == 0, r.stderr.strip().splitlines()[-1:] if r.returncode else []


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--code-dir", type=Path, help="folder with SME02.1.py, SME02.2.py, SME02.3.py")
    p.add_argument("--gold", action="store_true", help="test the gold solution (built in a temp folder)")
    p.add_argument("--log-dir", type=Path, help="optional folder for pass/fail files")
    args = p.parse_args()
    if args.gold == bool(args.code_dir):
        p.error("give exactly one of --gold or --code-dir")

    problem = json.loads((HERE / "problem.json").read_text(encoding="utf-8"))
    checks = [(s["step_number"], s["test_cases"], s["step_number"]) for s in problem["sub_steps"]]
    last = problem["sub_steps"][-1]["step_number"]
    checks.append((problem["problem_id"], problem["general_tests"], last))  # general test uses the full code
    if args.log_dir:
        args.log_dir.mkdir(parents=True, exist_ok=True)

    n_pass = 0
    with tempfile.TemporaryDirectory() as tmp:
        code_dir = args.code_dir
        if args.gold:
            code_dir = Path(tmp, "gold_code")
            write_gold_files(problem, code_dir)
        for target_id, tests, file_id in checks:
            code_file = code_dir / f"{file_id}.py"
            if not code_file.exists():
                print(f"{target_id}: missing {code_file.name}")
                continue
            script = Path(tmp, f"{target_id}.py")
            script.write_text(assemble(code_file.read_text(encoding="utf-8"), target_id, tests),
                              encoding="utf-8")
            ok, err = run(script)
            n_pass += ok
            result = "pass" if ok else "fail"
            if args.log_dir:
                (args.log_dir / f"{target_id}.txt").write_text(result, encoding="utf-8")
            print(f"{target_id}: {result}" + (f"  ({err[0]})" if err else ""))
    print(f"{n_pass}/{len(checks)} passed")


if __name__ == "__main__":
    main()
