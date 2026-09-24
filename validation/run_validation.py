#!/usr/bin/env python3
"""Grade the reference solution, an honest variant and plausible wrong answers with
the task's REAL verifier code (task/tests/test_outputs.py), without Docker.

  * the reference and the honest variant must score reward 1
  * every wrong answer must score reward 0
  * a shortcut the verifier is known NOT to catch is listed as such and asserted
    to pass, so the gap is recorded rather than hidden

The reference solution is executed unmodified except for its /app paths, which are
pointed at a scratch directory. Usage:  python3 validation/run_validation.py
(needs numpy and scipy). Exit status is non-zero if any expectation is violated.
"""
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
TASK = os.path.join(os.path.dirname(HERE), "task")
DATA = os.path.join(TASK, "environment", "data")
sys.path[:0] = [os.path.join(TASK, "tests"), HERE]


# The verifier's test functions are called directly, so pytest's decorators are
# replaced by pass-throughs. pytest.fail must still fail, so it raises AssertionError.
def _fail(msg):
    raise AssertionError(msg)


stub = types.ModuleType("pytest")
stub.fixture = lambda *a, **k: (lambda f: f)
stub.mark = types.SimpleNamespace(parametrize=lambda *a, **k: (lambda f: f))
stub.fail = _fail
sys.modules["pytest"] = stub

import test_outputs as V   # noqa: E402  (the verifier)
import variants            # noqa: E402

TESTS = [V.test_answer_has_all_57_rows, V.test_concentrations_positive,
         V.test_concentrations_within_tolerance]


def write(outdir, answers, F):
    with open(os.path.join(outdir, "answer.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["site", "compound", "cw_ng_per_L"])
        for (s, c), v in answers.items():
            w.writerow([s, c, f"{v:.6g}"])
    with open(os.path.join(outdir, "sampling_rates.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["site", "F"])
        for s, v in F.items():
            w.writerow([s, f"{v:.6g}"])


def run_reference(outdir):
    src = open(os.path.join(TASK, "solution", "reference_solution.py")).read()
    shutil.copytree(DATA, os.path.join(outdir, "data"))
    src = src.replace('"/app/', f'"{outdir}/').replace('"/app"', f'"{outdir}"')
    script = os.path.join(outdir, "ref.py")
    open(script, "w").write(src)
    subprocess.run([sys.executable, script], check=True, capture_output=True)


def grade(outdir):
    V.ANSWER = os.path.join(outdir, "answer.csv")
    V.RATES = os.path.join(outdir, "sampling_rates.csv")
    failures, within = [], None
    try:
        ans = V.answers()
        ok = sum(abs(ans[(s, c)] / t - 1) <= V.TOL["cw_rel"]
                 for s, d in V.KEY["cw_ng_per_L"].items() for c, t in d.items() if (s, c) in ans)
        within = f"{ok}/{V.TOL['cw_total']}"
        for t in TESTS:
            try:
                t(ans)
            except AssertionError:
                failures.append(t.__name__.replace("test_", ""))
    except AssertionError:
        failures.append("answer.csv unreadable")
    try:
        V.test_sampling_rates(V.rates())
    except AssertionError:
        failures.append("sampling_rates")
    return int(not failures), within, failures


def scaled(fn, conc_factor, f_factor):
    def inner(data):
        a, F = fn(data)
        return {k: v * conc_factor for k, v in a.items()}, {k: v * f_factor for k, v in F.items()}
    return inner


def constant(data):
    a, F = variants.honest_disk_with_binding_term(data)
    med = sorted(a.values())[len(a) // 2]
    return {k: med for k in a}, F


CASES = [
    ("reference (solution/reference_solution.py)", None, 1),
    ("honest: sheets + disk with fitted binding term", variants.honest_disk_with_binding_term, 1),
    ("wrong: one Speedisk G per site (textbook)", variants.wrong_textbook_one_G, 0),
    ("wrong: G from least hydrophobic compound", variants.wrong_least_hydrophobic_calibrant, 0),
    ("wrong: steady-state inversion of 30-day set", variants.wrong_steady_state_30day_set, 0),
    # Documented gap, asserted rather than hidden: the blanks are small next to the
    # accumulated masses, so skipping blank correction stays inside both bands.
    ("NOT CAUGHT: procedural blanks ignored", variants.wrong_no_blank_correction, 1),
    ("wrong: sheet mass in g, not kg (F x1000)",
     scaled(variants.honest_disk_with_binding_term, 1e-3, 1e3), 0),
    ("wrong: one constant for every value", constant, 0),
]


def main():
    ok_all = True
    print(f"{'submission':48s} {'reward':>6s} {'expect':>6s} {'within 20%':>10s}  failed")
    print("-" * 104)
    for name, fn, expected in CASES:
        with tempfile.TemporaryDirectory() as d:
            if fn is None:
                run_reference(d)
            else:
                write(d, *fn(DATA))
            reward, within, failures = grade(d)
        flag = "" if reward == expected else "   <-- UNEXPECTED"
        ok_all &= reward == expected
        print(f"{name:48s} {reward:>6d} {expected:>6d} {within or '-':>10s}  "
              f"{', '.join(failures) or '-'}{flag}")
    print("-" * 104)
    print("ALL EXPECTATIONS MET" if ok_all else "EXPECTATION VIOLATED")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
