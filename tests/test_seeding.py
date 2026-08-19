"""Seeding must be reproducible across processes.

The pre-specification fixes a seed so that a re-run draws the same random
gene sets. That promise is easy to break invisibly: Python randomizes string
hashing per process, so `hash(cohort)` in a seed expression yields different
draws on every invocation while the code still *looks* seeded.

Run:  pytest tests/test_seeding.py -q
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load_null_script():
    spec = importlib.util.spec_from_file_location(
        "null_calibration", ROOT / "scripts" / "04_null_calibration.py")
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ROOT / "src"))
    spec.loader.exec_module(mod)
    return mod


def test_offset_is_deterministic_within_a_process():
    m = _load_null_script()
    a = [m._stable_offset(c) for c in ("Mariathasan", "Kim", "Gide")]
    b = [m._stable_offset(c) for c in ("Mariathasan", "Kim", "Gide")]
    assert a == b


def test_offset_is_identical_across_processes():
    """The property that hash() fails. Run in fresh interpreters with
    different PYTHONHASHSEED values and require the same answer."""
    code = (
        "import hashlib;"
        "f=lambda n: int.from_bytes("
        "hashlib.blake2b(n.encode('utf-8'),digest_size=8).digest(),'big')%997;"
        "print([f(c) for c in ('Mariathasan','Kim','Gide','Braun')])"
    )
    outs = []
    for hashseed in ("0", "1", "12345"):
        r = subprocess.run([sys.executable, "-c", code], capture_output=True,
                           text=True, check=True,
                           env={"PYTHONHASHSEED": hashseed, "PATH": ""})
        outs.append(r.stdout.strip())
    assert len(set(outs)) == 1, f"offset varies across processes: {outs}"


def test_offset_matches_the_script_implementation():
    """Pins the test's reference formula to the script's own function."""
    import hashlib
    m = _load_null_script()
    for name in ("Mariathasan", "Kim", "Gide", "Braun", "VanDenEnde"):
        expected = int.from_bytes(
            hashlib.blake2b(name.encode("utf-8"), digest_size=8).digest(),
            "big") % 997
        assert m._stable_offset(name) == expected


def test_python_hash_would_have_been_unstable():
    """Documents why blake2b is used -- if this ever starts passing with
    hash(), the rationale in the script can be revisited."""
    code = "print(abs(hash('Mariathasan')) % 997)"
    outs = set()
    for hashseed in ("1", "2", "3"):
        r = subprocess.run([sys.executable, "-c", code], capture_output=True,
                           text=True, check=True,
                           env={"PYTHONHASHSEED": hashseed, "PATH": ""})
        outs.add(r.stdout.strip())
    assert len(outs) > 1, "hash() appears stable here; check PYTHONHASHSEED"


def test_distinct_cohorts_get_distinct_offsets():
    m = _load_null_script()
    cohorts = ["Mariathasan", "Braun", "Liu", "Limagne1", "Puch", "Van_Allen",
               "Gide", "VanDenEnde", "Kim", "Riaz", "Miao1", "Hugo", "Jung",
               "Nathanson", "Limagne2"]
    offs = [m._stable_offset(c) for c in cohorts]
    assert len(set(offs)) == len(offs), "two cohorts share a seed offset"
