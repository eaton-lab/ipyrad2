"""Exercise trimming from an installed wheel, where editable tests miss data loss."""
import gzip
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import zipfile


def test_installed_wheel_contains_adapters_and_trims_paired_reads(tmp_path):
    root = Path(__file__).resolve().parents[1]
    source = tmp_path / "source"
    source.mkdir()
    shutil.copytree(root / "ipyrad2", source / "ipyrad2",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    for name in ("pyproject.toml", "README.md", "LICENSE.txt"):
        if (root / name).exists():
            shutil.copy2(root / name, source / name)
    env = dict(os.environ, SETUPTOOLS_SCM_PRETEND_VERSION_FOR_IPYRAD2="0.0.0",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    subprocess.run([sys.executable, "-c",
                    "from setuptools.build_meta import build_wheel; build_wheel('dist')"],
                   cwd=source, env=env, check=True, capture_output=True, text=True)
    wheel, = (source / "dist").glob("*.whl")
    with zipfile.ZipFile(wheel) as archive:
        assert archive.read("ipyrad2/trimmer/adapters.fa") == (
            root / "ipyrad2/trimmer/adapters.fa").read_bytes()
        installed = tmp_path / "installed"
        archive.extractall(installed)
    # Start outside the source tree and force imports from the built artifact.
    env["PYTHONPATH"] = str(installed)
    rng = random.Random(23)
    inputs = []
    for mate, motif in ((1, "ATGCC"), (2, "CTGCAA")):
        path = tmp_path / f"sample_R{mate}.fastq.gz"
        with gzip.open(path, "wt") as out:
            for idx in range(20):
                seq = motif + "".join(rng.choices("ACGT", k=90))
                out.write(f"@read{idx}\n{seq}\n+\n{'I' * len(seq)}\n")
        inputs.append(str(path))
    output = tmp_path / "trimmed"
    command = (
        "import ipyrad2; from pathlib import Path; "
        "assert Path(ipyrad2.__file__).is_relative_to(Path('installed').resolve()); "
        "from ipyrad2.cli import main; main()"
    )
    result = subprocess.run(
        [sys.executable, "-c", command, "trim", "-d", *inputs, "-o", str(output),
         "-dx", "_R", "-di", "1", "-e1", "ATGCC,ATGCAT", "-e2", "CTGCAA",
         "-c", "1", "-t", "1"],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for mate in (1, 2):
        with gzip.open(output / f"sample.R{mate}.trimmed.fastq.gz", "rt") as inp:
            lines = inp.read().splitlines()
        assert len(lines) == 20 * 4
        assert all(len(lines[idx + 1]) == (89 if mate == 1 else 90) for idx in range(0, len(lines), 4))
