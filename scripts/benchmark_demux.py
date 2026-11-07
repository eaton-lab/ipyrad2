#!/usr/bin/env python
"""Compare demux repositories with full-record multiset fingerprints and resources.

Requires psutil in the benchmark environment (not a runtime dependency).
Input JSON: [{"name": "rad", "kwargs": {"fastqs": [...], "barcodes": "...",
"cutsite_1": "TGCAG", "cutsite_2": null, "max_reads": null}}].
--synthetic covers 96/384 samples, SE/PE, one/two lanes, balanced/skewed reads.
Timing excludes fixture generation and verification. RSS sums process RSS and
therefore double-counts shared pages. Disk peaks are sampled lower bounds.
"""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
import os
import platform
from pathlib import Path
import random
import sqlite3
import subprocess
import sys
import time


def worker(repo, config_path):
    sys.path.insert(0, str(Path(repo).resolve()))
    from ipyrad2.demuxer.demux import Demux

    kw = json.loads(Path(config_path).read_text())
    kw["fastqs"] = [Path(p) for p in kw["fastqs"]]
    for key in ("barcodes", "outdir", "tmpdir"):
        if kw.get(key) is not None:
            kw[key] = Path(kw[key])
    Demux(**kw).run()


def synthetic(root, nreads):
    rng = random.Random(914)
    barcodes = set()
    while len(barcodes) < 384:
        code = "".join(rng.choice("ACGT") for _ in range(10))
        if len(set(code)) > 1:
            barcodes.add(code)
    codes = sorted(barcodes)
    cases = []
    for nsamples in (96, 384):
        for paired in (False, True):
            for lanes in (1, 2):
                for skew in (False, True):
                    name = f'synth-{nsamples}-{"PE" if paired else "SE"}-{lanes}lane-{"skew" if skew else "balanced"}'
                    folder = root / name
                    folder.mkdir(parents=True, exist_ok=False)
                    bars = folder / "barcodes.tsv"
                    bars.write_text(
                        "".join(
                            f"s{i:03d}\t{code}\n"
                            for i, code in enumerate(codes[:nsamples])
                        )
                    )
                    paths = []
                    for lane in range(lanes):
                        r1 = folder / f"lane{lane}_R1.fastq.gz"
                        r2 = folder / f"lane{lane}_R2.fastq.gz"
                        paths.append(str(r1))
                        if paired:
                            paths.append(str(r2))
                        with gzip.open(r1, "wb", compresslevel=1) as a:
                            b = gzip.open(r2, "wb", compresslevel=1) if paired else None
                            try:
                                for idx in range(nreads):
                                    sample = 0 if skew and idx % 10 else idx % nsamples
                                    seq = (
                                        codes[sample]
                                        + "ATCGG"
                                        + "".join(rng.choice("ACGT") for _ in range(85))
                                    )
                                    header = f"@lane{lane}-read{idx}"
                                    a.write(
                                        f'{header}/1\n{seq}\n+\n{"I"*len(seq)}\n'.encode()
                                    )
                                    if b is not None:
                                        seq2 = "".join(
                                            rng.choice("ACGT") for _ in range(100)
                                        )
                                        b.write(
                                            f'{header}/2\n{seq2}\n+\n{"I"*100}\n'.encode()
                                        )
                            finally:
                                if b is not None:
                                    b.close()
                    cases.append(
                        {
                            "name": name,
                            "kwargs": dict(
                                fastqs=paths,
                                barcodes=str(bars),
                                cutsite_1="ATCGG",
                                cutsite_2=None,
                            ),
                        }
                    )
    return cases


def quartet(handle):
    first = handle.readline()
    if not first:
        return None
    lines = (first, handle.readline(), handle.readline(), handle.readline())
    if any(not line for line in lines):
        raise ValueError("incomplete output FASTQ")
    return b"".join(line.rstrip(b"\r\n") + b"\n" for line in lines)


def fingerprint(outdir):
    """Disk-backed sorted SHA256 record counts, including pair association/multiplicity."""
    database = outdir / "verification.sqlite"
    db = sqlite3.connect(database)
    db.execute(
        "create table records (digest blob primary key, n integer) without rowid"
    )
    summary = {}
    try:
        for r1 in sorted(outdir.glob("*_R1.fastq.gz")):
            r2 = r1.with_name(r1.name.replace("_R1.fastq.gz", "_R2.fastq.gz"))
            db.execute("delete from records")
            counts = Counter()
            total = 0
            with gzip.open(r1, "rb") as a:
                b = gzip.open(r2, "rb") if r2.exists() else None
                try:
                    while True:
                        left = quartet(a)
                        right = quartet(b) if b else None
                        if b and (left is None) != (right is None):
                            raise ValueError(f"unequal output mate counts: {r1}")
                        if left is None:
                            break
                        if right is not None:
                            ids = [
                                read.split(b"\n", 1)[0]
                                .split()[0]
                                .removesuffix(b"/1")
                                .removesuffix(b"/2")
                                for read in (left, right)
                            ]
                            if ids[0] != ids[1]:
                                raise ValueError(f"output mate IDs differ: {r1}")
                        payload = len(left).to_bytes(8, "big") + left + (right or b"")
                        counts[hashlib.sha256(payload).digest()] += 1
                        total += 1
                        if total % 10000 == 0:
                            db.executemany(
                                "insert into records values (?, ?) on conflict(digest) do update set n=n+excluded.n",
                                counts.items(),
                            )
                            counts.clear()
                    db.executemany(
                        "insert into records values (?, ?) on conflict(digest) do update set n=n+excluded.n",
                        counts.items(),
                    )
                    db.commit()
                finally:
                    if b:
                        b.close()
            digest = hashlib.sha256()
            for key, count in db.execute(
                "select digest,n from records order by digest"
            ):
                digest.update(key)
                digest.update(count.to_bytes(8, "big"))
            summary[r1.name] = {
                "records": total,
                "paired": r2.exists(),
                "multiset_sha256": digest.hexdigest(),
            }
    finally:
        db.close()
        database.unlink()
    return summary


def tree_sizes(folder):
    total = temporary = 0
    for base, _, files in os.walk(folder):
        for filename in files:
            try:
                size = (Path(base) / filename).stat().st_size
            except FileNotFoundError:
                continue
            total += size
            if filename.endswith((".fastq", ".bin")):
                temporary += size
    return total, temporary


def measure(repo, kw, folder, scratch, timeout):
    import psutil

    folder.mkdir(parents=True)
    config_path = folder / "config.json"
    config_path.write_text(json.dumps(kw, indent=2))
    cmd = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        str(repo),
        str(config_path),
    ]
    peak_rss = peak_disk = peak_temp = 0
    io_by_pid = {}
    started = time.monotonic()
    with (folder / "run.log").open("w") as log:
        proc = subprocess.Popen(
            cmd, cwd=repo, stdout=log, stderr=log, start_new_session=True
        )
        parent = psutil.Process(proc.pid)
        last_disk = 0
        try:
            while proc.poll() is None:
                if time.monotonic() - started > timeout:
                    raise TimeoutError(f"benchmark exceeded {timeout}s: {folder}")
                try:
                    processes = [parent] + parent.children(recursive=True)
                except psutil.NoSuchProcess:
                    processes = []
                rss = 0
                for item in processes:
                    try:
                        rss += item.memory_info().rss
                        io_counts = item.io_counters()
                        io_by_pid[item.pid] = (
                            io_counts.read_bytes,
                            io_counts.write_bytes,
                        )
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                peak_rss = max(peak_rss, rss)
                if time.monotonic() - last_disk >= 0.2:
                    size, temp = tree_sizes(Path(kw["outdir"]))
                    if scratch:
                        more_size, more_temp = tree_sizes(scratch)
                        size += more_size
                        temp += more_temp
                    peak_disk = max(peak_disk, size)
                    peak_temp = max(peak_temp, temp)
                    last_disk = time.monotonic()
                time.sleep(0.025)
        except BaseException:
            import signal

            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(3)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            raise
    elapsed = time.monotonic() - started
    peak_disk = max(peak_disk, tree_sizes(Path(kw["outdir"]))[0])
    result = dict(
        seconds=elapsed,
        exit_code=proc.returncode,
        peak_aggregate_rss_bytes=peak_rss,
        peak_output_disk_bytes=peak_disk,
        peak_temporary_fastq_bytes=peak_temp,
        sampled_read_bytes=sum(v[0] for v in io_by_pid.values()),
        sampled_write_bytes=sum(v[1] for v in io_by_pid.values()),
        command=cmd,
    )
    if proc.returncode == 0:
        result["outputs"] = fingerprint(Path(kw["outdir"]))
        reports = list(Path(kw["outdir"]).glob("ipyrad_demux_stats_*.json"))
        result["statistics"] = json.loads(reports[0].read_text())
    return result


def source_fingerprint(repo):
    return {
        str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted((repo / "ipyrad2" / "demuxer").glob("*.py"))
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", nargs=2, metavar=("REPO", "CONFIG"))
    parser.add_argument("--baseline", type=Path)
    parser.add_argument(
        "--candidate", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--out", type=Path)
    parser.add_argument("--cases", type=Path)
    parser.add_argument("--synthetic", action="store_true")
    parser.add_argument("--reads", type=int, default=100000)
    parser.add_argument("--cores", nargs="+", type=int, default=[1, 2, 4, 8])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--discard-fastqs", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(*args.worker)
        return
    if not args.baseline or not args.out or not (args.cases or args.synthetic):
        parser.error("--baseline, --out, and --cases or --synthetic are required")
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=False)
    metadata = dict(
        python=sys.version,
        platform=platform.platform(),
        cpu_count=os.cpu_count(),
        baseline=str(args.baseline.resolve()),
        candidate=str(args.candidate.resolve()),
        baseline_sources=source_fingerprint(args.baseline.resolve()),
        candidate_sources=source_fingerprint(args.candidate.resolve()),
        metrics_note="25 ms sampled summed process RSS (shared pages counted repeatedly); 200 ms sampled disk peaks; timing excludes verification. No cache dropping.",
    )
    if hasattr(os, "sched_getaffinity"):
        metadata["cpu_affinity"] = sorted(os.sched_getaffinity(0))
    (args.out / "metadata.json").write_text(json.dumps(metadata, indent=2))
    cases = json.loads(args.cases.read_text()) if args.cases else []
    if args.synthetic:
        cases += synthetic(args.out / "fixtures", args.reads)
    (args.out / "cases.json").write_text(json.dumps(cases, indent=2))
    results = []
    defaults = dict(
        cutsite_1=None,
        cutsite_2=None,
        max_mismatch=0,
        chunksize=10_000_000,
        merge_technical_replicates=False,
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=None,
        max_reads_kmer=100_000,
        log_level="WARNING",
    )
    for case in cases:
        for cores in args.cores:
            for repeat in range(args.repeats):
                pair = {}
                versions = [("baseline", args.baseline), ("candidate", args.candidate)]
                if repeat % 2:
                    versions.reverse()
                for version, repo in versions:
                    folder = args.out / case["name"] / f"c{cores}-r{repeat}" / version
                    kw = dict(defaults, **case["kwargs"])
                    kw.update(outdir=str(folder / "outputs"), cores=cores)
                    print(
                        f"{case['name']} c{cores} repeat {repeat}: {version}",
                        flush=True,
                    )
                    pair[version] = measure(
                        repo.resolve(), kw, folder, None, args.timeout
                    )
                old, new = pair["baseline"], pair["candidate"]
                result = dict(case=case["name"], cores=cores, repeat=repeat, **pair)
                result["outputs_equal"] = (
                    old.get("outputs") == new.get("outputs")
                    if old["exit_code"] == new["exit_code"] == 0
                    else False
                )
                result["statistics_equal"] = (
                    old.get("statistics") == new.get("statistics")
                    if old["exit_code"] == new["exit_code"] == 0
                    else False
                )
                result["runtime_ratio"] = new["seconds"] / old["seconds"]
                results.append(result)
                (args.out / "results.json").write_text(json.dumps(results, indent=2))
                print(
                    f"  output equality={result['outputs_equal']}, stats equality={result['statistics_equal']}, time ratio={result['runtime_ratio']:.3f}",
                    flush=True,
                )
                if args.discard_fastqs and result["outputs_equal"]:
                    for version in ("baseline", "candidate"):
                        for path in (
                            args.out
                            / case["name"]
                            / f"c{cores}-r{repeat}"
                            / version
                            / "outputs"
                        ).glob("*.fastq.gz"):
                            path.unlink()
    if source_fingerprint(args.candidate.resolve()) != metadata["candidate_sources"]:
        raise RuntimeError(
            "Candidate sources changed during benchmarking; rerun from a frozen snapshot"
        )
    if any(not r["outputs_equal"] or not r["statistics_equal"] for r in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
