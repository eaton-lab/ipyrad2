"""Reproducible, bounded-memory FASTQ sampling before denovo clustering."""

from contextlib import ExitStack, contextmanager
import gzip
from pathlib import Path
from random import Random
from tempfile import TemporaryDirectory

from ..demuxer.match import fastq_read_id, read_fastq_record
from ..utils.exceptions import IPyradError

SAMPLING_SEED = 42


def _open_input(path: Path):
    return (gzip.open if path.suffix == ".gz" else open)(path, "rb")


def _count_records(path: Path) -> int:
    """Count four-line records in a buffered pass, including concatenated gzip members."""
    lines = 0
    last = b"\n"
    with _open_input(path) as handle:
        for block in iter(lambda: handle.read(1024**2), b""):
            lines += block.count(b"\n")
            last = block[-1:]
    lines += last != b"\n"
    if lines % 4:
        raise IPyradError(f"FASTQ {path}: incomplete four-line record")
    return lines // 4


@contextmanager
def sampled_fastqs(sample: str, r1: Path, r2: Path | None, cap: int, directory: Path):
    """Yield the original inputs or temporary FASTQs sampled in lockstep."""
    total = _count_records(r1)
    counts = {"input": total, "retained": min(total, cap)}
    if total <= cap:
        if r2 is not None and _count_records(r2) != total:
            raise IPyradError(f"Sample {sample}: FASTQ mates have different record counts")
        yield r1, r2, counts
        return

    rng = Random(f"{SAMPLING_SEED}:{sample}")
    needed = cap
    with TemporaryDirectory(dir=directory) as scratch:
        selected = Path(scratch) / "R1.fastq.gz"
        selected_mate = Path(scratch) / "R2.fastq.gz" if r2 is not None else None
        with ExitStack() as stack:
            inp = stack.enter_context(_open_input(r1))
            out = stack.enter_context(gzip.open(selected, "wb", compresslevel=1))
            if r2 is not None:
                mate = stack.enter_context(_open_input(r2))
                out_mate = stack.enter_context(gzip.open(selected_mate, "wb", compresslevel=1))
            for offset in range(total):
                index = offset + 1
                read = read_fastq_record(inp, r1, index)
                if read is None:
                    raise IPyradError(f"FASTQ {r1}: input changed during sampling")
                if r2 is not None:
                    other = read_fastq_record(mate, r2, index)
                    if other is None:
                        raise IPyradError(f"Sample {sample}: FASTQ mates have different record counts")
                    if fastq_read_id(read[0], r1, index, 1) != fastq_read_id(other[0], r2, index, 2):
                        raise IPyradError(f"Sample {sample}, record {index}: FASTQ mate IDs differ")
                # With k remaining selections and n remaining records, select with probability k/n.
                if needed and rng.randrange(total - offset) < needed:
                    out.write(b"".join(read))
                    if r2 is not None:
                        out_mate.write(b"".join(other))
                    needed -= 1
            if inp.read(1):
                raise IPyradError(f"FASTQ {r1}: input changed during sampling")
            if r2 is not None and mate.read(1):
                raise IPyradError(f"Sample {sample}: FASTQ mates have different record counts")
        yield selected, selected_mate, counts
