"""Materialize explicitly grouped technical BAM replicates in assembly scratch."""

from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

from ..utils.exceptions import IPyradError
from ..utils.parallel import run_pipeline
from .read_filters import BIN_SAM


def _read_group_ids(header: str) -> list[str]:
    ids = []
    for line in header.splitlines():
        if line.startswith("@RG\t"):
            tags = dict(field.split(":", 1) for field in line.split("\t")[1:])
            if not tags.get("ID") or tags["ID"] in ids:
                raise IPyradError("Invalid or duplicate @RG ID within a grouped input BAM")
            ids.append(tags["ID"])
    return ids


def _member_configurations(headers, sample):
    """Preserve unique RG IDs, disambiguate collisions, and handle missing RG."""
    original = [_read_group_ids(header) for header in headers]
    counts = Counter(rg for ids in original for rg in ids)
    reserved = set(counts)

    def allocate(base):
        candidate = base
        i = 0
        while candidate in reserved:
            i += 1
            candidate = f"{base}.{i}"
        reserved.add(candidate)
        return candidate

    configs = []
    for i, ids in enumerate(original):
        mapping = {
            rg: rg if counts[rg] == 1 else allocate(f"ipyrad2.{i}.{rg}")
            for rg in ids
        }
        # An untagged read in a single-RG input has an unambiguous
        # library origin. Otherwise provide a distinct fallback group.
        synthetic = len(mapping) != 1
        configs.append({
            "sample": sample,
            "prefix": f"ipyrad2.{i}:",
            "read_groups": mapping,
            "untagged": allocate(f"ipyrad2.{i}.untagged") if synthetic else next(iter(mapping.values())),
            "synthetic_untagged": synthetic,
        })
    return configs


def rewrite_sam(source, target, config):
    """Normalize temporary RG sample names and isolate QNAMEs across BAMs."""
    groups = config["read_groups"]
    sample = config["sample"]
    prefix = config["prefix"]
    untagged = config["untagged"]
    header_finished = False
    for line in source:
        if not header_finished and line.startswith("@"):
            if line.startswith("@RG\t"):
                fields = line.rstrip("\n").split("\t")
                tags = dict(field.split(":", 1) for field in fields[1:])
                tags["ID"] = groups[tags["ID"]]
                tags["SM"] = sample
                line = "\t".join(["@RG", *(f"{k}:{v}" for k, v in tags.items())]) + "\n"
            target.write(line)
            continue
        if not header_finished:
            if config["synthetic_untagged"]:
                target.write(f"@RG\tID:{untagged}\tSM:{sample}\n")
            header_finished = True
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            raise IPyradError("Malformed SAM record while grouping BAMs")
        # HTSlib overlap correction uses query names; namespace templates by
        # source BAM even when two libraries reused identical QNAMEs.
        name = prefix + "q:" + fields[0]
        if len(name) > 254:
            name = prefix + "h:" + hashlib.sha256(fields[0].encode()).hexdigest()
        fields[0] = name
        tagged = False
        for i, field in enumerate(fields[11:], 11):
            if field.startswith("RG:"):
                if not field.startswith("RG:Z:") or tagged:
                    raise IPyradError("Invalid or duplicate RG tag in grouped input BAM")
                rg = field[5:]
                if rg not in groups:
                    raise IPyradError(f"Read refers to undeclared read group {rg!r}")
                fields[i] = "RG:Z:" + groups[rg]
                tagged = True
        if not tagged:
            fields.append("RG:Z:" + untagged)
        target.write("\t".join(fields) + "\n")
    if not header_finished and config["synthetic_untagged"]:
        target.write(f"@RG\tID:{untagged}\tSM:{sample}\n")


def merge_grouped_bams(*, sname: str, bam_paths: tuple[Path, ...],
                       tmpdir: Path, threads: int) -> Path:
    """Combine coordinate-sorted BAMs, without modifying or deduplicating inputs."""
    grouped = tmpdir / "grouped_bams"
    grouped.mkdir(exist_ok=True)
    key = hashlib.sha256(sname.encode()).hexdigest()[:24]
    out = grouped / f"{key}.bam"
    if out.exists():
        raise IPyradError(f"Grouped BAM already exists: {out}")
    headers = []
    for path in bam_paths:
        run_pipeline([[BIN_SAM, "quickcheck", str(path)]])
        _, header, _ = run_pipeline([[BIN_SAM, "view", "-H", str(path)]])
        headers.append(header.decode())
    configs = _member_configurations(headers, sname)
    try:
        with TemporaryDirectory(prefix=f".{key}.", dir=grouped) as scratch:
            scratch = Path(scratch)
            members = []
            for i, (path, config) in enumerate(zip(bam_paths, configs)):
                member = scratch / f"{i}.bam"
                config_path = scratch / f"{i}.json"
                config_path.write_text(json.dumps(config))
                run_pipeline([
                    [BIN_SAM, "view", "-h", str(path)],
                    [sys.executable, "-m", "ipyrad2.assembler.bam_grouping", str(config_path)],
                    [BIN_SAM, "view", "-b", "-o", str(member), "-"],
                ], ignore_sigpipe=False)
                members.append(member)
            run_pipeline([[BIN_SAM, "merge", "-@", str(max(0, threads - 1)), "-o",
                           str(out), *(str(path) for path in members)]])
        run_pipeline([[BIN_SAM, "quickcheck", str(out)]])
        # Indexing also verifies that records are coordinate sorted.
        run_pipeline([[BIN_SAM, "index", "-c", "-@", str(max(0, threads - 1)), str(out)]])
    except BaseException:
        out.unlink(missing_ok=True)
        Path(str(out) + ".csi").unlink(missing_ok=True)
        raise
    return out


if __name__ == "__main__":
    rewrite_sam(sys.stdin, sys.stdout, json.loads(Path(sys.argv[1]).read_text()))
