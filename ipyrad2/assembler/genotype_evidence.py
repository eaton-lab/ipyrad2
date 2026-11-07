"""Write native caller evidence for exactly matched, retained SNP records."""

from itertools import groupby
from pathlib import Path
from tempfile import TemporaryDirectory

from loguru import logger

from ..utils.exceptions import IPyradError
from ..utils.parallel import run_pipeline, stream_pipeline_lines
from .variants import BIN_BCF, _load_mask_entries_from_paths


_EVIDENCE_FIELDS = ("DP", "AD", "GQ", "PL")
_FT_HEADER = '##FORMAT=<ID=FT,Number=1,Type=String,Description="Evidence status: PASS, NoHardCall, NoRead or HardMask">'


def _read_header(lines):
    meta = []
    for line in lines:
        if line.startswith("##"):
            meta.append(line)
        elif line.startswith("#CHROM\t"):
            return meta, line
        else:
            break
    raise IPyradError("Genotype evidence requires a valid VCF header.")


def _coordinate_groups(lines, contigs):
    previous = None

    def coordinate(fields):
        return contigs[fields[0]], int(fields[1])

    for key, group in groupby((line.split("\t") for line in lines), coordinate):
        if previous is not None and key <= previous:
            raise IPyradError("Genotype evidence input records are out of coordinate order.")
        previous = key
        yield key, list(group)


def _masked_fields(entries, pos0):
    for entry in entries:
        intervals = entry["intervals"]
        cursor = entry["cursor"]
        while cursor < len(intervals) and intervals[cursor][1] <= pos0:
            cursor += 1
        entry["cursor"] = cursor
        if cursor < len(intervals) and intervals[cursor][0] <= pos0:
            yield entry["field_idx"]


def _evidence_sample(native, called, allele_count):
    dp, ad, gq, pl = (native.get(key, ".") for key in _EVIDENCE_FIELDS)
    if dp in (".", "0"):
        return ("./.", dp, ".", ".", ".", "NoRead")
    for key, value, width in (
        ("AD", ad, allele_count), ("PL", pl, allele_count * (allele_count + 1) // 2)
    ):
        if value != "." and len(value.split(",")) != width:
            raise IPyradError(f"Native {key} has the wrong cardinality for the matched SNP alleles.")
    gt = called.get("GT", "./.")
    return (gt, dp, ad, gq, pl, "NoHardCall" if "." in gt else "PASS")


def write_genotype_evidence(*, raw_vcf, final_vcf, out_bcf, hard_masks=(), shared_mask=None):
    """Join pre-normalization evidence to the canonical SNP calls by exact allele key."""
    out_bcf = Path(out_bcf)
    out_index = out_bcf.with_suffix(".bcf.csi")
    if out_bcf.exists() or out_index.exists():
        raise IPyradError(f"Genotype evidence output already exists: {out_bcf}")
    with TemporaryDirectory(prefix=".genotype-evidence-", dir=out_bcf.parent) as temp:
        plain = Path(temp) / "evidence.vcf"
        with (
            stream_pipeline_lines([[BIN_BCF, "view", "--no-version", str(raw_vcf)]]) as raw,
            stream_pipeline_lines([[BIN_BCF, "view", "--no-version", str(final_vcf)]]) as final,
            plain.open("w", encoding="utf-8") as out,
        ):
            raw_meta, raw_header = _read_header(raw)
            final_meta, final_header = _read_header(final)
            samples = final_header.split("\t")[9:]
            if raw_header.split("\t")[9:] != samples:
                raise IPyradError("Native and final VCF sample orders differ.")
            contigs = {
                line.split("=", 2)[2].split(",", 1)[0].rstrip(">"): i
                for i, line in enumerate(raw_meta) if line.startswith("##contig=<ID=")
            }
            sample_index = {name: i for i, name in enumerate(samples)}
            entries = {}
            masks = list(hard_masks)
            shared_field = 9 + len(samples)
            if shared_mask is not None:
                sample_index[""] = len(samples)
                masks.append({"": shared_mask})
            for layer in masks:
                existing = {name: path for name, path in layer.items() if path.is_file()}
                for chrom, intervals in _load_mask_entries_from_paths(existing, sample_index).items():
                    for entry in intervals:
                        entry["intervals"].sort()
                    entries.setdefault(chrom, []).extend(intervals)

            tags = ("GT", *_EVIDENCE_FIELDS)
            header = [line for line in final_meta if not line.startswith("##FORMAT=")]
            for tag in tags:
                lines = [line for line in raw_meta if line.startswith(f"##FORMAT=<ID={tag},")]
                if len(lines) != 1:
                    raise IPyradError(f"Native VCF is missing a unique FORMAT/{tag} declaration.")
                header.extend(lines)
            out.write("\n".join([*header, _FT_HEADER, final_header]) + "\n")
            raw_groups = _coordinate_groups(raw, contigs)
            current = next(raw_groups, None)
            written = skipped = 0
            for key, final_group in _coordinate_groups(final, contigs):
                while current is not None and current[0] < key:
                    current = next(raw_groups, None)
                originals = current[1] if current is not None and current[0] == key else []
                for called in final_group:
                    alleles = [called[3], *called[4].split(",")]
                    matches = [record for record in originals if record[3:5] == called[3:5]]
                    if (
                        len(matches) != 1
                        or sum(record[3:5] == called[3:5] for record in final_group) != 1
                        or not all(allele in ("A", "C", "G", "T") for allele in alleles)
                        or len(set(alleles)) != len(alleles)
                    ):
                        skipped += 1
                        continue
                    native = matches[0]
                    excluded = set(_masked_fields(entries.get(called[0], ()), key[1] - 1))
                    raw_tags, called_tags = native[8].split(":"), called[8].split(":")
                    fields = []
                    for index, (raw_sample, called_sample) in enumerate(zip(native[9:], called[9:], strict=True), 9):
                        if index in excluded or shared_field in excluded:
                            evidence = ("./.", ".", ".", ".", ".", "HardMask")
                        else:
                            evidence = _evidence_sample(
                                dict(zip(raw_tags, raw_sample.split(":"))),
                                dict(zip(called_tags, called_sample.split(":"))),
                                len(alleles),
                            )
                        fields.append(":".join(evidence))
                    out.write("\t".join([*called[:7], ".", ":".join([*tags, "FT"]), *fields]) + "\n")
                    written += 1
        staged = Path(temp) / out_bcf.name
        run_pipeline([[BIN_BCF, "view", "--no-version", "-Ob", "-o", str(staged), str(plain)]])
        run_pipeline([[BIN_BCF, "index", "-c", str(staged)]])
        staged.with_suffix(".bcf.csi").replace(out_index)
        try:
            staged.replace(out_bcf)
        except BaseException:
            out_index.unlink()
            raise
    logger.info("wrote genotype evidence: {} matched SNPs, {} skipped records", written, skipped)
    return out_bcf
