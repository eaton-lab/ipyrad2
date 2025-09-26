#!/usr/bin/env python

"""Group files into pairs based on file name matching patterns.

The `get_filenames_to_paired_fastqs` function is used in both demux
and step1 of assembly to return full Paths to fastq files as tuples
of either (R1, None) or (R1, R2).
"""

from __future__ import annotations
from typing import List, Tuple, Dict, Union, Iterable, Optional
from pathlib import Path
from collections import defaultdict
from loguru import logger

__all__ = [
    "get_paths_list_from_fastq_str",
    "get_name_to_fastq_dict",
]
PathPair = Tuple[Path, Path]
Buckets = Dict[str, List[Path]]


def get_paths_list_from_fastq_str(fastq_paths: Union[Path, List[Path]]) -> List[Path]:
    """Expand fastq_paths str argument into a list of Path objects.

    This is used within `get_fastq_tuples_dict_from_paths_list`.
    """
    expanded = []
    # ensure paths is a List[Path] but where the Path elements may be
    # regex path names that have not yet been expanded.

    # ensure it is a list
    if isinstance(fastq_paths, (str, Path)):
        fastq_paths = [fastq_paths]

    # ensure each is a Path object
    paths = []
    for path in fastq_paths:
        if isinstance(path, str):
            path = Path(path)
        paths.append(path)

    # for each Path in paths list expand into a list of Paths
    for path in paths:
        # expand a regex operator to possibly match multiple files
        # such as paired-end files.
        try:
            fastqs = list(path.parent.glob(path.name))
            assert fastqs
        except (ValueError, AssertionError):
            msg = f"No fastq data match input: {path}"
            logger.error(msg)
            raise ValueError(msg)
        expanded.extend([Path(i).expanduser().resolve() for i in fastqs])
    return expanded


def pair_or_single_by_name_right_trim(
    paths: Iterable[Path],
    delim_index: Optional[Tuple[str, int]] = None,  # e.g., ('_', 1), ('_', 2), ('_', -1)
    skip_paired: bool = False,
) -> Dict[str, PathPair]:
    """
    Pairing (first match wins):
      A) If delim_index is given, split at the Nth occurrence of `delim` (idx>0 from left,
         idx<0 from right), keep the left side; if all buckets have size==2, return pairs.
      B) Strip ALL file suffixes (e.g., '.fastq.gz' -> base) and check for perfect pairs.
      C) Right-trim character-by-character: for cut in [max_len..1], group by name[:cut],
         return the first cut that yields perfect pairs.

    Single-end fallback (if no pairing rule works):
      D) If delim_index produces unique keys (all buckets size==1), return {key: (p, -null-)}.
      E) Otherwise, find the minimal cut k (from 1..max_len) such that name[:k] is unique
         across all paths; then try to extend each key to the next '.' (excluding the dot)
         while keeping uniqueness. Return singles with whichever of these keeps uniqueness.

    Returns a dict with deterministic ordering; paired values are sorted by filename.
    """
    paths = list(paths)
    if not paths:
        return {}

    names = [p.name for p in paths]
    max_len = max(len(n) for n in names)

    # ---------- helpers ----------
    def key_left_of_nth_delim(name: str, delim: str, idx: int) -> str:
        """Split at the idx-th delimiter occurrence and keep the left part.
        idx>0: from left (1-based). idx<0: from right (-1 is last). If not found, return name.
        """
        pos: List[int] = []
        start, L = 0, len(delim)
        while True:
            j = name.find(delim, start)
            if j == -1:
                break
            pos.append(j)
            start = j + L
        if not pos or idx == 0:
            return name
        if idx > 0:
            return name[:pos[idx - 1]] if idx <= len(pos) else name
        k = len(pos) + idx  # negative index from right
        return name[:pos[k]] if 0 <= k < len(pos) else name

    def bucket_by_delim_idx(delim: str, idx: int) -> Buckets:
        b: Buckets = defaultdict(list)
        for p in paths:
            b[key_left_of_nth_delim(p.name, delim, idx)].append(p)
        return b

    def strip_all_suffixes(name: str) -> str:
        suf = Path(name).suffixes  # handles compound suffixes like .fastq.gz
        return name[: -len("".join(suf))] if suf else name

    def bucket_by_full_suffix_strip() -> Buckets:
        b: Buckets = defaultdict(list)
        for p in paths:
            b[strip_all_suffixes(p.name)].append(p)
        return b

    def bucket_by_cut(cut: int) -> Buckets:
        b: Buckets = defaultdict(list)
        for p in paths:
            b[p.name[:-cut]].append(p)
        return b

    def perfect_pairs(b: Buckets) -> bool:
        return all(len(v) == 2 for v in b.values()) and sum(len(v) for v in b.values()) == len(paths)

    def all_unique(b: Buckets) -> bool:
        return all(len(v) == 1 for v in b.values()) and len(b) == len(paths)

    def to_pairs(b: Buckets) -> Dict[str, PathPair]:
        out: Dict[str, PathPair] = {}
        for k, items in b.items():
            a, b_ = sorted(items, key=lambda x: x.name)
            out[k] = (a, b_)
        return dict(sorted(out.items()))

    def to_singles(b: Buckets) -> Dict[str, PathPair]:
        # Keep your sentinel Path("-null-") for R2
        out = {k: (v[0], Path("-null-")) for k, v in b.items()}
        return dict(sorted(out.items()))

    def extend_to_next_dot(name: str, k: int) -> str:
        """Extend name[:k] to the next '.' after k (dot excluded)."""
        j = name.find('.', k)
        return name[:j] if j != -1 else name[:k]
    # -----------------------------

    if not skip_paired:
        # ===== Pairing attempts =====
        # A) delimiter override
        if delim_index is not None:
            delim, idx = delim_index
            b = bucket_by_delim_idx(delim, idx)
            if perfect_pairs(b):
                return to_pairs(b)
            logger.debug("pairing files by name-delim failed, falling back to auto-detection")

        # B) strip all suffixes
        b_strip = bucket_by_full_suffix_strip()
        if perfect_pairs(b_strip):
            logger.debug("successfully paired samples")
            return to_pairs(b_strip)

        # C) right-trim (longest cut first)
        # for cut in range(max_len, 0, -1):
        for cut in range(0, max_len):
            b = bucket_by_cut(cut)
            if perfect_pairs(b):
                logger.debug("successfully paired samples")
                return to_pairs(b)
        logger.debug("pairing files by auto-detection failed, assuming data are R1 only")

    # ===== Single-end fallback =====
    # D) delimiter uniqueness
    if delim_index is not None:
        delim, idx = delim_index
        b = bucket_by_delim_idx(delim, idx)
        if all_unique(b):
            logger.debug("names parsed from non-paired files using name-delim")
            return to_singles(b)
        logger.debug("names could not be parsed from non-paired files using name-delim, falling back to auto-detect")

    # E) minimal unique prefix + optional dot-boundary prettify
    # Find minimal cut k with unique keys
    k_min_unique: Optional[int] = None
    b_min_unique: Optional[Buckets] = None
    for k in range(1, max_len + 1):
        b = bucket_by_cut(k)
        if all_unique(b):
            k_min_unique, b_min_unique = k, b
            break

    if b_min_unique is None:
        # Last resort: full names as keys
        return {p.name: (p, Path("-null-")) for p in sorted(paths, key=lambda x: x.name)}

    # Try to extend each key to next '.' while keeping uniqueness
    b_ext: Buckets = defaultdict(list)
    for (p,) in b_min_unique.values():
        b_ext[extend_to_next_dot(p.name, k_min_unique)].append(p)
    if all_unique(b_ext):
        return to_singles(b_ext)

    return to_singles(b_min_unique)


def get_name_to_fastq_dict(
    fastqs: List[Path | str],
    delim_index: Optional[Tuple[int, str]] = None,
    skip_paired: bool = False,
) -> Dict[str, PathPair]:
    """Return {name: (Path, Path)} from on or more str path args."""

    # parse paths and names-to-path-pairs
    paths = get_paths_list_from_fastq_str(fastqs)
    fastq_dict = pair_or_single_by_name_right_trim(paths, delim_index, skip_paired)

    # log to user
    logger.info("sample names parsed from file paths")
    max_len = max(len(i) for i in fastq_dict)
    for key, val in fastq_dict.items():
        key_padded = key + " " * (max_len - len(key))
        if val[1].name == "-null-":
            logger.info(f"{key_padded} <- {val[0].name}")
        else:
            logger.info(f"{key_padded} <- {(val[0].name, val[1].name)}")
    return fastq_dict


if __name__ == "__main__":

    pass
    # from ipyrad3.utils.logger import set_log_level
    # set_log_level("DEBUG")

    # path = Path("/home/deren/Documents/tools/ipyrad2/examples/Pedic-PE-ddRAD/*_R*")
    # get_fastq_tuples_dict_from_paths_list(path)

    fastqs = Path("/tmp/TRIM/*.R*.gz")
    fastqs = Path("/tmp/TRIM_SE/*.gz")
    fastqs = Path("/tmp/MAP_PE/*.bam")
    pdict = get_name_to_fastq_dict(fastqs)#, ("plate", 1))

    # get_fastq_tuples_dict_from_paths_list(path)


    # paths = get_paths_list_from_fastq_str(path)
    # ndict = get_paired_clean_names_dict(paths)
    # if not ndict:
    #     ndict = get_single_clean_names_dict(paths)
    # for k, v in ndict.items():
    #     print(k, v)

    # path = Path("/data/name.x._R1.fq")
    # logger.debug(path)
    # x = get_delimited_name_chunks(path)
    # logger.debug(x)

    # path = Path("/d/name.x._R1.fq.gz")
    # x = get_delimited_name_chunks(path)    
    # logger.warning(x)

    # name.x._R1.fq.gz  # full
    # name_R1.fq.gz  # 1
    # name.xR1.fq.gz # 2
    # name.x.fq.gz   # 3
    # name.x._R1gz   # 4
    # name.x._R1.fq  # 5                

    # all of these paired samples work
    # path1 = "/data/name.x.y.z_R1.suf.fix"
    # path2 = "/data/name.x.y.z_R2.suf.fix"
    # path1 = "/data/name.x.y.z_1.suf.fix"
    # path2 = "/data/name.x.y.z_2.suf.fix"
    # path1 = "/data/name.x.y.z_1_001.suf.fix"
    # path2 = "/data/name.x.y.z_2_001.suf.fix"
    # path1 = "/data/name_a_b_R1.suf.fix"
    # path2 = "/data/name_a_b_R2.suf.fix"
    # p1 = Path(path1)
    # p2 = Path(path2)
    # fastqs = [p1, p2]
    # res = get_paired_samples(fastqs)

    # se samples
    # path1 = "/data/name1_R1.suf.fix"
    # path2 = "/data/name2_R1.suf.fix"
    # path3 = "/data/name3.12_R1.suf.fix"
    # fastqs = [Path(i) for i in (path1, path2, path3)]
    # res = get_single_samples(fastqs)
    # logger.info(res)
    # # logger.info(get_delimited_name_chunks(p1))
    # # logger.info(get_delimited_name_chunks(p2))

    # assert groups

    # get_new_pairs([p1, p2])
    # print(p1)
    # path = p1
    # suffixes = path.suffixes
    # while path.suffix in suffixes:
    #     path = path.with_suffix('')
    # print(path)
    # print(drop_from_right(p1, ("_", "."), 0))  # drops _I.
    # print(drop_from_right(p1, ("_", "."), 1))  # drops sample_
    # print(drop_from_right(p1, ("_", "."), 2))  # drops nothing

    # FASTQ_PATH = Path("/tmp/tmp-test____/sample_I*")
    # FASTQS = get_paths_list_from_fastq_str(FASTQ_PATH)
    # pairs = get_fastq_tuples_dict_from_paths_list(FASTQS)

    # FASTQ_PATH = Path("../../pedtest/small_tmp_R*.gz")
    # FASTQS = get_paths_list_from_fastq_str(FASTQ_PATH)
    # pairs = get_fastq_tuples_dict_from_paths_list(FASTQS)

    # FASTQ_PATH = Path("../../pedtest/NEW_fastqs/*fastq.gz")
    # FASTQS = get_paths_list_from_fastq_str(FASTQ_PATH)
    # pairs = get_fastq_tuples_dict_from_paths_list(FASTQS)

    # FASTQ_PATH = Path("../../sra-fastqs/*fastq")
    # FASTQS = get_paths_list_from_fastq_str(FASTQ_PATH)
    # pairs = get_fastq_tuples_dict_from_paths_list(FASTQS)
