#!/usr/bin/env python

"""Parse names from PE or SE filenames.

PE
-----
- custom delim
- split on all '.' and keep longest left substring w/ proper pairing.
- split on all '_' and keep longest left substring w/ proper pairing.
- cut from right and longest with proper pairing.
- else assume SE

SE
----
- custom delim
- cut to find longest unique left substrings
"""


from typing import List, Tuple, Dict, Union
from pathlib import Path
from collections import defaultdict
from loguru import logger
from .exceptions import IPyradError


def expand_path(p: str | Path) -> Path:
    """Returns an absolute path after expanding ~ and env variables"""
    return Path(p).expanduser().absolute()


def get_paths_list_from_fastq_str(fastq_paths: Union[Path, List[Path]]) -> List[Path]:
    """Expand fastq_paths str (e.g., 'data/*.gz') into List[Path].
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
        # expand to a full path
        path = expand_path(path)

        # raise if path is a dir.
        if path.is_dir():
            raise IPyradError(f"{path} is a dir. Use regex to select files in the dir (e.g., './path/*.fastq.gz')")

        # expand a regex operator to possibly match multiple files
        # such as paired-end files.
        try:
            fastqs = list(path.parent.glob(path.name))
            assert fastqs
        except (ValueError, AssertionError):
            raise IPyradError(f"No fastq data match input: {path}")
        expanded.extend([expand_path(i) for i in fastqs])
    return list(set(expanded))


def get_name_to_fastq_dict(
    fastqs: List[Path | str],
    delim: str | None,
    delim_index: int | None,
    suffix: str | None = None,
) -> Dict[str, Tuple[Path, Path | None]]:
    """
    """
    # expand str to List[Path]
    paths_list = get_paths_list_from_fastq_str(fastqs)

    # parse List[Path] to {name: (Path, Path)} or {name: (Path, None)}
    fastq_dict = get_pairs_or_single_by_trim(paths_list, delim, delim_index)

    # add optional suffix
    if suffix is not None:
        fastq_dict = {f"{i}{suffix}": j for (i, j) in fastq_dict.items()}

    # report to logger
    total = len(fastq_dict)
    fmax = min(10, total)
    logger.info(f"showing first {fmax}/{total} names parsed from file paths")
    max_len = max(len(i) for i in fastq_dict)
    for fidx, name in enumerate(sorted(fastq_dict)):
        paths = fastq_dict[name]
        key_padded = name + " " * (max_len - len(name))
        if paths[1]:
            if fidx < fmax:
                logger.info(f"{key_padded} <- {(paths[0].name, paths[1].name)}")
        else:
            if fidx < fmax:
                logger.info(f"{key_padded} <- {paths[0].name}")
    return fastq_dict


def perfect_pairs(ndict: Dict[str, List[Path]], paths: List[Path]) -> bool:
    """valid PE name"""
    # All names have 2 fastq files
    a = all(len(v) == 2 for v in ndict.values())
    # Number of fq files detected == number passed in
    b = sum(len(v) for v in ndict.values()) == len(paths)

    # Check that last letters of sample names are unique
    if len(ndict) == 1:
        # In the case of raw undemuxed fastq data there will only be one PE pair
        c = True
    else:
        try:
            # Splitting pairs on '.' or '_' it can create ndict keys == "" so this
            # also handles that case by raising on [-1] index of an empty string
            c = len(set([i[-1] for i in ndict])) > 1
        except Exception as e:
            c = False
    return a & b & c


def all_unique(ndict: Dict[str, List[Path]], paths: List[Path]) -> bool:
    """valid SE name"""
    if len(ndict) == 1:
        # Only one fastq file, so unique by definition
        return True
    else:
        # All names have one and only one fastq file
        a = all(len(v) == 1 for v in ndict.values())
        # The number of kv pairs in ndict equals the number of fq passed in
        b = len(ndict) == len(paths)
        # Checks that the last letter is not the same in all sample names
        c = len(set([i[-1] for i in ndict])) > 1
        return a & b & c


def get_pairs_or_single_by_trim(
    fastqs: List[Path | str],
    delim: str | None,
    delim_index: int | None,
) -> Dict[str, Tuple[Path, Path | None]]:
    """..."""
    # try to pair sample by delim args
    if delim:
        names_to_paths = defaultdict(list)
        for path in fastqs:
            parts = path.name.split(delim)
            left = delim.join(parts[:delim_index])
            names_to_paths[left].append(path)
        if perfect_pairs(names_to_paths, fastqs):
            logger.info(f"paired files by user args: -dx={delim} -di={delim_index}")
            return {i: tuple(sorted(j)) for i, j in names_to_paths.items()}
        logger.warning("pairing files by name-delim failed. Falling back to auto-detection.")

    # try to pair samples by splitting on '.' from right
    hits = []
    idx = 1
    while 1:
        names_to_paths = defaultdict(list)
        for path in fastqs:
            stem = ".".join(path.name.split(".")[:-idx])
            names_to_paths[stem].append(path)
        if perfect_pairs(names_to_paths, fastqs):
            hits.append((idx, names_to_paths))
        if '' in names_to_paths:
            break
        idx += 1
    # keep the longest one
    if hits:
        idx, names_to_paths = max(hits, key=lambda x: len(list(x[1])[0]))
        logger.info(f"paired files by auto-splitting: -dx=. -di={-idx}")
        return {i: tuple(sorted(j)) for i, j in names_to_paths.items()}

    # try to pair samples by splitting on '_' from right
    hits = []
    idx = 1
    while 1:
        names_to_paths = defaultdict(list)
        for path in fastqs:
            stem = "_".join(path.name.split("_")[:-idx])
            names_to_paths[stem].append(path)
        if perfect_pairs(names_to_paths, fastqs):
            hits.append((idx, names_to_paths))
        if '' in names_to_paths:
            break
        idx += 1
    # keep the longest one
    if hits:
        idx, names_to_paths = max(hits, key=lambda x: len(list(x[1])[0]))
        logger.info(f"paired files by auto-splitting: -dx=_ -di={-idx}")
        return {i: tuple(sorted(j)) for i, j in names_to_paths.items()}

    # try to pair samples by splitting on '_' from left
    hits = []
    idx = 1
    while 1:
        names_to_paths = defaultdict(list)
        for path in fastqs:
            stem = "_".join(path.name.split("_")[:idx])
            names_to_paths[stem].append(path)
        if perfect_pairs(names_to_paths, fastqs):
            hits.append((idx, names_to_paths))
        if stem == path.name:
            break
        idx += 1
    # keep the longest one
    if hits:
        idx, names_to_paths = max(hits, key=lambda x: len(list(x[1])[0]))
        logger.info(f"paired files by auto-splitting: -dx=_ -di={idx}")
        return {i: tuple(sorted(j)) for i, j in names_to_paths.items()}

    # try to pair samples by stripping characters from right
    names = [i.name for i in fastqs]
    min_len = min(len(i) for i in names)
    hits = []
    for cut in range(1, min_len):
        names_to_paths = defaultdict(list)
        for path in fastqs:
            stem = path.name[:-cut]
            names_to_paths[stem].append(path)
        if perfect_pairs(names_to_paths, fastqs):
            logger.info("paired files by stripping characters")
            return {i: tuple(sorted(j)) for i, j in names_to_paths.items()}

    # --------------------------------------------------------------
    logger.info("failed to pair files, assuming data is single-end")
    # --------------------------------------------------------------

    # try to get unique SE names using delim args
    if delim:
        names_to_paths = defaultdict(list)
        for path in fastqs:
            parts = path.name.split(delim)
            left = delim.join(parts[:delim_index])
            names_to_paths[left].append(path)
        if all_unique(names_to_paths, fastqs):
            logger.info(f"parsed names by user args: -dx={delim} -di={delim_index}")
            return {i: (j[0], None) for i, j in names_to_paths.items()}
        logger.info("parsing names by user args failed. Falling back to auto-detection.")

    # If only 1 SE file then assume it is un-demux data and not a sample
    if len(fastqs) == 1:
        stem = fastqs[0].name.rsplit(".", 2)[0]
        return {stem: (fastqs[0], None)}

    # get unique SE name by stripping characters from right
    hits = []
    for cut in range(1, min_len):
        names_to_paths = defaultdict(list)
        for path in fastqs:
            stem = path.name[:-cut]
            names_to_paths[stem].append(path)
        if all_unique(names_to_paths, fastqs):
            logger.info("parsed names by stripping characters")
            return {i: (j[0], None) for i, j in names_to_paths.items()}

    raise IPyradError(
        "Cannot parse names from file names, likely because they do not share "
        "the same suffix. Try setting the delim args explicitly, or entering "
        "paths with different naming conventions as separate inputs, e.g., "
        "'-d *.fq.gz *.fastq.gz'."
    )


if __name__ == "__main__":

    pass
    # from ipyrad3.utils.logger import set_log_level
    # set_log_level("DEBUG")

    # path = Path("/home/deren/Documents/tools/ipyrad2/examples/Pedic-PE-ddRAD/*_R*")
    # get_fastq_tuples_dict_from_paths_list(path)

    a = Path("/home/deren/Documents/ipyrad-tests/examples/Pedic-PE-ddRAD/c*_R*")
    b = Path("/home/deren/Documents/ipyrad-tests/examples/Pedic-PE-ddRAD/m*_R*")
    fastqs = [a, b]
    # fastqs = Path("/tmp/TRIM_SE/*.gz")
    # fastqs = Path("/tmp/MAP_PE/*.bam")
    pdict = get_name_to_fastq_dict(fastqs, None, None)
    # print(pdict)
