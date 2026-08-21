#!/usr/bin/env python

"""Single-run BPP analysis helpers."""

from __future__ import annotations

import os
import platform
import re
import shutil
import subprocess as sps
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path
from types import SimpleNamespace

import requests
from loguru import logger

from ..extracters.seqex import SeqexEngine
from ..extracters.sequence_common import normalize_sequence_population_inputs
from ..extracters.sequence_common import resolve_sequence_sample_subset
from ...utils.exceptions import IPyradError


_MISSING_TOYTREE = """
You are missing required packages to use ipa.bpp().
First run the following conda install command:

conda install toytree -c conda-forge
"""

try:
    import toytree
except ImportError as exc:
    raise IPyradError(_MISSING_TOYTREE) from exc


DELIM = "___"
BPP_DOCS_VERSION = "4.8.7"
_MINIMUM_BPP_VERSION = (4, 8, 7)
_ANALYSES = {"a00", "a01", "a10", "a11", "msc-i", "msc-m"}
_SUBSTITUTION_MODELS = {"JC69", "K80", "F81", "HKY", "T92", "TN93", "F84", "GTR"}
_THETA_MODELS = {
    "linked-none",
    "linked-all",
    "linked-inner",
    "linked-msci",
    "linked-mscm",
}
_DEFAULT_LOCUS_RATE = ("0",)
_DEFAULT_CLOCK = ("1",)
_BPP_BINARY_SPECS = {
    ("linux", "x86_64"): (
        "4.8.7",
        "bpp-4.8.7-linux-x86_64.tar.gz",
        "bpp-4.8.7-linux-x86_64/bin/bpp",
    ),
    ("linux", "aarch64"): (
        "4.8.7",
        "bpp-4.8.7-linux-aarch64.tar.gz",
        "bpp-4.8.7-linux-aarch64/bin/bpp",
    ),
    ("darwin", "arm64"): (
        "4.8.7",
        "bpp-4.8.7-macos-aarch64.tar.gz",
        "bpp-4.8.7-macos-aarch64/bin/bpp",
    ),
    ("darwin", "x86_64"): (
        "4.8.7",
        "bpp-4.8.7-macos-x86_64.tar.gz",
        "bpp-4.8.7-macos-x86_64/bin/bpp",
    ),
    ("win32", "x86_64"): (
        "4.8.7",
        "bpp-4.8.7-win-x86_64.zip",
        "bpp-4.8.7-win-x86_64/bpp.exe",
    ),
}


def _normalize_machine(machine: str) -> str:
    """Normalize architecture strings to the values used by bundled BPP binaries."""
    machine = machine.lower()
    aliases = {
        "amd64": "x86_64",
        "x64": "x86_64",
        "arm64": "arm64",
        "aarch64": "aarch64",
    }
    return aliases.get(machine, machine)


def _bpp_target_key() -> tuple[str, str]:
    """Return the `(platform, machine)` key used for bundled BPP binaries."""
    plat = sys.platform
    if plat.startswith("linux"):
        plat = "linux"
    elif plat == "darwin":
        plat = "darwin"
    elif plat.startswith("win"):
        plat = "win32"
    return plat, _normalize_machine(platform.machine())


def _get_bpp_download_spec() -> SimpleNamespace:
    """Return the bundled BPP binary download metadata for the current target."""
    key = _bpp_target_key()
    if key not in _BPP_BINARY_SPECS:
        raise IPyradError(
            "No bundled BPP binary is available for platform={} arch={}.".format(*key)
        )
    version, archive_name, binary_relpath = _BPP_BINARY_SPECS[key]
    return SimpleNamespace(
        version=version,
        archive_name=archive_name,
        binary_relpath=binary_relpath,
        url=f"https://github.com/bpp/bpp/releases/download/v{version}/{archive_name}",
        archive_path=Path(tempfile.gettempdir()) / archive_name,
        extract_dir=Path(tempfile.gettempdir()),
        binary_path=Path(tempfile.gettempdir()) / binary_relpath,
    )


def _coerce_positive_int(value, label: str) -> int:
    """Parse one positive integer."""
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise IPyradError(f"{label} must be an integer.") from exc
    if parsed < 1:
        raise IPyradError(f"{label} must be >= 1.")
    return parsed


def _coerce_positive_number(value, label: str) -> float:
    """Parse one positive float."""
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise IPyradError(f"{label} must be numeric.") from exc
    if parsed <= 0:
        raise IPyradError(f"{label} must be > 0.")
    return parsed


def _normalize_prior_pair(name: str, values) -> tuple[float, float]:
    """Normalize a two-number prior."""
    if values is None or len(values) != 2:
        raise IPyradError(f"{name} must contain exactly two numeric values.")
    return (
        _coerce_positive_number(values[0], f"{name}[0]"),
        _coerce_positive_number(values[1], f"{name}[1]"),
    )


def _normalize_analysis(value) -> str:
    """Return one canonical, case-insensitive BPP analysis name."""
    analysis = str(value).strip().casefold()
    if analysis not in _ANALYSES:
        allowed = ", ".join(sorted(_ANALYSES))
        raise IPyradError(f"analysis must be one of: {allowed}.")
    return analysis


def _normalize_theta_prior(values) -> tuple[str, ...]:
    """Normalize a BPP theta prior."""
    if values is None:
        values = ("invgamma", 3, 0.002, "estimate")
    tokens = tuple(str(value) for value in values)
    if not tokens:
        raise IPyradError("theta_prior must not be empty.")
    distribution = tokens[0].casefold()
    if distribution in {"invgamma", "gamma"}:
        valid_lengths = {3, 4} if distribution == "invgamma" else {3}
        if len(tokens) not in valid_lengths:
            suffix = " [estimate|integrate]" if distribution == "invgamma" else ""
            raise IPyradError(f"theta_prior {distribution} requires A B{suffix}.")
        first = _coerce_positive_number(tokens[1], "theta_prior[1]")
        second = _coerce_positive_number(tokens[2], "theta_prior[2]")
        result = [distribution, f"{first:g}", f"{second:g}"]
        if len(tokens) == 4:
            estimate = tokens[3].casefold()
            if estimate not in {"estimate", "integrate"}:
                raise IPyradError(
                    "theta_prior inverse-gamma mode must be estimate or integrate."
                )
            if estimate == "integrate":
                result.append("int")
        return tuple(result)
    if distribution == "beta":
        if len(tokens) != 5:
            raise IPyradError("theta_prior beta requires P Q LOWER UPPER.")
        numbers = tuple(
            _coerce_positive_number(value, f"theta_prior[{idx}]")
            for idx, value in enumerate(tokens[1:], start=1)
        )
        if numbers[3] <= numbers[2]:
            raise IPyradError("theta_prior beta UPPER must be greater than LOWER.")
        return ("beta", *(f"{value:g}" for value in numbers))
    raise IPyradError("theta_prior distribution must be invgamma, gamma, or beta.")


def _normalize_tau_prior(values) -> tuple[str, ...]:
    """Normalize a BPP divergence-time prior."""
    if values is None:
        values = ("invgamma", 3, 0.03)
    tokens = tuple(str(value) for value in values)
    if len(tokens) != 3:
        raise IPyradError("tau_prior requires DIST A B.")
    distribution = tokens[0].casefold()
    if distribution not in {"invgamma", "gamma"}:
        raise IPyradError("tau_prior distribution must be invgamma or gamma.")
    first = _coerce_positive_number(tokens[1], "tau_prior[1]")
    second = _coerce_positive_number(tokens[2], "tau_prior[2]")
    return distribution, f"{first:g}", f"{second:g}"


def _normalize_alpha_prior(values) -> tuple[float, float, int] | None:
    """Normalize `alphaprior = alpha beta ncat`."""
    if values is None:
        return None
    if len(values) != 3:
        raise IPyradError("alpha_prior must contain exactly three values.")
    return (
        _coerce_positive_number(values[0], "alpha_prior[0]"),
        _coerce_positive_number(values[1], "alpha_prior[1]"),
        _coerce_positive_int(values[2], "alpha_prior[2]"),
    )


def _normalize_threads(value) -> tuple[int, ...] | None:
    """Normalize BPP thread specifications."""
    if value is None:
        return None
    if isinstance(value, int):
        return (_coerce_positive_int(value, "threads"),)
    parsed = tuple(_coerce_positive_int(i, "threads") for i in value)
    if len(parsed) not in (1, 3):
        raise IPyradError("threads must contain exactly 1 or 3 positive integers.")
    return parsed


def _normalize_seed(seed, label: str) -> int | None:
    """Normalize CLI/API seed inputs."""
    if seed is None:
        return None
    if isinstance(seed, str) and seed.casefold() == "none":
        return None
    return _coerce_positive_int(seed, label)


def _normalize_locus_rate(values) -> tuple[str, ...]:
    """Normalize supported BPP locus-rate modes."""
    if values is None:
        return _DEFAULT_LOCUS_RATE
    if not values:
        raise IPyradError("locus_rate must not be empty.")
    tokens = tuple(str(value) for value in values)
    if tokens == ("0",):
        return tokens
    if tokens[0] in {"2", "3"}:
        feature = "fixed rate files" if tokens[0] == "2" else "tip dating"
        raise IPyradError(f"locus_rate mode {tokens[0]} ({feature}) is not supported.")
    if len(tokens) != 5 or tokens[0] != "1":
        raise IPyradError("locus_rate must be `0` or `1 A B C iid|dir`.")
    numbers = [
        _coerce_positive_number(token, f"locus_rate[{idx}]")
        for idx, token in enumerate(tokens[1:4], start=1)
    ]
    dependence = tokens[4].casefold()
    if dependence not in {"iid", "dir"}:
        raise IPyradError("locus_rate mode 1 must end with iid or dir.")
    return ("1", *(f"{value:g}" for value in numbers), dependence)


def _normalize_clock(values) -> tuple[str, ...]:
    """Normalize supported BPP clock modes."""
    if values is None:
        return _DEFAULT_CLOCK
    if not values:
        raise IPyradError("clock must not be empty.")
    tokens = tuple(str(value) for value in values)
    if tokens == ("1",):
        return tokens
    if tokens[0] == "4":
        if len(tokens) != 2:
            raise IPyradError("clock mode 4 requires one positive parameter.")
        value = _coerce_positive_number(tokens[1], "clock[1]")
        return "4", f"{value:g}"
    if tokens[0] not in {"2", "3"} or len(tokens) != 6:
        raise IPyradError("clock must be `1`, `2|3 A B C iid|dir G|LN`, or `4 A`.")
    numbers = [
        _coerce_positive_number(token, f"clock[{idx}]")
        for idx, token in enumerate(tokens[1:4], start=1)
    ]
    dependence = tokens[4].casefold()
    distribution = tokens[5].upper()
    if dependence not in {"iid", "dir"}:
        raise IPyradError("clock modes 2 and 3 require iid or dir.")
    if distribution not in {"G", "LN"}:
        raise IPyradError("clock modes 2 and 3 require G or LN.")
    return (
        tokens[0],
        *(f"{value:g}" for value in numbers),
        dependence,
        distribution,
    )


def _normalize_delimitation_proposal(values) -> tuple[str, ...]:
    """Normalize an rjMCMC species-delimitation proposal."""
    if values is None:
        values = (0, 2)
    tokens = tuple(str(value) for value in values)
    if len(tokens) == 2 and tokens[0] == "0":
        epsilon = _coerce_positive_number(tokens[1], "delimitation_proposal[1]")
        return "0", f"{epsilon:g}"
    if len(tokens) == 3 and tokens[0] == "1":
        alpha = _coerce_positive_number(tokens[1], "delimitation_proposal[1]")
        mvalue = _coerce_positive_number(tokens[2], "delimitation_proposal[2]")
        return "1", f"{alpha:g}", f"{mvalue:g}"
    raise IPyradError(
        "delimitation_proposal must be `0 EPSILON` or `1 ALPHA M`."
    )


def _normalize_checkpoint(values) -> tuple[int, ...] | None:
    """Normalize checkpoint start and optional interval."""
    if values is None:
        return None
    if isinstance(values, int):
        values = (values,)
    parsed = tuple(_coerce_positive_int(value, "checkpoint") for value in values)
    if len(parsed) not in {1, 2}:
        raise IPyradError("checkpoint requires START and optional INTERVAL.")
    return parsed

def _split_top_level_commas(text: str) -> list[str]:
    """Split one migration token on commas, ignoring commas inside parentheses."""
    parts = []
    chunk = []
    depth = 0
    for char in text:
        if char == "(":
            depth += 1
        elif char == ")":
            depth = max(0, depth - 1)
        if char == "," and depth == 0:
            part = "".join(chunk).strip()
            if not part:
                raise IPyradError(f"Malformed migration token: {text}")
            parts.append(part)
            chunk = []
            continue
        chunk.append(char)
    tail = "".join(chunk).strip()
    if tail:
        parts.append(tail)
    if len(parts) != 2:
        raise IPyradError(
            "Each migration token must define one source,target pair. "
            f"Received: {text}"
        )
    return parts


def _normalize_migration_tokens(values) -> tuple[tuple[str, str], ...]:
    """Normalize `--msc-m` values into `(source, target)` pairs."""
    if not values:
        raise IPyradError("--msc-m requires one or more migration pairs.")
    pairs = []
    for value in values:
        source, target = _split_top_level_commas(str(value))
        pairs.append((source, target))
    return tuple(pairs)


def _read_tree_text(tree) -> str:
    """Return Newick text from an inline string or a file path."""
    if tree is None:
        raise IPyradError("A guide tree is required.")
    if isinstance(tree, Path):
        path = tree.expanduser().absolute()
        if not path.exists():
            raise IPyradError(f"guide tree file does not exist: {path}")
        text = path.read_text(encoding="utf-8")
    else:
        tree_text = str(tree).strip()
        path = Path(tree_text).expanduser()
        if path.exists() and path.is_file():
            text = path.absolute().read_text(encoding="utf-8")
        else:
            text = tree_text
    text = text.strip()
    if not text:
        raise IPyradError("guide tree text is empty.")
    if not text.endswith(";"):
        raise IPyradError("guide tree must end with ';'.")
    return text


def _parse_tree(tree_text: str):
    """Parse a guide tree with toytree."""
    try:
        return toytree.tree(tree_text)
    except Exception as exc:  # pragma: no cover - error type varies with parser
        raise IPyradError(f"failed to parse guide tree: {exc}") from exc


def _prefix_tip_labels(tree_text: str, tip_labels, prefix: str = DELIM) -> str:
    """Prefix only tip labels inside a Newick or extended Newick string."""
    prefixed = tree_text
    for label in sorted(tip_labels, key=len, reverse=True):
        pattern = rf"([,(])(\s*){re.escape(label)}(?=(?::|,|\)|\[))"
        prefixed = re.sub(pattern, rf"\1\2{prefix}{label}", prefixed)
    return prefixed


def _resolve_bpp_binary_version(binary: str) -> str | None:
    """Return a best-effort version string for the resolved BPP binary."""
    try:
        proc = sps.run(
            [binary, "--help"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, sps.SubprocessError):
        return None
    text = f"{proc.stdout}\n{proc.stderr}"
    match = re.search(r"\bv?(\d+\.\d+\.\d+)\b", text)
    return match.group(1) if match else None


def _resolve_bpp_binary(binary: str | None) -> str:
    """Resolve the BPP binary from an explicit path, PATH, or the bundled fallback."""
    if binary:
        resolved = Path(binary).expanduser().absolute()
        if not resolved.is_file():
            raise IPyradError(f"BPP binary does not exist: {resolved}")
        if not os.access(resolved, os.X_OK):
            raise IPyradError(f"BPP binary is not executable: {resolved}")
        return str(resolved)

    found = shutil.which("bpp")
    if found:
        return str(Path(found).resolve())

    spec = _get_bpp_download_spec()
    try:
        response = requests.get(spec.url, allow_redirects=True, stream=True, timeout=60)
        response.raise_for_status()
        with open(spec.archive_path, "wb") as archive:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    archive.write(chunk)
    except requests.RequestException as exc:
        raise IPyradError(
            f"Failed to download bundled BPP binary from {spec.url}: {exc}"
        ) from exc

    try:
        if str(spec.archive_path).endswith(".tar.gz"):
            with tarfile.open(spec.archive_path, "r:gz") as archive:
                archive.extractall(spec.extract_dir)
        elif str(spec.archive_path).endswith(".zip"):
            with zipfile.ZipFile(spec.archive_path, "r") as archive:
                archive.extractall(spec.extract_dir)
        else:
            raise IPyradError(f"Unsupported BPP archive format: {spec.archive_name}")
    except (tarfile.TarError, zipfile.BadZipFile, OSError) as exc:
        raise IPyradError(f"Failed to extract bundled BPP archive: {spec.archive_name}") from exc

    if not spec.binary_path.exists():
        raise IPyradError(
            f"Bundled BPP binary was not found after extraction: {spec.binary_path}"
        )
    return str(spec.binary_path.resolve())


def _extract_migration_labels(reference: str) -> tuple[str, ...]:
    """Extract plain labels from one migration reference expression."""
    labels = tuple(part for part in re.split(r"[(),\s]+", reference) if part)
    if not labels:
        raise IPyradError(f"Malformed migration reference: {reference}")
    return labels


class Bpp:
    """Prepare one BPP run, write its inputs, and optionally execute it."""

    def __init__(
        self,
        *,
        data,
        tree,
        imap,
        minmap,
        max_loci,
        min_length,
        name: str = "bpp",
        outdir: Path | str = "output-bpp",
        analysis: str = "a00",
        migration=None,
        theta_prior=None,
        tau_prior=None,
        species_model_prior: int | None = None,
        delimitation_proposal=None,
        phi_prior=(1.0, 1.0),
        w_prior=(2.0, 200.0),
        alpha_prior=None,
        substitution_model: str = "JC69",
        locus_rate=None,
        clock=None,
        theta_model: str = "linked-none",
        clean_data: bool = False,
        checkpoint=None,
        constraint_file=None,
        print_locus_rates: bool = False,
        print_gene_trees: bool = False,
        print_substitution_parameters: bool = False,
        burnin: int = 10000,
        samplefreq: int = 10,
        nsample: int = 100000,
        threads=None,
        locus_seed=None,
        mcmc_seed=None,
        binary: str | None = None,
        force: bool = False,
        log_level: str = "INFO",
    ) -> None:
        self.data = Path(data).expanduser().absolute()
        self.name = str(name)
        self.outdir = Path(outdir).expanduser().absolute()
        self.binary = None if binary is None else str(Path(binary).expanduser().absolute())
        self.force = bool(force)
        self.log_level = str(log_level)

        self.analysis = _normalize_analysis(analysis)
        self.max_loci = _coerce_positive_int(max_loci, "max_loci")
        self.min_length = _coerce_positive_int(min_length, "min_length")
        self.burnin = _coerce_positive_int(burnin, "burnin")
        self.samplefreq = _coerce_positive_int(samplefreq, "samplefreq")
        self.nsample = _coerce_positive_int(nsample, "nsample")
        self.theta_prior = _normalize_theta_prior(theta_prior)
        self.tau_prior = _normalize_tau_prior(tau_prior)
        self.phi_prior = _normalize_prior_pair("phi_prior", phi_prior)
        self.w_prior = _normalize_prior_pair("w_prior", w_prior)
        self.alpha_prior = _normalize_alpha_prior(alpha_prior)
        self.locus_rate = _normalize_locus_rate(locus_rate)
        self.clock = _normalize_clock(clock)
        self.delimitation_proposal = _normalize_delimitation_proposal(
            delimitation_proposal
        )
        self.threads = _normalize_threads(threads)
        self.locus_seed = _normalize_seed(locus_seed, "locus_seed")
        self.mcmc_seed = _normalize_seed(mcmc_seed, "mcmc_seed")
        self.clean_data = bool(clean_data)
        self.checkpoint = _normalize_checkpoint(checkpoint)
        self.print_locus_rates = bool(print_locus_rates)
        self.print_gene_trees = bool(print_gene_trees)
        self.print_substitution_parameters = bool(print_substitution_parameters)

        self.substitution_model = str(substitution_model).upper()
        self.theta_model = str(theta_model).casefold()
        self.species_model_prior = (
            None if species_model_prior is None else int(species_model_prior)
        )
        self.migration_pairs = (
            _normalize_migration_tokens(migration) if migration else ()
        )
        self.constraint_file = (
            None
            if constraint_file is None
            else Path(constraint_file).expanduser().absolute()
        )

        self.imap, self.minmap = normalize_sequence_population_inputs(imap, minmap)
        self.tree_text = _read_tree_text(tree)
        self.tree = _parse_tree(self.tree_text)
        self.tip_labels = tuple(self.tree.get_tip_labels())
        self.species_order = tuple(sorted(self.imap)) if self.imap else ()

        self.paths = SimpleNamespace(
            seqfile=self.outdir / f"{self.name}.phy",
            mapfile=self.outdir / f"{self.name}.imapfile.txt",
            ctlfile=self.outdir / f"{self.name}.ctl.txt",
            constraintfile=self.outdir / f"{self.name}.constraints.txt",
            statsfile=self.outdir / f"{self.name}.stats.txt",
            statsjson=self.outdir / f"{self.name}.stats.json",
            jobname=self.outdir / self.name,
            outfile=self.outdir / f"{self.name}.txt",
            mcmcfile=self.outdir / f"{self.name}.mcmc.txt",
            figtree=self.outdir / f"{self.name}.figtree.nex",
        )
        self.seqex = None
        self.nloci = self.max_loci
        self.effective_threads = self.threads

        self._validate()
        _, _, _, self.imap = resolve_sequence_sample_subset(
            self.data,
            exclude=None,
            include_reference=False,
            imap=self.imap,
        )

    @property
    def speciestree(self) -> bool:
        """Return whether BPP estimates species-tree topology."""
        return self.analysis in {"a01", "a11"}

    @property
    def speciesdelimitation(self) -> bool:
        """Return whether BPP estimates species delimitation."""
        return self.analysis in {"a10", "a11"}

    @property
    def geneflow(self) -> bool:
        """Return whether the selected model estimates gene flow."""
        return self.analysis in {"msc-i", "msc-m"}

    @property
    def algorithm(self) -> str:
        """Return the canonical analysis name used for logging and outputs."""
        return self.analysis

    def _validate(self) -> None:
        """Validate user inputs before writing files or launching BPP."""
        if not self.data.exists():
            raise IPyradError(f"data file does not exist: {self.data}")
        if self.data.suffix != ".hdf5":
            raise IPyradError("'data' must be an ipyrad2 .hdf5 file.")
        if not self.imap or not isinstance(self.imap, dict):
            raise IPyradError("imap is required and must resolve to a dictionary.")
        if self.minmap is not None and set(self.minmap) != set(self.imap):
            raise IPyradError("imap and minmap keys must match.")
        if set(self.species_order) != set(self.tip_labels):
            raise IPyradError(
                "IMAP keys must match guide tree tip names exactly.\n"
                f"imap={sorted(self.imap)}\n"
                f"tree={sorted(self.tip_labels)}"
            )
        if self.burnin <= 200:
            raise IPyradError("burnin must be greater than 200 with finetune enabled.")
        if self.threads and self.threads[0] > self.max_loci:
            raise IPyradError("threads cannot exceed max_loci.")

        model_sampling = self.analysis in {"a01", "a10", "a11"}
        if model_sampling:
            if self.species_model_prior is None:
                self.species_model_prior = 1
            allowed = {0, 1, 2, 3} if self.analysis == "a11" else {0, 1}
            if self.species_model_prior not in allowed:
                values = ", ".join(str(value) for value in sorted(allowed))
                raise IPyradError(
                    f"analysis {self.analysis.upper()} allows species_model_prior {values}."
                )
        elif self.species_model_prior is not None:
            raise IPyradError(
                "species_model_prior is only valid for A01, A10, and A11."
            )

        if self.analysis == "msc-i" and "&phi=" not in self.tree_text:
            raise IPyradError(
                "MSC-I requires an extended Newick tree containing &phi annotations."
            )
        if self.analysis == "msc-m" and not self.migration_pairs:
            raise IPyradError("MSC-M requires one or more migration pairs.")
        if self.analysis != "msc-m" and self.migration_pairs:
            raise IPyradError("migration pairs are only valid for MSC-M.")
        if self.migration_pairs:
            all_labels = set(self.tip_labels)
            for node in self.tree.treenode.traverse():
                if not node.is_leaf() and node.name:
                    all_labels.add(str(node.name))
            missing = sorted(
                {
                    label
                    for pair in self.migration_pairs
                    for reference in pair
                    for label in _extract_migration_labels(reference)
                    if label not in all_labels
                }
            )
            if missing:
                raise IPyradError(
                    "MSC-M migration labels are not present in the guide tree: "
                    + ", ".join(missing)
                )

        if self.substitution_model not in _SUBSTITUTION_MODELS:
            allowed = ", ".join(sorted(_SUBSTITUTION_MODELS))
            raise IPyradError(f"substitution_model must be one of: {allowed}.")
        if self.theta_model not in _THETA_MODELS:
            allowed = ", ".join(sorted(_THETA_MODELS))
            raise IPyradError(f"theta_model must be one of: {allowed}.")
        if self.theta_model == "linked-msci" and self.analysis != "msc-i":
            raise IPyradError("linked-msci theta_model is only valid for MSC-I.")
        if self.theta_model == "linked-mscm" and self.analysis != "msc-m":
            raise IPyradError("linked-mscm theta_model is only valid for MSC-M.")
        if self.analysis == "msc-i" and self.clock[0] == "3":
            raise IPyradError("clock mode 3 is not supported for MSC-I.")
        if self.locus_rate[0] == "1" and self.clock[0] in {"2", "3"}:
            if self.locus_rate[4] != self.clock[4]:
                raise IPyradError(
                    "locus_rate and clock must use the same iid or dir setting."
                )

        if not self.speciesdelimitation and self.delimitation_proposal != ("0", "2"):
            raise IPyradError(
                "delimitation_proposal is only valid for A10 and A11."
            )
        if self.constraint_file is not None:
            if self.analysis not in {"a01", "a11"}:
                raise IPyradError("constraint_file is only valid for A01 and A11.")
            if not self.constraint_file.is_file():
                raise IPyradError(
                    f"constraint file does not exist: {self.constraint_file}"
                )
            self._render_constraint_text()

    def _ensure_output_paths(self, *, include_results: bool) -> None:
        """Fail early when output files already exist and overwrite is disabled."""
        paths = [
            self.paths.seqfile,
            self.paths.mapfile,
            self.paths.ctlfile,
            self.paths.statsfile,
            self.paths.statsjson,
        ]
        if self.constraint_file is not None:
            paths.append(self.paths.constraintfile)
        if include_results:
            paths.extend([self.paths.outfile, self.paths.mcmcfile, self.paths.figtree])
        existing = next((path for path in paths if path.exists()), None)
        if existing is not None and not self.force:
            raise IPyradError(
                f"Output file already exists: {existing}. Use --force to overwrite."
            )

    def _prefixed_species_names(self) -> tuple[str, ...]:
        """Return sorted species labels with the BPP delimiter prefix."""
        return tuple(f"{DELIM}{name}" for name in self.species_order)

    def _render_tree_text(self) -> str:
        """Return guide-tree text with tip labels prefixed for BPP."""
        return _prefix_tip_labels(self.tree_text, self.tip_labels, prefix=DELIM)

    def _label_for_tree_reference(self, label: str) -> str:
        """Prefix tip labels inside one migration reference expression."""
        rendered = str(label)
        for tip in sorted(self.tip_labels, key=len, reverse=True):
            if rendered == tip:
                rendered = f"{DELIM}{tip}"
                continue
            pattern = rf"([,(])(\s*){re.escape(tip)}(?=(?:,|\)|$))"
            rendered = re.sub(pattern, rf"\1\2{DELIM}{tip}", rendered)
        return rendered

    def _render_constraint_text(self) -> str:
        """Validate and prefix species tokens in a BPP constraint file."""
        text = self.constraint_file.read_text(encoding="utf-8")
        aliases = set()
        statements = []
        for lineno, line in enumerate(text.splitlines(), start=1):
            statement = line.split("#", maxsplit=1)[0].strip()
            if not statement:
                continue
            keyword = statement.split(maxsplit=1)[0].casefold()
            if keyword not in {"define", "constraint", "outgroup"}:
                raise IPyradError(
                    f"constraint file line {lineno} must start with define, constraint, or outgroup."
                )
            if keyword == "define":
                match = re.match(r"(?i)^define\s+([^\s=]+)\s+(?:as|=)\s+", statement)
                if not match:
                    raise IPyradError(
                        f"constraint file line {lineno} has an invalid define statement."
                    )
                alias = match.group(1)
                if alias in self.species_order:
                    raise IPyradError(
                        f"constraint alias conflicts with population name: {alias}"
                    )
                aliases.add(alias)
            statements.append((lineno, statement))

        ignored = {"define", "constraint", "outgroup", "as", "not"}
        for lineno, statement in statements:
            tokens = re.findall(r"[A-Za-z0-9_.-]+", statement)
            unknown = sorted(
                token
                for token in tokens
                if token.casefold() not in ignored
                and token not in aliases
                and token not in self.species_order
            )
            if unknown:
                raise IPyradError(
                    f"constraint file line {lineno} contains unknown names: "
                    + ", ".join(unknown)
                )

        rendered = text
        for species in sorted(self.species_order, key=len, reverse=True):
            pattern = rf"(?<![A-Za-z0-9_.-]){re.escape(species)}(?![A-Za-z0-9_.-])"
            rendered = re.sub(pattern, f"{DELIM}{species}", rendered)
        return rendered

    def _write_seqfile(self) -> None:
        """Extract loci once and write the BPP-formatted sequence file."""
        self.seqex = SeqexEngine(
            data=self.data,
            name=self.name,
            outdir=self.outdir,
            out_format="phy",
            max_loci=self.max_loci,
            random_seed=self.locus_seed,
            min_length=self.min_length,
            windows=None,
            min_sample_coverage=len(self.imap),
            max_sample_missing=1.0,
            exclude=None,
            include_reference=False,
            imap=self.imap,
            minmap=self.minmap,
            concatenate=False,
            split=False,
            append_population=True,
            stdout=False,
            force=self.force,
            logged_command=None,
            cores=self.threads[0] if self.threads else 1,
            log_level=self.log_level,
        )
        loci = self.seqex.run()
        self.nloci = len(loci)
        self.effective_threads = self.threads
        if self.threads and self.threads[0] > self.nloci:
            self.effective_threads = (self.nloci, *self.threads[1:])
            logger.warning(
                "reducing BPP threads from {} to {} because only {} loci were written",
                self.threads[0],
                self.nloci,
                self.nloci,
            )

    def _write_mapfile(self) -> None:
        """Write the IMAP file expected by BPP."""
        longname = max(
            len(name)
            for species in self.species_order
            for name in self.imap[species]
        )
        formatstr = "{:<" + str(longname + 2) + "} {}"
        rows = [
            formatstr.format(sample, f"{DELIM}{species}")
            for species in self.species_order
            for sample in self.imap[species]
        ]
        self.paths.mapfile.write_text("\n".join(rows), encoding="utf-8")

    def _write_constraint_file(self) -> None:
        """Write a validated, namespace-adjusted constraint file when requested."""
        if self.constraint_file is not None:
            self.paths.constraintfile.write_text(
                self._render_constraint_text(), encoding="utf-8"
            )

    def _render_ctl(self) -> str:
        """Render one BPP 4.8.7 control file."""
        prefixed_species = self._prefixed_species_names()
        species_counts = " ".join(str(len(self.imap[name])) for name in self.species_order)
        phase = " ".join(["1"] * len(self.species_order))
        delimitation = (
            "1 " + " ".join(self.delimitation_proposal)
            if self.speciesdelimitation
            else "0"
        )
        lines = [
            "* I/O",
            f"seqfile = {self.paths.seqfile}",
            f"Imapfile = {self.paths.mapfile}",
            f"jobname = {self.paths.jobname}",
            "",
            "* DATA",
            f"nloci = {self.nloci}",
            "usedata = 1",
            f"cleandata = {int(self.clean_data)}",
            "",
            "* MODEL",
            f"speciesdelimitation = {delimitation}",
            f"speciestree = {int(self.speciestree)}",
        ]
        if self.species_model_prior is not None:
            lines.append(f"speciesmodelprior = {self.species_model_prior}")
        lines.extend(
            [
                f"species&tree = {len(prefixed_species)} {' '.join(prefixed_species)}",
                f"               {species_counts}",
                f"               {self._render_tree_text()}",
                f"phase = {phase}",
                f"model = {self.substitution_model}",
                "heredity = 0",
                f"thetamodel = {self.theta_model}",
                f"geneflow = {int(self.geneflow)}",
            ]
        )
        if self.constraint_file is not None:
            lines.append(f"constraintfile = {self.paths.constraintfile}")
        lines.extend(
            [
                "",
                "* PRIORS AND RATES",
                f"thetaprior = {' '.join(self.theta_prior)}",
                f"tauprior = {' '.join(self.tau_prior)}",
            ]
        )
        if self.alpha_prior is not None:
            lines.append(
                "alphaprior = {} {} {}".format(
                    f"{self.alpha_prior[0]:g}",
                    f"{self.alpha_prior[1]:g}",
                    self.alpha_prior[2],
                )
            )
        lines.extend(
            [
                f"locusrate = {' '.join(self.locus_rate)}",
                f"clock = {' '.join(self.clock)}",
            ]
        )
        if self.analysis == "msc-i":
            lines.append(f"phiprior = {self.phi_prior[0]:g} {self.phi_prior[1]:g}")
        if self.analysis == "msc-m":
            lines.append(f"wprior = {self.w_prior[0]:g} {self.w_prior[1]:g}")
            lines.append(f"migration = {len(self.migration_pairs)}")
            lines.extend(
                "  {} {}".format(
                    self._label_for_tree_reference(source),
                    self._label_for_tree_reference(target),
                )
                for source, target in self.migration_pairs
            )
        lines.extend(
            [
                "",
                "* MCMC PARAMS",
                f"seed = {-1 if self.mcmc_seed is None else self.mcmc_seed}",
                "finetune = 1",
                "print = 1 {} 0 {} {}".format(
                    int(self.print_locus_rates),
                    int(self.print_gene_trees),
                    int(self.print_substitution_parameters),
                ),
                f"burnin = {self.burnin}",
                f"sampfreq = {self.samplefreq}",
                f"nsample = {self.nsample}",
            ]
        )
        if self.checkpoint:
            lines.append(
                f"checkpoint = {' '.join(str(value) for value in self.checkpoint)}"
            )
        if self.effective_threads:
            lines.append(
                f"threads = {' '.join(str(value) for value in self.effective_threads)}"
            )
        return "\n".join(lines) + "\n"

    def _write_ctlfile(self) -> None:
        """Write the BPP control file."""
        self.paths.ctlfile.write_text(self._render_ctl(), encoding="utf-8")

    def write_inputs(self) -> SimpleNamespace:
        """Write the sequence, IMAP, constraint, and control files."""
        self._ensure_output_paths(include_results=False)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self._write_seqfile()
        self._write_mapfile()
        self._write_constraint_file()
        self._write_ctlfile()
        return self.paths

    def run(self) -> SimpleNamespace:
        """Write inputs, resolve a supported BPP binary, and execute the run."""
        self._ensure_output_paths(include_results=True)
        self.write_inputs()
        binary = _resolve_bpp_binary(self.binary)
        version = _resolve_bpp_binary_version(binary)
        if version is None:
            raise IPyradError(
                "Could not determine the BPP version; BPP >= 4.8.7 is required."
            )
        version_tuple = tuple(int(value) for value in version.split("."))
        if version_tuple < _MINIMUM_BPP_VERSION:
            raise IPyradError(
                f"BPP >= {BPP_DOCS_VERSION} is required; found {version}."
            )
        logger.info("bpp v{} ({})", version, binary)
        logger.info("running BPP and writing MCMC samples to {}", self.paths.mcmcfile)
        _call_bpp(binary, str(self.paths.ctlfile), self.algorithm)
        return self.paths


def run_bpp_method(
    *,
    data,
    name: str,
    outdir,
    tree,
    imap,
    minmap,
    max_loci: int,
    min_length: int,
    analysis: str,
    migration,
    theta_prior,
    tau_prior,
    species_model_prior,
    delimitation_proposal,
    phi_prior,
    w_prior,
    alpha_prior,
    substitution_model,
    locus_rate,
    clock,
    theta_model,
    clean_data,
    checkpoint,
    constraint_file,
    print_locus_rates,
    print_gene_trees,
    print_substitution_parameters,
    burnin: int,
    samplefreq: int,
    nsample: int,
    threads,
    locus_seed,
    mcmc_seed,
    binary,
    write_only: bool,
    force: bool,
    log_level: str = "INFO",
) -> None:
    """CLI entrypoint for one BPP run."""
    tool = Bpp(
        data=data,
        name=name,
        outdir=outdir,
        tree=tree,
        imap=imap,
        minmap=minmap,
        max_loci=max_loci,
        min_length=min_length,
        analysis=analysis,
        migration=migration,
        theta_prior=theta_prior,
        tau_prior=tau_prior,
        species_model_prior=species_model_prior,
        delimitation_proposal=delimitation_proposal,
        phi_prior=phi_prior,
        w_prior=w_prior,
        alpha_prior=alpha_prior,
        substitution_model=substitution_model,
        locus_rate=locus_rate,
        clock=clock,
        theta_model=theta_model,
        clean_data=clean_data,
        checkpoint=checkpoint,
        constraint_file=constraint_file,
        print_locus_rates=print_locus_rates,
        print_gene_trees=print_gene_trees,
        print_substitution_parameters=print_substitution_parameters,
        burnin=burnin,
        samplefreq=samplefreq,
        nsample=nsample,
        threads=threads,
        locus_seed=locus_seed,
        mcmc_seed=mcmc_seed,
        binary=binary,
        force=force,
        log_level=log_level,
    )
    if write_only:
        paths = tool.write_inputs()
        logger.info("wrote BPP inputs to {}", paths.ctlfile.parent)
        return
    paths = tool.run()
    logger.info("wrote BPP outputs to {}", paths.ctlfile.parent)

def _call_bpp(binary, ctlfile, alg):
    """Run one BPP job inside the control file directory and surface errors."""
    ctlpath = Path(ctlfile)
    workdir = ctlpath.parent
    cmd = [binary, "--cfile", ctlfile]
    proc = sps.run(
        cmd,
        cwd=workdir,
        stdout=sps.PIPE,
        stderr=sps.STDOUT,
        check=False,
    )
    if proc.returncode:
        message = proc.stdout.decode("utf-8", errors="replace")
        raise IPyradError(f"BPP failed with exit code {proc.returncode}:\n{message}")

    figfile = workdir / "FigTree.tre"
    if figfile.exists():
        os.replace(figfile, ctlpath.with_suffix("").with_suffix(".figtree.nex"))

    seed_used = workdir / "SeedUsed"
    if seed_used.exists():
        seed_used.unlink()
