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

from ..extracters.locus_extracter import LocusExtracter
from ..extracters.sequence_common import normalize_sequence_population_inputs
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
BPP_DOCS_VERSION = "4.8.6"
_DEFAULT_LOCUSRATE = ("1", "2", "3", "2", "iid")
_DEFAULT_CLOCK = ("2", "10.0", "100.0", "5.0", "dir", "LN")
_BPP_BINARY_SPECS = {
    ("linux", "x86_64"): (
        "4.8.6",
        "bpp-4.8.6-linux-x86_64.tar.gz",
        "bpp-4.8.6-linux-x86_64/bin/bpp",
    ),
    ("linux", "aarch64"): (
        "4.8.6",
        "bpp-4.8.6-linux-aarch64.tar.gz",
        "bpp-4.8.6-linux-aarch64/bin/bpp",
    ),
    ("darwin", "arm64"): (
        "4.8.6",
        "bpp-4.8.6-macos-aarch64.tar.gz",
        "bpp-4.8.6-macos-aarch64/bin/bpp",
    ),
    ("darwin", "x86_64"): (
        "4.8.4",
        "bpp-4.8.4-macos-x86_64.tar.gz",
        "bpp-4.8.4-macos-x86_64/bin/bpp",
    ),
    ("win32", "x86_64"): (
        "4.8.6",
        "bpp-4.8.6-win-x86_64.zip",
        "bpp-4.8.6-win-x86_64/bpp.exe",
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
    """Normalize two-number priors used by the CLI."""
    if values is None or len(values) != 2:
        raise IPyradError(f"{name} must contain exactly two numeric values.")
    return (
        _coerce_positive_number(values[0], f"{name}[0]"),
        _coerce_positive_number(values[1], f"{name}[1]"),
    )


def _normalize_alpha_prior(values) -> tuple[float, float, int]:
    """Normalize `alphaprior = alpha beta ncat`."""
    if values is None or len(values) != 3:
        raise IPyradError("alphaprior must contain exactly three values.")
    return (
        _coerce_positive_number(values[0], "alphaprior[0]"),
        _coerce_positive_number(values[1], "alphaprior[1]"),
        _coerce_positive_int(values[2], "alphaprior[2]"),
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


def _normalize_seed(seed) -> int | None:
    """Normalize CLI/API seed inputs."""
    if seed is None:
        return None
    if isinstance(seed, str):
        if seed.lower() == "none":
            return None
    return _coerce_positive_int(seed, "seed")


def _normalize_string_tokens(name: str, values, default: tuple[str, ...]) -> tuple[str, ...]:
    """Normalize one free-form ctl token list."""
    if values is None:
        return default
    if not values:
        raise IPyradError(f"{name} must not be empty.")
    return tuple(str(value) for value in values)


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
    """Prepare one BPP run, write one control-file set, and optionally execute it."""

    def __init__(
        self,
        *,
        data,
        name: str = "bpp",
        outdir: Path | str = "output-bpp",
        tree=None,
        imap=None,
        minmap=None,
        max_loci: int = 100,
        min_length: int = 100,
        msc_i: bool = False,
        msc_m=None,
        speciestree: bool = False,
        speciesdelimitation: bool = False,
        thetaprior=(3.0, 0.03),
        tauprior=(3.0, 0.03),
        speciesmodelprior: int = 0,
        phiprior=(1.0, 1.0),
        wprior=(2.0, 200.0),
        alphaprior=(1.0, 1.0, 4),
        locusrate=None,
        clock=None,
        burnin: int = 1000,
        samplefreq: int = 2,
        nsample: int = 10000,
        threads=None,
        seed=None,
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

        self.max_loci = _coerce_positive_int(max_loci, "max_loci")
        self.min_length = _coerce_positive_int(min_length, "min_length")
        self.burnin = _coerce_positive_int(burnin, "burnin")
        self.samplefreq = _coerce_positive_int(samplefreq, "samplefreq")
        self.nsample = _coerce_positive_int(nsample, "nsample")
        self.speciesmodelprior = int(speciesmodelprior)
        self.thetaprior = _normalize_prior_pair("thetaprior", thetaprior)
        self.tauprior = _normalize_prior_pair("tauprior", tauprior)
        self.phiprior = _normalize_prior_pair("phiprior", phiprior)
        self.wprior = _normalize_prior_pair("wprior", wprior)
        self.alphaprior = _normalize_alpha_prior(alphaprior)
        self.locusrate = _normalize_string_tokens("locusrate", locusrate, _DEFAULT_LOCUSRATE)
        self.clock = _normalize_string_tokens("clock", clock, _DEFAULT_CLOCK)
        self.threads = _normalize_threads(threads)
        self.seed = _normalize_seed(seed)

        self.msc_i = bool(msc_i)
        self.msc_m_pairs = _normalize_migration_tokens(msc_m) if msc_m else ()
        self.speciestree = bool(speciestree)
        self.speciesdelimitation = bool(speciesdelimitation)

        self.imap, self.minmap = normalize_sequence_population_inputs(imap, minmap)
        self.tree_text = _read_tree_text(tree)
        self.tree = _parse_tree(self.tree_text)
        self.tip_labels = tuple(self.tree.get_tip_labels())
        self.species_order = tuple(sorted(self.imap)) if self.imap else ()

        self.paths = SimpleNamespace(
            seqfile=self.outdir / f"{self.name}.phy",
            mapfile=self.outdir / f"{self.name}.imapfile.txt",
            ctlfile=self.outdir / f"{self.name}.ctl.txt",
            statsfile=self.outdir / f"{self.name}.stats.txt",
            jobname=self.outdir / self.name,
            outfile=self.outdir / f"{self.name}.txt",
            mcmcfile=self.outdir / f"{self.name}.mcmc.txt",
            figtree=self.outdir / f"{self.name}.figtree.nex",
        )
        self.lex = None

        self._validate()

    @property
    def algorithm(self) -> str:
        """Return the BPP analysis code for the selected model class."""
        if self.speciestree:
            return "10"
        if self.speciesdelimitation:
            return "01"
        return "00"

    def _validate(self) -> None:
        """Validate user inputs before writing files or launching BPP."""
        if not self.data.exists():
            raise IPyradError(f"data file does not exist: {self.data}")
        if self.data.suffix != ".hdf5":
            raise IPyradError("'data' must be an ipyrad2 .hdf5 file.")
        if not self.imap:
            raise IPyradError("imap is required for BPP analyses.")
        if not isinstance(self.imap, dict):
            raise IPyradError("imap must resolve to a dictionary.")
        if self.minmap is not None and set(self.minmap) != set(self.imap):
            raise IPyradError("imap and minmap keys must match.")
        if set(self.species_order) != set(self.tip_labels):
            raise IPyradError(
                "IMAP keys must match guide tree tip names exactly.\n"
                f"imap={sorted(self.imap)}\n"
                f"tree={sorted(self.tip_labels)}"
            )
        if self.speciesmodelprior not in {0, 1, 2, 3}:
            raise IPyradError("speciesmodelprior must be one of 0, 1, 2, or 3.")
        if self.speciestree and self.speciesmodelprior not in {0, 1}:
            raise IPyradError(
                "speciestree analyses require speciesmodelprior 0 or 1."
            )
        selected = sum(
            [
                self.msc_i,
                bool(self.msc_m_pairs),
                self.speciestree,
                self.speciesdelimitation,
            ]
        )
        if selected > 1:
            raise IPyradError(
                "Choose only one of --msc-i, --msc-m, --speciestree, or --speciesdelimitation."
            )
        if self.threads and self.threads[0] > self.max_loci:
            raise IPyradError("threads cannot exceed max_loci.")
        if self.msc_i and "&phi=" not in self.tree_text:
            raise IPyradError(
                "MSC-I requires --tree to contain extended Newick introgression annotations."
            )
        if self.msc_m_pairs:
            all_labels = set(self.tip_labels)
            for node in self.tree.treenode.traverse():
                if node.is_leaf():
                    continue
                if node.name:
                    all_labels.add(str(node.name))
            missing = sorted(
                {
                    label
                    for pair in self.msc_m_pairs
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

    def _ensure_output_paths(self, *, include_results: bool) -> None:
        """Fail early when output files already exist and overwrite is disabled."""
        paths = [self.paths.seqfile, self.paths.mapfile, self.paths.ctlfile, self.paths.statsfile]
        if include_results:
            paths.extend([self.paths.outfile, self.paths.mcmcfile, self.paths.figtree])
        existing = next((path for path in paths if path.exists()), None)
        if existing is not None and not self.force:
            raise IPyradError(
                f"Output file already exists: {existing}. Use --force to overwrite."
            )

    def _prefixed_species_names(self) -> tuple[str, ...]:
        """Return the sorted species labels with the BPP delimiter prefix."""
        return tuple(f"{DELIM}{name}" for name in self.species_order)

    def _render_tree_text(self) -> str:
        """Return the guide tree text with tip labels prefixed for BPP."""
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

    def _write_seqfile(self) -> None:
        """Extract loci once and write the BPP-formatted sequence file."""
        self.lex = LocusExtracter(
            data=self.data,
            name=self.name,
            outdir=self.outdir,
            out_format="bpp",
            nloci=self.max_loci,
            min_length=self.min_length,
            windows=None,
            min_sample_coverage=len(self.imap),
            max_sample_missing=1.0,
            exclude=None,
            include_reference=False,
            imap=self.imap,
            minmap=self.minmap,
            stdout=False,
            force=self.force,
            random_seed=self.seed,
        )
        self.lex._DELIM = "^" + DELIM
        self.lex._run()

    def _write_mapfile(self) -> None:
        """Write the IMAP file expected by BPP."""
        longname = 0
        for key in self.species_order:
            for name in self.imap[key]:
                longname = max(longname, len(name))
        formatstr = "{:<" + str(longname + len(DELIM) + 2) + "} {}"
        rows = [
            formatstr.format(f"{DELIM}{sample}", f"{DELIM}{species}")
            for species in self.species_order
            for sample in self.imap[species]
        ]
        with open(self.paths.mapfile, "w", encoding="utf-8") as out:
            out.write("\n".join(rows))

    def _render_ctl(self) -> str:
        """Render one BPP control file."""
        prefixed_species = self._prefixed_species_names()
        species_counts = " ".join(str(len(self.imap[name])) for name in self.species_order)
        phase = " ".join(["1"] * len(self.species_order))
        lines = [
            "* I/O",
            f"seqfile = {self.paths.seqfile}",
            f"Imapfile = {self.paths.mapfile}",
            f"jobname = {self.paths.jobname}",
            "",
            "* DATA",
            f"nloci = {self.max_loci}",
            "usedata = 1",
            "cleandata = 0",
            "",
            "* MODEL",
            f"speciestree = {1 if self.speciestree else 0}",
            "speciesdelimitation = 1 0 2" if self.speciesdelimitation else "speciesdelimitation = 0",
            f"speciesmodelprior = {self.speciesmodelprior}",
            f"species&tree = {len(prefixed_species)} {' '.join(prefixed_species)}",
            f"               {species_counts}",
            f"               {self._render_tree_text()}",
            f"phase = {phase}",
            "heredity = 0",
            "thetamodel = linked-none",
            f"geneflow = {1 if (self.msc_i or self.msc_m_pairs) else 0}",
            "",
            "* PRIORS",
            f"thetaprior = invgamma {self.thetaprior[0]:g} {self.thetaprior[1]:g} E",
            f"tauprior = invgamma {self.tauprior[0]:g} {self.tauprior[1]:g}",
            f"alphaprior = {self.alphaprior[0]:g} {self.alphaprior[1]:g} {self.alphaprior[2]}",
            f"locusrate = {' '.join(self.locusrate)}",
            f"clock = {' '.join(self.clock)}",
        ]
        if self.msc_i:
            lines.append(f"phiprior = {self.phiprior[0]:g} {self.phiprior[1]:g}")
        if self.msc_m_pairs:
            lines.append(f"wprior = {self.wprior[0]:g} {self.wprior[1]:g}")
            lines.append(f"migration = {len(self.msc_m_pairs)}")
            lines.extend(
                "  {} {}".format(
                    self._label_for_tree_reference(source),
                    self._label_for_tree_reference(target),
                )
                for source, target in self.msc_m_pairs
            )
        lines.extend(
            [
                "",
                "* MCMC PARAMS",
                f"seed = {-1 if self.seed is None else self.seed}",
                "finetune = 1",
                "print = 1 0 0 1 0",
                f"burnin = {self.burnin}",
                f"sampfreq = {self.samplefreq}",
                f"nsample = {self.nsample}",
            ]
        )
        if self.threads:
            lines.append(f"threads = {' '.join(str(i) for i in self.threads)}")
        return "\n".join(lines) + "\n"

    def _write_ctlfile(self) -> None:
        """Write the BPP control file."""
        with open(self.paths.ctlfile, "w", encoding="utf-8") as out:
            out.write(self._render_ctl())

    def write_inputs(self) -> SimpleNamespace:
        """Write the seqfile, IMAP file, and ctl file for one BPP run."""
        self._ensure_output_paths(include_results=False)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self._write_seqfile()
        self._write_mapfile()
        self._write_ctlfile()
        return self.paths

    def run(self) -> SimpleNamespace:
        """Write inputs, resolve the BPP binary, and execute one BPP run."""
        self._ensure_output_paths(include_results=True)
        self.write_inputs()
        binary = _resolve_bpp_binary(self.binary)
        version = _resolve_bpp_binary_version(binary)
        if version:
            logger.info("bpp v{} ({})", version, binary)
        else:
            logger.info("bpp ({})", binary)
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
    msc_i: bool,
    msc_m,
    speciestree: bool,
    speciesdelimitation: bool,
    thetaprior,
    tauprior,
    speciesmodelprior: int,
    phiprior,
    wprior,
    alphaprior,
    locusrate,
    clock,
    burnin: int,
    samplefreq: int,
    nsample: int,
    threads,
    seed,
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
        msc_i=msc_i,
        msc_m=msc_m,
        speciestree=speciestree,
        speciesdelimitation=speciesdelimitation,
        thetaprior=thetaprior,
        tauprior=tauprior,
        speciesmodelprior=speciesmodelprior,
        phiprior=phiprior,
        wprior=wprior,
        alphaprior=alphaprior,
        locusrate=locusrate,
        clock=clock,
        burnin=burnin,
        samplefreq=samplefreq,
        nsample=nsample,
        threads=threads,
        seed=seed,
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

    if alg == "00":
        figfile = workdir / "FigTree.tre"
        if figfile.exists():
            os.replace(figfile, ctlpath.with_suffix("").with_suffix(".figtree.nex"))

    seed_used = workdir / "SeedUsed"
    if seed_used.exists():
        seed_used.unlink()
