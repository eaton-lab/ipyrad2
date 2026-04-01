#!/usr/bin/env python

"""BPP analysis helpers for writing input files, running jobs, and summarizing results."""

import os
import sys
import glob
import time
import copy
import tempfile
import re
import itertools
import platform
import shutil
import subprocess as sps
import tarfile
import zipfile
from pathlib import Path
import requests

import numpy as np
import pandas as pd

from types import SimpleNamespace
from ..extractors.locus_extractor import LocusExtractor
from ...utils.exceptions import IPyradError
from ...utils.parallel import run_with_pool
from ...utils.progress import ProgressBar

_MISSING_SCIPY = """
You are missing required packages to use ipa.bpp().
First run the following conda install command:

conda install scipy -c conda-forge
"""
_MISSING_TOYTREE = """
You are missing required packages to use ipa.bpp().
First run the following conda install command:

conda install toytree -c conda-forge
"""

try:
    import scipy.stats as ss
except ImportError:

    raise IPyradError(_MISSING_SCIPY)

try:
    import toytree
    # from toytree.utils import bpp2newick
except ImportError:
    raise IPyradError(_MISSING_TOYTREE)


DELIM = "___"
BPP_DOCS_VERSION = "4.8.6"
_BPP_BINARY_SPECS = {
    ("linux", "x86_64"): ("4.8.6", "bpp-4.8.6-linux-x86_64.tar.gz", "bpp-4.8.6-linux-x86_64/bin/bpp"),
    ("linux", "aarch64"): ("4.8.6", "bpp-4.8.6-linux-aarch64.tar.gz", "bpp-4.8.6-linux-aarch64/bin/bpp"),
    ("darwin", "arm64"): ("4.8.6", "bpp-4.8.6-macos-aarch64.tar.gz", "bpp-4.8.6-macos-aarch64/bin/bpp"),
    ("darwin", "x86_64"): ("4.8.4", "bpp-4.8.4-macos-x86_64.tar.gz", "bpp-4.8.4-macos-x86_64/bin/bpp"),
    ("win32", "x86_64"): ("4.8.6", "bpp-4.8.6-win-x86_64.zip", "bpp-4.8.6-win-x86_64/bpp.exe"),
}
_THETA_TAU_PRIOR_FAMILIES = {"gamma", "invgamma"}
_PHI_PRIOR_FAMILIES = {"beta"}


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


def _coerce_positive_number(value, label: str) -> float:
    """Parse one prior parameter as a positive float."""
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise IPyradError(f"{label} must be numeric.") from exc
    if parsed <= 0:
        raise IPyradError(f"{label} must be > 0.")
    return parsed


def _normalize_prior_spec(name: str, value) -> tuple[str, float, float]:
    """Normalize legacy and explicit prior tuples to `(family, a, b)`."""
    if not isinstance(value, tuple):
        raise IPyradError(f"{name} must be provided as a tuple.")

    if name in {"thetaprior", "tauprior"}:
        valid_families = _THETA_TAU_PRIOR_FAMILIES
        default_family = "invgamma"
    elif name == "phiprior":
        valid_families = _PHI_PRIOR_FAMILIES
        default_family = "beta"
    else:
        raise IPyradError(f"Unknown prior name: {name}")

    if len(value) == 2:
        family = default_family
        a, b = value
    elif len(value) == 3 and isinstance(value[0], str):
        family = value[0].lower()
        a, b = value[1:]
    else:
        raise IPyradError(
            f"{name} must be a 2-tuple `(a, b)` or a 3-tuple `(family, a, b)`."
        )

    if family not in valid_families:
        valid = ", ".join(sorted(valid_families))
        raise IPyradError(f"{name} prior family must be one of: {valid}")
    return family, _coerce_positive_number(a, f"{name}[0]"), _coerce_positive_number(b, f"{name}[1]")


def _format_prior_spec(prior: tuple[str, float, float], *, estimate_theta: bool = False) -> str:
    """Render a normalized prior tuple for the BPP control file."""
    family, a, b = prior
    rendered = f"{family} {a:g} {b:g}"
    if estimate_theta:
        rendered = f"{rendered} E"
    return rendered


def _draw_gamma_from_range(min_value: float, max_value: float) -> tuple[float, float]:
    """Approximate a bounded uncertainty range with a gamma(a, b-rate) prior."""
    min_value = _coerce_positive_number(min_value, "minimum value")
    max_value = _coerce_positive_number(max_value, "maximum value")
    if max_value <= min_value:
        raise IPyradError("maximum value must be greater than minimum value.")
    mean = (max_value + min_value) / 2.0
    var = ((max_value - min_value) ** 2) / 16.0
    return mean ** 2 / var, mean / var


def _sample_bpp_prior(prior: tuple[str, float, float], *, size: int, random_state) -> np.ndarray:
    """Sample from a normalized BPP prior specification."""
    family, a, b = prior
    if family == "gamma":
        return ss.gamma.rvs(a, scale=1 / b, size=size, random_state=random_state)
    if family == "invgamma":
        return ss.invgamma.rvs(a, scale=b, size=size, random_state=random_state)
    if family == "beta":
        return ss.beta.rvs(a, b, size=size, random_state=random_state)
    raise IPyradError(f"Unsupported prior family: {family}")


def _bpp_prior_pdf(prior: tuple[str, float, float], xvals: np.ndarray) -> np.ndarray:
    """Evaluate the BPP prior density across x-values."""
    family, a, b = prior
    if family == "gamma":
        return ss.gamma.pdf(xvals, a, scale=1 / b)
    if family == "invgamma":
        return ss.invgamma.pdf(xvals, a, scale=b)
    if family == "beta":
        return ss.beta.pdf(xvals, a, b)
    raise IPyradError(f"Unsupported prior family: {family}")


def _bpp_prior_xvals(prior: tuple[str, float, float], *, lower: float = 0.005, upper: float = 0.995, nvalues: int = 100) -> np.ndarray:
    """Return plotting x-values for a normalized BPP prior."""
    family, a, b = prior
    if family == "gamma":
        return np.linspace(ss.gamma.ppf(lower, a, scale=1 / b), ss.gamma.ppf(upper, a, scale=1 / b), nvalues)
    if family == "invgamma":
        return np.linspace(ss.invgamma.ppf(lower, a, scale=b), ss.invgamma.ppf(upper, a, scale=b), nvalues)
    if family == "beta":
        return np.linspace(ss.beta.ppf(lower, a, b), ss.beta.ppf(upper, a, b), nvalues)
    raise IPyradError(f"Unsupported prior family: {family}")


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


def _map_bpp_columns_to_node_ids(tree, columns) -> dict[str, int]:
    """Map BPP parameter column names onto MRCA node indexes in the guide tree."""
    mapped = {}
    for col in columns:
        tips = col.split(DELIM)[1:]
        if not tips:
            raise IPyradError(f"Cannot map BPP parameter without tip labels: {col}")
        mapped[col] = tree.get_mrca_node(*tips).idx
    return mapped


def _assign_tree_distances_from_divergence_row(tree, divergence_row: pd.Series) -> None:
    """Set branch lengths on a tree from absolute divergence times at each node."""
    for node in tree.treenode.traverse("postorder"):
        if node.is_leaf():
            # toytree nodes expose `.dist` as read-only; update the private
            # branch-length field on the copied tree we are annotating.
            node._dist = float(divergence_row.loc[node.up.idx])
        elif node.up:
            parent_time = float(divergence_row.loc[node.up.idx])
            node_time = float(divergence_row.loc[node.idx])
            node._dist = max(0.0, parent_time - node_time)


class Bpp(object):
    """
    Prepare BPP inputs, run replicate jobs, and summarize BPP outputs.

    Parameters:
    -----------
    name: str
        A name for this analysis object.

    data: str
        The path to a .hdf5 file produced by ipyrad.

    imap: dict
        A Python dictionary with 'species' names as keys, and lists of sample
        names for the values. Any sample that is not included in the imap
        dictionary will be filtered out of the data when converting the .loci
        file into the bpp formatted sequence file. Each species in the imap
        dictionary must also be present in the input 'guidetree'.

    guidetree: str
        A newick string species tree hypothesis [e.g., (((a,b),(c,d)),e);]
        All taxa in the imap dictionary must also be present in the guidetree.
        Tree can also be a filename of a newick string.

    load_existing_results: bool
        If True then any existing results files saved in the working
        directory and with the entered name will be loaded and attached 
        to this object. This is useful if returning to a notebook later and 
        you want to summarize results. 

    reps_sample_loci (bool):
        if True then when nloci is set this will randomly sample a 
        different set of N loci in each replicate, rather than sampling
        just the first N loci < nloci. 

    infer_sptree:
        Default=0, only infer parameters on a fixed species tree. If 1, then the
        input tree is treated as a guidetree and tree search is employed to find
        the best tree. The results will include support values for the inferred
        topology.

    infer_delimit:
        Default=0, no delimitation. If 1 then splits in the tree that separate
        'species' will be collapsed to test whether fewer species are a better
        fit to the data than the number in the input guidetree.

    infer_delimit_args:
        Species delimitation algorithm is a two-part tuple. The first value
        is the algorithm (0 or 1) and the second value is a tuple of arguments
        for the given algorithm. See other ctl files for examples of what the
        delimitation line looks like. This is where you can enter the params
        (e.g., alpha, migration) for the two different algorithms.
        For example, the following args would produce the following ctl lines:
         alg=0, epsilon=5
         > delimit_alg = (0, 5)
         speciesdelimitation = 1 0 5

         alg=1, alpha=2, migration=1
         > delimit_alg = (1, 2, 1)
         speciesdelimitation = 1 1 2 1

         alg=1, alpha=2, migration=1, diagnosis=0, ?=1
         > delimit_alg = (1, 2, 1, 0, 1)
         speciesdelimitation = 1 1 2 1 0 1

    nloci (int):
        The max number of loci that will be used in an analysis.
    seed:
        A random number seed at start of analysis.
    burnin:
        Number of burnin generations in mcmc
    nsample:
        Number of mcmc generations to run.
    sampfreq:
        How often to sample from the mcmc chain.
    thetaprior:
        Prior on theta (4Neu). Legacy `(a, b)` tuples are treated as inverse-gamma;
        explicit tuples like `("invgamma", a, b)` and `("gamma", a, b)` are also
        accepted.
    tauprior
        Prior on root tau. Legacy `(a, b)` tuples are treated as inverse-gamma;
        explicit tuples like `("invgamma", a, b)` and `("gamma", a, b)` are also
        accepted.
    usedata:
        If false inference proceeds without sequence data (can be used to test
        the effect of priors on the tree distributions).
    cleandata:
        If 1 then sites with missing or hetero characters are removed.
    finetune:
        See bpp documentation.

    """    

    # init object for params
    def __init__(
        self,
        name,
        data=None,
        workdir="analysis-bpp", 
        guidetree=None, 
        imap=None, 
        minmap=None,
        nloci=100,
        length=1000,
        windows=None,
        max_sample_missing=1.0,
        reps_resample_loci=False,
        # load_existing_results=False,
        cores=4,
        log_level="DEBUG",
        *args, 
        **kwargs):

        # results files
        #self.files = Params()
        self.files = {"mcmcfiles":[], "outfiles":[], "treefiles":[]}
        self.files = SimpleNamespace(**self.files)
        self.files.mcmcfiles = []
        self.files.outfiles = []
        self.files.treefiles = []

        # store args
        self.data = data
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.guidetree = guidetree
        self.imap = imap
        self.minmap = minmap
        self.nloci = nloci
        self.length = length
        self.windows = windows
        self.max_sample_missing = max_sample_missing
        self.reps_resample_loci = reps_resample_loci
        # self.load_existing_results = load_existing_results
        self.cores = cores
        self.log_level = log_level

        # update kwargs 
        self.asyncs = []
        self.kwargs = {
            "binary": None,
            "nloci": 10,
            "infer_sptree": 0,
            "infer_delimit": 0,
            "infer_delimit_args": (0, 2),
            "speciesmodelprior": 1,
            "seed": 12345,
            "burnin": 10000,
            "nsample": 100000,
            "sampfreq": 2,
            "thetaprior": ("invgamma", 3, 0.002),
            "tauprior": ("invgamma", 3, 0.002),
            "phiprior": ("beta", 1, 1),
            "usedata": 1,
            "cleandata": 0,
            "finetune": 1,
            "copied": False,
        }

        # can set prior for visual plotting and/or for real analysis
        self._check_kwargs(kwargs)
        self.kwargs.update(kwargs)

        # binary is needed for running or loading and combining results
        self._check_binary()

        # update and check kwargs if data else do nothing which allows 
        # loading dummy objects for summarizing existing results.
        if data:
            self._check_args()


    def _check_binary(self):
        """Resolve the BPP binary from kwargs, PATH, or the bundled fallback."""
        if not sys.modules.get("toytree"):
            raise ImportError(_MISSING_TOYTREE)
        binary = self.kwargs.get("binary")
        if binary:
            resolved = os.path.realpath(os.path.expanduser(str(binary)))
            if not os.path.isfile(resolved):
                raise IPyradError(f"BPP binary does not exist: {resolved}")
            if not os.access(resolved, os.X_OK):
                raise IPyradError(f"BPP binary is not executable: {resolved}")
            self.kwargs["binary"] = resolved
            return

        found = shutil.which("bpp")
        if found:
            self.kwargs["binary"] = os.path.realpath(found)
            return

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
        self.kwargs["binary"] = os.path.realpath(spec.binary_path)


    def _check_kwargs(self, kwargs):
        """Validate and normalize user-supplied kwargs for the BPP runner."""
        for kwarg in list(kwargs):
            if kwarg not in self.kwargs:
                raise IPyradError(
                    "argument {} is either incorrect or no longer supported; "
                    "please check the latest documentation".format(kwarg)
                )

            if kwarg in ["nloci", "burnin", "sampfreq", "nsample", "seed"]:
                kwargs[kwarg] = int(kwargs[kwarg])

            if kwarg in ["thetaprior", "tauprior", "phiprior"]:
                kwargs[kwarg] = _normalize_prior_spec(kwarg, kwargs[kwarg])


    def _check_args(self):
        """Validate the HDF5 input, guide tree, IMAP mapping, and working directory."""
        # expand path to data
        self.data = os.path.realpath(os.path.expanduser(self.data))

        # check for data input
        if not os.path.exists(self.data):
            raise IPyradError(f"data file does not exist: {self.data}")
        if '.hdf5' not in self.data:
            raise IPyradError(
                "'data' argument must be an ipyrad2 .hdf5 file.")

        # set the guidetree
        if not self.guidetree:
            raise IPyradError(
                "must enter a 'guidetree' argument (a newick file or string).")
        self.tree = toytree.tree(self.guidetree)

        # check workdir
        if not self.workdir:
            self.workdir = os.path.join(os.path.curdir, "analysis-bpp")
        self.workdir = os.path.realpath(os.path.expanduser(self.workdir))           
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # parsing imap dictionary, or create simple 1-1 mapping
        if not self.imap:
            raise IPyradError("no IMAP dictionary.")
            # self.imap = {i: [i] for i in self.tree.get_tip_labels()}

        # check that all tree tip labels are represented in imap
        itips = set(self.imap.keys())
        ttips = set(self.tree.get_tip_labels())
        if itips.difference(ttips):
            raise IPyradError(
                "guidetree tips not in IMAP dictionary: {}"
                .format(itips.difference(ttips)))
        if ttips.difference(itips):
            raise IPyradError(
                "IMAP keys not in guidtree: {}"
                .format(ttips.difference(itips)))

        # check that minmap is OK or set it.
        # ...

        # checks
        if not isinstance(self.imap, dict):
            raise IPyradError("you must enter an IMAP dictionary")
        if set(self.imap.keys()) != set(self.tree.get_tip_labels()):
            raise IPyradError(
                "IMAP keys must match guidetree names: \n{}\n{}"
                .format(self.imap.keys(), self.tree.get_tip_labels())
            )


    def _load_existing_results(self, name, workdir, quiet=False):
        """
        Load existing results files for an object with this workdir and name. 
        This does NOT reload the parameter settings for the object...
        """
        # get mcmcs
        path = os.path.realpath(os.path.join(workdir, name))
        mcmcs = sorted(glob.glob("{}_r*.mcmc.txt".format(path)))
        # bpp outfile needs to be weeded out from the other *.txt files it creates
        outs = sorted(glob.glob("{}_r*.txt".format(path)))
        pattern = re.compile(rf"{name}_r(\d+)\.txt$")
        outs = [f for f in outs if pattern.search(f)]
        # get trees
        trees = sorted(glob.glob("{}_r*.tre".format(path)))

        for mcmcfile in mcmcs:
            if mcmcfile not in self.files.mcmcfiles:
                self.files.mcmcfiles.append(mcmcfile)
        for outfile in outs:
            if outfile not in self.files.outfiles:
                self.files.outfiles.append(outfile)
        for tree in trees:
            if tree not in self.files.treefiles:
                self.files.treefiles.append(tree)        

        if not quiet:
            print("[ipa.bpp] found {} existing result files".format(len(mcmcs)))


    @property
    def _algorithm(self):
        if self.kwargs["infer_sptree"]:
            if self.kwargs["infer_delimit"]:
                return "11"
            return "10"
        else:
            if self.kwargs["infer_delimit"]:
                return "01"
            return "00"


    def _run(self, force, nreps, dry_run):
        """Distribute BPP replicate jobs in parallel."""

        # clear out pre-existing files for this object
        self.files.mcmcfiles = []
        self.files.outfiles = []
        self.files.treefiles = []
        self.asyncs = []

        # apply locus extracter filtering
        self.lex = LocusExtractor(
            data=self.data, 
            name=self.name,
            outdir=self.workdir,
            imap=self.imap,
            out_format="bpp",
            nloci=self.nloci,
            windows=self.windows,
            min_length=self.length,
            minmap=self.minmap,
            min_sample_coverage=len(self.imap),  # ENFORCE at least 1 per spp.
            max_sample_missing=self.max_sample_missing,
            exclude=None,
            stdout=False,
            force=force,
        )
        self.lex._DELIM = "^" + DELIM

        version = _resolve_bpp_binary_version(self.kwargs["binary"])
        if version:
            print(f"[ipa bpp] bpp v{version} ({self.kwargs['binary']})")
        else:
            print(f"[ipa bpp] bpp ({self.kwargs['binary']})")

        # initiate random seed 
        np.random.seed(self.kwargs["seed"])

        printstr = "Writing {} bpp sequence files".format(nreps)
        prog = ProgressBar(nreps, None, printstr)
        prog.update()

        # replicate jobs
        jobs = {}
        for job in range(nreps):
            
            self.lex._run(postfix=str(job))
            prog.finished += 1
            prog.update()

            # make repname and make ctl filename
            self._name = "{}_r{}".format(self.name, job)
            ctlhandle = os.path.realpath(
                os.path.join(self.workdir, "{}.ctl.txt".format(self._name)))

            # skip if ctlfile exists
            if (not force) and (os.path.exists(ctlhandle)):
                print("Named ctl file already exists. Use force=True to"
                      " overwrite\nFilename:{}".format(ctlhandle))
                continue

            # submit job to run
            else:
                # write imap groupings to the imapfile
                self._write_mapfile()

                # change seed for each rep. CTL has other file paths.
                self._seed = np.random.randint(0, 1e9)
                ctlfile = self._write_ctlfile()

            kwargs = dict(binary=self.kwargs["binary"],
                          ctlfile=ctlfile,
                          alg=self._algorithm)
            jobs[job] = (_call_bpp, kwargs)

        if not dry_run:
            print("\n[ipa.bpp] distributing {} bpp jobs (name={}, nloci={})"
                  .format(nreps, self.name, self.nloci))
            if jobs:
                run_with_pool(jobs, self.log_level, self.cores, msg="Running bpp")

        else:
            # report on the files written
            print("\n[ipa.bpp] wrote {} bpp ctl files (name={}, nloci={})"
                  .format(nreps, self.name, self.nloci))


    def _write_mapfile(self):
        """
        Writes the IMAP formatted file for bpp from the ipa IMAP
        """
        # get outfile path
        self.mapfile = os.path.realpath(
            os.path.join(
                self.workdir, self._name + ".imapfile.txt"
            )
        )

        # get longest name in the file to create fmt string: e.g., '{:<10} {}'
        longname = 0
        for key in sorted(self.imap.keys()):
            for name in self.imap[key]:
                longname = max(longname, len(name))
        formatstr = "{:<" + str(longname + 2) + "} {}"

        # open handle for writing and add delimiter to names
        with open(self.mapfile, 'w') as mapfile:
            data = [
                formatstr.format(DELIM + val, DELIM + key) 
                for key in sorted(self.imap) 
                for val in self.imap[key]
            ]
            mapfile.write("\n".join(data))


    def _write_ctlfile(self):
        """ write outfile with any args in argdict """

        # get full path to out files for this repname
        jobname = os.path.realpath(os.path.join(self.workdir, self._name))
        mcmcfile = "{}.mcmc.txt".format(jobname)
        outfile = "{}.txt".format(jobname)

        # store files for this rep
        if mcmcfile not in self.files.mcmcfiles:
            self.files.mcmcfiles.append(mcmcfile)
        if outfile not in self.files.outfiles:
            self.files.outfiles.append(outfile)

        # get tree with delimiters on names
        tmptre = self.tree.copy()
        tmptre = tmptre.set_node_data(
            "name", 
            {i: DELIM + str(tmptre[i].name) for i in range(tmptre.ntips)},
        )

        # expand options to fill ctl file
        ctlstring = CTLFILE.format(**{
            "seqfile": self.lex.outfile,
            "mapfile": self.mapfile,
            "jobname": jobname,

            "nloci": self.nloci,
            "usedata": self.kwargs["usedata"],
            "cleandata": self.kwargs["cleandata"],

            "infer_sptree": int(self.kwargs["infer_sptree"]),
            "infer_delimit": int(self.kwargs["infer_delimit"]),
            "infer_delimit_args": (
                " ".join([str(i) for i in self.kwargs["infer_delimit_args"]])
                if self.kwargs["infer_delimit"]
                else ""),
            "nsp": len(self.imap),
            "spnames": " ".join([DELIM + i for i in sorted(self.imap)]),
            "spcounts": " ".join([str(len(self.imap[i])) for i in sorted(self.imap)]),
            "spnewick": tmptre.write(dist_formatter=None),
            "speciesmodelprior": self.kwargs["speciesmodelprior"],

            "thetaprior": _format_prior_spec(
                self.kwargs["thetaprior"],
                estimate_theta=self._algorithm == "00",
            ),
            "tauprior": _format_prior_spec(self.kwargs["tauprior"]),
            "phiprior": _format_prior_spec(self.kwargs["phiprior"]),

            "seed": self._seed,
            "finetune": self.kwargs["finetune"],
            "burnin": self.kwargs["burnin"],
            "sampfreq": self.kwargs["sampfreq"],
            "nsample": self.kwargs["nsample"],
        })

        # write out the ctl file
        ctlhandle = os.path.realpath(
            "{}.ctl.txt".format(os.path.join(self.workdir, self._name)))
        with open(ctlhandle, 'w') as out:
            out.write(ctlstring)

        return ctlhandle



    def copy(self, name):
        """ 
        Returns a copy of the bpp object with the same parameter settings
        but with the files.mcmcfiles and files.outfiles attributes cleared, 
        and with a new 'name' attribute. 

        Parameters
        ----------
        name (str):
            A name for the new copied bpp object that will be used for the 
            output files created by the object. 
        """
        if name == self.name:
            raise Exception(
                "new object must have a different 'name' than its parent")
        asyncs = self.asyncs
        self.asyncs = []
        newself = copy.deepcopy(self)
        newself.name = name
        newself.files.mcmcfiles = []
        newself.files.outfiles = []
        newself.asyncs = []
        self.asyncs = asyncs
        return newself



    def summarize_results(self, algorithm, individual_results=False, quiet=False):
        """ 
        Prints a summarized table of results from replicate runs, or,
        if individual_result=True, then returns a list of separate
        dataframes for each replicate run. 
        """
        # reports number of results found
        self._load_existing_results(self.name, self.workdir, quiet)
        if algorithm not in ["00", "01", "10", "11"]:
            raise IPyradError(f"Unsupported BPP algorithm: {algorithm}")
        if not quiet:
            print(
                "[ipa.bpp] summarizing algorithm '{}' results"
                .format(algorithm))

        # algorithms supported
        if algorithm == "00":
            return self._summarize_00(individual_results)
        if algorithm == "10":
            return self._summarize_10(individual_results)
        if algorithm == "01":
            return self._summarize_01(individual_results)
        raise IPyradError("Summary support for algorithm '11' is not yet implemented.")



    def _summarize_00(self, individual_results, ):
        """
        Combines MCMC files together and writes a ctl file then calls
        bpp with the print=-1 option which means "read in" so that it 
        will compute a new posterior table...
        """

        if not self.files.mcmcfiles:
            raise IPyradError("No result files found.")

        # load mcmc tables of posteriors
        dfs = [
            pd.read_csv(i, sep='\t', index_col=0)
            for i in self.files.mcmcfiles
        ]

        # return a list of parsed CSV results
        if individual_results:
            # load out tables of summarized posteriors
            if not self.files.outfiles:
                raise IPyradError("No BPP summary output files were found.")
            try:
                tables = []
                for ofile in self.files.outfiles:
                    tables.append(self._parse_A00_out(ofile))
            except (IndexError, FileNotFoundError, OSError, ValueError) as exc:
                raise IPyradError(
                    "BPP job cannot be summarized because it did not finish."
                ) from exc

            return tables, dfs

        # concatenate each CSV and then get stats w/ describe
        else:
            print('[ipa.bpp] combining mcmc files')

            # new file handles
            handle = self.files.mcmcfiles[0].rsplit("_", 1)[0] + "_concat"
            cf = handle + ".mcmc.txt"
            of = handle + ".txt"
            af = handle + ".conditional_a1b1.txt"

            # existing and new ctl files
            ctlfile = os.path.join(self.workdir, self.name + "_r0.ctl.txt")
            newctl = os.path.join(self.workdir, self.name + "_tmp.ctl.txt")

            # write a concatenated mcmc file
            concat_mcmc = pd.concat(dfs, ignore_index=True)
            # bpp expects the index col to be labeled 'Gen' and `ignore_index`
            # wipes that out here, so reset it
            concat_mcmc.index.name = "Gen"
            concat_mcmc.to_csv(cf, sep="\t", float_format="%.6f")

            # write a concatenated a1b1 file
            cdfs = []
            for mcmcfile in self.files.mcmcfiles:
                conditional = mcmcfile.replace("mcmc", "conditional_a1b1")
                if not os.path.exists(conditional):
                    raise IPyradError(f"Missing BPP conditional_a1b1 file: {conditional}")
                cdfs.append(pd.read_csv(conditional, sep='\t', index_col=0))
            # write a concatenated mcmc file
            concat_conditional = pd.concat(cdfs, ignore_index=True)
            # bpp expects the index col to be labeled 'Gen' and `ignore_index`
            # wipes that out here, so reset it
            concat_conditional.index.name = "Gen"
            concat_conditional.to_csv(af, sep="\t", float_format="%.6f")

            # write a tmp ctl file with print=-1 and mcmcfile=cf
            with open(ctlfile, 'r') as infile:
                with open(newctl, 'w') as outfile:
                    cdat = infile.readlines()
                    for line in cdat:
                        if 'jobname' in line:
                            line = "jobname = {}\n".format(handle)
                        if 'print' in line:
                            line = "print = -1 0 0 0 0\n"
                        outfile.write(line)

            # run bpp on the new ctlfile
            _call_bpp(self.kwargs["binary"], newctl, "00")

            # cleanup tmp file
            os.remove(newctl)

            # load the new table
            table = self._parse_A00_out(of)
            return table, concat_mcmc


    def _parse_A00_out(self, ofile):
        nodes = []
        rows = []
        in_nodes = False
        in_table = False
        with open(ofile, 'r') as infile:
            for line in infile:
                line = line.strip().split()

                # detect node label block
                if "(+1)" in line:
                    in_nodes = True
                elif in_nodes and line == []:
                    in_nodes = False
                elif in_nodes:
                    # Just save the column with the node name
                    nodes.append(line[3])

                # detect parameter estimate block
                if "param" in line and "rho1" in line:
                    rows.append(line)
                    in_table = True
                elif "lnL" in line:
                    rows.append(line)
                    in_table = False
                elif in_table:
                    rows.append(line)
            if len(rows) < 4 or not nodes:
                raise IPyradError(f"Failed to parse BPP A00 output: {ofile}")
            # remove blank line and '---' separator
            rows.pop(1)
            rows.pop(-2)

        # Increment index because bpp uses 1-based population indexing
        nodes = pd.DataFrame(nodes)
        nodes.index = [x+1 for x in nodes.index]

        # Load params data and fix the column labels
        params = pd.DataFrame(rows[1:], columns=rows[0]).T
        params.columns = params.iloc[0]
        params = params.iloc[1:]

        def _get_new_header(val):
            if val == "lnL":
                return val
            else:
                p, idx = val.split(":")
                return f"{p}_{idx}{nodes.loc[int(idx)][0]}"
        columns = [_get_new_header(x) for x in params.columns]
        params.columns = columns
        return params.astype(float)


    def _summarize_01(self, individual_results):
        dfs = []
        for ofile in self.files.outfiles:
            with open(ofile, 'r') as infile:
                dat = infile.read().split("posterior\n")[1]
                table, dat = dat.split("Order of ancestral nodes:")
                data = [i.strip().split() for i in table.strip().split("\n")]
                df = pd.DataFrame(
                    data=data,
                    columns=["x", "delim", "prior", "posterior"]
                )
                df = df.drop(columns=['x'])
                df["nspecies"] = [i.count("1") + 1 for i in df["delim"]]
                df["posterior"] = df["posterior"].astype(float)
                dfs.append(df)

        if individual_results:
            return dfs
        else:
            df = dfs[0]
            for odf in dfs[1:]:
                df["posterior"] += odf.posterior
            df["posterior"] /= len(dfs)
            return df



    def _summarize_10(self, individual_results):
        """
        Returns a tuple with (tree, treedist) where tree is a Toytree 
        containing the MJrule tree with support values on the edges and 
        treedist is a multitree containing the posterior distribution of trees. 
        """
        # store results 
        trees = []
        treelists = []

        # get best trees
        for treefile in self.files.outfiles:
            with open(treefile, 'r') as infile:

                # jump to end of file to get besttree
                line = None
                for line in infile:
                    pass
                if line is None:
                    raise IPyradError(f"BPP tree summary file was empty: {treefile}")
                newick = line.split(";")[0] + ";"

                # get majority-rule tree
                # while 1:
                #     line = next(infile)
                #     if line.startswith("(C) Majority"):
                #         newick = next(infile).strip()
                #         break

                # extract proper newick from bpp last line
                newick = newick.replace(" #", "")
                tree = toytree.tree(newick)

                # convert support values to ints
                for node in tree.treenode.traverse():
                    node.support = int(round(node.support * 100))
                trees.append(tree)

        # get posteriors
        for treefile in self.files.mcmcfiles:
            post = []
            with open(treefile, 'r') as infile:
                btrees = infile.readlines()
                for newick in btrees:
                    # newick = bpp2newick(tre)
                    post.append(newick)
                treelists.append(post)

        # return results
        if individual_results:
            return trees, [toytree.mtree(i) for i in treelists]
        else:
            return trees, toytree.mtree(list(itertools.chain(*treelists)))



    def draw_priors(self, gentime_min, gentime_max, mutrate_min, mutrate_max, invgamma=True, seed=123):
        """Draw the configured priors and their derived Ne/divergence distributions."""
        import toyplot
        del invgamma
        rng = np.random.default_rng(seed)

        # setup canvas
        canvas0 = toyplot.Canvas(width=800, height=250)
        ax0 = canvas0.cartesian(
            bounds=("10%", "30%", "20%", "80%"), 
            xlabel="prior on u (mut. rate x10^-8)",
        )
        ax1 = canvas0.cartesian(
            bounds=("40%", "60%", "20%", "80%"),
            xlabel="prior on theta (4Neu)",
        )
        ax2 = canvas0.cartesian(
            bounds=("70%", "90%", "20%", "80%"),
            xlabel="prior on Ne (theta/4u)",
        )
        for ax in (ax0, ax1, ax2):
            ax.y.ticks.labels.show = False
            ax.x.ticks.show = True
        ax0.y.label.text = "density"

        # distribution of mutation_rates ---------------------------------
        a, b = _draw_gamma_from_range(mutrate_min, mutrate_max)
        muts_rvs = ss.gamma.rvs(a, scale=1 / b, random_state=rng, size=1000)

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, scale=1 / b),
                ss.gamma.ppf(1 - edge, a, scale=1 / b),
                100)
            y = ss.gamma.pdf(x, a, scale=1 / b)
            ax0.fill(x * 1e8, y, opacity=0.33, color=toyplot.color.Palette()[0])

            if cix == 0.95:
                ax0.label.text = (
                    "95% CI: {:.1f} - {:.1f}"
                    .format(round(x[0] * 1e8, 1), round(x[-1] * 1e8, 1)))


        # distribution of prior on theta ---------------------------------
        theta_prior = self.kwargs["thetaprior"]
        theta_rvs = _sample_bpp_prior(theta_prior, size=1000, random_state=rng)

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = _bpp_prior_xvals(theta_prior, lower=edge, upper=1 - edge)
            y = _bpp_prior_pdf(theta_prior, x)
            ax1.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[1])

            if cix == 0.95:
                ax1.label.text = (
                    "95% CI: {:.3f} - {:.3f}"
                    .format(x[0], x[-1]))


        # distribution of effective population sizes --------------------
        ne_rvs = (theta_rvs / (muts_rvs * 4))
        mean, var, std = ss.bayes_mvs(ne_rvs)
        a = mean.statistic ** 2 / var.statistic
        b = mean.statistic / var.statistic

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, **{"scale": 1 / b}),
                ss.gamma.ppf(1 - edge, a, **{"scale": 1 / b}),
                100)
            y = ss.gamma.pdf(x, a, **{"scale": 1 / b})
            ax2.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[2])

            if cix == 0.95:
                ax2.label.text = (
                    "95% CI: {:.0f} - {:.0f}"
                    .format(x[0], x[-1]))


        # plot  ---------------------------------------------------------
        canvas1 = toyplot.Canvas(width=800, height=250)
        ax3 = canvas1.cartesian(
            bounds=("10%", "30%", "20%", "80%"), 
            xlabel="prior on generation times",
        )
        ax4 = canvas1.cartesian(
            bounds=("40%", "60%", "20%", "80%"),
            xlabel="prior on root tau (%)",
        )
        ax5 = canvas1.cartesian(
            bounds=("70%", "90%", "20%", "80%"),
            xlabel="prior on root Ma div (tau*gen/u)",            
        )
        for ax in (ax3, ax4, ax5):
            ax.y.ticks.labels.show = False
            ax.x.ticks.show = True
        ax3.y.label.text = "density"

        # distribution of generation times -------------------------------
        a, b = _draw_gamma_from_range(gentime_min, gentime_max)
        gens_rvs = ss.gamma.rvs(a, scale=1 / b, random_state=rng, size=1000)

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, scale=1 / b),
                ss.gamma.ppf(1 - edge, a, scale=1 / b),
                100)
            y = ss.gamma.pdf(x, a, scale=1 / b)
            ax3.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[0])

            if cix == 0.95:
                ax3.label.text = (
                    "95% CI: {:.1f} - {:.1f}"
                    .format(round(x[0], 1), round(x[-1], 1)))


        # distribution of prior on tau ----------------------------------
        tau_prior = self.kwargs["tauprior"]
        tau_rvs = _sample_bpp_prior(tau_prior, size=1000, random_state=rng)

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = _bpp_prior_xvals(tau_prior, lower=edge, upper=1 - edge)
            y = _bpp_prior_pdf(tau_prior, x)
            ax4.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[1])

            if cix == 0.95:
                ax4.label.text = (
                    "95% CI: {:.4f} - {:.4f}"
                    .format(x[0], x[-1]))


        # distribution of divergence times in years ----------------------
        div_rvs = (gens_rvs * tau_rvs) / muts_rvs
        div_rvs /= 1e6
        mean, var, std = ss.bayes_mvs(div_rvs)
        a = mean.statistic ** 2 / var.statistic
        b = mean.statistic / var.statistic

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, **{"scale": 1 / b}),
                ss.gamma.ppf(1 - edge, a, **{"scale": 1 / b}),
                100)
            y = ss.gamma.pdf(x, a, **{"scale": 1 / b})
            ax5.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[2])

            if cix == 0.95:
                ax5.label.text = (
                    "95% CI: {:.1f} - {:.1f}"
                    .format(x[0], x[-1]))

        return canvas0, canvas1#, (ax0, ax1, ax2, ax3, ax4, ax5)



    def get_transformed_values(self, mcmc, param, gentime_min, gentime_max, mutrate_min, mutrate_max):
        """
        Transform one posterior column into real units using empirical samples.
        """
        # init transformer tool
        tx = Transformer(mcmc, gentime_min, gentime_max, mutrate_min, mutrate_max)
        return tx.transform(param)



    def transform(self, mcmc, gentime_min, gentime_max, mutrate_min, mutrate_max, nsamp=1000):
        """
        Transform posterior theta/tau samples into Ne and divergence-time units.

        The BPP posterior samples are transformed directly rather than being
        re-fit to another parametric distribution first. Returned dataframes are
        keyed by guide-tree node index, and the returned trees use transformed
        divergence times as branch lengths.
        """
        # check that use supplied a tree to the bpp object
        if not hasattr(self, 'tree'):
            if not hasattr(self, 'guidetree'):
                raise IPyradError("requires a 'guidetree' to be set for the bpp object.")
            else:
                self.tree = toytree.tree(self.guidetree)

        tx = Transformer(mcmc, gentime_min, gentime_max, mutrate_min, mutrate_max)
        summary = pd.DataFrame(
            index=["mean", "median", "std", "min", "max", "2.5%", "97.5%"],
            columns=mcmc.columns,
            dtype=float,
        )
        transformed = {}
        for col in mcmc.columns:
            if col == "lnL":
                continue
            values = tx.transform(col)
            transformed[col] = values
            summary.loc["mean", col] = float(np.mean(values))
            summary.loc["median", col] = float(np.median(values))
            summary.loc["std", col] = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
            summary.loc["min", col] = float(np.min(values))
            summary.loc["max", col] = float(np.max(values))
            pc0, pc1 = np.percentile(values, [2.5, 97.5])
            summary.loc["2.5%", col] = float(pc0)
            summary.loc["97.5%", col] = float(pc1)

        tau_cols = [col for col in transformed if col.startswith("tau_")]
        theta_cols = [col for col in transformed if col.startswith("theta_")]
        tau_map = _map_bpp_columns_to_node_ids(self.tree, tau_cols)
        theta_map = _map_bpp_columns_to_node_ids(self.tree, theta_cols)

        divs = summary.loc[:, tau_cols].copy()
        divs.columns = [tau_map[col] for col in tau_cols]
        divs = divs.reindex(sorted(divs.columns), axis=1)

        popsize = summary.loc[:, theta_cols].copy()
        popsize.columns = [theta_map[col] for col in theta_cols]
        popsize = popsize.reindex(sorted(popsize.columns), axis=1)

        newtree = self.tree.copy()
        for node in newtree.treenode.traverse("postorder"):
            node.Ne = popsize.loc["median", node.idx] if node.idx in popsize.columns else 10000
        _assign_tree_distances_from_divergence_row(newtree, divs.loc["median"])

        mtree = toytree.mtree([newtree.write()] * nsamp)
        tau_frame = pd.DataFrame({tau_map[col]: transformed[col] for col in tau_cols})
        tau_frame = tau_frame.reindex(sorted(tau_frame.columns), axis=1)
        if tau_frame.shape[0] >= nsamp:
            sampled_tau = tau_frame.iloc[:nsamp].reset_index(drop=True)
        else:
            sampled_tau = tau_frame.sample(
                n=nsamp,
                replace=True,
                random_state=tx.seed,
            ).reset_index(drop=True)
        for tidx, tree in enumerate(mtree.treelist):
            _assign_tree_distances_from_divergence_row(tree, sampled_tau.loc[tidx])
        return divs, popsize, newtree, mtree



    def draw_posteriors(self, mcmc, gentime_min, gentime_max, mutrate_min, mutrate_max, invgamma=True, seed=123):
        """Draw posterior densities on top of the configured priors."""
        import toyplot
        del invgamma
        rng = np.random.default_rng(seed)
        tx = Transformer(mcmc, gentime_min, gentime_max, mutrate_min, mutrate_max, seed=seed)

        # setup canvas
        canvas0 = toyplot.Canvas(width=800, height=250)
        ax0 = canvas0.cartesian(
            bounds=("10%", "30%", "20%", "80%"), 
            xlabel="prior on u (mut. rate x10^-8)",
        )
        ax1 = canvas0.cartesian(
            bounds=("40%", "60%", "20%", "80%"),
            xlabel="prior on theta (4Neu)",
        )
        ax2 = canvas0.cartesian(
            bounds=("70%", "90%", "20%", "80%"),
            xlabel="prior on Ne (theta/4u)",
        )
        for ax in (ax0, ax1, ax2):
            ax.y.ticks.labels.show = False
            ax.x.ticks.show = True
        ax0.y.label.text = "density"

        # distribution of mutation_rates ---------------------------------
        a, b = _draw_gamma_from_range(mutrate_min, mutrate_max)
        muts_rvs = ss.gamma.rvs(a, scale=1 / b, random_state=rng, size=10000)

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, scale=1 / b),
                ss.gamma.ppf(1 - edge, a, scale=1 / b),
                100)
            y = ss.gamma.pdf(x, a, scale=1 / b)
            ax0.fill(x * 1e8, y, opacity=0.33, color=toyplot.color.Palette()[0])

            if cix == 0.95:
                ax0.label.text = (
                    "95% CI: {:.1f} - {:.1f}"
                    .format(round(x[0] * 1e8, 1), round(x[-1] * 1e8, 1)))


        # distribution of prior on theta ---------------------------------
        theta_prior = self.kwargs["thetaprior"]
        theta_rvs = _sample_bpp_prior(theta_prior, size=10000, random_state=rng)

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = _bpp_prior_xvals(theta_prior, lower=edge, upper=1 - edge)
            y = _bpp_prior_pdf(theta_prior, x)
            ax1.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[1])

            if cix == 0.95:
                ax1.label.text = (
                    "95% CI: {:.3f} - {:.3f}"
                    .format(x[0], x[-1]))


        # distribution of effective population sizes --------------------
        ne_rvs = (theta_rvs / (muts_rvs * 4))
        mean, var, std = ss.bayes_mvs(ne_rvs)
        a = mean.statistic ** 2 / var.statistic
        b = mean.statistic / var.statistic

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, **{"scale": 1 / b}),
                ss.gamma.ppf(1 - edge, a, **{"scale": 1 / b}),
                100)
            y = ss.gamma.pdf(x, a, **{"scale": 1 / b})
            ax2.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[2])

            if cix == 0.95:
                ax2.label.text = (
                    "95% CI: {:.0f} - {:.0f}"
                    .format(x[0], x[-1]))


        # -----------
        thetacols = [i for i in mcmc.columns if "theta_" in i] 
        for col in thetacols:
            theta = mcmc.loc[:, col]
            mags, edges = np.histogram(
                theta,
                bins=100, 
                range=(0, 0.03),
                density=True,
            )
            
            ax1.plot(edges[:-1][mags>0], mags[mags>0], color='black', opacity=0.5);
            ne_vals = tx.transform(col)
            nes = pd.Series(ne_vals).sample(
                min(ne_vals.size, 10000),
                replace=ne_vals.size < 10000,
                random_state=seed,
            )
            mags, edges = np.histogram(
                nes,
                bins=100, 
                range=(0, 1200000),
                density=True,
            )
            ax2.plot(edges[:-1][mags>0], mags[mags>0], color='black', opacity=0.5);


        # plot  ---------------------------------------------------------
        canvas1 = toyplot.Canvas(width=800, height=250)
        ax3 = canvas1.cartesian(
            bounds=("10%", "30%", "20%", "80%"), 
            xlabel="prior on generation times",
        )
        ax4 = canvas1.cartesian(
            bounds=("40%", "60%", "20%", "80%"),
            xlabel="prior on root tau (%)",
        )
        ax5 = canvas1.cartesian(
            bounds=("70%", "90%", "20%", "80%"),
            xlabel="prior on root Ma div (tau*gen/u)",            
        )
        for ax in (ax3, ax4, ax5):
            ax.y.ticks.labels.show = False
            ax.x.ticks.show = True
        ax3.y.label.text = "density"

        # distribution of generation times -------------------------------
        a, b = _draw_gamma_from_range(gentime_min, gentime_max)
        gens_rvs = ss.gamma.rvs(a, scale=1 / b, random_state=rng, size=10000)

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, scale=1 / b),
                ss.gamma.ppf(1 - edge, a, scale=1 / b),
                100)
            y = ss.gamma.pdf(x, a, scale=1 / b)
            ax3.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[0])

            if cix == 0.95:
                ax3.label.text = (
                    "95% CI: {:.1f} - {:.1f}"
                    .format(round(x[0], 1), round(x[-1], 1)))


        # distribution of prior on tau ----------------------------------
        tau_prior = self.kwargs["tauprior"]
        tau_rvs = _sample_bpp_prior(tau_prior, size=10000, random_state=rng)

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = _bpp_prior_xvals(tau_prior, lower=edge, upper=1 - edge)
            y = _bpp_prior_pdf(tau_prior, x)
            ax4.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[1])

            if cix == 0.95:
                ax4.label.text = (
                    "95% CI: {:.4f} - {:.4f}"
                    .format(x[0], x[-1]))


        # distribution of divergence times in years ----------------------
        div_rvs = (gens_rvs * tau_rvs) / muts_rvs
        div_rvs /= 1e6
        mean, var, std = ss.bayes_mvs(div_rvs)
        a = mean.statistic ** 2 / var.statistic
        b = mean.statistic / var.statistic

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, **{"scale": 1 / b}),
                ss.gamma.ppf(1 - edge, a, **{"scale": 1 / b}),
                100)
            y = ss.gamma.pdf(x, a, **{"scale": 1 / b})
            ax5.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[2])

            if cix == 0.95:
                ax5.label.text = (
                    "95% CI: {:.1f} - {:.1f}"
                    .format(x[0], x[-1]))

        # --------------------------------------------------------
        # POSTERIOR
        # 
        # get the deepest divergence time
        taus = sorted([i for i in mcmc.columns if "tau_" in i], key=lambda x: len(x))[-1]
        taus = mcmc.loc[:, taus]

        # had range(0, 0.03) as an arg, removed on 2020-11-25
        mags, edges = np.histogram(taus, bins=100, density=True)

        ax4.plot(edges[:-1][mags>0], mags[mags>0], color='black', opacity=0.5);
        div_vals = tx.transform(taus.name) / 1e6
        divs = pd.Series(div_vals).sample(
            min(div_vals.size, 10000),
            replace=div_vals.size < 10000,
            random_state=seed,
        )

        # had range(0, 75) as an arg, removed on 2020-11-25
        mags, edges = np.histogram(divs.to_list(), bins=100, density=True)
        ax5.plot(edges[:-1][mags>0], mags[mags>0], color='black', opacity=0.5);

        # get all Ne values
        return canvas0, canvas1



    def draw_posteriors_stacked(self, gamma_tuples, labels, **kwargs):
        """
        ...

        """
        import toyplot
        c = toyplot.Canvas(225, 350)
        ax = c.cartesian()

        vals = gamma_tuples
        vals = [(i[0] ** 2 / i[1], i[0] / i[1]) for i in vals]

        marks = []
        ax.hlines(len(vals))
        for vidx, val in enumerate(vals):

            idx = len(vals) - vidx - 1
            a, b = val
            x = np.linspace(
                ss.gamma.ppf(0.01, a, scale=1/b),
                ss.gamma.ppf(0.99, a, scale=1/b),
                100,
            )

            mark = ax.fill(
                x, 
                (ss.gamma.pdf(x, a, scale=1/b) / ss.gamma.pdf(x, a, scale=1/b).max()) * 1.25,
                style={"fill": toyplot.color.Palette()[0], "stroke": "white", "stroke-width": 1.5},
                baseline=np.repeat(idx, x.size),
                annotation=True,
            )
            ax.hlines(idx)
            marks.append(mark) 

        ax.hlines([-0.3, len(vals) + 0.35])
        minm = min([i.domain('x')[0] for i in marks])
        maxm = max([i.domain('x')[1] for i in marks])
        ax.vlines([minm - maxm * 0.05, maxm + maxm * 0.05])
        ax.x.ticks.locator = toyplot.locator.Extended(count=3, only_inside=True)
        ax.x.ticks.show = True
        ax.y.ticks.locator = toyplot.locator.Explicit(
            locations=np.arange(len(vals)) + 0.5,
            labels=labels,
        )
        return c, ax, marks



    def draw_posterior_tree(
        self,
        mcmc, 
        gentime_min, 
        gentime_max, 
        mutrate_min, 
        mutrate_max, 
        node_dists=None,
        nbins=50,
        ydrop=2,
        **kwargs):
        """
        

        Parameters
        ----------
        node_dists: list
            A list of node indices to show histogram distributions under.
        ydrop: int
            The amount of y space to leave for histograms.
        """
        import toyplot
        import toytree
        icolors = copy.deepcopy(toytree.icolors2)

        # get results as transformed dataframes and trees
        dfdiv, dfne, ttre, tmtre = self.transform(
            mcmc, 
            gentime_min, gentime_max,
            mutrate_min, mutrate_max,
        )

        # do not allow any tips in node_dists:
        if node_dists:
            for nidx in node_dists:
                if ttre.idx_dict[nidx].is_leaf():
                    raise IPyradError(
                        "error in node_dists: cannot plot div time for tip nodes")

        # setup plot dims
        height = (275 if "height" not in kwargs else kwargs["height"])
        width = (450 if "width" not in kwargs else kwargs["width"])
        canvas = toyplot.Canvas(height=height, width=width)
        axes = canvas.cartesian(yshow=False) 

        # draw the tree on canvas
        mark = ttre.draw(
            axes=axes,
            node_sizes=0,
            edge_type='c', 
            # node_labels="idx",
            layout='r',
            scalebar=True,
        );

        # add error bars at nodes
        colors = {}
        for node in ttre.treenode.traverse():
            if not node.is_leaf():
                if node.idx in node_dists:
                    colors[node.idx] = next(icolors)
                else:
                    colors[node.idx] = "grey"               
                axes.rectangle(
                    -dfdiv.loc["2.5%", node.idx],
                    -dfdiv.loc["97.5%", node.idx],
                    node.x - 0.15, 
                    node.x + 0.15, 
                    color=colors[node.idx],
                    opacity=0.45,
                )

                if node.idx in node_dists:
                    axes.scatterplot(
                        -dfdiv.loc["median", node.idx],
                        node.x,
                        marker="o",
                        size=8,
                        mstyle={"stroke": "black", "fill": colors[node.idx]},
                        opacity=0.45
                    )

        # add vertical lines at select nodes
        for nidx in node_dists:
            axes.plot(
                [-dfdiv.loc["median", nidx], -dfdiv.loc["median", nidx]],
                [-ydrop, ttre.get_node_coordinates()[nidx, 1]],
                color=colors[nidx], #'black',  #toytree.colors[1],
                opacity=0.45,
                style={"stroke-dasharray": "2,2"}
            )

        # add histograms under nodes
        hists = {}
        for nidx in node_dists:  #, nodename in enumerate(node_dists):

            # get tips desc from this node
            tips = self.tree.get_tip_labels(nidx)
            taus = [i for i in mcmc.columns if "tau_" in i]
            matches = []
            for i in taus:
                if all([t in i for t in tips]):
                    matches.append(i)
            match = sorted(matches, key=lambda x: len(x))[0]

            # get tranformed values for this specific node
            vals = self.get_transformed_values(
                mcmc,
                match,
                gentime_min, gentime_max,
                mutrate_min, mutrate_max,
            )

            # plot histogram of these values
            p0, p1 = np.percentile(vals, [2.5, 97.5])
            mags, edges = np.histogram(
                vals, 
                bins=nbins,
                range=(p0, p1),
                density=True,
            )
            hists[nidx] = (mags, edges)

        # get max mags
        maxmags = max([i[0].max() for i in hists.values()])
        maxmags *= 1.10
        
        # plot normalized histos maxed at maxmags
        for nidx in node_dists:        

            # reload data
            mags, edges = hists[nidx]

            # plot bars at midpoints
            nudge = ((edges[1] - edges[0]) / 2)
            axes.fill(
                -(edges[:-1] + nudge),
                (mags / maxmags) * ydrop,
                color=colors[nidx],
                opacity=0.25,
                baseline=[-ydrop] * nbins
            )

            axes.plot(
                a=-edges[1:], 
                b=((mags / maxmags) * ydrop) - ydrop,
                color=colors[nidx],  #toytree.colors[nidx],
                opacity=0.5,
            )

        # get tick spacers
        ndigits = len(str(dfdiv.loc["97.5%"].max().astype(int)))
        raised = ndigits - 2

        # style the axes
        axes.x.label.text = "time (x 10^{})".format(raised)
        marks = np.linspace(0, dfdiv.loc["97.5%"].max().astype(int), raised)
        ticks = marks / (10**raised)
        axes.x.ticks.locator = toyplot.locator.Explicit(
            [-1 * i for i in marks], 
            ["{:.0f}".format(i) for i in ticks.round(0)]
        )
        axes.padding = 10

        # SAVE THE FIGURE AND DISPLAY IT
        #toyplot.html.render(canvas, "/tmp/test.html")
        return canvas, axes
                





class Transformer(object):
    """Transform posterior theta/tau samples into Ne and divergence-time units."""
    def __init__(self, df, gentime_min, gentime_max, mutrate_min, mutrate_max, seed=123):

        self.df = df
        if self.df is None or self.df.empty:
            raise IPyradError("Cannot transform an empty BPP posterior table.")

        self.seed = int(seed)
        self._rng = np.random.default_rng(self.seed)
        self.gentime_a, self.gentime_b = _draw_gamma_from_range(gentime_min, gentime_max)
        self.mutrate_a, self.mutrate_b = _draw_gamma_from_range(mutrate_min, mutrate_max)
        self._sample_gentime_rvs()
        self._sample_mutrate_rvs()


    def _sample_gentime_rvs(self):
        """Sample generation times from the user-provided uncertainty range."""
        self.gentime_rvs = ss.gamma.rvs(
            self.gentime_a,
            scale=1 / self.gentime_b,
            random_state=self._rng,
            size=self.df.shape[0],
        )


    def _sample_mutrate_rvs(self):
        """Sample mutation rates from the user-provided uncertainty range."""
        self.mutrate_rvs = ss.gamma.rvs(
            self.mutrate_a,
            scale=1 / self.mutrate_b,
            random_state=self._rng,
            size=self.df.shape[0],
        )


    def _get_gentime_x(self, nvalues=100):
        xvals = np.linspace(
            ss.gamma.ppf(0.0001, self.gentime_a, scale=1 / self.gentime_b),
            ss.gamma.ppf(0.9999, self.gentime_a, scale=1 / self.gentime_b),
            nvalues,
        )
        return xvals


    def _get_mutrate_x(self, nvalues=100):
        xvals = np.linspace(
            ss.gamma.ppf(0.0001, self.mutrate_a, scale=1 / self.mutrate_b),
            ss.gamma.ppf(0.9999, self.mutrate_a, scale=1 / self.mutrate_b),
            nvalues,
        )
        return xvals


    def _get_parameter_values(self, colname):
        """Return one posterior parameter column plus the aligned uncertainty draws."""
        if colname not in self.df.columns:
            raise IPyradError(f"posterior column not found: {colname}")
        series = pd.to_numeric(self.df[colname], errors="coerce")
        valid = series.notna().to_numpy()
        if not np.any(valid):
            raise IPyradError(f"posterior column has no numeric values: {colname}")
        values = series.to_numpy(dtype=float, na_value=np.nan)[valid]
        return values, self.gentime_rvs[valid], self.mutrate_rvs[valid]


    def transform(self, colname):
        """Transform one posterior column to divergence-time or Ne units."""
        values, gentime, mutrate = self._get_parameter_values(colname)
        if "tau" in colname:
            return (values * gentime) / mutrate
        if "theta" in colname:
            return values / (mutrate * 4)
        raise IPyradError(f"Unsupported BPP posterior parameter: {colname}")


    # def transform_plot(self, axes=None, **kwargs):
    #     """

    #     """
    #     # get mean, var, std of the mcmc posterior values
    #     mean, var, std = ss.bayes_mvs(self.div_rvs)

    #     a = mean.statistic ** 2 / var.statistic
    #     b = mean.statistic / var.statistic
    #     x = np.linspace(
    #         ss.gamma.ppf(0.025, a, **{'scale': 1 / b}),
    #         ss.gamma.ppf(0.975, a, **{'scale': 1 / b})
    #     )
    #     pdf = ss.gamma.pdf(x, a, scale=1 / b)
    #     print("mean: {:.0f}".format(mean.statistic))
    #     print("95% CI: {:.0f}-{:.0f}".format(pdf[0], pdf[-1]))

    #     canvas, axes, mark = draw_dist(
    #         mean.statistic, 
    #         var.statistic, 
    #         xlabel="Divergence time", 
    #         axes=axes,
    #         **kwargs)


    #     if "theta" in colname:
    #         self._transform_theta(colname)
    #         mean, var, std = ss.bayes_mvs(self.ne_rvs)
    #         print("mean: {:.0f}".format(mean.statistic))
    #         print("95% CI: {:.0f}-{:.0f}".format(mean.minmax[0], mean.minmax[1]))
    #         canvas, axes, mark = draw_dist(
    #             mean.statistic, var.statistic, 
    #             "Effective population size", **kwargs)
    #     return canvas, axes, mark




    # def draw_dists(self, mcmcs):
    #     """

    #     """
    #     mcmcs = self.summarize_results("00")




def draw_dists(mcmcs, **kwargs):
    """
    Draw several posterior distributions on a shared axis.
    """
    import toyplot
    canvas = toyplot.Canvas(225, 50 * len(mcmcs))
    axes = canvas.cartesian()

    # store marks for each drawn distribution
    marks = []

    # draw each distribution
    for idx, dist in enumerate(mcmcs):

        # get mean, var, and a,b of gamma distribution
        mean, var, std = ss.bayes_mvs(dist)
        a = mean.statistic ** 2 / var.statistic
        b = mean.statistic / var.statistic

        # how far down from top
        tidx = len(mcmcs) - idx

        # x-axis points along distribution density
        x = np.linspace(
            ss.gamma.ppf(0.01, a, scale=1 / b),
            ss.gamma.ppf(0.99, a, scale=1 / b),
            100,
        )

        # get the density
        density = ss.gamma.pdf(x, a, scale=1 / b)

        # draw it
        mark = axes.fill(
            x, 
            (density / density.max()) * 1.25,
            baseline=np.repeat(tidx, x.size),
            style={
                "fill": "grey", 
                "stroke": "white", 
                "stroke-width": 2,
            }
        )
        marks.append(mark)  
        axes.hlines(idx)

    # style the axes
    axes.hlines([-0.25, len(mcmcs) + 0.35])
    minm = min([i.domain('x')[0] for i in marks])
    maxm = max([i.domain('x')[1] for i in marks])    
    axes.vlines([minm - minm * 0.05, maxm + maxm * 0.05])
    axes.x.ticks.locator = toyplot.locator.Extended(4, only_inside=True)
    axes.x.ticks.show = True
    axes.y.ticks.locator = toyplot.locator.Explicit(
        locations=np.arange(len(mcmcs)) + 0.5,
        labels=range(len(mcmcs)),
    )
    return canvas, axes, marks


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

    # BPP writes side-effect files next to the working directory rather than to
    # the configured jobname path, so keep the subprocess cwd local to the
    # control file and then rename the tree artifact into the normal prefix.
    if alg == "00":
        figfile = workdir / "FigTree.tre"
        newfigpath = ctlpath.with_suffix("").with_suffix(".figtree.nex")
        if figfile.exists():
            os.replace(figfile, newfigpath)

    seed_used = workdir / "SeedUsed"
    if seed_used.exists():
        seed_used.unlink()


def draw_dist(mean, var, xlabel=None, axes=None, **kwargs):
    """
    Tranformer class subfunction to draw posterior density.
    """
    import toyplot
    a = mean ** 2 / var
    b = mean / var

    # set up plots
    canvas = None
    if axes is None:
        canvas = toyplot.Canvas(width=400, height=300)
        axes = canvas.cartesian(ylabel="density", xlabel=xlabel)

    # default kwargs
    if "color" not in kwargs:
        kwargs["color"] = toyplot.color.Palette()[0]
    if "opacity" not in kwargs:
        kwargs["opacity"] = 0.25


    # confidence intervals shaded
    for civ in [0.99, 0.95, 0.50]:

        # stroke the wide interval only
        style = ({"stroke": "black", "stroke-width": 2} if civ == 0.99 else {})

        # get 100 values evenly spaced across 99% 
        edge = (1 - civ) / 2.
        x = np.linspace(
            ss.gamma.ppf(edge, a, **{'scale': 1 / b}),
            ss.gamma.ppf(1 - edge, a, **{'scale': 1 / b}), 
            100)

        # plot values across range of gamma
        mark = axes.fill(
            x,  # / 1e6,
            ss.gamma.pdf(x, a, **{'scale': 1 / b}),
            style=style,
            **kwargs
        )

    axes.y.ticks.labels.show = False
    axes.x.ticks.show = True
    return canvas, axes, mark


CTLFILE = """
* I/O 
seqfile = {seqfile}
Imapfile = {mapfile}
jobname = {jobname}

* DATA
nloci = {nloci}
usedata = {usedata}
cleandata = {cleandata}

* MODEL
speciestree = {infer_sptree}
speciesdelimitation = {infer_delimit} {infer_delimit_args}
speciesmodelprior = {speciesmodelprior}
species&tree = {nsp} {spnames}
               {spcounts}
               {spnewick}

* PRIORS
thetaprior = {thetaprior}
tauprior = {tauprior}
phiprior = {phiprior}

* MCMC PARAMS
seed = {seed}
finetune = {finetune}
print = 1 0 0 0 0
burnin = {burnin}
sampfreq = {sampfreq}
nsample = {nsample}
"""

