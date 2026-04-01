# Installation

The recommended way to install ipyrad2 is to use ``conda``, which will not only install ``ipyrad2`` but also all of its required external dependencies: `fastp`, `samtools`, `bcftools`, `bedtools`, `bwa-mem2`, `vsearch`, and `mafft`.

## Recommended Install

```bash
conda install ipyrad2 -c conda-forge -c bioconda
```

## For developers: installing a development version

Set up an environment with all dependencies
```bash
git clone https://github.com/eaton-lab/ipyrad2.git
cd ipyrad2
conda env create -f environment.yml -n ipyrad2
conda activate ipyrad2
```

Then install ipyrad2 itself:

```bash
pip install -e . --no-deps
```

## Optional Analysis Dependencies

The core environment is enough for the main assembly workflow, but some downstream analysis commands require extra Python packages that are listed as optional dependencies in `pyproject.toml`.

The main optional analysis packages are:

- `scikit-learn`
- `umap-learn`
- `toyplot`
- `toytree`

You can install them either with conda:

```bash
conda install scikit-learn umap-learn toyplot toytree -c conda-forge
```
```bash
conda install ipyrad2[analysis] -c conda-forge
```

## Verify the Install

A minimal check is:

```bash
ipyrad2 -h
```

If you are a returning ipyrad user, note that the source install also exposes a legacy-style entrypoint:

```bash
ipyrad2-classic -h
```

That command exists for the older parameter-file workflow, but the main documentation in this site focuses on the newer `ipyrad2` subcommand interface.

## Common Pitfalls

### Forgetting to activate the environment

If `ipyrad2` or one of its external binaries is “not found”, first check that the conda environment is active:

```bash
conda activate ipyrad2
```

### Installing only the Python package

`pip install .` alone does not install the external tools used by trimming, mapping, and assembly. That is the main reason the docs recommend creating the conda environment first.

### Missing analysis extras

Some analysis methods import optional packages lazily. That means the main install can succeed, but a method such as PCA or UMAP may later fail with a message about missing `scikit-learn` or `umap-learn`. Install those extras only if you need the corresponding methods.

### Using the wrong shell or Python

If `ipyrad2` starts but behaves like a different installation than expected, check which executable is being used:

```bash
which ipyrad2
```

That is especially useful if you have multiple conda environments or both editable and non-editable installs on the same machine.

## Where to Go Next

- Read [What Is ipyrad?](./what-is-ipyrad.md) for the short project overview.
- Read [Files and Data Types](./files-and-data-types.md) for the workflow file model.
- Read [Assembly](../assembly/index.md) for the end-to-end workflow.
- Read the [Analysis Guide](../analyses/index.md) if you plan to use the downstream analysis commands.
