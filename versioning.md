# Versioning and Releases

`ipyrad2` uses `setuptools-scm`, so the package version comes from Git tags.
Do not edit a version string in `ipyrad2/__init__.py`, `pyproject.toml`, or
the generated `ipyrad2/_version.py` file when preparing a release.

## Tag Naming

Use these tag formats:

- Stable release: `vX.Y.Z`
- Alpha release: `vX.Y.ZaN`
- Beta release: `vX.Y.ZbN`
- Release candidate: `vX.Y.ZrcN`

Examples:

- `v0.2.0`
- `v0.3.0a1`
- `v0.3.0a2`
- `v0.3.0b1`
- `v0.3.0rc1`

Use PEP 440 style suffixes only. Do not use forms like `v0.3.0-alpha1` or
`v0.3.0-beta1`.

## What Changes the Version

The version changes when you create a new Git tag.

- Stable release: create a new `vX.Y.Z` tag.
- Pre-release: create a new `vX.Y.ZaN`, `vX.Y.ZbN`, or `vX.Y.ZrcN` tag.
- Final release after pre-releases: create a new stable tag such as `v0.3.0`.

Do not move or reuse existing tags. Each new alpha, beta, release candidate,
or final release gets its own new tag on a newer commit.

## Before You Tag

Release from `main` and from a clean worktree.

```bash
git checkout main
git pull origin main
git status --short
```

`git status --short` should print nothing before you tag. Commit any release
changes first.

This repository currently has no tags yet, so the first release can be either:

- `v0.2.0` for the first stable release
- `v0.2.0a1` for the first alpha release

## Stable Release Steps

Choose the next stable version number, then create and push an annotated tag.

```bash
git checkout main
git pull origin main
git status --short
git tag -a v0.2.0 -m "v0.2.0"
git push origin main
git push origin v0.2.0
```

## Pre-Release Steps

Choose the next alpha, beta, or release-candidate tag, then create and push it.

```bash
git checkout main
git pull origin main
git status --short
git tag -a v0.3.0a1 -m "v0.3.0a1"
git push origin main
git push origin v0.3.0a1
```

Example follow-up pre-releases:

```bash
git tag -a v0.3.0a2 -m "v0.3.0a2"
git push origin v0.3.0a2

git tag -a v0.3.0b1 -m "v0.3.0b1"
git push origin v0.3.0b1

git tag -a v0.3.0rc1 -m "v0.3.0rc1"
git push origin v0.3.0rc1
```

When the prerelease cycle is complete, create the final stable tag separately:

```bash
git tag -a v0.3.0 -m "v0.3.0"
git push origin v0.3.0
```

## After Pushing

After you push the tag:

- confirm the tag exists on GitHub
- confirm it points to the intended commit
- if you rebuild or reinstall from that tagged commit, the runtime version
  should match the tag

Example checks:

```bash
git show v0.2.0 --no-patch
python -m ipyrad2.cli.cli_main -v
```

## GitHub Action Release Behavior

The repository workflow at `.github/workflows/publish-pypi.yml` handles package
builds for release tags.

- Any pushed tag matching `v*` triggers a build.
- Stable tags matching exactly `vX.Y.Z` publish to PyPI automatically.
- Pre-release tags such as `v0.3.0a1`, `v0.3.0b1`, and `v0.3.0rc1` build but do
  not publish to PyPI.

The workflow is intended to use PyPI trusted publishing. Configure the PyPI
`ipyrad2` project to trust this GitHub repository and the publish workflow. No
PyPI API token secret is needed when trusted publishing is set up.
