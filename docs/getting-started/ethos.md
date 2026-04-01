# The ipyrad Ethos

ipyrad2 keeps the core goals of the original ipyrad project: reduced-representation genomic workflows should be simple, resourceful, reproducible, flexible, and transparent. Those ideas matter as much in ipyrad2 as they did in the earlier project. What changed is the interface: instead of centering the workflow on a parameter-file-driven step system, ipyrad2 exposes explicit command stages with clearer inputs, outputs, logs, and stats.

## Core Values

- **Simple**: the workflow should be easy to install, easy to run from the command line, and easy to understand one stage at a time.
- **Resourceful**: users should get useful logs, stats files, and concrete outputs that make it possible to diagnose mistakes instead of starting over blindly.
- **Reproducible**: the same commands on the same data should produce stable, inspectable results that can be revisited later.
- **Flexible**: the tool should support multiple reduced-representation library designs, multiple entry points into the workflow, and multiple downstream paths once loci and SNPs are assembled.
- **Transparent**: intermediate and final outputs should stay visible and inspectable rather than disappearing behind opaque internal state.

## What That Means in Practice

In ipyrad2, these values show up as a staged workflow: reads can be demultiplexed, trimmed, mapped or used to build a denovo pseudoreference, assembled into loci and variants, and then exported or analyzed directly. Logs, stats summaries, HDF5-backed data products, and method-specific result tables are all part of that design.

Users still need to choose sensible trimming, coverage, mapping, and filtering settings. The point is to make those decisions visible and manageable, while providing reasonable default values.

## Contact us
If you encounter problems running ipyrad2 the best way to get a quick response is to [raise an issue on GitHub](https://github.com/eaton-lab/ipyrad2/issues).


## Where to Go Next

- Read [What Is ipyrad2?](./what-is-ipyrad.md) for the short project overview.
- Read [Installation](./installation.md) to set up the environment.
- Read [Files and Data Types](./files-and-data-types.md) for the file model and workflow entry points.
- Read [Assembly](../assembly/index.md) for the end-to-end workflow.
