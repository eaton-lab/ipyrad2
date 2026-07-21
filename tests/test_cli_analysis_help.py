import argparse

import pytest

from ipyrad2.cli.cli_main import setup_parsers


def _get_top_level_subparsers() -> argparse._SubParsersAction:
    parser = setup_parsers()
    return next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )


def _get_tool_parser(tool: str) -> argparse.ArgumentParser:
    return _get_top_level_subparsers().choices[tool]


def test_top_level_help_groups_commands_and_examples_are_updated() -> None:
    parser = setup_parsers()
    help_text = parser.format_help()
    subparsers = _get_top_level_subparsers()

    expected_sections = [
        "assembly subcommands",
        "data export/conversion subcommands",
        "analysis subcommands",
        "options:",
        "Note",
        "Assembly pipeline",
        "Export/Conversion examples",
        "Analysis examples",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    expected_commands = [
        "demux",
        "trim",
        "denovo",
        "map",
        "assemble",
        "inspect",
        "wex",
        "lex",
        "snpex",
        "vcf2hdf5",
        "pca",
        "dapc",
        "snmf",
        "admixture",
        "popgen",
        "bpp",
        "baba",
        "treeslider",
    ]
    start = help_text.index("assembly subcommands")
    indices = []
    for item in expected_commands:
        idx = help_text.index(f"    {item}", start)
        indices.append(idx)
        start = idx + 1
    assert indices == sorted(indices)

    assert "analysis" not in subparsers.choices
    assert "vcf2hdf5" in subparsers.choices
    assert "baba" in subparsers.choices
    assert "treeslider" in subparsers.choices
    assert "$ ipyrad2 demux -h" in help_text
    assert "$ ipyrad2 wex -d OUT/HDF5 -m 10" in help_text
    assert "$ ipyrad2 pca -d OUT/HDF5 -i IMAP -g MINMAP -I sample --plot" in help_text
    assert "ipyrad2 analysis" not in help_text
    assert "vcf-to-hdf5" not in help_text


def test_wex_help_groups_examples_and_formats_are_updated() -> None:
    help_text = _get_tool_parser("wex").format_help()

    expected_sections = [
        "Core inputs:",
        "Locus sampling:",
        "Filtering and samples:",
        "Output control:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    expected_order = [
        "-d, --data",
        "-n, --name",
        "-o, --out",
        "-O, --out-format",
        "-w, --windows",
        "-m, --min-sample-coverage",
        "-r, --max-sample-missing",
        "-e, --exclude",
        "-R, --include-reference",
        "-i, --imap",
        "-g, --minmap",
        "-P, --print-scaffold-table",
        "-x, --stdout",
        "-f, --force",
        "-l, --log-level",
        "-h, --help",
    ]
    start = help_text.index("Core inputs:")
    indices = []
    for item in expected_order:
        idx = help_text.index(item, start)
        indices.append(idx)
        start = idx + 1
    assert indices == sorted(indices)

    assert "ipyrad2 wex: extract one alignment from selected genomic windows" in help_text
    assert "$ ipyrad2 wex -d HDF5 -o OUT/" in help_text
    assert "$ ipyrad2 wex -d HDF5 -o OUT/ -w windows.bed -O fa" in help_text
    assert "If omitted, wex selects" in help_text
    assert "the full length of all scaffolds." in help_text
    assert "BED uses standard" in help_text
    assert "0-based half-open coordinates" in help_text
    assert "--include-reference" in help_text
    assert "assembly_reference_sequence" in help_text
    assert "ipyrad2 analysis wex" not in help_text
    assert "inex" not in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_lex_help_groups_examples_and_logging_are_updated() -> None:
    help_text = _get_tool_parser("lex").format_help()

    expected_sections = [
        "Core inputs:",
        "Locus sampling:",
        "Filtering and samples:",
        "Output control:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    expected_order = [
        "-d, --data",
        "-n, --name",
        "-o, --out",
        "-O, --out-format",
        "-w, --windows",
        "-N, --max-loci",
        "-s, --random-seed",
        "-L, --min-length",
        "-m, --min-sample-coverage",
        "-r, --max-sample-missing",
        "-e, --exclude",
        "-R, --include-reference",
        "-i, --imap",
        "-g, --minmap",
        "-C, --concatenate",
        "-P, --print-scaffold-table",
        "-x, --stdout",
        "-f, --force",
        "-l, --log-level",
        "-h, --help",
    ]
    start = help_text.index("Core inputs:")
    indices = []
    for item in expected_order:
        idx = help_text.index(item, start)
        indices.append(idx)
        start = idx + 1
    assert indices == sorted(indices)

    assert "ipyrad2 lex: extract delimited loci from HDF5 database" in help_text
    assert "$ ipyrad2 lex -d assembly.hdf5 -o OUT/ -w Chr1:1-50000 -N 25 -L 300 -O bpp" in help_text
    assert "-s 123 -C -O phy" in help_text
    assert "Append selected loci end to end" in help_text
    assert "Minimum whole-locus length in bp before and after filtering" in help_text
    assert "--include-reference" in help_text
    assert "assembly_reference_sequence" in help_text
    assert "-l, --log-level" in help_text
    assert "--log-level" in help_text
    assert "ipyrad2 analysis lex" not in help_text
    assert "--nloci" not in help_text
    assert "--length" not in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_treeslider_help_groups_examples_and_tree_options_are_present() -> None:
    help_text = _get_tool_parser("treeslider").format_help()

    expected_sections = [
        "Core inputs:",
        "Window planning:",
        "Filtering and samples:",
        "Tree inference:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    expected_order = [
        "-d, --data",
        "-n, --name",
        "-o, --out",
        "--window-size",
        "--slide-size",
        "--scaffolds",
        "-P, --print-scaffold-table",
        "-m, --min-sample-coverage",
        "--min-sample-alignment-length",
        "--min-alignment-length",
        "-e, --exclude",
        "-R, --include-reference",
        "-i, --imap",
        "-g, --minmap",
        "-j, --jobs",
        "--threads",
        "--workers",
        "--bs-trees",
        "--model",
        "--raxml-ng-binary",
        "--seed",
        "--redo",
        "-f, --force",
        "-l, --log-level",
        "-h, --help",
    ]
    start = help_text.index("Core inputs:")
    indices = []
    for item in expected_order:
        idx = help_text.index(item, start)
        indices.append(idx)
        start = idx + 1
    assert indices == sorted(indices)

    assert "ipyrad2 treeslider: extract filtered windows and infer one tree per window" in help_text
    assert "$ ipyrad2 treeslider -d assembly.hdf5 --print-scaffold-table" in help_text
    assert "$ ipyrad2 treeslider -d assembly.hdf5 -o OUT/ -n windows --window-size 100000 --slide-size 50000" in help_text
    assert "one-tree-per-locus mode" in help_text
    assert "shell-style wildcard patterns" in help_text
    assert "assembly_reference_sequence" in help_text
    assert "filtering and alignment-writing jobs" in help_text
    assert "--raxml-ng-binary" in help_text
    assert "--redo" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_snpex_help_groups_examples_and_reference_controls_are_present() -> None:
    help_text = _get_tool_parser("snpex").format_help()

    expected_sections = [
        "Core inputs:",
        "Filtering and samples:",
        "Linkage and export:",
        "Performance and overwrite:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    expected_order = [
        "-d, --data",
        "-n, --name",
        "-o, --out",
        "-m, --min-sample-coverage",
        "-r, --max-sample-missing",
        "-a, --min-minor-allele-frequency",
        "--min-genotype-depth",
        "--min-site-qual",
        "-e, --exclude",
        "-R, --include-reference",
        "-i, --imap",
        "-g, --minmap",
        "--no-subsample",
        "--seed",
        "--plink",
        "--phylip",
        "--nexus",
        "--fasta",
        "--treemix",
        "--eems",
        "--impute-method",
        "-c, --cores",
        "-f, --force",
        "-l, --log-level",
        "-h, --help",
    ]
    start = help_text.index("Core inputs:")
    indices = []
    for item in expected_order:
        idx = help_text.index(item, start)
        indices.append(idx)
        start = idx + 1
    assert indices == sorted(indices)

    assert "ipyrad2 snpex: extract filtered SNP matrices from HDF5 database" in help_text
    assert "$ ipyrad2 snpex -d assembly.hdf5 -o OUT/" in help_text
    assert "$ ipyrad2 snpex -d assembly.hdf5 -o OUT/ -n SNPSET --seed 123 --phylip" in help_text
    assert "$ ipyrad2 snpex -d assembly.hdf5 -o OUT/ --no-subsample --plink --impute-method sample" in help_text
    assert "$ ipyrad2 snpex -d assembly.hdf5 -o OUT/ -i POPs.txt -g MINs.txt" in help_text
    assert "-I" in help_text
    assert "Population-to-minimum-coverage mapping file" in help_text
    assert "`population<TAB>min`" in help_text
    assert "--min-genotype-depth" in help_text
    assert "--min-site-qual" in help_text
    assert "--include-reference" in help_text
    assert "assembly_reference_sequence" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_vcf2hdf5_help_groups_examples_and_logging_are_present() -> None:
    help_text = _get_tool_parser("vcf2hdf5").format_help()

    expected_sections = [
        "Core inputs:",
        "Conversion:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    expected_order = [
        "-d, --data",
        "-n, --name",
        "-o, --out",
        "-b, --ld-block-size",
        "-f, --force",
        "-l, --log-level",
        "-h, --help",
    ]
    start = help_text.index("Core inputs:")
    indices = []
    for item in expected_order:
        idx = help_text.index(item, start)
        indices.append(idx)
        start = idx + 1
    assert indices == sorted(indices)

    assert "ipyrad2 vcf2hdf5: convert VCF to SNP-capable HDF5 database" in help_text
    assert "$ ipyrad2 vcf2hdf5 -d variants.vcf.gz -o OUT/" in help_text
    assert "output-vcf2hdf5" in help_text
    assert "vcf-to-hdf5" not in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_pca_help_groups_examples_and_imputation_are_present() -> None:
    help_text = _get_tool_parser("pca").format_help()

    expected_sections = [
        "Core inputs:",
        "Method and linkage:",
        "Plotting:",
        "Filtering and samples:",
        "Performance and overwrite:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    assert "ipyrad2 pca: run PCA, t-SNE, or UMAP on SNP HDF5 data" in help_text
    assert "-I" in help_text
    assert "--impute-method" in help_text
    assert "sample` or `zero" in help_text
    assert "zero-fill" not in help_text
    assert "deprecated alias" not in help_text
    assert "--no-subsample" in help_text
    assert "-M str, --method str" in help_text
    assert "Method to run: pca, tsne, or umap." in help_text
    assert "serial method initialization" in help_text
    assert "--plot" in help_text
    assert "--plot-width" in help_text
    assert "--plot-height" in help_text
    assert "--plot-marker-size" in help_text
    assert "--plot-colors" in help_text
    assert "population color file" in help_text
    assert "Population-to-minimum-coverage mapping file" in help_text
    assert "`population<TAB>min`" in help_text
    assert "--min-genotype-depth" in help_text
    assert "--min-site-qual" in help_text
    assert "chunked SNP filtering and UMAP embedding" in help_text
    assert "$ ipyrad2 pca -d snps.hdf5 -o OUT/ --plot" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_snmf_help_groups_examples_and_imputation_are_present() -> None:
    help_text = _get_tool_parser("snmf").format_help()

    expected_sections = [
        "Core inputs:",
        "Clustering:",
        "Regularization and scoring:",
        "Filtering and samples:",
        "Performance and overwrite:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    assert "ipyrad2 snmf: run sNMF-style clustering on SNP HDF5 data" in help_text
    assert "-I" in help_text
    assert "--impute-method" in help_text
    assert "--k-range" in help_text
    assert "cross-entropy" in help_text
    assert "--alpha-w" in help_text
    assert "--alpha-h" in help_text
    assert "--l1-ratio" in help_text
    assert "--n-init" in help_text
    assert "--cv-replicates" in help_text
    assert "--cv-holdout" in help_text
    assert "Population-to-minimum-coverage mapping file" in help_text
    assert "`population<TAB>min`" in help_text
    assert "--min-genotype-depth" in help_text
    assert "--min-site-qual" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_dapc_help_groups_examples_and_imputation_are_present() -> None:
    help_text = _get_tool_parser("dapc").format_help()

    expected_sections = [
        "Core inputs:",
        "Clustering:",
        "Filtering and samples:",
        "Performance and overwrite:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    assert "ipyrad2 dapc: run DAPC-style clustering on SNP HDF5 data" in help_text
    assert "-I" in help_text
    assert "--impute-method" in help_text
    assert "--n-pcs" in help_text
    assert "Population-to-minimum-coverage mapping file" in help_text
    assert "`population<TAB>min`" in help_text
    assert "--min-genotype-depth" in help_text
    assert "--min-site-qual" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_admixture_help_groups_examples_and_imputation_are_present() -> None:
    help_text = _get_tool_parser("admixture").format_help()

    expected_sections = [
        "Core inputs:",
        "Clustering:",
        "Filtering and samples:",
        "Performance and overwrite:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    assert "ipyrad2 admixture: run external ADMIXTURE on SNP HDF5 data" in help_text
    assert "-I" in help_text
    assert "--impute-method" in help_text
    assert "--k-range" in help_text
    assert "--binary" in help_text
    assert "--keep-intermediates" in help_text
    assert "Population-to-minimum-coverage mapping file" in help_text
    assert "`population<TAB>min`" in help_text
    assert "--min-genotype-depth" in help_text
    assert "--min-site-qual" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_popgen_help_groups_examples_and_backend_controls_are_present() -> None:
    help_text = _get_tool_parser("popgen").format_help()

    expected_sections = [
        "Core inputs:",
        "Filtering and samples:",
        "SNP-backed options:",
        "Windowing:",
        "Performance and overwrite:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    assert "ipyrad2 popgen: compute genome-wide population-genetic statistics" in help_text
    assert "--stats" in help_text
    assert "--subsample-unlinked" in help_text
    assert "--seed" in help_text
    assert "--include-reference" in help_text
    assert "fis" in help_text
    assert "fit" in help_text
    assert "--window-size" in help_text
    assert "--loci-per-window" in help_text
    assert "Sequence HDF5 supports the full panel" in help_text
    assert "population-aware coverage filtering" in help_text
    assert "global `-m` filter" in help_text
    assert "--min-genotype-depth" in help_text
    assert "--min-site-qual" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_baba_help_groups_examples_and_tree_controls_are_present() -> None:
    help_text = _get_tool_parser("baba").format_help()

    expected_sections = [
        "Core inputs:",
        "Filtering and samples:",
        "Statistics and resampling:",
        "Optional outputs:",
        "Performance and overwrite:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    assert "ipyrad2 baba: compute ABBA/BABA admixture metrics from SNP HDF5 data" in help_text
    assert "--tests" in help_text
    assert "--tree" in help_text
    assert "--resampling" in help_text
    assert "--bootstrap-replicates" in help_text
    assert "--jackknife-block-bp" in help_text
    assert "--f-branch" in help_text
    assert "--f-branch-p-threshold" in help_text
    assert "--write-block-table" in help_text
    assert "--clustering-stats" in help_text
    assert "requires `--imap`" in help_text
    assert "Input SNP HDF5 file containing `genos` and `snpsmap`." in help_text
    assert "Population-to-minimum-coverage file; applied per quartet and requires `--imap`." in help_text
    assert "--min-genotype-depth" in help_text
    assert "--min-site-qual" in help_text
    assert "Significance method; `auto` prefers physical-block jackknife." in help_text
    assert "P-value threshold for zeroing non-significant tree f_G values." in help_text
    assert "$ ipyrad2 baba -d assembly.hdf5 -o OUT/ --tests quartets.tsv" in help_text
    assert "$ ipyrad2 baba -d assembly.hdf5 -o OUT/ --tree species.nwk --f-branch" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_bpp_help_groups_examples_and_runtime_controls_are_present() -> None:
    help_text = _get_tool_parser("bpp").format_help()

    expected_sections = [
        "Core inputs:",
        "Locus sampling:",
        "Model selection:",
        "Priors:",
        "Rate variation:",
        "MCMC sampler:",
        "Runtime:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    expected_order = [
        "-d, --data",
        "-o, --out",
        "-n, --name",
        "--tree",
        "-i, --imap",
        "-g, --minmap",
        "-N, --max-loci",
        "-L, --min-length",
        "--msc-i",
        "--msc-m",
        "--speciestree",
        "--speciesdelimitation",
        "--thetaprior",
        "--tauprior",
        "--speciesmodelprior",
        "--phiprior",
        "--wprior",
        "--alphaprior",
        "--locusrate",
        "--clock",
        "--burnin",
        "--samplefreq",
        "--nsample",
        "--threads",
        "--seed",
        "--write-only",
        "-f, --force",
        "-l, --log-level",
        "-h, --help",
    ]
    start = help_text.index("Core inputs:")
    indices = []
    for item in expected_order:
        idx = help_text.index(item, start)
        indices.append(idx)
        start = idx + 1
    assert indices == sorted(indices)

    assert "ipyrad2 bpp: stage one BPP analysis from sequence HDF5 data" in help_text
    assert "$ ipyrad2 bpp -d assembly.hdf5 -o OUT/ -n demo --tree species.nwk -i IMAP.txt -g MINMAP.txt -N 1000 -L 100 --threads 50 --seed 123 --tauprior 3 0.05" in help_text
    assert "--msc-m" in help_text
    assert "--speciesdelimitation" in help_text
    assert "--alphaprior" in help_text
    assert "--locusrate" in help_text
    assert "--clock" in help_text
    assert "--write-only" in help_text
    assert "sample<TAB>population" in help_text
    assert "glob<TAB>population" in help_text
    assert "population<TAB>min" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_top_level_parser_accepts_flattened_export_and_analysis_subcommands() -> None:
    parser = setup_parsers()

    lex_args = parser.parse_args(
        [
            "lex",
            "-d",
            "assembly.hdf5",
            "-o",
            "OUT",
            "-N",
            "10",
            "-L",
            "150",
        ]
    )
    assert lex_args.subcommand == "lex"
    assert lex_args.max_loci == 10
    assert lex_args.min_length == 150

    snpex_args = parser.parse_args(
        [
            "snpex",
            "-d",
            "assembly.hdf5",
            "-o",
            "OUT",
            "-n",
            "snpset",
        ]
    )
    assert snpex_args.subcommand == "snpex"
    assert snpex_args.name == "snpset"
    assert snpex_args.impute_method is None

    snpex_impute_args = parser.parse_args(
        [
            "snpex",
            "-d",
            "assembly.hdf5",
            "-o",
            "OUT",
            "-I",
            "sample",
        ]
    )
    assert snpex_impute_args.impute_method == "sample"

    vcf_args = parser.parse_args(
        [
            "vcf2hdf5",
            "-d",
            "variants.vcf.gz",
            "-o",
            "OUT",
        ]
    )
    assert vcf_args.subcommand == "vcf2hdf5"
    assert str(vcf_args.out) == "OUT"

    pca_args = parser.parse_args(
        [
            "pca",
            "-d",
            "snps.hdf5",
            "-o",
            "OUT",
        ]
    )
    assert pca_args.subcommand == "pca"
    assert pca_args.impute_method == "sample"
    assert pca_args.plot is False
    assert pca_args.plot_width == 400
    assert pca_args.plot_height == 300
    assert pca_args.plot_marker_size == 10

    pca_plot_args = parser.parse_args(
        [
            "pca",
            "-d",
            "snps.hdf5",
            "-o",
            "OUT",
            "--plot",
            "-I",
            "zero",
        ]
    )
    assert pca_plot_args.subcommand == "pca"
    assert pca_plot_args.plot is True
    assert pca_plot_args.impute_method == "zero"

    snmf_args = parser.parse_args(
        [
            "snmf",
            "-d",
            "snps.hdf5",
            "-o",
            "OUT",
            "-k",
            "2",
            "-I",
            "none",
        ]
    )
    assert snmf_args.subcommand == "snmf"
    assert snmf_args.k == 2
    assert snmf_args.impute_method == "none"
    assert snmf_args.alpha_w == 1e-4
    assert snmf_args.alpha_h == "same"
    assert snmf_args.l1_ratio == 1.0
    assert snmf_args.n_init == 10
    assert snmf_args.cv_replicates == 5
    assert snmf_args.cv_holdout == 0.1

    dapc_args = parser.parse_args(
        [
            "dapc",
            "-d",
            "snps.hdf5",
            "-o",
            "OUT",
            "--k-range",
            "2:4",
            "-I",
            "none",
        ]
    )
    assert dapc_args.subcommand == "dapc"
    assert dapc_args.k_range == "2:4"
    assert dapc_args.impute_method == "none"

    admixture_args = parser.parse_args(
        [
            "admixture",
            "-d",
            "snps.hdf5",
            "-o",
            "OUT",
            "-k",
            "2",
            "-I",
            "none",
        ]
    )
    assert admixture_args.subcommand == "admixture"
    assert admixture_args.k == 2
    assert admixture_args.impute_method == "none"

    popgen_args = parser.parse_args(
        [
            "popgen",
            "-d",
            "assembly.hdf5",
            "-o",
            "OUT",
            "--stats",
            "pi,dxy,fst",
        ]
    )
    assert popgen_args.subcommand == "popgen"
    assert popgen_args.stats == "pi,dxy,fst"

    popgen_window_args = parser.parse_args(
        [
            "popgen",
            "-d",
            "assembly.hdf5",
            "-o",
            "OUT",
            "--stats",
            "pi,fis",
            "--window-size",
            "1000",
            "--step-size",
            "500",
        ]
    )
    assert popgen_window_args.subcommand == "popgen"
    assert popgen_window_args.window_size == 1000
    assert popgen_window_args.step_size == 500

    baba_args = parser.parse_args(
        [
            "baba",
            "-d",
            "snps.hdf5",
            "-o",
            "OUT",
            "--tests",
            "quartets.tsv",
            "--resampling",
            "bootstrap",
            "--bootstrap-replicates",
            "100",
            "--f-branch-p-threshold",
            "0.02",
            "--seed",
            "9",
        ]
    )
    assert baba_args.subcommand == "baba"
    assert str(baba_args.tests) == "quartets.tsv"
    assert baba_args.resampling == "bootstrap"
    assert baba_args.bootstrap_replicates == 100
    assert baba_args.f_branch_p_threshold == 0.02
    assert baba_args.seed == 9

    bpp_args = parser.parse_args(
        [
            "bpp",
            "-d",
            "assembly.hdf5",
            "-o",
            "OUT",
            "--tree",
            "(a,b);",
            "-i",
            "IMAP.txt",
            "-g",
            "MINMAP.txt",
            "-N",
            "25",
            "-L",
            "100",
            "--threads",
            "8",
            "--seed",
            "none",
            "--write-only",
        ]
    )
    assert bpp_args.subcommand == "bpp"
    assert bpp_args.max_loci == 25
    assert bpp_args.min_length == 100
    assert bpp_args.threads == [8]
    assert bpp_args.seed == "none"
    assert bpp_args.write_only is True


def test_removed_nested_analysis_command_no_longer_parses() -> None:
    parser = setup_parsers()

    with pytest.raises(SystemExit):
        parser.parse_args(["analysis", "wex", "-d", "assembly.hdf5"])


def test_treeslider_is_a_real_top_level_command_parser() -> None:
    parser = setup_parsers()

    tool = "treeslider"
    help_text = _get_tool_parser(tool).format_help()
    assert f"ipyrad2 {tool}:" in help_text
    assert "Reserved command placeholder." not in help_text
    assert "not implemented yet" not in help_text

    args = parser.parse_args([tool, "-d", "assembly.hdf5"])
    assert args.subcommand == tool
    assert args.data.name == "assembly.hdf5"
    assert args.jobs == 4

    parallel_args = parser.parse_args([tool, "-d", "assembly.hdf5", "-j", "7"])
    assert parallel_args.jobs == 7

    with pytest.raises(SystemExit):
        parser.parse_args(
            [tool, "-d", "assembly.hdf5", "--redo", "--force"]
        )
