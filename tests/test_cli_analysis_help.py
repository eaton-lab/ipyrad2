import argparse

from ipyrad2.cli.cli_main import setup_parsers


def _get_analysis_tool_parser(tool: str) -> argparse.ArgumentParser:
    parser = setup_parsers()
    subparsers = next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    analysis_parser = subparsers.choices["analysis"]
    analysis_subparsers = next(
        action for action in analysis_parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    return analysis_subparsers.choices[tool]


def test_wex_help_groups_examples_and_formats_are_updated() -> None:
    help_text = _get_analysis_tool_parser("wex").format_help()

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

    assert "ipyrad2 analysis wex: extract one alignment from selected genomic windows" in help_text
    assert "$ ipyrad2 analysis wex -d HDF5 -o OUT/" in help_text
    assert "$ ipyrad2 analysis wex -d HDF5 -o OUT/ -w windows.bed -O fa" in help_text
    assert "If omitted, wex selects" in help_text
    assert "the full length of all scaffolds." in help_text
    assert "BED uses standard" in help_text
    assert "0-based half-open coordinates" in help_text
    assert "--include-reference" in help_text
    assert "assembly_reference_sequence" in help_text
    assert "ipyrad analysis wex" not in help_text
    assert "inex" not in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_lex_help_groups_examples_and_logging_are_updated() -> None:
    help_text = _get_analysis_tool_parser("lex").format_help()

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
        "-L, --min-length",
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

    assert "ipyrad2 analysis lex: extract delimited loci from HDF5 database" in help_text
    assert "$ ipyrad2 analysis lex -d assembly.hdf5 -o OUT/ -w Chr1:1-50000 -N 25 -L 300 -O bpp" in help_text
    assert "Minimum whole-locus length in bp before and after filtering" in help_text
    assert "--include-reference" in help_text
    assert "assembly_reference_sequence" in help_text
    assert "-l, --log-level" in help_text
    assert "--log-level" in help_text
    assert "ipyrad analysis lex" not in help_text
    assert "--nloci" not in help_text
    assert "--length" not in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_snpex_help_groups_examples_and_reference_controls_are_present() -> None:
    help_text = _get_analysis_tool_parser("snpex").format_help()

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

    assert "ipyrad2 analysis snpex: extract filtered SNP matrices from HDF5 database" in help_text
    assert "$ ipyrad2 analysis snpex -d assembly.hdf5 -o OUT/" in help_text
    assert "$ ipyrad2 analysis snpex -d assembly.hdf5 -o OUT/ -n SNPSET --seed 123 --phylip" in help_text
    assert "$ ipyrad2 analysis snpex -d assembly.hdf5 -o OUT/ --no-subsample --plink --impute-method sample" in help_text
    assert "$ ipyrad2 analysis snpex -d assembly.hdf5 -o OUT/ -i POPs.txt -g MINs.txt" in help_text
    assert "-I" in help_text
    assert "Population-to-minimum-coverage mapping file" in help_text
    assert "`population<TAB>min`" in help_text
    assert "--include-reference" in help_text
    assert "assembly_reference_sequence" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_vcf_to_hdf5_help_groups_examples_and_logging_are_present() -> None:
    help_text = _get_analysis_tool_parser("vcf-to-hdf5").format_help()

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

    assert "ipyrad2 analysis vcf-to-hdf5: convert VCF to SNP-capable HDF5 database" in help_text
    assert "$ ipyrad2 analysis vcf-to-hdf5 -d variants.vcf.gz -o OUT/" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_pca_help_groups_examples_and_imputation_are_present() -> None:
    help_text = _get_analysis_tool_parser("pca").format_help()

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

    assert "ipyrad2 analysis pca: run PCA, t-SNE, or UMAP on SNP HDF5 data" in help_text
    assert "-I" in help_text
    assert "--impute-method" in help_text
    assert "sample` or `zero" in help_text
    assert "zero-fill" not in help_text
    assert "deprecated alias" not in help_text
    assert "--no-subsample" in help_text
    assert "-M, --method" in help_text
    assert "Method to run: pca, tsne, or umap." in help_text
    assert "--plot" in help_text
    assert "--plot-width" in help_text
    assert "--plot-height" in help_text
    assert "--plot-marker-size" in help_text
    assert "Population-to-minimum-coverage mapping file" in help_text
    assert "`population<TAB>min`" in help_text
    assert "$ ipyrad2 analysis pca -d snps.hdf5 -o OUT/ --plot" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_snmf_help_groups_examples_and_imputation_are_present() -> None:
    help_text = _get_analysis_tool_parser("snmf").format_help()

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

    assert "ipyrad2 analysis snmf: run sNMF-style clustering on SNP HDF5 data" in help_text
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
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_dapc_help_groups_examples_and_imputation_are_present() -> None:
    help_text = _get_analysis_tool_parser("dapc").format_help()

    expected_sections = [
        "Core inputs:",
        "Clustering:",
        "Filtering and samples:",
        "Performance and overwrite:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    assert "ipyrad2 analysis dapc: run DAPC-style clustering on SNP HDF5 data" in help_text
    assert "-I" in help_text
    assert "--impute-method" in help_text
    assert "--n-pcs" in help_text
    assert "Population-to-minimum-coverage mapping file" in help_text
    assert "`population<TAB>min`" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_admixture_help_groups_examples_and_imputation_are_present() -> None:
    help_text = _get_analysis_tool_parser("admixture").format_help()

    expected_sections = [
        "Core inputs:",
        "Clustering:",
        "Filtering and samples:",
        "Performance and overwrite:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    assert "ipyrad2 analysis admixture: run external ADMIXTURE on SNP HDF5 data" in help_text
    assert "-I" in help_text
    assert "--impute-method" in help_text
    assert "--k-range" in help_text
    assert "--binary" in help_text
    assert "--keep-intermediates" in help_text
    assert "Population-to-minimum-coverage mapping file" in help_text
    assert "`population<TAB>min`" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_popgen_help_groups_examples_and_backend_controls_are_present() -> None:
    help_text = _get_analysis_tool_parser("popgen").format_help()

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

    assert "ipyrad2 analysis popgen: compute genome-wide population-genetic statistics" in help_text
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
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_bpp_help_groups_examples_and_runtime_controls_are_present() -> None:
    help_text = _get_analysis_tool_parser("bpp").format_help()

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

    assert "ipyrad2 analysis bpp: stage one BPP analysis from sequence HDF5 data" in help_text
    assert "$ ipyrad2 analysis bpp -d assembly.hdf5 -o OUT/ -n demo --tree species.nwk -i IMAP.txt -g MINMAP.txt -N 1000 -L 100 --threads 50 --seed 123 --tauprior 3 0.05" in help_text
    assert "--msc-m" in help_text
    assert "--speciesdelimitation" in help_text
    assert "--alphaprior" in help_text
    assert "--locusrate" in help_text
    assert "--clock" in help_text
    assert "--write-only" in help_text
    assert "sample<TAB>population" in help_text
    assert "population<TAB>min" in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_analysis_parser_accepts_phase2_subcommands() -> None:
    parser = setup_parsers()

    lex_args = parser.parse_args(
        [
            "analysis",
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
    assert lex_args.tool == "lex"
    assert lex_args.max_loci == 10
    assert lex_args.min_length == 150

    snpex_args = parser.parse_args(
        [
            "analysis",
            "snpex",
            "-d",
            "assembly.hdf5",
            "-o",
            "OUT",
            "-n",
            "snpset",
        ]
    )
    assert snpex_args.subcommand == "analysis"
    assert snpex_args.tool == "snpex"
    assert snpex_args.name == "snpset"
    assert snpex_args.impute_method is None

    snpex_impute_args = parser.parse_args(
        [
            "analysis",
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
            "analysis",
            "vcf-to-hdf5",
            "-d",
            "variants.vcf.gz",
            "-o",
            "OUT",
        ]
    )
    assert vcf_args.subcommand == "analysis"
    assert vcf_args.tool == "vcf-to-hdf5"
    assert str(vcf_args.out) == "OUT"

    pca_args = parser.parse_args(
        [
            "analysis",
            "pca",
            "-d",
            "snps.hdf5",
            "-o",
            "OUT",
        ]
    )
    assert pca_args.tool == "pca"
    assert pca_args.impute_method == "sample"
    assert pca_args.plot is False
    assert pca_args.plot_width == 400
    assert pca_args.plot_height == 300
    assert pca_args.plot_marker_size == 10

    pca_plot_args = parser.parse_args(
        [
            "analysis",
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
    assert pca_plot_args.tool == "pca"
    assert pca_plot_args.plot is True
    assert pca_plot_args.impute_method == "zero"

    snmf_args = parser.parse_args(
        [
            "analysis",
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
    assert snmf_args.tool == "snmf"
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
            "analysis",
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
    assert dapc_args.tool == "dapc"
    assert dapc_args.k_range == "2:4"
    assert dapc_args.impute_method == "none"

    admixture_args = parser.parse_args(
        [
            "analysis",
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
    assert admixture_args.tool == "admixture"
    assert admixture_args.k == 2
    assert admixture_args.impute_method == "none"

    popgen_args = parser.parse_args(
        [
            "analysis",
            "popgen",
            "-d",
            "assembly.hdf5",
            "-o",
            "OUT",
            "--stats",
            "pi,dxy,fst",
        ]
    )
    assert popgen_args.tool == "popgen"
    assert popgen_args.stats == "pi,dxy,fst"

    popgen_window_args = parser.parse_args(
        [
            "analysis",
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
    assert popgen_window_args.tool == "popgen"
    assert popgen_window_args.window_size == 1000
    assert popgen_window_args.step_size == 500

    bpp_args = parser.parse_args(
        [
            "analysis",
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
    assert bpp_args.tool == "bpp"
    assert bpp_args.max_loci == 25
    assert bpp_args.min_length == 100
    assert bpp_args.threads == [8]
    assert bpp_args.seed == "none"
    assert bpp_args.write_only is True
