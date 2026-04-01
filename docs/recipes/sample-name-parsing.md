# Using -dx and -di to pair and name samples

This recipe shows how to use `-dx` and `-di` when FASTQ filenames contain extra provider, lane, or run tokens that you do not want in the final sample name.

The example command below uses `ipyrad2 trim`, but the same delimiter logic applies to other ipyrad2 commands that accept `-dx` and `-di`, including `map` and `denovo`.

## The Rule

`-dx` and `-di` work by splitting the filename on a literal delimiter string and then keeping the text left of the Nth occurrence of that delimiter.

In practical terms, ipyrad2 does this:

```python
parts = filename.split(delim)
sample_name = delim.join(parts[:delim_index])
```

That means:

- `-dx` chooses the literal split string
- `-di 1` keeps the text left of the first occurrence
- `-di 2` keeps the text left of the second occurrence
- and so on

This affects two things at once:

- which files are grouped together as one sample
- what sample name is written into downstream outputs

## Example 1: Drop the provider sample-sheet token

Suppose your paired-end files look like this:

```text
AMA_21040XD-01-07_S39_L002_R1.fastq.gz
AMA_21040XD-01-07_S39_L002_R2.fastq.gz
AMA_21040XD-01-08_S40_L002_R1.fastq.gz
AMA_21040XD-01-08_S40_L002_R2.fastq.gz
```

Here, the biological sample names you want are:

```text
AMA_21040XD-01-07
AMA_21040XD-01-08
```

The `_S39` and `_S40` tokens come from the sequencing sample sheet, and `_L002` is a lane token. They are useful metadata, but not part of the sample identity you want in the trimmed outputs.

Use:

```bash
ipyrad2 trim -d FASTQs/*.fastq.gz -o TRIMMED/ -dx _S -di 1
```

Why this works:

```python
"AMA_21040XD-01-07_S39_L002_R1.fastq.gz".split("_S")
# ["AMA_21040XD-01-07", "39_L002_R1.fastq.gz"]
```

With `-di 1`, ipyrad2 keeps only the text left of the first `_S`, so the parsed sample name becomes:

```text
AMA_21040XD-01-07
```

The same logic applies to the matching R2 file, so those two files are treated as one paired sample. The outputs will be written as:

```text
AMA_21040XD-01-07.R1.trimmed.fastq.gz
AMA_21040XD-01-07.R2.trimmed.fastq.gz
AMA_21040XD-01-08.R1.trimmed.fastq.gz
AMA_21040XD-01-08.R2.trimmed.fastq.gz
```

## Example 2: Keep more text by changing -di

Now suppose your filenames look like this:

```text
Canis.lupus.popA.rep1.R1.fastq.gz
Canis.lupus.popA.rep1.R2.fastq.gz
```

If you split on `.`:

```python
"Canis.lupus.popA.rep1.R1.fastq.gz".split(".")
# ["Canis", "lupus", "popA", "rep1", "R1", "fastq", "gz"]
```

Then:

- `-dx . -di 2` gives `Canis.lupus`
- `-dx . -di 3` gives `Canis.lupus.popA`

So these two commands would keep different sample names:

```bash
ipyrad2 trim -d FASTQs/*.fastq.gz -o TRIMMED/ -dx . -di 2
```

```bash
ipyrad2 trim -d FASTQs/*.fastq.gz -o TRIMMED/ -dx . -di 3
```

This is the main use of `-di`: you can decide exactly how much of the left side of the filename should count as the sample name.

## What This Recipe Solves

This recipe is for:

- making sure R1 and R2 files group together the way you expect
- removing provider or run tokens from output sample names
- keeping sample names short, consistent, and biologically meaningful

It is not a recipe for merging multiple lanes or several separate paired files into one sample during `trim`. If one biological sample appears in several independent R1/R2 file pairs, do not assume that `-dx` and `-di` will merge them into one trimmed sample automatically.

## Troubleshooting

If files are grouped incorrectly or sample names are not what you expected:

- print a few filenames and test the split logic manually
- choose a delimiter that appears exactly where you want the sample name to stop
- adjust `-di` until the kept left-hand substring matches the sample identity
- prefer a specific delimiter such as `_S` over a more general one like `_` when possible

If you are unsure what ipyrad2 parsed, check the log output at the start of the run. Commands that use `-dx` and `-di` report the first sample names they parsed from file paths before launching the main work.

## Related Pages

- [trim](../assembly/trim.md)
- [map](../assembly/map.md)
- [Recipes](./index.md)
